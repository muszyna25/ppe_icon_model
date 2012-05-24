!>
!! Scheduler for internal post-processing.
!!
!! This module manages a "job queue" for internal post-processing
!! tasks on the compute PEs.
!!
!! For example, the interpolation of model variables onto lon-lat
!! fields constitutes such a task. Allocating and computing only those
!! fields which are required for output (or as intermediate results)
!! can save memory and computing time. Jobs are processed according to
!! user-defined priority levels.
!!
!! @author F. Prill, DWD
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2012-03-01)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
!! TODO[FP] To increase performance, one should allocate/deallocate
!!          the vertical interpolation coefficient tables "vcoeff"
!!          only when necessary.
!!
!! TODO[FP] Interpolation tasks are performed more often than
!!          necessary: The activity flag must be adjusted to the output
!!          intervals!
!!
!! TODO[FP] Do not insert post-processing tasks for patches with 0 cells.
!!
!! -----------------------------------------------------------------------------------
MODULE mo_pp_scheduler

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message, message_text, finish
  USE mo_impl_constants,          ONLY: SUCCESS,                      &
    & VINTP_TYPE_Z, VINTP_TYPE_P_OR_Z, VINTP_TYPE_NONE,               &
    & VINTP_METHOD_UV, VINTP_METHOD_LIN, HINTP_TYPE_NONE,             &     
    & VINTP_METHOD_QV, HINTP_TYPE_LONLAT, VINTP_METHOD_LIN_NLEVP1,    &
    & max_dom, max_var_ml, max_var_pl, max_var_hl
  USE mo_model_domain,            ONLY: t_patch, p_patch
  USE mo_var_list,                ONLY: add_var, get_all_var_names,         &
    &                                   create_hor_interp_metadata,         &
    &                                   nvar_lists, var_lists
  USE mo_var_list_element,        ONLY: t_var_list_element, level_type_ml,  &
    &                                   level_type_pl, level_type_hl
  USE mo_var_metadata,            ONLY: t_var_metadata, t_vert_interp_meta
  USE mo_intp,                    ONLY: verts2cells_scalar
  USE mo_intp_data_strc,          ONLY: t_int_state, lonlat_grid_list, &
    &                                   t_lon_lat_intp, p_int_state
  USE mo_nh_vert_interp,          ONLY: prepare_vert_interp, &
    &                                   lin_intp, uv_intp, qv_intp 
  USE mo_nonhydro_types,          ONLY: t_nh_state, t_nh_prog, t_nh_diag, &
    &                                   t_nh_metrics
  USE mo_nonhydro_state,          ONLY: p_nh_state
  USE mo_opt_diagnostics,         ONLY: t_nh_diag_pz, t_nh_opt_diag, t_vcoeff, &
    &                                   vcoeff_deallocate, p_nh_opt_diag
  USE mo_nwp_phy_state,           ONLY: t_nwp_phy_diag, prm_diag
  USE mo_nh_pzlev_config,         ONLY: t_nh_pzlev_config, nh_pzlev_config
  USE mo_name_list_output_config, ONLY: t_output_name_list, &
    &                                   vname_len, first_output_name_list
  USE mo_parallel_config,         ONLY: nproma
  USE mo_dynamics_config,         ONLY: nnow
  USE mo_util_string,             ONLY: int2string, remove_duplicates, &
    &                                   difference, toupper
  USE mo_cf_convention,           ONLY: t_cf_var
  USE mo_grib2,                   ONLY: t_grib2_var
  USE mo_cdi_constants,           ONLY: GRID_CELL, GRID_REFERENCE,               &
    &                                   GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,  &
    &                                   ZAXIS_PRESSURE, GRID_REGULAR_LONLAT,     &
    &                                   GRID_UNSTRUCTURED_EDGE,                  &
    &                                   GRID_UNSTRUCTURED_VERT, ZAXIS_SURFACE
  USE mo_linked_list,             ONLY: t_var_list, t_list_element
  USE mo_lonlat_grid,             ONLY: t_lon_lat_grid
  USE mo_intp_lonlat,             ONLY: rbf_interpol_lonlat_nl, &
    &                                   rbf_vec_interpol_lonlat_nl
  USE mo_sync,                    ONLY: sync_patch_array,                        &
    &                                   SYNC_C, SYNC_E, SYNC_V

  IMPLICIT NONE

  ! max. name string length
  INTEGER, PARAMETER :: MAX_NAME_LENGTH   =   64
  INTEGER, PARAMETER :: HIGH_PRIORITY     =    0  
  INTEGER, PARAMETER :: LOW_PRIORITY      =  100  
  INTEGER, PARAMETER :: DEFAULT_PRIORITY1 =    9   
  INTEGER, PARAMETER :: DEFAULT_PRIORITY2 =   10  

  ! level of output verbosity
  INTEGER :: dbg_level = 0
  
  !--- Available post-processing tasks
  !------ setup tasks (coefficients,...)
  INTEGER, PARAMETER :: TASK_INIT_VER_PZ       = 1  !< task: setup pz-interpolation
  INTEGER, PARAMETER :: TASK_INIT_VER_Z        = 2  !< task: setup only z-interpolation
  INTEGER, PARAMETER :: TASK_FINALIZE_PZ       = 3  !< task: deallocate pz-interpolation
  !------ interpolation tasks:
  INTEGER, PARAMETER :: TASK_INTP_HOR_LONLAT   = 4  !< task: lon-lat
  INTEGER, PARAMETER :: TASK_INTP_VER_PZLEV    = 5  !< task: vertical p or z-levels
  INTEGER, PARAMETER :: TASK_INTP_SYNC         = 6  !< task: synchronizes halo regions

  ! interface definition
  PRIVATE

  ! functions and subroutines
  PUBLIC :: pp_scheduler_init
  PUBLIC :: pp_scheduler_process
  PUBLIC :: pp_scheduler_finalize
  PUBLIC :: new_simulation_status
  ! data types
  PUBLIC :: t_simulation_status

  !--- JOB QUEUE DEFINITION ----------------------------------------------------------


  !> data necessary for job input.
  !
  !  This type is likely to contain more data than really needed for
  !  your specific post-processing job (we are lacking polymorphic
  !  data structures).
  !
  !  @note Please avoid using any non-static data from other places
  !  for your post-processing tasks!  The only exceptions are the use of
  !  the global, volatile "nnow" value and "nproma"!
  !
  !  @note Elements of this type are COPIED. Therefore avoid large
  !  data structures or use POINTERs.
  TYPE t_data_input
    ! pointer for model variable (array)
    TYPE (t_var_list_element), POINTER :: var

    INTEGER                            :: jg ! domain ID

    TYPE(t_patch),             POINTER :: p_patch         
    TYPE(t_int_state),         POINTER :: p_int_state     
    TYPE(t_nh_state),          POINTER :: p_nh_state      
    TYPE(t_nwp_phy_diag),      POINTER :: prm_diag        
    TYPE(t_nh_opt_diag),       POINTER :: p_nh_opt_diag   
    TYPE(t_nh_pzlev_config),   POINTER :: nh_pzlev_config 
  END TYPE t_data_input


  !> data necessary for job output.
  !  See also @p t_data_input.
  TYPE t_data_output
    ! pointer for model variable (array)
    TYPE (t_var_list_element), POINTER :: var    => NULL()

    ! (optional) pointer for second component of model variable.
    ! necessary for lon-lat interpolation of edge-based fields.
    TYPE (t_var_list_element), POINTER :: var_2  => NULL()
  END TYPE t_data_output


  !> Definition of simulation status, a list of LOGICAL flags like
  !  "first_step", "last_step", "output_time"
  !
  !  Based on these values, we determine if a post-processing
  !  task is "active" and will be processed.
  !
  !  Flags are stored in a contiguous LOGICAL array to make
  !  Boolean comparisons more convenient.
  !
  !  @note There might be better places in the code for such a
  !  variable!
  TYPE t_simulation_status
    ! l_output_step, l_first_step, l_last_step
    LOGICAL :: status_flags(3)
  END TYPE t_simulation_status


  !> A variable of type @p t_job_queue defines a single post-processing task.
  !
  !  Jobs with smaller priority values are processed first.
  TYPE t_job_queue
    
    INTEGER                         :: job_priority   !< Task priority.
    CHARACTER(len=MAX_NAME_LENGTH)  :: job_name       !< job name string (for status output)
    INTEGER                         :: job_type       !< task type (quasi function pointer)
    TYPE(t_simulation_status)       :: activity       !< "under which conditions does this task run?"

    TYPE(t_data_input)              :: data_input     !< input of post-processing task
    TYPE(t_data_output)             :: data_output    !< result of post-processing task

    TYPE(t_job_queue), POINTER      :: next => NULL() !< pointer to next element in list

  END TYPE t_job_queue


  INTERFACE assign_if_present  ! purely internal
    MODULE PROCEDURE assign_if_present_logical
  END INTERFACE


  !--- MODULE DATA -------------------------------------------------------------------
  TYPE(t_job_queue), POINTER   :: job_queue  =>  NULL() !< head of (ordered) job queue


CONTAINS

  !--- SCHEDULER ---------------------------------------------------------------------

  !---------------------------------------------------------------
  !> Setup of post-processing job queue.
  !
  ! E.g. setup of optional diagnostic quantities like pz-level
  ! interpolation. Computation of these additional variables is
  ! registered as post-processing tasks.
  !
  ! @note This subroutine must be called after the namelists have been
  !       read and the nh_state has been constructed and _before_
  !       initialization of the name list output.
  !
  SUBROUTINE pp_scheduler_init(l_init_prm_diag)
    LOGICAL, INTENT(IN) :: l_init_prm_diag

    ! local variables
    CHARACTER(*), PARAMETER :: routine =  &
      &  TRIM("mo_pp_scheduler:pp_scheduler_init")
    INTEGER                               :: &
      &  jg, ndom, ierrstat, ivar, i, j, idx, nvars_ll, nlev, &
      &  nblks_lonlat, ilev_type, max_var, ilev
    LOGICAL                               :: &
      &  l_jg_active, found, l_horintp
    TYPE (t_output_name_list), POINTER    :: p_onl
    TYPE(t_job_queue),         POINTER    :: task
    TYPE(t_var_list),          POINTER    :: p_opt_diag_list
    REAL(wp), POINTER                     :: p_opt_field_r3d(:,:,:)
    TYPE(t_list_element),      POINTER    :: element, new_element, new_element_2
    CHARACTER(LEN=vname_len),  POINTER    :: varlist(:)
    INTEGER, ALLOCATABLE                  :: ll_vargrid(:)
    CHARACTER(LEN=vname_len), ALLOCATABLE :: ll_varlist(:)
    INTEGER, ALLOCATABLE                  :: ll_varlevs(:)
    CHARACTER(LEN=vname_len)              :: vname
    TYPE(t_var_metadata),      POINTER    :: info
    INTEGER                               :: var_shape(5)
    TYPE (t_lon_lat_intp),     POINTER    :: ptr_int_lonlat

    if (dbg_level > 5)  CALL message(routine, "Enter")

    !-------------------------------------------------------------
    !--- setup of vertical interpolation onto p/z-levels

    CALL pp_scheduler_init_pz(l_init_prm_diag)

    !-------------------------------------------------------------
    !-- horizontal interpolation regular lon-lat grids

    ndom = SIZE(p_nh_opt_diag)
    ALLOCATE(ll_varlist(ndom*max_var_ml), ll_vargrid(ndom*max_var_ml), &
      &      ll_varlevs(ndom*max_var_ml), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! flag. will be set if horizontal interpolation tasks have been
    ! created
    l_horintp = .FALSE.

    ! loop over the output definitions, collect pairs of variable
    ! names/lon-lat interpolation requests

    p_onl => first_output_name_list
    nvars_ll = 0

    NML_LOOP : DO
      IF (.NOT.ASSOCIATED(p_onl)) EXIT NML_LOOP

      DO ilev=1,3
   
        SELECT CASE(ilev)
        CASE (1) 
          varlist => p_onl%ml_varlist
          ilev_type  =  level_type_ml
          max_var    =  max_var_ml
        CASE (2)
          varlist => p_onl%pl_varlist
          ilev_type  =  level_type_pl
          max_var    =  max_var_pl
        CASE (3)
          varlist => p_onl%hl_varlist
          ilev_type  =  level_type_hl
          max_var    =  max_var_hl
        END SELECT
        IF (varlist(1) == ' ') CYCLE
      
        ! Selection criterion: 
        ! - lon-lat interpolation is requested
        IF (p_onl%remap == 1) THEN

          ! check, if "all" variables are desired:
          IF (toupper(TRIM(varlist(1))) == 'ALL') THEN
            IF (dbg_level > 8)  CALL message(routine, "ALL vars for model levels")
            j = nvars_ll
            CALL get_all_var_names(ll_varlist, nvars_ll,                &
              &                    opt_hor_intp_type=HINTP_TYPE_LONLAT, &
              &                    opt_vlevel_type=ilev_type,           &
              &                    opt_loutput=.TRUE.)
            ll_vargrid((j+1):nvars_ll) = p_onl%lonlat_id
            ll_varlevs((j+1):nvars_ll) = ilev_type
            EXIT NML_LOOP
          END IF

          DO ivar=1,max_var
            IF (varlist(ivar) == ' ') CYCLE
            nvars_ll=nvars_ll+1
            ll_varlist(nvars_ll) = varlist(ivar)
            ll_vargrid(nvars_ll) = p_onl%lonlat_id
            ll_varlevs(nvars_ll) = ilev_type
          END DO
        END IF
      END DO

      p_onl => p_onl%next
    END DO NML_LOOP

    !-- loop over requested variables, add variables ("add_var") and
    !-- register interpolation tasks for all requested domains:
      
    DO ivar=1,nvars_ll
      IF (dbg_level > 8) &
        CALL message(routine, "horizontal interpolation: "&
        &     //"Looking for input var '"//TRIM(ll_varlist(ivar))//"'")
      found = .FALSE.
        
      !- loop over model level variables
      ! Note that there may be several variables with different time levels,
      ! we just add unconditionally all
      LIST_LOOP : DO i = 1,nvar_lists
        ! Do not inspect lists which are disabled for output
        IF (.NOT. var_lists(i)%p%loutput) CYCLE
        ! Do not inspect lists if vertical level type does not
        ! match (p/z/model levels):
        IF (var_lists(i)%p%vlevel_type/=ll_varlevs(ivar)) CYCLE LIST_LOOP
        ! loop only over variables on requested domains:
        jg = var_lists(i)%p%patch_id
        IF (.NOT. lonlat_grid_list(ll_vargrid(nvars_ll))%l_dom(jg)) &
          &  CYCLE LIST_LOOP
          
        element => NULL()
        VAR_LOOP : DO
          IF(.NOT.ASSOCIATED(element)) THEN
            element => var_lists(i)%p%first_list_element
          ELSE
            element => element%next_list_element
          ENDIF
          IF(.NOT.ASSOCIATED(element)) EXIT
            
          info => element%field%info
          ! Do not inspect element if it is a container
          IF (info%lcontainer) CYCLE VAR_LOOP
          ! Do not inspect element if "loutput=.false."
          IF (.NOT. info%loutput) CYCLE VAR_LOOP
          ! Do not inspect element if it does not support horizontal
          ! interpolation
          IF (info%hor_interp%hor_intp_type==HINTP_TYPE_NONE) CYCLE VAR_LOOP
          ! Check for matching name
          vname = ll_varlist(ivar)

          ! Check for suffix of time-dependent variables:
          idx = INDEX(info%name,'.TL')
          IF(idx == 0) THEN
            IF (TRIM(vname) /= TRIM(info%name)) CYCLE VAR_LOOP
          ELSE
            IF (TRIM(vname) /= TRIM(info%name(1:idx-1))) CYCLE VAR_LOOP
          ENDIF

          IF ((info%hgrid /= GRID_UNSTRUCTURED_CELL) .AND.  &
            & (info%hgrid /= GRID_UNSTRUCTURED_EDGE) .AND.  &
            & (info%hgrid /= GRID_UNSTRUCTURED_VERT)) &
            CYCLE VAR_LOOP

          ! Found it, add it to the variable list of optional
          ! diagnostics
         
          SELECT CASE(ll_varlevs(ivar))
          CASE (level_type_ml)
            p_opt_diag_list => p_nh_opt_diag(jg)%opt_diag_list         
          CASE (level_type_pl)
            p_opt_diag_list => p_nh_opt_diag(jg)%opt_diag_list_p
          CASE (level_type_hl)
            p_opt_diag_list => p_nh_opt_diag(jg)%opt_diag_list_z
          END SELECT

          nlev   = info%used_dimensions(2)
          ! set local values for "nblks" and "npromz"
          ptr_int_lonlat => lonlat_grid_list(ll_vargrid(ivar))%intp(jg)
          nblks_lonlat   =  (ptr_int_lonlat%nthis_local_pts - 1)/nproma + 1
          var_shape         =  info%used_dimensions(:)
          IF (info%vgrid == ZAXIS_SURFACE) THEN
            var_shape(2:3)   =  (/ 1, nblks_lonlat /)
          ELSE
            var_shape(3)     =  nblks_lonlat
          END IF

          SELECT CASE (info%hgrid)
          CASE (GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT)
            CALL add_var( p_opt_diag_list, info%name, p_opt_field_r3d,          &
              &           GRID_REGULAR_LONLAT, info%vgrid, info%cf, info%grib2, &
              &           ldims=var_shape, lrestart=.FALSE.,                    &
              &           loutput=.TRUE., new_element=new_element,              &
              &           hor_interp=create_hor_interp_metadata(                &
              &             hor_intp_type=HINTP_TYPE_NONE ) )
          CASE (GRID_UNSTRUCTURED_EDGE)
            IF(idx == 0) THEN
              vname = TRIM(ll_varlist(ivar))
            ELSE
              vname = TRIM(ll_varlist(ivar))//TRIM(info%name(idx:LEN(info%name)))
            END IF
            CALL add_var( p_opt_diag_list, TRIM(vname)//".X", p_opt_field_r3d,  &
              &           GRID_REGULAR_LONLAT, info%vgrid, info%cf, info%grib2, &
              &           ldims=var_shape, lrestart=.FALSE.,                    &
              &           loutput=.TRUE., new_element=new_element,              &
              &           hor_interp=create_hor_interp_metadata(                &
              &             hor_intp_type=HINTP_TYPE_NONE ))
            CALL add_var( p_opt_diag_list, TRIM(vname)//".Y", p_opt_field_r3d,  &
              &           GRID_REGULAR_LONLAT, info%vgrid, info%cf, info%grib2, &
              &           ldims=var_shape, lrestart=.FALSE.,                    &
              &           loutput=.TRUE., new_element=new_element_2,            &
              &           hor_interp=create_hor_interp_metadata(                &
              &             hor_intp_type=HINTP_TYPE_NONE ))
          CASE DEFAULT
            CALL finish(routine, "Unsupported grid type!")
          END SELECT
  
          ! link this new variable to the lon-lat grid:
          new_element%field%info%hor_interp%lonlat_id = ll_vargrid(ivar)
            
          !-- create and add post-processing task
          task => pp_task_insert(DEFAULT_PRIORITY2)
          WRITE (task%job_name, *) "horizontal interp. ",TRIM(info%name),", DOM ",jg
          task%data_input%p_nh_state      => NULL()
          task%data_input%prm_diag        => NULL()
          task%data_input%nh_pzlev_config => NULL()
          task%data_input%jg            =  jg           
          task%data_input%p_patch       => p_patch(jg)
          task%data_input%p_nh_opt_diag => p_nh_opt_diag(jg)
          task%data_input%p_int_state   => p_int_state(jg)
          task%job_type                 =  TASK_INTP_HOR_LONLAT
          task%activity                 =  new_simulation_status(l_output_step=.TRUE.)
          task%data_input%var           => element%field       ! set input variable
          task%data_output%var          => new_element%field   ! set output variable
          IF (info%hgrid == GRID_UNSTRUCTURED_EDGE) THEN
            new_element_2%field%info%hor_interp%lonlat_id = ll_vargrid(ivar)
            task%data_output%var_2      => new_element_2%field ! set Y-component
          END IF
          
          ! Flag. Denotes that at least one interpolation task has
          ! been created.
          l_horintp = .TRUE.
            
          found = .TRUE.
          
        ENDDO VAR_LOOP ! loop over vlist "i"
      ENDDO LIST_LOOP ! i = 1,nvar_lists
        
      ! Check that at least one element with this name has been found
      IF(.NOT. found) &
        CALL finish(routine,'No feasible variable found: ' &
        &   //TRIM(ll_varlist(ivar)))
    END DO ! ivar

    DEALLOCATE(ll_varlist, ll_vargrid, ll_varlevs, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    ! If at least one interpolation task has been created, we add a
    ! setup task which synchronizes the halo regions:
    IF (l_horintp) THEN
      task => pp_task_insert(HIGH_PRIORITY)
      WRITE (task%job_name, *) "horizontal interp. SYNC"
      task%job_type = TASK_INTP_SYNC
      task%activity = new_simulation_status(l_output_step=.TRUE.)
    END IF

    IF (dbg_level > 5)  CALL message(routine, "Done")
    
  END SUBROUTINE pp_scheduler_init


  !---------------------------------------------------------------
  !> Setup of pz-level interpolation tasks.
  !
  ! See SUBROUTINE pp_scheduler_init for further details.
  !
  SUBROUTINE pp_scheduler_init_pz(l_init_prm_diag)
    LOGICAL, INTENT(IN) :: l_init_prm_diag

    ! local variables
    CHARACTER(*), PARAMETER :: routine =  &
      &  TRIM("mo_pp_scheduler:pp_scheduler_init_pz")
    INTEGER                            :: &
      &  jg, ndom, ientr, nblks_c, ierrstat, ivar, i, idx, &
      &  iaxis, vgrid, nlev
    LOGICAL                            :: &
      &  l_jg_active, l_intp_p, l_intp_z, found
    TYPE (t_output_name_list), POINTER :: p_onl
    TYPE(t_job_queue),         POINTER :: task
    INTEGER                            :: shape3d(3)
    TYPE(t_cf_var)                     :: cf_desc
    TYPE(t_grib2_var)                  :: grib2_desc
    TYPE(t_nh_diag_pz),        POINTER :: p_diag_pz
    TYPE(t_var_list),          POINTER :: p_opt_diag_list_p, p_opt_diag_list_z, &
      &                                   p_opt_diag_list
    REAL(wp), POINTER                  :: p_opt_field_r3d(:,:,:)
    TYPE(t_list_element),      POINTER :: element, new_element
    ! variable lists (for all domains + output name lists):
    INTEGER                            :: nvars_pl, nvars_hl, nvars, nvars_predef
    CHARACTER(LEN=vname_len), TARGET, ALLOCATABLE  :: pl_varlist(:), hl_varlist(:)
    CHARACTER(LEN=vname_len), POINTER  :: varlist(:)
    CHARACTER(LEN=vname_len)           :: varlist_predef(5)
    CHARACTER(LEN=10)                  :: prefix
    TYPE(t_var_metadata), POINTER      :: info

    if (dbg_level > 5)  CALL message(routine, "Enter")
    ndom = SIZE(p_nh_opt_diag)
   
    ! loop over the domains and output definitions
    ! - check for which domain p- or z-level interpolation has been
    !   requested
    ! - register setup routine for pz-level interpolation (must be
    !   executed ahead of interpolation tasks)
    ! - for each variable: register post-processing task

    ALLOCATE(pl_varlist(ndom*max_var_pl), hl_varlist(ndom*max_var_hl), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    DOM_LOOP : DO jg=1,ndom
      IF (dbg_level > 8)  CALL message(routine, "DOM "//int2string(jg))

      !-- check if any output name list requests p- or z-level
      !-- interpolation for this domain, collect the list of variables

      ! loop in search of pressure-level interpolation      
      p_onl => first_output_name_list
      nvars_pl = 0
      l_intp_p = .FALSE.
      NML_LOOP_P : DO
        IF (.NOT.ASSOCIATED(p_onl)) EXIT NML_LOOP_P
      
        ! If dom(:) was not specified in namelist input, it is set
        ! completely to -1.  In this case all domains are wanted in
        ! the output
        l_jg_active = (ANY(p_onl%dom(:) == jg) .OR. p_onl%dom(1) <= 0)
        ! Selection criteria: 
        ! - domain is requested
        ! - p-level interpolation is requested
        IF (l_jg_active .AND. (p_onl%pl_varlist(1) /= ' ')) THEN
          ! check, if "all" variables are desired:
          IF (toupper(TRIM(p_onl%pl_varlist(1))) == 'ALL') THEN
            IF (dbg_level > 8)  CALL message(routine, "ALL vars for p-level")
            CALL get_all_var_names(pl_varlist, nvars_pl, &
              &                    opt_vert_intp_type=VINTP_TYPE_P_OR_Z, &
              &                    opt_loutput=.TRUE.)
            IF (nvars_pl > 0) l_intp_p = .TRUE.
            EXIT NML_LOOP_P
          END IF

          l_intp_p = .TRUE.
          DO ivar=1,max_var_pl
            IF (p_onl%pl_varlist(ivar) == ' ') CYCLE
            nvars_pl=nvars_pl+1
            pl_varlist(nvars_pl) = p_onl%pl_varlist(ivar)
          END DO
        END IF
        p_onl => p_onl%next
      END DO NML_LOOP_P

      ! loop in search of height-level interpolation
      p_onl => first_output_name_list
      nvars_hl = 0
      l_intp_z = .FALSE.
      NML_LOOP_H : DO
        IF (.NOT.ASSOCIATED(p_onl)) EXIT NML_LOOP_H
      
        ! If dom(:) was not specified in namelist input, it is set
        ! completely to -1.  In this case all domains are wanted in
        ! the output
        l_jg_active = (ANY(p_onl%dom(:) == jg) .OR. p_onl%dom(1) <= 0)
        ! Selection criteria: 
        ! - domain is requested
        ! - z-level interpolation is requested
        IF (l_jg_active .AND. (p_onl%hl_varlist(1) /= ' ')) THEN
          ! check, if "all" variables are desired:
          IF (toupper(TRIM(p_onl%hl_varlist(1))) == 'ALL') THEN
            IF (dbg_level > 8)  CALL message(routine, "ALL vars for z-level")
            CALL get_all_var_names(hl_varlist, nvars_hl, &
              &                    opt_vert_intp_type=VINTP_TYPE_Z,       &
              &                    opt_vert_intp_type2=VINTP_TYPE_P_OR_Z, &
              &                    opt_loutput=.TRUE.)
            IF (nvars_hl > 0) l_intp_z = .TRUE.
            EXIT NML_LOOP_H
          END IF

          l_intp_z = .TRUE.
          DO ivar=1,max_var_hl
            IF (p_onl%hl_varlist(ivar) == ' ') CYCLE
            nvars_hl=nvars_hl+1
            hl_varlist(nvars_hl) = p_onl%hl_varlist(ivar)
          END DO
        END IF
        p_onl => p_onl%next
      END DO NML_LOOP_H

      IF (dbg_level > 8)  THEN
        DO i=1,nvars_pl
          WRITE (message_text,*) "p var list: ", TRIM(pl_varlist(i))
          CALL message(routine, message_text)
        END DO
        DO i=1,nvars_hl
          WRITE (message_text,*) "h var list: ", TRIM(hl_varlist(i))
          CALL message(routine, message_text)
        END DO
      END IF
      ! now, we have total variables lists "hl_varlist(1:nvars_hl)"
      ! and "pl_varlist(1:nvars_pl)"
      
      ! skip domain if no p/z-interpolation requested:
      IF (.NOT. (l_intp_z .OR. l_intp_p)) CYCLE DOM_LOOP

      ! remove duplicates from variable lists
      CALL remove_duplicates(pl_varlist, nvars_pl)
      CALL remove_duplicates(hl_varlist, nvars_hl)

      !-- First, add some diagnostic variables which are essential for
      !-- p/z-level interpolation, e.g. temp_z, pres_z:
      p_diag_pz         => p_nh_opt_diag(jg)%diag_pz
      p_opt_diag_list_z => p_nh_opt_diag(jg)%opt_diag_list_z
      p_opt_diag_list_p => p_nh_opt_diag(jg)%opt_diag_list_p
      ientr     = 16   ! "entropy" of horizontal slice

      ! predefined array shapes
      nblks_c   = p_patch(jg)%nblks_c
      shape3d = (/ nproma, nh_pzlev_config(jg)%nzlev, nblks_c /)
      nvars_predef = 0
        
      ! temp         (nproma,nzlev,nblks_c)        
      cf_desc    = t_cf_var('temperature', 'K', 'temperature')
      grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
      nvars_predef = nvars_predef + 1
      varlist_predef(nvars_predef) = "temp"
      CALL add_var( p_opt_diag_list_z, varlist_predef(nvars_predef), p_diag_pz%z_temp, &
        & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE, cf_desc, grib2_desc, &
        & ldims=shape3d )
      
      ! pres         (nproma,nzlev,nblks_c)
      nvars_predef = nvars_predef + 1
      varlist_predef(nvars_predef) = "pres"
      cf_desc    = t_cf_var('pressure', 'Pa', 'pressure')
      grib2_desc = t_grib2_var(0, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_opt_diag_list_z, varlist_predef(nvars_predef), p_diag_pz%z_pres, &
        & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE, cf_desc, grib2_desc, &
        & ldims=shape3d )

      ! tracer_qv
      nvars_predef = nvars_predef + 1
      varlist_predef(nvars_predef) = "qv"
      cf_desc    = t_cf_var('tracer_qv', 'kg kg-1', 'specific_humidity')
      grib2_desc = t_grib2_var(0, 1, 201, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_opt_diag_list_z, varlist_predef(nvars_predef),          &
        & p_diag_pz%z_tracer_iqv, GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,     &
        & cf_desc, grib2_desc, ldims=shape3d)

      ! tot_qv
      nvars_predef = nvars_predef + 1
      varlist_predef(nvars_predef) = "tot_qv"
      cf_desc    = t_cf_var('tot_qv', '','total_specific_humidity')
      grib2_desc = t_grib2_var(0, 6, 6, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_opt_diag_list_z, varlist_predef(nvars_predef),          &
        & p_diag_pz%z_tot_cld_iqv, GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,    &
        & cf_desc, grib2_desc, ldims=shape3d)

      IF (l_intp_p) THEN
        shape3d = (/ nproma, nh_pzlev_config(jg)%nplev, nblks_c /)
        ! GEOPOT
        cf_desc    = t_cf_var('z', 'm2 s-2', 'geopotential')
        grib2_desc = t_grib2_var(0, 3, 4, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_opt_diag_list_p, 'z', p_diag_pz%p_geopot,             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, grib2_desc,      &
          & ldims=shape3d )
        ! temp
        cf_desc    = t_cf_var('temperature', 'K', 'temperature')
        grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_opt_diag_list_p, 'temp', p_diag_pz%p_temp,            &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, grib2_desc,      &
          & ldims=shape3d )
      END IF

      !-- register interpolation setup as a post-processing task

      task => pp_task_insert(HIGH_PRIORITY)
      task%data_input%p_int_state      => NULL()
      task%data_input%jg               =  jg           
      task%data_input%p_patch          => p_patch(jg)
      task%data_input%p_nh_state       => p_nh_state(jg)
      task%data_input%p_nh_opt_diag    => p_nh_opt_diag(jg)
      IF (l_init_prm_diag) THEN
        task%data_input%prm_diag       => prm_diag(jg)
      ELSE
        task%data_input%prm_diag       => NULL() 
      END IF
      task%data_input%nh_pzlev_config  => nh_pzlev_config(jg)
      task%activity     = new_simulation_status(l_output_step=.TRUE.)
      IF (l_intp_p) THEN
        task%job_name = "Init: pz-level interpolation, DOM "//TRIM(int2string(jg))
        task%job_type = TASK_INIT_VER_PZ
      ELSE
        task%job_name = "Init:  z-level interpolation, DOM "//TRIM(int2string(jg))
        task%job_type = TASK_INIT_VER_Z
      END IF

      !-- register clean-up routine as a post-processing task

      task => pp_task_insert(LOW_PRIORITY)
      task%activity     = new_simulation_status(l_output_step=.TRUE.)
      task%data_input%p_nh_opt_diag => p_nh_opt_diag(jg)
      task%job_name     = "Clean-up: pz-level interpolation, level "//TRIM(int2string(jg))
      task%job_type     = TASK_FINALIZE_PZ
      task%data_input%p_int_state      => NULL()
      task%data_input%jg               =  jg           
      task%data_input%p_patch          => p_patch(jg)
      task%data_input%p_nh_state       => p_nh_state(jg)
      task%data_input%prm_diag         => NULL() 
      task%data_input%nh_pzlev_config  => nh_pzlev_config(jg)

      ! remove already defined variables from list of requested output
      ! fields:
      CALL difference(hl_varlist, nvars_hl, varlist_predef, nvars_predef)
      IF (l_intp_p) &
        CALL difference(pl_varlist, nvars_pl, (/ "z   ", "temp" /), 2)

      !-- loop over requested p- and z-level variables, add variables
      !-- ("add_var") and register interpolation tasks:
      DO iaxis=1,2
        IF (iaxis == 1) THEN
          prefix  =  "z-level"
          varlist => hl_varlist
          nvars   =  nvars_hl
          nlev    =  nh_pzlev_config(jg)%nzlev
          vgrid   =  ZAXIS_ALTITUDE
          p_opt_diag_list => p_opt_diag_list_z
        END IF
        IF (iaxis == 2) THEN
          prefix  =  "p-level"
          varlist => pl_varlist
          nvars   =  nvars_pl
          nlev    =  nh_pzlev_config(jg)%nplev
          vgrid   =  ZAXIS_PRESSURE
          p_opt_diag_list => p_opt_diag_list_p
        END IF
        
        DO ivar=1,nvars
          IF (dbg_level > 8) &
            CALL message(routine, &
            &      TRIM(prefix)//": Looking for input var '"//TRIM(varlist(ivar))//"'")
          found = .FALSE.
        
          !- loop over model level variables
          ! Note that there may be several variables with different time levels,
          ! we just add unconditionally all
          DO i = 1,nvar_lists
            ! Do not inspect lists which are disabled for output
            IF (.NOT. var_lists(i)%p%loutput) CYCLE
            ! loop only over model level variables
            IF (var_lists(i)%p%vlevel_type /= level_type_ml) CYCLE         
            ! loop only over variables of current domain
            IF (var_lists(i)%p%patch_id /= jg) CYCLE

            element => NULL()
            DO
              IF(.NOT.ASSOCIATED(element)) THEN
                element => var_lists(i)%p%first_list_element
              ELSE
                element => element%next_list_element
              ENDIF
              IF(.NOT.ASSOCIATED(element)) EXIT

              info => element%field%info
              ! Do not inspect element if it is a container
              IF (info%lcontainer) CYCLE
              ! Do not inspect element if "loutput=.false."
              IF (.NOT. info%loutput) CYCLE
              ! Inspect element only if vertical interpolation matches
              IF (iaxis == 1) THEN
                IF ((info%vert_interp%vert_intp_type/=VINTP_TYPE_Z) .AND.  &
                  & (info%vert_interp%vert_intp_type/=VINTP_TYPE_P_OR_Z)) CYCLE
              END IF
              IF (iaxis == 2) THEN
                IF (info%vert_interp%vert_intp_type/=VINTP_TYPE_P_OR_Z) CYCLE
              END IF
              IF (info%vert_interp%vert_intp_type==VINTP_TYPE_NONE) CYCLE

              ! Check for matching name
              idx = INDEX(info%name,'.TL')
              IF(idx == 0) THEN
                IF (varlist(ivar) /= info%name) CYCLE
              ELSE
                IF (varlist(ivar) /= info%name(1:idx-1)) CYCLE
              ENDIF

              ! Found it, add it to the variable list of optional
              ! diagnostics
              IF ( (info%used_dimensions(1) /= nproma) .OR.  &
                &  (info%used_dimensions(3) /= nblks_c) ) THEN
                CALL finish(routine, "Unexpected field size!")
              END IF
              shape3d  = (/ nproma, nlev, nblks_c /)

              CALL add_var( p_opt_diag_list, info%name, p_opt_field_r3d, &
                &           info%hgrid, vgrid, info%cf, info%grib2,      &
                &           ldims=shape3d, lrestart=.FALSE.,           &
                &           loutput=.TRUE., new_element=new_element)

              !-- add post-processing task

              task => pp_task_insert(DEFAULT_PRIORITY1)
              task%job_name        =  &
                &  TRIM(prefix)//" interp. "//TRIM(info%name)  &
                &  //", DOM "//TRIM(int2string(jg))
              task%job_type        =  TASK_INTP_VER_PZLEV
              task%activity        =  new_simulation_status(l_output_step=.TRUE.)
              task%data_input%jg               =  jg 
              task%data_input%p_patch          => p_patch(jg)          
              task%data_input%p_int_state      => NULL()
              task%data_input%p_nh_state       => p_nh_state(jg)
              task%data_input%nh_pzlev_config  => nh_pzlev_config(jg)
              task%data_input%p_nh_opt_diag    => p_nh_opt_diag(jg)
              IF (l_init_prm_diag) THEN
                task%data_input%prm_diag       => prm_diag(jg)
              ELSE
                task%data_input%prm_diag       => NULL() 
              END IF
              task%data_input%var  => element%field       ! set input variable
              task%data_output%var => new_element%field   ! set output variable

              found = .TRUE.

            ENDDO ! loop over vlist "i"
          ENDDO ! i = 1,nvar_lists
        
          ! Check that at least one element with this name has been found
          IF(.NOT. found) &
            CALL finish(routine,TRIM(prefix)//" interpolation: No feasible variable found: "&
            &                   //TRIM(varlist(ivar)))
        END DO ! ivar
      END DO ! height/pressure axis
    END DO DOM_LOOP  ! jg
    
    DEALLOCATE(pl_varlist, hl_varlist, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    IF (dbg_level > 5)  CALL message(routine, "Done")
    
  END SUBROUTINE pp_scheduler_init_pz


  !---------------------------------------------------------------
  !> Loop over job queue, call active tasks.
  !
  SUBROUTINE pp_scheduler_process(simulation_status)
    TYPE(t_simulation_status), INTENT(IN) :: simulation_status
    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_pp_scheduler:pp_scheduler_process")
    TYPE(t_job_queue), POINTER :: ptr_task

    ptr_task => job_queue
    ! loop over job queue
    LOOP_JOB : DO
      IF (.NOT. ASSOCIATED(ptr_task)) EXIT
      IF (.NOT. pp_task_is_active(ptr_task, simulation_status)) THEN
        ptr_task => ptr_task%next
        CYCLE LOOP_JOB
      END IF

      IF (dbg_level > 5) THEN
        WRITE(message_text,*) "Staging task '", TRIM(ptr_task%job_name), "'"
        CALL message(routine, TRIM(message_text))
      END IF

      SELECT CASE ( ptr_task%job_type )
      CASE ( TASK_INIT_VER_PZ, TASK_INIT_VER_Z, TASK_FINALIZE_PZ)
        CALL pp_task_pzlev_setup(ptr_task)
      CASE ( TASK_INTP_HOR_LONLAT )
        CALL pp_task_lonlat(ptr_task)
      CASE ( TASK_INTP_VER_PZLEV )
        CALL pp_task_pzlev(ptr_task)
      CASE ( TASK_INTP_SYNC )
        call pp_task_sync()
      CASE DEFAULT
        CALL finish(routine, "Unknown post-processing job.")
      END SELECT

      ptr_task => ptr_task%next
    END DO LOOP_JOB
    
  END SUBROUTINE pp_scheduler_process


  !---------------------------------------------------------------
  !> Destruction of post-processing job queue
  !
  SUBROUTINE pp_scheduler_finalize()
    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_pp_scheduler:pp_scheduler_finalize")
    INTEGER                    :: ierrstat
    TYPE(t_job_queue), POINTER :: tmp
    
    CALL message(routine, "")
    ! destroy linked list
    DO
      IF (.NOT. ASSOCIATED(job_queue)) EXIT
      ! remove list head
      tmp => job_queue%next
      DEALLOCATE(job_queue, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
      job_queue => tmp
    END DO
    
  END SUBROUTINE pp_scheduler_finalize


  !--- POST-PROCESSING TASKS ---------------------------------------------------------

  !---------------------------------------------------------------
  !> Performs interpolation of a variable onto a regular grid.
  !
  !  This is only a wrapper for the corresponding routines from the
  !  interpolation module.
  SUBROUTINE pp_task_lonlat(ptr_task)
    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_pp_scheduler:pp_task_lonlat")
    INTEGER                            :: &
      &  nblks_ll, npromz_ll, lonlat_id, jg,     &
      &  in_var_idx, out_var_idx, out_var_idx_2, &
      &  ierrstat, dim1, dim2
    TYPE (t_var_list_element), POINTER :: in_var, out_var, out_var_2
    TYPE (t_var_metadata),     POINTER :: p_info
    TYPE (t_lon_lat_intp),     POINTER :: ptr_int_lonlat
    REAL(wp), ALLOCATABLE              :: tmp_var(:,:,:)
    TYPE(t_patch),             POINTER :: p_patch

    p_patch        => ptr_task%data_input%p_patch      ! patch
    p_info         => ptr_task%data_input%var%info
    in_var         => ptr_task%data_input%var
    out_var        => ptr_task%data_output%var

    lonlat_id      =  ptr_task%data_output%var%info%hor_interp%lonlat_id
    jg             =  ptr_task%data_input%jg
    ptr_int_lonlat => lonlat_grid_list(lonlat_id)%intp(jg)

    in_var_idx        = 1
    IF (in_var%info%lcontained)  in_var_idx  = in_var%info%ncontained
    out_var_idx       = 1
    IF (out_var%info%lcontained) out_var_idx = out_var%info%ncontained

    ! For edge-based interpolation: retrieve data on Y-component:
    IF (ASSOCIATED(ptr_task%data_output%var_2)) THEN
      out_var_2      => ptr_task%data_output%var_2
      out_var_idx_2  =  1
      IF (out_var_2%info%lcontained) out_var_idx_2 = out_var_2%info%ncontained
    END IF

    nblks_ll  = (ptr_int_lonlat%nthis_local_pts - 1)/nproma + 1
    npromz_ll =  ptr_int_lonlat%nthis_local_pts - (nblks_ll-1)*nproma

    SELECT CASE (p_info%hgrid)
    CASE (GRID_UNSTRUCTURED_CELL)
      IF (p_info%vgrid == ZAXIS_SURFACE) THEN
        ! For 2D variables (nproma, nblks) we first copy this to 1-level
        ! 3D variable (nproma, nlevs, nblks). This requires a temporary
        ! variable:
        dim1 = p_info%used_dimensions(1)
        dim2 = p_info%used_dimensions(2)
        ALLOCATE(tmp_var(dim1, 1, dim2), STAT=ierrstat)
        tmp_var(:,1,:) = in_var%r_ptr(:,:,1,1,1)

        ! for cell-based variables: interpolate gradients (finite
        ! differences) and reconstruct
        CALL rbf_interpol_lonlat_nl(              &
          &   tmp_var(:,:,:),                     &
          &   ptr_int_lonlat,                     &
          &   out_var%r_ptr(:,:,:,out_var_idx,1), &
          &   nblks_ll, npromz_ll)

        ! clean up:
        DEALLOCATE(tmp_var, STAT=ierrstat)
        IF (ierrstat /= SUCCESS)  CALL finish (routine, 'deallocation failed')
      ELSE
        ! for cell-based variables: interpolate gradients (finite
        ! differences) and reconstruct
        CALL rbf_interpol_lonlat_nl(              &
          &   in_var%r_ptr(:,:,:,in_var_idx,1),   &
          &   ptr_int_lonlat,                     &
          &   out_var%r_ptr(:,:,:,out_var_idx,1), &
          &   nblks_ll, npromz_ll)
      END IF ! 2D
      ! --------------------------------------------------------------
      !
    CASE (GRID_UNSTRUCTURED_EDGE)
      IF (p_info%vgrid == ZAXIS_SURFACE) THEN
        ! For 2D variables (nproma, nblks) we first copy this to 1-level
        ! 3D variable (nproma, nlevs, nblks). This requires a temporary
        ! variable:
        dim1 = p_info%used_dimensions(1)
        dim2 = p_info%used_dimensions(2)
        ALLOCATE(tmp_var(dim1, 1, dim2), STAT=ierrstat)
        tmp_var(:,1,:) = in_var%r_ptr(:,:,1,1,1)

        ! for edge-based variables: simple interpolation
        CALL rbf_vec_interpol_lonlat_nl(              &
          &   tmp_var(:,:,:),                         &
          &   ptr_int_lonlat,                         &
          &   out_var%r_ptr(:,:,:,out_var_idx,1),     &
          &   out_var_2%r_ptr(:,:,:,out_var_idx_2,1), &
          &   nblks_ll, npromz_ll)
        ! clean up:
        DEALLOCATE(tmp_var, STAT=ierrstat)
        IF (ierrstat /= SUCCESS)  CALL finish (routine, 'deallocation failed')
      ELSE
        ! for edge-based variables: simple interpolation
        CALL rbf_vec_interpol_lonlat_nl(              &
          &   in_var%r_ptr(:,:,:,in_var_idx,1),       &
          &   ptr_int_lonlat,                         &
          &   out_var%r_ptr(:,:,:,out_var_idx,1),     &
          &   out_var_2%r_ptr(:,:,:,out_var_idx_2,1), &
          &   nblks_ll, npromz_ll)
      END IF ! 2D
      ! --------------------------------------------------------------
      !
    CASE (GRID_UNSTRUCTURED_VERT)
      ! vertex-based variables (e.g. vorticity "VOR") are treated in a
      ! special way: They are first interpolated onto the cell centers
      ! and afterwards treated as scalar variables of type
      ! "GRID_UNSTRUCTURED_CELL"

      ! Note: The whole process is overly expensive! For the case that
      ! there are other vertex-based variables for output, this
      ! routine must be optimized!

      ! this requires a temporary variable:
      dim1 = p_info%used_dimensions(1)
      dim2 = p_info%used_dimensions(2)
      ALLOCATE(tmp_var(dim1, dim2, p_patch%nblks_c), STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation failed')
      CALL verts2cells_scalar( in_var%r_ptr(:,:,:,in_var_idx,1), p_patch,      &
        &                      ptr_task%data_input%p_int_state%verts_aw_cells, &
        &                      tmp_var(:,:,:), 1, dim2 )
      CALL rbf_interpol_lonlat_nl( tmp_var(:,:,:), ptr_int_lonlat,     &
        &                          out_var%r_ptr(:,:,:,out_var_idx,1), &
        &                          nblks_ll, npromz_ll )
      ! clean up:
      DEALLOCATE(tmp_var, STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'deallocation failed')
    CASE DEFAULT
      CALL finish(routine, 'Unknown grid type.')
    END SELECT

  END SUBROUTINE pp_task_lonlat


  !---------------------------------------------------------------
  !> Performs synchronization of halo regions.
  !  This is necessary before starting the lon-lat interpolation.
  !  All variables that are part of an interpolation task are
  !  synchronized.
  !
  ! To avoid unnecessary overhead in the synchronization,
  ! several improvements are possible:
  !   - Use "*_mult" synchronization routines

  !   - Copy 2D fields into a 3D field which is
  !     synchronized. Afterwards, 2D fields are extracted again from
  !     the temporary 3D field.
  !
  !   - Introduce some kind of meta information of each variable list.
  !     For example, prognostic fields must not be synchronized,
  !     therefore the corresponding variable lists can be marked as
  !     "skip_sync".
  !
  SUBROUTINE pp_task_sync()
    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_pp_scheduler:pp_task_sync")
    TYPE(t_job_queue), POINTER :: ptr_task
    INTEGER                            :: in_var_idx
    TYPE (t_var_list_element), POINTER :: in_var
    TYPE (t_var_metadata),     POINTER :: p_info
    TYPE(t_patch),             POINTER :: p_patch

    ptr_task => job_queue
    ! loop over job queue
    LOOP_JOB : DO
      IF (.NOT. ASSOCIATED(ptr_task)) EXIT
      IF (ptr_task%job_type == TASK_INTP_HOR_LONLAT) THEN
        p_info      => ptr_task%data_input%var%info
        p_patch     => ptr_task%data_input%p_patch
        in_var      => ptr_task%data_input%var
        in_var_idx  =  1
        IF (in_var%info%lcontained) &
          &  in_var_idx = in_var%info%ncontained
        
        SELECT CASE (p_info%hgrid)
        CASE (GRID_UNSTRUCTURED_CELL)
          IF (p_info%vgrid == ZAXIS_SURFACE) THEN
            CALL sync_patch_array(SYNC_C, p_patch, in_var%r_ptr(:,:,1,1,1) )
          ELSE
            CALL sync_patch_array(SYNC_C, p_patch, in_var%r_ptr(:,:,:,in_var_idx,1) )
          END IF
        CASE (GRID_UNSTRUCTURED_EDGE)
          IF (p_info%vgrid == ZAXIS_SURFACE) THEN
            CALL sync_patch_array(SYNC_E, p_patch, in_var%r_ptr(:,:,1,1,1) )
          ELSE
            CALL sync_patch_array(SYNC_E, p_patch, in_var%r_ptr(:,:,:,in_var_idx,1) )
          END IF
        CASE (GRID_UNSTRUCTURED_VERT)
          CALL sync_patch_array(SYNC_V, p_patch, in_var%r_ptr(:,:,:,in_var_idx,1) )
        CASE DEFAULT
          CALL finish(routine, 'Unknown grid type.')
        END SELECT
      END IF
      !
      ptr_task => ptr_task%next
    END DO LOOP_JOB

  END SUBROUTINE pp_task_sync


  !---------------------------------------------------------------
  !> Performs setup of vertical interpolation.
  !
  !  This is only a wrapper for the corresponding routines from the
  !  interpolation module.
  SUBROUTINE pp_task_pzlev_setup(ptr_task)
    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_pp_scheduler:pp_task_pzlev_setup")
    INTEGER                            :: &
      &  jg, nzlev, nplev
    TYPE(t_patch),             POINTER :: p_patch
    TYPE(t_nh_metrics),        POINTER :: p_metrics    

    ! prognostic state: note that we only use p_prog(nnow(jg))
    TYPE(t_nh_prog),           POINTER :: p_prog
    TYPE(t_nh_diag),           POINTER :: p_diag
    TYPE(t_nh_diag_pz),        POINTER :: p_diag_pz
    TYPE(t_nwp_phy_diag),      POINTER :: prm_diag

    TYPE(t_nh_pzlev_config),   POINTER :: nh_pzlev_config

    ! patch, state, and metrics
    jg             =  ptr_task%data_input%jg
    p_patch        => ptr_task%data_input%p_patch
    p_metrics      => ptr_task%data_input%p_nh_state%metrics
    p_prog         => ptr_task%data_input%p_nh_state%prog(nnow(jg))
    p_diag         => ptr_task%data_input%p_nh_state%diag
    p_diag_pz      => ptr_task%data_input%p_nh_opt_diag%diag_pz
    prm_diag       => ptr_task%data_input%prm_diag

    ! pz-level interpolation data
    nh_pzlev_config   => ptr_task%data_input%nh_pzlev_config

    nzlev             =  nh_pzlev_config%nzlev
    nplev             =  nh_pzlev_config%nplev
                      
    SELECT CASE ( ptr_task%job_type )
    CASE ( TASK_INIT_VER_PZ )
      ! build data structure "vcoeff" containing coefficient tables
      CALL prepare_vert_interp(p_patch, p_prog, p_diag, prm_diag, nzlev, nplev,  & ! in
        &                      p_diag_pz%z_temp, p_diag_pz%z_tracer_iqv,         & ! inout
        &                      p_diag_pz%z_tot_cld_iqv,                          & ! inout
        &                      p_diag_pz%z_pres, p_diag_pz%p_geopot,             & ! inout
        &                      p_diag_pz%p_temp,                                 & ! inout
        &                      nh_pzlev_config, p_metrics,                       & ! in
        &                      p_diag_pz%vcoeff_z, p_diag_pz%vcoeff_p )            ! inout
      !
    CASE ( TASK_INIT_VER_Z )
      ! build data structure "vcoeff" containing coefficient tables
      CALL prepare_vert_interp(p_patch, p_prog, p_diag, prm_diag, nzlev, nplev,  & ! in
        &                      p_diag_pz%z_temp, p_diag_pz%z_tracer_iqv,         & ! inout
        &                      p_diag_pz%z_tot_cld_iqv,                          & ! inout
        &                      p_diag_pz%z_pres, p_diag_pz%p_geopot,             & ! inout
        &                      p_diag_pz%p_temp,                                 & ! inout
        &                      nh_pzlev_config, p_metrics,                       & ! in
        &                      p_diag_pz%vcoeff_z )                                ! inout
      !
    CASE ( TASK_FINALIZE_PZ )
      ! deallocate coefficient tables:
      CALL vcoeff_deallocate(p_diag_pz%vcoeff_z)
      CALL vcoeff_deallocate(p_diag_pz%vcoeff_p)
      !
    CASE DEFAULT
      CALL finish(routine, "Internal error.")
    END SELECT ! vert_intp_method

  END SUBROUTINE pp_task_pzlev_setup


  !---------------------------------------------------------------
  !> Performs vertical interpolation of a variable onto p/z-levels.
  !
  !  This is only a wrapper for the corresponding routines from the
  !  interpolation module.
  SUBROUTINE pp_task_pzlev(ptr_task)
    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_pp_scheduler:pp_task_pzlev")
    INTEGER                            :: &
      &  vert_intp_type, vert_intp_method, jg,  &
      &  in_var_idx, out_var_idx, nlev, nlevp1, &
      &  nzlev, nplev, npzlev
    TYPE(t_patch),             POINTER :: p_patch
    TYPE(t_nh_metrics),        POINTER :: p_metrics    
    TYPE(t_nh_prog),           POINTER :: p_prog
    TYPE(t_nh_diag),           POINTER :: p_diag
    TYPE(t_nh_diag_pz),        POINTER :: p_diag_pz
    TYPE(t_nwp_phy_diag),      POINTER :: prm_diag
    TYPE(t_vert_interp_meta),  POINTER :: pzlev_flags

    TYPE (t_var_list_element), POINTER :: in_var, out_var
    TYPE(t_var_metadata),      POINTER :: p_info
    TYPE(t_vcoeff),            POINTER :: vcoeff
    TYPE(t_nh_pzlev_config),   POINTER :: nh_pzlev_config
    REAL(wp),                  POINTER :: p_z3d(:,:,:)

    LOGICAL                            :: &
      &  l_hires_intp, l_restore_fricred, l_loglin, &
      &  l_extrapol, l_satlimit, l_restore_pbldev,  &
      &  l_pd_limit, l_restore_sfcinv, l_hires_corr
    REAL(wp)                           :: &
      &  lower_limit, extrapol_dist

    ! input/output field for this task
    p_info            => ptr_task%data_input%var%info
    in_var            => ptr_task%data_input%var
    out_var           => ptr_task%data_output%var

    in_var_idx        = 1
    if (ptr_task%data_input%var%info%lcontained) &
      in_var_idx = ptr_task%data_input%var%info%ncontained
    out_var_idx       = 1
    if (ptr_task%data_output%var%info%lcontained) &
      out_var_idx = ptr_task%data_output%var%info%ncontained

    !--- load some items from input/output data structures
    vert_intp_type   = p_info%vert_interp%vert_intp_type
    vert_intp_method = p_info%vert_interp%vert_intp_method

    ! patch, state, and metrics
    jg                =  ptr_task%data_input%jg
    p_patch           => ptr_task%data_input%p_patch
    p_metrics         => ptr_task%data_input%p_nh_state%metrics
    p_prog            => ptr_task%data_input%p_nh_state%prog(nnow(jg))
    p_diag            => ptr_task%data_input%p_nh_state%diag
    p_diag_pz         => ptr_task%data_input%p_nh_opt_diag%diag_pz
    prm_diag          => ptr_task%data_input%prm_diag

    nh_pzlev_config   => ptr_task%data_input%nh_pzlev_config
    nlev              = p_patch%nlev
    nlevp1            = p_patch%nlevp1
    nzlev             = nh_pzlev_config%nzlev
    nplev             = nh_pzlev_config%nplev

    ! pz-level interpolation data
    IF ( vert_intp_type == VINTP_TYPE_Z) THEN
      ! vertical levels for z-level interpolation
      npzlev  =   nzlev
      vcoeff  =>  p_diag_pz%vcoeff_z
      p_z3d   =>  nh_pzlev_config%z3d
    ELSE
      ! vertical levels for p-level interpolation
      npzlev  =   nplev
      vcoeff  =>  p_diag_pz%vcoeff_p
      p_z3d   =>  p_diag_pz%p_geopot
    END IF
                     
    ! interpolation flags + parameters
    pzlev_flags => in_var%info%vert_interp
    l_hires_intp      = pzlev_flags%l_hires_intp      
    l_restore_fricred = pzlev_flags%l_restore_fricred 
    l_loglin          = pzlev_flags%l_loglin          
    l_extrapol        = pzlev_flags%l_extrapol        
    l_satlimit        = pzlev_flags%l_satlimit        
    l_restore_pbldev  = pzlev_flags%l_restore_pbldev  
    l_pd_limit        = pzlev_flags%l_pd_limit
    l_restore_sfcinv  = pzlev_flags%l_restore_sfcinv 
    l_hires_corr      = pzlev_flags%l_hires_corr     
    lower_limit       = pzlev_flags%lower_limit       
    extrapol_dist     = pzlev_flags%extrapol_dist    

    !-- perform some consistency checks
    IF (p_info%ndims /= 3) &
      & CALL finish(routine, "Wrong number of variables dimensions!")
    IF (.NOT. vcoeff%l_initialized) &
      CALL finish(routine, "Interpolation coefficients not yet initialized!")

    !--- actually perform vertical interpolation task
    SELECT CASE ( vert_intp_method )
    CASE ( VINTP_METHOD_UV )
      IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_UV")
      CALL uv_intp(in_var%r_ptr(:,:,:,in_var_idx,1),                     & !in
        &          out_var%r_ptr(:,:,:,out_var_idx,1),                   & !out
        &          p_metrics%z_mc, p_z3d,                                & !in
        &          p_patch%nblks_c, p_patch%npromz_c, nlev, npzlev,      & !in
        &          vcoeff%coef1, vcoeff%coef2,                           & !in
        &          vcoeff%coef3, vcoeff%wfac_lin,                        & !in
        &          vcoeff%idx0_cub, vcoeff%idx0_lin,                     & !in
        &          vcoeff%bot_idx_cub, vcoeff%bot_idx_lin,               & !in
        &          vcoeff%wfacpbl1, vcoeff%kpbl1, vcoeff%wfacpbl2,       & !in
        &          vcoeff%kpbl2,                                         & !in
        &          l_hires_intp=l_hires_intp,                            & !in
        &          l_restore_fricred=l_restore_fricred )                   !in
      !
    CASE ( VINTP_METHOD_LIN )        
      IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_LIN")
      CALL lin_intp(in_var%r_ptr(:,:,:,in_var_idx,1),                    & !inout
        &           out_var%r_ptr(:,:,:,out_var_idx,1),                  & !out
        &           p_patch%nblks_c, p_patch%npromz_c, nlev, npzlev,     & !in
        &           vcoeff%wfac_lin, vcoeff%idx0_lin,                    & !in
        &           vcoeff%bot_idx_lin, vcoeff%wfacpbl1,                 & !in
        &           vcoeff%kpbl1, vcoeff%wfacpbl2, vcoeff%kpbl2,         & !in
        &           l_loglin=l_loglin,                                   & !in
        &           l_extrapol=l_extrapol, l_pd_limit=l_pd_limit,        & !in
        &           lower_limit=lower_limit )                              !in
      !
    CASE ( VINTP_METHOD_LIN_NLEVP1 )        
      IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_LIN_NLEVP1")
      CALL lin_intp(in_var%r_ptr(:,:,:,in_var_idx,1),                    & !inout
        &           out_var%r_ptr(:,:,:,out_var_idx,1),                  & !out
        &           p_patch%nblks_c, p_patch%npromz_c, nlevp1, npzlev,   & !in
        &           vcoeff%wfac_lin_nlevp1, vcoeff%idx0_lin_nlevp1,      & !in
        &           vcoeff%bot_idx_lin_nlevp1,                           & !in
        &           vcoeff%wfacpbl1_nlevp1, vcoeff%kpbl1_nlevp1,         & !in
        &           vcoeff%wfacpbl2_nlevp1, vcoeff%kpbl2_nlevp1,         & !in
        &           l_loglin=l_loglin,                                   & !in
        &           l_extrapol=l_extrapol, l_pd_limit=l_pd_limit,        & !in
        &           lower_limit=lower_limit )                              !in
      !
    CASE (VINTP_METHOD_QV )
      IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_QV")
      CALL qv_intp(in_var%r_ptr(:,:,:,in_var_idx,1),                   & !in
        &          out_var%r_ptr(:,:,:,out_var_idx,1),                 & !out
        &          p_metrics%z_mc, p_z3d, p_diag%temp,                 & !in
        &          p_diag%pres, p_diag_pz%p_temp, nh_pzlev_config%p3d, & !in
        &          p_patch%nblks_c, p_patch%npromz_c, nlev, npzlev,    & !in
        &          vcoeff%coef1, vcoeff%coef2, vcoeff%coef3,           & !in
        &          vcoeff%wfac_lin, vcoeff%idx0_cub, vcoeff%idx0_lin,  & !in
        &          vcoeff%bot_idx_cub, vcoeff%bot_idx_lin,             & !in
        &          vcoeff%wfacpbl1, vcoeff%kpbl1,                      & !in
        &          vcoeff%wfacpbl2, vcoeff%kpbl2,                      & !in
        &          l_satlimit=l_satlimit, lower_limit=lower_limit,     & !in
        &          l_restore_pbldev=l_restore_pbldev )                   !in
    END SELECT ! vert_intp_method

  END SUBROUTINE pp_task_pzlev


  !--- UTILITY ROUTINES --------------------------------------------------------------

  !---------------------------------------------------------------
  !> @return .TRUE. if given post-processing task is in active state.
  !
  ! Tasks may be inactive, e.g. outside the output intervals.
  FUNCTION pp_task_is_active(ptr_task, simulation_status)
    LOGICAL :: pp_task_is_active
    TYPE(t_job_queue), POINTER :: ptr_task
    TYPE(t_simulation_status),  INTENT(IN) :: simulation_status

    ! compare simulation status to post-processing tasks activity
    ! flags, then check if any of the activity conditions is
    ! fulfilled:
    pp_task_is_active = .FALSE.
    IF (ANY(ptr_task%activity%status_flags(:)  .AND.  &
      &     simulation_status%status_flags(:))) pp_task_is_active = .TRUE.

  END FUNCTION pp_task_is_active  


  !---------------------------------------------------------------
  !> Insert task into job queue.
  FUNCTION pp_task_insert(job_priority) RESULT(element)
    INTEGER, INTENT(IN) :: job_priority
    TYPE(t_job_queue), POINTER :: element
    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_pp_scheduler:pp_task_insert")
    TYPE(t_job_queue), POINTER :: tmp, nb_left
    INTEGER                    :: ierrstat

    IF (dbg_level > 5) &
         & CALL message(routine, "Inserting pp task")
    
    ! find the correct position in list:
    tmp     => job_queue
    nb_left => NULL()
    DO
      IF (.NOT. ASSOCIATED(tmp))                EXIT
      IF (tmp%job_priority > job_priority) EXIT
      nb_left => tmp
      tmp     => tmp%next
    END DO
    ! insert element into linked list
    IF (ASSOCIATED(nb_left)) THEN
      ALLOCATE(nb_left%next, stat=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      element => nb_left%next
    ELSE
      ALLOCATE(job_queue, stat=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      element => job_queue
    END IF
    ! fill new list element with data
    element%job_priority =  job_priority
    element%next         => tmp
  END FUNCTION pp_task_insert


  !------------------------------------------------------------------------------------------------


  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for "t_simulation_status" variables
  ! 
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION new_simulation_status(l_output_step, l_first_step, l_last_step)  &
    RESULT(sim_status)

    TYPE(t_simulation_status) :: sim_status
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  l_output_step, l_first_step, l_last_step

    ! set default values
    sim_status%status_flags(:) = (/ .FALSE., .FALSE., .FALSE. /)
    ! supersede with user definitions
    CALL assign_if_present(sim_status%status_flags(1), l_output_step)
    CALL assign_if_present(sim_status%status_flags(2), l_first_step)
    CALL assign_if_present(sim_status%status_flags(3), l_last_step)

  END FUNCTION new_simulation_status


  !------------------------------------------------------------------------------------------------
  !
  ! private routine to assign values if actual parameters are present
  !
  SUBROUTINE assign_if_present_logical (y,x)
    LOGICAL, INTENT(inout)        :: y
    LOGICAL, INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_logical

END MODULE mo_pp_scheduler
