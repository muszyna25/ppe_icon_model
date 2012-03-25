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
!!
!! TODO[FP] So far, pz-level + lon-lat interpolation is not yet implemented!
!!
!! TODO[FP] To increase performance, one should allocate/deallocate
!!          the vertical interpolation coefficient tables "vcoeff"
!!          only when necessary.
!!
!! -----------------------------------------------------------------------------------
MODULE mo_pp_scheduler

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message, message_text, finish
  USE mo_impl_constants,          ONLY: SUCCESS,                      &
    & VINTP_TYPE_Z, VINTP_TYPE_P, VINTP_TYPE_P_OR_Z, VINTP_TYPE_NONE, &
    & VINTP_METHOD_UV, VINTP_METHOD_LIN,                              &     
    & VINTP_METHOD_QV, max_dom    
  USE mo_model_domain,            ONLY: t_patch
  USE mo_var_list,                ONLY: t_var_list, new_var_list, &
  &                                     default_var_list_settings, add_var, &
  &                                     get_all_var_names
  USE mo_var_list_element,        ONLY: level_type_ml
  USE mo_var_list_element,        ONLY: t_var_list_element
  USE mo_var_metadata,            ONLY: t_var_metadata, t_vert_interp_meta
  USE mo_intp,                    ONLY: verts2cells_scalar
  USE mo_intp_data_strc,          ONLY: t_int_state
  USE mo_nh_vert_interp,          ONLY: prepare_vert_interp, &
    &                                   lin_intp, uv_intp, qv_intp 
  USE mo_nonhydro_types,          ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_opt_diagnostics,         ONLY: t_nh_diag_pz, t_nh_opt_diag, t_vcoeff, &
    &                                   vcoeff_deallocate
  USE mo_nwp_phy_state,           ONLY: t_nwp_phy_diag
  USE mo_run_config,              ONLY: iqv
  USE mo_nh_pzlev_config,         ONLY: t_nh_pzlev_config
  USE mo_name_list_output_config, ONLY: t_output_name_list, max_var_pl, max_var_hl, &
    &                                   vname_len
  USE mo_parallel_config,         ONLY: nproma
  USE mo_dynamics_config,         ONLY: nnow
  USE mo_util_string,             ONLY: int2string, remove_duplicates, &
    &                                   difference, toupper
  USE mo_cf_convention,           ONLY: t_cf_var
  USE mo_grib2,                   ONLY: t_grib2_var
  USE mo_cdi_constants,           ONLY: GRID_CELL, GRID_REFERENCE,               &
    &                                   GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,  &
    &                                   ZAXIS_PRESSURE
  USE mo_linked_list,             ONLY: t_list_element

  IMPLICIT NONE

  ! max. name string length
  INTEGER, PARAMETER :: MAX_NAME_LENGTH   =   64
  INTEGER, PARAMETER :: HIGH_PRIORITY     =    0  
  INTEGER, PARAMETER :: LOW_PRIORITY      =  100  
  INTEGER, PARAMETER :: DEFAULT_PRIORITY  =   10  

  ! level of output verbosity
  INTEGER :: dbg_level = 0
  
  !--- Available post-processing tasks
  !------ setup tasks (coefficients,...)
  INTEGER, PARAMETER :: TASK_INIT_VER_PZ       = 1  !< task: setup pz-interpolation
  INTEGER, PARAMETER :: TASK_INIT_VER_Z        = 2  !< task: setup only z-interpolation
  INTEGER, PARAMETER :: TASK_FINALIZE_PZ       = 3  !< task: deallocate pz-interpolation
  !------ interpolation tasks:
  INTEGER, PARAMETER :: TASK_INTP_HOR_LONLAT   = 4  !< task: lon-lat
  INTEGER, PARAMETER :: TASK_INTP_HOR_VERTEX   = 5  !< task: vertex to cell
  INTEGER, PARAMETER :: TASK_INTP_VER_PZLEV    = 6  !< task: vertical p or z-levels
  INTEGER, PARAMETER :: TASK_INTP_VERHOR_PZ_LL = 7  !< task: vertical&horizontal

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
  !  This type is likely to contain more data than really needed four
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

    TYPE(t_patch),             POINTER :: p_patch         => NULL()
    TYPE(t_int_state),         POINTER :: p_int_state     => NULL()
    TYPE(t_nh_state),          POINTER :: p_nh_state      => NULL()
    TYPE(t_nwp_phy_diag),      POINTER :: prm_diag        => NULL()
    TYPE(t_nh_opt_diag),       POINTER :: p_nh_opt_diag   => NULL()
    TYPE(t_nh_pzlev_config),   POINTER :: nh_pzlev_config => NULL()
  END TYPE t_data_input


  !> data necessary for job output.
  !  See also @p t_data_input.
  TYPE t_data_output
    ! pointer for model variable (array)
    TYPE (t_var_list_element), POINTER :: var
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


  !>
  !! Generic interface for registering new post-processing jobs of
  !! different types.
  INTERFACE add_pp_task
    MODULE PROCEDURE add_task_lonlat_intp
    MODULE PROCEDURE add_task_verts2cells
  END INTERFACE

  INTERFACE assign_if_present  ! purely internal
    MODULE PROCEDURE assign_if_present_logical
  END INTERFACE


  !--- MODULE DATA -------------------------------------------------------------------
  TYPE(t_job_queue), POINTER   :: job_queue  !< head of (ordered) job queue


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
  SUBROUTINE pp_scheduler_init(p_patch, p_nh_state, prm_diag, p_nh_opt_diag, &
    &                          nh_pzlev_config, first_output_name_list,      &
    &                          var_lists, nvar_lists)

    TYPE(t_patch),              TARGET,  INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),           TARGET,  INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_nwp_phy_diag),       TARGET,  INTENT(IN)    :: prm_diag(:)
    TYPE(t_nh_opt_diag),        TARGET,  INTENT(INOUT) :: p_nh_opt_diag(:)
    TYPE(t_nh_pzlev_config),    TARGET,  INTENT(IN)    :: nh_pzlev_config(0:max_dom)
    TYPE (t_output_name_list),  POINTER :: first_output_name_list
    TYPE(t_var_list),                    INTENT(IN)    :: var_lists(:)
    INTEGER,                             INTENT(IN)    :: nvar_lists

    ! local variables
    CHARACTER(*), PARAMETER :: routine =  &
      &  TRIM("mo_pp_scheduler:pp_scheduler_init")
    INTEGER                            :: &
      &  jg, ndom, ientr, nblks_c, ierrstat, ivar, i, idx, &
      &  iaxis, vgrid, nlev
    LOGICAL                            :: &
      &  l_jg_active, l_intp_p, l_intp_z, found
    TYPE (t_output_name_list), POINTER :: p_onl
    TYPE(t_job_queue)                  :: task
    INTEGER                            :: shape3d_c(3)
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
   
    !-------------------------------------------------------------
    !-- horizontal interpolation vertices->cells

    ! DEVELOPMENT

    !-------------------------------------------------------------
    !-- horizontal interpolation regular lon-lat grids

    ! DEVELOPMENT

    !-------------------------------------------------------------
    !--- setup of vertical interpolation onto p/z-levels

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

      !-- create data structure for post-processing task
      task%data_input%jg               =  jg           
      task%data_input%p_patch          => p_patch(jg)
      task%data_input%p_nh_state       => p_nh_state(jg)
      task%data_input%p_nh_opt_diag    => p_nh_opt_diag(jg)
      task%data_input%prm_diag         => prm_diag(jg)
      task%data_input%nh_pzlev_config  => nh_pzlev_config(jg)

      !-- check if any name list requests p- or z-level interpolation
      !-- for this domain, collect the list of variables

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
      p_diag_pz         => task%data_input%p_nh_opt_diag%diag_pz
      p_opt_diag_list_z => task%data_input%p_nh_opt_diag%opt_diag_list_z
      p_opt_diag_list_p => task%data_input%p_nh_opt_diag%opt_diag_list_p
      ientr     = 16   ! "entropy" of horizontal slice

      ! predefined array shapes
      nblks_c   = p_patch(jg)%nblks_c
      shape3d_c = (/ nproma, task%data_input%nh_pzlev_config%nzlev, nblks_c /)
      nvars_predef = 0
        
      ! temp         (nproma,nzlev,nblks_c)        
      cf_desc    = t_cf_var('temperature', 'K', 'temperature')
      grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
      nvars_predef = nvars_predef + 1
      varlist_predef(nvars_predef) = "temp"
      CALL add_var( p_opt_diag_list_z, varlist_predef(nvars_predef), p_diag_pz%z_temp, &
        & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE, cf_desc, grib2_desc, &
        & ldims=shape3d_c )
      
      ! pres         (nproma,nzlev,nblks_c)
      nvars_predef = nvars_predef + 1
      varlist_predef(nvars_predef) = "pres"
      cf_desc    = t_cf_var('pressure', 'Pa', 'pressure')
      grib2_desc = t_grib2_var(0, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_opt_diag_list_z, varlist_predef(nvars_predef), p_diag_pz%z_pres, &
        & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE, cf_desc, grib2_desc, &
        & ldims=shape3d_c )

      ! tracer_qv
      nvars_predef = nvars_predef + 1
      varlist_predef(nvars_predef) = "qv"
      cf_desc    = t_cf_var('tracer_qv', 'kg kg-1', 'specific_humidity')
      grib2_desc = t_grib2_var(0, 1, 201, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_opt_diag_list_z, varlist_predef(nvars_predef),          &
        & p_diag_pz%z_tracer_iqv, GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,     &
        & cf_desc, grib2_desc, ldims=shape3d_c)

      ! tot_qv
      nvars_predef = nvars_predef + 1
      varlist_predef(nvars_predef) = "tot_qv"
      cf_desc    = t_cf_var('tot_qv', '','total_specific_humidity')
      grib2_desc = t_grib2_var(0, 6, 6, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_opt_diag_list_z, varlist_predef(nvars_predef),          &
        & p_diag_pz%z_tot_cld_iqv, GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,    &
        & cf_desc, grib2_desc, ldims=shape3d_c)

      IF (l_intp_p) THEN
        shape3d_c = (/ nproma, task%data_input%nh_pzlev_config%nplev, nblks_c /)
        ! GEOPOT
        cf_desc    = t_cf_var('z', 'm2 s-2', 'geopotential')
        grib2_desc = t_grib2_var(0, 3, 4, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_opt_diag_list_p, 'z', p_diag_pz%p_geopot,             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, grib2_desc,      &
          & ldims=shape3d_c )
        ! temp
        cf_desc    = t_cf_var('temperature', 'K', 'temperature')
        grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_opt_diag_list_p, 'temp', p_diag_pz%p_temp,            &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, grib2_desc,      &
          & ldims=shape3d_c )
      END IF

      !-- register interpolation setup as a post-processing task
      task%job_priority = HIGH_PRIORITY
      task%activity     = new_simulation_status(l_output_step=.TRUE.)
      IF (l_intp_p) THEN
        task%job_name = "Init: pz-level interpolation, DOM "//TRIM(int2string(jg))
        task%job_type = TASK_INIT_VER_PZ
      ELSE
        task%job_name = "Init:  z-level interpolation, DOM "//TRIM(int2string(jg))
        task%job_type = TASK_INIT_VER_Z
      END IF
      CALL pp_task_insert(task)

      !-- register clean-up routine as a post-processing task
      task%activity     = new_simulation_status(l_output_step=.TRUE.)
      task%job_name     = "Clean-up: pz-level interpolation, level "//TRIM(int2string(jg))
      task%job_priority = LOW_PRIORITY
      task%job_type     = TASK_FINALIZE_PZ
      CALL pp_task_insert(task)

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
          nlev    =  task%data_input%nh_pzlev_config%nzlev
          vgrid   =  ZAXIS_ALTITUDE
          p_opt_diag_list => p_opt_diag_list_z
        END IF
        IF (iaxis == 2) THEN
          prefix  =  "p-level"
          varlist => pl_varlist
          nvars   =  nvars_pl
          nlev    =  task%data_input%nh_pzlev_config%nplev
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
            ! loop only over model level variables
            IF (var_lists(i)%p%level_type /= level_type_ml) CYCLE         
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
              ! Do not inspect element if it does not contain info for
              ! vertical interpolation
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
              shape3d_c  = (/ nproma, nlev, nblks_c /)

              CALL add_var( p_opt_diag_list, info%name, p_opt_field_r3d, &
                &           info%hgrid, vgrid, info%cf, info%grib2,      &
                &           ldims=shape3d_c, lrestart=.FALSE.,           &
                &           loutput=.TRUE., new_element=new_element)

              !-- add post-processing task
              task%job_name        =  &
                &  TRIM(prefix)//" interp. "//TRIM(info%name)  &
                &  //", DOM "//TRIM(int2string(jg))
              task%job_priority    =  DEFAULT_PRIORITY
              task%job_type        =  TASK_INTP_VER_PZLEV
              task%activity        =  new_simulation_status(l_output_step=.TRUE.)
              task%data_input%var  => element%field       ! set input variable
              task%data_output%var => new_element%field   ! set output variable
              CALL pp_task_insert(task)
            
              found = .TRUE.

            ENDDO ! loop over vlist "i"
          ENDDO ! i = 1,nvar_lists
        
          ! Check that at least one element with this name has been found
          IF(.NOT. found) &
            CALL finish(routine,'No feasible variable found: '//TRIM(varlist(ivar)))
        END DO ! ivar
      END DO ! height/pressure axis
    END DO DOM_LOOP  ! jg
    
    DEALLOCATE(pl_varlist, hl_varlist, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    IF (dbg_level > 5)  CALL message(routine, "Done")
    
  END SUBROUTINE pp_scheduler_init


  !---------------------------------------------------------------
  !> Register interpolation task: lon-lat interpolation
  !
  SUBROUTINE add_task_lonlat_intp()
    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_pp_scheduler:add_task_lonlat_intp")
    TYPE(t_job_queue) :: task

    ! DEVELOPMENT
    CALL message(routine, "end")
    CALL pp_task_insert(task)
  END SUBROUTINE add_task_lonlat_intp


  !---------------------------------------------------------------
  !> Register interpolation task: vertex -> cells
  !
  SUBROUTINE add_task_verts2cells(in_var, out_var, p_patch, p_int_state, job_priority)
    TYPE (t_var_list_element), POINTER :: in_var, out_var
    TYPE(t_patch),             POINTER :: p_patch
    TYPE(t_int_state),         POINTER :: p_int_state
    INTEGER,                  OPTIONAL, INTENT(IN) :: job_priority
    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_pp_scheduler:add_task_verts2cells")
    TYPE(t_job_queue) :: task

    ! meta data 
    IF (PRESENT(job_priority)) THEN
      task%job_priority = job_priority
    ELSE
      task%job_priority = DEFAULT_PRIORITY
    END IF
    task%job_name     = TRIM("verts2cells: "//in_var%info%name)
    task%job_type     = TASK_INTP_HOR_VERTEX
    CALL message(routine, TRIM(task%job_name))

    ! input
    task%data_input%var          => in_var
    task%data_input%p_patch      => p_patch
    task%data_input%p_int_state  => p_int_state
    ! output
    task%data_output%var         => out_var

    CALL pp_task_insert(task)
  END SUBROUTINE add_task_verts2cells


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
      CASE ( TASK_INTP_HOR_VERTEX )
        CALL pp_task_vertex2cell(ptr_task)
      CASE ( TASK_INTP_VERHOR_PZ_LL )
        CALL pp_task_pz_lonlat(ptr_task)
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

    ! DEVELOPMENT

  END SUBROUTINE pp_task_lonlat


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
      &  vert_intp_type, vert_intp_method, jg, &
      &  in_var_idx, out_var_idx, nlev,        &
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


  !---------------------------------------------------------------
  !> Performs interpolation of a vertex-based variable onto cells.
  !
  !  This is only a wrapper for the corresponding routines from the
  !  interpolation module.
  SUBROUTINE pp_task_vertex2cell(ptr_task)
    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_pp_scheduler:pp_task_vertex2cell")
    INTEGER                            :: nlev, &
      &  in_var_idx, out_var_idx
    TYPE(t_patch),             POINTER :: p_patch
    TYPE(t_int_state),         POINTER :: p_int_state
    TYPE (t_var_list_element), POINTER :: in_var, out_var
    TYPE(t_var_metadata),      POINTER :: p_info

    p_patch     => ptr_task%data_input%p_patch      ! patch
    nlev        =  p_patch%nlev                     ! no. of levels
    p_int_state => ptr_task%data_input%p_int_state  ! interpolation coeffs

    p_info      => ptr_task%data_input%var%info
    in_var      => ptr_task%data_input%var
    out_var     => ptr_task%data_output%var

    in_var_idx        = 1
    if (ptr_task%data_input%var%info%lcontained) &
      in_var_idx = ptr_task%data_input%var%info%ncontained
    out_var_idx       = 1
    if (ptr_task%data_output%var%info%lcontained) &
      out_var_idx = ptr_task%data_output%var%info%ncontained

    IF (p_info%ndims /= 3) &
      & CALL finish(routine, "Wrong number of variables dimensions!")
      
    CALL verts2cells_scalar( in_var%r_ptr(:,:,:,in_var_idx,1), p_patch,  &
      &                      p_int_state%verts_aw_cells,        &
      &                      out_var%r_ptr(:,:,:,out_var_idx,1), 1, nlev )

  END SUBROUTINE pp_task_vertex2cell


  !---------------------------------------------------------------
  !> Performs interpolation onto pz-levels _and_ lon-lat grid.
  !
  !  This is only a wrapper for the corresponding routines from the
  !  interpolation module.
  SUBROUTINE pp_task_pz_lonlat(ptr_task)
    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_pp_scheduler:pp_task_p_lonlat")

  END SUBROUTINE pp_task_pz_lonlat


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
  SUBROUTINE pp_task_insert(task)
    TYPE(t_job_queue), INTENT(IN) :: task
    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_pp_scheduler:pp_task_insert")
    TYPE(t_job_queue), POINTER :: tmp, nb_left, element
    INTEGER                    :: ierrstat

    IF (dbg_level > 5) THEN
      WRITE(message_text,*) "Inserting pp task '", TRIM(task%job_name), "'"
      CALL message(routine, TRIM(message_text))
    END IF

    ! find the correct position in list:
    tmp     => job_queue
    nb_left => NULL()
    DO
      IF (.NOT. ASSOCIATED(tmp))                EXIT
      IF (tmp%job_priority > task%job_priority) EXIT
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
    element      = task
    element%next => tmp

  END SUBROUTINE pp_task_insert


  !------------------------------------------------------------------------------------------------


  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for "t_simulation" variables
  ! 
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION new_simulation_status(l_output_step, l_first_step, l_last_step)  &
    RESULT(sim_status)

    TYPE(t_simulation_status) :: sim_status
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  l_output_step, l_first_step, l_last_step

    ! set default values
    sim_status%status_flags(1) = .FALSE.
    sim_status%status_flags(2) = .FALSE.
    sim_status%status_flags(3) = .FALSE.
    ! supersede with user definitions
    CALL assign_if_present(sim_status%status_flags(1), l_output_step)
    CALL assign_if_present(sim_status%status_flags(2), l_first_step)
    CALL assign_if_present(sim_status%status_flags(3), l_last_step)

  END FUNCTION new_simulation_status


  !
  ! private routine to assign values if actual parameters are present
  !
  SUBROUTINE assign_if_present_logical (y,x)
    LOGICAL, INTENT(inout)        :: y
    LOGICAL, INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_logical

END MODULE mo_pp_scheduler
