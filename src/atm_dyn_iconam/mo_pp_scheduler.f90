! Option directive to avoid a possible SX compiler bug with OpenMP + MPI;
! assumes that pointer references are overlapped in optimization.
! disable moving invariant expressions outside of a loop by the compiler.
! [2012-11-08, F. Prill, DWD / 2012-11-27, J. Beismann, NEC]
!option! -O overlap
!option! -O nomove
!>
!! Scheduler for internal post-processing.
!!
!! This module manages a "job queue" for internal post-processing
!! tasks on the compute PEs.
!!
!! For example, the interpolation of model variables onto lon-lat
!! fields constitutes such a task. Allocating and computing only those
!! fields which are required for output (or as intermediate results)
!! can save memory and computing time. 

!! Jobs are processed according to user-defined priority levels. In
!! principle, all tasks with the same job priority could be processed
!! simultaneously (with OpenMP).
!!
!! List of post-processing tasks implemented so far:
!!
!! Type of task                        Execution priority
!! -------------------------------------------------------------------
!!
!! z/p/i interpolation setup           HIGH_PRIORITY
!! compute vertical velocity           DEFAULT_PRIORITY0
!! compute rel. humidity               DEFAULT_PRIORITY0
!! compute mean sea level pressure     DEFAULT_PRIORITY0
!! vertical interpolation              DEFAULT_PRIORITY1
!! horizontal interpolation edge->cell DEFAULT_PRIORITY2
!! horizontal synchronization          DEFAULT_PRIORITY3
!! horizontal interpolation            DEFAULT_PRIORITY4
!! z/p/i interpolation clean-up        LOW_PRIORITY
!!
!!
!! @author F. Prill, DWD
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2012-03-01)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
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
  USE mo_impl_constants,          ONLY: SUCCESS, HINTP_TYPE_NONE, max_dom, max_var_ml,      &
    &                                   max_var_pl, max_var_hl, max_var_il, TASK_NONE,      &
    &                                   TASK_INIT_VER_Z, TASK_INIT_VER_P, TASK_INIT_VER_I,  &
    &                                   TASK_FINALIZE_IPZ, TASK_INTP_HOR_LONLAT,            &
    &                                   TASK_INTP_VER_PLEV, TASK_INTP_SYNC, TASK_INTP_MSL,  &
    &                                   TASK_COMPUTE_RH, TASK_INTP_VER_ZLEV,                &
    &                                   TASK_INTP_VER_ILEV, TASK_INTP_EDGE2CELL,            &
    &                                   max_phys_dom, UNDEF_TIMELEVEL, ALL_TIMELEVELS,      &
    &                                   vname_len, TASK_COMPUTE_OMEGA
  USE mo_model_domain,            ONLY: t_patch, p_patch, p_phys_patch
  USE mo_var_list,                ONLY: add_var, get_all_var_names,                         &
    &                                   nvar_lists, var_lists, get_var_name,                &
    &                                   get_var_timelevel
  USE mo_var_list_element,        ONLY: t_var_list_element, level_type_ml,                  &
    &                                   level_type_pl, level_type_hl, level_type_il
  USE mo_var_metadata_types,      ONLY: t_var_metadata, t_vert_interp_meta,                 &
    &                                   VINTP_TYPE_LIST, VARNAME_LEN, t_post_op_meta
  USE mo_var_metadata,            ONLY: create_hor_interp_metadata, vintp_types,            &
    &                                   vintp_type_id
  USE mo_intp_data_strc,          ONLY: t_int_state, lonlat_grid_list,                      &
    &                                   t_lon_lat_intp, p_int_state,                        &
    &                                   MAX_LONLAT_GRIDS
  USE mo_nonhydro_state,          ONLY: p_nh_state
  USE mo_opt_diagnostics,         ONLY: t_nh_diag_pz, t_nh_opt_diag, p_nh_opt_diag
  USE mo_nwp_phy_types,           ONLY: t_nwp_phy_diag
  USE mo_nwp_phy_state,           ONLY: prm_diag
  USE mo_nh_pzlev_config,         ONLY: t_nh_pzlev_config, nh_pzlev_config
  USE mo_name_list_output_config, ONLY: first_output_name_list
  USE mo_name_list_output_types,  ONLY: t_output_name_list
  USE mo_parallel_config,         ONLY: nproma
  USE mo_cf_convention,           ONLY: t_cf_var
  USE mo_grib2,                   ONLY: t_grib2_var
  USE mo_util_string,             ONLY: int2string, remove_duplicates,                      &
    &                                   difference, toupper, tolower
  USE mo_cdi_constants,           ONLY: GRID_CELL, GRID_REFERENCE,                          &
    &                                   GRID_UNSTRUCTURED_CELL, ZA_ALTITUDE,                &
    &                                   ZA_PRESSURE, GRID_REGULAR_LONLAT,                   &
    &                                   GRID_UNSTRUCTURED_EDGE,                             &
    &                                   DATATYPE_FLT32, DATATYPE_PACK16, ZA_ISENTROPIC,     &
    &                                   is_2d_field
  USE mo_linked_list,             ONLY: t_var_list, t_list_element, find_list_element
  USE mo_grid_config,             ONLY: n_dom
  USE mo_pp_tasks,                ONLY: pp_task_lonlat, pp_task_sync, pp_task_ipzlev_setup, &
    &                                   pp_task_ipzlev, pp_task_compute_field,              &
    &                                   pp_task_intp_msl, pp_task_edge2cell,                & 
    &                                   t_data_input, t_data_output,                        &
    &                                   t_simulation_status, t_job_queue, job_queue,        &
    &                                   MAX_NAME_LENGTH, HIGH_PRIORITY,                     &
    &                                   DEFAULT_PRIORITY0, DEFAULT_PRIORITY1,               &
    &                                   DEFAULT_PRIORITY2, DEFAULT_PRIORITY3,               &
    &                                   DEFAULT_PRIORITY4, LOW_PRIORITY, dbg_level
  USE mo_fortran_tools,           ONLY: assign_if_present


  IMPLICIT NONE

  ! interface definition
  PRIVATE

  CHARACTER(len=*), PARAMETER :: &
    &  version = '$Id$'

  ! functions and subroutines
  PUBLIC :: pp_scheduler_init
  PUBLIC :: pp_scheduler_process
  PUBLIC :: pp_scheduler_finalize
  PUBLIC :: new_simulation_status

  ! some constants (for better readability):
  CHARACTER(*), PARAMETER :: vn_name = TRIM("vn")

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
    CHARACTER(*), PARAMETER :: routine =  TRIM("mo_pp_scheduler:pp_scheduler_init")
    INTEGER                          :: jg, i
    TYPE(t_list_element), POINTER    :: element, element_pres

    if (dbg_level > 5)  CALL message(routine, "Enter")

    !-------------------------------------------------------------
    !--- setup of optional diagnostic fields updated by the 
    !    post-processing scheduler

    !- loop over model level variables
    DO i = 1,nvar_lists
      jg = var_lists(i)%p%patch_id         

      element => NULL()
      DO
        IF(.NOT.ASSOCIATED(element)) THEN
          element => var_lists(i)%p%first_list_element
        ELSE
          element => element%next_list_element
        ENDIF
        IF(.NOT.ASSOCIATED(element)) EXIT

        IF (element%field%info%l_pp_scheduler_task /= TASK_NONE) THEN

          IF (dbg_level > 5)  &
            & CALL message(routine, "Inserting pp task: "//TRIM(element%field%info%name))

          SELECT CASE(element%field%info%l_pp_scheduler_task)
          CASE (TASK_COMPUTE_RH) 
            ! relative humidity
            CALL pp_scheduler_register( name=element%field%info%name, jg=jg, p_out_var=element, &
              &                         l_init_prm_diag=l_init_prm_diag, job_type=TASK_COMPUTE_RH )
            !
          CASE (TASK_COMPUTE_OMEGA)
            ! vertical velocity
            CALL pp_scheduler_register( name=element%field%info%name, jg=jg, p_out_var=element, &
              &                         l_init_prm_diag=l_init_prm_diag, job_type=TASK_COMPUTE_OMEGA )
            !
          CASE (TASK_INTP_MSL)   
            ! mean sea level pressure
            !
            ! find the standard pressure field:
            element_pres => find_list_element (p_nh_state(jg)%diag_list, 'pres')
            IF (ASSOCIATED (element)) THEN
              ! register task for interpolation to z=0:
              CALL pp_scheduler_register( name=element%field%info%name, jg=jg, p_out_var=element, &
                &                         job_type=TASK_INTP_MSL,                                 &
                &                         l_init_prm_diag=l_init_prm_diag,                        &
                &                         opt_p_in_var=element_pres)
            END IF
            CASE DEFAULT
            CALL finish(routine, "Unknown pp task type!")
          END SELECT
        END IF
        
      ENDDO ! loop over vlist "i"
    ENDDO ! i = 1,nvar_lists

    !-------------------------------------------------------------
    !--- setup of vertical interpolation onto i/p/z-levels

    CALL pp_scheduler_init_ipz(l_init_prm_diag)

    !-------------------------------------------------------------
    !-- horizontal interpolation regular lon-lat grids

    CALL pp_scheduler_init_lonlat

    IF (dbg_level > 5)  CALL message(routine, "Done")
    
  END SUBROUTINE pp_scheduler_init


  !---------------------------------------------------------------
  !> (Internal) Utility routine, add a new "vn" field for a given
  !  axis, based on the meta-data of the standard "vn". - This is done
  !  for all available time levels.
  !
  SUBROUTINE init_vn_horizontal(ll_grid_id, lev_type)
    INTEGER,          INTENT(IN)  :: ll_grid_id      !< lon-lat grid number
    INTEGER,          INTENT(IN)  :: lev_type        !< level type: p/z/i/m
    ! local variables
    CHARACTER(*), PARAMETER :: routine =  TRIM("mo_pp_scheduler:init_vn_horizontal")
    TYPE(t_list_element), POINTER :: element_u, element_v, element, new_element, new_element_2
    INTEGER                       :: i, shape3d_ll(3), shape3d_e(3), nblks_lonlat, &
      &                              nblks_e, nlev, jg, tl
    TYPE(t_job_queue),    POINTER :: task
    TYPE(t_var_metadata), POINTER :: info
    REAL(wp),             POINTER :: p_opt_field_r3d(:,:,:)
    CHARACTER(len=varname_len)    :: name
    CHARACTER(len=4)              :: suffix
    TYPE(t_cf_var)                :: cf
    TYPE(t_grib2_var)             :: grib2
    TYPE(t_post_op_meta)          :: post_op
    TYPE(t_var_list), POINTER     :: dst_varlist     !< destination variable list
    TYPE (t_lon_lat_intp), POINTER:: ptr_int_lonlat
    CHARACTER(LEN=1)              :: prefix

    IF (dbg_level > 5)  CALL message(routine, "Enter")
    
    ! initialize "new_element" pointer (cf. NEC compiler bugs DWD0121
    ! and DWD0123 for hybrid parallelization)
    new_element   => NULL()
    new_element_2 => NULL()

    !- loop over model level variables
    ! Note that there are several "vn" variables with different time
    ! levels, we just add unconditionally all
    DO i = 1,nvar_lists
      jg = var_lists(i)%p%patch_id

      SELECT CASE(lev_type)
      CASE (level_type_ml)
        dst_varlist => p_nh_opt_diag(jg)%opt_diag_list         
        prefix = "m"
      CASE (level_type_pl)
        dst_varlist => p_nh_opt_diag(jg)%opt_diag_list_p
        prefix = "p"
      CASE (level_type_hl)
        dst_varlist => p_nh_opt_diag(jg)%opt_diag_list_z
        prefix = "z"
      CASE (level_type_il)
        dst_varlist => p_nh_opt_diag(jg)%opt_diag_list_i
        prefix = "i"
      END SELECT

      ! Do not inspect lists which are disabled for output
      IF (.NOT. var_lists(i)%p%loutput) CYCLE
      ! loop only over model level variables
      IF (var_lists(i)%p%vlevel_type /= lev_type) CYCLE         
      ! loop only over variables of where domain was requested
      IF (.NOT. lonlat_grid_list(ll_grid_id)%l_dom(jg)) CYCLE

      ! now, search for "vn" in the variable list:
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

        ! Check for matching name
        IF (vn_name /= TRIM(tolower(get_var_name(element%field)))) CYCLE

        ! get time level
        tl = get_var_timelevel(element%field)
        suffix = ""
        IF (tl /= -1)  WRITE (suffix,'(".TL",i1)') tl

        !- find existing variables "u", "v" (for copying the meta-data):
        element_u => find_list_element (p_nh_state(jg)%diag_list, "u")
        element_v => find_list_element (p_nh_state(jg)%diag_list, "v")
        
        !- predefined array shapes
        nlev = element%field%info%used_dimensions(2)
        nblks_e   = p_patch(jg)%nblks_e
        ptr_int_lonlat => lonlat_grid_list(ll_grid_id)%intp(jg)
        nblks_lonlat   =  (ptr_int_lonlat%nthis_local_pts - 1)/nproma + 1
        shape3d_ll = (/ nproma, nlev, nblks_lonlat /)
        shape3d_e  = (/ nproma, nlev, nblks_e /)

        !-- create new cell-based variables "u", "v" on lon-lat grid
        !   for the same time level
        IF (dbg_level > 8) &
          CALL message(routine, "horizontal interpolation: create u/v variables on lon-lat grid")
        name    = TRIM(get_var_name(element_u%field))//suffix
        cf      = element_u%field%info%cf
        grib2   = element_u%field%info%grib2
        post_op = element_u%field%info%post_op
        CALL add_var( dst_varlist, TRIM(name), p_opt_field_r3d,                           &
          & GRID_REGULAR_LONLAT, info%vgrid, cf, grib2,                                   &
          & ldims=shape3d_ll, lrestart=.FALSE., in_group=element_u%field%info%in_group,   &
          & new_element=new_element, loutput=.TRUE., post_op=post_op )

        name    = TRIM(get_var_name(element_v%field))//suffix
        cf      = element_v%field%info%cf
        grib2   = element_v%field%info%grib2
        post_op = element_v%field%info%post_op
        CALL add_var( dst_varlist, TRIM(name), p_opt_field_r3d,                           &
          & GRID_REGULAR_LONLAT, info%vgrid, cf, grib2,                                   &
          & ldims=shape3d_ll, lrestart=.FALSE., in_group=element_v%field%info%in_group,   &
          & new_element=new_element_2, loutput=.TRUE., post_op=post_op )

        ! link these new variables to the lon-lat grid:
        new_element%field%info%hor_interp%lonlat_id   = ll_grid_id
        new_element_2%field%info%hor_interp%lonlat_id = ll_grid_id

        !-- create and add post-processing task
        task => pp_task_insert(DEFAULT_PRIORITY4)
        WRITE (task%job_name, *) "horizontal interp. ",TRIM(info%name),", ",prefix//"-levels", ", DOM ",jg
        IF (dbg_level > 8) CALL message(routine, task%job_name)
        task%data_input%p_nh_state          => NULL()
        task%data_input%prm_diag            => NULL()
        task%data_input%nh_pzlev_config     => NULL()
        task%data_input%jg                  =  jg           
        task%data_input%p_patch             => p_patch(jg)
        task%data_input%p_nh_opt_diag       => p_nh_opt_diag(jg)
        task%data_input%p_int_state         => p_int_state(jg)
        task%job_type                       =  TASK_INTP_HOR_LONLAT
        task%activity                       =  new_simulation_status(l_output_step=.TRUE.)
        task%activity%ldom_active(jg)       =  .TRUE.
        task%activity%i_timelevel(jg)       =  get_var_timelevel(element%field)
        task%data_input%var                 => element%field       ! set input variable
        task%data_output%var                => new_element%field   ! set output variable "u"
        task%data_output%var_2              => new_element_2%field ! set output variable "v"
      END DO
    END DO
    if (dbg_level > 5)  CALL message(routine, "Done")

  END SUBROUTINE init_vn_horizontal


  !---------------------------------------------------------------
  !> Setup of lon-lat interpolation tasks.
  !
  SUBROUTINE pp_scheduler_init_lonlat

    ! local variables
    CHARACTER(*), PARAMETER :: routine =  TRIM("mo_pp_scheduler:pp_scheduler_init_lonlat")
    INTEGER                               :: &
      &  jg, ndom, ierrstat, ivar, i, j, nvars_ll, &
      &  nblks_lonlat, ilev_type, max_var, ilev, n_uv_hrz_intp
    LOGICAL                               :: found, l_horintp
    TYPE (t_output_name_list), POINTER    :: p_onl
    TYPE(t_job_queue),         POINTER    :: task
    TYPE(t_var_list),          POINTER    :: p_opt_diag_list
    REAL(wp), POINTER                     :: p_opt_field_r3d(:,:,:)
    INTEGER,  POINTER                     :: p_opt_field_i3d(:,:,:)
    TYPE(t_list_element),      POINTER    :: element, new_element
    CHARACTER(LEN=vname_len),  POINTER    :: varlist(:)
    INTEGER, ALLOCATABLE                  :: ll_vargrid(:)
    CHARACTER(LEN=vname_len), ALLOCATABLE :: ll_varlist(:)
    INTEGER, ALLOCATABLE                  :: ll_varlevs(:)
    CHARACTER(LEN=vname_len)              :: vname
    TYPE(t_var_metadata),      POINTER    :: info
    INTEGER                               :: var_shape(5)
    TYPE (t_lon_lat_intp),     POINTER    :: ptr_int_lonlat
    INTEGER                               :: uv_hrz_intp_grid(4*MAX_LONLAT_GRIDS), &
      &                                      uv_hrz_intp_levs(4*MAX_LONLAT_GRIDS)

    if (dbg_level > 5)  CALL message(routine, "Enter")

    ! initialize "new_element" pointer (cf. NEC compiler bugs DWD0121
    ! and DWD0123 for hybrid parallelization)
    new_element   => NULL()

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
    nvars_ll            =  0
    n_uv_hrz_intp       =  0
    uv_hrz_intp_grid(:) = -1
    uv_hrz_intp_levs(:) = -1

    NML_LOOP : DO
      IF (.NOT.ASSOCIATED(p_onl)) EXIT NML_LOOP

      ! do the same for the three level types: model levels (ml),
      ! height levels (hl), pressure levels (pl) and isentropic levels (il):
      DO ilev=1,4
   
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
        CASE (4)
          varlist => p_onl%il_varlist
          ilev_type  =  level_type_il
          max_var    =  max_var_il
        END SELECT
        IF (varlist(1) == ' ') CYCLE
      
        ! Selection criterion: 
        ! - lon-lat interpolation is requested
        IF (p_onl%remap == 1) THEN

          DO ivar=1,max_var
            IF (varlist(ivar) == ' ') CYCLE
            nvars_ll=nvars_ll+1
            ll_varlist(nvars_ll) = varlist(ivar)
            ll_vargrid(nvars_ll) = p_onl%lonlat_id
            ll_varlevs(nvars_ll) = ilev_type
            ! return a special flag, if var name matches "u" or "v":
            IF ((TRIM(ll_varlist(nvars_ll)) == "u") .OR. &
              & (TRIM(ll_varlist(nvars_ll)) == "v")) THEN
              ! check if this lon-lat grid has not yet been
              ! registered (because we may specify "u" AND "v"):
              found = .FALSE.
              DO j=1,n_uv_hrz_intp
                IF ((uv_hrz_intp_grid(n_uv_hrz_intp) == p_onl%lonlat_id) .AND. &
                  & (uv_hrz_intp_levs(n_uv_hrz_intp) == ilev_type)) THEN
                  found = .TRUE.
                END IF
              END DO
              IF (.NOT. found) THEN
                n_uv_hrz_intp = n_uv_hrz_intp + 1
                uv_hrz_intp_grid(n_uv_hrz_intp) = p_onl%lonlat_id
                uv_hrz_intp_levs(n_uv_hrz_intp) = ilev_type
              END IF
            END IF
          END DO

        END IF
      END DO
      p_onl => p_onl%next
    END DO NML_LOOP

    !-- take care of the special case of "u", "v", where we take "vn"
    !   as the source for horizontal interpolation:
    IF (n_uv_hrz_intp > 0) THEN
      DO i=1,n_uv_hrz_intp
        ! insert post-processing tasks for "vn":
        CALL init_vn_horizontal(uv_hrz_intp_grid(i), uv_hrz_intp_levs(i))
      END DO
    END IF

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
        ! match (p/z/i/model levels):
        IF (var_lists(i)%p%vlevel_type/=ll_varlevs(ivar)) CYCLE LIST_LOOP
        ! loop only over variables on requested domains:
        jg = var_lists(i)%p%patch_id
        IF (.NOT. lonlat_grid_list(ll_vargrid(ivar))%l_dom(jg)) CYCLE LIST_LOOP
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

          ! "u", "v" are processed separately, see above.
          IF ((TRIM(vname) == "u") .OR. (TRIM(vname) == "v"))  CYCLE VAR_LOOP

          ! Check for matching name (take care of suffix of
          ! time-dependent variables):
          IF (TRIM(vname) /= TRIM(tolower(get_var_name(element%field)))) CYCLE VAR_LOOP

          IF (info%hgrid /= GRID_UNSTRUCTURED_CELL)  CYCLE VAR_LOOP

          ! Found it, add it to the variable list of optional
          ! diagnostics       
          SELECT CASE(ll_varlevs(ivar))
          CASE (level_type_ml)
            p_opt_diag_list => p_nh_opt_diag(jg)%opt_diag_list         
          CASE (level_type_pl)
            p_opt_diag_list => p_nh_opt_diag(jg)%opt_diag_list_p
          CASE (level_type_hl)
            p_opt_diag_list => p_nh_opt_diag(jg)%opt_diag_list_z
          CASE (level_type_il)
            p_opt_diag_list => p_nh_opt_diag(jg)%opt_diag_list_i
          END SELECT

          ! set local values for "nblks" and "npromz"
          ptr_int_lonlat => lonlat_grid_list(ll_vargrid(ivar))%intp(jg)
          nblks_lonlat   =  (ptr_int_lonlat%nthis_local_pts - 1)/nproma + 1
          var_shape      =  info%used_dimensions(:)
          IF (is_2d_field(info%vgrid) .AND. (info%ndims /= 2)) THEN
            CALL finish(routine, "Inconsistent dimension info: "//TRIM(info%name)//"!")
          END IF
          IF (is_2d_field(info%vgrid)) THEN
            var_shape(2:3)   =  (/ 1, nblks_lonlat /)
          ELSE
            var_shape(3)     =  nblks_lonlat
          END IF

          SELECT CASE (info%hgrid)
          CASE (GRID_UNSTRUCTURED_CELL)
            !--- REAL fields
            IF (ASSOCIATED(element%field%r_ptr)) THEN
              CALL add_var( p_opt_diag_list, info%name, p_opt_field_r3d,          &
                &           GRID_REGULAR_LONLAT, info%vgrid, info%cf, info%grib2, &
                &           ldims=var_shape, lrestart=.FALSE.,                    &
                &           loutput=.TRUE., new_element=new_element,              &
                &           isteptype=info%isteptype,                             &
                &           hor_interp=create_hor_interp_metadata(                &
                &               hor_intp_type=HINTP_TYPE_NONE ),                  &
                &           post_op=info%post_op )
            END IF
            !--- INTEGER fields
            IF (ASSOCIATED(element%field%i_ptr)) THEN
              CALL add_var( p_opt_diag_list, info%name, p_opt_field_i3d,          &
                &           GRID_REGULAR_LONLAT, info%vgrid, info%cf, info%grib2, &
                &           ldims=var_shape, lrestart=.FALSE.,                    &
                &           loutput=.TRUE., new_element=new_element,              &
                &           isteptype=info%isteptype,                             &
                &           hor_interp=create_hor_interp_metadata(                &
                &               hor_intp_type=HINTP_TYPE_NONE ),                  &
                &           post_op=info%post_op )
            END IF
          CASE DEFAULT
            CALL finish(routine, "Unsupported grid type!")
          END SELECT
  
          ! link this new variable to the lon-lat grid:
          new_element%field%info%hor_interp%lonlat_id = ll_vargrid(ivar)
            
          !-- create and add post-processing task
          task => pp_task_insert(DEFAULT_PRIORITY4)
          WRITE (task%job_name, *) "horizontal interp. ",TRIM(info%name),", DOM ",jg
          task%data_input%p_nh_state      => NULL()
          task%data_input%prm_diag        => NULL()
          task%data_input%nh_pzlev_config => NULL()
          task%data_input%jg              =  jg           
          task%data_input%p_patch         => p_patch(jg)
          task%data_input%p_nh_opt_diag   => p_nh_opt_diag(jg)
          task%data_input%p_int_state     => p_int_state(jg)
          task%job_type                   =  TASK_INTP_HOR_LONLAT
          task%activity                   =  new_simulation_status(l_output_step=.TRUE.)
          task%activity%ldom_active(jg)   =  .TRUE.
          task%activity%i_timelevel(jg)   =  get_var_timelevel(element%field)
          task%data_input%var             => element%field       ! set input variable
          task%data_output%var            => new_element%field   ! set output variable
          
          ! Flag. Denotes that at least one interpolation task has
          ! been created.
          l_horintp = .TRUE.
          found     = .TRUE.
          
        ENDDO VAR_LOOP ! loop over vlist "i"
      ENDDO LIST_LOOP ! i = 1,nvar_lists
    END DO ! ivar

    DEALLOCATE(ll_varlist, ll_vargrid, ll_varlevs, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    ! If at least one interpolation task has been created, we add a
    ! setup task which synchronizes the halo regions:
    IF (l_horintp) THEN
      IF (dbg_level >= 10) &
        CALL message(routine, "Creating synchronization task for horizontal interpolation.")
      task => pp_task_insert(DEFAULT_PRIORITY3)
      WRITE (task%job_name, *) "horizontal interp. SYNC"
      task%job_type = TASK_INTP_SYNC
      task%activity = new_simulation_status(l_output_step=.TRUE.)
      task%activity%ldom_active(:)  = .FALSE. ! i.e. no domain-wise (in-)activity 
      task%activity%i_timelevel(:)  = ALL_TIMELEVELS
    END IF
    IF (dbg_level > 5)  CALL message(routine, "Done")
    
  END SUBROUTINE pp_scheduler_init_lonlat


  !---------------------------------------------------------------
  !> (Internal) Utility routine, collecting variable names from output
  !  namelist.
  !
  SUBROUTINE collect_output_variables(jg, vintp_name, max_var, &
    &                                 nvars, l_intp, var_names, l_uv_vertical_intp)
    INTEGER, INTENT(IN)      :: jg                              !< current domain
    CHARACTER(LEN=*), INTENT(IN) :: vintp_name                  !< "P", "Z", "I"
    INTEGER, INTENT(IN)      :: max_var                         !< maximum no. of variables
    INTEGER, INTENT(OUT)     :: nvars                           !< actual no. of variables
    LOGICAL, INTENT(OUT)     :: l_intp                          !< Flag. .FALSE. if there is no variable
    CHARACTER(LEN=vname_len), INTENT(INOUT) :: var_names(:)     !< list of variable names (strings)
    LOGICAL, INTENT(OUT)     :: l_uv_vertical_intp              !< Flag. .TRUE., if "u" or "v" contained
    ! local variables
    CHARACTER(*), PARAMETER :: routine =  &
      &  TRIM("mo_pp_scheduler::collect_output_variables")
    TYPE (t_output_name_list), POINTER :: p_onl
    LOGICAL :: l_jg_active, vert_intp_type(SIZE(VINTP_TYPE_LIST))
    INTEGER :: ivar, iphys_dom
    CHARACTER(LEN=vname_len), POINTER :: nml_varlist(:)         !< varlist (hl/ml/pl/il) in output_nml namelist

    vert_intp_type(:) = vintp_types(TRIM(vintp_name))

    l_uv_vertical_intp = .FALSE.
    p_onl => first_output_name_list
    nvars = 0
    var_names(:) = " "
    l_intp = .FALSE.
    NML_LOOP : DO
      IF (.NOT.ASSOCIATED(p_onl)) EXIT NML_LOOP
      SELECT CASE (TRIM(toupper(vintp_name)))
      CASE ("Z")
        nml_varlist => p_onl%hl_varlist
      CASE ("P")
        nml_varlist => p_onl%pl_varlist
      CASE ("I")
        nml_varlist => p_onl%il_varlist
      CASE DEFAULT
        CALL finish(routine, "Internal error!")
      END SELECT

      IF (dbg_level >= 15)  WRITE (0,*) nml_varlist 

      ! If dom(:) was not specified in namelist input, it is set
      ! completely to -1.  In this case all domains are searched:
      l_jg_active = (p_onl%dom(1) <= 0)
      DO iphys_dom=1,max_phys_dom
        IF (p_onl%dom(iphys_dom) > 0) THEN
          l_jg_active = l_jg_active .OR. (jg == p_phys_patch(p_onl%dom(iphys_dom))%logical_id)
        END IF
      END DO

      ! Selection criteria: 
      ! - domain is requested
      ! - "Z"/"P"/"I"-level interpolation is requested
      IF (l_jg_active .AND. (nml_varlist(1) /= ' ')) THEN
        l_intp = .TRUE.
        DO ivar=1,max_var
          IF (nml_varlist(ivar) == ' ') CYCLE
          nvars=nvars+1
          var_names(nvars) = nml_varlist(ivar)
          ! return a special flag, if var name matches "u" or "v":
          l_uv_vertical_intp = l_uv_vertical_intp .OR. &
            &                  (TRIM(var_names(nvars)) == "u") .OR. &
            &                  (TRIM(var_names(nvars)) == "v")
        END DO
      END IF
      p_onl => p_onl%next
    END DO NML_LOOP
  END SUBROUTINE collect_output_variables


  !---------------------------------------------------------------
  !> (Internal) Utility routine, add a new variable field for a given
  !  axis, based on the meta-data of an existing variable field (which
  !  is defined on model/half levels).
  !
  SUBROUTINE copy_variable(name, src_varlist, dst_axis, shape3d, ptr, dst_varlist)
    CHARACTER(len=*),          INTENT(IN)    :: name        !< name of variable
    TYPE(t_var_list),          INTENT(IN)    :: src_varlist !< source variable list
    INTEGER,                   INTENT(IN)    :: dst_axis    !< destination axis
    INTEGER,                   INTENT(IN)    :: shape3d(3)  !< shape of variable field
    REAL(wp),         POINTER                :: ptr(:,:,:)  !< reference to field
    TYPE(t_var_list), POINTER                :: dst_varlist !< destination variable list
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_pp_scheduler:copy_variable")
    TYPE(t_list_element), POINTER :: element
     
    ! find existing variable
    element => find_list_element (src_varlist, TRIM(name))
    IF (.NOT. ASSOCIATED (element)) CALL finish(routine, "Variable not found!")
    ! add new variable, copy the meta-data from the existing variable
    CALL add_var( dst_varlist, TRIM(name), ptr, element%field%info%hgrid, dst_axis,  &
      &           element%field%info%cf, element%field%info%grib2, ldims=shape3d,    &
      &           post_op=element%field%info%post_op, loutput=.TRUE., lrestart=.FALSE. )
  END SUBROUTINE copy_variable


  !---------------------------------------------------------------
  !> (Internal) Utility routine, add a new "vn" field for a given
  !  axis, based on the meta-data of the standard "vn". - This is done
  !  for all available time levels.
  !
  SUBROUTINE init_vn_vertical(jg, job_type, prefix, l_init_prm_diag, nlev, dst_axis, dst_varlist)
    INTEGER,          INTENT(IN)  :: jg              !< domain number
    INTEGER,          INTENT(IN)  :: job_type        !< vertical interpolation type
    CHARACTER(LEN=*), INTENT(IN)  :: prefix          !< job name prefix
    LOGICAL,          INTENT(IN)  :: l_init_prm_diag !< Flag. If .TRUE., then prm_diag data structure is available
    INTEGER,          INTENT(IN)  :: nlev            !< number of vertical levels
    INTEGER,          INTENT(IN)  :: dst_axis        !< destination axis
    TYPE(t_var_list), POINTER     :: dst_varlist     !< destination variable list
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_pp_scheduler:init_vn_vertical")
    TYPE(t_list_element), POINTER :: element_u, element_v, element, vn_element, new_element, new_element_2
    INTEGER                       :: i, shape3d_c(3), shape3d_e(3), nblks_c, nblks_e, tl
    TYPE(t_job_queue),    POINTER :: task
    TYPE(t_var_metadata), POINTER :: info
    REAL(wp),             POINTER :: p_opt_field_r3d(:,:,:)
    CHARACTER(len=varname_len)    :: name
    CHARACTER(len=4)              :: suffix
    TYPE(t_cf_var)                :: cf
    TYPE(t_grib2_var)             :: grib2
    TYPE(t_post_op_meta)          :: post_op

    if (dbg_level > 5)  CALL message(routine, "Enter")

    ! initialize "new_element" pointer (cf. NEC compiler bugs DWD0121
    ! and DWD0123 for hybrid parallelization)
    new_element   => NULL()
    new_element_2 => NULL()
     
    !- find existing variables "u", "v" (for copying the meta-data):
    element_u => find_list_element (p_nh_state(jg)%diag_list, "u")
    element_v => find_list_element (p_nh_state(jg)%diag_list, "v")

    !- predefined array shapes
    nblks_c   = p_patch(jg)%nblks_c
    nblks_e   = p_patch(jg)%nblks_e
    shape3d_c = (/ nproma, nlev, nblks_c /)
    shape3d_e = (/ nproma, nlev, nblks_e /)

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

        ! Check for matching name (take care of suffix of
        ! time-dependent variables):
        IF (TRIM(vn_name) /= TRIM(tolower(get_var_name(element%field)))) CYCLE

        ! get time level
        tl = get_var_timelevel(element%field)
        suffix = ""
        IF (tl /= -1)  WRITE (suffix,'(".TL",i1)') tl

        !-- create a new z/p/i-variable "vn":
        CALL add_var( dst_varlist, TRIM(info%name), p_opt_field_r3d, element%field%info%hgrid,    &
          &           dst_axis, info%cf, info%grib2, ldims=shape3d_e,                             &
          &           vert_interp=info%vert_interp, new_element=vn_element,                       &
          &           post_op=info%post_op, lrestart=.FALSE. )
         
        !-- create a post-processing task for vertical interpolation of "vn"
        task => pp_task_insert(DEFAULT_PRIORITY1)
        task%job_name        =  &
          &  TRIM(prefix)//" interp. "//TRIM(info%name)//", DOM "//TRIM(int2string(jg))
        IF (dbg_level > 8) CALL message(routine, task%job_name)

        task%job_type                    =  job_type
        task%activity                    =  new_simulation_status(l_output_step=.TRUE.)
        task%activity%ldom_active(jg)    =  .TRUE.
        task%activity%i_timelevel(jg)    =  get_var_timelevel(element%field)
        task%data_input%jg               =  jg 
        task%data_input%p_patch          => p_patch(jg)          
        task%data_input%p_int_state      => p_int_state(jg)
        task%data_input%p_nh_state       => p_nh_state(jg)
        task%data_input%nh_pzlev_config  => nh_pzlev_config(jg)
        task%data_input%p_nh_opt_diag    => p_nh_opt_diag(jg)
        IF (l_init_prm_diag) THEN
          task%data_input%prm_diag       => prm_diag(jg)
        ELSE
          task%data_input%prm_diag       => NULL() 
        END IF
        task%data_input%var  => element%field      ! set input variable
        task%data_output%var => vn_element%field   ! set output variable

        !-- create new cell-based variables "u", "v" on the same vertical axis

        name    = TRIM(get_var_name(element_u%field))//suffix
        cf      = element_u%field%info%cf
        grib2   = element_u%field%info%grib2
        post_op = element_u%field%info%post_op
        CALL add_var( dst_varlist, TRIM(name), p_opt_field_r3d,                           &
          & GRID_UNSTRUCTURED_CELL, dst_axis, cf, grib2,                                  &
          & ldims=shape3d_c, lrestart=.FALSE., in_group=element_u%field%info%in_group,    &
          & new_element=new_element, post_op=post_op )

        name    = TRIM(get_var_name(element_v%field))//suffix
        cf      = element_v%field%info%cf
        grib2   = element_v%field%info%grib2
        post_op = element_v%field%info%post_op
        CALL add_var( dst_varlist, TRIM(name), p_opt_field_r3d,                           &
          & GRID_UNSTRUCTURED_CELL, dst_axis, cf, grib2,                                  &
          & ldims=shape3d_c, lrestart=.FALSE., in_group=element_v%field%info%in_group,    &
          & new_element=new_element_2, post_op=post_op )

        !-- create a post-processing task for edge2cell interpolation "vn" -> "u","v"
        task => pp_task_insert(DEFAULT_PRIORITY2)
        WRITE (task%job_name, *) "edge2cell interp. ",TRIM(info%name),", DOM ",jg
        task%data_input%p_nh_state      => NULL()
        task%data_input%prm_diag        => NULL()
        task%data_input%nh_pzlev_config => NULL()
        task%data_input%jg              =  jg           
        task%data_input%p_patch         => p_patch(jg)
        task%data_input%p_nh_opt_diag   => p_nh_opt_diag(jg)
        task%data_input%p_int_state     => p_int_state(jg)
        task%job_type                   =  TASK_INTP_EDGE2CELL
        task%activity                   =  new_simulation_status(l_output_step=.TRUE.)
        task%activity%ldom_active(jg)   =  .TRUE.
        task%activity%i_timelevel(jg)   =  get_var_timelevel(element%field)
        task%data_input%var             => vn_element%field    ! set input variable
        task%data_output%var            => new_element%field   ! set output variable
        task%data_output%var_2          => new_element_2%field ! set Y-component
      END DO
    END DO

    if (dbg_level > 5)  CALL message(routine, "Done")

  END SUBROUTINE init_vn_vertical


  !---------------------------------------------------------------
  !> Setup of i/p/z-level interpolation tasks.
  !  - collects lists of variables for i/p/z interpolation
  !  - adds variable fields where interpolation results will be stored
  !  - creates "post-processing tasks", i.e. entries in a list 
  !    which is regularly traversed during the model run.
  !
  ! See SUBROUTINE pp_scheduler_init for further details.
  !
  SUBROUTINE pp_scheduler_init_ipz(l_init_prm_diag)
    LOGICAL, INTENT(IN) :: l_init_prm_diag !< Flag. If .TRUE., then prm_diag data structure is available

    ! local variables
    CHARACTER(*), PARAMETER :: routine =  &
      &  TRIM("mo_pp_scheduler:pp_scheduler_init_ipz")
    INTEGER,      PARAMETER :: init_tasks(3) = &
      &  (/ TASK_INIT_VER_Z, TASK_INIT_VER_P, TASK_INIT_VER_I /)
    CHARACTER,    PARAMETER :: init_names(3) = &
      &  (/ 'z', 'p', 'i' /)
    INTEGER                            :: &
      &  jg, ndom, ibits, nblks_c, nblks_v, ierrstat, ivar, i,      &
      &  iaxis, vgrid, nlev, nvars_pl, nvars_hl, nvars_il, nvars,   &
      &  job_type, z_id, p_id, i_id, shape3d(3)
    LOGICAL                            :: &
      &  l_intp_p, l_intp_z, l_intp_i, found, &
      &  l_uv_vertical_intp_z, l_uv_vertical_intp_p, l_uv_vertical_intp_i, &
      &  l_uv_vertical_intp
    TYPE(t_job_queue),         POINTER :: task
    TYPE(t_nh_diag_pz),        POINTER :: p_diag_pz
    TYPE(t_var_list),          POINTER :: p_opt_diag_list_p, p_opt_diag_list_z, &
      &                                   p_opt_diag_list_i, p_opt_diag_list
    REAL(wp),                  POINTER :: p_opt_field_r3d(:,:,:)
    TYPE(t_list_element),      POINTER :: element, new_element
    ! variable lists (for all domains + output name lists):
    CHARACTER(LEN=vname_len), TARGET, ALLOCATABLE  :: &
         &                                pl_varlist(:), hl_varlist(:), il_varlist(:)
    CHARACTER(LEN=vname_len),  POINTER :: varlist(:)
    CHARACTER(LEN=10)                  :: prefix
    TYPE(t_var_metadata),      POINTER :: info
    TYPE(t_cf_var)                     :: cf_desc
    TYPE(t_grib2_var)                  :: grib2_desc

    ! initialize "new_element" pointer (cf. NEC compiler bugs DWD0121
    ! and DWD0123 for hybrid parallelization)
    new_element   => NULL()

    if (dbg_level > 5)  CALL message(routine, "Enter")
    ndom = SIZE(p_nh_opt_diag)
   
    ! loop over the domains and output definitions
    ! - check for which domain p- or z-level interpolation has been
    !   requested
    ! - register setup routine for pz-level interpolation (must be
    !   executed ahead of interpolation tasks)
    ! - for each variable: register post-processing task

    ALLOCATE(pl_varlist(ndom*max_var_pl), hl_varlist(ndom*max_var_hl), &
         &   il_varlist(ndom*max_var_il), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! list indices for the vertical interpolation types:
!CDIR NOIEXPAND
    z_id = vintp_type_id("Z")
!CDIR NOIEXPAND
    p_id = vintp_type_id("P")
!CDIR NOIEXPAND
    i_id = vintp_type_id("I")

    DOM_LOOP : DO jg=1,ndom
      IF (dbg_level > 8)  CALL message(routine, "DOM "//int2string(jg))

      !-- check if any output name list requests p- or z- or i-level 
      !-- interpolation for this domain, collect the list of variables

      ! loop in search of pressure-level interpolation      
      CALL collect_output_variables(jg, "P", max_var_pl, &                                ! in
        &                           nvars_pl, l_intp_p, pl_varlist, l_uv_vertical_intp_p) ! out

      ! loop in search of height-level interpolation
      CALL collect_output_variables(jg, "Z", max_var_hl, &                                ! in
        &                           nvars_hl, l_intp_z, hl_varlist, l_uv_vertical_intp_z) ! out

      ! loop in search of i-level interpolation
      CALL collect_output_variables(jg, "I", max_var_il, &                                ! in
        &                           nvars_il, l_intp_i, il_varlist, l_uv_vertical_intp_i) ! out

      ! now, we have total variables lists "hl_varlist(1:nvars_hl)"
      ! and "pl_varlist(1:nvars_pl)" and "il_varlist(1:nvars_il)"

      ! some debugging output...
      IF (dbg_level > 8)  THEN
        DO i=1,nvars_pl
          WRITE (message_text,*) "p var list: ", TRIM(pl_varlist(i))
          CALL message(routine, message_text)
        END DO
        DO i=1,nvars_hl
          WRITE (message_text,*) "h var list: ", TRIM(hl_varlist(i))
          CALL message(routine, message_text)
        END DO
        DO i=1,nvars_il
          WRITE (message_text,*) "i var list: ", TRIM(il_varlist(i))
          CALL message(routine, message_text)
        END DO
      END IF
      
      ! skip domain if no p/z-interpolation requested:
      IF (.NOT. (l_intp_z .OR. l_intp_p .OR. l_intp_i)) CYCLE DOM_LOOP

      ! remove duplicates from variable lists
      IF (l_intp_z) CALL remove_duplicates(hl_varlist, nvars_hl)
      IF (l_intp_p) CALL remove_duplicates(pl_varlist, nvars_pl)
      IF (l_intp_i) CALL remove_duplicates(il_varlist, nvars_il)

      !-- First, add some diagnostic variables which are essential for
      !-- p/z-level interpolation, e.g. temp_z, pres_z:
      p_diag_pz         => p_nh_opt_diag(jg)%diag_pz
      p_opt_diag_list_z => p_nh_opt_diag(jg)%opt_diag_list_z
      p_opt_diag_list_p => p_nh_opt_diag(jg)%opt_diag_list_p
      p_opt_diag_list_i => p_nh_opt_diag(jg)%opt_diag_list_i
      ibits     = DATATYPE_PACK16   ! "entropy" of horizontal slice

      ! predefined array shapes
      nblks_c   = p_patch(jg)%nblks_c
      nblks_v   = p_patch(jg)%nblks_v

      ! add new variable fields for the z/p/i-axis, based on the
      ! meta-data of an existing variable field (which is defined on
      ! model/half levels):
      IF (l_intp_z) THEN
        shape3d = (/ nproma, nh_pzlev_config(jg)%nzlev, nblks_c /)
        CALL copy_variable("temp", p_nh_state(jg)%diag_list, ZA_ALTITUDE, shape3d, &
          &                p_diag_pz%z_temp, p_opt_diag_list_z)
        CALL copy_variable("pres", p_nh_state(jg)%diag_list, ZA_ALTITUDE, shape3d, &
          &                p_diag_pz%z_pres, p_opt_diag_list_z)
      END IF
      IF (l_intp_p) THEN
        shape3d = (/ nproma, nh_pzlev_config(jg)%nplev, nblks_c /)
        cf_desc    = t_cf_var('gh', 'm', 'geopotential height', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 3, 5, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_opt_diag_list_p, 'gh', p_diag_pz%p_gh,                  &
          & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE, cf_desc, grib2_desc,           &
          & ldims=shape3d, lrestart=.FALSE. )
        CALL copy_variable("temp",   p_nh_state(jg)%diag_list,    ZA_PRESSURE, shape3d, &
          &                p_diag_pz%p_temp, p_opt_diag_list_p)
      END IF
      IF (l_intp_i) THEN
        shape3d = (/ nproma, nh_pzlev_config(jg)%nilev, nblks_c /)
        cf_desc    = t_cf_var('gh', 'm', 'geopotential height', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 3, 5, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_opt_diag_list_i, 'gh', p_diag_pz%i_gh,                  &
          & GRID_UNSTRUCTURED_CELL, ZA_ISENTROPIC, cf_desc, grib2_desc,         &
          & ldims=shape3d, lrestart=.FALSE. )
        CALL copy_variable("temp",   p_nh_state(jg)%diag_list,    ZA_ISENTROPIC, shape3d, &
          &                p_diag_pz%i_temp, p_opt_diag_list_i)
      END IF

      !-- register interpolation setup as post-processing task(s)
      DO i=1,SIZE(init_tasks)
        IF ((init_tasks(i) == TASK_INIT_VER_Z) .AND. .NOT. l_intp_z) CYCLE
        IF ((init_tasks(i) == TASK_INIT_VER_P) .AND. .NOT. l_intp_p) CYCLE
        IF ((init_tasks(i) == TASK_INIT_VER_I) .AND. .NOT. l_intp_i) CYCLE

        task => pp_task_insert(HIGH_PRIORITY)
        task%data_input%p_int_state      => p_int_state(jg)
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
        task%activity%ldom_active(jg)    = .TRUE.
        task%activity%i_timelevel(jg)    = ALL_TIMELEVELS
        task%job_type                    = init_tasks(i)
        task%job_name                    = "Init: "//init_names(i)//"-level interpolation, DOM "//TRIM(int2string(jg))
      END DO

      !-- register clean-up routine as a post-processing task
      task => pp_task_insert(LOW_PRIORITY)
      task%activity     = new_simulation_status(l_output_step=.TRUE.)
      task%activity%ldom_active(jg)    =  .TRUE.
      task%activity%i_timelevel(jg)    =  ALL_TIMELEVELS
      task%data_input%p_nh_opt_diag => p_nh_opt_diag(jg)
      task%job_name     = "Clean-up: ipz-level interpolation, level "//TRIM(int2string(jg))
      task%job_type     = TASK_FINALIZE_IPZ
      task%data_input%p_int_state      => NULL()
      task%data_input%jg               =  jg           
      task%data_input%p_patch          => p_patch(jg)
      task%data_input%p_nh_state       => p_nh_state(jg)
      task%data_input%prm_diag         => NULL() 
      task%data_input%nh_pzlev_config  => nh_pzlev_config(jg)

      ! remove already defined variables from list of requested output
      ! fields:
      IF (l_intp_z) CALL difference(hl_varlist, nvars_hl, &
        &                           (/ "temp  ", "pres  ", "u     ", "v     " /), 4)
      IF (l_intp_p) CALL difference(pl_varlist, nvars_pl, &
        &                           (/ "gh    ", "temp  ", "u     ", "v     " /), 4)
      IF (l_intp_i) CALL difference(il_varlist, nvars_il, &
        &                           (/ "gh    ", "temp  ", "u     ", "v     " /), 4)

      !-- loop over requested p-, z-, and i-level variables, add variables
      !-- ("add_var") and register interpolation tasks:
      DO iaxis=1,3
        IF (iaxis == 1) THEN
          prefix  =  "z-level"
          varlist => hl_varlist
          nvars   =  nvars_hl
          nlev    =  nh_pzlev_config(jg)%nzlev
          vgrid   =  ZA_ALTITUDE
          p_opt_diag_list => p_opt_diag_list_z
          job_type = TASK_INTP_VER_ZLEV
          l_uv_vertical_intp = l_uv_vertical_intp_z
        END IF
        IF (iaxis == 2) THEN
          prefix  =  "p-level"
          varlist => pl_varlist
          nvars   =  nvars_pl
          nlev    =  nh_pzlev_config(jg)%nplev
          vgrid   =  ZA_PRESSURE
          p_opt_diag_list => p_opt_diag_list_p
          job_type = TASK_INTP_VER_PLEV
          l_uv_vertical_intp = l_uv_vertical_intp_p
        END IF
        IF (iaxis == 3) THEN
          prefix  =  "i-level"
          varlist => il_varlist
          nvars   =  nvars_il
          nlev    =  nh_pzlev_config(jg)%nilev
          vgrid   =  ZA_ISENTROPIC
          p_opt_diag_list => p_opt_diag_list_i
          job_type = TASK_INTP_VER_ILEV
          l_uv_vertical_intp = l_uv_vertical_intp_i
        END IF

        !-- if "u", "v" appear in the variable list...
        IF (l_uv_vertical_intp) THEN
          ! for each time level, create a new z/p/i-variable "vn",
          ! create new cell-based variables "u", "v" on the same
          ! vertical axis, create a post-processing task for vertical
          ! interpolation of "vn", create a post-processing task for
          ! edge2cell interpolation "vn" -> "u","v":
          CALL init_vn_vertical(jg, job_type, prefix, l_init_prm_diag, &
            &                   nlev, vgrid, p_opt_diag_list)
        END IF

        DO ivar=1,nvars
          IF (dbg_level > 8) &
            CALL message(routine, TRIM(prefix)//": Looking for input var '"//TRIM(varlist(ivar))//"'")
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
                IF (.NOT. info%vert_interp%vert_intp_type(z_id)) CYCLE
              END IF
              IF (iaxis == 2) THEN
                IF (.NOT. info%vert_interp%vert_intp_type(p_id)) CYCLE
              END IF
              IF (iaxis == 3) THEN
                IF (.NOT. info%vert_interp%vert_intp_type(i_id)) CYCLE
              END IF

              ! Check for matching name (take care of suffix of
              ! time-dependent variables):
              IF (TRIM(varlist(ivar)) /= TRIM(tolower(get_var_name(element%field)))) CYCLE

              ! throw error message, if this variable is not a REAL field:
              IF (.NOT. ASSOCIATED(element%field%r_ptr)) THEN
                CALL finish(routine, TRIM(info%name)//": i/p/z interpolation implemented for REAL fields only.")
              END IF

              ! Found it, add it to the variable list of optional
              ! diagnostics
              IF ( (info%used_dimensions(1) /= nproma)  .OR.   &
                &  ((info%used_dimensions(3) /= nblks_c) .AND. &
                &   (info%used_dimensions(3) /= nblks_v)) ) THEN
                CALL finish(routine, "Unexpected field size!")
              END IF
              ! Note: Even vertex-based variables are interpolated
              ! onto a cell-based variable, since we interpolate the
              ! vertex-based vars to cell-based vars first:
              shape3d  = (/ info%used_dimensions(1), nlev, nblks_c /)

              CALL add_var( p_opt_diag_list, info%name, p_opt_field_r3d, &
                &           info%hgrid, vgrid, info%cf, info%grib2,      &
                &           ldims=shape3d, lrestart=.FALSE.,             &
                &           isteptype=info%isteptype,                    &
                &           loutput=.TRUE., new_element=new_element,     &
                &           post_op=info%post_op)

              !-- add post-processing task for interpolation

              task => pp_task_insert(DEFAULT_PRIORITY1)
              task%job_name        =  &
                &  TRIM(prefix)//" interp. "//TRIM(info%name)  &
                &  //", DOM "//TRIM(int2string(jg))
              IF (dbg_level > 8) CALL message(routine, task%job_name)

              task%job_type                    =  job_type
              task%activity                    =  new_simulation_status(l_output_step=.TRUE.)
              task%activity%ldom_active(jg)    =  .TRUE.
              task%activity%i_timelevel(jg)    =  get_var_timelevel(element%field)
              task%data_input%jg               =  jg 
              task%data_input%p_patch          => p_patch(jg)          
              task%data_input%p_int_state      => p_int_state(jg)
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
      END DO ! height/pressure/isentropic axis
    END DO DOM_LOOP  ! jg
    
    DEALLOCATE(pl_varlist, hl_varlist, il_varlist,STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    IF (dbg_level > 5)  CALL message(routine, "Done")
    
  END SUBROUTINE pp_scheduler_init_ipz


  !---------------------------------------------------------------
  !> Register a new post-processing task for computing additional 
  !  diagnostic fields (not for interpolation).
  !
  !  The new task will be added to the dynamic list of post-processing
  !  jobs and called on a regular basis.
  !
  SUBROUTINE pp_scheduler_register(name, jg, p_out_var,             &
    &                              l_init_prm_diag, job_type,       &
    &                              opt_priority, opt_l_output_step, &
    &                              opt_p_in_var)

    CHARACTER(LEN=*)                   , INTENT(IN) :: name
    INTEGER                            , INTENT(IN) :: jg
    TYPE (t_list_element), POINTER                  :: p_out_var
    LOGICAL                            , INTENT(IN) :: l_init_prm_diag
    INTEGER                            , INTENT(IN) :: job_type
    INTEGER, OPTIONAL                  , INTENT(IN) :: opt_priority
    LOGICAL, OPTIONAL                  , INTENT(IN) :: opt_l_output_step
    TYPE (t_list_element), POINTER, OPTIONAL        :: opt_p_in_var
    ! local variables
    LOGICAL                    :: l_output_step
    INTEGER                    :: priority
    TYPE(t_job_queue), POINTER :: task
    
    ! set default values
    l_output_step = .TRUE.
    priority      = DEFAULT_PRIORITY0
    ! assign optional parameters:
    CALL assign_if_present(priority, opt_priority)
    CALL assign_if_present(l_output_step, opt_l_output_step)
    ! create a post-processing task and fill its input/output data:
    task => pp_task_insert(priority)
    WRITE (task%job_name, *) TRIM(name),", DOM ",jg
    IF (PRESENT(opt_p_in_var)) THEN
      task%data_input%var            => opt_p_in_var%field
    END IF
    task%data_input%jg               =  jg           
    task%data_input%p_nh_state       => p_nh_state(jg)
    task%data_input%p_patch          => p_patch(jg)
    task%data_input%nh_pzlev_config  => nh_pzlev_config(jg)
    IF (l_init_prm_diag) THEN
      task%data_input%prm_diag       => prm_diag(jg)
    END IF
    task%data_output%var             => p_out_var%field
    task%job_type                    =  job_type
    task%activity                    =  new_simulation_status(l_output_step=l_output_step)
    task%activity%ldom_active(jg)    =  .TRUE.
    task%activity%i_timelevel(jg)    =  ALL_TIMELEVELS

  END SUBROUTINE pp_scheduler_register


  !---------------------------------------------------------------
  !> Loop over job queue, call active tasks.
  !
  SUBROUTINE pp_scheduler_process(simulation_status)
    TYPE(t_simulation_status), INTENT(IN) :: simulation_status
    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_pp_scheduler:pp_scheduler_process")
    TYPE(t_job_queue), POINTER :: ptr_task

    IF (dbg_level >= 10) THEN
      CALL message(routine,"Processing task list...")
    END IF
    ptr_task => job_queue
    ! loop over job queue
    LOOP_JOB : DO
      IF (.NOT. ASSOCIATED(ptr_task)) EXIT
      IF (.NOT. pp_task_is_active(ptr_task, simulation_status)) THEN
        IF (dbg_level > 5) THEN
          WRITE(message_text,*) "Skipping task '", TRIM(ptr_task%job_name), "'"
          CALL message(routine, TRIM(message_text))
        END IF
        ptr_task => ptr_task%next
        CYCLE LOOP_JOB
      END IF

      IF (dbg_level > 5) THEN
        WRITE(message_text,*) "Staging task '", TRIM(ptr_task%job_name), "'"
        CALL message(routine, TRIM(message_text))
      END IF

      SELECT CASE ( ptr_task%job_type )

        ! initialize vertical interpolation:
      CASE ( TASK_INIT_VER_Z, TASK_INIT_VER_P, TASK_INIT_VER_I, &
           & TASK_FINALIZE_IPZ)
        CALL pp_task_ipzlev_setup(ptr_task)

        ! perform horizontal interpolation:
      CASE ( TASK_INTP_HOR_LONLAT )
        CALL pp_task_lonlat(ptr_task)

        ! perform vertical interpolation:
      CASE ( TASK_INTP_VER_PLEV, TASK_INTP_VER_ZLEV, TASK_INTP_VER_ILEV )
        CALL pp_task_ipzlev(ptr_task)

        ! synchronize halo regions:
      CASE ( TASK_INTP_SYNC )
        CALL pp_task_sync(simulation_status)

        ! compute mean sea level pressure:
      CASE ( TASK_INTP_MSL )
        CALL pp_task_intp_msl(ptr_task)

        ! compute relative humidty, vertical velocity
      CASE ( TASK_COMPUTE_RH, TASK_COMPUTE_OMEGA )
        CALL pp_task_compute_field(ptr_task)

        ! vector reconstruction on cell centers:
      CASE ( TASK_INTP_EDGE2CELL )
        CALL pp_task_edge2cell(ptr_task)

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


  !--- UTILITY ROUTINES --------------------------------------------------------------

  !---------------------------------------------------------------
  !> @return .TRUE. if given post-processing task is in active state.
  ! 
  ! Tasks may be inactive, e.g. outside the output intervals.
  FUNCTION pp_task_is_active(ptr_task, sim_status)
    LOGICAL :: pp_task_is_active
    TYPE(t_job_queue), POINTER :: ptr_task
    TYPE(t_simulation_status),  INTENT(IN) :: sim_status

    INTEGER :: i

    ! compare simulation status to post-processing tasks activity
    ! flags, then check if any of the activity conditions is
    ! fulfilled:
    pp_task_is_active = .FALSE.
    IF (ANY(ptr_task%activity%status_flags(:)  .AND.  &
      &     sim_status%status_flags(:))) pp_task_is_active = .TRUE.

    ! check, if current task applies only to domains which are
    ! "active":
    DO i=1,n_dom
      IF (ptr_task%activity%ldom_active(i) .AND. &
        & .NOT. sim_status%ldom_active(i)) THEN
        pp_task_is_active = .FALSE.
      END IF
    END DO

    ! check, if current task matches the variable time level (TL1,
    ! TL2, ...) of the simulation status:
    DO i=1,n_dom
      IF (ptr_task%activity%ldom_active(i)) THEN
        IF  ((ptr_task%activity%i_timelevel(i) /= ALL_TIMELEVELS) .AND.  &
          &  (ptr_task%activity%i_timelevel(i) /= sim_status%i_timelevel(i))) THEN
          pp_task_is_active = .FALSE.
        END IF
      END IF
    END DO

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

    IF (dbg_level > 5)  CALL message(routine, "Inserting pp task")
    
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
    element%job_priority            =  job_priority
    element%next                    => tmp
    element%activity%i_timelevel(:) =  UNDEF_TIMELEVEL
  END FUNCTION pp_task_insert


  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for "t_simulation_status" variables
  ! 
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION new_simulation_status(l_output_step, l_first_step, l_last_step, &
    &                            l_dom_active, i_timelevel)  &
    RESULT(sim_status)

    TYPE(t_simulation_status) :: sim_status
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  l_output_step, l_first_step, l_last_step
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  l_dom_active(:)
    INTEGER, INTENT(IN), OPTIONAL      :: &
      &  i_timelevel(:)
    ! local variables
    INTEGER :: ndom

    ! set default values
    sim_status%status_flags(:) = (/ .FALSE., .FALSE., .FALSE. /)

    ! supersede with user definitions
    CALL assign_if_present(sim_status%status_flags(1), l_output_step)
    CALL assign_if_present(sim_status%status_flags(2), l_first_step)
    CALL assign_if_present(sim_status%status_flags(3), l_last_step)

    ! as a default, all domains are "inactive", i.e. the activity
    ! flags are not considered:
    sim_status%ldom_active(:) = .FALSE.
    IF  (PRESENT(l_dom_active)) THEN
      ndom = SIZE(l_dom_active)
      sim_status%ldom_active(1:ndom) = l_dom_active(1:ndom)
    END IF

    ! as a default, no special timelevel is set for (in-)activity:
    sim_status%i_timelevel(:) = ALL_TIMELEVELS
    CALL assign_if_present(sim_status%i_timelevel, i_timelevel)

  END FUNCTION new_simulation_status

END MODULE mo_pp_scheduler
