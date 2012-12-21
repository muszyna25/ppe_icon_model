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
  USE mo_impl_constants,          ONLY: SUCCESS,                       &
    & HINTP_TYPE_NONE, HINTP_TYPE_LONLAT,                              &
    & max_dom, max_var_ml, max_var_pl, max_var_hl, max_var_il,         &
    & TASK_NONE, TASK_INIT_VER_Z, TASK_INIT_VER_P, TASK_INIT_VER_I,    &
    & TASK_FINALIZE_IPZ,                                               &
    & TASK_INTP_HOR_LONLAT, TASK_INTP_VER_PLEV, TASK_INTP_SYNC,        &
    & TASK_INTP_MSL, TASK_COMPUTE_RH, TASK_INTP_VER_ZLEV,              &
    & TASK_INTP_VER_ILEV, TASK_INTP_EDGE2CELL, max_phys_dom
  USE mo_model_domain,            ONLY: t_patch, p_patch, p_phys_patch
  USE mo_var_list,                ONLY: add_var, get_all_var_names,         &
    &                                   create_hor_interp_metadata,         &
    &                                   nvar_lists, var_lists, vintp_types, &
    &                                   vintp_type_id
  USE mo_var_list_element,        ONLY: t_var_list_element, level_type_ml,  &
    &                                   level_type_pl, level_type_hl, level_type_il
  USE mo_var_metadata,            ONLY: t_var_metadata, t_vert_interp_meta, &
    &                                   VINTP_TYPE_LIST
  USE mo_intp_data_strc,          ONLY: t_int_state, lonlat_grid_list, &
    &                                   t_lon_lat_intp, p_int_state
  USE mo_nonhydro_types,          ONLY: t_nh_state, t_nh_prog, t_nh_diag, &
    &                                   t_nh_metrics
  USE mo_nonhydro_state,          ONLY: p_nh_state
  USE mo_opt_diagnostics,         ONLY: t_nh_diag_pz, t_nh_opt_diag, p_nh_opt_diag
  USE mo_nwp_phy_types,           ONLY: t_nwp_phy_diag
  USE mo_nwp_phy_state,           ONLY: prm_diag
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
    &                                   GRID_UNSTRUCTURED_CELL, ZA_ALTITUDE,     &
    &                                   ZA_PRESSURE, GRID_REGULAR_LONLAT,        &
    &                                   GRID_UNSTRUCTURED_EDGE,                  &
    &                                   GRID_UNSTRUCTURED_VERT,                  &
    &                                   DATATYPE_FLT32, DATATYPE_PACK16, ZA_ISENTROPIC, &
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
    &                                   DEFAULT_PRIORITY2, LOW_PRIORITY, dbg_level
  USE mo_fortran_tools,           ONLY: assign_if_present

  IMPLICIT NONE

  ! interface definition
  PRIVATE

  ! functions and subroutines
  PUBLIC :: pp_scheduler_init
  PUBLIC :: pp_scheduler_register
  PUBLIC :: pp_scheduler_process
  PUBLIC :: pp_scheduler_finalize
  PUBLIC :: new_simulation_status


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
          CASE (TASK_COMPUTE_RH) ! relative humidity
            CALL pp_scheduler_register( name=element%field%info%name, jg=jg, p_out_var=element, &
              &                         l_init_prm_diag=l_init_prm_diag, job_type=TASK_COMPUTE_RH )
          CASE (TASK_INTP_MSL)   ! mean sea level pressure
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
  !> Setup of lon-lat interpolation tasks.
  !
  SUBROUTINE pp_scheduler_init_lonlat

    ! local variables
    CHARACTER(*), PARAMETER :: routine =  &
      &  TRIM("mo_pp_scheduler:pp_scheduler_init_lonlat")
    INTEGER                               :: &
      &  jg, ndom, ierrstat, ivar, i, j, idx, nvars_ll, &
      &  nblks_lonlat, ilev_type, max_var, ilev
    LOGICAL                               :: &
      &  found, l_horintp
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
        IF (.NOT. lonlat_grid_list(ll_vargrid(ivar))%l_dom(jg)) &
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
          CASE (level_type_il)
            p_opt_diag_list => p_nh_opt_diag(jg)%opt_diag_list_i
          END SELECT

          ! set local values for "nblks" and "npromz"
          ptr_int_lonlat => lonlat_grid_list(ll_vargrid(ivar))%intp(jg)
          nblks_lonlat   =  (ptr_int_lonlat%nthis_local_pts - 1)/nproma + 1
          var_shape         =  info%used_dimensions(:)
          IF (is_2d_field(info%vgrid) .AND. (info%ndims /= 2)) THEN
            CALL finish(routine, "Inconsistent dimension info!")
          END IF
          IF (is_2d_field(info%vgrid)) THEN
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
              &           isteptype=info%isteptype,                             &
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
              &           isteptype=info%isteptype,                             &
              &           hor_interp=create_hor_interp_metadata(                &
              &             hor_intp_type=HINTP_TYPE_NONE ))
            CALL add_var( p_opt_diag_list, TRIM(vname)//".Y", p_opt_field_r3d,  &
              &           GRID_REGULAR_LONLAT, info%vgrid, info%cf, info%grib2, &
              &           ldims=var_shape, lrestart=.FALSE.,                    &
              &           loutput=.TRUE., new_element=new_element_2,            &
              &           isteptype=info%isteptype,                             &
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
          task%activity%ldom_active(jg) =  .TRUE.
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

!      ! Check that at least one element with this name has been found
!      IF(.NOT. found) &
!        CALL finish(routine,'No feasible variable found: ' &
!        &   //TRIM(ll_varlist(ivar)))
    END DO ! ivar

    DEALLOCATE(ll_varlist, ll_vargrid, ll_varlevs, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    ! If at least one interpolation task has been created, we add a
    ! setup task which synchronizes the halo regions:
    IF (l_horintp) THEN
      IF (dbg_level >= 10) &
        CALL message(routine, "Creating synchronization task for horizontal interpolation.")
      task => pp_task_insert(DEFAULT_PRIORITY1)
      WRITE (task%job_name, *) "horizontal interp. SYNC"
      task%job_type = TASK_INTP_SYNC
      task%activity = new_simulation_status(l_output_step=.TRUE.)
      task%activity%ldom_active(:) = .TRUE.
    END IF

    IF (dbg_level > 5)  CALL message(routine, "Done")
    
  END SUBROUTINE pp_scheduler_init_lonlat


  !---------------------------------------------------------------
  !> Setup of i/p/z-level interpolation tasks.
  !
  ! See SUBROUTINE pp_scheduler_init for further details.
  !
  SUBROUTINE pp_scheduler_init_ipz(l_init_prm_diag)
    LOGICAL, INTENT(IN) :: l_init_prm_diag

    ! local variables
    CHARACTER(*), PARAMETER :: routine =  &
      &  TRIM("mo_pp_scheduler:pp_scheduler_init_ipz")
    INTEGER,      PARAMETER :: init_tasks(3) = &
      &  (/ TASK_INIT_VER_Z, TASK_INIT_VER_P, TASK_INIT_VER_I /)
    CHARACTER,    PARAMETER :: init_names(3) = &
      &  (/ 'z', 'p', 'i' /)
    INTEGER                            :: &
      &  jg, ndom, ibits, nblks_c, nblks_v, ierrstat, ivar, i, idx, &
      &  iaxis, vgrid, nlev, iphys_dom
    LOGICAL                            :: &
      &  l_jg_active, l_intp_p, l_intp_z, l_intp_i, found
    TYPE (t_output_name_list), POINTER :: p_onl
    TYPE(t_job_queue),         POINTER :: task
    INTEGER                            :: shape3d(3)
    TYPE(t_cf_var)                     :: cf_desc
    TYPE(t_grib2_var)                  :: grib2_desc
    TYPE(t_nh_diag_pz),        POINTER :: p_diag_pz
    TYPE(t_var_list),          POINTER :: p_opt_diag_list_p, p_opt_diag_list_z, &
      &                                   p_opt_diag_list_i, p_opt_diag_list
    REAL(wp),                  POINTER :: p_opt_field_r3d(:,:,:)
    TYPE(t_list_element),      POINTER :: element, new_element
    ! variable lists (for all domains + output name lists):
    INTEGER                            :: nvars_pl, nvars_hl, nvars_il, nvars, &
         &                                job_type
    CHARACTER(LEN=vname_len), TARGET, ALLOCATABLE  :: &
         &                                pl_varlist(:), hl_varlist(:), il_varlist(:)
    CHARACTER(LEN=vname_len),  POINTER :: varlist(:)
    CHARACTER(LEN=10)                  :: prefix
    TYPE(t_var_metadata),      POINTER :: info
    INTEGER :: z_id, p_id, i_id
    LOGICAL :: vert_intp_type(SIZE(VINTP_TYPE_LIST))

    if (dbg_level > 5)  CALL message(routine, "Enter")
    ndom = SIZE(p_nh_opt_diag)
   
    ! loop over the domains and output definitions
    ! - check for which domain p- or z-level interpolation has been
    !   requested
    ! - register setup routine for pz-level interpolation (must be
    !   executed ahead of interpolation tasks)
    ! - for each variable: register post-processing task

    ALLOCATE(pl_varlist(ndom*max_var_pl), &
         &   hl_varlist(ndom*max_var_hl), &
         &   il_varlist(ndom*max_var_il), &
         &   STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

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
      p_onl => first_output_name_list
      nvars_pl = 0
      l_intp_p = .FALSE.
      NML_LOOP_P : DO
        IF (.NOT.ASSOCIATED(p_onl)) EXIT NML_LOOP_P

        IF (dbg_level >= 15)  WRITE (0,*) p_onl%pl_varlist 

        ! If dom(:) was not specified in namelist input, it is set
        ! completely to -1.  In this case all domains are wanted in
        ! the output
        l_jg_active = (p_onl%dom(1) <= 0)
        DO iphys_dom=1,max_phys_dom
          IF (p_onl%dom(iphys_dom) > 0) THEN
            l_jg_active = l_jg_active .OR. (jg == p_phys_patch(p_onl%dom(iphys_dom))%logical_id)
          END IF
        END DO

        ! Selection criteria: 
        ! - domain is requested
        ! - p-level interpolation is requested
        IF (l_jg_active .AND. (p_onl%pl_varlist(1) /= ' ')) THEN
          ! check, if "all" variables are desired:
          IF (toupper(TRIM(p_onl%pl_varlist(1))) == 'ALL') THEN
            IF (dbg_level > 8)  CALL message(routine, "ALL vars for p-level")
!CDIR NOIEXPAND
            vert_intp_type(:) = vintp_types("P")
            CALL get_all_var_names(pl_varlist, nvars_pl,              &
              &                    opt_vert_intp_type=vert_intp_type, &
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
        l_jg_active = (p_onl%dom(1) <= 0)
        DO iphys_dom=1,max_phys_dom
          IF (p_onl%dom(iphys_dom) > 0) THEN
            l_jg_active = l_jg_active .OR. (jg == p_phys_patch(p_onl%dom(iphys_dom))%logical_id)
          END IF
        END DO
        ! Selection criteria: 
        ! - domain is requested
        ! - z-level interpolation is requested
        IF (l_jg_active .AND. (p_onl%hl_varlist(1) /= ' ')) THEN
          ! check, if "all" variables are desired:
          IF (toupper(TRIM(p_onl%hl_varlist(1))) == 'ALL') THEN
            IF (dbg_level > 8)  CALL message(routine, "ALL vars for z-level")
!CDIR NOIEXPAND
            vert_intp_type(:) = vintp_types("Z")
            CALL get_all_var_names(hl_varlist, nvars_hl,               &
              &                    opt_vert_intp_type=vert_intp_type,  &
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

      ! loop in search of i-level interpolation
      p_onl => first_output_name_list
      nvars_il = 0
      l_intp_i = .FALSE.
      NML_LOOP_I : DO
        IF (.NOT.ASSOCIATED(p_onl)) EXIT NML_LOOP_I
      
        ! If dom(:) was not specified in namelist input, it is set
        ! completely to -1.  In this case all domains are wanted in
        ! the output
        l_jg_active = (p_onl%dom(1) <= 0)
        DO iphys_dom=1,max_phys_dom
          IF (p_onl%dom(iphys_dom) > 0) THEN
            l_jg_active = l_jg_active .OR. (jg == p_phys_patch(p_onl%dom(iphys_dom))%logical_id)
          END IF
        END DO
        ! Selection criteria: 
        ! - domain is requested
        ! - i-level interpolation is requested
        IF (l_jg_active .AND. (p_onl%il_varlist(1) /= ' ')) THEN
          ! check, if "all" variables are desired:
          IF (toupper(TRIM(p_onl%il_varlist(1))) == 'ALL') THEN
            IF (dbg_level > 8)  CALL message(routine, "ALL vars for i-level")
!CDIR NOIEXPAND
            vert_intp_type(:) = vintp_types("I")
            CALL get_all_var_names(il_varlist, nvars_il, &
              &                    opt_vert_intp_type=vert_intp_type,   &
              &                    opt_loutput=.TRUE.)
            IF (nvars_il > 0) l_intp_i = .TRUE.
            EXIT NML_LOOP_I
          END IF

          l_intp_i = .TRUE.
          DO ivar=1,max_var_il
            IF (p_onl%il_varlist(ivar) == ' ') CYCLE
            nvars_il=nvars_il+1
            il_varlist(nvars_il) = p_onl%il_varlist(ivar)
          END DO
        END IF
        p_onl => p_onl%next
      END DO NML_LOOP_I

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
      ! now, we have total variables lists "hl_varlist(1:nvars_hl)"
      ! and "pl_varlist(1:nvars_pl)" and "il_varlist(1:nvars_il)"
      
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

      shape3d = (/ nproma, nh_pzlev_config(jg)%nzlev, nblks_c /)
        
      IF (l_intp_z) THEN
        ! temp         (nproma,nzlev,nblks_c)        
        cf_desc    = t_cf_var('temperature', 'K', 'temperature', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_opt_diag_list_z, "temp", p_diag_pz%z_temp,              &
          & GRID_UNSTRUCTURED_CELL, ZA_ALTITUDE, cf_desc, grib2_desc,           &
          & ldims=shape3d )

        ! pres         (nproma,nzlev,nblks_c)
        cf_desc    = t_cf_var('pressure', 'Pa', 'pressure', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_opt_diag_list_z, "pres", p_diag_pz%z_pres,              &
          & GRID_UNSTRUCTURED_CELL, ZA_ALTITUDE, cf_desc, grib2_desc,           &
          & ldims=shape3d )
      END IF
      IF (l_intp_p) THEN
        shape3d = (/ nproma, nh_pzlev_config(jg)%nplev, nblks_c /)
        ! GEOPOT
        cf_desc    = t_cf_var('gh', 'm', 'geopotential height', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 3, 5, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_opt_diag_list_p, 'gh', p_diag_pz%p_geopot,              &
          & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE, cf_desc, grib2_desc,           &
          & ldims=shape3d )
        ! temp
        cf_desc    = t_cf_var('temperature', 'K', 'temperature', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_opt_diag_list_p, 'temp', p_diag_pz%p_temp,              &
          & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE, cf_desc, grib2_desc,           &
          & ldims=shape3d )
      END IF
      IF (l_intp_i) THEN
        shape3d = (/ nproma, nh_pzlev_config(jg)%nilev, nblks_c /)
        ! GEOPOT
        cf_desc    = t_cf_var('gh', 'm', 'geopotential height', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 3, 5, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_opt_diag_list_i, 'gh', p_diag_pz%i_geopot,              &
          & GRID_UNSTRUCTURED_CELL, ZA_ISENTROPIC, cf_desc, grib2_desc,         &
          & ldims=shape3d )
        ! temp
        cf_desc    = t_cf_var('temperature', 'K', 'temperature', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_opt_diag_list_i, 'temp', p_diag_pz%i_temp,              &
          & GRID_UNSTRUCTURED_CELL, ZA_ISENTROPIC, cf_desc, grib2_desc,         &
          & ldims=shape3d )
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
        task%job_type                    = init_tasks(i)
        task%job_name                    = "Init: "//init_names(i)//"-level interpolation, DOM "//TRIM(int2string(jg))
      END DO

      !-- register clean-up routine as a post-processing task

      task => pp_task_insert(LOW_PRIORITY)
      task%activity     = new_simulation_status(l_output_step=.TRUE.)
      task%activity%ldom_active(jg) =  .TRUE.
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
      IF (l_intp_z) &
        CALL difference(hl_varlist, nvars_hl, (/ "temp  ", "pres  " /), 2)
      IF (l_intp_p) &
        CALL difference(pl_varlist, nvars_pl, (/ "gh    ", "temp  " /), 2)
      IF (l_intp_i) &
        CALL difference(il_varlist, nvars_il, (/ "gh    ", "temp  " /), 2)

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
        END IF
        IF (iaxis == 2) THEN
          prefix  =  "p-level"
          varlist => pl_varlist
          nvars   =  nvars_pl
          nlev    =  nh_pzlev_config(jg)%nplev
          vgrid   =  ZA_PRESSURE
          p_opt_diag_list => p_opt_diag_list_p
          job_type = TASK_INTP_VER_PLEV
        END IF
        IF (iaxis == 3) THEN
          prefix  =  "i-level"
          varlist => il_varlist
          nvars   =  nvars_il
          nlev    =  nh_pzlev_config(jg)%nilev
          vgrid   =  ZA_ISENTROPIC
          p_opt_diag_list => p_opt_diag_list_i
          job_type = TASK_INTP_VER_ILEV
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
                IF (.NOT. info%vert_interp%vert_intp_type(z_id)) CYCLE
              END IF
              IF (iaxis == 2) THEN
                IF (.NOT. info%vert_interp%vert_intp_type(p_id)) CYCLE
              END IF
              IF (iaxis == 3) THEN
                IF (.NOT. info%vert_interp%vert_intp_type(i_id)) CYCLE
              END IF

              ! Check for matching name
              idx = INDEX(info%name,'.TL')
              IF(idx == 0) THEN
                IF (varlist(ivar) /= info%name) CYCLE
              ELSE
                IF (varlist(ivar) /= info%name(1:idx-1)) CYCLE
              ENDIF

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
                &           loutput=.TRUE., new_element=new_element)

              !-- add post-processing task for interpolation

              task => pp_task_insert(DEFAULT_PRIORITY0)
              task%job_name        =  &
                &  TRIM(prefix)//" interp. "//TRIM(info%name)  &
                &  //", DOM "//TRIM(int2string(jg))
              IF (dbg_level > 8) CALL message(routine, task%job_name)

              task%job_type                    =  job_type
              task%activity                    =  new_simulation_status(l_output_step=.TRUE.)
              task%activity%ldom_active(jg)    =  .TRUE.
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
        CALL pp_task_sync()
        ! compute mean sea level pressure:
      CASE ( TASK_INTP_MSL )
        CALL pp_task_intp_msl(ptr_task)
        ! compute relative humidty
      CASE ( TASK_COMPUTE_RH )
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
  FUNCTION pp_task_is_active(ptr_task, simulation_status)
    LOGICAL :: pp_task_is_active
    TYPE(t_job_queue), POINTER :: ptr_task
    TYPE(t_simulation_status),  INTENT(IN) :: simulation_status

    INTEGER :: i

    ! compare simulation status to post-processing tasks activity
    ! flags, then check if any of the activity conditions is
    ! fulfilled:
    pp_task_is_active = .FALSE.
    IF (ANY(ptr_task%activity%status_flags(:)  .AND.  &
      &     simulation_status%status_flags(:))) pp_task_is_active = .TRUE.

    ! check, if current task applies only to domains which are
    ! "active":
    DO i=1,n_dom
      IF (ptr_task%activity%ldom_active(i) .AND. &
        & .NOT. simulation_status%ldom_active(i)) THEN
        pp_task_is_active = .FALSE.
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
  !
  ! Quasi-constructor for "t_simulation_status" variables
  ! 
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION new_simulation_status(l_output_step, l_first_step, l_last_step, &
    &                            l_dom_active)  &
    RESULT(sim_status)

    TYPE(t_simulation_status) :: sim_status
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  l_output_step, l_first_step, l_last_step
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  l_dom_active(:)
    ! local variables
    INTEGER :: ndom

    ! set default values
    sim_status%status_flags(:) = (/ .FALSE., .FALSE., .FALSE. /)

    ! supersede with user definitions
    CALL assign_if_present(sim_status%status_flags(1), l_output_step)
    CALL assign_if_present(sim_status%status_flags(2), l_first_step)
    CALL assign_if_present(sim_status%status_flags(3), l_last_step)

    ! as a default, all domains are "inactive":
    sim_status%ldom_active(:) = .FALSE.
    IF  (PRESENT(l_dom_active)) THEN
      ndom = SIZE(l_dom_active)
      sim_status%ldom_active(1:ndom) = l_dom_active(1:ndom)
    END IF

  END FUNCTION new_simulation_status

END MODULE mo_pp_scheduler
