!>
!! Tasks for internal post-processing.
!!
!! The subroutines in this module can be inserted into a dynamic "job queue".
!! See module "mo_pp_scheduler" for detailed info.
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
!! -----------------------------------------------------------------------------------
MODULE mo_pp_tasks

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message, message_text, finish
  USE mo_impl_constants,          ONLY: SUCCESS,                      &
    & VINTP_TYPE_Z, VINTP_TYPE_P_OR_Z, VINTP_TYPE_NONE,               &
    & VINTP_METHOD_UV, VINTP_METHOD_LIN, HINTP_TYPE_NONE,             &     
    & VINTP_METHOD_QV, HINTP_TYPE_LONLAT, VINTP_METHOD_LIN_NLEVP1,    &
    & TASK_NONE, TASK_INIT_VER_PZ, TASK_INIT_VER_Z, TASK_FINALIZE_PZ, &
    & TASK_INTP_HOR_LONLAT, TASK_INTP_VER_PLEV, TASK_INTP_SYNC,       &
    & TASK_COMPUTE_RH, TASK_INTP_VER_ZLEV
  USE mo_model_domain,            ONLY: t_patch, p_patch
  USE mo_var_list_element,        ONLY: t_var_list_element, level_type_ml,  &
    &                                   level_type_pl, level_type_hl
  USE mo_var_metadata,            ONLY: t_var_metadata, t_vert_interp_meta
  USE mo_intp,                    ONLY: verts2cells_scalar
  USE mo_intp_data_strc,          ONLY: t_int_state, lonlat_grid_list, &
    &                                   t_lon_lat_intp, p_int_state
  USE mo_nh_vert_interp,          ONLY: prepare_vert_interp, &
    &                                   lin_intp, uv_intp, qv_intp, &
    &                                   pressure_intp_msl
  USE mo_nonhydro_types,          ONLY: t_nh_state, t_nh_prog, t_nh_diag, &
    &                                   t_nh_metrics
  USE mo_nonhydro_state,          ONLY: p_nh_state
  USE mo_opt_diagnostics,         ONLY: t_nh_diag_pz, t_nh_opt_diag, t_vcoeff, &
    &                                   vcoeff_deallocate, p_nh_opt_diag
  USE mo_nwp_phy_types,           ONLY: t_nwp_phy_diag
  USE mo_nwp_phy_state,           ONLY: prm_diag
  USE mo_nh_pzlev_config,         ONLY: t_nh_pzlev_config, nh_pzlev_config
  USE mo_name_list_output_config, ONLY: t_output_name_list, &
    &                                   vname_len, first_output_name_list
  USE mo_parallel_config,         ONLY: nproma
  USE mo_dynamics_config,         ONLY: nnow
  USE mo_cdi_constants,           ONLY: GRID_CELL, GRID_REFERENCE,               &
    &                                   GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,  &
    &                                   ZAXIS_PRESSURE, GRID_REGULAR_LONLAT,     &
    &                                   GRID_UNSTRUCTURED_EDGE,                  &
    &                                   GRID_UNSTRUCTURED_VERT, ZAXIS_SURFACE,   &
    &                                   DATATYPE_FLT32, DATATYPE_PACK16
  USE mo_linked_list,             ONLY: t_var_list, t_list_element
  USE mo_lonlat_grid,             ONLY: t_lon_lat_grid
  USE mo_intp_lonlat,             ONLY: rbf_interpol_lonlat_nl, &
    &                                   rbf_vec_interpol_lonlat_nl
  USE mo_sync,                    ONLY: sync_patch_array,                        &
    &                                   SYNC_C, SYNC_E, SYNC_V
  USE mo_util_phys,               ONLY: compute_field_rel_hum

  IMPLICIT NONE

  ! interface definition
  PRIVATE

  ! max. name string length
  INTEGER, PARAMETER, PUBLIC :: MAX_NAME_LENGTH   =   64

  ! priority levels for tasks (smaller is earlier):
  INTEGER, PARAMETER, PUBLIC  :: HIGH_PRIORITY     =    0  
  INTEGER, PARAMETER, PUBLIC  :: DEFAULT_PRIORITY0 =    8  
  INTEGER, PARAMETER, PUBLIC  :: DEFAULT_PRIORITY1 =    9   
  INTEGER, PARAMETER, PUBLIC  :: DEFAULT_PRIORITY2 =   10  
  INTEGER, PARAMETER, PUBLIC  :: LOW_PRIORITY      =  100  

  ! level of output verbosity
  INTEGER, PUBLIC :: dbg_level = 0

  ! functions and subroutines
  PUBLIC :: pp_task_lonlat
  PUBLIC :: pp_task_sync
  PUBLIC :: pp_task_pzlev_setup
  PUBLIC :: pp_task_pzlev
  PUBLIC :: pp_task_intp_msl
  PUBLIC :: pp_task_compute_field
  ! variables
  PUBLIC :: job_queue
  ! data types
  PUBLIC :: t_data_input
  PUBLIC :: t_data_output
  PUBLIC :: t_simulation_status
  PUBLIC :: t_job_queue

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


  !--- MODULE DATA -------------------------------------------------------------------
  TYPE(t_job_queue), POINTER   :: job_queue  =>  NULL() !< head of (ordered) job queue


CONTAINS

  !--- POST-PROCESSING TASKS ---------------------------------------------------------

  !---------------------------------------------------------------
  !> Performs interpolation of a variable onto a regular grid.
  !
  !  This is only a wrapper for the corresponding routines from the
  !  interpolation module.
  SUBROUTINE pp_task_lonlat(ptr_task)
    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_pp_tasks:pp_task_lonlat")
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
        tmp_var(:,1,:) = in_var%r_ptr(:,:,in_var_idx,1,1)

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
        tmp_var(:,1,:) = in_var%r_ptr(:,:,in_var_idx,1,1)

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
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_pp_tasks:pp_task_sync")
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
            CALL sync_patch_array(SYNC_C, p_patch, in_var%r_ptr(:,:,in_var_idx,1,1) )
          ELSE
            CALL sync_patch_array(SYNC_C, p_patch, in_var%r_ptr(:,:,:,in_var_idx,1) )
          END IF
        CASE (GRID_UNSTRUCTURED_EDGE)
          IF (p_info%vgrid == ZAXIS_SURFACE) THEN
            CALL sync_patch_array(SYNC_E, p_patch, in_var%r_ptr(:,:,in_var_idx,1,1) )
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
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_pp_tasks:pp_task_pzlev_setup")
    INTEGER                            :: jg, nzlev, nplev
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
      IF (dbg_level >= 10)  CALL message(routine, "TASK_INIT_VER_PZ")
      CALL prepare_vert_interp(p_patch, p_prog, p_diag, prm_diag, nzlev, nplev,  & ! in
        &                      p_diag_pz%z_temp, p_diag_pz%z_tracer_iqv,         & ! inout
        &                      p_diag_pz%z_tot_cld_iqv,                          & ! inout
        &                      p_diag_pz%z_pres, p_diag_pz%p_geopot,             & ! inout
        &                      p_diag_pz%p_temp,                                 & ! inout
        &                      nh_pzlev_config%p3d, nh_pzlev_config%z3d,         & ! in
        &                      p_metrics,                                        & ! in
        &                      vcoeff_z=p_diag_pz%vcoeff_z, vcoeff_p=p_diag_pz%vcoeff_p )            ! inout
      !
    CASE ( TASK_INIT_VER_Z )
      ! build data structure "vcoeff" containing coefficient tables
      IF (dbg_level >= 10)  CALL message(routine, "TASK_INIT_VER_PZ")
      CALL prepare_vert_interp(p_patch, p_prog, p_diag, prm_diag, nzlev, nplev,  & ! in
        &                      p_diag_pz%z_temp, p_diag_pz%z_tracer_iqv,         & ! inout
        &                      p_diag_pz%z_tot_cld_iqv,                          & ! inout
        &                      p_diag_pz%z_pres, p_diag_pz%p_geopot,             & ! inout
        &                      p_diag_pz%p_temp,                                 & ! inout
        &                      nh_pzlev_config%p3d, nh_pzlev_config%z3d,         & ! in
        &                      p_metrics,                                        & ! in
        &                      vcoeff_z=p_diag_pz%vcoeff_z )                                ! inout
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
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_pp_tasks:pp_task_pzlev")
    INTEGER                            :: &
      &  vert_intp_method, jg,                  &
      &  in_var_idx, out_var_idx, nlev, nlevp1, &
      &  nzlev, nplev, npzlev, npromz, nblks,   &
      &  dim2, ierrstat
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
    REAL(wp), ALLOCATABLE              :: tmp_var(:,:,:)

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
    SELECT CASE ( ptr_task%job_type )
    CASE ( TASK_INTP_VER_ZLEV )
      ! vertical levels for z-level interpolation
      npzlev  =   nzlev
      vcoeff  =>  p_diag_pz%vcoeff_z
      p_z3d   =>  nh_pzlev_config%z3d
    CASE ( TASK_INTP_VER_PLEV )
      ! vertical levels for p-level interpolation
      npzlev  =   nplev
      vcoeff  =>  p_diag_pz%vcoeff_p
      p_z3d   =>  p_diag_pz%p_geopot
    CASE DEFAULT
      CALL finish(routine, "Unknown post-processing job.")
    END SELECT
                     
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

    nblks  = p_patch%nblks_c
    npromz = p_patch%npromz_c

    dim2 = UBOUND(in_var%r_ptr, 2)
    ALLOCATE(tmp_var(nproma, dim2, nblks), STAT=ierrstat)
    IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation failed')

    SELECT CASE ( p_info%hgrid )
    CASE (GRID_UNSTRUCTURED_CELL) 
      tmp_var(:,:,:) = in_var%r_ptr(:,:,:,in_var_idx,1)
    CASE (GRID_UNSTRUCTURED_VERT) 
      CALL verts2cells_scalar( in_var%r_ptr(:,:,:,in_var_idx,1), p_patch,      &
        &                      ptr_task%data_input%p_int_state%verts_aw_cells, &
        &                      tmp_var(:,:,:), 1, dim2 )
    END SELECT


    !--- actually perform vertical interpolation task
    SELECT CASE ( vert_intp_method )
    CASE ( VINTP_METHOD_UV )
      IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_UV")
      CALL uv_intp(tmp_var(:,:,:),                                       & !in
        &          out_var%r_ptr(:,:,:,out_var_idx,1),                   & !out
        &          p_metrics%z_mc, p_z3d,                                & !in
        &          nblks, npromz, nlev, npzlev,                          & !in
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
      CALL lin_intp(tmp_var(:,:,:),                                      & !inout
        &           out_var%r_ptr(:,:,:,out_var_idx,1),                  & !out
        &           nblks, npromz, nlev, npzlev,                         & !in
        &           vcoeff%wfac_lin, vcoeff%idx0_lin,                    & !in
        &           vcoeff%bot_idx_lin, vcoeff%wfacpbl1,                 & !in
        &           vcoeff%kpbl1, vcoeff%wfacpbl2, vcoeff%kpbl2,         & !in
        &           l_loglin=l_loglin,                                   & !in
        &           l_extrapol=l_extrapol, l_pd_limit=l_pd_limit,        & !in
        &           lower_limit=lower_limit )                              !in
      !
    CASE ( VINTP_METHOD_LIN_NLEVP1 )        
      IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_LIN_NLEVP1")
      CALL lin_intp(tmp_var(:,:,:),                                      & !inout
        &           out_var%r_ptr(:,:,:,out_var_idx,1),                  & !out
        &           nblks, npromz, nlevp1, npzlev,                       & !in
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
      CALL qv_intp(tmp_var(:,:,:),                                     & !in
        &          out_var%r_ptr(:,:,:,out_var_idx,1),                 & !out
        &          p_metrics%z_mc, p_z3d, p_diag%temp,                 & !in
        &          p_diag%pres, p_diag_pz%p_temp, nh_pzlev_config%p3d, & !in
        &          nblks, npromz, nlev, npzlev,                        & !in
        &          vcoeff%coef1, vcoeff%coef2, vcoeff%coef3,           & !in
        &          vcoeff%wfac_lin, vcoeff%idx0_cub, vcoeff%idx0_lin,  & !in
        &          vcoeff%bot_idx_cub, vcoeff%bot_idx_lin,             & !in
        &          vcoeff%wfacpbl1, vcoeff%kpbl1,                      & !in
        &          vcoeff%wfacpbl2, vcoeff%kpbl2,                      & !in
        &          l_satlimit=l_satlimit, lower_limit=lower_limit,     & !in
        &          l_restore_pbldev=l_restore_pbldev )                   !in
    END SELECT ! vert_intp_method

    ! clean up
    DEALLOCATE(tmp_var, STAT=ierrstat)
    IF (ierrstat /= SUCCESS)  CALL finish (routine, 'deallocation failed')

  END SUBROUTINE pp_task_pzlev


  !---------------------------------------------------------------
  !> Performs interpolation of a 2D field onto mean sea level, z=0.
  !
  !  This routine is completely independent from the data structures
  !  used for pz-level interpolation.
  !
  !  @note It could save computational capacity, if the temporary
  !        fields for temperature etc. would be allocated once and
  !        destroyed once for each task.
  !
  SUBROUTINE pp_task_intp_msl(ptr_task)
    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables    
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_pp_tasks:pp_task_intp_msl")
    LOGICAL, PARAMETER :: lpres_msl_gme = .TRUE.

    REAL(wp), ALLOCATABLE              :: z_temp(:,:,:), z_tracer_iqv(:,:,:),  &
      &                                   z_tot_cld_iqv(:,:,:), z_pres(:,:,:)
    REAL(wp), POINTER                  :: z_0(:,:,:)
    REAL(wp)                           :: lower_limit
    INTEGER                            :: nzlev, nplev, nblks, npromz, jg,     &
      &                                   in_var_idx, out_var_idx, dim1, dim2, &
      &                                   ierrstat, sea_lev
    LOGICAL                            :: l_extrapol, l_pd_limit, l_loglin,    &
      &                                   l_only_p2z
    TYPE(t_nh_pzlev_config),   POINTER :: nh_pzlev_config
    TYPE(t_vcoeff)                     :: vcoeff
    TYPE(t_vert_interp_meta),  POINTER :: pzlev_flags
    TYPE (t_var_list_element), POINTER :: in_var, out_var
    TYPE(t_var_metadata),      POINTER :: p_info
    REAL(wp), ALLOCATABLE              :: tmp_var_out(:,:,:)
    TYPE(t_patch),             POINTER :: p_patch
    TYPE(t_nh_metrics),        POINTER :: p_metrics    
    TYPE(t_nh_prog),           POINTER :: p_prog
    TYPE(t_nh_diag),           POINTER :: p_diag
    TYPE(t_nwp_phy_diag),      POINTER :: prm_diag

    ! just a single z-level... 
    nzlev   = 1
    sea_lev = 1
    nplev   = 0

    ! patch, state, and metrics
    jg             =  ptr_task%data_input%jg
    p_patch        => ptr_task%data_input%p_patch
    p_metrics      => ptr_task%data_input%p_nh_state%metrics
    p_prog         => ptr_task%data_input%p_nh_state%prog(nnow(jg))
    p_diag         => ptr_task%data_input%p_nh_state%diag
    prm_diag       => ptr_task%data_input%prm_diag

    ! pz-level interpolation data
    nh_pzlev_config   => ptr_task%data_input%nh_pzlev_config

    ! input/output field for this task
    p_info            => ptr_task%data_input%var%info
    in_var            => ptr_task%data_input%var
    out_var           => ptr_task%data_output%var

    IF (TRIM(p_info%name) /= "pres")  CALL message(routine, "Invalid input field!")

    ! interpolation flags + parameters
    pzlev_flags => in_var%info%vert_interp
    l_loglin       = pzlev_flags%l_loglin          
    l_extrapol     = pzlev_flags%l_extrapol        
    l_pd_limit     = pzlev_flags%l_pd_limit
    lower_limit    = pzlev_flags%lower_limit       

    SELECT CASE ( p_info%hgrid )
    CASE (GRID_UNSTRUCTURED_CELL) 
      nblks  = p_patch%nblks_c
      npromz = p_patch%npromz_c
    CASE (GRID_UNSTRUCTURED_VERT) 
      nblks  = p_patch%nblks_v
      npromz = p_patch%npromz_v
    END SELECT

    in_var_idx  = 1
    IF (in_var%info%lcontained)  in_var_idx  = in_var%info%ncontained
    out_var_idx = 1
    IF (out_var%info%lcontained) out_var_idx = out_var%info%ncontained

    IF (.NOT. lpres_msl_gme) THEN
      ! First create a 1-level 3D variable for output  (nproma, nlevs, nblks):
      dim1 = p_info%used_dimensions(1)
      dim2 = p_info%used_dimensions(3)
      ALLOCATE(tmp_var_out(dim1, nzlev, dim2), STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation failed')

      ! allocate temporary variables for temperature, etc.
      ALLOCATE(z_temp(dim1, nzlev, dim2), z_tracer_iqv(dim1, nzlev, dim2),  &
        &      z_tot_cld_iqv(dim1, nzlev, dim2), z_pres(dim1, nzlev, dim2), &
        &      z_0(dim1, nzlev, dim2), STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation failed')

      ! initialize trivial height field:
      z_0(:,sea_lev,:) = 0._wp

      ! build data structure "vcoeff" containing coefficient tables
      CALL prepare_vert_interp(p_patch, p_prog, p_diag, prm_diag, nzlev, nplev, & ! in
        &           z_temp, z_tracer_iqv, z_tot_cld_iqv, z_pres,                & ! inout
        &           p_z3d_out=z_0, p_metrics=p_metrics,                         & ! in
        &           vcoeff_z=vcoeff, l_only_p2z=.TRUE. )                          ! inout, in

      ! copy back pressure:
      out_var%r_ptr(:,:,out_var_idx,1,1) = z_pres(:,sea_lev,:)

      ! clean up:
      DEALLOCATE(tmp_var_out, STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'deallocation failed')
      DEALLOCATE(z_temp, z_tracer_iqv, z_tot_cld_iqv, z_pres, z_0, STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'deallocation failed')
      ! deallocate coefficient tables:
      CALL vcoeff_deallocate(vcoeff)

    ELSE

      ! Interpolate pressure on z-levels
      CALL pressure_intp_msl(p_diag%pres, p_diag%pres_sfc, p_diag%temp, &  ! in
        &                    p_metrics%z_ifc,                           &  ! in
        &                    out_var%r_ptr(:,:,out_var_idx,1,1),        &  ! out
        &                    nblks, npromz, p_patch%nlev )                 ! in
      
    END IF

  END SUBROUTINE pp_task_intp_msl


  !---------------------------------------------------------------
  !> Performs computation of optional diagnostic fields.
  !
  !  Selects subroutines for field computation based on the variable
  !  name.
  !
  !  @note This could be easily replaced by a procedure pointer,
  !        alas, this is an F2003 feature.
  !
  !  @todo Change order of processing: First, interpolate input fields
  !        onto z-levels, then compute rel_hum.
  !
  SUBROUTINE pp_task_compute_field(ptr_task)
    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_pp_tasks:pp_task_compute_field")
    INTEGER                            :: jg, out_var_idx
    TYPE (t_var_list_element), POINTER :: out_var
    TYPE(t_var_metadata),      POINTER :: p_info
    TYPE(t_patch),             POINTER :: p_patch
    TYPE(t_nh_prog),           POINTER :: p_prog
    TYPE(t_nh_diag),           POINTER :: p_diag
    
    ! output field for this task
    out_var   => ptr_task%data_output%var
    p_info    => out_var%info    
    out_var   => out_var
    out_var_idx = 1
    if (out_var%info%lcontained)  out_var_idx = out_var%info%ncontained

    ! input data required for computation:
    jg        =  ptr_task%data_input%jg
    p_patch   => ptr_task%data_input%p_patch
    p_prog    => ptr_task%data_input%p_nh_state%prog(nnow(jg))
    p_diag    => ptr_task%data_input%p_nh_state%diag

    SELECT CASE(ptr_task%job_type)
    CASE (TASK_COMPUTE_RH)
      CALL compute_field_rel_hum(p_patch, p_prog, p_diag, &
        &                        out_var%r_ptr(:,:,:,out_var_idx,1))
    END SELECT

  END SUBROUTINE pp_task_compute_field

END MODULE mo_pp_tasks
