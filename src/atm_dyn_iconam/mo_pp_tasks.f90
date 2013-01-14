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
    & VINTP_METHOD_UV, VINTP_METHOD_LIN, VINTP_METHOD_QV,             &
    & VINTP_METHOD_LIN_NLEVP1,                                        &
    & TASK_NONE, TASK_INIT_VER_Z, TASK_INIT_VER_P, TASK_INIT_VER_I,   &
    & TASK_FINALIZE_IPZ,                                              &
    & TASK_INTP_HOR_LONLAT, TASK_INTP_VER_PLEV, TASK_INTP_SYNC,       &
    & TASK_COMPUTE_RH, TASK_INTP_VER_ZLEV, TASK_INTP_VER_ILEV,        &
    & PRES_MSL_METHOD_SAI, PRES_MSL_METHOD_GME, max_dom
  USE mo_model_domain,            ONLY: t_patch, p_patch
  USE mo_var_list_element,        ONLY: t_var_list_element, level_type_ml,  &
    &                                   level_type_pl, level_type_hl
  USE mo_var_metadata,            ONLY: t_var_metadata, t_vert_interp_meta
  USE mo_intp,                    ONLY: verts2cells_scalar, cell_avg,       &
    &                                   cells2edges_scalar
  USE mo_intp_data_strc,          ONLY: t_int_state, lonlat_grid_list,      &
    &                                   t_lon_lat_intp, p_int_state
  USE mo_intp_rbf,                ONLY: rbf_vec_interpol_cell
  USE mo_nh_vert_interp,          ONLY: prepare_vert_interp_z,              &
    &                                   prepare_vert_interp_p,              &
    &                                   prepare_vert_interp_i,              &
    &                                   lin_intp, uv_intp, qv_intp,         &
    &                                   diagnose_pmsl, diagnose_pmsl_gme,   &
    &                                   prepare_extrap
  USE mo_nonhydro_types,          ONLY: t_nh_state, t_nh_prog, t_nh_diag,   &
    &                                   t_nh_metrics
  USE mo_nonhydro_state,          ONLY: p_nh_state
  USE mo_opt_diagnostics,         ONLY: t_nh_diag_pz, t_nh_opt_diag, t_vcoeff, &
    &                                   vcoeff_allocate, vcoeff_deallocate,    &
    &                                   p_nh_opt_diag, t_vcoeff_lin,           &
    &                                   t_vcoeff_cub
  USE mo_nwp_phy_types,           ONLY: t_nwp_phy_diag
  USE mo_nwp_phy_state,           ONLY: prm_diag
  USE mo_nh_pzlev_config,         ONLY: t_nh_pzlev_config, nh_pzlev_config
  USE mo_name_list_output_config, ONLY: t_output_name_list, &
    &                                   vname_len, first_output_name_list
  USE mo_parallel_config,         ONLY: nproma
  USE mo_dynamics_config,         ONLY: nnow
  USE mo_cdi_constants,           ONLY: GRID_CELL, GRID_REFERENCE,               &
    &                                   GRID_UNSTRUCTURED_CELL, GRID_REGULAR_LONLAT, &
    &                                   GRID_UNSTRUCTURED_EDGE,                  &
    &                                   DATATYPE_FLT32,                          &
    &                                   DATATYPE_PACK16, is_2d_field
  USE mo_linked_list,             ONLY: t_var_list, t_list_element
  USE mo_lonlat_grid,             ONLY: t_lon_lat_grid
  USE mo_intp_lonlat,             ONLY: rbf_interpol_lonlat_nl,                  &
    &                                   rbf_vec_interpol_lonlat_nl
  USE mo_sync,                    ONLY: sync_patch_array,                        &
    &                                   SYNC_C, SYNC_E, SYNC_V
  USE mo_util_phys,               ONLY: compute_field_rel_hum
  USE mo_io_config,               ONLY: itype_pres_msl

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
  INTEGER, PARAMETER, PUBLIC  :: DEFAULT_PRIORITY3 =   11  
  INTEGER, PARAMETER, PUBLIC  :: LOW_PRIORITY      =  100  

  ! level of output verbosity
  INTEGER, PUBLIC :: dbg_level = 0

  ! functions and subroutines
  PUBLIC :: pp_task_lonlat
  PUBLIC :: pp_task_sync
  PUBLIC :: pp_task_ipzlev_setup
  PUBLIC :: pp_task_ipzlev
  PUBLIC :: pp_task_intp_msl
  PUBLIC :: pp_task_compute_field
  PUBLIC :: pp_task_edge2cell
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
    ! active domains
    LOGICAL :: ldom_active(max_dom)
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

    IF (is_2d_field(p_info%vgrid) .AND. (p_info%ndims /= 2)) THEN
      CALL finish(routine, "Inconsistent dimension info!")
    END IF

    SELECT CASE (p_info%hgrid)
    CASE (GRID_UNSTRUCTURED_CELL)
      IF (is_2d_field(p_info%vgrid)) THEN
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
      IF (is_2d_field(p_info%vgrid)) THEN
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
        
        IF (is_2d_field(p_info%vgrid) .AND. (p_info%ndims /= 2)) THEN
          CALL finish(routine, "Inconsistent dimension info!")
        END IF

        SELECT CASE (p_info%hgrid)
        CASE (GRID_UNSTRUCTURED_CELL)
          IF (is_2d_field(p_info%vgrid)) THEN
            CALL sync_patch_array(SYNC_C, p_patch, in_var%r_ptr(:,:,in_var_idx,1,1) )
          ELSE
            CALL sync_patch_array(SYNC_C, p_patch, in_var%r_ptr(:,:,:,in_var_idx,1) )
          END IF
        CASE (GRID_UNSTRUCTURED_EDGE)
          IF (is_2d_field(p_info%vgrid)) THEN
            CALL sync_patch_array(SYNC_E, p_patch, in_var%r_ptr(:,:,in_var_idx,1,1) )
          ELSE
            CALL sync_patch_array(SYNC_E, p_patch, in_var%r_ptr(:,:,:,in_var_idx,1) )
          END IF
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
  SUBROUTINE pp_task_ipzlev_setup(ptr_task)
    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_pp_tasks:pp_task_ipzlev_setup")
    INTEGER                            :: jg, nzlev, nplev, nilev
    TYPE(t_patch),             POINTER :: p_patch
    TYPE(t_nh_metrics),        POINTER :: p_metrics    

    ! prognostic state: note that we only use p_prog(nnow(jg))
    TYPE(t_nh_prog),           POINTER :: p_prog
    TYPE(t_nh_diag),           POINTER :: p_diag
    TYPE(t_nh_diag_pz),        POINTER :: p_diag_pz
    TYPE(t_nwp_phy_diag),      POINTER :: prm_diag
    TYPE(t_nh_pzlev_config),   POINTER :: nh_pzlev_config
    TYPE(t_int_state),         POINTER :: intp_hrz

    ! patch, state, and metrics
    jg             =  ptr_task%data_input%jg
    p_patch        => ptr_task%data_input%p_patch
    p_metrics      => ptr_task%data_input%p_nh_state%metrics
    p_prog         => ptr_task%data_input%p_nh_state%prog(nnow(jg))
    p_diag         => ptr_task%data_input%p_nh_state%diag
    p_diag_pz      => ptr_task%data_input%p_nh_opt_diag%diag_pz
    prm_diag       => ptr_task%data_input%prm_diag
    intp_hrz       => ptr_task%data_input%p_int_state

    ! ipz-level interpolation data
    nh_pzlev_config   => ptr_task%data_input%nh_pzlev_config

    nzlev          =  nh_pzlev_config%nzlev
    nplev          =  nh_pzlev_config%nplev
    nilev          =  nh_pzlev_config%nilev

    ! build data structure "vcoeff" containing coefficient tables                      
    SELECT CASE ( ptr_task%job_type )
    CASE ( TASK_INIT_VER_Z )
      IF (dbg_level >= 10)  CALL message(routine, "TASK_INIT_VER_Z")
      CALL prepare_vert_interp_z(p_patch, p_diag, p_metrics, intp_hrz, nzlev,          &
        &                        p_diag_pz%z_temp, p_diag_pz%z_pres,                   &
        &                        nh_pzlev_config%z3d, p_diag_pz%vcoeff_z)
      !
    CASE ( TASK_INIT_VER_P )
      IF (dbg_level >= 10)  CALL message(routine, "TASK_INIT_VER_P")
      CALL prepare_vert_interp_p(p_patch, p_diag, p_metrics, intp_hrz, nplev,          &
        &                        p_diag_pz%p_gh, p_diag_pz%p_temp,                     &
        &                        nh_pzlev_config%p3d, p_diag_pz%vcoeff_p)
      !
    CASE ( TASK_INIT_VER_I )
      IF (dbg_level >= 10)  CALL message(routine, "TASK_INIT_VER_I")
      CALL prepare_vert_interp_i(p_patch, p_prog, p_diag, p_metrics, intp_hrz, nilev,  &
        &                        p_diag_pz%i_gh, p_diag_pz%i_temp,                     &
        &                        nh_pzlev_config%i3d, p_diag_pz%vcoeff_i)
      !
    CASE ( TASK_FINALIZE_IPZ )
      ! deallocate coefficient tables:
      CALL vcoeff_deallocate(p_diag_pz%vcoeff_z)
      CALL vcoeff_deallocate(p_diag_pz%vcoeff_p)
      CALL vcoeff_deallocate(p_diag_pz%vcoeff_i)
      !
    CASE DEFAULT
      CALL finish(routine, "Internal error.")
    END SELECT ! vert_intp_method

  END SUBROUTINE pp_task_ipzlev_setup


  !---------------------------------------------------------------
  !> Performs vertical interpolation of a variable onto i/p/z-levels.
  !
  !  This is only a wrapper for the corresponding routines from the
  !  interpolation module.
  !
  SUBROUTINE pp_task_ipzlev(ptr_task)
    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_pp_tasks:pp_task_ipzlev")
    INTEGER                            :: &
      &  vert_intp_method, jg,                    &
      &  in_var_idx, out_var_idx, nlev, nlevp1,   &
      &  n_ipzlev, npromz, nblks, dim2, ierrstat
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
    REAL(wp),                  POINTER :: p_z3d(:,:,:), p_pres(:,:,:), p_temp(:,:,:)
    REAL(wp), ALLOCATABLE              :: tmp_var(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET      :: z_me(:,:,:), p_z3d_edge(:,:,:)
    REAL(wp),                  POINTER :: in_z3d(:,:,:), in_z_mc(:,:,:)
    TYPE(t_int_state),         POINTER :: intp_hrz

    LOGICAL                            :: &
      &  l_hires_intp, l_restore_fricred, l_loglin, &
      &  l_extrapol, l_satlimit, l_restore_pbldev,  &
      &  l_pd_limit, l_restore_sfcinv, l_hires_corr
    REAL(wp)                           :: &
      &  lower_limit, extrapol_dist
    TYPE (t_vcoeff_lin), POINTER       :: &
      &  vcoeff_lin, vcoeff_lin_nlevp1
    TYPE (t_vcoeff_cub), POINTER       :: &
      &  vcoeff_cub

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
    intp_hrz          => ptr_task%data_input%p_int_state

    nh_pzlev_config   => ptr_task%data_input%nh_pzlev_config
    nlev              = p_patch%nlev
    nlevp1            = p_patch%nlevp1

    ! pz-level interpolation data
    SELECT CASE ( ptr_task%job_type )
    CASE ( TASK_INTP_VER_ZLEV )
      ! vertical levels for z-level interpolation
      n_ipzlev  =  nh_pzlev_config%nzlev
      vcoeff  =>  p_diag_pz%vcoeff_z
      p_z3d   =>  nh_pzlev_config%z3d
      p_pres  =>  p_diag_pz%z_pres
      p_temp  =>  p_diag_pz%z_temp
    CASE ( TASK_INTP_VER_PLEV )
      ! vertical levels for p-level interpolation
      n_ipzlev  =  nh_pzlev_config%nplev
      vcoeff  =>  p_diag_pz%vcoeff_p
      p_z3d   =>  p_diag_pz%p_gh
      p_pres  =>  nh_pzlev_config%p3d
      p_temp  =>  p_diag_pz%p_temp
    CASE ( TASK_INTP_VER_ILEV )
      ! vertical levels for isentropic-level interpolation
      n_ipzlev  =   nh_pzlev_config%nilev
      vcoeff  =>  p_diag_pz%vcoeff_i
      p_z3d   =>  p_diag_pz%i_gh
      p_pres  =>  nh_pzlev_config%p3d ! ** this still needs to be fixed! we either need i_pres here
      p_temp  =>  p_diag_pz%i_temp    !    or have to turn off the saturation adjustment for theta levels
                                      !    or have to use log-linear interpolation in this case **
      IF (vert_intp_method == VINTP_METHOD_QV) & ! Let's stop with an error message for the time being
        CALL finish(routine, "QV interpolation to isentropic levels not available.")
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

    SELECT CASE ( p_info%hgrid )
    CASE (GRID_UNSTRUCTURED_CELL) 
      nblks  = p_patch%nblks_c
      npromz = p_patch%npromz_c
      ! 
      vcoeff_lin        => vcoeff%lin_cell
      vcoeff_lin_nlevp1 => vcoeff%lin_cell_nlevp1
      vcoeff_cub        => vcoeff%cub_cell
      in_z3d            => p_z3d
      in_z_mc           => p_metrics%z_mc

    CASE (GRID_UNSTRUCTURED_EDGE) 
      nblks  = p_patch%nblks_e
      npromz = p_patch%npromz_e
      !
      vcoeff_lin        => vcoeff%lin_edge
      vcoeff_lin_nlevp1 => NULL()
      vcoeff_cub        => vcoeff%cub_edge

      ! Compute geometric height at edge points (temporary variable)
      ALLOCATE(p_z3d_edge(nproma,n_ipzlev,nblks), z_me(nproma,p_patch%nlev,nblks), STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed')

      CALL cells2edges_scalar(p_metrics%z_mc, p_patch, intp_hrz%c_lin_e,    &
        &                     z_me, opt_fill_latbc=.TRUE.)
      CALL cells2edges_scalar(p_z3d, p_patch, intp_hrz%c_lin_e, p_z3d_edge, &
        &                     opt_fill_latbc=.TRUE.)
      in_z3d            => p_z3d_edge
      in_z_mc           => z_me
    END SELECT

    dim2 = UBOUND(in_var%r_ptr, 2)
    ALLOCATE(tmp_var(nproma, dim2, nblks), STAT=ierrstat)
    IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation failed')

    SELECT CASE ( p_info%hgrid )
    CASE (GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE) 
      ! consistency check:
      IF ((UBOUND(in_var%r_ptr,1) > nproma) .OR.  &
        & (UBOUND(in_var%r_ptr,2) > dim2)   .OR.  &
        & (UBOUND(in_var%r_ptr,3) > nblks)) THEN
        CALL finish(routine, "Inconsistent array dimensions")
      END IF
      tmp_var(:,:,:) = in_var%r_ptr(:,:,:,in_var_idx,1)
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

    !--- actually perform vertical interpolation task
    IF (.NOT. ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0)))) THEN

      SELECT CASE ( vert_intp_method )
      CASE ( VINTP_METHOD_UV )
        IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_UV")
        IF (.NOT. ASSOCIATED(vcoeff_lin)) CALL finish(routine, "Internal error!")
        IF (.NOT. ASSOCIATED(vcoeff_cub)) CALL finish(routine, "Internal error!")
        CALL uv_intp(tmp_var(:,:,:),                                                & !in
          &          out_var%r_ptr(:,:,:,out_var_idx,1),                            & !out
          &          in_z_mc, in_z3d,                                               & !in
          &          nblks, npromz, nlev, n_ipzlev,                                 & !in
          &          vcoeff_cub%coef1, vcoeff_cub%coef2,                            & !in
          &          vcoeff_cub%coef3, vcoeff_lin%wfac_lin,                         & !in
          &          vcoeff_cub%idx0_cub, vcoeff_lin%idx0_lin,                      & !in
          &          vcoeff_cub%bot_idx_cub, vcoeff_lin%bot_idx_lin,                & !in
          &          vcoeff_lin%wfacpbl1, vcoeff_lin%kpbl1,                         & !in
          &          vcoeff_lin%wfacpbl2, vcoeff_lin%kpbl2,                         & !in
          &          l_hires_intp=l_hires_intp,                                     & !in
          &          l_restore_fricred=l_restore_fricred )                            !in
        !
      CASE ( VINTP_METHOD_LIN )        
        IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_LIN")
        IF (.NOT. ASSOCIATED(vcoeff_lin)) CALL finish(routine, "Internal error!")
        CALL lin_intp(tmp_var(:,:,:),                                               & !inout
          &           out_var%r_ptr(:,:,:,out_var_idx,1),                           & !out
          &           nblks, npromz, nlev, n_ipzlev,                                & !in
          &           vcoeff_lin%wfac_lin, vcoeff_lin%idx0_lin,                     & !in
          &           vcoeff_lin%bot_idx_lin,                                       & !in
          &           vcoeff_lin%wfacpbl1, vcoeff_lin%kpbl1,                        & !in
          &           vcoeff_lin%wfacpbl2, vcoeff_lin%kpbl2,                        & !in
          &           l_loglin=l_loglin,                                            & !in
          &           l_extrapol=l_extrapol, l_pd_limit=l_pd_limit,                 & !in
          &           lower_limit=lower_limit )                                       !in
        !
      CASE ( VINTP_METHOD_LIN_NLEVP1 )        
        IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_LIN_NLEVP1")
        IF (.NOT. ASSOCIATED(vcoeff_lin_nlevp1)) CALL finish(routine, "Internal error!")
        CALL lin_intp(tmp_var(:,:,:),                                               & !inout
          &           out_var%r_ptr(:,:,:,out_var_idx,1),                           & !out
          &           nblks, npromz, nlevp1, n_ipzlev,                              & !in
          &           vcoeff_lin_nlevp1%wfac_lin,                                   & !in
          &           vcoeff_lin_nlevp1%idx0_lin,                                   & !in
          &           vcoeff_lin_nlevp1%bot_idx_lin,                                & !in
          &           vcoeff_lin_nlevp1%wfacpbl1, vcoeff_lin_nlevp1%kpbl1,          & !in
          &           vcoeff_lin_nlevp1%wfacpbl2, vcoeff_lin_nlevp1%kpbl2,          & !in
          &           l_loglin=l_loglin,                                            & !in
          &           l_extrapol=l_extrapol, l_pd_limit=l_pd_limit,                 & !in
          &           lower_limit=lower_limit )                                       !in
        !
      CASE (VINTP_METHOD_QV )
        IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_QV")
        IF (.NOT. ASSOCIATED(vcoeff_lin)) CALL finish(routine, "Internal error!")
        IF (.NOT. ASSOCIATED(vcoeff_cub)) CALL finish(routine, "Internal error!")
        CALL qv_intp(tmp_var(:,:,:),                                                & !in
          &          out_var%r_ptr(:,:,:,out_var_idx,1),                            & !out
          &          in_z_mc, in_z3d, p_diag%temp,                                  & !in
          &          p_diag%pres, p_temp, p_pres,                                   & !in
          &          nblks, npromz, nlev, n_ipzlev,                                 & !in
          &          vcoeff_cub%coef1, vcoeff_cub%coef2,                            & !in
          &          vcoeff_cub%coef3,                                              & !in
          &          vcoeff_lin%wfac_lin, vcoeff_cub%idx0_cub,                      & !in
          &          vcoeff_lin%idx0_lin,                                           & !in
          &          vcoeff_cub%bot_idx_cub, vcoeff_lin%bot_idx_lin,                & !in
          &          vcoeff_lin%wfacpbl1, vcoeff_lin%kpbl1,                         & !in
          &          vcoeff_lin%wfacpbl2, vcoeff_lin%kpbl2,                         & !in
          &          l_satlimit=l_satlimit, lower_limit=lower_limit,                & !in
          &          l_restore_pbldev=l_restore_pbldev )                              !in
      END SELECT ! vert_intp_method

    END IF

    ! clean up
    DEALLOCATE(tmp_var, STAT=ierrstat)
    IF (ierrstat /= SUCCESS)  CALL finish (routine, 'deallocation failed')
    IF (p_info%hgrid == GRID_UNSTRUCTURED_EDGE) THEN
      DEALLOCATE(p_z3d_edge, z_me, STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'DEALLOCATE failed')
    END IF

  END SUBROUTINE pp_task_ipzlev


  !---------------------------------------------------------------
  !> Performs interpolation of a 2D field onto mean sea level, z=0.
  !
  !  This routine is completely independent from the data structures
  !  used for pz-level interpolation.
  !
  SUBROUTINE pp_task_intp_msl(ptr_task)
    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables    
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_pp_tasks:pp_task_intp_msl")
    INTEGER,  PARAMETER :: nzlev         =        1     ! just a single z-level... 
    REAL(wp), PARAMETER :: ZERO_HEIGHT   =    0._wp, &
      &                    EXTRAPOL_DIST = -500._wp

    INTEGER                            :: nblks_c, npromz_c, nblks_e, jg,          &
      &                                   in_var_idx, out_var_idx, nlev, i_endblk
    TYPE(t_nh_pzlev_config),   POINTER :: nh_pzlev_config
    TYPE(t_vcoeff)                     :: vcoeff
    TYPE (t_var_list_element), POINTER :: in_var, out_var
    TYPE(t_var_metadata),      POINTER :: p_info
    TYPE(t_patch),             POINTER :: p_patch
    TYPE(t_nh_metrics),        POINTER :: p_metrics    
    TYPE(t_nh_prog),           POINTER :: p_prog
    TYPE(t_nh_diag),           POINTER :: p_diag

    REAL(wp) :: pmsl_aux(nproma,1,ptr_task%data_input%p_patch%nblks_c), &
                pmsl_avg(nproma,1,ptr_task%data_input%p_patch%nblks_c)

    ! patch, state, and metrics
    jg             =  ptr_task%data_input%jg
    p_patch        => ptr_task%data_input%p_patch
    p_metrics      => ptr_task%data_input%p_nh_state%metrics
    p_prog         => ptr_task%data_input%p_nh_state%prog(nnow(jg))
    p_diag         => ptr_task%data_input%p_nh_state%diag

    ! pz-level interpolation data
    nh_pzlev_config   => ptr_task%data_input%nh_pzlev_config

    ! input/output field for this task
    p_info            => ptr_task%data_input%var%info
    in_var            => ptr_task%data_input%var
    out_var           => ptr_task%data_output%var

    IF (TRIM(p_info%name) /= "pres")  CALL message(routine, "Invalid input field!")

    nlev     = p_patch%nlev
    nblks_c  = p_patch%nblks_c
    npromz_c = p_patch%npromz_c
    nblks_e  = p_patch%nblks_e
    
    in_var_idx  = 1
    IF (in_var%info%lcontained)  in_var_idx  = in_var%info%ncontained
    out_var_idx = 1
    IF (out_var%info%lcontained) out_var_idx = out_var%info%ncontained

    SELECT CASE (itype_pres_msl)
    CASE (PRES_MSL_METHOD_SAI) ! stepwise analytical integration 

      IF (dbg_level >= 10)  CALL message(routine, "PRES_MSL_METHOD_SAI: stepwise analytical integration")
      ! allocate coefficient table:
      CALL vcoeff_allocate(nblks_c, nblks_e, NZLEV, vcoeff)
      ! compute extrapolation coefficients:
      CALL prepare_extrap(p_metrics%z_mc,                                     & !in
        &                 nblks_c, npromz_c, nlev,                            & !in
        &                 vcoeff%lin_cell%kpbl1, vcoeff%lin_cell%wfacpbl1,    & !out
        &                 vcoeff%lin_cell%kpbl2, vcoeff%lin_cell%wfacpbl2   )   !out
      ! Interpolate pressure on z-level "0": 
      CALL diagnose_pmsl(p_diag%pres, p_diag%tempv, p_metrics%z_mc,           &
        &                pmsl_aux(:,1,:),                                     &
        &                nblks_c, npromz_c, p_patch%nlev,                       &
        &                vcoeff%lin_cell%wfacpbl1, vcoeff%lin_cell%kpbl1,     &
        &                vcoeff%lin_cell%wfacpbl2, vcoeff%lin_cell%kpbl2,     &
        &                ZERO_HEIGHT, EXTRAPOL_DIST)
      ! deallocate coefficient tables:
      CALL vcoeff_deallocate(vcoeff)

    CASE (PRES_MSL_METHOD_GME) ! GME-type extrapolation

      IF (dbg_level >= 10)  CALL message(routine, "PRES_MSL_METHOD_GME")
      ! Interpolate pressure on z-level "0":
      CALL diagnose_pmsl_gme(p_diag%pres, p_diag%pres_sfc, p_diag%temp, &  ! in
        &                    p_metrics%z_ifc,                           &  ! in
        &                    pmsl_aux(:,1,:),                           &  ! out
        &                    nblks_c, npromz_c, p_patch%nlev )             ! in

    CASE DEFAULT
      CALL finish(routine, 'Internal error!')
    END SELECT

    IF (jg > 1) THEN ! copy outermost nest boundary row in order to avoid missing values
      i_endblk = p_patch%cells%end_blk(1,1)
      pmsl_avg(:,1,1:i_endblk) = pmsl_aux(:,1,1:i_endblk)
    ENDIF

    CALL cell_avg(pmsl_aux, p_patch, p_int_state(jg)%c_bln_avg, pmsl_avg)
    out_var%r_ptr(:,:,out_var_idx,1,1) = pmsl_avg(:,1,:)

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


  !---------------------------------------------------------------
  !> Performs interpolation of a edge-based variable onto cell centers.
  !
  !  This is only a wrapper for the corresponding routines from the
  !  interpolation module.
  SUBROUTINE pp_task_edge2cell(ptr_task)
    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_pp_tasks:pp_task_edge2cell")
    INTEGER :: &
      &  in_var_idx, out_var_idx_1, out_var_idx_2
    TYPE (t_var_list_element), POINTER :: in_var, out_var_1, out_var_2
    TYPE (t_var_metadata),     POINTER :: p_info
    TYPE(t_patch),             POINTER :: p_patch
    TYPE(t_int_state),         POINTER :: intp_hrz

    p_patch        => ptr_task%data_input%p_patch      ! patch
    intp_hrz       => ptr_task%data_input%p_int_state
    in_var         => ptr_task%data_input%var
    p_info         => ptr_task%data_input%var%info

    ! Consistency check: We make the following assumptions:
    ! - This is a 3D variable
    IF (is_2d_field(p_info%vgrid))  CALL finish(routine, "Internal error!")
    ! - We have two output components:
    IF (.NOT. ASSOCIATED(ptr_task%data_output%var) .OR.  &
      & .NOT. ASSOCIATED(ptr_task%data_output%var_2)) THEN
      CALL finish(routine, "Internal error!")
    END IF

    in_var_idx        = 1
    IF (in_var%info%lcontained)  in_var_idx  = in_var%info%ncontained

    out_var_1      => ptr_task%data_output%var
    out_var_idx_1  = 1
    IF (out_var_1%info%lcontained) out_var_idx_1 = out_var_1%info%ncontained
    out_var_2      => ptr_task%data_output%var_2
    out_var_idx_2  =  1
    IF (out_var_2%info%lcontained) out_var_idx_2 = out_var_2%info%ncontained

    CALL rbf_vec_interpol_cell(in_var%r_ptr(:,:,:,in_var_idx,1),        &   !< normal wind comp.
      &                        p_patch, intp_hrz,                       &   !< patch, interpolation state
      &                        out_var_1%r_ptr(:,:,:,out_var_idx_1,1),  &
      &                        out_var_2%r_ptr(:,:,:,out_var_idx_2,1) )     !< reconstr. u,v wind

  END SUBROUTINE pp_task_edge2cell

END MODULE mo_pp_tasks
