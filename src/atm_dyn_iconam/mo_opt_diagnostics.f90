!>
!! Definition of optional (diagnostic) model variables
!!
!! In the metadata of each model variable ("add_var") it is specified
!! *how* certain post-processing tasks, e.g. vertical interpolation
!! onto p/z-levels, are treated. In the namelist, users can then
!! specify *if* computations for these variables are performed.
!!
!! If so, the resulting model fields are appended to the list of
!! internal post-processing tasks (each field forms its own task). As
!! we do not know in advance the contents of this list, we call them
!! "optional diagnostics".
!!
!! @author F. Prill (DWD)
!!
!! @par Revision History
!! Initial implementation by F. Prill, DWD (2012-03-07)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_opt_diagnostics

  USE mo_kind,                 ONLY: wp
  USE mo_parallel_config,      ONLY: nproma
  USE mo_linked_list,          ONLY: t_var_list
  USE mo_model_domain,         ONLY: t_patch, t_subset_range
  USE mo_nonhydro_types,       ONLY: t_nh_diag,t_nh_prog
  USE mo_echam_phy_memory,     ONLY: prm_field, prm_tend
  USE mo_impl_constants,       ONLY: SUCCESS, MAX_CHAR_LENGTH,           &
    &                                VINTP_METHOD_QV,                    &
    &                                VINTP_METHOD_PRES,                  &
    &                                VINTP_METHOD_LIN,                   &
    &                                VINTP_METHOD_LIN_NLEVP1
  USE mo_exception,            ONLY: finish!!$, message, message_text
  USE mo_grid_config,          ONLY: n_dom
  USE mo_run_config,           ONLY: ntracer,iqv,iqc,iqi
  USE mo_advection_config,     ONLY: t_advection_config, advection_config
  USE mo_cdi,                  ONLY: DATATYPE_FLT32, DATATYPE_PACK16,                &
    &                                DATATYPE_PACK24, TSTEP_INSTANT,                 &
    &                                DATATYPE_FLT64
  USE mo_cdi_constants,        ONLY: GRID_UNSTRUCTURED_CELL, GRID_REFERENCE,          &
    &                                GRID_CELL, ZA_HYBRID, ZA_HYBRID_HALF, ZA_SURFACE, &
    &                                ZA_MEANSEA
  USE mo_var_list,             ONLY: default_var_list_settings,     &
    &                                new_var_list, delete_var_list, add_var, add_ref
  USE mo_var_list_element,     ONLY: level_type_ml, level_type_pl,  &
    &                                level_type_hl, level_type_il
  USE mo_name_list_output_config,ONLY: first_output_name_list, is_variable_in_output
  USE mo_io_config,            ONLY: lnetcdf_flt64_output
  USE mo_gribout_config,       ONLY: gribout_config
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_grib2,                ONLY: t_grib2_var, grib2_var
  USE mo_var_metadata,         ONLY: create_vert_interp_metadata,            &
    &                                groups, vintp_types
  USE mo_tracer_metadata,      ONLY: create_tracer_metadata
  USE mo_statistics,           ONLY: add_fields
  USE mo_util_dbg_prnt,        ONLY: dbg_print

  IMPLICIT NONE

  PRIVATE


  ! data types
  PUBLIC :: t_nh_opt_diag         ! optional diagnostic variables (data type)
  PUBLIC :: t_nh_acc
  PUBLIC :: p_nh_opt_diag         ! state vector of optional diagnostic variables
                                  ! e.g. variables on p- and/or z-levels
  PUBLIC :: t_nh_diag_pz
  PUBLIC :: t_vcoeff, t_vcoeff_lin, t_vcoeff_cub
  ! subroutines
  PUBLIC :: vcoeff_allocate, vcoeff_deallocate
  PUBLIC :: construct_opt_diag
  PUBLIC :: destruct_opt_diag
  PUBLIC :: update_opt_acc, reset_opt_acc, calc_mean_opt_acc

  ! Sub-type of "t_vcoeff" containing linear interpolation
  ! coefficients
  TYPE t_vcoeff_lin
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: wfac_lin            ! (nproma, nlev, nblks)
    INTEGER,  ALLOCATABLE, DIMENSION(:,:,:) :: idx0_lin            ! (nproma, nlev, nblks)
    INTEGER,  ALLOCATABLE, DIMENSION(:,:)   :: bot_idx_lin         ! (nproma, nblks)
    INTEGER,  ALLOCATABLE, DIMENSION(:,:)   :: kpbl1, kpbl2        ! (nproma, nblks)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: wfacpbl1, wfacpbl2  ! (nproma, nblks)
  END TYPE t_vcoeff_lin

  ! Sub-type of "t_vcoeff" containing cubic interpolation
  ! coefficients
  TYPE t_vcoeff_cub
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: coef1, coef2, coef3 ! (nproma, nlev, nblks)
    INTEGER,  ALLOCATABLE, DIMENSION(:,:,:) :: idx0_cub            ! (nproma, nlev, nblks)
    INTEGER,  ALLOCATABLE, DIMENSION(:,:)   :: bot_idx_cub         ! (nproma, nblks)
  END TYPE t_vcoeff_cub

  ! Derived type containing coefficient tables for vertical
  ! interpolation.
  TYPE t_vcoeff
    LOGICAL :: &
      & l_initialized = .FALSE.,  &
      & l_allocated   = .FALSE.

    ! LINEAR interpolation data
    TYPE (t_vcoeff_lin) ::    &
      &    lin_cell,          &  !< cell centers: interpolation data for the model levels
      &    lin_cell_nlevp1,   &  !< cell centers: interpolation data for the vertical interface of cells, "nlevp1"
      &    lin_edge              !< edge midpts:  interpolation data for the model levels

    ! CUBIC interpolation (model levels)
    TYPE (t_vcoeff_cub) ::    &
      &    cub_cell,          &  !< cell centers: interpolation data for the model levels
      &    cub_edge

  END TYPE t_vcoeff

  TYPE t_pointer_3d_wp
    REAL(wp),POINTER :: p(:,:,:)  ! pointer to 3D array
  END TYPE t_pointer_3d_wp
  ! variable to be accumulated manually
  TYPE t_nh_acc
    REAL(wp), POINTER   &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    , CONTIGUOUS        &
#endif
    &  ::               &
    !
    ! dynamics
    &  rho(:,:,:),      &
    &  qv(:,:,:),       &
    &  qc(:,:,:),       &
    &  qi(:,:,:),       &
    &  temp(:,:,:),     &
    &  pres_sfc(:,:),   &
    &  pres_msl(:,:),   &
    &  pres(:,:,:),     &
    &  pres_ifc(:,:,:), &
    &  u(:,:,:),        &
    &  v(:,:,:),        &
    &  w(:,:,:),        &
    &  omega(:,:,:),    &
    !
    ! tracers container
    &  tracer(:,:,:,:), &

    !
    ! echam physics
    &  cosmu0(:,:),     &
    &  rsdt(:,:),       &
    &  relhum(:,:,:),   &
    &  aclc(:,:,:),     &
    &  aclcov(:,:),     &
    &  rsfl(:,:),       &
    &  rsfc(:,:),       &
    &  ssfl(:,:),       &
    &  ssfc(:,:),       &
    &  totprec(:,:),    &
    &  qvi(:,:),        &
    &  xlvi(:,:),       &
    &  xivi(:,:),       &
    &  swflxsfc(:,:),   &
    &  swflxtoa(:,:),   &
    &  lwflxsfc(:,:),   &
    &  lwflxtoa(:,:),   &
    &  tsfc(:,:),       &
    &  evap(:,:),       &
    &  lhflx(:,:),      &
    &  shflx(:,:),      &
    &  u_stress(:,:),   &
    &  v_stress(:,:),   &
    &  u_stress_sso(:,:),    &
    &  v_stress_sso(:,:),    &
    &  dissipation_sso(:,:), &
    &  seaice(:,:),     &
    &  siced(:,:),      &
    &  albedo(:,:),     &
    &  sfcWind(:,:),    &
    &  uas(:,:),        &
    &  vas(:,:),        &
    &  tas(:,:),        &
    &  dew2(:,:),       &
    !
    ! tendencies
    ! - temperature:
    &  tend_ta(:,:,:)     ,&
    &  tend_ta_dyn(:,:,:) ,&
    &  tend_ta_phy(:,:,:) ,&
    &  tend_ta_rlw(:,:,:) ,&
    &  tend_ta_rlw_impl(:,:),&
    &  tend_ta_rsw(:,:,:) ,&
    &  tend_ta_cld(:,:,:) ,&
    &  tend_ta_cnv(:,:,:) ,&
    &  tend_ta_vdf(:,:,:) ,&
    &  tend_ta_gwh(:,:,:) ,&
    &  tend_ta_sso(:,:,:) ,&
    !
    !  - u-wind:
    &  tend_ua(:,:,:)     ,&
    &  tend_ua_dyn(:,:,:) ,&
    &  tend_ua_phy(:,:,:) ,&
    &  tend_ua_cnv(:,:,:) ,&
    &  tend_ua_vdf(:,:,:) ,&
    &  tend_ua_gwh(:,:,:) ,&
    &  tend_ua_sso(:,:,:) ,&
    !
    !  - v-wind:
    &  tend_va(:,:,:)     ,&
    &  tend_va_dyn(:,:,:) ,&
    &  tend_va_phy(:,:,:) ,&
    &  tend_va_cnv(:,:,:) ,&
    &  tend_va_vdf(:,:,:) ,&
    &  tend_va_gwh(:,:,:) ,&
    &  tend_va_sso(:,:,:) !!$,&
!!$    !
!!$    !  - specific humidity
!!$    &  tend_hus(:,:,:)    ,&
!!$    &  tend_hus_dyn(:,:,:),&
!!$    &  tend_hus_phy(:,:,:),&
!!$    &  tend_hus_cld(:,:,:),&
!!$    &  tend_hus_cnv(:,:,:),&
!!$    &  tend_hus_vdf(:,:,:),&
!!$    !
!!$    !  - xl and xi
!!$    &  tend_clw_dtr(:,:,:),&
!!$    &  tend_cli_dtr(:,:,:)

    TYPE(t_pointer_3d_wp),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer

    ! Internal counter for accumulation operations
    INTEGER :: numberOfAccumulations

    ! logicals for presence of time mean output variables in the output name lists
    !
    !  inidcate if any time averaged variable is requested for the output
    LOGICAL :: l_any_m
    !
    !  dynamics
    LOGICAL :: l_ua_m
    LOGICAL :: l_va_m
    LOGICAL :: l_wa_m
    LOGICAL :: l_rho_m
    LOGICAL :: l_ta_m
    LOGICAL :: l_ps_m
    LOGICAL :: l_psl_m
    LOGICAL :: l_pfull_m
    LOGICAL :: l_phalf_m
    LOGICAL :: l_wap_m
    !
    !  tracers
    LOGICAL :: l_tracer_m
    !
    !  physics
    LOGICAL :: l_cosmu0_m
    LOGICAL :: l_rsdt_m
    LOGICAL :: l_hur_m
    LOGICAL :: l_cl_m
    LOGICAL :: l_clt_m
    LOGICAL :: l_prlr_m
    LOGICAL :: l_prls_m
    LOGICAL :: l_prcr_m
    LOGICAL :: l_prcs_m
    LOGICAL :: l_pr_m
    LOGICAL :: l_prw_m
    LOGICAL :: l_cllvi_m
    LOGICAL :: l_clivi_m
    LOGICAL :: l_rsns_m
    LOGICAL :: l_rsnt_m
    LOGICAL :: l_rlns_m
    LOGICAL :: l_rlnt_m
    LOGICAL :: l_ts_m
    LOGICAL :: l_evspsbl_m
    LOGICAL :: l_hfls_m
    LOGICAL :: l_hfss_m
    LOGICAL :: l_tauu_m
    LOGICAL :: l_tauv_m
    LOGICAL :: l_tauu_sso_m
    LOGICAL :: l_tauv_sso_m
    LOGICAL :: l_diss_sso_m
    LOGICAL :: l_sic_m
    LOGICAL :: l_sit_m
    LOGICAL :: l_albedo_m
    LOGICAL :: l_sfcWind_m
    LOGICAL :: l_uas_m
    LOGICAL :: l_vas_m
    LOGICAL :: l_tas_m
    LOGICAL :: l_dew2_m
    !
    !  tendencies
    !  of temperature:
    LOGICAL :: l_tend_ta_m
    LOGICAL :: l_tend_ta_dyn_m
    LOGICAL :: l_tend_ta_phy_m
    LOGICAL :: l_tend_ta_rsw_m
    LOGICAL :: l_tend_ta_rlw_m
    LOGICAL :: l_tend_ta_rlw_impl_m
    LOGICAL :: l_tend_ta_cld_m
    LOGICAL :: l_tend_ta_cnv_m
    LOGICAL :: l_tend_ta_vdf_m
    LOGICAL :: l_tend_ta_gwh_m
    LOGICAL :: l_tend_ta_sso_m
    !
    !  of u-wind:
    LOGICAL :: l_tend_ua_m
    LOGICAL :: l_tend_ua_dyn_m
    LOGICAL :: l_tend_ua_phy_m
    LOGICAL :: l_tend_ua_cnv_m
    LOGICAL :: l_tend_ua_vdf_m
    LOGICAL :: l_tend_ua_gwh_m
    LOGICAL :: l_tend_ua_sso_m
    !
    !  of v-wind:
    LOGICAL :: l_tend_va_m
    LOGICAL :: l_tend_va_dyn_m
    LOGICAL :: l_tend_va_phy_m
    LOGICAL :: l_tend_va_cnv_m
    LOGICAL :: l_tend_va_vdf_m
    LOGICAL :: l_tend_va_gwh_m
    LOGICAL :: l_tend_va_sso_m
    !
!!$    !  of specific humidity
!!$    LOGICAL :: l_tend_hus_m
!!$    LOGICAL :: l_tend_hus_dyn_m
!!$    LOGICAL :: l_tend_hus_phy_m
!!$    LOGICAL :: l_tend_hus_cld_m
!!$    LOGICAL :: l_tend_hus_cnv_m
!!$    LOGICAL :: l_tend_hus_vdf_m
!!$    !
!!$    !  of xl and xi
!!$    LOGICAL :: l_tend_clw_dtr_m
!!$    LOGICAL :: l_tend_cli_dtr_m

  END TYPE t_nh_acc


  ! State vector for diagnostic variables on p-, z- and/or i-levels
  !
  ! @note The pointers which are collected in this derived type
  !       constitute only the minimum set of fields that are required
  !       for i/p/z-level interpolation. All other variables are
  !       stored inside the "opt_diag_list_p", "opt_diag_list_z"
  !       variable lists.
  TYPE t_nh_diag_pz

    REAL(wp), POINTER ::    &
      !--- cells (nproma,nlev,nblks)
      ! fields that are essential for z-level interpolation:
      &  z_temp(:,:,:),        & ! temperature                  [K]
      &  z_pres(:,:,:),        & ! pressure                     [Pa]
      ! fields that are essential for p-level interpolation only:
      &  p_gh    (:,:,:),      & ! geopotential height          [m]
      &  p_temp(:,:,:),        & ! temperature                  [K]
      ! fields that are essential for interpolation on isentropes only:
      &  i_gh    (:,:,:),      & ! geopotential height          [m]
      &  i_temp(:,:,:)           ! temperature                  [K]

    ! coefficient tables for vertical interpolation. There exist
    ! different kinds of coefficients: For p-, z-,and for
    ! i-level-interpolation; for cells and for edges.
    TYPE(t_vcoeff) :: vcoeff_z, vcoeff_p, vcoeff_i

  END TYPE t_nh_diag_pz


  ! List of optional diagnostics + necessary meta data
  TYPE t_nh_opt_diag

    ! diag_pz: data structure containing coefficient tables and
    ! pointers to a few number of fields which are required for
    ! interpolation of model variables to p/z-levels
    TYPE(t_nh_diag_pz) :: diag_pz

    TYPE(t_nh_acc)     :: acc

    ! opt_diag_list: List of optional diagnostics variables.
    !
    ! The "opt_diag_list_*" lists contain all variables that have been
    ! interpolated onto p/z-levels
    TYPE(t_var_list)   :: opt_diag_list,   opt_diag_list_p, &
      &                   opt_diag_list_z, opt_diag_list_i, &
      &                   opt_acc_list

  END TYPE t_nh_opt_diag


  ! Actual instantiation of optional diagnostics type "t_nh_opt_diag"
  TYPE(t_nh_opt_diag), TARGET, ALLOCATABLE :: p_nh_opt_diag(:)


CONTAINS

  ! setup of accumulation variables
  SUBROUTINE construct_opt_acc(p_patch,list,p_acc)
    TYPE(t_patch),        INTENT(IN) :: p_patch
    TYPE(t_var_list)                 :: list
    TYPE(t_nh_acc)                   :: p_acc

    ! LOCAL ===================================================================
    INTEGER :: nblks_c       !< number of cell blocks to allocate
!!$    INTEGER :: nblks_e       !< number of edge blocks to allocate
!!$    INTEGER :: nblks_v       !< number of vertex blocks to allocate

    INTEGER :: nlev
    INTEGER :: nlevp1

    INTEGER :: jt

    INTEGER :: shape2d  (2)
    INTEGER :: shape2d_c(2), shape3d_c(3), shape3d_chalf(3), shape4d_c(4)
!!$    INTEGER :: shape2d_e(2), shape3d_e(3)
!!$    INTEGER ::               shape3d_v(3)

    INTEGER :: ibits,iextbits     !< "entropy" of horizontal slice
    INTEGER :: DATATYPE_PACK_VAR  !< variable "entropy" for some thermodynamic fields
    INTEGER :: datatype_flt       !< floating point accuracy in NetCDF output

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc
    TYPE(t_advection_config), POINTER :: advconf
    ! =========================================================================

    !determine size of arrays
    nblks_c = p_patch%nblks_c
!!$    nblks_e = p_patch%nblks_e
!!$    nblks_v = p_patch%nblks_v

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice
    iextbits = DATATYPE_PACK24

    IF (gribout_config(p_patch%id)%lgribout_24bit) THEN  ! analysis
      ! higher accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK24
    ELSE
      ! standard accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK16
    ENDIF

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ! pointer to advection_config(jg) to save some paperwork
    advconf => advection_config(p_patch%id)

    ! predefined array shapes
    shape2d_c     = (/nproma,          nblks_c    /)
    shape2d       = shape2d_c
    shape3d_c     = (/nproma, nlev   , nblks_c    /)
    shape3d_chalf = (/nproma, nlevp1 , nblks_c    /)
    shape4d_c     = (/nproma, nlev   , nblks_c, ntracer     /)
!!$    shape2d_e     = (/nproma,          nblks_e    /)
!!$    shape3d_e     = (/nproma, nlev   , nblks_e    /)
!!$    shape3d_v     = (/nproma, nlev   , nblks_v    /)

    p_acc%l_any_m = .FALSE.

    ! PROGS {{{
    p_acc%l_ua_m  = is_variable_in_output(first_output_name_list, var_name="ua_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_ua_m
    IF (p_acc%l_ua_m) THEN
       cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'Zonal wind (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 2, 2, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'ua_m', p_acc%u,                                        &
                   & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
                   & ldims=shape3d_c,                                              &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &   vert_intp_type=vintp_types("P","Z","I") ),                  &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_va_m  = is_variable_in_output(first_output_name_list, var_name="va_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_va_m
    IF (p_acc%l_va_m) THEN
       cf_desc    = t_cf_var('northward_wind', 'm s-1', 'Meridional wind (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 2, 3, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'va_m', p_acc%v,                                        &
                   & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
                   & ldims=shape3d_c,                                              &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &   vert_intp_type=vintp_types("P","Z","I") ),                  &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_wa_m  = is_variable_in_output(first_output_name_list, var_name="wa_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_wa_m
    IF (p_acc%l_wa_m) THEN
       cf_desc    = t_cf_var('upward_air_velocity', 'm s-1', 'Vertical velocity (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 2, 9, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'wa_m', p_acc%w,                                        &
                   & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,  &
                   & ldims=shape3d_chalf,                                          &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &   vert_intp_type=vintp_types("P","Z","I"),                    &
                   &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),                 &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_rho_m = is_variable_in_output(first_output_name_list, var_name="rho_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_rho_m
    IF (p_acc%l_rho_m) THEN
       cf_desc    = t_cf_var('air_density', 'kg m-3', 'density (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 3, 10, DATATYPE_PACK_VAR, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'rho_m', p_acc%rho,                                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
                   & ldims=shape3d_c,                                              &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &   vert_intp_type=vintp_types("P","Z","I"),                    &
                   &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_ta_m  = is_variable_in_output(first_output_name_list, var_name="ta_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_ta_m
    IF (p_acc%l_ta_m) THEN
       cf_desc    = t_cf_var('air temperature', 'K', 'Temperature', datatype_flt)
       grib2_desc = grib2_var(0, 0, 0, DATATYPE_PACK_VAR, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'ta_m', p_acc%temp,                                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
                   & ldims=shape3d_c,                                              &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &             vert_intp_type=vintp_types("P","Z","I"),          &
                   &             vert_intp_method=VINTP_METHOD_LIN ),              &
                   & in_group=groups("prog_timemean","atmo_timemean"))
    END IF

    p_acc%l_ps_m  = is_variable_in_output(first_output_name_list, var_name="ps_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_ps_m
    IF (p_acc%l_ps_m) THEN
       cf_desc    = t_cf_var('surface_air_pressure', 'Pa', 'surface pressure (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'ps_m', p_acc%pres_sfc,                                 &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                   & ldims=shape2d_c,                                              &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_psl_m = is_variable_in_output(first_output_name_list, var_name="psl_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_psl_m
    IF (p_acc%l_psl_m) THEN
       cf_desc    = t_cf_var('mean sea level pressure', 'Pa',                      &
         &                   'mean sea level pressure (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 3, 1, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'psl_m', p_acc%pres_msl,                                &
                   & GRID_UNSTRUCTURED_CELL, ZA_MEANSEA, cf_desc, grib2_desc,      &
                   & ldims=shape2d_c,                                              &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_pfull_m = is_variable_in_output(first_output_name_list, var_name="pfull_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_pfull_m
    IF (p_acc%l_pfull_m) THEN
       cf_desc    = t_cf_var('air_pressure', 'Pa', 'pressure at full level (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 3, 0, DATATYPE_PACK_VAR, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'pfull_m', p_acc%pres,                                  &
                   & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
                   & ldims=shape3d_c, lrestart=.FALSE. ,                           &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &             vert_intp_type=vintp_types("P","Z","I"),          &
                   &             vert_intp_method=VINTP_METHOD_PRES ),             &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_phalf_m = is_variable_in_output(first_output_name_list, var_name="phalf_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_phalf_m
    IF (p_acc%l_phalf_m) THEN
       cf_desc    = t_cf_var('air_pressure', 'Pa', 'pressure at half level (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'phalf_m', p_acc%pres_ifc,                              &
                   & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,  &
                   & ldims=shape3d_chalf, lrestart=.FALSE.,                        &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &             vert_intp_type=vintp_types("P","Z","I"),          &
                   &             vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),       &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_wap_m = is_variable_in_output(first_output_name_list, var_name="wap_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_wap_m
    IF (p_acc%l_wap_m) THEN
       cf_desc    = t_cf_var('omega', 'Pa/s', 'vertical velocity (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 2, 8, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list,"wap_m", p_acc%omega,                                    &
                   & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                            &
                   & cf_desc, grib2_desc,                                          &
                   & ldims=shape3d_c,                                              &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &             vert_intp_type=vintp_types("P","Z","I"),          &
                   &             vert_intp_method=VINTP_METHOD_LIN,                &
                   &             l_loglin=.FALSE., l_extrapol=.FALSE.),            &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF
    ! }}}

    ! TRACERS {{{
    ! support qv,qc,qi because they are always there
    IF (ntracer > 0) THEN
       p_acc%l_tracer_m = is_variable_in_output(first_output_name_list, var_name="hus_m") .OR. &
                        & is_variable_in_output(first_output_name_list, var_name="clw_m") .OR. &
                        & is_variable_in_output(first_output_name_list, var_name="cli_m")
       p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tracer_m
       IF (p_acc%l_tracer_m) THEN
          cf_desc    = t_cf_var('tracer', 'kg kg-1', 'air tracer (time mean)', datatype_flt)
          grib2_desc = grib2_var(0,20,2, ibits, GRID_REFERENCE, GRID_CELL)
          CALL add_var( list, 'tracer_m', p_acc%tracer,                         &
                      & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, &
                      & ldims=shape4d_c ,                                       &
                      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

          ALLOCATE(p_acc%tracer_ptr(ntracer))
          DO jt=1,ntracer
             IF (jt == iqv ) CALL add_ref(                                          &
                  &  list, 'tracer_m', 'hus_m', p_acc%tracer_ptr(jt)%p,             &
                  &  GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &
                  &  t_cf_var('specific_humidity', 'kg kg-1',                       &
                  &           'specific_humidity (time mean)', datatype_flt),       &
                  &  grib2_var( 0, 1, 0, ibits, GRID_REFERENCE, GRID_CELL),         &
                  &  ldims=shape3d_c,                                               &
                  &  tlev_source=1,                                                 &
                  &  tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                  &              ihadv_tracer=advconf%ihadv_tracer(iqv),            &
                  &              ivadv_tracer=advconf%ivadv_tracer(iqv)),           &
                  &  vert_interp=create_vert_interp_metadata(                       &
                  &              vert_intp_type=vintp_types("P","Z","I"),           &
                  &              vert_intp_method=VINTP_METHOD_QV,                  &
                  &              l_satlimit=.FALSE.,                                &
                  &              lower_limit=2.5e-6_wp, l_restore_pbldev=.FALSE. ), &
                  &  in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
                  &                  "tracer_timemean","atmo_timemean"))

             IF ( jt == iqc )  CALL add_ref(                                        &
                  &  list, 'tracer_m', 'clw_m', p_acc%tracer_ptr(jt)%p,             &
                  &  GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &
                  &  t_cf_var('specific_cloud_water_content', 'kg kg-1',            &
                  &           'specific_cloud_water_content (time mean)',datatype_flt), &
                  &  grib2_var(0, 1, 22, ibits, GRID_REFERENCE, GRID_CELL),         &
                  &  ldims=shape3d_c,                                               &
                  &  tlev_source=1,                                                 &
                  &  tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                  &              ihadv_tracer=advconf%ihadv_tracer(iqc),            &
                  &              ivadv_tracer=advconf%ivadv_tracer(iqc)),           &
                  &  vert_interp=create_vert_interp_metadata(                       &
                  &              vert_intp_type=vintp_types("P","Z","I"),           &
                  &              vert_intp_method=VINTP_METHOD_LIN,                 &
                  &              l_loglin=.FALSE.,                                  &
                  &              l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                  &              lower_limit=0._wp  ),                              &
                  &  in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
                  &                  "tracer_timemean","atmo_timemean"))

             IF ( jt == iqi ) CALL add_ref(                                         &
                  &  list, 'tracer_m', 'cli_m', p_acc%tracer_ptr(jt)%p,             &
                  &  GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &
                  &  t_cf_var('specific_cloud_ice_content', 'kg kg-1',              &
                  &           'specific_cloud_ice_content (time mean)', datatype_flt),  &
                  &  grib2_var(0, 1, 82, ibits, GRID_REFERENCE, GRID_CELL),         &
                  &  ldims=shape3d_c,                                               &
                  &  tlev_source=1,                                                 &
                  &  tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                  &              ihadv_tracer=advconf%ihadv_tracer(iqi),            &
                  &              ivadv_tracer=advconf%ivadv_tracer(iqi)),           &
                  &  vert_interp=create_vert_interp_metadata(                       &
                  &              vert_intp_type=vintp_types("P","Z","I"),           &
                  &              vert_intp_method=VINTP_METHOD_LIN,                 &
                  &              l_loglin=.FALSE.,                                  &
                  &              l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                  &              lower_limit=0._wp  ),                              &
                  &  in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
                  &                  "tracer_timemean","atmo_timemean"))
          END DO
       END IF
    END IF
    ! }}}

    ! ECHAM {{{
    p_acc%l_cosmu0_m = is_variable_in_output(first_output_name_list, var_name="cosmu0_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_cosmu0_m
    IF (p_acc%l_cosmu0_m) THEN
       cf_desc    = t_cf_var('cosmu0', '', 'cosine of the zenith angle (time mean)', datatype_flt)
       grib2_desc = grib2_var(192,214,1, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'cosmu0_m', p_acc%cosmu0,                                  &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"))
    END IF

    p_acc%l_rsdt_m = is_variable_in_output(first_output_name_list, var_name="rsdt_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_rsdt_m
    IF (p_acc%l_rsdt_m) THEN
       cf_desc    = t_cf_var('rsdt', 'W m-2',                                                    &
                   &         'downward shortwave flux at the top of the atmosphere (time mean)', &
                   &         datatype_flt)
       grib2_desc = grib2_var(0,4,7, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'rsdt_m', p_acc%rsdt,                                &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean") )
    END IF

    p_acc%l_hur_m = is_variable_in_output(first_output_name_list, var_name="hur_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_hur_m
    IF (p_acc%l_hur_m) THEN
       cf_desc    = t_cf_var('hur', '',                       &
                   &         'relative humidity (time mean)', &
                   &         datatype_flt)
       grib2_desc = grib2_var(0,1,1, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'hur_m', p_acc%relhum,                                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                               &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape3d_c,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_cl_m  = is_variable_in_output(first_output_name_list, var_name="cl_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_cl_m
    IF (p_acc%l_cl_m) THEN
       cf_desc    = t_cf_var('cl', 'm2 m-2',                    &
                   &         'cloud area fraction (time mean)', &
                   &         datatype_flt)
       grib2_desc = grib2_var(0,6,22, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'cl_m', p_acc%aclc,                                        &
                   & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                               &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape3d_c,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_clt_m = is_variable_in_output(first_output_name_list, var_name="clt_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_clt_m
    IF (p_acc%l_clt_m) THEN
       cf_desc    = t_cf_var('clt', 'm2 m-2',                 &
                   &         'total cloud cover (time mean)', &
                   &         datatype_flt)
       grib2_desc = grib2_var(0,6,1, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'clt_m', p_acc%aclcov,                                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_prlr_m = is_variable_in_output(first_output_name_list, var_name="prlr_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_prlr_m
    IF (p_acc%l_prlr_m) THEN
       cf_desc    = t_cf_var('prlr', 'kg m-2 s-1',                                 &
                   &         'large-scale precipitation flux (water) (time mean)', &
                   &         datatype_flt)
       grib2_desc = grib2_var(0,1,77, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'prlr_m', p_acc%rsfl,                                      &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_prls_m = is_variable_in_output(first_output_name_list, var_name="prls_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_prls_m
    IF (p_acc%l_prls_m) THEN
       cf_desc    = t_cf_var('prls', 'kg m-2 s-1',                                &
                   &         'large-scale precipitation flux (snow) (time mean)', &
                   &         datatype_flt)
       grib2_desc = grib2_var(0,1,59, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'prls_m', p_acc%ssfl,                                      &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_prcr_m = is_variable_in_output(first_output_name_list, var_name="prcr_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_prcr_m
    IF (p_acc%l_prcr_m) THEN
       cf_desc    = t_cf_var('prcr', 'kg m-2 s-1',                                &
                   &         'convective precipitation flux (water) (time mean)', &
                   &         datatype_flt)
       grib2_desc = grib2_var(0,1,76, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'prcr_m', p_acc%rsfc,                                      &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_prcs_m = is_variable_in_output(first_output_name_list, var_name="prcs_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_prcs_m
    IF (p_acc%l_prcs_m) THEN
       cf_desc    = t_cf_var('prcs', 'kg m-2 s-1',                               &
                   &         'convective precipitation flux (snow) (time mean)', &
                   &         datatype_flt)
       grib2_desc = grib2_var(0,1,58, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'prcs_m', p_acc%ssfc,                                      &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_pr_m = is_variable_in_output(first_output_name_list, var_name="pr_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_pr_m
    IF (p_acc%l_pr_m) THEN
       cf_desc    = t_cf_var('pr', 'kg m-2 s-1',               &
                   &         'precipitation flux (time mean)', &
                   &         datatype_flt)
       grib2_desc = grib2_var(0, 1, 52, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'pr_m', p_acc%totprec,                                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_prw_m = is_variable_in_output(first_output_name_list, var_name="prw_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_prw_m
    IF (p_acc%l_prw_m) THEN
       cf_desc    = t_cf_var('total_vapour', 'kg m-2',                         &
                   &         'vertically integrated water vapour (time mean)', &
                   &         datatype_flt)
       grib2_desc = grib2_var(0,1,64, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'prw_m', p_acc%qvi,                                        &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_cllvi_m = is_variable_in_output(first_output_name_list, var_name="cllvi_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_cllvi_m
    IF (p_acc%l_cllvi_m) THEN
       cf_desc    = t_cf_var('total_cloud_water', 'kg m-2', &
                   & 'vertically integrated cloud water (time mean)', datatype_flt)
       grib2_desc = grib2_var(0,1,69, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'cllvi_m', p_acc%xlvi,                                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_clivi_m = is_variable_in_output(first_output_name_list, var_name="clivi_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_clivi_m
    IF (p_acc%l_clivi_m) THEN
       cf_desc    = t_cf_var('total_cloud_ice', 'kg m-2', &
                   & 'vertically integrated cloud ice (time mean)', datatype_flt)
       grib2_desc = grib2_var(0,1,70, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'clivi_m', p_acc%xivi,                                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_rsns_m = is_variable_in_output(first_output_name_list, var_name="rsns_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_rsns_m
    IF (p_acc%l_rsns_m) THEN
       cf_desc    = t_cf_var('rsns', 'W m-2', ' shortwave net flux at surface (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'rsns_m', p_acc%swflxsfc,                                  &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_rsnt_m = is_variable_in_output(first_output_name_list, var_name="rsnt_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_rsnt_m
    IF (p_acc%l_rsnt_m) THEN
       cf_desc    = t_cf_var('rsnt', 'W m-2', ' shortwave net flux at TOA (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'rsnt_m', p_acc%swflxtoa,                                  &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_rlns_m = is_variable_in_output(first_output_name_list, var_name="rlns_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_rlns_m
    IF (p_acc%l_rlns_m) THEN
       cf_desc    = t_cf_var('rlns', 'W m-2', 'longwave net flux at surface (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'rlns_m', p_acc%lwflxsfc,                                  &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_rlnt_m = is_variable_in_output(first_output_name_list, var_name="rlnt_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_rlnt_m
    IF (p_acc%l_rlnt_m) THEN
       cf_desc    = t_cf_var('rlnt', 'W m-2', 'longwave net flux at TOA (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'rlnt_m', p_acc%lwflxtoa,&
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & cf_desc, grib2_desc,                                             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_ts_m = is_variable_in_output(first_output_name_list, var_name="ts_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_ts_m
    IF (p_acc%l_ts_m) THEN
       cf_desc    = t_cf_var('surface_temperature', '', 'surface temperature (time mean)', datatype_flt)
       grib2_desc = grib2_var(0,0,0, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'ts_m', p_acc%tsfc,                                        &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"))
    END IF

    p_acc%l_evspsbl_m = is_variable_in_output(first_output_name_list, var_name="evspsbl_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_evspsbl_m
    IF (p_acc%l_evspsbl_m) THEN
       CALL add_var( list, 'evspsbl_m', p_acc%evap,                                   &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & t_cf_var('evap', 'kg m-2 s-1', 'evaporation (time mean)',        &
                   & datatype_flt),                                                   &
                   & grib2_var(0,1,6,iextbits, GRID_REFERENCE, GRID_CELL),            &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT                                 )
    END IF

    p_acc%l_hfls_m = is_variable_in_output(first_output_name_list, var_name="hfls_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_hfls_m
    IF (p_acc%l_hfls_m) THEN
       CALL add_var( list, 'hfls_m', p_acc%lhflx,                                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & t_cf_var('hfls', 'W m-2 ', 'latent heat flux (time mean)',       &
                   & datatype_flt),                                                   &
                   & grib2_var(0,0,10, ibits, GRID_REFERENCE, GRID_CELL),             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT                                 )
    END IF

    p_acc%l_hfss_m = is_variable_in_output(first_output_name_list, var_name="hfss_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_hfss_m
    IF (p_acc%l_hfss_m) THEN
       CALL add_var( list, 'hfss_m', p_acc%shflx,                                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                   & t_cf_var('hfss', 'W m-2 ', 'sensible heat flux (time mean)',     &
                   & datatype_flt),                                                   &
                   & grib2_var(0,0,11, ibits, GRID_REFERENCE, GRID_CELL),             &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
                   & isteptype=TSTEP_INSTANT                                 )
    END IF

    p_acc%l_tauu_m = is_variable_in_output(first_output_name_list, var_name="tauu_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tauu_m
    IF (p_acc%l_tauu_m) THEN
       CALL add_var( list, 'tauu_m', p_acc%u_stress,                                             &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                         &
                   & t_cf_var('u_stress', 'N m-2', 'u-momentum flux at the surface (time mean)', &
                   &          datatype_flt),                                                     &
                   & grib2_var(0,2,17, ibits, GRID_REFERENCE, GRID_CELL),                        &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),            &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_tauv_m = is_variable_in_output(first_output_name_list, var_name="tauv_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tauv_m
    IF (p_acc%l_tauv_m) THEN
       CALL add_var( list, 'tauv_m', p_acc%v_stress,                                             &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                         &
                   & t_cf_var('v_stress', 'N m-2', 'v-momentum flux at the surface (time mean)', &
                   &          datatype_flt),                                                     &
                   & grib2_var(0,2,18, ibits, GRID_REFERENCE, GRID_CELL),                        &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),            &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_tauu_sso_m = is_variable_in_output(first_output_name_list, var_name="tauu_sso_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tauu_sso_m
    IF (p_acc%l_tauu_sso_m) THEN
       CALL add_var( list, 'tauu_sso_m', p_acc%u_stress_sso,                                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                         &
                   & t_cf_var('u_stress', 'N m-2',                                               &
                   &          'zonal stress from subgrid scale orographic drag (time mean)',     &
                   &          datatype_flt),                                                     &
                   & grib2_var(0,2,17, ibits, GRID_REFERENCE, GRID_CELL),                        &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),            &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_tauv_sso_m = is_variable_in_output(first_output_name_list, var_name="tauv_sso_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tauv_sso_m
    IF (p_acc%l_tauv_sso_m) THEN
       CALL add_var( list, 'tauv_sso_m', p_acc%v_stress_sso,                                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                         &
                   & t_cf_var('v_stress', 'N m-2',                                               &
                   &          'meridional stress from subgrid scale orographic drag (time mean)',&
                   &          datatype_flt),                                                     &
                   & grib2_var(0,2,18, ibits, GRID_REFERENCE, GRID_CELL),                        &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),            &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_diss_sso_m = is_variable_in_output(first_output_name_list, var_name="diss_sso_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_diss_sso_m
    IF (p_acc%l_diss_sso_m) THEN
       CALL add_var( list, 'diss_sso_m', p_acc%dissipation_sso,                                  &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                         &
                   & t_cf_var('dissipation_sso', '',                                             &
                   &          'dissipation of orographic waves (time mean)',                     &
                   &          datatype_flt),                                                     &
                   & grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL),                       &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),            &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_sic_m = is_variable_in_output(first_output_name_list, var_name="sic_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_sic_m
    IF (p_acc%l_sic_m) THEN
       CALL add_var( list, 'sic_m', p_acc%seaice,                                                &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                         &
                   & t_cf_var('sea_ice_cover', '',                                               &
                   &          'fraction of ocean covered by sea ice (time mean)',                &
                   &          datatype_flt),                                                     &
                   & grib2_var(10,2,0, ibits, GRID_REFERENCE, GRID_CELL),                        &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),            &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_sit_m = is_variable_in_output(first_output_name_list, var_name="sit_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_sit_m
    IF (p_acc%l_sit_m) THEN
       CALL add_var( list, 'sit_m', p_acc%siced,                                                 &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                         &
                   & t_cf_var('sea_ice_thickness', 'm',                                          &
                   &          'sea ice thickness (time mean)',                                   &
                   &          datatype_flt),                                                     &
                   & grib2_var(10,2,1, ibits, GRID_REFERENCE, GRID_CELL),                        &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),            &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_albedo_m = is_variable_in_output(first_output_name_list, var_name="albedo_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_albedo_m
    IF (p_acc%l_albedo_m) THEN
       CALL add_var( list, 'albedo_m', p_acc%albedo,                                             &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                         &
                   & t_cf_var('albedo', '',                                                      &
                   &          'surface albedo (time mean)',                                      &
                   &          datatype_flt),                                                     &
                   & grib2_var(0,19,1, ibits, GRID_REFERENCE, GRID_CELL),                        &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),            &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_sfcWind_m = is_variable_in_output(first_output_name_list, var_name="sfcWind_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_sfcWind_m
    IF (p_acc%l_sfcWind_m) THEN
       CALL add_var( list, 'sfcWind_m', p_acc%sfcWind,                                           &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                         &
                   & t_cf_var('sfcWind', 'm s-1',                                                &
                   &          '10m windspeed (time mean)',                                       &
                   &          datatype_flt),                                                     &
                   & grib2_var(0,2,1, ibits, GRID_REFERENCE, GRID_CELL),                         &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),            &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_uas_m = is_variable_in_output(first_output_name_list, var_name="uas_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_uas_m
    IF (p_acc%l_uas_m) THEN
       CALL add_var( list, 'uas_m', p_acc%uas,                                                   &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                         &
                   & t_cf_var('uas', 'm s-1',                                                    &
                   &          'zonal wind in 10m (time mean)',                                   &
                   &          datatype_flt),                                                     &
                   & grib2_var(0,2,2, ibits, GRID_REFERENCE, GRID_CELL),                         &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),            &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_vas_m = is_variable_in_output(first_output_name_list, var_name="vas_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_vas_m
    IF (p_acc%l_vas_m) THEN
       CALL add_var( list, 'vas_m', p_acc%vas,                                                   &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                         &
                   & t_cf_var('vas', 'm s-1',                                                    &
                   &          'meridional wind in 10m (time mean)',                              &
                   &          datatype_flt),                                                     &
                   & grib2_var(0,2,3, ibits, GRID_REFERENCE, GRID_CELL),                         &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),            &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_tas_m = is_variable_in_output(first_output_name_list, var_name="tas_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tas_m
    IF (p_acc%l_tas_m) THEN
       CALL add_var( list, 'tas_m', p_acc%tas,                                                   &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                         &
                   & t_cf_var('tas', 'K',                                                        &
                   &          'temperature in 2m (time mean)',                                   &
                   &          datatype_flt),                                                     &
                   & grib2_var(0,0,0, ibits, GRID_REFERENCE, GRID_CELL),                         &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),            &
                   & isteptype=TSTEP_INSTANT )
    END IF

    p_acc%l_dew2_m = is_variable_in_output(first_output_name_list, var_name="dew2_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_dew2_m
    IF (p_acc%l_dew2_m) THEN
       CALL add_var( list, 'dew2_m', p_acc%dew2,                                                 &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                         &
                   & t_cf_var('dew2', 'K',                                                       &
                   &          'dew point temperature in 2m (time mean)',                         &
                   &          datatype_flt),                                                     &
                   & grib2_var(0,0,6, ibits, GRID_REFERENCE, GRID_CELL),                         &
                   & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),            &
                   & isteptype=TSTEP_INSTANT )
    END IF

    !------------------------------
    ! Temperature tendencies
    !------------------------------
    p_acc%l_tend_ta_m     = is_variable_in_output(first_output_name_list, var_name="tend_ta_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ta_m
    IF (p_acc%l_tend_ta_m) THEN
       cf_desc = t_cf_var('temperature_tendency', 'K s-1',                                    &
            &             'temperature tendency (time mean)',                                 &
            &             datatype_flt)
       grib2_desc = grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ta_m', p_acc%tend_ta,                                        &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ta_dyn_m = is_variable_in_output(first_output_name_list, var_name="tend_ta_dyn_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ta_dyn_m
    IF (p_acc%l_tend_ta_dyn_m) THEN
       cf_desc = t_cf_var('temperature_tendency_dyn', 'K s-1',                                &
            &             'temperature tendency due to resolved dynamics (time mean)',        &
            &             datatype_flt)
       grib2_desc = grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ta_dyn_m', p_acc%tend_ta_dyn,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ta_phy_m = is_variable_in_output(first_output_name_list, var_name="tend_ta_phy_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ta_phy_m
    IF (p_acc%l_tend_ta_phy_m) THEN
       cf_desc = t_cf_var('temperature_tendency_phy', 'K s-1',                                &
            &             'temperature tendency due to param. processes (time mean)',         &
            &             datatype_flt)
       grib2_desc = grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ta_phy_m', p_acc%tend_ta_phy,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ta_rsw_m = is_variable_in_output(first_output_name_list, var_name="tend_ta_rsw_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ta_rsw_m
    IF (p_acc%l_tend_ta_rsw_m) THEN
       cf_desc = t_cf_var('temperature_tendency_rsw', 'K s-1',                                &
            &             'temperature tendency due to shortwave radiation (time mean)',      &
            &             datatype_flt)
       grib2_desc = grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ta_rsw_m', p_acc%tend_ta_rsw,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ta_rlw_m = is_variable_in_output(first_output_name_list, var_name="tend_ta_rlw_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ta_rlw_m
    IF (p_acc%l_tend_ta_rlw_m) THEN
       cf_desc = t_cf_var('temperature_tendency_rlw', 'K s-1',                                &
            &             'temperature tendency due to longwave radiation (time mean)',       &
            &             datatype_flt)
       grib2_desc = grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ta_rlw_m', p_acc%tend_ta_rlw,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ta_rlw_impl_m = is_variable_in_output(first_output_name_list, var_name="tend_ta_rlw_impl_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ta_rlw_impl_m
    IF (p_acc%l_tend_ta_rlw_impl_m) THEN
       cf_desc = t_cf_var('temperature_tendency_rlw_impl', 'K s-1',                           &
            &             'temperature tendency due to LW rad. due to implicit land surface temperature change (time mean)', &
            &             datatype_flt)
       grib2_desc = grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ta_rlw_impl_m', p_acc%tend_ta_rlw_impl,                      &
            &        GRID_UNSTRUCTURED_CELL, ZA_surface, cf_desc, grib2_desc, ldims=shape2d )
    END IF

    p_acc%l_tend_ta_cld_m = is_variable_in_output(first_output_name_list, var_name="tend_ta_cld_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ta_cld_m
    IF (p_acc%l_tend_ta_cld_m) THEN
       cf_desc = t_cf_var('temperature_tendency_cld', 'K s-1',                                &
            &             'temperature tendency due large scale cloud processes (time mean)', &
            &             datatype_flt)
       grib2_desc = grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ta_cld_m', p_acc%tend_ta_cld,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ta_cnv_m = is_variable_in_output(first_output_name_list, var_name="tend_ta_cnv_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ta_cnv_m
    IF (p_acc%l_tend_ta_cnv_m) THEN
       cf_desc = t_cf_var('temperature_tendency_cnv', 'K s-1',                                &
            &             'temperature tendency due convective cloud processes (time mean)',  &
            &             datatype_flt)
       grib2_desc = grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ta_cnv_m', p_acc%tend_ta_cnv,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ta_vdf_m = is_variable_in_output(first_output_name_list, var_name="tend_ta_vdf_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ta_vdf_m
    IF (p_acc%l_tend_ta_vdf_m) THEN
       cf_desc = t_cf_var('temperature_tendency_vdf', 'K s-1',                                &
            &             'temperature tendency due vertical diffusion (time mean)',          &
            &             datatype_flt)
       grib2_desc = grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ta_vdf_m', p_acc%tend_ta_vdf,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ta_gwh_m = is_variable_in_output(first_output_name_list, var_name="tend_ta_gwh_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ta_gwh_m
    IF (p_acc%l_tend_ta_gwh_m) THEN
       cf_desc = t_cf_var('temperature_tendency_gwh', 'K s-1',                                &
            &             'temperature tendency due non-orographic gravity waves (time mean)',&
            &             datatype_flt)
       grib2_desc = grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ta_gwh_m', p_acc%tend_ta_gwh,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ta_sso_m = is_variable_in_output(first_output_name_list, var_name="tend_ta_sso_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ta_sso_m
    IF (p_acc%l_tend_ta_sso_m) THEN
       cf_desc = t_cf_var('temperature_tendency_sso', 'K s-1',                                &
            &             'temperature tendency due sub grid scale orography (time mean)',    &
            &             datatype_flt)
       grib2_desc = grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ta_sso_m', p_acc%tend_ta_sso,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    !------------------------------
    ! U-wind tendencies
    !------------------------------
    p_acc%l_tend_ua_m     = is_variable_in_output(first_output_name_list, var_name="tend_ua_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ua_m
    IF (p_acc%l_tend_ua_m) THEN
       cf_desc = t_cf_var('u_wind_tendency', 'm s-2',                                         &
            &             'u-wind tendency (time mean)',                                      &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ua_m', p_acc%tend_ua,                                        &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ua_dyn_m = is_variable_in_output(first_output_name_list, var_name="tend_ua_dyn_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ua_dyn_m
    IF (p_acc%l_tend_ua_dyn_m) THEN
       cf_desc = t_cf_var('u_wind_tendency_dyn', 'm s-2',                                     &
            &             'u-wind tendency due to resolved dynamics (time mean)',             &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ua_dyn_m', p_acc%tend_ua_dyn,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ua_phy_m = is_variable_in_output(first_output_name_list, var_name="tend_ua_phy_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ua_phy_m
    IF (p_acc%l_tend_ua_phy_m) THEN
       cf_desc = t_cf_var('u_wind_tendency_phy', 'm s-2',                                     &
            &             'u-wind tendency due to param. processes (time mean)',              &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ua_phy_m', p_acc%tend_ua_phy,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ua_cnv_m = is_variable_in_output(first_output_name_list, var_name="tend_ua_cnv_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ua_cnv_m
    IF (p_acc%l_tend_ua_cnv_m) THEN
       cf_desc = t_cf_var('u_wind_tendency_cnv', 'm s-2',                                     &
            &             'u-wind tendency due to convective cloud precesses (time mean)',    &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ua_cnv_m', p_acc%tend_ua_cnv,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ua_vdf_m = is_variable_in_output(first_output_name_list, var_name="tend_ua_vdf_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ua_vdf_m
    IF (p_acc%l_tend_ua_vdf_m) THEN
       cf_desc = t_cf_var('u_wind_tendency_vdf', 'm s-2',                                     &
            &             'u-wind tendency due to vertical diffusion (time mean)',            &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ua_vdf_m', p_acc%tend_ua_vdf,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ua_gwh_m = is_variable_in_output(first_output_name_list, var_name="tend_ua_gwh_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ua_gwh_m
    IF (p_acc%l_tend_ua_gwh_m) THEN
       cf_desc = t_cf_var('u_wind_tendency_gwh', 'm s-2',                                     &
            &             'u-wind tendency due to non-orographic gravity waves (time mean)',  &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ua_gwh_m', p_acc%tend_ua_gwh,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_ua_sso_m = is_variable_in_output(first_output_name_list, var_name="tend_ua_sso_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_ua_sso_m
    IF (p_acc%l_tend_ua_sso_m) THEN
       cf_desc = t_cf_var('u_wind_tendency_sso', 'm s-2',                                     &
            &             'u-wind tendency due to sub grid scale orography (time mean)',      &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_ua_sso_m', p_acc%tend_ua_sso,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    !------------------------------
    ! V-wind tendencies
    !------------------------------
    p_acc%l_tend_va_m     = is_variable_in_output(first_output_name_list, var_name="tend_va_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_va_m
    IF (p_acc%l_tend_va_m) THEN
       cf_desc = t_cf_var('v_wind_tendency', 'm s-2',                                         &
            &             'v-wind tendency (time mean)',                                      &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_va_m', p_acc%tend_va,                                        &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_va_dyn_m = is_variable_in_output(first_output_name_list, var_name="tend_va_dyn_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_va_dyn_m
    IF (p_acc%l_tend_va_dyn_m) THEN
       cf_desc = t_cf_var('v_wind_tendency_dyn', 'm s-2',                                     &
            &             'v-wind tendency due to resolved dynamics (time mean)',             &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_va_dyn_m', p_acc%tend_va_dyn,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_va_phy_m = is_variable_in_output(first_output_name_list, var_name="tend_va_phy_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_va_phy_m
    IF (p_acc%l_tend_va_phy_m) THEN
       cf_desc = t_cf_var('v_wind_tendency_phy', 'm s-2',                                     &
            &             'v-wind tendency due to param. processes (time mean)',              &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_va_phy_m', p_acc%tend_va_phy,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_va_cnv_m = is_variable_in_output(first_output_name_list, var_name="tend_va_cnv_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_va_cnv_m
    IF (p_acc%l_tend_va_cnv_m) THEN
       cf_desc = t_cf_var('v_wind_tendency_cnv', 'm s-2',                                     &
            &             'v-wind tendency due to convective cloud precesses (time mean)',    &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_va_cnv_m', p_acc%tend_va_cnv,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_va_vdf_m = is_variable_in_output(first_output_name_list, var_name="tend_va_vdf_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_va_vdf_m
    IF (p_acc%l_tend_va_vdf_m) THEN
       cf_desc = t_cf_var('v_wind_tendency_vdf', 'm s-2',                                     &
            &             'v-wind tendency due to vertical diffusion (time mean)',            &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_va_vdf_m', p_acc%tend_va_vdf,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_va_gwh_m = is_variable_in_output(first_output_name_list, var_name="tend_va_gwh_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_va_gwh_m
    IF (p_acc%l_tend_va_gwh_m) THEN
       cf_desc = t_cf_var('v_wind_tendency_gwh', 'm s-2',                                     &
            &             'v-wind tendency due to non-orographic gravity waves (time mean)',  &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_va_gwh_m', p_acc%tend_va_gwh,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF

    p_acc%l_tend_va_sso_m = is_variable_in_output(first_output_name_list, var_name="tend_va_sso_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_va_sso_m
    IF (p_acc%l_tend_va_sso_m) THEN
       cf_desc = t_cf_var('v_wind_tendency_sso', 'm s-2',                                     &
            &             'v-wind tendency due to sub grid scale orography (time mean)',      &
            &             datatype_flt)
       grib2_desc = grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( list, 'tend_va_sso_m', p_acc%tend_va_sso,                                &
            &        GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d_c, &
            &        vert_interp=create_vert_interp_metadata(                                 &
            &        vert_intp_type=vintp_types("P","Z","I"),                                 &
            &        vert_intp_method=VINTP_METHOD_LIN,                                       &
            &        l_extrapol=.FALSE. ) )
    END IF


!!$    p_acc%l_tend_hus_m     = is_variable_in_output(first_output_name_list, var_name="tend_hus_m")
!!$    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_hus_m
!!$    IF (p_acc%l_tend_hus_m) THEN
!!$    END IF
!!$
!!$    p_acc%l_tend_hus_dyn_m = is_variable_in_output(first_output_name_list, var_name="tend_hus_dyn_m")
!!$    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_hus_dyn_m
!!$    IF (p_acc%l_tend_hus_dyn_m) THEN
!!$    END IF
!!$
!!$    p_acc%l_tend_hus_phy_m = is_variable_in_output(first_output_name_list, var_name="tend_hus_phy_m")
!!$    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_hus_phy_m
!!$    IF (p_acc%l_tend_hus_phy_m) THEN
!!$    END IF
!!$
!!$    p_acc%l_tend_hus_cld_m = is_variable_in_output(first_output_name_list, var_name="tend_hus_cld_m")
!!$    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_hus_cld_m
!!$    IF (p_acc%l_tend_hus_cld_m) THEN
!!$    END IF
!!$
!!$    p_acc%l_tend_hus_cnv_m = is_variable_in_output(first_output_name_list, var_name="tend_hus_cnv_m")
!!$    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_hus_cnv_m
!!$    IF (p_acc%l_tend_hus_cnv_m) THEN
!!$    END IF
!!$
!!$    p_acc%l_tend_hus_vdf_m = is_variable_in_output(first_output_name_list, var_name="tend_hus_vdf_m")
!!$    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_hus_vdf_m
!!$    IF (p_acc%l_tend_hus_vdf_m) THEN
!!$    END IF
!!$
!!$
!!$    p_acc%l_tend_clw_dtr_m = is_variable_in_output(first_output_name_list, var_name="tend_clw_dtr_m")
!!$    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_clw_dtr_m
!!$    IF (p_acc%l_tend_clw_dtr_m) THEN
!!$    END IF
!!$
!!$    p_acc%l_tend_cli_dtr_m = is_variable_in_output(first_output_name_list, var_name="tend_cli_dtr_m")
!!$    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tend_cli_dtr_m
!!$    IF (p_acc%l_tend_cli_dtr_m) THEN
!!$    END IF
    ! }}}

    p_acc%numberOfAccumulations = 0
  END SUBROUTINE construct_opt_acc


  SUBROUTINE update_opt_acc(acc, nh_prog, rho, nh_diag, subset, levels)
    TYPE(t_nh_acc),  INTENT(INOUT)   :: acc
    TYPE(t_nh_prog), INTENT(IN)      :: nh_prog     ! for jg=1
    REAL(wp), INTENT(IN)             :: rho(:,:,:)
    TYPE(t_nh_diag), INTENT(IN)      :: nh_diag     ! for jg=1
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER , INTENT(IN)             :: levels

    INTEGER :: jt

    INTEGER, PARAMETER :: jg = 1 ! used for prm_field(jg)%...

    !WRITE(message_text,'(a,i2)') '(pre ): numberOfAccumulations:',acc%numberOfAccumulations
    !CALL message('update_opt_nh_acc', TRIM(message_text))
    IF (acc%l_ua_m)    CALL add_fields(acc%u       , nh_diag%u       , subset, levels=levels)
    IF (acc%l_va_m)    CALL add_fields(acc%v       , nh_diag%v       , subset, levels=levels)
    IF (acc%l_wa_m)    CALL add_fields(acc%w       , nh_prog%w       , subset, levels=levels+1)
    IF (acc%l_rho_m) THEN
      CALL add_fields(acc%rho     , rho             , subset, levels=levels)
      CALL dbg_print('RHO Update FROM',nh_prog%rho  ,'opt_diag',5, in_subset=subset)
      CALL dbg_print('RHO Update TO  ',acc%rho      ,'opt_diag',5, in_subset=subset)
    ENDIF
    IF (acc%l_ta_m)    CALL add_fields(acc%temp    , nh_diag%temp    , subset, levels=levels)
    IF (acc%l_ps_m)    CALL add_fields(acc%pres_sfc, nh_diag%pres_sfc, subset)
    IF (acc%l_psl_m)   CALL add_fields(acc%pres_msl, nh_diag%pres_msl, subset)
    IF (acc%l_pfull_m) CALL add_fields(acc%pres    , nh_diag%pres    , subset, levels=levels)
    IF (acc%l_phalf_m) CALL add_fields(acc%pres_ifc, nh_diag%pres_ifc, subset, levels=levels+1)
    IF (acc%l_wap_m)   CALL add_fields(acc%omega   , nh_diag%omega   , subset, levels=levels)

    IF (ntracer > 0) THEN
       IF (acc%l_tracer_m) THEN
          DO jt=1,ntracer
             CALL add_fields(acc%tracer(:,:,:,jt) ,nh_prog%tracer(:,:,:,jt),subset,levels=levels)
             CALL dbg_print('Tracer Update FROM',nh_prog%tracer(:,:,:,jt),'opt_diag',5, in_subset=subset)
             CALL dbg_print('Tracer Update TO  ',acc%tracer(:,:,:,jt),    'opt_diag',5, in_subset=subset)
          END DO
       END IF
    END IF

    IF (acc%l_cosmu0_m)   CALL add_fields(acc%cosmu0         , prm_field(jg)%cosmu0         , subset)
    IF (acc%l_rsdt_m)     CALL add_fields(acc%rsdt           , prm_field(jg)%rsdt           , subset)
    IF (acc%l_hur_m)      CALL add_fields(acc%relhum         , prm_field(jg)%relhum         , subset, levels=levels)
    IF (acc%l_cl_m)       CALL add_fields(acc%aclc           , prm_field(jg)%aclc           , subset, levels=levels)
    IF (acc%l_clt_m)      CALL add_fields(acc%aclcov         , prm_field(jg)%aclcov         , subset)
    IF (acc%l_prlr_m)     CALL add_fields(acc%rsfl           , prm_field(jg)%rsfl           , subset)
    IF (acc%l_prcr_m)     CALL add_fields(acc%rsfc           , prm_field(jg)%rsfc           , subset)
    IF (acc%l_prls_m)     CALL add_fields(acc%ssfl           , prm_field(jg)%ssfl           , subset)
    IF (acc%l_prcs_m)     CALL add_fields(acc%ssfc           , prm_field(jg)%ssfc           , subset)
    IF (acc%l_pr_m)       CALL add_fields(acc%totprec        , prm_field(jg)%totprec        , subset)
    IF (acc%l_prw_m)      CALL add_fields(acc%qvi            , prm_field(jg)%qvi            , subset)
    IF (acc%l_cllvi_m)    CALL add_fields(acc%xlvi           , prm_field(jg)%xlvi           , subset)
    IF (acc%l_clivi_m)    CALL add_fields(acc%xivi           , prm_field(jg)%xivi           , subset)
    IF (acc%l_rsns_m)     CALL add_fields(acc%swflxsfc       , prm_field(jg)%swflxsfc       , subset)
    IF (acc%l_rsnt_m)     CALL add_fields(acc%swflxtoa       , prm_field(jg)%swflxtoa       , subset)
    IF (acc%l_rlns_m)     CALL add_fields(acc%lwflxsfc       , prm_field(jg)%lwflxsfc       , subset)
    IF (acc%l_rlnt_m)     CALL add_fields(acc%lwflxtoa       , prm_field(jg)%lwflxtoa       , subset)
    IF (acc%l_ts_m)       CALL add_fields(acc%tsfc           , prm_field(jg)%tsfc           , subset)
    IF (acc%l_evspsbl_m)  CALL add_fields(acc%evap           , prm_field(jg)%evap           , subset)
    IF (acc%l_hfls_m)     CALL add_fields(acc%lhflx          , prm_field(jg)%lhflx          , subset)
    IF (acc%l_hfss_m)     CALL add_fields(acc%shflx          , prm_field(jg)%shflx          , subset)
    IF (acc%l_tauu_m)     CALL add_fields(acc%u_stress       , prm_field(jg)%u_stress       , subset)
    IF (acc%l_tauv_m)     CALL add_fields(acc%v_stress       , prm_field(jg)%v_stress       , subset)
    IF (acc%l_tauu_sso_m) CALL add_fields(acc%u_stress_sso   , prm_field(jg)%u_stress_sso   , subset)
    IF (acc%l_tauv_sso_m) CALL add_fields(acc%v_stress_sso   , prm_field(jg)%v_stress_sso   , subset)
    IF (acc%l_diss_sso_m) CALL add_fields(acc%dissipation_sso, prm_field(jg)%dissipation_sso, subset)
    IF (acc%l_sic_m)      CALL add_fields(acc%seaice         , prm_field(jg)%seaice         , subset)
    IF (acc%l_sit_m)      CALL add_fields(acc%siced          , prm_field(jg)%siced          , subset)
    IF (acc%l_albedo_m)   CALL add_fields(acc%albedo         , prm_field(jg)%albedo         , subset)
    IF (acc%l_sfcWind_m)  CALL add_fields(acc%sfcWind        , prm_field(jg)%sfcWind        , subset)
    IF (acc%l_uas_m)      CALL add_fields(acc%uas            , prm_field(jg)%uas            , subset)
    IF (acc%l_vas_m)      CALL add_fields(acc%vas            , prm_field(jg)%vas            , subset)
    IF (acc%l_tas_m)      CALL add_fields(acc%tas            , prm_field(jg)%tas            , subset)
    IF (acc%l_dew2_m)     CALL add_fields(acc%dew2           , prm_field(jg)%dew2           , subset)

    IF (acc%l_tend_ta_m    )      CALL add_fields(acc%tend_ta         , prm_tend(jg)%temp         , subset, levels=levels)
    IF (acc%l_tend_ta_dyn_m)      CALL add_fields(acc%tend_ta_dyn     , prm_tend(jg)%temp_dyn     , subset, levels=levels)
    IF (acc%l_tend_ta_phy_m)      CALL add_fields(acc%tend_ta_phy     , prm_tend(jg)%temp_phy     , subset, levels=levels)
    IF (acc%l_tend_ta_rsw_m)      CALL add_fields(acc%tend_ta_rsw     , prm_tend(jg)%temp_rsw     , subset, levels=levels)
    IF (acc%l_tend_ta_rlw_m)      CALL add_fields(acc%tend_ta_rlw     , prm_tend(jg)%temp_rlw     , subset, levels=levels)
    IF (acc%l_tend_ta_rlw_impl_m) CALL add_fields(acc%tend_ta_rlw_impl, prm_tend(jg)%temp_rlw_impl, subset)
    IF (acc%l_tend_ta_cld_m)      CALL add_fields(acc%tend_ta_cld     , prm_tend(jg)%temp_cld     , subset, levels=levels)
    IF (acc%l_tend_ta_cnv_m)      CALL add_fields(acc%tend_ta_cnv     , prm_tend(jg)%temp_cnv     , subset, levels=levels)
    IF (acc%l_tend_ta_vdf_m)      CALL add_fields(acc%tend_ta_vdf     , prm_tend(jg)%temp_vdf     , subset, levels=levels)
    IF (acc%l_tend_ta_gwh_m)      CALL add_fields(acc%tend_ta_gwh     , prm_tend(jg)%temp_gwh     , subset, levels=levels)
    IF (acc%l_tend_ta_sso_m)      CALL add_fields(acc%tend_ta_sso     , prm_tend(jg)%temp_sso     , subset, levels=levels)

    IF (acc%l_tend_ua_m    ) CALL add_fields(acc%tend_ua    , prm_tend(jg)%u    , subset, levels=levels)
    IF (acc%l_tend_ua_dyn_m) CALL add_fields(acc%tend_ua_dyn, prm_tend(jg)%u_dyn, subset, levels=levels)
    IF (acc%l_tend_ua_phy_m) CALL add_fields(acc%tend_ua_phy, prm_tend(jg)%u_phy, subset, levels=levels)
    IF (acc%l_tend_ua_cnv_m) CALL add_fields(acc%tend_ua_cnv, prm_tend(jg)%u_cnv, subset, levels=levels)
    IF (acc%l_tend_ua_vdf_m) CALL add_fields(acc%tend_ua_vdf, prm_tend(jg)%u_vdf, subset, levels=levels)
    IF (acc%l_tend_ua_gwh_m) CALL add_fields(acc%tend_ua_gwh, prm_tend(jg)%u_gwh, subset, levels=levels)
    IF (acc%l_tend_ua_sso_m) CALL add_fields(acc%tend_ua_sso, prm_tend(jg)%u_sso, subset, levels=levels)

    IF (acc%l_tend_va_m    ) CALL add_fields(acc%tend_va    , prm_tend(jg)%v    , subset, levels=levels)
    IF (acc%l_tend_va_dyn_m) CALL add_fields(acc%tend_va_dyn, prm_tend(jg)%v_dyn, subset, levels=levels)
    IF (acc%l_tend_va_phy_m) CALL add_fields(acc%tend_va_phy, prm_tend(jg)%v_phy, subset, levels=levels)
    IF (acc%l_tend_va_cnv_m) CALL add_fields(acc%tend_va_cnv, prm_tend(jg)%v_cnv, subset, levels=levels)
    IF (acc%l_tend_va_vdf_m) CALL add_fields(acc%tend_va_vdf, prm_tend(jg)%v_vdf, subset, levels=levels)
    IF (acc%l_tend_va_gwh_m) CALL add_fields(acc%tend_va_gwh, prm_tend(jg)%v_gwh, subset, levels=levels)
    IF (acc%l_tend_va_sso_m) CALL add_fields(acc%tend_va_sso, prm_tend(jg)%v_sso, subset, levels=levels)

    acc%numberOfAccumulations = acc%numberOfAccumulations + 1
    !WRITE(message_text,'(a,i2)') '(post): numberOfAccumulations:',acc%numberOfAccumulations
    !CALL message('update_opt_nh_acc', TRIM(message_text))

  END SUBROUTINE update_opt_acc


  SUBROUTINE reset_opt_acc(acc)
    TYPE(t_nh_acc) :: acc!(n_dom)

    INTEGER :: jt

    IF (acc%l_ua_m)    acc%u        = 0.0_wp
    IF (acc%l_va_m)    acc%v        = 0.0_wp
    IF (acc%l_wa_m)    acc%w        = 0.0_wp
    IF (acc%l_rho_m)   acc%rho      = 0.0_wp
    IF (acc%l_ta_m)    acc%temp     = 0.0_wp
    IF (acc%l_ps_m)    acc%pres_sfc = 0.0_wp
    IF (acc%l_psl_m)   acc%pres_msl = 0.0_wp
    IF (acc%l_pfull_m) acc%pres     = 0.0_wp
    IF (acc%l_phalf_m) acc%pres_ifc = 0.0_wp
    IF (acc%l_wap_m)   acc%omega    = 0.0_wp

    IF (ntracer > 0) THEN
       IF (acc%l_tracer_m) THEN
          DO jt=1,ntracer
             acc%tracer(:,:,:,jt) = 0.0_wp
          END DO
       END IF
    END IF

    IF (acc%l_cosmu0_m)   acc%cosmu0          = 0.0_wp
    IF (acc%l_rsdt_m)     acc%rsdt            = 0.0_wp
    IF (acc%l_hur_m)      acc%relhum          = 0.0_wp
    IF (acc%l_cl_m)       acc%aclc            = 0.0_wp
    IF (acc%l_clt_m)      acc%aclcov          = 0.0_wp
    IF (acc%l_prlr_m)     acc%rsfl            = 0.0_wp
    IF (acc%l_prls_m)     acc%ssfl            = 0.0_wp
    IF (acc%l_prcr_m)     acc%rsfc            = 0.0_wp
    IF (acc%l_prcs_m)     acc%ssfc            = 0.0_wp
    IF (acc%l_pr_m)       acc%totprec         = 0.0_wp
    IF (acc%l_prw_m)      acc%qvi             = 0.0_wp
    IF (acc%l_cllvi_m)    acc%xlvi            = 0.0_wp
    IF (acc%l_clivi_m)    acc%xivi            = 0.0_wp
    IF (acc%l_rsns_m)     acc%swflxsfc        = 0.0_wp
    IF (acc%l_rsnt_m)     acc%swflxtoa        = 0.0_wp
    IF (acc%l_rlns_m)     acc%lwflxsfc        = 0.0_wp
    IF (acc%l_rlnt_m)     acc%lwflxtoa        = 0.0_wp
    IF (acc%l_ts_m)       acc%tsfc            = 0.0_wp
    IF (acc%l_evspsbl_m)  acc%evap            = 0.0_wp
    IF (acc%l_hfls_m)     acc%lhflx           = 0.0_wp
    IF (acc%l_hfss_m)     acc%shflx           = 0.0_wp
    IF (acc%l_tauu_m)     acc%u_stress        = 0.0_wp
    IF (acc%l_tauv_m)     acc%v_stress        = 0.0_wp
    IF (acc%l_tauu_sso_m) acc%u_stress_sso    = 0.0_wp
    IF (acc%l_tauv_sso_m) acc%v_stress_sso    = 0.0_wp
    IF (acc%l_diss_sso_m) acc%dissipation_sso = 0.0_wp
    IF (acc%l_sic_m)      acc%seaice          = 0.0_wp
    IF (acc%l_sit_m)      acc%siced           = 0.0_wp
    IF (acc%l_albedo_m)   acc%albedo          = 0.0_wp
    IF (acc%l_sfcWind_m)  acc%sfcWind         = 0.0_wp
    IF (acc%l_uas_m)      acc%uas             = 0.0_wp
    IF (acc%l_vas_m)      acc%vas             = 0.0_wp
    IF (acc%l_tas_m)      acc%tas             = 0.0_wp
    IF (acc%l_dew2_m)     acc%dew2            = 0.0_wp

    IF (acc%l_tend_ta_m    )      acc%tend_ta          = 0.0_wp
    IF (acc%l_tend_ta_dyn_m)      acc%tend_ta_dyn      = 0.0_wp
    IF (acc%l_tend_ta_phy_m)      acc%tend_ta_phy      = 0.0_wp
    IF (acc%l_tend_ta_rsw_m)      acc%tend_ta_rsw      = 0.0_wp
    IF (acc%l_tend_ta_rlw_m)      acc%tend_ta_rlw      = 0.0_wp
    IF (acc%l_tend_ta_rlw_impl_m) acc%tend_ta_rlw_impl = 0.0_wp
    IF (acc%l_tend_ta_cld_m)      acc%tend_ta_cld      = 0.0_wp
    IF (acc%l_tend_ta_cnv_m)      acc%tend_ta_cnv      = 0.0_wp
    IF (acc%l_tend_ta_vdf_m)      acc%tend_ta_vdf      = 0.0_wp
    IF (acc%l_tend_ta_gwh_m)      acc%tend_ta_gwh      = 0.0_wp
    IF (acc%l_tend_ta_sso_m)      acc%tend_ta_sso      = 0.0_wp

    IF (acc%l_tend_ua_m    ) acc%tend_ua     = 0.0_wp
    IF (acc%l_tend_ua_dyn_m) acc%tend_ua_dyn = 0.0_wp
    IF (acc%l_tend_ua_phy_m) acc%tend_ua_phy = 0.0_wp
    IF (acc%l_tend_ua_cnv_m) acc%tend_ua_cnv = 0.0_wp
    IF (acc%l_tend_ua_vdf_m) acc%tend_ua_vdf = 0.0_wp
    IF (acc%l_tend_ua_gwh_m) acc%tend_ua_gwh = 0.0_wp
    IF (acc%l_tend_ua_sso_m) acc%tend_ua_sso = 0.0_wp

    IF (acc%l_tend_va_m    ) acc%tend_va     = 0.0_wp
    IF (acc%l_tend_va_dyn_m) acc%tend_va_dyn = 0.0_wp
    IF (acc%l_tend_va_phy_m) acc%tend_va_phy = 0.0_wp
    IF (acc%l_tend_va_cnv_m) acc%tend_va_cnv = 0.0_wp
    IF (acc%l_tend_va_vdf_m) acc%tend_va_vdf = 0.0_wp
    IF (acc%l_tend_va_gwh_m) acc%tend_va_gwh = 0.0_wp
    IF (acc%l_tend_va_sso_m) acc%tend_va_sso = 0.0_wp

    acc%numberOfAccumulations = 0

  END SUBROUTINE reset_opt_acc


  SUBROUTINE calc_mean_opt_acc(acc)
    TYPE(t_nh_acc) :: acc!(n_dom)

    INTEGER  :: jt
    REAL(wp) :: xfactor

    xfactor = 1._wp/REAL(acc%numberOfAccumulations,wp)

    IF (acc%l_ua_m)    acc%u        = acc%u        *xfactor
    IF (acc%l_va_m)    acc%v        = acc%v        *xfactor
    IF (acc%l_wa_m)    acc%w        = acc%w        *xfactor
    IF (acc%l_rho_m)   acc%rho      = acc%rho      *xfactor
    IF (acc%l_ta_m)    acc%temp     = acc%temp     *xfactor
    IF (acc%l_ps_m)    acc%pres_sfc = acc%pres_sfc *xfactor
    IF (acc%l_psl_m)   acc%pres_msl = acc%pres_msl *xfactor
    IF (acc%l_pfull_m) acc%pres     = acc%pres     *xfactor
    IF (acc%l_phalf_m) acc%pres_ifc = acc%pres_ifc *xfactor
    IF (acc%l_wap_m)   acc%omega    = acc%omega    *xfactor

    IF (ntracer > 0) THEN
       IF (acc%l_tracer_m) THEN
          DO jt=1,ntracer
             acc%tracer(:,:,:,jt) = acc%tracer(:,:,:,jt) *xfactor
          END DO
       END IF
    END IF

    IF (acc%l_cosmu0_m)   acc%cosmu0          = acc%cosmu0          *xfactor
    IF (acc%l_rsdt_m)     acc%rsdt            = acc%rsdt            *xfactor
    IF (acc%l_hur_m)      acc%relhum          = acc%relhum          *xfactor
    IF (acc%l_cl_m)       acc%aclc            = acc%aclc            *xfactor
    IF (acc%l_clt_m)      acc%aclcov          = acc%aclcov          *xfactor
    IF (acc%l_prlr_m)     acc%rsfl            = acc%rsfl            *xfactor
    IF (acc%l_prcr_m)     acc%rsfc            = acc%rsfc            *xfactor
    IF (acc%l_prls_m)     acc%ssfl            = acc%ssfl            *xfactor
    IF (acc%l_prcs_m)     acc%ssfc            = acc%ssfc            *xfactor
    IF (acc%l_pr_m)       acc%totprec         = acc%totprec         *xfactor
    IF (acc%l_prw_m)      acc%qvi             = acc%qvi             *xfactor
    IF (acc%l_cllvi_m)    acc%xlvi            = acc%xlvi            *xfactor
    IF (acc%l_clivi_m)    acc%xivi            = acc%xivi            *xfactor
    IF (acc%l_rsns_m)     acc%swflxsfc        = acc%swflxsfc        *xfactor
    IF (acc%l_rsnt_m)     acc%swflxtoa        = acc%swflxtoa        *xfactor
    IF (acc%l_rlns_m)     acc%lwflxsfc        = acc%lwflxsfc        *xfactor
    IF (acc%l_rlnt_m)     acc%lwflxtoa        = acc%lwflxtoa        *xfactor
    IF (acc%l_ts_m)       acc%tsfc            = acc%tsfc            *xfactor
    IF (acc%l_evspsbl_m)  acc%evap            = acc%evap            *xfactor
    IF (acc%l_hfls_m)     acc%lhflx           = acc%lhflx           *xfactor
    IF (acc%l_hfss_m)     acc%shflx           = acc%shflx           *xfactor
    IF (acc%l_tauu_m)     acc%u_stress        = acc%u_stress        *xfactor
    IF (acc%l_tauv_m)     acc%v_stress        = acc%v_stress        *xfactor
    IF (acc%l_tauu_sso_m) acc%u_stress_sso    = acc%u_stress_sso    *xfactor
    IF (acc%l_tauv_sso_m) acc%v_stress_sso    = acc%v_stress_sso    *xfactor
    IF (acc%l_diss_sso_m) acc%dissipation_sso = acc%dissipation_sso *xfactor
    IF (acc%l_sic_m)      acc%seaice          = acc%seaice          *xfactor
    IF (acc%l_sit_m)      acc%siced           = acc%siced           *xfactor
    IF (acc%l_albedo_m)   acc%albedo          = acc%albedo          *xfactor
    IF (acc%l_sfcWind_m)  acc%sfcWind         = acc%sfcWind         *xfactor
    IF (acc%l_uas_m)      acc%uas             = acc%uas             *xfactor
    IF (acc%l_vas_m)      acc%vas             = acc%vas             *xfactor
    IF (acc%l_tas_m)      acc%tas             = acc%tas             *xfactor
    IF (acc%l_dew2_m)     acc%dew2            = acc%dew2            *xfactor

    IF (acc%l_tend_ta_m    )      acc%tend_ta          = acc%tend_ta          *xfactor
    IF (acc%l_tend_ta_dyn_m)      acc%tend_ta_dyn      = acc%tend_ta_dyn      *xfactor
    IF (acc%l_tend_ta_phy_m)      acc%tend_ta_phy      = acc%tend_ta_phy      *xfactor
    IF (acc%l_tend_ta_rsw_m)      acc%tend_ta_rsw      = acc%tend_ta_rsw      *xfactor
    IF (acc%l_tend_ta_rlw_m)      acc%tend_ta_rlw      = acc%tend_ta_rlw      *xfactor
    IF (acc%l_tend_ta_rlw_impl_m) acc%tend_ta_rlw_impl = acc%tend_ta_rlw_impl *xfactor
    IF (acc%l_tend_ta_cld_m)      acc%tend_ta_cld      = acc%tend_ta_cld      *xfactor
    IF (acc%l_tend_ta_cnv_m)      acc%tend_ta_cnv      = acc%tend_ta_cnv      *xfactor
    IF (acc%l_tend_ta_vdf_m)      acc%tend_ta_vdf      = acc%tend_ta_vdf      *xfactor
    IF (acc%l_tend_ta_gwh_m)      acc%tend_ta_gwh      = acc%tend_ta_gwh      *xfactor
    IF (acc%l_tend_ta_sso_m)      acc%tend_ta_sso      = acc%tend_ta_sso      *xfactor

    IF (acc%l_tend_ua_m    ) acc%tend_ua     = acc%tend_ua     *xfactor
    IF (acc%l_tend_ua_dyn_m) acc%tend_ua_dyn = acc%tend_ua_dyn *xfactor
    IF (acc%l_tend_ua_phy_m) acc%tend_ua_phy = acc%tend_ua_phy *xfactor
    IF (acc%l_tend_ua_cnv_m) acc%tend_ua_cnv = acc%tend_ua_cnv *xfactor
    IF (acc%l_tend_ua_vdf_m) acc%tend_ua_vdf = acc%tend_ua_vdf *xfactor
    IF (acc%l_tend_ua_gwh_m) acc%tend_ua_gwh = acc%tend_ua_gwh *xfactor
    IF (acc%l_tend_ua_sso_m) acc%tend_ua_sso = acc%tend_ua_sso *xfactor

    IF (acc%l_tend_va_m    ) acc%tend_va     = acc%tend_va     *xfactor
    IF (acc%l_tend_va_dyn_m) acc%tend_va_dyn = acc%tend_va_dyn *xfactor
    IF (acc%l_tend_va_phy_m) acc%tend_va_phy = acc%tend_va_phy *xfactor
    IF (acc%l_tend_va_cnv_m) acc%tend_va_cnv = acc%tend_va_cnv *xfactor
    IF (acc%l_tend_va_vdf_m) acc%tend_va_vdf = acc%tend_va_vdf *xfactor
    IF (acc%l_tend_va_gwh_m) acc%tend_va_gwh = acc%tend_va_gwh *xfactor
    IF (acc%l_tend_va_sso_m) acc%tend_va_sso = acc%tend_va_sso *xfactor

  END SUBROUTINE calc_mean_opt_acc

  !-------------
  !
  !> Add optional diagnostic variable lists (might remain empty)
  !
  SUBROUTINE construct_opt_diag(p_patch, l_init_pz)
    TYPE(t_patch),        INTENT(IN)   :: p_patch(n_dom)
    LOGICAL,              INTENT(IN)   :: l_init_pz

    ! local variables
    CHARACTER(*), PARAMETER :: routine =  &
      &  TRIM("mo_opt_diagnostics:construct_opt_diag")
    INTEGER                            :: jg, ist
    CHARACTER(len=MAX_CHAR_LENGTH)     :: listname

    ! initialize data structure for optional diagnostics
    ALLOCATE(p_nh_opt_diag(n_dom), STAT=ist)
    IF (ist /= SUCCESS) &
      CALL finish (TRIM(routine), 'Allocation of optional diagnostics failed')

    DO jg = 1, n_dom

      WRITE(listname,'(a,i2.2)') 'nh_state_opt_diag_of_domain_',jg
      CALL new_var_list( p_nh_opt_diag(jg)%opt_diag_list, TRIM(listname), &
        & patch_id=p_patch(jg)%id, vlevel_type=level_type_ml )
      CALL default_var_list_settings( p_nh_opt_diag(jg)%opt_diag_list,    &
        & lrestart=.FALSE. )

      IF (.NOT. l_init_pz) CYCLE

      WRITE(listname,'(a,i2.2)') 'nh_state_opt_diag_z_of_domain_',jg
      CALL new_var_list( p_nh_opt_diag(jg)%opt_diag_list_z, TRIM(listname), &
        & patch_id=p_patch(jg)%id, vlevel_type=level_type_hl )
      CALL default_var_list_settings( p_nh_opt_diag(jg)%opt_diag_list_z,    &
        & lrestart=.FALSE. )

      WRITE(listname,'(a,i2.2)') 'nh_state_opt_diag_p_of_domain_',jg
      CALL new_var_list( p_nh_opt_diag(jg)%opt_diag_list_p, TRIM(listname), &
        & patch_id=p_patch(jg)%id, vlevel_type=level_type_pl )
      CALL default_var_list_settings( p_nh_opt_diag(jg)%opt_diag_list_p,    &
        & lrestart=.FALSE. )

      WRITE(listname,'(a,i2.2)') 'nh_state_opt_diag_i_of_domain_',jg
      CALL new_var_list( p_nh_opt_diag(jg)%opt_diag_list_i, TRIM(listname), &
        & patch_id=p_patch(jg)%id, vlevel_type=level_type_il )
      CALL default_var_list_settings( p_nh_opt_diag(jg)%opt_diag_list_i,    &
        & lrestart=.FALSE. )

      WRITE(listname,'(a,i2.2)') 'nh_accumulation_for_ProgAndDiag_of_domain_',jg
      CALL new_var_list( p_nh_opt_diag(jg)%opt_acc_list, TRIM(listname), &
        & patch_id=p_patch(jg)%id, vlevel_type=level_type_ml )
      CALL default_var_list_settings( p_nh_opt_diag(jg)%opt_acc_list,    &
        & lrestart=.FALSE.,loutput=.TRUE. )
    ENDDO ! jg

    ! provisional construction of memory for a hardwired set of variables on domain 1
    CALL construct_opt_acc( p_patch(1),                    &
                          & p_nh_opt_diag(1)%opt_acc_list, &
                          & p_nh_opt_diag(1)%acc           )

  END SUBROUTINE construct_opt_diag


  !-------------
  !
  !> Clear optional diagnostic variable lists
  !
  SUBROUTINE destruct_opt_diag()
    ! local variables
    CHARACTER(*), PARAMETER :: routine =  &
      &  TRIM("mo_opt_diagnostics:destruct_opt_diag")
    INTEGER :: jg, ist

    DO jg = 1, n_dom
      CALL delete_var_list( p_nh_opt_diag(jg)%opt_diag_list_z )
      CALL delete_var_list( p_nh_opt_diag(jg)%opt_diag_list_p )
      CALL delete_var_list( p_nh_opt_diag(jg)%opt_diag_list_i )
      CALL delete_var_list( p_nh_opt_diag(jg)%opt_diag_list   )
      CALL delete_var_list( p_nh_opt_diag(jg)%opt_acc_list    )
    ENDDO ! jg

    ! Delete optional diagnostics
    DEALLOCATE(p_nh_opt_diag, STAT=ist)
    IF (ist /= SUCCESS) &
      CALL finish(TRIM(routine),'Deallocation for optional diagnostics failed.')

  END SUBROUTINE destruct_opt_diag


  !-------------
  !>
  ! Initialize a variable containing coefficient tables for vertical
  ! interpolation. There exist to different kinds of coefficients: For
  ! p- and for z-level-interpolation.
  SUBROUTINE vcoeff_lin_allocate(nblks, nlev, vcoeff_lin)
    INTEGER,                   INTENT(IN)    :: nblks
    INTEGER,                   INTENT(IN)    :: nlev
    TYPE(t_vcoeff_lin),        INTENT(INOUT) :: vcoeff_lin

    CHARACTER(*), PARAMETER :: routine = TRIM("mo_opt_diagnostics:vcoeff_lin_allocate")
    INTEGER :: ierrstat

    ! real(wp)
    ALLOCATE(vcoeff_lin%wfac_lin(nproma,nlev,nblks), STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! integer
    ALLOCATE( vcoeff_lin%idx0_lin(nproma,nlev,nblks), vcoeff_lin%bot_idx_lin(nproma,nblks),   &
      &       STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ALLOCATE( vcoeff_lin%wfacpbl1(nproma,nblks), vcoeff_lin%wfacpbl2(nproma,nblks),           &
      &       STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ALLOCATE( vcoeff_lin%kpbl1(nproma,nblks), vcoeff_lin%kpbl2(nproma,nblks),                 &
      &       STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! Initialization
    vcoeff_lin%wfac_lin    = 0._wp
    vcoeff_lin%idx0_lin    = 0
    vcoeff_lin%bot_idx_lin = 0
    vcoeff_lin%wfacpbl1    = 0._wp
    vcoeff_lin%wfacpbl2    = 0._wp
    vcoeff_lin%kpbl1       = 0
    vcoeff_lin%kpbl2       = 0
  END SUBROUTINE vcoeff_lin_allocate


  !-------------
  !>
  ! Initialize a variable containing coefficient tables for cubic
  ! vertical interpolation.
  !
  SUBROUTINE vcoeff_cub_allocate(nblks, nlev, vcoeff_cub)
    INTEGER,                   INTENT(IN)    :: nblks
    INTEGER,                   INTENT(IN)    :: nlev
    TYPE(t_vcoeff_cub),        INTENT(INOUT) :: vcoeff_cub

    CHARACTER(*), PARAMETER :: routine = TRIM("mo_opt_diagnostics:vcoeff_cub_allocate")
    INTEGER :: ierrstat

    ! real(wp)
    ALLOCATE( vcoeff_cub%coef1(nproma,nlev,nblks), vcoeff_cub%coef2(nproma,nlev,nblks),     &
      &       vcoeff_cub%coef3(nproma,nlev,nblks), STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! integer
    ALLOCATE( vcoeff_cub%idx0_cub(nproma,nlev,nblks), vcoeff_cub%bot_idx_cub(nproma,nblks), &
      &       STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! Initialization
    vcoeff_cub%coef1       = 0._wp
    vcoeff_cub%coef2       = 0._wp
    vcoeff_cub%coef3       = 0._wp
    vcoeff_cub%idx0_cub    = 0
    vcoeff_cub%bot_idx_cub = 0
  END SUBROUTINE vcoeff_cub_allocate


  !-------------
  !>
  ! Initialize a variable containing coefficient tables for vertical
  ! interpolation. There exist to different kinds of coefficients: For
  ! p- and for z-level-interpolation.
  SUBROUTINE vcoeff_allocate(nblks_c, nblks_e, nlev, vcoeff)
    INTEGER,                           INTENT(IN)    :: nblks_c, nblks_e
    INTEGER,                           INTENT(IN)    :: nlev
    TYPE(t_vcoeff),                    INTENT(INOUT) :: vcoeff

!!$    CHARACTER(*), PARAMETER :: routine = TRIM("mo_opt_diagnostics:vcoeff_allocate")

    IF (.NOT. vcoeff%l_allocated) THEN
      CALL vcoeff_lin_allocate(nblks_c, nlev, vcoeff%lin_cell)
      CALL vcoeff_lin_allocate(nblks_c, nlev, vcoeff%lin_cell_nlevp1)
      CALL vcoeff_lin_allocate(nblks_e, nlev, vcoeff%lin_edge)

      ! CUBIC interpolation coefficients:
      CALL vcoeff_cub_allocate(nblks_c, nlev, vcoeff%cub_cell)
      CALL vcoeff_cub_allocate(nblks_e, nlev, vcoeff%cub_edge)

      vcoeff%l_allocated = .TRUE.
    END IF
  END SUBROUTINE vcoeff_allocate


  !-------------
  !>
  ! Clear a coefficient tables for linear vertical interpolation.
  SUBROUTINE vcoeff_lin_deallocate(vcoeff_lin)
    TYPE(t_vcoeff_lin), INTENT(INOUT) :: vcoeff_lin

    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_opt_diagnostics:vcoeff_lin_deallocate")
    INTEGER :: ierrstat

    ! real(wp)
    DEALLOCATE( vcoeff_lin%wfac_lin, STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    ! integer
    DEALLOCATE( vcoeff_lin%idx0_lin, vcoeff_lin%bot_idx_lin, STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    DEALLOCATE( vcoeff_lin%wfacpbl1, vcoeff_lin%wfacpbl2, STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    DEALLOCATE( vcoeff_lin%kpbl1, vcoeff_lin%kpbl2, STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  END SUBROUTINE vcoeff_lin_deallocate


  !-------------
  !>
  ! Clear a coefficient tables for cubic vertical interpolation.
  SUBROUTINE vcoeff_cub_deallocate(vcoeff_cub)
    TYPE(t_vcoeff_cub), INTENT(INOUT) :: vcoeff_cub

    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_opt_diagnostics:vcoeff_cub_deallocate")
    INTEGER :: ierrstat

    ! CUBIC interpolation coefficients:
    ! real(wp)
    DEALLOCATE( vcoeff_cub%coef1, vcoeff_cub%coef2, vcoeff_cub%coef3, STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    ! integer
    DEALLOCATE( vcoeff_cub%idx0_cub, vcoeff_cub%bot_idx_cub, STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  END SUBROUTINE vcoeff_cub_deallocate


  !-------------
  !>
  ! Clear a variable containing coefficient tables for vertical
  ! interpolation. There exist to different kinds of coefficients: For
  ! p-, z- and for i-level-interpolation.
  SUBROUTINE vcoeff_deallocate(vcoeff)
    TYPE(t_vcoeff), INTENT(INOUT) :: vcoeff

!!$    CHARACTER(*), PARAMETER :: routine = &
!!$      &  TRIM("mo_opt_diagnostics:vcoeff_deallocate")

    ! deallocate coefficient tables:
    IF (vcoeff%l_allocated) THEN
      CALL vcoeff_lin_deallocate(vcoeff%lin_cell)
      CALL vcoeff_lin_deallocate(vcoeff%lin_cell_nlevp1)
      CALL vcoeff_lin_deallocate(vcoeff%lin_edge)

      call vcoeff_cub_deallocate(vcoeff%cub_cell)
      call vcoeff_cub_deallocate(vcoeff%cub_edge)

      vcoeff%l_allocated = .FALSE.
    END IF

    vcoeff%l_initialized = .FALSE.
  END SUBROUTINE vcoeff_deallocate

END MODULE mo_opt_diagnostics



