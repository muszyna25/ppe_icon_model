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
  USE mo_parallel_config,      ONLY: nproma, use_dp_mpi2io
  USE mo_linked_list,          ONLY: t_var_list
  USE mo_model_domain,         ONLY: t_patch, t_subset_range
  USE mo_nonhydro_types,       ONLY: t_nh_diag,t_nh_prog
  USE mo_echam_phy_memory,     ONLY: prm_field
  USE mo_impl_constants,       ONLY: SUCCESS, MAX_CHAR_LENGTH,           &
    &                                VINTP_METHOD_QV,                    &
    &                                VINTP_METHOD_LIN
  USE mo_exception,            ONLY: finish
  USE mo_grid_config,          ONLY: n_dom
  USE mo_run_config,           ONLY: ntracer,iqv,iqc,iqi
  USE mo_advection_config,     ONLY: t_advection_config, advection_config
  USE mo_cdi_constants,        ONLY: GRID_UNSTRUCTURED_CELL, GRID_REFERENCE,         &
    &                                GRID_CELL, ZA_HYBRID, ZA_SURFACE,               &
    &                                ZA_MEANSEA, DATATYPE_FLT32, DATATYPE_PACK16,    &
    &                                DATATYPE_PACK24, TSTEP_INSTANT,                 &
    &                                DATATYPE_FLT64
  USE mo_var_list,             ONLY: default_var_list_settings,     &
    &                                new_var_list, delete_var_list, add_var, add_ref
  USE mo_var_list_element,     ONLY: level_type_ml, level_type_pl,  &
    &                                level_type_hl, level_type_il
  USE mo_gribout_config,       ONLY: gribout_config
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_grib2,                ONLY: t_grib2_var
  USE mo_var_metadata,         ONLY: create_tracer_metadata,                 &
    &                                create_vert_interp_metadata,            &
    &                                groups, vintp_types
  USE mo_statistics,           ONLY: add_fields
  USE mo_util_dbg_prnt,        ONLY: dbg_print
  USE mo_exception,            ONLY: message, message_text

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
#ifdef _CRAYFTN
    , CONTIGUOUS        &
#endif
    &  ::               &
    &  rho(:,:,:),      &
    &  qv(:,:,:),       &
    &  qc(:,:,:),       &
    &  qi(:,:,:),       &
    &  temp(:,:,:),     &
    &  pres_sfc(:,:),   &
    &  pres_msl(:,:),   &
    &  u(:,:,:),        &
    &  v(:,:,:),        &
    &  cosmu0(:,:),     &
    &  flxdwswtoa(:,:), &
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
    &  tracer(:,:,:,:)

    TYPE(t_pointer_3d_wp),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer

    ! Internal counter for accumulation operations
    INTEGER :: numberOfAccumulations

    LOGICAL :: l_pres_msl
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
  SUBROUTINE construct_opt_acc(p_patch,list,p_acc,echam_forcing_active, l_pres_msl)
    TYPE(t_patch),        INTENT(IN) :: p_patch
    TYPE(t_var_list)                 :: list
    TYPE(t_nh_acc)                   :: p_acc
    LOGICAL , INTENT(IN)             :: echam_forcing_active
    LOGICAL , OPTIONAL, INTENT(IN)   :: l_pres_msl

    ! LOCAL ===================================================================
    INTEGER :: nblks_c       !< number of cell blocks to allocate
!!$    INTEGER :: nblks_e       !< number of edge blocks to allocate
!!$    INTEGER :: nblks_v       !< number of vertex blocks to allocate

    INTEGER :: nlev
!!$    INTEGER :: nlevp1

    INTEGER :: jt

    INTEGER :: shape2d  (2)
    INTEGER :: shape2d_c(2), shape3d_c(3), shape4d_c(4)
!!$    INTEGER :: shape2d_e(2), shape3d_e(3)
!!$    INTEGER ::               shape3d_v(3)
 
    INTEGER :: ibits,iextbits     !< "entropy" of horizontal slice
    INTEGER :: DATATYPE_PACK_VAR, dataType  !< variable "entropy" for some thermodynamic fields

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
!!$    nlevp1 = p_patch%nlevp1

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice
    iextbits = DATATYPE_PACK24

    IF (gribout_config(p_patch%id)%lgribout_24bit) THEN  ! analysis
      ! higher accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK24
    ELSE
      ! standard accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK16
    ENDIF

    ! pointer to advection_config(jg) to save some paperwork
    advconf => advection_config(p_patch%id)

    ! predefined array shapes
    shape2d_c     = (/nproma,          nblks_c    /)
    shape2d       = shape2d_c
    shape3d_c     = (/nproma, nlev   , nblks_c    /)
    shape4d_c     = (/nproma, nlev   , nblks_c, ntracer     /)
!!$    shape2d_e     = (/nproma,          nblks_e    /)
!!$    shape3d_e     = (/nproma, nlev   , nblks_e    /)
!!$    shape3d_v     = (/nproma, nlev   , nblks_v    /)

    IF ( use_dp_mpi2io ) THEN
      dataType = DATATYPE_FLT64
    ELSE
      dataType = DATATYPE_FLT32
    ENDIF

    ! PROGS {{{
    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'Zonal wind (time mean)', dataType)
    grib2_desc = t_grib2_var(0, 2, 2, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'ua_m', p_acc%u,                                        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
                & ldims=shape3d_c,                                              &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I") ),                  &
                & in_group=groups("prog_timemean","atmo_timemean") )

    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'Meridional wind (time mean)', dataType)
    grib2_desc = t_grib2_var(0, 2, 3, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'va_m', p_acc%v,                                        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
                & ldims=shape3d_c,                                              &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I") ),                  &
                & in_group=groups("prog_timemean","atmo_timemean") )

    cf_desc    = t_cf_var('air_temperature', 'K', 'Temperature', dataType)
    grib2_desc = t_grib2_var(0, 0, 0, DATATYPE_PACK_VAR, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'ta_m', p_acc%temp,                                     &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
                & ldims=shape3d_c,                                              &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("P","Z","I"),          &
                &             vert_intp_method=VINTP_METHOD_LIN ),              &
                & in_group=groups("prog_timemean","atmo_timemean"))

    cf_desc    = t_cf_var('surface_air_pressure', 'Pa', 'surface pressure (time mean)', dataType)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'ps_m', p_acc%pres_sfc,                                 &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                &    in_group=groups("prog_timemean","atmo_timemean"),          &
                & ldims=shape2d_c)

    p_acc%l_pres_msl = .FALSE.
    IF (PRESENT(l_pres_msl)) THEN
      IF (l_pres_msl) THEN
        cf_desc    = t_cf_var('mean sea level pressure', 'Pa',                  &
          &                   'mean sea level pressure (time mean)', dataType)
        grib2_desc = t_grib2_var(0, 3, 1, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( list, 'psl_m', p_acc%pres_msl,                            &
          &           GRID_UNSTRUCTURED_CELL, ZA_MEANSEA, cf_desc, grib2_desc,  &
          &              in_group=groups("prog_timemean","atmo_timemean"),      &
          &           ldims=shape2d_c)
        p_acc%l_pres_msl = .TRUE.
      END IF
    END IF

    cf_desc    = t_cf_var('air_density', 'kg m-3', 'air density (time mean)', dataType)
    grib2_desc = t_grib2_var(0, 3, 10, DATATYPE_PACK_VAR, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'rho_m', p_acc%rho,  &
      &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
      &           ldims=shape3d_c,                                             &
      &           vert_interp=create_vert_interp_metadata(                     &
      &              vert_intp_type=vintp_types("P","Z","I"),                  &
      &              vert_intp_method=VINTP_METHOD_LIN ),                      &
      &           in_group=groups("prog_timemean","atmo_timemean")) 
    ! }}}
    ! TRACERS {{{
    ! support qv,qc,qi because they area always there
    IF (ntracer > 0) THEN
      cf_desc    = t_cf_var('tracer', 'kg kg-1', 'air tracer (time mean)', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(0,20,2, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( list, 'tracer', p_acc%tracer,                           &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, &
        &           ldims=shape4d_c ,                                       &
        &           lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ALLOCATE(p_acc%tracer_ptr(ntracer))
      DO jt=1,ntracer
        IF (jt == iqv ) CALL add_ref( &
          &  list, 'tracer', 'hus_m', p_acc%tracer_ptr(jt)%p,               &
          &  GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &
          &  t_cf_var('specific_humidity', 'kg kg-1',                       &
          &           'specific_cloud_ice_content (time mean)', dataType),  &
          &  t_grib2_var( 0, 1, 0, ibits, GRID_REFERENCE, GRID_CELL),       &
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

        IF ( jt == iqc )  CALL add_ref(                                     &
          &  list, 'tracer', 'clw_m', p_acc%tracer_ptr(jt)%p,               &
          &  GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &
          &  t_cf_var('specific_cloud_water_content', 'kg kg-1',            &
          &           'specific_cloud_ice_content (time mean)', dataType),  &
          &  t_grib2_var(0, 1, 22, ibits, GRID_REFERENCE, GRID_CELL),       &
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

        IF ( jt == iqi ) CALL add_ref(                                      &
          &  list, 'tracer', 'cli_m', p_acc%tracer_ptr(jt)%p,               &
          &  GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &
          &  t_cf_var('specific_cloud_ice_content', 'kg kg-1',              &
          &           'specific_cloud_ice_content (time mean)', dataType),  &
          &  t_grib2_var(0, 1, 82, ibits, GRID_REFERENCE, GRID_CELL),       &
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
    ! }}}

    ! ECHAM {{{
    IF (echam_forcing_active) THEN
    cf_desc    = t_cf_var('cosmu0', '', 'cosine of the zenith angle (time mean)', dataType)
    grib2_desc = t_grib2_var(192,214,1, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'cosmu0_m', p_acc%cosmu0,                                         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d, &
                & in_group=groups("echam_timemean","atmo_timemean"))

    cf_desc    = t_cf_var('rsdt', 'W m-2',                                                    &
                &         'downward shortwave flux at the top of the atmosphere (time mean)', &
                &         dataType)
    grib2_desc = t_grib2_var(0,4,7, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'rsdt_m', p_acc%flxdwswtoa,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,&
                & in_group=groups("echam_timemean","atmo_timemean"), ldims=shape2d )

    cf_desc    = t_cf_var('clt', 'm2 m-2', &
               & 'total cloud cover (time mean)', dataType)
    grib2_desc = t_grib2_var(0,6,1, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'clt_m', p_acc%aclcov,                   &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d, in_group=groups("echam_timemean","atmo_timemean"),&
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('prlr', 'kg m-2 s-1',    &
               & 'large-scale precipitation flux (water) (time mean)', dataType)
    grib2_desc = t_grib2_var(0,1,77, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'prlr_m', p_acc%rsfl,                    &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"), &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('prcr', 'kg m-2 s-1',    &
               & 'convective precipitation flux (water) (time mean)', dataType)
    grib2_desc = t_grib2_var(0,1,76, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'prcr_m', p_acc%rsfc,                    &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),&
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('prls', 'kg m-2 s-1',    &
               & 'large-scale precipitation flux (snow) (time mean)', dataType)
    grib2_desc = t_grib2_var(0,1,59, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'prls_m', p_acc%ssfl,                    &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d, in_group=groups("echam_timemean","atmo_timemean"),&
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('prcs', 'kg m-2 s-1',    &
               & 'convective precipitation flux (snow) (time mean)', dataType)
    grib2_desc = t_grib2_var(0,1,58, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'prcs_m', p_acc%ssfc,                    &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d, in_group=groups("echam_timemean","atmo_timemean"),&
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('pr', 'kg m-2 s-1',               &
         &                'precipitation flux (time mean)', &
         &                dataType)
    grib2_desc = t_grib2_var(0, 1, 52, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'pr_m', p_acc%totprec,                   &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d, in_group=groups("echam_timemean","atmo_timemean"),&
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('total_vapour', 'kg m-2', 'vertically integrated water vapour (time mean)', &
         &                dataType)
    grib2_desc = t_grib2_var(0,1,64, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'prw_m', p_acc%qvi,                      &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),&
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('total_cloud_water', 'kg m-2',&
               & 'vertically integrated cloud water (time mean)', dataType)
    grib2_desc = t_grib2_var(0,1,69, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'cllvi_m', p_acc%xlvi,                   &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),&
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('total_cloud_ice', 'kg m-2',&
               & 'vertically integrated cloud ice (time mean)', dataType)
    grib2_desc = t_grib2_var(0,1,70, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'clivi_m', p_acc%xivi,                   &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),&
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('rsns', 'W m-2', ' shortwave net flux at surface (time mean)', dataType)
    grib2_desc = t_grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'rsns_m', p_acc%swflxsfc,                &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),&
         &        isteptype=TSTEP_INSTANT )
        
    cf_desc    = t_cf_var('rsnt', 'W m-2', ' shortwave net flux at TOA (time mean)', dataType)
    grib2_desc = t_grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'rsnt_m', p_acc%swflxtoa,                &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),&
         &        isteptype=TSTEP_INSTANT )
        
    cf_desc    = t_cf_var('rlns', 'W m-2', 'longwave net flux at surface (time mean)', dataType)
    grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'rlns_m', p_acc%lwflxsfc,                &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),&
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('rlnt', 'W m-2', 'longwave net flux at TOA (time mean)', dataType)
    grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'rlnt_m', p_acc%lwflxtoa,&
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),&
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('surface_temperature', '', 'surface temperature (time mean)', dataType)
    grib2_desc = t_grib2_var(0,0,0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( list, 'ts_m', p_acc%tsfc,                      &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
              & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"))

    CALL add_var( list, 'evspsbl_m', p_acc%evap,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
                & t_cf_var('evap', 'kg m-2 s-1', 'evaporation (time mean)', &
                & dataType),                                     &
                & t_grib2_var(0,1,6,iextbits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),&
                & isteptype=TSTEP_INSTANT                                 )

    CALL add_var( list, 'hfls_m', p_acc%lhflx,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('hfls', 'W m-2 ', 'latent heat flux (time mean)', &
                & dataType),                                        &
                & t_grib2_var(0,0,10, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),&
                & isteptype=TSTEP_INSTANT                                 )

    CALL add_var( list, 'hfss_m', p_acc%shflx,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('hfss', 'W m-2 ', 'sensible heat flux (time mean)', &
                & dataType),                                        &
                & t_grib2_var(0,0,11, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),&
                & isteptype=TSTEP_INSTANT                                 )

    CALL add_var( list, 'tauu_m', p_acc%u_stress        ,       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('u_stress', 'N m-2', 'u-momentum flux at the surface (time mean)',          &
                &          dataType),                                     &
                & t_grib2_var(0,2,17, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),&
                & isteptype=TSTEP_INSTANT )

    CALL add_var( list, 'tauv_m', p_acc%v_stress,               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('v_stress', 'N m-2', 'v-momentum flux at the surface (time mean)',          &
                &          dataType),                                     &
                & t_grib2_var(0,2,18, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d,in_group=groups("echam_timemean","atmo_timemean"),&
                & isteptype=TSTEP_INSTANT )
    END IF
    ! }}}

    p_acc%numberOfAccumulations = 0
  END SUBROUTINE construct_opt_acc


  SUBROUTINE update_opt_acc(acc, nh_prog, rho, nh_diag, subset, levels, echam_forcing_active)
    TYPE(t_nh_acc),  INTENT(INOUT)   :: acc
    TYPE(t_nh_prog), INTENT(IN)      :: nh_prog
    REAL(wp), INTENT(IN)             :: rho(:,:,:)
    TYPE(t_nh_diag), INTENT(IN)      :: nh_diag
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER , INTENT(IN)             :: levels
    LOGICAL , INTENT(IN)             :: echam_forcing_active

    INTEGER :: jt

    !WRITE(message_text,'(a,i2)') '(pre ): numberOfAccumulations:',acc%numberOfAccumulations
    CALL message('update_opt_nh_acc', TRIM(message_text))
    CALL add_fields(acc%u       , nh_diag%u       , subset, levels=levels)
    CALL add_fields(acc%v       , nh_diag%v       , subset, levels=levels)
    CALL add_fields(acc%temp    , nh_diag%temp    , subset, levels=levels)
    CALL add_fields(acc%pres_sfc, nh_diag%pres_sfc, subset)
    IF (acc%l_pres_msl) &
      & CALL add_fields(acc%pres_msl, nh_diag%pres_msl, subset)
    CALL add_fields(acc%rho     , rho             , subset, levels=levels)
    CALL dbg_print('RHO Update FROM',nh_prog%rho  ,'opt_diag',5, in_subset=subset)
    CALL dbg_print('RHO Update TO  ',acc%rho      ,'opt_diag',5, in_subset=subset)
    IF (ntracer > 0) THEN
      DO jt=1,ntracer
        CALL add_fields(acc%tracer(:,:,:,jt) ,nh_prog%tracer(:,:,:,jt),subset,levels=levels)
        CALL dbg_print('Tracer Update FROM',nh_prog%tracer(:,:,:,jt),'opt_diag',5, in_subset=subset)
        CALL dbg_print('Tracer Update TO  ',acc%tracer(:,:,:,jt),    'opt_diag',5, in_subset=subset)
      END DO
    ENDIF

    IF (echam_forcing_active) CALL update_opt_acc_echam(acc, 1, subset)

    acc%numberOfAccumulations = acc%numberOfAccumulations + 1
    !WRITE(message_text,'(a,i2)') '(post): numberOfAccumulations:',acc%numberOfAccumulations
    !CALL message('update_opt_nh_acc', TRIM(message_text))
  END SUBROUTINE update_opt_acc


  SUBROUTINE reset_opt_acc(acc,echam_forcing_active)
    TYPE(t_nh_acc) :: acc!(n_dom)
    LOGICAL , INTENT(IN)             :: echam_forcing_active

    INTEGER :: jt

    acc%u                            = 0.0_wp
    acc%v                            = 0.0_wp
    acc%temp                         = 0.0_wp
    acc%pres_sfc                     = 0.0_wp
    IF (acc%l_pres_msl) acc%pres_msl = 0.0_wp
    acc%rho                          = 0.0_wp
    !WRITE(message_text,'(a,i2)') '(    ): numberOfAccumulations:',acc%numberOfAccumulations
    !CALL message('reset_opt_nh_acc', TRIM(message_text))
    IF (ntracer > 0) THEN
      DO jt=1,ntracer
        acc%tracer(:,:,:,jt) = 0.0_wp
      END DO
    ENDIF
    IF (echam_forcing_active) THEN
      acc%cosmu0     = 0.0_wp
      acc%flxdwswtoa = 0.0_wp
      acc%aclcov     = 0.0_wp
      acc%rsfl       = 0.0_wp
      acc%rsfc       = 0.0_wp
      acc%ssfl       = 0.0_wp
      acc%ssfc       = 0.0_wp
      acc%totprec    = 0.0_wp
      acc%qvi        = 0.0_wp
      acc%xlvi       = 0.0_wp
      acc%xivi       = 0.0_wp
      acc%swflxsfc   = 0.0_wp
      acc%swflxtoa   = 0.0_wp
      acc%lwflxsfc   = 0.0_wp
      acc%lwflxtoa   = 0.0_wp
      acc%tsfc       = 0.0_wp
      acc%evap       = 0.0_wp
      acc%lhflx      = 0.0_wp
      acc%shflx      = 0.0_wp
      acc%u_stress   = 0.0_wp
      acc%v_stress   = 0.0_wp
    ENDIF
    acc%numberOfAccumulations = 0
  END SUBROUTINE reset_opt_acc


  SUBROUTINE calc_mean_opt_acc(acc,echam_forcing_active)
    TYPE(t_nh_acc) :: acc!(n_dom)
    LOGICAL , INTENT(IN)             :: echam_forcing_active

    INTEGER :: jt

    acc%u        = acc%u        /REAL(acc%numberOfAccumulations,wp)
    acc%v        = acc%v        /REAL(acc%numberOfAccumulations,wp)
    acc%temp     = acc%temp     /REAL(acc%numberOfAccumulations,wp)
    acc%pres_sfc = acc%pres_sfc /REAL(acc%numberOfAccumulations,wp)
    IF (acc%l_pres_msl) acc%pres_msl = acc%pres_msl /REAL(acc%numberOfAccumulations,wp)
    acc%rho      = acc%rho      /REAL(acc%numberOfAccumulations,wp)
    !WRITE(message_text,'(a,i2)') '(    ): numberOfAccumulations:',acc%numberOfAccumulations
    !CALL message('calc_mean_opt_nh_acc', TRIM(message_text))
    IF (ntracer > 0) THEN
      DO jt=1,ntracer
        acc%tracer(:,:,:,jt) = acc%tracer(:,:,:,jt) /REAL(acc%numberOfAccumulations,wp)
      END DO
    ENDIF
    IF (echam_forcing_active) THEN
      acc%cosmu0     = acc%cosmu0     /REAL(acc%numberOfAccumulations,wp)
      acc%flxdwswtoa = acc%flxdwswtoa /REAL(acc%numberOfAccumulations,wp)
      acc%aclcov     = acc%aclcov     /REAL(acc%numberOfAccumulations,wp)
      acc%rsfl       = acc%rsfl       /REAL(acc%numberOfAccumulations,wp)
      acc%rsfc       = acc%rsfc       /REAL(acc%numberOfAccumulations,wp)
      acc%ssfl       = acc%ssfl       /REAL(acc%numberOfAccumulations,wp)
      acc%ssfc       = acc%ssfc       /REAL(acc%numberOfAccumulations,wp)
      acc%totprec    = acc%totprec    /REAL(acc%numberOfAccumulations,wp)
      acc%qvi        = acc%qvi        /REAL(acc%numberOfAccumulations,wp)
      acc%xlvi       = acc%xlvi       /REAL(acc%numberOfAccumulations,wp)
      acc%xivi       = acc%xivi       /REAL(acc%numberOfAccumulations,wp)
      acc%swflxsfc   = acc%swflxsfc   /REAL(acc%numberOfAccumulations,wp)
      acc%swflxtoa   = acc%swflxtoa   /REAL(acc%numberOfAccumulations,wp)
      acc%lwflxsfc   = acc%lwflxsfc   /REAL(acc%numberOfAccumulations,wp)
      acc%lwflxtoa   = acc%lwflxtoa   /REAL(acc%numberOfAccumulations,wp)
      acc%tsfc       = acc%tsfc       /REAL(acc%numberOfAccumulations,wp)
      acc%evap       = acc%evap       /REAL(acc%numberOfAccumulations,wp)
      acc%lhflx      = acc%lhflx      /REAL(acc%numberOfAccumulations,wp)
      acc%shflx      = acc%shflx      /REAL(acc%numberOfAccumulations,wp)
      acc%u_stress   = acc%u_stress   /REAL(acc%numberOfAccumulations,wp)
      acc%v_stress   = acc%v_stress   /REAL(acc%numberOfAccumulations,wp)
    ENDIF
  END SUBROUTINE calc_mean_opt_acc

  SUBROUTINE update_opt_acc_echam(acc,jg,subset)
    TYPE(t_nh_acc) :: acc
    INTEGER :: jg
    TYPE(t_subset_range) , INTENT(IN):: subset

    CALL add_fields(acc%cosmu0          , prm_field(jg)%cosmu0    , subset)
    CALL add_fields(acc%flxdwswtoa      , prm_field(jg)%flxdwswtoa, subset)
    CALL add_fields(acc%aclcov          , prm_field(jg)%aclcov    , subset)
    CALL add_fields(acc%rsfl            , prm_field(jg)%rsfl      , subset)
    CALL add_fields(acc%rsfc            , prm_field(jg)%rsfc      , subset)
    CALL add_fields(acc%ssfl            , prm_field(jg)%ssfl      , subset)
    CALL add_fields(acc%ssfc            , prm_field(jg)%ssfc      , subset)
    CALL add_fields(acc%totprec         , prm_field(jg)%totprec   , subset)
    CALL add_fields(acc%qvi             , prm_field(jg)%qvi       , subset)
    CALL add_fields(acc%xlvi            , prm_field(jg)%xlvi      , subset)
    CALL add_fields(acc%xivi            , prm_field(jg)%xivi      , subset)
    CALL add_fields(acc%swflxsfc        , prm_field(jg)%swflxsfc  , subset)
    CALL add_fields(acc%swflxtoa        , prm_field(jg)%swflxtoa  , subset)
    CALL add_fields(acc%lwflxsfc        , prm_field(jg)%lwflxsfc  , subset)
    CALL add_fields(acc%lwflxtoa        , prm_field(jg)%lwflxtoa  , subset)
    CALL add_fields(acc%tsfc            , prm_field(jg)%tsfc      , subset)
    CALL add_fields(acc%evap            , prm_field(jg)%evap      , subset)
    CALL add_fields(acc%lhflx           , prm_field(jg)%lhflx     , subset)
    CALL add_fields(acc%shflx           , prm_field(jg)%shflx     , subset)
    CALL add_fields(acc%u_stress        , prm_field(jg)%u_stress  , subset)
    CALL add_fields(acc%v_stress        , prm_field(jg)%v_stress  , subset)
  END SUBROUTINE update_opt_acc_echam

  !-------------
  !
  !> Add optional diagnostic variable lists (might remain empty)
  !
  SUBROUTINE construct_opt_diag(p_patch, l_init_pz, echam_forcing_active,l_pres_msl)
    TYPE(t_patch),        INTENT(IN)   :: p_patch(n_dom)
    LOGICAL,              INTENT(IN)   :: l_init_pz
    LOGICAL, INTENT(IN)                :: echam_forcing_active
    LOGICAL , OPTIONAL, INTENT(IN)     :: l_pres_msl
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
    CALL construct_opt_acc(p_patch(1), &
      &                    p_nh_opt_diag(1)%opt_acc_list, &
      &                    p_nh_opt_diag(1)%acc,echam_forcing_active, &
      &                    l_pres_msl=l_pres_msl)
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

    CHARACTER(*), PARAMETER :: routine = TRIM("mo_opt_diagnostics:vcoeff_allocate")

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



