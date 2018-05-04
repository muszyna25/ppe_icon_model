!>
!! Contains basic diagnostics for ICON ocean model.
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/02)
!!  Extended   by Stephan Lorenz,   MPI-M (2012)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_diagnostics
  USE mo_kind,               ONLY: wp, dp, i8
#ifdef _OPENMP
  USE omp_lib
#endif
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range, t_subset_indexed
  USE mo_grid_tools,         ONLY: get_oriented_edges_from_global_vertices, check_global_indexes
  USE mo_mpi,                ONLY: my_process_is_stdio, p_field_sum, get_my_mpi_work_id, &
    & p_comm_work_test, p_comm_work, p_io, p_bcast, my_process_is_mpi_workroot
  USE mo_sync,               ONLY: global_sum_array, disable_sync_checks, enable_sync_checks, &
    &                              sync_c, sync_e, sync_patch_array
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates, cvec2gvec
  USE mo_advection_utils,    ONLY: laxfr_upflux
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  USE mo_dbg_nml,            ONLY: idbg_val
  USE mo_math_constants,     ONLY: rad2deg, dbl_eps
  USE mo_impl_constants,     ONLY: sea_boundary,sea, &
    & min_rlcell, min_rledge, min_rlcell, &
    & max_char_length, min_dolic
  USE mo_timer,              ONLY: timer_calc_moc, timer_start, timer_stop
  USE mo_cdi_constants,      ONLY: GRID_EDGE, GRID_CELL, GRID_UNSTRUCTURED_EDGE, &
    & GRID_UNSTRUCTURED_CELL
  USE mo_ocean_nml,          ONLY: n_zlev, no_tracer, &
    & gibraltar, &
    & denmark_strait, &
    & drake_passage, &
    & florida_strait, &
    & indonesian_throughflow,&
    & scotland_iceland, &
    & mozambique, &
    & framStrait, &
    & beringStrait, &
    & barentsOpening, &
    & agulhas, &
    & agulhas_long, &
    & agulhas_longer, &
    & ab_const, ab_beta, ab_gam, iswm_oce, discretization_scheme, &
    & iforc_oce, No_Forcing, i_sea_ice, diagnostics_level, &
    & diagnose_for_horizontalVelocity, OceanReferenceDensity
  USE mo_sea_ice_nml,        ONLY: kice
  USE mo_dynamics_config,    ONLY: nold,nnew
  USE mo_parallel_config,    ONLY: nproma, p_test_run
  USE mo_run_config,         ONLY: dtime, nsteps
  USE mo_physical_constants, ONLY: grav, rhos, rhoi,sice
  USE mo_model_domain,       ONLY: t_patch, t_patch_3d,t_patch_vert, t_grid_edges
  USE mo_ocean_types,        ONLY: t_hydro_ocean_state, t_hydro_ocean_diag
  USE mo_ocean_diagnostics_types,  ONLY: t_ocean_regions, t_ocean_region_volumes, &
    &  t_ocean_region_areas, t_ocean_monitor
  USE mo_hamocc_types,       ONLY: t_hamocc_state
  USE mo_hamocc_diagnostics, ONLY: get_monitoring 
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_exception,          ONLY: message, finish, message_text, warning
  USE mo_sea_ice_types,      ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_ocean_surface_types,ONLY: t_ocean_surface
  USE mo_linked_list,        ONLY: t_var_list
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_scalar_product,     ONLY: map_edges2cell_3d
  USE mo_io_units,           ONLY: find_next_free_unit
  USE mo_util_file,          ONLY: util_symlink, util_rename, util_islink, util_unlink
  USE mo_statistics,         ONLY: subset_sum, levels_horizontal_mean, total_mean, gather_sums
  USE mo_fortran_tools,      ONLY: assign_if_present

  USE mo_linked_list,        ONLY: t_var_list
  USE mo_var_list,           ONLY: add_var,                  &
    &                              new_var_list,             &
    &                              delete_var_list,          &
    &                              default_var_list_settings,&
    &                              add_ref
  USE mo_var_metadata,       ONLY: groups
  USE mo_cf_convention
  USE mo_grib2,              ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, GRID_UNSTRUCTURED
  USE mo_cdi_constants,      ONLY: GRID_EDGE, GRID_CELL, GRID_UNSTRUCTURED_EDGE, &
    &                              GRID_UNSTRUCTURED_CELL
  USE mo_zaxis_type,         ONLY: ZA_DEPTH_BELOW_SEA
  USE mo_mpi,                ONLY: my_process_is_mpi_parallel, p_sum
  USE mo_io_config,          ONLY: lnetcdf_flt64_output
  USE mo_name_list_output_init, ONLY: isRegistered

  USE mtime,                 ONLY: datetime, MAX_DATETIME_STR_LEN, datetimeToPosixString

  IMPLICIT NONE

  !PRIVATE

  CHARACTER(LEN=12)           :: str_module    = 'oceDiag     '  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

  INTEGER :: moc_unit  = -1 ! file handle for the global timeseries output
  CHARACTER(LEN=max_char_length) :: diag_fname, moc_fname
  INTEGER, PARAMETER :: linecharacters  = 2048

  !
  ! PUBLIC INTERFACE
  !
  PUBLIC :: calc_slow_oce_diagnostics, calc_fast_oce_diagnostics
  PUBLIC :: construct_oce_diagnostics
  PUBLIC :: destruct_oce_diagnostics
  PUBLIC :: calc_moc
  PUBLIC :: calc_psi

  INTERFACE calc_moc
    MODULE PROCEDURE calc_moc_acc
    MODULE PROCEDURE calc_moc_internal
  END INTERFACE

  TYPE t_oce_section
    TYPE(t_subset_indexed) :: subset
    REAL(wp), POINTER :: orientation(:)
  END TYPE t_oce_section

  INTEGER, PARAMETER  :: oce_section_count = 13
  PRIVATE             :: oce_section_count
  TYPE(t_oce_section) :: oce_sections(oce_section_count)
  PRIVATE             :: oce_sections

  TYPE(t_ocean_region_volumes),SAVE :: ocean_region_volumes
  PRIVATE                           :: ocean_region_volumes
  TYPE(t_ocean_region_areas),SAVE   :: ocean_region_areas
  PRIVATE                           :: ocean_region_areas


  TYPE(t_var_list) :: horizontal_velocity_diagnostics
  ! addtional diagnostics
  REAL(wp), POINTER :: veloc_adv_horz_u(:,:,:),  veloc_adv_horz_v(:,:,:), &
    & laplacian_horz_u(:,:,:), laplacian_horz_v(:,:,:), vn_u(:,:,:), vn_v(:,:,:), &
    & mass_flx_e_u(:,:,:), mass_flx_e_v(:,:,:), pressure_grad_u(:,:,:), pressure_grad_v(:,:,:), &
    & potential_vort_e(:,:,:), potential_vort_c(:,:,:)

 CHARACTER(LEN=*), PARAMETER :: module_name="mo_ocean_statistics"

CONTAINS

  !-------------------------------------------------------------------------
  !  The constructor of the types related to ocean diagnostics
  !>
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !!
!<Optimize:inUse>
  SUBROUTINE construct_oce_diagnostics( patch_3D, ocean_state, datestring )
    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3D
    TYPE(t_hydro_ocean_state), TARGET      :: ocean_state
    CHARACTER(LEN=32)                      :: datestring

    !local variable
    INTEGER :: i,ist
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_ocean_diagnostics:construct_oce_diagnostics')
    !-----------------------------------------------------------------------
    INTEGER  :: nblks_e,blockNo,jc,jk, region_index,start_index,end_index
    REAL(wp) :: surface_area, surface_height, prism_vol, prism_area, column_volume

    TYPE(t_patch), POINTER        :: patch_2d
    TYPE(t_subset_range), POINTER :: owned_cells
    INTEGER, POINTER              :: regions(:,:)
    TYPE(t_ocean_regions)         :: ocean_regions
    CHARACTER(LEN=max_char_length) :: listname
    INTEGER                       :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    !-----------------------------------------------------------------------
    patch_2d => patch_3D%p_patch_2d(1)
    regions => patch_3D%regio_c
    !-----------------------------------------------------------------------
    owned_cells => patch_2d%cells%owned
    nblks_e = patch_2d%nblks_e
    !-----------------------------------------------------------------------
    WRITE(listname,'(a)')  'horizontal_velocity_diagnostics'
    CALL new_var_list(horizontal_velocity_diagnostics, listname, patch_id=patch_2d%id)
    CALL default_var_list_settings( horizontal_velocity_diagnostics,            &
      & lrestart=.FALSE.,model_type='oce',loutput=.TRUE. )
    !-----------------------------------------------------------------------
    IF (diagnose_for_horizontalVelocity) THEN
      CALL add_var(horizontal_velocity_diagnostics, 'veloc_adv_horz_u', veloc_adv_horz_u, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('veloc_adv_horz_u','m/s','velocity advection zonal', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'veloc_adv_horz_v', veloc_adv_horz_v, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('veloc_adv_horz_v','m/s','velocity advection meridional', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'laplacian_horz_u', laplacian_horz_u, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('laplacian_horz_u','m/s','velocity laplacian zonal', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'laplacian_horz_v', laplacian_horz_v, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('laplacian_horz_v','m/s','velocity laplacian meridional', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'vn_u', vn_u, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('vn_u','m/s','edge velocity zonal', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'vn_v', vn_v, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('vn_v','m/s','edge velocity meridional', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'mass_flx_e_u', mass_flx_e_u, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('mass_flx_e_u','m*m/s','mass flux zonal', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'mass_flx_e_v', mass_flx_e_v, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('mass_flx_e_v','m*m/s','mass flux meridional', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'pressure_grad_u', pressure_grad_u, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('pressure_grad_u','N','pressure gradient zonal', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'pressure_grad_v', pressure_grad_v, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('pressure_grad_v','N','pressure gradient meridional', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'potential_vort_e', potential_vort_e, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('vn_v','1/s','potential vorticity at edges', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'potential_vort_c', potential_vort_c, &
        & grid_unstructured_cell, za_depth_below_sea, &
        & t_cf_var('vn_v','1/s','potential vorticity at cells', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,patch_2d%alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

    ENDIF

    IF (diagnostics_level < 1) RETURN

    CALL message (TRIM(routine), 'start')

    ! CALL check_global_indexes(patch_2d)

    ! compute subsets for given sections path allong edges
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(1)%subset,      &
      & orientation = oce_sections(1)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array = gibraltar,            &
      & subset_name = 'gibraltar')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(2)%subset,      &
      & orientation = oce_sections(2)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =denmark_strait,        &
      & subset_name = 'denmark_strait')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(3)%subset,      &
      & orientation = oce_sections(3)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =drake_passage,         &
      & subset_name = 'drake_passage')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(4)%subset,      &
      & orientation = oce_sections(4)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =indonesian_throughflow,&
      & subset_name = 'indonesian_throughflow')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(5)%subset,      &
      & orientation = oce_sections(5)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =scotland_iceland,      &
      & subset_name = 'scotland_iceland')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(6)%subset,      &
      & orientation = oce_sections(6)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =mozambique,      &
      & subset_name = 'mozambique')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(7)%subset,      &
      & orientation = oce_sections(7)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =framStrait,      &
      & subset_name = 'framStrait')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(8)%subset,      &
      & orientation = oce_sections(8)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =beringStrait,      &
      & subset_name = 'beringStrait')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(9)%subset,      &
      & orientation = oce_sections(9)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =barentsOpening,      &
      & subset_name = 'barentsOpening')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(10)%subset,      &
      & orientation = oce_sections(10)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =agulhas,      &
      & subset_name = 'agulhas')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(11)%subset,      &
      & orientation = oce_sections(11)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =agulhas_long,      &
      & subset_name = 'agulhas_long')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(12)%subset,      &
      & orientation = oce_sections(12)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =agulhas_longer,      &
      & subset_name = 'agulhas_longer')

    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(13)%subset,       &
      & orientation = oce_sections(13)%orientation, &
      & patch_3d = patch_3D,                         &
      & global_vertex_array =florida_strait,         &
      & subset_name = 'florida_strait')
!     CALL finish("","")

    surface_area   = 0.0_wp
    surface_height = 0.0_wp
    prism_vol      = 0.0_wp
    prism_area     = 0.0_wp
    ocean_region_areas%total = 0.0_wp
    ! compute regional ocean volumes
    DO blockNo = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index

        ! area
        prism_area               = patch_2d%cells%area(jc,blockNo)
        ocean_region_areas%total = ocean_region_areas%total + prism_area

        ! volume
        CALL compute_vertical_volume(blockNo,jc, &
          & prism_area, &
          & ocean_state%p_prog(nnew(1))%h(jc,blockNo), &
          & patch_3D%p_patch_1d(1)%prism_thick_c(jc,:,blockNo), &
          & patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo), &
          & column_volume)
        ocean_region_volumes%total = ocean_region_volumes%total + column_volume

        region_index = regions(jc,blockNo)
        IF (ocean_regions%greenland_iceland_norwegian_sea == region_index) THEN
          ocean_region_volumes%greenland_iceland_norwegian_sea = &
            & ocean_region_volumes%greenland_iceland_norwegian_sea + column_volume
          ocean_region_areas%greenland_iceland_norwegian_sea   = &
            & ocean_region_areas%greenland_iceland_norwegian_sea + prism_area

        ELSEIF (ocean_regions%arctic_ocean == region_index) THEN
          ocean_region_volumes%arctic_ocean                    = ocean_region_volumes%arctic_ocean + column_volume
          ocean_region_areas%arctic_ocean                      = ocean_region_areas%arctic_ocean + prism_area

        ELSEIF (ocean_regions%labrador_sea == region_index) THEN
          ocean_region_volumes%labrador_sea                    = ocean_region_volumes%labrador_sea + column_volume
          ocean_region_areas%labrador_sea                      = ocean_region_areas%labrador_sea + prism_area

        ELSEIF (ocean_regions%north_atlantic == region_index) THEN
          ocean_region_volumes%north_atlantic                  = ocean_region_volumes%north_atlantic + column_volume
          ocean_region_areas%north_atlantic                    = ocean_region_areas%north_atlantic + prism_area

        ELSEIF (ocean_regions%tropical_atlantic == region_index) THEN
          ocean_region_volumes%tropical_atlantic               = ocean_region_volumes%tropical_atlantic + column_volume
          ocean_region_areas%tropical_atlantic                 = ocean_region_areas%tropical_atlantic + prism_area

        ELSEIF (ocean_regions%southern_ocean == region_index) THEN
          ocean_region_volumes%southern_ocean                  = ocean_region_volumes%southern_ocean + column_volume
          ocean_region_areas%southern_ocean                    = ocean_region_areas%southern_ocean + prism_area

        ELSEIF (ocean_regions%indian_ocean == region_index) THEN
          ocean_region_volumes%indian_ocean                    = ocean_region_volumes%indian_ocean + column_volume
          ocean_region_areas%indian_ocean                      = ocean_region_areas%indian_ocean + prism_area

        ELSEIF (ocean_regions%tropical_pacific == region_index) THEN
          ocean_region_volumes%tropical_pacific                = ocean_region_volumes%tropical_pacific + column_volume
          ocean_region_areas%tropical_pacific                  = ocean_region_areas%tropical_pacific + prism_area

        ELSEIF (ocean_regions%north_pacific == region_index) THEN
          ocean_region_volumes%north_pacific                   = ocean_region_volumes%north_pacific + column_volume
          ocean_region_areas%north_pacific                     = ocean_region_areas%north_pacific + prism_area

        ELSEIF (ocean_regions%caribbean == region_index) THEN
          ocean_region_volumes%caribbean                       = ocean_region_volumes%caribbean + column_volume
          ocean_region_areas%caribbean                         = ocean_region_areas%caribbean + prism_area
        END IF

      END DO
    END DO
    ! compute global values

    CALL disable_sync_checks()
    ocean_region_volumes%total = global_sum_array(ocean_region_volumes%total)
    ocean_region_areas%total   = global_sum_array(ocean_region_areas%total)
    CALL enable_sync_checks()

    CALL message (TRIM(routine), 'end')
  END SUBROUTINE construct_oce_diagnostics
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !  !The destructor of the types related to ocean diagnostics
  !>
  !!
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_oce_diagnostics()
    !
    !
    !local variables
    INTEGER :: i,iret
    CHARACTER(LEN=max_char_length)  :: linkname
    CHARACTER(LEN=max_char_length)  :: message_text

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_ocean_diagnostics:destruct_oce_diagnostics')
    !-----------------------------------------------------------------------
    CALL delete_var_list(horizontal_velocity_diagnostics)

    IF (diagnostics_level <= 0) RETURN

    CALL message (TRIM(routine), 'end')
  END SUBROUTINE destruct_oce_diagnostics
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE compute_vertical_volume(blockNo,jc,prism_area,surface_height,thicknesses,max_vertical_level,volume)
    INTEGER,  INTENT(in)  :: blockNo,jc,max_vertical_level
    REAL(wp), INTENT(in)  :: prism_area, surface_height, thicknesses(:)
    REAL(wp), INTENT(inout) :: volume

    INTEGER :: jk
    REAL(wp) :: surface_height_,prism_vol_

    volume  = 0.0_wp
    DO jk = 1,max_vertical_level
      !local volume
      surface_height_ = MERGE(surface_height,0.0_wp, 1 == jk)
      prism_vol_      = prism_area * (thicknesses(jk) + surface_height_)
      !Fluid volume wrt lsm
      volume          = volume + prism_vol_
    END DO
  END SUBROUTINE compute_vertical_volume
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  ! !  calculate_oce_diagnostics
  !
  ! @par Revision History
  ! Developed  by  Peter Korn, MPI-M (2010).
  !
  SUBROUTINE calc_slow_oce_diagnostics(patch_3D, ocean_state, surface_fluxes, ice, &
    & timestep, this_datetime)
    TYPE(t_patch_3D ),TARGET, INTENT(in)    :: patch_3D
    TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_ocean_surface),    INTENT(in)    :: surface_fluxes
    TYPE (t_sea_ice),   INTENT(in)          :: ice
    INTEGER, INTENT(in) :: timestep
    TYPE(datetime), POINTER                 :: this_datetime

    !Local variables
    INTEGER :: start_cell_index, end_cell_index!,i_startblk_c, i_endblk_c,
    INTEGER :: start_edge_index, end_edge_index !,i_startblk_c, i_endblk_c,
    INTEGER :: jk,jc,je,blockNo, level !,je
    INTEGER :: edge_1_of_cell_idx, edge_1_of_cell_blk
    INTEGER :: edge_2_of_cell_idx, edge_2_of_cell_blk
    INTEGER :: edge_3_of_cell_idx, edge_3_of_cell_blk
    INTEGER :: i_no_t, i
    REAL(wp):: prism_vol, surface_height, prism_area, surface_area, z_w
!   REAL(wp):: ssh_global_mean !TODO
    INTEGER :: reference_timestep
    TYPE(t_patch), POINTER :: patch_2d
    REAL(wp) :: sflux

    TYPE(t_subset_range), POINTER :: owned_cells, edges_inDomain
    TYPE(t_ocean_monitor),  POINTER :: monitor
    CHARACTER(LEN=linecharacters) :: line, nvars
    CHARACTER(LEN=linecharacters) :: fmt_string, real_fmt
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: datestring
    REAL(wp), PARAMETER :: equator = 0.00001_wp
    TYPE(t_ocean_regions)         :: ocean_regions
    CHARACTER(len=32) :: fmtstr

    !-----------------------------------------------------------------------
    patch_2d       => patch_3D%p_patch_2d(1)
    owned_cells    => patch_2d%cells%owned
    edges_inDomain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------

    IF (diagnose_for_horizontalVelocity) THEN
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_edge_index,end_edge_index, je, level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_inDomain%start_block,edges_inDomain%end_block
        CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
        DO je =  start_edge_index, end_edge_index
          DO level=1,patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

            veloc_adv_horz_u(je, level, blockNo) = &
              & ocean_state%p_diag%veloc_adv_horz(je, level, blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v1
            veloc_adv_horz_v(je, level, blockNo) = &
              & ocean_state%p_diag%veloc_adv_horz(je, level, blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v2
            laplacian_horz_u(je, level, blockNo) = &
              & ocean_state%p_diag%laplacian_horz(je, level, blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v1
            laplacian_horz_v(je, level, blockNo) = &
              & ocean_state%p_diag%laplacian_horz(je, level, blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v2
            vn_u(je, level, blockNo) = &
              & ocean_state%p_prog(nnew(1))%vn(je, level, blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v1
            vn_v(je, level, blockNo) = &
              & ocean_state%p_prog(nnew(1))%vn(je, level, blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v2
            mass_flx_e_u(je, level, blockNo) = &
              & ocean_state%p_diag%mass_flx_e(je, level, blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v1
            mass_flx_e_v(je, level, blockNo) = &
              & ocean_state%p_diag%mass_flx_e(je, level, blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v2
            pressure_grad_u(je, level, blockNo) = &
              & ocean_state%p_diag%press_grad(je, level, blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v1
            pressure_grad_v(je, level, blockNo) = &
              & ocean_state%p_diag%press_grad(je, level, blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v2

          ENDDO
        ENDDO
      ENDDO
!ICON_OMP_END_DO

!ICON_OMP_DO PRIVATE(start_edge_index,end_edge_index, je, level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_inDomain%start_block,edges_inDomain%end_block
        CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
        DO je =  start_edge_index, end_edge_index
          DO level=1,patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)
            ! this is when it is calculated by the veloc_adv_horz_mimetic_rot
            potential_vort_e(je, level, blockNo) = ocean_state%p_diag%veloc_adv_horz(je, level, blockNo)
          ENDDO
        ENDDO
      ENDDO
!ICON_OMP_END_DO


!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, jc, level, &
!ICON_OMP edge_1_of_cell_idx, edge_1_of_cell_blk, edge_2_of_cell_idx, edge_2_of_cell_blk, &
!ICON_OMP edge_3_of_cell_idx, edge_3_of_cell_blk) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = owned_cells%start_block, owned_cells%end_block
        CALL get_index_range(owned_cells, blockNo, start_cell_index, end_cell_index)
        DO jc =  start_cell_index, end_cell_index
          edge_1_of_cell_idx  = patch_2d%cells%edge_idx(jc,blockNo,1)
          edge_1_of_cell_blk  = patch_2d%cells%edge_blk(jc,blockNo,1)
          edge_2_of_cell_idx  = patch_2d%cells%edge_idx(jc,blockNo,2)
          edge_2_of_cell_blk  = patch_2d%cells%edge_blk(jc,blockNo,2)
          edge_3_of_cell_idx  = patch_2d%cells%edge_idx(jc,blockNo,3)
          edge_3_of_cell_blk  = patch_2d%cells%edge_blk(jc,blockNo,3)

          DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)

            potential_vort_c(jc, level, blockNo) = &
              & (potential_vort_e(edge_1_of_cell_idx,level,edge_1_of_cell_blk)  &
              & +potential_vort_e(edge_2_of_cell_idx,level,edge_2_of_cell_blk)  &
              & +potential_vort_e(edge_3_of_cell_idx,level,edge_3_of_cell_blk))/3.0_wp

          END DO
        END DO
      END DO
!ICON_OMP_END_PARALLEL

    ENDIF


    !-----------------------------------------------------------------------
    IF (diagnostics_level < 1) RETURN

    monitor        => ocean_state%p_diag%monitor
    surface_area   = 0.0_wp
    surface_height = 0.0_wp
    prism_vol      = 0.0_wp
    prism_area     = 0.0_wp
    z_w            = 0.0_wp

    fmtstr = '%Y-%m-%d %H:%M:%S'
    call datetimeToPosixString(this_datetime, datestring, fmtstr)
    !CALL reset_ocean_monitor(monitor)

    !cell loop to calculate cell based monitored fields volume, kinetic energy and tracer content
    SELECT CASE (iswm_oce)
    CASE (1) ! shallow water mode
      !Potential energy in SW-casep_patch%patch_oce%del_zlev_m(1)
      DO blockNo = owned_cells%start_block,owned_cells%end_block
        CALL get_index_range(owned_cells, blockNo, start_cell_index, end_cell_index)

        DO jc =  start_cell_index, end_cell_index
          IF ( patch_3D%lsm_c(jc,1,blockNo) <= sea_boundary ) THEN
            surface_area = surface_area + patch_2d%cells%area(jc,blockNo)

            prism_vol      = patch_2d%cells%area(jc,blockNo)*patch_3D%p_patch_1d(1)%prism_thick_c(jc,1,blockNo)
            monitor%volume = monitor%volume + prism_vol


!TODO          monitor%pot_energy = monitor%pot_energy &
!TODO            & + 0.5_wp*grav* prism_vol*patch_3D%p_patch_1d(1)%prism_thick_c(jc,1,blockNo)

!TODO          monitor%kin_energy = monitor%kin_energy + ocean_state%p_diag%kin(jc,1,blockNo)*prism_vol

!TODO          monitor%total_energy=monitor%kin_energy+monitor%pot_energy
!             DO i_no_t=1, no_tracer
!               monitor%tracer_content(i_no_t) = monitor%tracer_content(i_no_t)&
!                 & + prism_vol*ocean_state%p_prog(nold(1))%tracer(jc,1,blockNo,i_no_t)
!             END DO
          ENDIF
        END DO
      END DO

    CASE default !3D model
      DO blockNo = owned_cells%start_block, owned_cells%end_block
        CALL get_index_range(owned_cells, blockNo, start_cell_index, end_cell_index)
        !We are dealing with the surface layer first
        DO jc =  start_cell_index, end_cell_index

          ! area
          prism_area   = patch_2d%cells%area(jc,blockNo)
          surface_area = surface_area + prism_area
          ! sum of top layer vertical velocities abolsute values
          monitor%absolute_vertical_velocity = &
            & monitor%absolute_vertical_velocity + ABS(ocean_state%p_diag%w(jc,1,blockNo))*prism_area

          monitor%HeatFlux_ShortWave   = monitor%HeatFlux_ShortWave   + surface_fluxes%HeatFlux_ShortWave(jc,blockNo)*prism_area
          monitor%HeatFlux_LongWave    = monitor%HeatFlux_LongWave    + surface_fluxes%HeatFlux_LongWave (jc,blockNo)*prism_area
          monitor%HeatFlux_Sensible    = monitor%HeatFlux_Sensible    + surface_fluxes%HeatFlux_Sensible (jc,blockNo)*prism_area
          monitor%HeatFlux_Latent      = monitor%HeatFlux_Latent      + surface_fluxes%HeatFlux_Latent   (jc,blockNo)*prism_area
          monitor%FrshFlux_SnowFall    = monitor%FrshFlux_SnowFall    + surface_fluxes%FrshFlux_SnowFall(jc,blockNo)*prism_area
          monitor%FrshFlux_TotalSalt   = monitor%FrshFlux_TotalSalt   + surface_fluxes%FrshFlux_TotalSalt(jc,blockNo)*prism_area
          monitor%FrshFlux_TotalOcean    = monitor%FrshFlux_TotalOcean    + &
            & surface_fluxes%FrshFlux_TotalOcean(jc,blockNo)*prism_area
          monitor%FrshFlux_TotalIce    = monitor%FrshFlux_TotalIce    + surface_fluxes%FrshFlux_TotalIce(jc,blockNo)*prism_area
          monitor%FrshFlux_VolumeIce   = monitor%FrshFlux_VolumeIce   + surface_fluxes%FrshFlux_VolumeIce (jc,blockNo)*prism_area
          monitor%FrshFlux_VolumeTotal  = monitor%FrshFlux_VolumeTotal  + &
            & surface_fluxes%FrshFlux_VolumeTotal(jc,blockNo)*prism_area
          monitor%HeatFlux_Relax = monitor%HeatFlux_Relax + surface_fluxes%HeatFlux_Relax(jc,blockNo)*prism_area
          monitor%FrshFlux_Relax = monitor%FrshFlux_Relax + surface_fluxes%FrshFlux_Relax(jc,blockNo)*prism_area
          monitor%TempFlux_Relax = monitor%TempFlux_Relax + surface_fluxes%TempFlux_Relax(jc,blockNo)*prism_area
          monitor%SaltFlux_Relax = monitor%SaltFlux_Relax + surface_fluxes%SaltFlux_Relax(jc,blockNo)*prism_area

          ! ice budgets
          ! heat
          !
          ! salinity
          ! volume

          DO jk = 1,patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)

            !local volume
            surface_height = MERGE(ocean_state%p_prog(nold(1))%h(jc,blockNo),0.0_wp, 1 == jk)
            prism_vol      = prism_area * (patch_3D%p_patch_1d(1)%prism_thick_c(jc,jk,blockNo) + surface_height)

            !Fluid volume
            monitor%volume = monitor%volume + prism_vol

            !kinetic energy
!TODO           monitor%kin_energy = monitor%kin_energy + ocean_state%p_diag%kin(jc,jk,blockNo)*prism_vol
!TODO
!TODO           !Potential energy
!TODO           IF(jk==1)THEN
!TODO             z_w = (ocean_state%p_diag%w(jc,jk,blockNo)*ocean_state%p_prog(nold(1))%h(jc,blockNo)&
!TODO               & +ocean_state%p_diag%w(jc,jk+1,blockNo)*0.5_wp*patch_3D%p_patch_1d(1)%del_zlev_i(jk))&
!TODO               & /(0.5_wp*patch_3D%p_patch_1d(1)%del_zlev_i(jk)+ocean_state%p_prog(nold(1))%h(jc,blockNo))
!TODO           ELSEIF(jk>1.AND.jk<n_zlev)THEN
!TODO             z_w = (ocean_state%p_diag%w(jc,jk,blockNo)*patch_3D%p_patch_1d(1)%del_zlev_i(jk)&
!TODO               & +ocean_state%p_diag%w(jc,jk+1,blockNo)*patch_3D%p_patch_1d(1)%del_zlev_i(jk+1))&
!TODO               & /(patch_3D%p_patch_1d(1)%del_zlev_i(jk)+patch_3D%p_patch_1d(1)%del_zlev_i(jk+1))
!TODO           ENDIF
!TODO           monitor%pot_energy = monitor%pot_energy + grav*z_w* ocean_state%p_diag%rho(jc,jk,blockNo)* prism_vol

            !Tracer content
!             DO i_no_t=1, no_tracer
!               monitor%tracer_content(i_no_t) = &
!                 & monitor%tracer_content(i_no_t) + prism_vol*ocean_state%p_prog(nold(1))%tracer(jc,jk,blockNo,i_no_t)
!             END DO
          END DO
        END DO
      END DO

    END SELECT

    ! compute global sums {
    CALL disable_sync_checks()
    monitor%volume                     = global_sum_array(monitor%volume)
    surface_area                       = global_sum_array(surface_area)
!TODO    monitor%kin_energy                 = global_sum_array(monitor%kin_energy)/monitor%volume
!TODO    monitor%pot_energy                 = global_sum_array(monitor%pot_energy)/monitor%volume
!TODO    monitor%total_energy               = global_sum_array(monitor%total_energy)/monitor%volume
    monitor%total_salt                 = calc_total_salt_content(patch_2d, &
      &                                                          patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
      &                                                          ice, ocean_state,surface_fluxes,ice%zUnderIce)
!TODO    monitor%vorticity                  = global_sum_array(monitor%vorticity)



!TODO    monitor%potential_enstrophy        = global_sum_array(monitor%potential_enstrophy)
    monitor%absolute_vertical_velocity = global_sum_array(monitor%absolute_vertical_velocity)/surface_area

    IF (iforc_oce > No_Forcing) THEN
      monitor%HeatFlux_ShortWave         = global_sum_array(monitor%HeatFlux_ShortWave)/surface_area
      monitor%HeatFlux_LongWave          = global_sum_array(monitor%HeatFlux_LongWave )/surface_area
      monitor%HeatFlux_Sensible          = global_sum_array(monitor%HeatFlux_Sensible )/surface_area
      monitor%HeatFlux_Latent            = global_sum_array(monitor%HeatFlux_Latent   )/surface_area
      monitor%FrshFlux_SnowFall          = global_sum_array(monitor%FrshFlux_SnowFall)/surface_area
      monitor%FrshFlux_TotalSalt         = global_sum_array(monitor%FrshFlux_TotalSalt)/surface_area
      monitor%FrshFlux_TotalOcean        = global_sum_array(monitor%FrshFlux_TotalOcean)/surface_area
      monitor%FrshFlux_TotalIce          = global_sum_array(monitor%FrshFlux_TotalIce)/surface_area
      monitor%FrshFlux_VolumeIce         = global_sum_array(monitor%FrshFlux_VolumeIce)/surface_area
      monitor%FrshFlux_VolumeTotal       = global_sum_array(monitor%FrshFlux_VolumeTotal)/surface_area
      monitor%HeatFlux_Relax             = global_sum_array(monitor%HeatFlux_Relax)/surface_area
      monitor%FrshFlux_Relax             = global_sum_array(monitor%FrshFlux_Relax)/surface_area
      monitor%SaltFlux_Relax             = global_sum_array(monitor%SaltFlux_Relax)/surface_area
      monitor%TempFlux_Relax             = global_sum_array(monitor%TempFlux_Relax)/surface_area
    ENDIF

    CALL enable_sync_checks()

    ! fluxes through given paths
    IF (my_process_is_stdio() .AND. idbg_val > 0) &
      & WRITE(0,*) "---------------  fluxes --------------------------------"
    DO i=1,oce_section_count
      sflux = section_flux(oce_sections(i), ocean_state%p_prog(nnew(1))%vn)
      !
      ! #slo# disabled since subset%block is not allocated (#3759, HPC_sun_debug)
      ! #ifdef NOMPI
      !     IF (my_process_is_stdio()) &
      !       & write(0,*) oce_sections(i)%subset%name, ":", sflux, 'at edges:',oce_sections(i)%subset%block
      ! #else
      IF (my_process_is_stdio() .AND. idbg_val > 0) &
        & WRITE(0,*) oce_sections(i)%subset%name, ":", sflux
      ! #endif

      SELECT CASE (i)
      CASE (1)
        monitor%gibraltar              = sflux*OceanReferenceDensity
      CASE (2)
        monitor%denmark_strait         = sflux*OceanReferenceDensity
      CASE (3)
        monitor%drake_passage          = sflux*OceanReferenceDensity
      CASE (4)
        monitor%indonesian_throughflow = sflux*OceanReferenceDensity
      CASE (5)
        monitor%scotland_iceland       = sflux*OceanReferenceDensity
      CASE (6)
        monitor%mozambique             = sflux*OceanReferenceDensity
      CASE (7)
        monitor%framStrait             = sflux*OceanReferenceDensity
      CASE (8)
        monitor%beringStrait           = sflux*OceanReferenceDensity
      CASE (9)
        monitor%barentsOpening         = sflux*OceanReferenceDensity
      CASE (10)
        monitor%agulhas                = sflux*OceanReferenceDensity
      CASE (11)
        monitor%agulhas_long           = sflux*OceanReferenceDensity
      CASE (12)
        monitor%agulhas_longer         = sflux*OceanReferenceDensity
      CASE (13)
        monitor%florida_strait         = sflux*OceanReferenceDensity
      END SELECT
    ENDDO

    IF (i_sea_ice >= 1) THEN
      IF (my_process_is_stdio() .AND. idbg_val > 0) &
        & WRITE(0,*) "--------------- ice fluxes -----------------------------"
          ! ice transport through given paths
      sflux = section_ice_flux(oce_sections(7), ice%hi*ice%conc, ice%vn_e) ! TODO: Replace hi/conc to himean

      IF (my_process_is_stdio() .AND. idbg_val > 0) &
        & WRITE(0,*) oce_sections(7)%subset%name, ":", sflux

      monitor%ice_framStrait             = sflux
    ENDIF

    IF (my_process_is_stdio() .AND. idbg_val > 0) &
      & WRITE(0,*) "---------------  end fluxes ----------------------------"

  END SUBROUTINE calc_slow_oce_diagnostics

!<Optimize:inUse>
  SUBROUTINE calc_fast_oce_diagnostics(patch_2d, patch_3d, dolic, prism_thickness, depths, &
          &  p_diag, sea_surface_height, tracers, p_atm_f, hamocc, ice, lhamocc)
    TYPE(t_patch ),TARGET :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(inout)     :: patch_3d
    INTEGER,  POINTER                           :: dolic(:,:)
    REAL(wp), POINTER                           :: prism_thickness(:,:,:)
    REAL(wp), INTENT(in)                        :: depths(:)
    TYPE(t_hydro_ocean_diag), TARGET            :: p_diag
    REAL(wp), POINTER                           :: sea_surface_height(:,:)
    REAL(wp), POINTER                           :: tracers(:,:,:,:)
    TYPE(t_atmos_fluxes ),    INTENT(IN)        :: p_atm_f
    TYPE(t_hamocc_state), TARGET, INTENT(inout) :: hamocc
    TYPE(t_sea_ice),          INTENT(inout)     :: ice
    LOGICAL, INTENT(IN)                         :: lhamocc

    !Local variables
    INTEGER :: start_cell_index, end_cell_index!,i_startblk_c, i_endblk_c,
    INTEGER :: jk,jc,blockNo!,je
    REAL(wp):: ssh_global_mean,total_runoff_flux,total_heat_flux, &
      &        total_fresh_water_flux,total_evaporation_flux, &
      &        ice_volume_nh, ice_volume_sh, ice_extent_nh, ice_extent_sh, &
      &        global_mean_potEnergy, global_mean_kinEnergy, global_mean_totalEnergy, &
      &        global_mean_potEnstrophy

    TYPE(t_subset_range), POINTER :: owned_cells
    TYPE(t_ocean_monitor),  POINTER :: monitor
    !-----------------------------------------------------------------------
    owned_cells    => patch_2d%cells%owned
    monitor        => p_diag%monitor

    !cell loop to calculate cell based monitored fields volume, kinetic energy and tracer content
    SELECT CASE (iswm_oce)
    CASE (1) ! shallow water mode

    CASE default !3D model
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index) SCHEDULE(dynamic)
      DO blockNo = owned_cells%start_block, owned_cells%end_block
        CALL get_index_range(owned_cells, blockNo, start_cell_index, end_cell_index)
        !We are dealing with the surface layer first
        DO jc =  start_cell_index, end_cell_index

         !p_diag%mld(jc,blockNo) = calc_mixed_layer_depth(p_diag%zgrad_rho(jc,:,blockNo),&
         !  & 0.125_wp, &
         !  & dolic(jc,blockNo), &
         !  & prism_thickness(jc,:,blockNo), &
         !  & depths(1))
         !
         ! compute mixed layer depth
          ! save the maximal mixed layer depth  - RESET TO 0 ECHT OUTPUT TIMESTEP
          p_diag%mld(jc,blockNo) = MAX(calc_mixed_layer_depth(p_diag%zgrad_rho(jc,:,blockNo),&
            & 0.125_wp, &
            & dolic(jc,blockNo), &
            & prism_thickness(jc,:,blockNo), &
            & depths(1)),p_diag%mld(jc,blockNo))
          p_diag%condep(jc,blockNo) = REAL(calc_condep(p_diag%zgrad_rho(jc,:,blockNo), dolic(jc,blockNo)),KIND=wp)

        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

      ! {{{ compute global mean values of:
      ! sea surface height
      ssh_global_mean = 0.0_wp
      IF (isRegistered('ssh_global')) THEN
        CALL levels_horizontal_mean( sea_surface_height, &
            & patch_2d%cells%area(:,:), &
            & owned_cells, &
            & ssh_global_mean)
      END IF
      monitor%ssh_global = ssh_global_mean

      ! total heat flux
      total_heat_flux = 0.0_wp
      IF (isRegistered('HeatFlux_Total_Global')) THEN
      call levels_horizontal_mean( p_atm_f%HeatFlux_Total, &
          & patch_2d%cells%area(:,:), &
          & owned_cells, &
          & total_heat_flux)
      END IF
      monitor%HeatFlux_Total = total_heat_flux

      ! total fresh water flux
      total_fresh_water_flux = 0.0_wp
      IF (isRegistered('FrshFlux_Precipitation_Global')) THEN
      call levels_horizontal_mean( p_atm_f%FrshFlux_Precipitation, &
          & patch_2d%cells%area(:,:), &
          & owned_cells, &
          & total_fresh_water_flux)
      END IF
      monitor%FrshFlux_Precipitation = total_fresh_water_flux

      ! total evaporation
      total_evaporation_flux = 0.0_wp
      IF (isRegistered('FrshFlux_Evaporation_Global')) THEN
      call levels_horizontal_mean( p_atm_f%FrshFlux_Evaporation, &
          & patch_2d%cells%area(:,:), &
          & owned_cells, &
          & total_evaporation_flux)
      END IF
      monitor%FrshFlux_Evaporation = total_evaporation_flux

      ! total runoff
      total_runoff_flux = 0.0_wp
      IF (isRegistered('FrshFlux_Runoff_Global')) THEN
      call levels_horizontal_mean( p_atm_f%FrshFlux_Runoff, &
          & patch_2d%cells%area(:,:), &
          & owned_cells, &
          & total_runoff_flux)
      END IF
      monitor%FrshFlux_Runoff = total_runoff_flux

      ! ice volume and extend
      ice_volume_nh = 0.0_wp
      IF (isRegistered('ice_volume_nh')) THEN
      ice_volume_nh = subset_sum( ice%vol(:,1,:)*p_diag%northernHemisphere(:,:), &
          & owned_cells )
      END IF
      monitor%ice_volume_nh = ice_volume_nh/1.0e9_wp !scaling to km^3

      ice_volume_sh = 0.0_wp
      IF (isRegistered('ice_volume_sh')) THEN
      ice_volume_sh = subset_sum( ice%vol(:,1,:)*p_diag%southernHemisphere(:,:), &
          & owned_cells)
      END IF
      monitor%ice_volume_sh = ice_volume_sh/1.0e9_wp !scaling to km^3

      ice_extent_nh = 0.0_wp
      IF (isRegistered('ice_extent_nh')) THEN
      ice_extent_nh = subset_sum( ice%concsum*p_diag%northernHemisphere*patch_2d%cells%area, &
          & owned_cells)
      END IF
      monitor%ice_extent_nh = ice_extent_nh/1.0e6_wp !scaling to km^2

      ice_extent_sh = 0.0_wp
      IF (isRegistered('ice_extent_sh')) THEN
      ice_extent_sh = subset_sum( ice%concsum*p_diag%southernHemisphere*patch_2d%cells%area, &
          & owned_cells)
      END IF
      monitor%ice_extent_sh = ice_extent_sh/1.0e6_wp !scaling to km^2

      ! energy/enstrophy
      global_mean_potEnergy = 0.0_wp
      IF (isRegistered('pot_energy_global')) THEN
        global_mean_potEnergy = potential_energy(& 
            & p_diag%w, &
!TODO       & p_prog(nold(1))%h,&
            & sea_surface_height , & ! this is h_new, the old implementation used h_old
            & p_diag%rho, &
            & patch_3D%p_patch_1d(1)%del_zlev_i, &
            & patch_3D%p_patch_1d(1)%prism_volume, &
            & owned_cells)
      END IF
      monitor%pot_energy = global_mean_potEnergy

      global_mean_kinEnergy = 0.0_wp
      IF (isRegistered('kin_energy_global')) THEN
         global_mean_kinEnergy = total_mean( p_diag%kin, &
          & patch_3d%p_patch_1d(1)%prism_volume, &
          & owned_cells )
      END IF
      monitor%kin_energy = global_mean_kinEnergy

      global_mean_totalEnergy = 0.0_wp
      IF (isRegistered('total_energy_global')) THEN
        ! use precomputed variables
        IF (isRegistered('kin_energy_global') .AND. isRegistered('pot_energy_global')) THEN
          global_mean_totalEnergy = global_mean_kinEnergy + global_mean_potEnergy
        END IF
      END IF
      monitor%total_energy = global_mean_totalEnergy

      global_mean_potEnstrophy = 0.0_wp
      IF (isRegistered('potential_enstrophy_global')) THEN
      END IF
      monitor%potential_enstrophy = global_mean_potEnstrophy
      !}}}

      ! calc moc each timestep from non-accumulated vertical veloc
      if ( isRegistered('global_moc') .or. isRegistered('atlant_moc') .or. isRegistered('pacind_moc') ) then
        CALL timer_start(timer_calc_moc)
        CALL calc_moc(patch_2d, patch_3d, &
           & p_diag%w, &
           & p_diag%global_moc, &
           & p_diag%atlantic_moc, &
           & p_diag%pacific_moc)
        CALL timer_stop(timer_calc_moc)
      endif

      CALL dbg_print('Diag: mld',p_diag%mld,str_module,4,in_subset=owned_cells)
      
      ! hamocc global diagnostics
      IF (lhamocc) CALL get_monitoring( hamocc, sea_surface_height , tracers, patch_3d)

    END SELECT
  END SUBROUTINE calc_fast_oce_diagnostics
  !-------------------------------------------------------------------------

  !TODO potential_energy(&
  !TODO     & p_diag%w,p_prog(nold(1))%h,&
  !TODO     & p_diag%rho, &
  !TODO     & patch_3D%p_patch_1d(1))%del_zlev_i, &
  !TODO     & patch_3D%p_patch_1d(1)%prism_volume, &
  !TODO     & owned_cells)
  !-------------------------------------------------------------------------
  FUNCTION potential_energy(w,h,rho,del_zlev_i,weights,in_subset)
    REAL(wp), INTENT(IN) :: w(:,:,:)
    REAL(wp), INTENT(IN) :: h(:,:)
    REAL(wp), INTENT(IN) :: rho(:,:,:)
    REAL(wp), INTENT(IN) :: del_zlev_i(:)
    REAL(wp), INTENT(IN) :: weights(:,:,:)
    TYPE(t_subset_range), INTENT(IN) :: in_subset

    REAL(wp) :: potential_energy

#define VerticalDim_Position 2

    REAL(wp), ALLOCATABLE :: sum_value(:,:), sum_weight(:,:), total_weight(:), total_sum(:)
    INTEGER :: block, level, start_index, end_index, idx, start_vertical, end_vertical
    INTEGER :: allocated_levels, no_of_threads, myThreadNo
    REAL(wp) :: z_w, totalSum, totalWeight

    CHARACTER(LEN=*), PARAMETER :: method_name=module_name//':potential_energy'


    IF (in_subset%no_of_holes > 0) CALL warning(method_name, "there are holes in the subset")

    no_of_threads = 1
    myThreadNo = 0
#ifdef _OPENMP
    no_of_threads = omp_get_max_threads()
#endif

    allocated_levels = SIZE(w,VerticalDim_Position)
    ALLOCATE( sum_value(allocated_levels, 0:no_of_threads-1), &
      & sum_weight(allocated_levels, 0:no_of_threads-1), &
      & total_weight(allocated_levels), total_sum(allocated_levels) )

    start_vertical = 1
    end_vertical = SIZE(w, VerticalDim_Position)

    IF (start_vertical > end_vertical) &
      & CALL finish(method_name, "start_vertical > end_vertical")
    IF ( allocated_levels < end_vertical) &
      & CALL finish(method_name, "allocated_levels < end_vertical")

!ICON_OMP_PARALLEL PRIVATE(myThreadNo)
!$  myThreadNo = omp_get_thread_num()
!ICON_OMP_SINGLE
!$  no_of_threads = OMP_GET_NUM_THREADS()
!ICON_OMP_END_SINGLE NOWAIT
    sum_value(:,  myThreadNo) = 0.0_wp
    sum_weight(:,  myThreadNo) = 0.0_wp
    IF (ASSOCIATED(in_subset%vertical_levels)) THEN
!ICON_OMP_DO PRIVATE(block, start_index, end_index, idx)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        DO idx = start_index, end_index
          DO level = start_vertical, MIN(end_vertical, in_subset%vertical_levels(idx,block)) - 1
            z_w = MERGE( &
              & (w(idx,level,block)*h(idx,block) &
              &  + w(idx,level+1,block)*0.5_wp*del_zlev_i(level)) &
              & /(0.5_wp*del_zlev_i(level)+h(idx,block)) &
              & , &
              & (w(idx,level,block)*del_zlev_i(level) &
              &  + w(idx,level+1,block)*del_zlev_i(level+1)) &
              & /(del_zlev_i(level)+del_zlev_i(level+1)) &
              & , &
              & 1 .EQ. level)

            sum_value(level, myThreadNo)  = sum_value(level, myThreadNo) + &
              & grav*z_w*rho(idx, level, block) * weights(idx, level, block)

            sum_weight(level, myThreadNo)  = sum_weight(level, myThreadNo) + weights(idx, level, block)

          ENDDO
        ENDDO
      ENDDO
!ICON_OMP_END_DO

    ELSE ! no in_subset%vertical_levels

!ICON_OMP_DO PRIVATE(block, start_index, end_index)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        DO idx = start_index, end_index
          ! since we have the same numbder of vertical layers, the weight is the same
          ! for all levels. Compute it only for the first level, and then copy it
          DO level = start_vertical, end_vertical - 1
            z_w = MERGE( &
              & (w(idx,level,block)*h(idx,block) &
              &  + w(idx,level+1,block)*0.5_wp*del_zlev_i(level)) &
              & /(0.5_wp*del_zlev_i(level)+h(idx,block)) &
              & , &
              & (w(idx,level,block)*del_zlev_i(level) &
              &  + w(idx,level+1,block)*del_zlev_i(level+1)) &
              & /(del_zlev_i(level)+del_zlev_i(level+1)) &
              & , &
              & 1 .EQ. level)

            sum_value(level, myThreadNo)  = sum_value(level, myThreadNo) + &
              & grav*z_w*rho(idx, level, block) * weights(idx,level, block)
            sum_weight(level, myThreadNo)  = sum_weight(start_vertical, myThreadNo) + weights(idx, level, block)
          ENDDO
        ENDDO
      ENDDO
!ICON_OMP_END_DO

    ENDIF
!ICON_OMP_END_PARALLEL

    ! gather the total level sum of this process in total_sum(level)
    total_sum(:)     = 0.0_wp
    total_weight(:) = 0.0_wp
    DO myThreadNo=0, no_of_threads-1
      DO level = start_vertical, end_vertical - 1
        ! write(0,*) myThreadNo, level, " sum=", sum_value(level, myThreadNo), sum_weight(level, myThreadNo)
        total_sum(level)    = total_sum(level)    + sum_value(level, myThreadNo)
        total_weight(level) = total_weight(level) + sum_weight(level, myThreadNo)
      ENDDO
    ENDDO
    DEALLOCATE(sum_value, sum_weight)

    ! Collect the value and weight sums (at all procs)
    CALL gather_sums(total_sum, total_weight)


    totalSum = 0.0_wp
    totalWeight = 0.0_wp
    DO level = start_vertical, end_vertical
      totalSum    = totalSum    + total_sum(level)
      totalWeight = totalWeight + total_weight(level)
    ENDDO
    DEALLOCATE(total_weight)
    DEALLOCATE(total_sum)

    potential_energy = totalSum / totalWeight
  END FUNCTION potential_energy

  !-------------------------------------------------------------------------
  REAL(wp) FUNCTION section_flux(in_oce_section, velocity_values)
    TYPE(t_oce_section) :: in_oce_section
    REAL(wp), POINTER :: velocity_values(:,:,:)

    INTEGER :: i, k, edge_idx, edge_block
    REAL(wp) :: oriented_length
    REAL(wp), ALLOCATABLE :: flux_weights(:,:)
    TYPE(t_grid_edges), POINTER ::  edges
    TYPE(t_patch_vert),POINTER :: patch_vertical

    CHARACTER(LEN=*), PARAMETER :: method_name='mo_ocean_diagnostics:section_flux'

    edges          => in_oce_section%subset%patch%edges
    patch_vertical => in_oce_section%subset%patch_3d%p_patch_1d(1)

    ! calculate weights
    ! flux_weights can also be preallocated
    ALLOCATE(flux_weights(n_zlev, MAX(in_oce_section%subset%SIZE, 1)))
    flux_weights(:,:) = 0.0_wp
    DO i=1, in_oce_section%subset%SIZE

      edge_idx   = in_oce_section%subset%idx(i)
      edge_block = in_oce_section%subset%BLOCK(i)
      oriented_length = edges%primal_edge_length(edge_idx, edge_block) * &
        & in_oce_section%orientation(i) ! this can also be pre-calculated and stored in in_oce_section%orientation

      !write(0,*) "oriented_length:",  oriented_length

      DO k=1, n_zlev
        flux_weights(k, i) = patch_vertical%prism_thick_e(edge_idx, k, edge_block) * oriented_length ! maybe also use slm
        !write(0,*) i, k, in_oce_section%subset%name, " flux_weights:",  flux_weights(k, i), &
        !  & patch_vertical%prism_thick_e(edge_idx, k, edge_block)
        !write(0,*) i, k, in_oce_section%subset%name, " velocity_value:", velocity_values(edge_idx, k, edge_block)
      ENDDO

    ENDDO


    section_flux = subset_sum(                           &
      & values                 = velocity_values,        &
      & indexed_subset         = in_oce_section%subset,  &
      & subset_indexed_weights = flux_weights)

    DEALLOCATE(flux_weights)

    !write(0,*) get_my_mpi_work_id(), ": section_flux on subset ", in_oce_section%subset%name, ":", &
    !  & section_flux, in_oce_section%subset%size

  END FUNCTION section_flux
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  REAL(wp) FUNCTION section_ice_flux(in_oce_section, ice_hmean, ice_vel)
    TYPE(t_oce_section)  :: in_oce_section
    REAL(wp), INTENT(IN) :: ice_hmean(:,:,:)
    REAL(wp), POINTER    :: ice_vel(:,:)

    TYPE(t_grid_edges), POINTER ::  edges
    INTEGER :: i, k, edge_idx, edge_block, cell_idx(2), cell_block(2)
    REAL(wp) :: vel_vals, sum_value, oriented_length
    REAL(wp), ALLOCATABLE :: fluxes(:,:)
    INTEGER :: communicator

    CHARACTER(LEN=*), PARAMETER :: method_name='mo_ocean_diagnostics:section_ice_flux'

    edges          => in_oce_section%subset%patch%edges

    ! calculate weights
    ! fluxes can also be preallocated
    ALLOCATE(fluxes(kice, MAX(in_oce_section%subset%SIZE, 1)))
    fluxes(:,:) = 0.0_wp
    sum_value  = 0.0_wp

    DO i=1, in_oce_section%subset%SIZE

      edge_idx   = in_oce_section%subset%idx(i)
      edge_block = in_oce_section%subset%BLOCK(i)

      cell_idx   = edges%cell_idx(edge_idx, edge_block,:)
      cell_block = edges%cell_blk(edge_idx, edge_block,:)

      oriented_length = edges%primal_edge_length(edge_idx, edge_block) * &
        & in_oce_section%orientation(i) ! this can also be pre-calculated and stored in in_oce_section%orientation

      vel_vals = ice_vel(edge_idx, edge_block)

      ! compute the first order upwind flux using cell-centered ice_hmean vals and edge-centered vel_vals
      ! Same as in upwind_hflux_ice which is used to calculate ice advection
      DO k=1,kice
          fluxes(k, i) = laxfr_upflux( vel_vals, ice_hmean(cell_idx(1),k,cell_block(1)), &
          &                    ice_hmean(cell_idx(2),k,cell_block(2)) ) * oriented_length
          !write(0,*) i, k, in_oce_section%subset%name, " fluxes:",  fluxes(k, i), &
          !  & patch_vertical%prism_thick_e(edge_idx, k, edge_block)
          !write(0,*) i, k, in_oce_section%subset%name, " velocity_value:", velocity_values(edge_idx, k, edge_block)
      ENDDO

    ENDDO

    sum_value = sum(fluxes(:,:))

    DEALLOCATE(fluxes)

    ! the global min, max is avaliable only to stdio process
    IF (my_process_is_mpi_parallel()) THEN

      communicator = in_oce_section%subset%patch%work_communicator
      ! these are avaliable to all processes
      section_ice_flux = p_sum( sum_value,  comm=communicator)

    ELSE

      section_ice_flux = sum_value

    ENDIF

!    write(0,*) get_my_mpi_work_id(), ": section_ice_flux on subset ", in_oce_section%subset%name, ":", &
!      & section_ice_flux, in_oce_section%subset%size

  END FUNCTION section_ice_flux
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !
  !!  Calculation of meridional overturning circulation (MOC)
  !
  !   Calculation of meridional overturning circulation for different basins
  !   (Atlantic, Pacific, Indian, global)
  !>
  !!
  !! @par Revision History
  !! Developed  by  Stephan Lorenz, MPI-M (2012).
  !!  based on code from MPIOM
  !
  ! TODO: implement variable output dimension (1 deg resolution) and smoothing extent
  ! TODO: calculate the 1 deg resolution meridional distance
  !!
  SUBROUTINE calc_moc_acc (patch_2d, patch_3D, w, this_datetime)

    TYPE(t_patch), TARGET, INTENT(in)  :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3D
    REAL(wp), INTENT(in)               :: w(:,:,:)   ! vertical velocity at cell centers
    ! dims: (nproma,nlev+1,alloc_cell_blocks)
    TYPE(datetime), POINTER            :: this_datetime
    !
    ! local variables
    ! INTEGER :: i
    INTEGER, PARAMETER ::  jbrei=3   !  latitudinal smoothing area is 2*jbrei-1 rows of 1 deg
    INTEGER :: blockNo, jc, jk, start_index, end_index !, il_e, ib_e
    INTEGER :: lbrei, lbr, idate, itime
    INTEGER :: mpi_comm
    INTEGER(i8) :: i1,i2,i3,i4

    REAL(wp) :: z_lat, z_lat_deg, z_lat_dim
    REAL(wp) :: global_moc(180,n_zlev), atlant_moc(180,n_zlev), pacind_moc(180,n_zlev)
    REAL(dp) :: local_moc(180), res_moc(180)

    TYPE(t_subset_range), POINTER :: dom_cells

    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: routine = ('mo_ocean_diagnostics:calc_moc')

    !-----------------------------------------------------------------------

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    global_moc(:,:) = 0.0_wp
    pacind_moc(:,:) = 0.0_wp
    atlant_moc(:,:) = 0.0_wp

    ! set barrier:
    ! CALL MPI_BARRIER(0)

    ! with all cells no sync is necessary
    !owned_cells => patch_2d%cells%owned
    dom_cells   => patch_2d%cells%in_domain

    !write(81,*) 'MOC: datetime:',datetime

    DO jk = 1, n_zlev   !  not yet on intermediate levels
      DO blockNo = dom_cells%start_block, dom_cells%end_block
        CALL get_index_range(dom_cells, blockNo, start_index, end_index)
        DO jc = start_index, end_index

          !  could be replaced by vertical loop to bottom
          IF ( patch_3D%lsm_c(jc,jk,blockNo) <= sea_boundary ) THEN

            ! lbrei: corresponding latitude row of 1 deg extension
            !       1 south pole
            !     180 north pole
            z_lat = patch_2d%cells%center(jc,blockNo)%lat
            z_lat_deg = z_lat*rad2deg
            lbrei = NINT(90.0_wp + z_lat_deg)
            lbrei = MAX(lbrei,1)
            lbrei = MIN(lbrei,180)

            ! get neighbor edge for scaling
            !   il_e = patch_2d%cells%edge_idx(jc,blockNo,1)
            !   ib_e = patch_2d%cells%edge_blk(jc,blockNo,1)

            ! z_lat_dim: scale to 1 deg resolution
            ! z_lat_dim: latitudinal extent of triangle divided by latitudinal smoothing extent
            !   z_lat_dim = patch_2d%edges%primal_edge_length(il_e,ib_e) / &
            !     & (REAL(2*jbrei, wp) * 111111._wp*1.3_wp)
            z_lat_dim = 1.0_wp

            ! distribute MOC over (2*jbrei)+1 latitude rows
            !  - no weighting with latitudes done
            !  - lbrei: index of 180 X 1 deg meridional resolution
            DO lbr = -jbrei, jbrei
              lbrei = NINT(90.0_wp + z_lat_deg + REAL(lbr, wp) * z_lat_dim)
              lbrei = MAX(lbrei,1)
              lbrei = MIN(lbrei,180)

              global_moc(lbrei,jk) = global_moc(lbrei,jk) - &
              !  multiply with wet (or loop to bottom)
                & patch_2d%cells%area(jc,blockNo) * OceanReferenceDensity * w(jc,jk,blockNo) * &
                & patch_3D%wet_c(jc,jk,blockNo) / &
                & REAL(2*jbrei + 1, wp)

              IF (patch_3D%basin_c(jc,blockNo) == 1) THEN         !  1: Atlantic; 0: Land

                atlant_moc(lbrei,jk) = atlant_moc(lbrei,jk) - &
                  & patch_2d%cells%area(jc,blockNo) * OceanReferenceDensity * w(jc,jk,blockNo) * &
                  & patch_3D%wet_c(jc,jk,blockNo) / &
                  & REAL(2*jbrei + 1, wp)
              ELSE IF (patch_3D%basin_c(jc,blockNo) >= 2) THEN   !  2: Indian; 4: Pacific
                pacind_moc(lbrei,jk) = pacind_moc(lbrei,jk) - &
                  & patch_2d%cells%area(jc,blockNo) * OceanReferenceDensity * w(jc,jk,blockNo) * &
                  & patch_3D%wet_c(jc,jk,blockNo) / &
                  & REAL(2*jbrei + 1, wp)
              END IF

            END DO

          END IF
        END DO
      END DO

      ! test parallelization:
      ! function field_sum_all using mpi_allreduce and working precisions wp does not exist
      ! res_moc(:) = p_field_sum_all_wp(global_moc(:,jk))
      ! res_moc(:) = p_field_sum_all_wp(atlant_moc(:,jk))
      ! res_moc(:) = p_field_sum_all_wp(pacind_moc(:,jk))

      ! function field_sum using mpi_reduce, then broadcast
      local_moc(:)     = REAL(global_moc(:,jk),dp)
      res_moc(:)       = p_field_sum(local_moc, mpi_comm)
      CALL p_bcast(res_moc(:), p_io, mpi_comm)
      global_moc(:,jk) = REAL(res_moc(:),wp)

      local_moc(:)     = REAL(atlant_moc(:,jk),dp)
      res_moc(:)       = p_field_sum(local_moc, mpi_comm)
      CALL p_bcast(res_moc(:), p_io, mpi_comm)
      atlant_moc(:,jk) = REAL(res_moc(:),wp)

      local_moc(:)     = REAL(pacind_moc(:,jk),dp)
      res_moc(:)       = p_field_sum(local_moc, mpi_comm)
      CALL p_bcast(res_moc(:), p_io, mpi_comm)
      pacind_moc(:,jk) = REAL(res_moc(:),wp)

    END DO  ! n_zlev-loop

    IF (my_process_is_stdio()) THEN
      DO lbr=179,1,-1   ! fixed to 1 deg meridional resolution

        global_moc(lbr,:)=global_moc(lbr+1,:)+global_moc(lbr,:)
        atlant_moc(lbr,:)=atlant_moc(lbr+1,:)+atlant_moc(lbr,:)
        pacind_moc(lbr,:)=pacind_moc(lbr+1,:)+pacind_moc(lbr,:)

      END DO

      ! write out MOC in extra format, file opened in mo_hydro_ocean_run  - integer*8
      !  - correct date in extra format - i.e YYYYMMDD - no time info
      !idate=datetime%month*1000000+datetime%day*10000+datetime%hour*100+datetime%minute
      idate = this_datetime%date%year*10000+this_datetime%date%month*100+this_datetime%date%day
      itime = this_datetime%time%hour*100+this_datetime%time%minute
      WRITE(message_text,*) 'Write MOC at year =',this_datetime%date%year,', date =',idate,' time =', itime
      CALL message (TRIM(routine), message_text)

      DO jk = 1,n_zlev
        i1=INT(idate,i8)
        i2 = INT(777,i8)
        i3 = INT(patch_3D%p_patch_1d(1)%zlev_i(jk),i8)
        i4 = INT(180,i8)
        WRITE(moc_unit) i1,i2,i3,i4
        WRITE(moc_unit) (global_moc(lbr,jk),lbr=1,180)
        i2 = INT(778,i8)
        WRITE(moc_unit) i1,i2,i3,i4
        WRITE(moc_unit) (atlant_moc(lbr,jk),lbr=1,180)
        i2 = INT(779,i8)
        WRITE(moc_unit) i1,i2,i3,i4
        WRITE(moc_unit) (pacind_moc(lbr,jk),lbr=1,180)

      END DO
    END IF

  END SUBROUTINE calc_moc_acc
  SUBROUTINE calc_moc_internal (patch_2d, patch_3D, w, global_moc, atlant_moc, pacind_moc)

    TYPE(t_patch), TARGET, INTENT(in)  :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3D
    REAL(wp), INTENT(in)               :: w(:,:,:)   ! vertical velocity at cell centers
    ! dims: (nproma,nlev+1,alloc_cell_blocks)
    REAL(wp), INTENT(OUT) :: global_moc(:,:), atlant_moc(:,:), pacind_moc(:,:) !dims(n_zlev,180)
    !
    ! local variables
    ! INTEGER :: i
    INTEGER, PARAMETER ::  latSmooth = 3   !  latitudinal smoothing area is 2*jbrei-1 rows of 1 deg
    INTEGER :: block, level, start_index, end_index, idx
    INTEGER :: ilat, l, idate, itime
    INTEGER :: mpi_comm

    REAL(wp) :: lat, ilat_inv, deltaMoc
    REAL(dp) :: local_moc(180), res_moc(180)

    TYPE(t_subset_range), POINTER :: cells

    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: routine = ('mo_ocean_diagnostics:calc_moc')

    !-----------------------------------------------------------------------
    mpi_comm = MERGE(p_comm_work_test, p_comm_work, p_test_run)

    global_moc(:,:) = 0.0_wp
    pacind_moc(:,:) = 0.0_wp
    atlant_moc(:,:) = 0.0_wp

    ! with all cells no sync is necessary
    !owned_cells => patch_2d%cells%owned
    cells   => patch_2d%cells%in_domain

!    DO jk = 1, n_zlev   !  not yet on intermediate levels
      DO block = cells%start_block, cells%end_block
        CALL get_index_range(cells, block, start_index, end_index)
        DO idx = start_index, end_index
          lat = patch_2d%cells%center(idx,block)%lat*rad2deg
          DO level = 1, cells%vertical_levels(idx,block)

            deltaMoc = patch_2d%cells%area(idx,block) * OceanReferenceDensity * w(idx,level,block)

            ! lat: corresponding latitude row of 1 deg extension
            !       1 south pole
            !     180 north pole
            ilat = NINT(90.0_wp + lat)
            ilat = MAX(ilat,1)
            ilat = MIN(ilat,180)

            ! distribute MOC over (2*jbrei)+1 latitude rows
            !  - no weighting with latitudes done
            !  - lat: index of 180 X 1 deg meridional resolution
            DO l = -latSmooth, latSmooth
              ilat = NINT(90.0_wp + lat + REAL(l, wp))
              ilat = MAX(ilat,1)
              ilat = MIN(ilat,180)
              ilat_inv = 1.0_wp / REAL(2*latSmooth + 1, wp)

              global_moc(level,ilat) =       global_moc(level,ilat) - deltaMoc*ilat_inv
              atlant_moc(level,ilat) = MERGE(atlant_moc(level,ilat) - deltaMoc*ilat_inv, 0.0_wp, patch_3D%basin_c(idx,block) == 1)
              pacind_moc(level,ilat) = MERGE(pacind_moc(level,ilat) - deltaMoc*ilat_inv, 0.0_wp, patch_3D%basin_c(idx,block) >= 2)
            END DO
          END DO
        END DO
      END DO
!   END DO

      ! test parallelization:
      ! function field_sum_all using mpi_allreduce and working precisions wp does not exist
      ! res_moc(:) = p_field_sum_all_wp(global_moc(:,jk))
      ! res_moc(:) = p_field_sum_all_wp(atlant_moc(:,jk))
      ! res_moc(:) = p_field_sum_all_wp(pacind_moc(:,level))

      global_moc = p_field_sum(global_moc,mpi_comm)
      atlant_moc = p_field_sum(atlant_moc,mpi_comm)
      pacind_moc = p_field_sum(pacind_moc,mpi_comm)

,     ! function field_sum using mpi_reduce, then broadcast
!     local_moc(:)     = REAL(global_moc(level,:),dp)
!     res_moc(:)       = p_field_sum(local_moc,mpi_comm)
!     CALL p_bcast(res_moc(:), p_io, mpi_comm)
!     global_moc(level,:) = REAL(res_moc(:),wp)
!
!     local_moc(:)     = REAL(atlant_moc(level,:),dp)
!     res_moc(:)       = p_field_sum(local_moc, mpi_comm)
!     CALL p_bcast(res_moc(:), p_io, mpi_comm)
!     atlant_moc(level,:) = REAL(res_moc(:),wp)
!
!     local_moc(:)     = REAL(pacind_moc(level,:),dp)
!     res_moc(:)       = p_field_sum(local_moc, mpi_comm)
!     CALL p_bcast(res_moc(:), p_io, mpi_comm)
!     pacind_moc(level,:) = REAL(res_moc(:),wp)


    IF (my_process_is_mpi_workroot()) THEN
      DO l=179,1,-1   ! fixed to 1 deg meridional resolution

        global_moc(:,l)=global_moc(:,l+1)+global_moc(:,l)
        atlant_moc(:,l)=atlant_moc(:,l+1)+atlant_moc(:,l)
        pacind_moc(:,l)=pacind_moc(:,l+1)+pacind_moc(:,l)

      END DO
    END IF

  END SUBROUTINE calc_moc_internal
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !
  !!  Calculation of horizontal stream function
  !
  !>
  !!
  !! @par Revision History
  !! Developed  by  Stephan Lorenz, MPI-M (2012).
  !!  based on code from MPIOM
  !
  ! TODO: implement variable output dimension (1 deg resolution) and smoothing extent
  !!
!<Optimize:inUse>
  SUBROUTINE calc_psi (patch_3D, u, prism_thickness, u_vint, this_datetime)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3D
    REAL(wp), INTENT(in)               :: u(:,:,:)     ! zonal velocity at cell centers
    REAL(wp), INTENT(in)               :: prism_thickness(:,:,:)       ! elevation on cell centers
    ! dims: (nproma,nlev,alloc_cell_blocks)
    REAL(wp), INTENT(inout)            :: u_vint(:,:)  ! barotropic zonal velocity on icon grid
    TYPE(datetime), POINTER            :: this_datetime
    !
    ! local variables
    ! INTEGER :: i

    ! switch for writing stream function (not yet in namelist); 1: icon-grid; 2: regular grid output
    INTEGER, PARAMETER ::  idiag_psi = 1

    INTEGER, PARAMETER ::  nlat = 180                    ! meridional dimension of regular grid
    INTEGER, PARAMETER ::  nlon = 360                    ! zonal dimension of regular grid

    ! smoothing area is 2*jsmth-1 lat/lon areas of 1 deg
    INTEGER, PARAMETER ::  jsmth = 3
    INTEGER :: blockNo, jc, jk, start_index, end_index
    INTEGER :: jlat, jlon, jlt, jln, jltx, jlnx, jsmth2
    INTEGER(i8)        :: idate, iextra(4)


    REAL(wp) :: z_lat_deg, z_lon_deg, z_lat_dist, delta_z, rsmth
    REAL(wp) :: z_uint_reg(nlon,nlat)                     ! vertical integral on regular grid
    REAL(wp) :: psi_reg(nlon,nlat)                        ! horizontal stream function

    TYPE(t_patch), POINTER  :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells, dom_cells

    !CHARACTER(len=max_char_length), PARAMETER :: routine = ('mo_ocean_diagnostics:calc_psi')

    !-----------------------------------------------------------------------

    psi_reg(:,:)    = 0.0_wp
    z_uint_reg(:,:) = 0.0_wp

    jsmth2          = 2*jsmth + 1
    rsmth           = REAL(jsmth2*jsmth2, wp)

    ! with all cells no sync is necessary
    patch_2d  => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    dom_cells => patch_2d%cells%in_domain

    ! (1) barotropic system:
    !     vertical integration of zonal velocity times vertical layer thickness [m/s*m]
!ICON_OMP_PARALLEL_DO PRIVATE(jc, jk, start_index, end_index) SCHEDULE(dynamic)
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      u_vint(:,blockNo)     = 0.0_wp
      DO jc = start_index, end_index

        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)

          u_vint(jc,blockNo) = u_vint(jc,blockNo) - u(jc,jk,blockNo) * prism_thickness(jc,jk,blockNo)

        END DO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('calc_psi:    u_vint    ', u_vint        , str_module, idt_src, in_subset=patch_2d%cells%owned)
    !---------------------------------------------------------------------

    IF (idiag_psi == 1) RETURN

    ! (2) distribute integrated zonal velocity (u*dz) on 1x1 deg grid
    !     this code is not mature yet

    ! in domain: count all cells only once
    DO blockNo = dom_cells%start_block, dom_cells%end_block
      CALL get_index_range(dom_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        z_lat_deg = patch_2d%cells%center(jc,blockNo)%lat * rad2deg
        z_lon_deg = patch_2d%cells%center(jc,blockNo)%lon * rad2deg

        !  ! 0 <= lon <= 360 deg
        !  z_lon_deg = z_lon_deg + 180.0_wp

        ! jlat/jlon: corresponding latitude/longitude coordinates of 1 deg extension
        ! jlat: 1 = south of 89.0S; 89 = 1S-Eq.; 90 = Eq-1N;  180 = north of 89N
        ! jlon: 1 = 180W-179W; 180 = 1-0 deg west; 360 = 179E-180E

        jlat = NINT(91.0_wp + z_lat_deg)
        jlon = NINT(z_lon_deg + 180.5_wp)

        ! distribute stream function over rsmth=(2*jsmth+1)**2 lat/lon regular grid points
        !  - no weighting with latitudes done
        !  - no correction with regular lsm done
        DO jltx = jlat-jsmth, jlat+jsmth

          jlt = jltx
          IF (jlt <    1) jlt =      1-jlt  ! apply equatorwards
          IF (jlt > nlat) jlt = 2*nlat-jlt  ! apply equatorwards
          DO jlnx = jlon-jsmth, jlon+jsmth

            jln = jlnx
            IF (jln <    1) jln = jln+nlon  ! circular boundary
            IF (jln > nlon) jln = jln-nlon  ! circular boundary

            z_uint_reg(jln,jlt) = z_uint_reg(jln,jlt) + u_vint(jc,blockNo) / rsmth


          END DO
        END DO

      END DO
    END DO

    ! (3) calculate meridional integral on regular grid starting from south pole:

    DO jlt = nlat-1, 1, -1
      z_uint_reg(:,jlt) = z_uint_reg(:,jlt) + z_uint_reg(:,jlt+1)
    END DO

    ! (4) calculate stream function: scale with length of 1 deg*rho [m/s*m*m*kg/m3=kg/s]

    ! meridional distance of 1 deg
    ! ATTENTION - fixed 1 deg resolution should be related to icon-resolution
    z_lat_dist = 111111.0_wp  ! * 1.3_wp ??

    psi_reg(:,:) = z_uint_reg(:,:) * z_lat_dist * OceanReferenceDensity

    ! stream function on icon grid without calculation of meridional integral
    !  - tbd after interpolation to regular grid externally
    !  psi    (:,:) = u_vint    (:,:)              * OceanReferenceDensity


    ! write out in extra format - integer*8
    idate = INT(this_datetime%date%month*1000000+this_datetime%date%day*10000 &
         &     +this_datetime%time%hour*100+this_datetime%time%minute,i8)
    WRITE(0,*) 'write global PSI at iyear, idate:',this_datetime%date%year, idate

    iextra(1) = INT(idate,i8)
    iextra(2) = INT(780,i8)
    iextra(3) = INT(0,i8)
    iextra(4) = INT(nlon*nlat,i8)

    WRITE(80) (iextra(blockNo),blockNo=1,4)
    WRITE(80) ((psi_reg(jln,jlt),jln=1,nlon),jlt=1,nlat)

    DO jlat=1,nlat
      WRITE(82,*) 'jlat=',jlat
      WRITE(82,'(1p10e12.3)') (psi_reg(jlon,jlat),jlon=1,nlon)
    ENDDO

    !---------DEBUG DIAGNOSTICS-------------------------------------------
 !  idt_src=3  ! output print level (1-5, fix)
 !  CALL dbg_print('calc_psi:    psi_reg   ', psi_reg       , str_module, idt_src, in_subset=patch_2d%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE calc_psi
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !
  !!  Calculation of horizontal stream function using normal velocity on edges
  !
  !>
  !>  Calculation of horizontal stream function using normal velocity on edges
  !!
  !! @par Revision History
  !! Developed  by  Stephan Lorenz, MPI-M (2014).
  !!
!<Optimize:inUse>
  SUBROUTINE calc_psi_vn (patch_3D, vn, prism_thickness_e, op_coeff, u_vint, v_vint, this_datetime)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3D
    REAL(wp), INTENT(in)               :: vn(:,:,:)                 ! normal velocity at cell edges
    REAL(wp), INTENT(in)               :: prism_thickness_e(:,:,:)  ! elevation on cell edges
    ! dims: (nproma,nlev,alloc_edge_blocks)
    REAL(wp), INTENT(inout)            :: u_vint(:,:)               ! barotropic zonal velocity on cell centers
    REAL(wp), INTENT(inout)            :: v_vint(:,:)               ! barotropic meridional velocity on cell centers
    TYPE(t_operator_coeff),INTENT(in)  :: op_coeff
    TYPE(datetime), POINTER            :: this_datetime
    !
    INTEGER  :: blockNo, jc, je, jk, start_index, end_index
    ! vertical integral vn on edges and in cartesian coordinates - no 2-dim mapping is available
    REAL(wp) :: vn_vint(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: u_2d(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)   !  scratch arrays for test
    REAL(wp) :: v_2d(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)   !  scratch arrays for test
    TYPE(t_cartesian_coordinates) :: vint_cc(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    TYPE(t_patch), POINTER  :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges, all_cells

    !CHARACTER(len=max_char_length), PARAMETER :: routine = ('mo_ocean_diagnostics:calc_psi_vn')

    !-----------------------------------------------------------------------

    patch_2d  => patch_3d%p_patch_2d(1)
    all_edges => patch_2d%edges%ALL
    all_cells => patch_2d%cells%ALL

    vn_vint  (:,:,:)    = 0.0_wp
    vint_cc(:,:,:)%x(1) = 0.0_wp
    vint_cc(:,:,:)%x(2) = 0.0_wp
    vint_cc(:,:,:)%x(3) = 0.0_wp
    u_2d       (:,:)    = 0.0_wp
    v_2d       (:,:)    = 0.0_wp

    ! (1) barotropic system:
    !     vertical integration of normal velocity times vertical layer thickness [m/s*m]
!ICON_OMP_PARALLEL_DO PRIVATE(je, jk, start_index, end_index) SCHEDULE(dynamic)
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_index, end_index)
      vn_vint(:,1,blockNo) = 0.0_wp
      DO je = start_index, end_index

        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

          vn_vint(je,1,blockNo) = vn_vint(je,1,blockNo) - vn(je,jk,blockNo) * prism_thickness_e(je,jk,blockNo)

        END DO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

    ! (2) remapping normal velocity to zonal and meridional velocity at cell centers
    CALL sync_patch_array(sync_e, patch_2d, vn_vint)

    CALL map_edges2cell_3d(patch_3D, vn_vint, op_coeff, vint_cc)

    CALL sync_patch_array(sync_c, patch_2d, vint_cc(:,:,:)%x(1))
    CALL sync_patch_array(sync_c, patch_2d, vint_cc(:,:,:)%x(2))
    CALL sync_patch_array(sync_c, patch_2d, vint_cc(:,:,:)%x(3))

    ! calculate zonal and meridional velocity:
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        CALL cvec2gvec(vint_cc(jc,1,blockNo)%x(1), vint_cc(jc,1,blockNo)%x(2), vint_cc(jc,1,blockNo)%x(3), &
          &            patch_2d%cells%center(jc,blockNo)%lon, patch_2d%cells%center(jc,blockNo)%lat,  &
          &            u_2d(jc,blockNo), v_2d(jc,blockNo))
!         &            u_vint(jc,blockNo), v_vint(jc,blockNo))
      END DO
    END DO
    !CALL sync_patch_array(sync_c, patch_2d, u_vint)
    !CALL sync_patch_array(sync_c, patch_2d, v_vint)
    CALL sync_patch_array(sync_c, patch_2d, u_2d)
    CALL sync_patch_array(sync_c, patch_2d, v_2d)

    ! hack for test: calc_psy for u_vint, calc_psi_vn for v_vint - accumulated and written out
    v_vint(:,:) = u_2d(:,:)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('calc_psi_vn: u_2d        ', u_2d        , str_module, idt_src, in_subset=patch_2d%cells%owned)
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('calc_psi_vn: v_2d        ', v_2d        , str_module, idt_src, in_subset=patch_2d%cells%owned)
    CALL dbg_print('calc_psi_vn: vint_cc%x(1)' ,vint_cc%x(1), str_module, idt_src, in_subset=patch_2d%cells%owned)
    !---------------------------------------------------------------------


  END SUBROUTINE calc_psi_vn
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!taken from MPIOM
!<Optimize:inUse>
  FUNCTION calc_condep(vertical_density_gradient,max_lev) result(condep)
    REAL(wp),INTENT(in)  :: vertical_density_gradient(n_zlev)
    INTEGER, INTENT(in)  :: max_lev
    INTEGER :: condep

    INTEGER :: jk
    INTEGER :: maxcondep     !< maximum convective penetration level
    REAL(wp) :: masked_vertical_density_gradient(n_zlev)

    condep = 1

    ! remove dbl_eps, which  is added in the vertical gradient computation
    masked_vertical_density_gradient = MAX(vertical_density_gradient - dbl_eps,0.0_wp)

    !! diagnose maximum convection level
    !! condep = maximum model level penetrated by vertically continous
    !! convection from the surface downward
    !! calculated over integration period ; it should be written out
    !! as snapshot at the end of the run
   maxcondep=1
    DO jk=2,max_lev
      IF (masked_vertical_density_gradient(jk) .ne. 0.0_wp) THEN
        maxcondep = jk
        EXIT
      ENDIF
    ENDDO

    condep = max0(maxcondep,condep)
  END FUNCTION calc_condep
!<Optimize:inUse>
  FUNCTION calc_mixed_layer_depth(vertical_density_gradient,critical_value,max_lev,thickness, depth_of_first_layer) &
    & result(mixed_layer_depth)
    REAL(wp), TARGET :: vertical_density_gradient(n_zlev)
    REAL(wp), INTENT(in)  :: critical_value
    INTEGER,  INTENT(in)  :: max_lev
    REAL(wp), INTENT(in)  :: thickness(n_zlev)
    REAL(wp), INTENT(in)  :: depth_of_first_layer

    REAL(wp) :: sigh        ,zzz
    REAL(wp) :: mixed_layer_depth
    REAL(wp) :: masked_vertical_density_gradient(n_zlev)
    INTEGER :: jk

    sigh              = critical_value
    mixed_layer_depth = depth_of_first_layer
    masked_vertical_density_gradient = MAX(vertical_density_gradient,0.0_wp)

    ! This diagnostic calculates the mixed layer depth.
    ! It uses the incremental density increase between two
    ! levels and substracts it from the initial density criterion (sigcrit)
    ! and adds the level thickness (zzz) to zmld. This is done till
    ! the accumulated density increase between the surface and
    ! layer k is sigcrit or sigh = O, respectively.

    ! stabio(k) = insitu density gradient
    ! sigh = remaining density difference

    DO jk = 2, max_lev
      IF (sigh .GT. 1.e-6_wp) THEN
        zzz               = MIN(sigh/(ABS(masked_vertical_density_gradient(jk))+1.0E-19_wp),thickness(jk))
        sigh              = MAX(0._wp, sigh-zzz*masked_vertical_density_gradient(jk))
        mixed_layer_depth = mixed_layer_depth + zzz
      ELSE
        sigh = 0._wp
      ENDIF
    ENDDO

  END FUNCTION calc_mixed_layer_depth

  FUNCTION calc_total_salt_content(patch_2d, thickness, ice, ocean_state, surface_fluxes, zUnderIce, &
      & computation_type) RESULT(total_salt_content)
    TYPE(t_patch),POINTER                                 :: patch_2d
    REAL(wp),DIMENSION(nproma,patch_2d%alloc_cell_blocks),INTENT(IN) :: zUnderIce
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: thickness
    TYPE (t_sea_ice),       INTENT(IN)                    :: ice
    TYPE(t_hydro_ocean_state)                             :: ocean_state
    TYPE(t_ocean_surface)                                 :: surface_fluxes
    INTEGER,INTENT(IN), OPTIONAL                          :: computation_type
    REAL(wp)                                              :: total_salt_content

    REAL(wp), DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks) :: salt

    salt               = calc_salt_content(patch_2d, thickness, ice, ocean_state, surface_fluxes, zUnderIce, computation_type)
    total_salt_content = global_sum_array(salt)
  END FUNCTION calc_total_salt_content

  FUNCTION calc_salt_content(patch_2d, thickness, ice, ocean_state, surface_fluxes, zUnderIce, &
      & computation_type) &
      & RESULT(salt)
    TYPE(t_patch),POINTER                                 :: patch_2d
    REAL(wp),DIMENSION(nproma,patch_2d%alloc_cell_blocks),INTENT(IN) :: zUnderIce
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: thickness
    TYPE (t_sea_ice),       INTENT(IN)                    :: ice
    TYPE(t_hydro_ocean_state)                             :: ocean_state
    TYPE(t_ocean_surface)                                 :: surface_fluxes
    INTEGER,INTENT(IN), OPTIONAL                          :: computation_type

    ! locals
    REAL(wp), DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks) :: salt

    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: saltInSeaice, saltInLiquidWater
    TYPE(t_subset_range), POINTER                         :: subset
    INTEGER                                               :: block, cell, cellStart,cellEnd, level
    INTEGER                                               :: my_computation_type
    IF(no_tracer<=1)RETURN
    my_computation_type = 0
    salt         = 0.0_wp

    CALL assign_if_present(my_computation_type, computation_type)
    subset => patch_2d%cells%owned
    DO block = subset%start_block, subset%end_block
      CALL get_index_range(subset, block, cellStart, cellEnd)
      DO cell = cellStart, cellEnd
        IF (subset%vertical_levels(cell,block) < 1) CYCLE
        SELECT CASE (my_computation_type)
        CASE (0) ! use zunderIce for volume in tracer change, multiply flux with top layer salinity
        ! surface:
          saltInSeaice(cell,block)    = sice*rhoi &
            &                         * SUM(ice%hi(cell,:,block)*ice%conc(cell,:,block)) &
            &                         * patch_2d%cells%area(cell,block)
        !!DN This is no longer needed since we now update surface salinity
        !directly
        !!DN   ocean_state%p_prog(nold(1))%tracer(cell,1,block,2) = (ocean_state%p_prog(nold(1))%tracer(cell,1,block,2)*zUnderIce(cell,block) &
        !!DN   &                                            -   dtime &
        !!DN   &                                              * surface_fluxes%FrshFlux_TotalSalt(cell,block) &
        !!DN   &                                              * ocean_state%p_prog(nold(1))%tracer(cell,1,block,2)) &
        !!DN   &                                           /ice%zUnderIce(cell,block)
          saltInLiquidWater(cell,block) = ocean_state%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * zUnderIce(cell,block)*OceanReferenceDensity &
            &                    * patch_2d%cells%area(cell,block)
        CASE (1) ! use zunderIce for volume in tracer change, multiply flux with top layer salinity
        ! surface:
          saltInSeaice(cell,block)    = sice*rhoi &
            &                         * SUM(ice%hi(cell,:,block)*ice%conc(cell,:,block)) &
            &                         * patch_2d%cells%area(cell,block)
        !!DN This is no longer needed since we now update surface salinity directly
          ocean_state%p_prog(nold(1))%tracer(cell,1,block,2) = &
            & (ocean_state%p_prog(nold(1))%tracer(cell,1,block,2)*zUnderIce(cell,block) &
            &   -  dtime  * surface_fluxes%FrshFlux_TotalSalt(cell,block) &
            &      * ocean_state%p_prog(nold(1))%tracer(cell,1,block,2)) &
            &  /zUnderIce(cell,block)
          saltInLiquidWater(cell,block) = ocean_state%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * zUnderIce(cell,block)*OceanReferenceDensity &
            &                    * patch_2d%cells%area(cell,block)
        END SELECT

        salt(cell,1,block) = saltInSeaice(cell,block) + saltInLiquidWater(cell,block)
        DO level=2,subset%vertical_levels(cell,block)
          salt(cell,level,block) = ocean_state%p_prog(nold(1))%tracer(cell,level,block,2) &
            &                    * thickness(cell,level,block)*OceanReferenceDensity &
            &                    * patch_2d%cells%area(cell,block)
        END DO

        ! rest of the underwater world
      END DO ! cell
    END DO !block
  END FUNCTION calc_salt_content

  SUBROUTINE reset_ocean_monitor(monitor)
    TYPE(t_ocean_monitor) :: monitor
    monitor%volume(:)                     = 0.0_wp
!   monitor%kin_energy(:)                 = 0.0_wp
!   monitor%pot_energy(:)                 = 0.0_wp
!   monitor%total_energy(:)               = 0.0_wp
    monitor%total_salt(:)                 = 0.0_wp
!   monitor%vorticity(:)                  = 0.0_wp
!   monitor%enstrophy(:)                  = 0.0_wp
!   monitor%potential_enstrophy(:)        = 0.0_wp
    monitor%absolute_vertical_velocity(:) = 0.0_wp
    monitor%HeatFlux_ShortWave(:)         = 0.0_wp
    monitor%HeatFlux_LongWave(:)          = 0.0_wp
    monitor%HeatFlux_Sensible(:)          = 0.0_wp
    monitor%HeatFlux_Latent(:)            = 0.0_wp
    monitor%FrshFlux_SnowFall(:)          = 0.0_wp
    monitor%FrshFlux_TotalSalt(:)         = 0.0_wp
    monitor%FrshFlux_TotalOcean(:)        = 0.0_wp
    monitor%FrshFlux_TotalIce(:)          = 0.0_wp
    monitor%FrshFlux_VolumeIce(:)         = 0.0_wp
    monitor%FrshFlux_VolumeTotal(:)       = 0.0_wp
    monitor%HeatFlux_Relax(:)             = 0.0_wp
    monitor%FrshFlux_Relax(:)             = 0.0_wp
    monitor%TempFlux_Relax(:)             = 0.0_wp
    monitor%SaltFlux_Relax(:)             = 0.0_wp
    monitor%ice_framStrait(:)             = 0.0_wp
    monitor%florida_strait(:)             = 0.0_wp
    monitor%gibraltar(:)                  = 0.0_wp
    monitor%denmark_strait(:)             = 0.0_wp
    monitor%drake_passage(:)              = 0.0_wp
    monitor%indonesian_throughflow(:)     = 0.0_wp
    monitor%scotland_iceland(:)           = 0.0_wp
    monitor%mozambique(:)                 = 0.0_wp
    monitor%framStrait(:)                 = 0.0_wp
    monitor%beringStrait(:)               = 0.0_wp
    monitor%barentsOpening(:)             = 0.0_wp
    monitor%agulhas(:)                    = 0.0_wp
    monitor%agulhas_long(:)               = 0.0_wp
    monitor%agulhas_longer(:)             = 0.0_wp
    monitor%t_mean_na_200m(:)             = 0.0_wp
    monitor%t_mean_na_800m(:)             = 0.0_wp
    monitor%ice_ocean_heat_budget(:)      = 0.0_wp
    monitor%ice_ocean_salinity_budget(:)  = 0.0_wp
    monitor%ice_ocean_volume_budget(:)    = 0.0_wp
  END SUBROUTINE reset_ocean_monitor
END MODULE mo_ocean_diagnostics
