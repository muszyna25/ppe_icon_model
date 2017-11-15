!>
!! Provide an implementation of the sea-ice model.
!!
!! Provide an implementation of the parameters of the surface module (sea ice)
!! used between the atmopshere and the hydrostatic ocean model.
!!
!! @author Peter Korn, MPI
!! @author Dirk Notz, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!----------------------------
#include "omp_definitions.inc"
#ifndef _OPENMP
#include "consistent_fma.inc"
#endif
!----------------------------
MODULE mo_sea_ice
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2007
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime, ltimer
  USE mo_coupling_config,     ONLY: is_coupled_run
  USE mo_dynamics_config,     ONLY: nold, nnew
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D, t_patch_vert
  USE mo_exception,           ONLY: finish, message
  USE mo_impl_constants,      ONLY: success, max_char_length, sea_boundary
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL,      &
    &                               GRID_UNSTRUCTURED_VERT, GRID_VERTEX,    &
    &                               GRID_UNSTRUCTURED_EDGE, GRID_EDGE
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref, ki, ks, Tf, albi, albim, albsm, albs,   &
    &                               fr_fac, mu, alf, alv, albedoW_sim, clw, cpd, zemiss_def, rd, &
    &                               stbo, tmelt, ci, Cd_ia, sice, alb_sno_vis, alb_sno_nir,      &
    &                               alb_ice_vis, alb_ice_nir
  USE mo_math_constants,      ONLY: rad2deg
  USE mo_statistics,          ONLY: add_fields
  USE mo_ocean_nml,           ONLY: no_tracer, limit_seaice, seaice_limit
  USE mo_sea_ice_nml,         ONLY: i_ice_therm, i_ice_dyn, hnull, hmin, hci_layer, &
    &                               i_ice_albedo, leadclose_1, leadclose_2n, use_IceInitialization_fromTemperature, &
    &                               use_constant_tfreez, use_calculated_ocean_stress, use_no_flux_gradients, t_heat_base, &
    &                               init_analytic_conc_param, init_analytic_hi_param, &
    &                               init_analytic_hs_param, init_analytic_temp_under_ice
  USE mo_ocean_types,           ONLY: t_hydro_ocean_state
  USE mo_ocean_state,           ONLY: v_base, ocean_restart_list, ocean_default_list
  USE mo_var_list,            ONLY: add_var
  USE mo_var_metadata,        ONLY: groups
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_cf_convention
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, GRID_UNSTRUCTURED
  USE mo_zaxis_type,          ONLY: ZA_GENERIC_ICE, ZA_SURFACE
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, &
    &                               t_atmos_for_ocean, t_sea_ice_acc, t_sea_ice_budgets
  USE mo_sea_ice_winton,      ONLY: ice_growth_winton, set_ice_temp_winton
  USE mo_sea_ice_zerolayer,   ONLY: ice_growth_zerolayer, set_ice_temp_zerolayer, set_ice_temp_zero_nogradients
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_dbg_nml,             ONLY: idbg_mxmn, idbg_val
  USE mo_ice_fem_utils,       ONLY: fem_ice_wrap, ice_advection, ice_advection_vla, ice_ocean_stress
  USE mo_ice_fem_init,        ONLY: init_fem_wgts, destruct_fem_wgts, ice_fem_grid_init, ice_fem_grid_post
  USE mo_ice_init,            ONLY: ice_init_fem
!  USE mo_grid_config,         ONLY: n_dom   ! restrict sea-ice model to the global domain for the time being
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_ice_fast, timer_ice_slow
  USE mo_fortran_tools,       ONLY: assign_if_present
  USE mo_io_config,           ONLY: lnetcdf_flt64_output

  IMPLICIT NONE

  PRIVATE

  ! Public interface

  ! Definition of forcing types
  ! public types
  ! contained in mo_sea_ice_types

  ! public subroutines
  PUBLIC :: construct_sea_ice
  PUBLIC :: destruct_sea_ice
  PUBLIC :: construct_atmos_for_ocean
  PUBLIC :: construct_atmos_fluxes
  PUBLIC :: destruct_atmos_for_ocean

  PUBLIC :: ice_init
  PUBLIC :: ice_zero
!  PUBLIC :: set_ice_albedo
!  PUBLIC :: sum_fluxes
!  PUBLIC :: ave_fluxes
  PUBLIC :: ice_fast
  PUBLIC :: ice_slow
  PUBLIC :: upper_ocean_TS
  PUBLIC :: ice_conc_change
  PUBLIC :: ice_clean_up_thd, ice_clean_up_dyn
  PUBLIC :: calc_bulk_flux_ice
  PUBLIC :: calc_bulk_flux_oce
  PUBLIC :: update_ice_statistic, compute_mean_ice_statistics, reset_ice_statistics
  PUBLIC :: salt_content_in_surface
  PUBLIC :: energy_content_in_surface

  !to be put into namelist
  !  INTEGER :: i_no_ice_thick_class = 1

  CHARACTER(len=12)           :: str_module    = 'SeaIce'  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1         ! Level of detail for 1 line debug



CONTAINS

  !-------------------------------------------------------------------------
  !
  !> Constructor of sea-ice model, allocates all components and assigns zero.
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !
  SUBROUTINE construct_sea_ice(p_patch_3D, p_ice, i_no_ice_thick_class)
    TYPE(t_patch_3D),TARGET,INTENT(IN)    :: p_patch_3D
    TYPE (t_sea_ice),       INTENT(INOUT) :: p_ice
    INTEGER,                INTENT(IN)    :: i_no_ice_thick_class

    !Local variables
    !INTEGER i
    INTEGER :: ibits = DATATYPE_PACK16

    INTEGER :: alloc_cell_blocks, nblks_v, nblks_e, ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_sea_ice'

    TYPE(t_patch),POINTER    :: p_patch
    INTEGER                  :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    p_patch           => p_patch_3D%p_patch_2D(1)
    alloc_cell_blocks =  p_patch%alloc_cell_blocks
    nblks_v           =  p_patch%nblks_v
    nblks_e           =  p_patch%nblks_e

    p_ice%kice = i_no_ice_thick_class

    CALL add_var(ocean_restart_list, 'alb', p_ice%alb ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('alb', '', 'albedo of snow-ice system', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'Tsurf', p_ice%Tsurf ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('Tsurf', '', 'surface temperature', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'T1', p_ice%T1 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('T1', 'C', 'Temperature upper layer', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'T2', p_ice%T2 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('T2', 'C', 'Temperature lower layer', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'E1', p_ice%E1 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('E1', 'Jm/kg', 'Energy content upper layer', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'E2', p_ice%E2 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('E2', 'Jm/kg', 'Energy content lower layer', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'vol', p_ice%vol ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('vol', 'm^3', 'ice volume', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'vols', p_ice%vols ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('vols', 'm^3', 'snow volume', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'hi', p_ice%hi ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hi', 'm', 'ice thickness', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'hs', p_ice%hs ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hs', 'm', 'snow thickness', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'hiold', p_ice%hiold ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hiold', 'm', 'ice thickness (last timstep)', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'hsold', p_ice%hsold ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hsold', 'm', 'snow thickness (last timstep)', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'Qtop', p_ice%Qtop ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('Qtop', 'W/m^2', 'Energy flux available for surface melting', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'Qbot', p_ice%Qbot ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('Qbot', 'W/m^2', 'Energy flux at ice-ocean interface', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
 !  CALL add_var(ocean_restart_list, 'Qcond', p_ice%Qcond ,&
 !    &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
 !    &          t_cf_var('Qcond', 'W/m^2', 'Conductive heat flux through ice at bottom', datatype_flt),&
 !    &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
 !    &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
 !    &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'heatocei', p_ice%heatocei ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('heatocei', 'J', 'Energy to ocean when all ice is melted', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'snow_to_ice', p_ice%snow_to_ice ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('snow_to_ice', 'm', 'amount of snow that is transformed to ice', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'surfmelt', p_ice%surfmelt ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('surfmelt', 'm', 'surface melt water running into ocean', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'surfmeltT', p_ice%surfmeltT ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('surfmeltT', 'C', 'Mean temperature of surface melt water', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'evapwi', p_ice%evapwi ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('evapwi', 'kg/m^2', 'amount of evaporated water if no ice left', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'conc', p_ice%conc ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('conc', '', 'ice concentration in each ice class', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'draft', p_ice%draft ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('draft', 'm', 'water equivalent of ice and snow on ice covered area', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'ice_u_prog', p_ice%u_prog ,&
      &          GRID_UNSTRUCTURED_VERT, ZA_SURFACE, &
      &          t_cf_var('ice_u_prog', 'm/s', 'zonal velocity', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_VERTEX),&
      &          ldims=(/nproma,nblks_v/), lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'ice_v_prog', p_ice%v_prog ,&
      &          GRID_UNSTRUCTURED_VERT, ZA_SURFACE, &
      &          t_cf_var('ice_v', 'm/s', 'meridional velocity', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_VERTEX),&
      &          ldims=(/nproma,nblks_v/), lrestart_cont=.TRUE.)
      
    CALL add_var(ocean_restart_list, 'ice_u', p_ice%u ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('ice_u', 'm/s', 'zonal velocity', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_diag"),&
      &          lrestart_cont=.FALSE.)
    CALL add_var(ocean_restart_list, 'ice_v', p_ice%v ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('ice_v', 'm/s', 'meridional velocity', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_diag"),&
      &          lrestart_cont=.FALSE.)

    CALL add_var(ocean_restart_list, 'ice_vn', p_ice%vn_e ,&
      &          GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, &
      &          t_cf_var('ice_vn', 'm/s', 'zonal velocity', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE),&
      &          ldims=(/nproma,nblks_e/),&
      &          lrestart_cont=.FALSE.)

    CALL add_var(ocean_restart_list, 'concSum', p_ice%concSum ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('concSum', '', 'total ice concentration', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'newice', p_ice%newice ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('newice', 'm', 'new ice growth in open water', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'totalsnowfall', p_ice%totalsnowfall ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('totalsnowfall', 'm', 'total snow fall on sea ice', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    
      CALL add_var(ocean_restart_list, 'zUnderIce', p_ice%zUnderIce ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('zUnderIce', 'm', 'water in upper ocean grid cell below ice', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
      CALL add_var(ocean_restart_list, 'draftave', p_ice%draftave ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('draftave', 'm', 'averaged water equivalent of ice and snow on grid area', &
      &                   datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    ALLOCATE(p_ice%hi_lim(i_no_ice_thick_class), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for hi_lim failed')
    END IF


    IF(p_ice%kice==1)THEN
      p_ice%hi_lim = 0.0_wp
    ELSEIF(p_ice%kice==8)THEN
      p_ice%hi_lim(:)=(/ 0.0_wp, 0.1_wp, 0.3_wp, 0.7_wp, 1.1_wp, 1.5_wp, 2.0_wp, 2.5_wp /)
    ENDIF

    IF ( i_ice_dyn == 1 ) THEN ! AWI dynamics
      CALL init_fem_wgts(p_patch_3D)
    ENDIF

    CALL add_var(ocean_default_list, 'zHeatOceI', p_ice%zHeatOceI,GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('zHeatOceI', 'W/m^2', 'Oceanic Heat flux', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! add accumulated fields
    CALL add_var(ocean_default_list, 'hi_acc', p_ice%acc%hi ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hi_acc', 'm', 'ice thickness', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_default"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_default_list, 'hs_acc', p_ice%acc%hs ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hs_acc', 'm', 'snow thickness', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_default"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_default_list, 'conc_acc', p_ice%acc%conc ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('conc_acc', '', 'ice concentration in each ice class', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_default"),&
      &          lrestart_cont=.TRUE.)
    IF ( i_ice_dyn == 1 ) THEN ! AWI dynamics
      CALL add_var(ocean_default_list, 'ice_u_acc', p_ice%acc%u ,&
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('ice_u_acc', 'm/s', 'zonal velocity', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_default"),&
        &          lrestart_cont=.FALSE.)
      CALL add_var(ocean_default_list, 'ice_v_acc', p_ice%acc%v ,&
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('ice_v_acc', 'm/s', 'meridional velocity', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_default"),&
        &          lrestart_cont=.FALSE.)
    ELSE  !  do not write into ice_default
      CALL add_var(ocean_default_list, 'ice_u_acc', p_ice%acc%u ,&
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('ice_u_acc', 'm/s', 'zonal velocity', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_diag"),&
        &          lrestart_cont=.FALSE.)
      CALL add_var(ocean_default_list, 'ice_v_acc', p_ice%acc%v ,&
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('ice_v_acc', 'm/s', 'meridional velocity', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_diag"),&
        &          lrestart_cont=.FALSE.)
    ENDIF

    CALL message(TRIM(routine), 'end' )

    CALL construct_sea_ice_budgets(p_patch_3D,p_ice%budgets, ocean_default_list)
  END SUBROUTINE construct_sea_ice

  
  SUBROUTINE construct_sea_ice_budgets(patch_3d,budgets, varlist)
    TYPE(t_patch_3d) :: patch_3d
    TYPE(t_sea_ice_budgets) :: budgets
    TYPE(t_var_list)        :: varlist

    INTEGER :: alloc_cell_blocks
    TYPE(t_patch), POINTER :: patch
    INTEGER :: ibits = DATATYPE_PACK16
    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    patch             => patch_3D%p_patch_2D(1)
    alloc_cell_blocks =  patch%alloc_cell_blocks

    CALL add_var(varlist, 'salt_00', budgets%salt_00 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('salt_00', 'kg', '', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_budgets"))
!   CALL add_var(varlist, 'salt_01', budgets%salt_00 ,&
!     &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
!     &          t_cf_var('salt_01', 'kg', '', datatype_flt),&
!     &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
!     &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_budgets"))
!   CALL add_var(varlist, 'salt_02', budgets%salt_00 ,&
!     &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
!     &          t_cf_var('salt_02', 'kg', '', datatype_flt),&
!     &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
!     &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_budgets"))
!   CALL add_var(varlist, 'salt_03', budgets%salt_00 ,&
!     &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
!     &          t_cf_var('salt_03', 'kg', '', datatype_flt),&
!     &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
!     &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_budgets"))
  END SUBROUTINE construct_sea_ice_budgets
  
  !-------------------------------------------------------------------------
  !
  !> Destructor of sea-ice model, deallocates all components.
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !
  SUBROUTINE destruct_sea_ice(p_ice)
    TYPE (t_sea_ice),  INTENT (INOUT) :: p_ice
    !Local variables
    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_sea_ice'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    DEALLOCATE(p_ice%hi_lim, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for hi_lim failed')
    END IF

    IF ( i_ice_dyn == 1 ) THEN ! AWI dynamics
      CALL destruct_fem_wgts
    ENDIF

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE destruct_sea_ice

  !-------------------------------------------------------------------------
  !
  !>
  !! Constructor of atmospheric reprsentation  in ocean.
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-07)
  !
  SUBROUTINE construct_atmos_for_ocean(p_patch, p_as)
    !
    TYPE(t_patch),                INTENT(IN):: p_patch
    TYPE(t_atmos_for_ocean ), INTENT(INOUT) :: p_as

    ! Local variables
    INTEGER :: alloc_cell_blocks, ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_atmos_for_ocean'

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    alloc_cell_blocks = p_patch%alloc_cell_blocks

    ALLOCATE(p_as%tafo(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'tafo')

    ALLOCATE(p_as%ftdew(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'ftdew')

    ALLOCATE(p_as%fclou(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'fclou')

    ALLOCATE(p_as%fu10(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'fu10')

    ALLOCATE(p_as%fswr(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'fswr')

    ALLOCATE(p_as%pao(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'pao')

    ALLOCATE(p_as%u(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'u')

    ALLOCATE(p_as%v(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'v')

    ALLOCATE(p_as%precip(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'precip')

    ALLOCATE(p_as%evap(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'evap')

    ALLOCATE(p_as%runoff(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'runoff')

    ALLOCATE(p_as%topBoundCond_windStress_u(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'topBoundCond_windStress_u')

    ALLOCATE(p_as%topBoundCond_windStress_v(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'topBoundCond_windStress_v')

    ALLOCATE(p_as%FrshFlux_Precipitation(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'FrshFlux_Precipitation')

    ALLOCATE(p_as%FrshFlux_Runoff(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'FrshFlux_Runoff')

    ALLOCATE(p_as%data_surfRelax_Temp(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'data_surfRelax_Temp')

    ALLOCATE(p_as%data_surfRelax_Salt(nproma,alloc_cell_blocks), STAT=ist)
    CALL finish_unless_allocate(ist,routine,'data_surfRelax_Salt')

    p_as%tafo  (:,:)                    = 0.0_wp
    p_as%ftdew (:,:)                    = 0.0_wp
    p_as%fclou (:,:)                    = 0.0_wp
    p_as%fu10  (:,:)                    = 0.0_wp
    p_as%fswr  (:,:)                    = 0.0_wp
    p_as%pao   (:,:)                    = 0.0_wp
    p_as%u     (:,:)                    = 0.0_wp
    p_as%v     (:,:)                    = 0.0_wp
    p_as%precip(:,:)                    = 0.0_wp
    p_as%evap  (:,:)                    = 0.0_wp
    p_as%runoff(:,:)                    = 0.0_wp
    p_as%topBoundCond_windStress_u(:,:) = 0.0_wp
    p_as%topBoundCond_windStress_v(:,:) = 0.0_wp
    p_as%FrshFlux_Precipitation(:,:)    = 0.0_wp
    p_as%FrshFlux_Runoff(:,:)           = 0.0_wp
    p_as%data_surfRelax_Temp(:,:)       = 0.0_wp
    p_as%data_surfRelax_Salt(:,:)       = 0.0_wp

    CALL message(TRIM(routine), 'end')

  END SUBROUTINE construct_atmos_for_ocean
  !-------------------------------------------------------------------------
  !
  !>
  !!  Destructor of atmospheric reprsentation  in ocean.
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011)
  !
  SUBROUTINE destruct_atmos_for_ocean(p_as)
    !
    TYPE(t_atmos_for_ocean ), INTENT(INOUT) :: p_as

    ! Local variables
    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_atmos_for_ocean'

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )


    DEALLOCATE(p_as%tafo, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for tafo failed')
    END IF
    DEALLOCATE(p_as%ftdew, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for ftdew failed')
    END IF
    DEALLOCATE(p_as%fclou, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for fclou failed')
    END IF

    DEALLOCATE(p_as%fu10, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for fu10 failed')
    END IF

    DEALLOCATE(p_as%fswr, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for fswr failed')
    END IF

    DEALLOCATE(p_as%pao, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for pao failed')
    END IF

    DEALLOCATE(p_as%u, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for u failed')
    END IF
    DEALLOCATE(p_as%v, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for v failed')
    END IF

    DEALLOCATE(p_as%precip, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for precip failed')
    END IF
    DEALLOCATE(p_as%evap, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for evap failed')
    END IF
    DEALLOCATE(p_as%runoff, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for runoff failed')
    END IF

    CALL message(TRIM(routine), 'end')

  END SUBROUTINE destruct_atmos_for_ocean
  !-------------------------------------------------------------------------
  !
  !>
  !! Constructor of atmos fluxes for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011)
  !
  SUBROUTINE construct_atmos_fluxes(p_patch, atmos_fluxes, i_no_ice_thick_class)
    !
    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_atmos_fluxes ), INTENT(INOUT) :: atmos_fluxes
    INTEGER,               INTENT(IN)    :: i_no_ice_thick_class
    ! Local variables
    INTEGER :: ibits = DATATYPE_PACK16
    INTEGER :: alloc_cell_blocks, ist

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_atmos_fluxes'
    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    alloc_cell_blocks = p_patch%alloc_cell_blocks

    ! TODO: cleanup unused fluxes, sort in ice_diag
    ! Variables with 3 dimensions: fluxes over ice-covered part, second dimension is ice class

    CALL add_var(ocean_default_list, 'atmos_fluxes_lat', atmos_fluxes%lat,                            &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_lat', 'W/m2', 'atmos_fluxes_lat', datatype_flt),            &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_sens', atmos_fluxes%sens,                          &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_sens', 'W/m2', 'atmos_fluxes_sens', datatype_flt),          &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_LWout', atmos_fluxes%LWout,                        &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_LWout', 'W/m2', 'atmos_fluxes_LWout', datatype_flt),        &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_LWnet', atmos_fluxes%LWnet,                        &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_LWnet', 'W/m2', 'atmos_fluxes_LWnet', datatype_flt),        &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_SWnet', atmos_fluxes%SWnet,                        &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_SWnet', 'W/m2', 'atmos_fluxes_SWnet', datatype_flt),        &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_bot', atmos_fluxes%bot,                            &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_bot', 'W/m2', 'atmos_fluxes_bot', datatype_flt),            &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_dsensdT', atmos_fluxes%dsensdT,                    &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_dsensdT', 'W/m2/K', 'atmos_fluxes_dsensdT', datatype_flt),  &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_dlatdT', atmos_fluxes%dlatdT,                      &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_dlatdT', 'W/m2/K', 'atmos_fluxes_dlatdT', datatype_flt),    &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_dLWdT', atmos_fluxes%dLWdT,                        &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_dLWdT', 'W/m2/K', 'atmos_fluxes_dLWdT', datatype_flt),      &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! Variables with 2 dimensions: fluxes over icefree part in sea ice model

    CALL add_var(ocean_default_list, 'atmos_fluxes_latw', atmos_fluxes%latw,                          &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_latw', 'W/m2', 'atmos_fluxes_latw', datatype_flt),          &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_sensw', atmos_fluxes%sensw,                        &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_sensw', 'W/m2', 'atmos_fluxes_sensw', datatype_flt),        &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_LWoutw', atmos_fluxes%LWoutw,                      &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_LWoutw', 'W/m2', 'atmos_fluxes_LWoutw', datatype_flt),      &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_LWnetw', atmos_fluxes%LWnetw,                      &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_LWnetw', 'W/m2', 'atmos_fluxes_LWnetw', datatype_flt),      &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_SWnetw', atmos_fluxes%SWnetw,                      &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_SWnetw', 'W/m2', 'atmos_fluxes_SWnetw', datatype_flt),      &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_LWin', atmos_fluxes%LWin,                          &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_LWin', 'W/m2', 'atmos_fluxes_LWin', datatype_flt),          &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_rprecw', atmos_fluxes%rprecw,                      &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_rprecw', 'm/s', 'atmos_fluxes_rprecw', datatype_flt),       &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL add_var(ocean_default_list, 'atmos_fluxes_rpreci', atmos_fluxes%rpreci,                      &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
      &          t_cf_var('atmos_fluxes_rpreci', 'm/s', 'atmos_fluxes_rpreci', datatype_flt),       &
      &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
      
    ! Coupling fluxes must go into restart file:
    IF (is_coupled_run()) THEN

      CALL add_var(ocean_restart_list, 'atmos_fluxes_stress_x', atmos_fluxes%stress_x,                  &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
        &          t_cf_var('atmos_fluxes_stress_x', 'Pa',   'atmos_fluxes_stress_x', datatype_flt),  &
        &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),                      &
        &          lrestart_cont=.TRUE.)
     
      CALL add_var(ocean_restart_list, 'atmos_fluxes_stress_y', atmos_fluxes%stress_y,                  &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
        &          t_cf_var('atmos_fluxes_stress_y', 'Pa',   'atmos_fluxes_stress_y', datatype_flt),  &
        &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),                      &
        &          lrestart_cont=.TRUE.)
     
      CALL add_var(ocean_restart_list, 'atmos_fluxes_stress_xw', atmos_fluxes%stress_xw,                &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
        &          t_cf_var('atmos_fluxes_stress_xw', 'Pa',   'atmos_fluxes_stress_xw', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),                      &
        &          lrestart_cont=.TRUE.)
     
      CALL add_var(ocean_restart_list, 'atmos_fluxes_stress_yw', atmos_fluxes%stress_yw,                &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
        &          t_cf_var('atmos_fluxes_stress_yw', 'Pa',   'atmos_fluxes_stress_yw', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),                      &
        &          lrestart_cont=.TRUE.)

    ELSE

      CALL add_var(ocean_default_list, 'atmos_fluxes_stress_x', atmos_fluxes%stress_x,                  &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
        &          t_cf_var('atmos_fluxes_stress_x', 'Pa',   'atmos_fluxes_stress_x', datatype_flt),  &
        &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
     
      CALL add_var(ocean_default_list, 'atmos_fluxes_stress_y', atmos_fluxes%stress_y,                  &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
        &          t_cf_var('atmos_fluxes_stress_y', 'Pa',   'atmos_fluxes_stress_y', datatype_flt),  &
        &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
     
      CALL add_var(ocean_default_list, 'atmos_fluxes_stress_xw', atmos_fluxes%stress_xw,                &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
        &          t_cf_var('atmos_fluxes_stress_xw', 'Pa',   'atmos_fluxes_stress_xw', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
     
      CALL add_var(ocean_default_list, 'atmos_fluxes_stress_yw', atmos_fluxes%stress_yw,                &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                  &
        &          t_cf_var('atmos_fluxes_stress_yw', 'Pa',   'atmos_fluxes_stress_yw', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),             &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ENDIF  !  coupled

    !albedos need to go into the restart
    CALL add_var(ocean_restart_list, 'albvisdirw', atmos_fluxes%albvisdirw ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('albvisdirw', '', 'albvisdirw', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albvisdifw', atmos_fluxes%albvisdifw ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('albvisdifw', '', 'albvisdifw', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albnirdirw', atmos_fluxes%albnirdirw ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('albnirdirw', '', 'albnirdirw', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albnirdifw', atmos_fluxes%albnirdifw ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('albnirdifw', '', 'albnirdifw', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albvisdir', atmos_fluxes%albvisdir ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('albvisdir', '', 'albvisdir', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albvisdif', atmos_fluxes%albvisdif ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('albvisdif', '', 'albvisdif', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albnirdir', atmos_fluxes%albnirdir ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('albnirdir', '', 'albnirdir', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albnirdif', atmos_fluxes%albnirdif ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('albnirdif', '', 'albnirdif', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    ! Initialize with zero
    atmos_fluxes%counter = 0

    ! Initialise the albedos sensibly
    atmos_fluxes%albvisdir (:,:,:) = albi
    atmos_fluxes%albvisdif (:,:,:) = albi
    atmos_fluxes%albnirdir (:,:,:) = albi
    atmos_fluxes%albnirdif (:,:,:) = albi
    atmos_fluxes%albvisdirw(:,:) = albedoW_sim
    atmos_fluxes%albvisdifw(:,:) = albedoW_sim
    atmos_fluxes%albnirdirw(:,:) = albedoW_sim
    atmos_fluxes%albnirdifw(:,:) = albedoW_sim

    ! forcing of zonal component of velocity equation,
    CALL add_var(ocean_default_list,'atmos_fluxes_topBC_windStress_u', atmos_fluxes%topBoundCond_windStress_u,     &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_topBC_windStress_u', '', 'atmos_fluxes_topBoundCond_windStress_u', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! forcing of meridional component of velocity equation,
    CALL add_var(ocean_default_list,'atmos_fluxes_topBC_windStress_v', atmos_fluxes%topBoundCond_windStress_v,     &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_topBC_windStress_v', '', 'atmos_fluxes_topBoundCond_windStress_v', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! relaxation fields {{{
    !
    ! forcing of temperature in vertical diffusion equation     [K*m/s]
    CALL add_var(ocean_default_list,'topBoundCond_Temp_vdiff  ',atmos_fluxes%topBoundCond_Temp_vdiff ,     &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_topBC_Temp_vdiff', '', 'atmos_fluxes_topBoundCond_Temp_vdiff', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! forcing of salinity in vertical diffusion equation        [psu*m/s]
    CALL add_var(ocean_default_list,'topBoundCond_Salt_vdiff  ',atmos_fluxes%topBoundCond_Salt_vdiff ,     & 
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_topBC_Salt_vdiff', '', 'atmos_fluxes_topBoundCond_Salt_vdiff', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! contains data to which salinity is relaxed                [psu]
    CALL add_var(ocean_default_list,'data_surfRelax_Salt      ',atmos_fluxes%data_surfRelax_Salt     ,     & 
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_data_surfRelax_Salt', '', 'atmos_fluxes_data_surfRelax_Salt', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    ! contains data to which temperature is relaxed             [K]
    CALL add_var(ocean_default_list,'atmos_fluxes_data_surfRelax_Temp', atmos_fluxes%data_surfRelax_Temp,     &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_data_surfRelax_Temp', '', 'atmos_fluxes_data_surfRelax_Temp', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
    CALL add_var(ocean_default_list, 'atmos_fluxes_HeatFlux_Relax', atmos_fluxes%HeatFlux_Relax , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('atmos_fluxes_HeatFlux_Relax', 'W/m2', 'atmos_fluxes_HeatFlux_Relax', datatype_flt),&
    &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
    CALL add_var(ocean_default_list, 'atmos_fluxes_FrshFlux_Relax', atmos_fluxes%FrshFlux_Relax , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('atmos_fluxes_FrshFlux_Relax', 'm/s', 'atmos_fluxes_FrshFlux_Relax', datatype_flt),&
    &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
    CALL add_var(ocean_default_list, 'atmos_fluxes_TempFlux_Relax', atmos_fluxes%TempFlux_Relax , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('atmos_fluxes_TempFlux_Relax', 'K/s', 'atmos_fluxes_TempFlux_Relax', datatype_flt),&
    &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
    CALL add_var(ocean_default_list, 'atmos_fluxes_SaltFlux_Relax', atmos_fluxes%SaltFlux_Relax , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('atmos_fluxes_SaltFlux_Relax', 'psu/s', 'atmos_fluxes_SaltFlux_Relax', datatype_flt),&
    &          grib2_var(255, 255, 255, ibits          , GRID_UNSTRUCTURED, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
    ! }}}
      
    ! Coupling fluxes must go into restart file:
    IF (is_coupled_run()) THEN

      ! surface short wave heat flux                              [W/m2]
      CALL add_var(ocean_restart_list,'atmos_fluxes_HeatFlux_ShortWave', atmos_fluxes%HeatFlux_ShortWave,         &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_HeatFlux_ShortWave', '[W/m2]', 'atmos_fluxes_HeatFlux_ShortWave', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),lrestart_cont=.TRUE.)
   
      ! surface long wave heat flux                               [W/m2]
      CALL add_var(ocean_restart_list,'atmos_fluxes_HeatFlux_LongWave', atmos_fluxes%HeatFlux_LongWave,           &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_HeatFlux_LongWave', '[W/m2]', 'atmos_fluxes_HeatFlux_LongWave', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),lrestart_cont=.TRUE.)
  
      ! surface sensible heat flux                                [W/m2]
      CALL add_var(ocean_restart_list,'atmos_fluxes_HeatFlux_Sensible', atmos_fluxes%HeatFlux_Sensible,           &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_HeatFlux_Sensible', '[W/m2]', 'atmos_fluxes_HeatFlux_Sensible', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),lrestart_cont=.TRUE.)
  
      ! surface latent heat flux                                  [W/m2]
      CALL add_var(ocean_restart_list,'atmos_fluxes_HeatFlux_Latent', atmos_fluxes%HeatFlux_Latent,               &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_HeatFlux_Latent', '[W/m2]', 'atmos_fluxes_HeatFlux_Latent', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),lrestart_cont=.TRUE.)
      ! total heat flux                                  [W/m2]
      CALL add_var(ocean_restart_list,'atmos_fluxes_HeatFlux_Total', atmos_fluxes%HeatFlux_Total,               &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_HeatFlux_Total', '[m/s]', 'atmos_fluxes_HeatFlux_Total', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),lrestart_cont=.TRUE.)
  
      ! total precipitation flux                                  [m/s]
      CALL add_var(ocean_restart_list,'atmos_fluxes_FrshFlux_Precipitation', atmos_fluxes%FrshFlux_Precipitation, &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_FrshFlux_Precipitation','[m/s]','atmos_fluxes_FrshFlux_Precipitation',datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),lrestart_cont=.TRUE.)
  
      ! total snow flux                                           [m/s]
      CALL add_var(ocean_restart_list,'atmos_fluxes_FrshFlux_SnowFall', atmos_fluxes%FrshFlux_SnowFall,           &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_FrshFlux_SnowFall', '[m/s]', 'atmos_fluxes_FrshFlux_SnowFall', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),lrestart_cont=.TRUE.)
  
      ! evaporation flux                                          [m/s]
      CALL add_var(ocean_restart_list,'atmos_fluxes_FrshFlux_Evaporation', atmos_fluxes%FrshFlux_Evaporation,     &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_FrshFlux_Evaporation', '[m/s]', 'atmos_fluxes_FrshFlux_Evaporation', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),lrestart_cont=.TRUE.)
  
      ! river runoff flux                                         [m/s]
      CALL add_var(ocean_restart_list,'atmos_fluxes_FrshFlux_Runoff', atmos_fluxes%FrshFlux_Runoff, &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_FrshFlux_Runoff', '[m/s]', 'atmos_fluxes_FrshFlux_Runoff', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),lrestart_cont=.TRUE.)

    ELSE

      ! surface short wave heat flux                              [W/m2]
      CALL add_var(ocean_default_list,'atmos_fluxes_HeatFlux_ShortWave', atmos_fluxes%HeatFlux_ShortWave,         &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_HeatFlux_ShortWave', '[W/m2]', 'atmos_fluxes_HeatFlux_ShortWave', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))  
   
      ! surface long wave heat flux                               [W/m2]
      CALL add_var(ocean_default_list,'atmos_fluxes_HeatFlux_LongWave', atmos_fluxes%HeatFlux_LongWave,           &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_HeatFlux_LongWave', '[W/m2]', 'atmos_fluxes_HeatFlux_LongWave', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))  

      ! surface sensible heat flux                                [W/m2]
      CALL add_var(ocean_default_list,'atmos_fluxes_HeatFlux_Sensible', atmos_fluxes%HeatFlux_Sensible,           &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_HeatFlux_Sensible', '[W/m2]', 'atmos_fluxes_HeatFlux_Sensible', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))  

      ! surface latent heat flux                                  [W/m2]
      CALL add_var(ocean_default_list,'atmos_fluxes_HeatFlux_Latent', atmos_fluxes%HeatFlux_Latent,               &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_HeatFlux_Latent', '[W/m2]', 'atmos_fluxes_HeatFlux_Latent', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))  
      ! total heat flux                                  [W/m2]
      CALL add_var(ocean_default_list,'atmos_fluxes_HeatFlux_Total', atmos_fluxes%HeatFlux_Total,               &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_HeatFlux_Total', '[W/m2]', 'atmos_fluxes_HeatFlux_Total', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))  

      ! total precipitation flux                                  [m/s]
      CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_Precipitation', atmos_fluxes%FrshFlux_Precipitation, &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_FrshFlux_Precipitation','[m/s]','atmos_fluxes_FrshFlux_Precipitation',datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))  

      ! total snow flux                                           [m/s]
      CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_SnowFall', atmos_fluxes%FrshFlux_SnowFall,           &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_FrshFlux_SnowFall', '[m/s]', 'atmos_fluxes_FrshFlux_SnowFall', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))  

      ! evaporation flux                                          [m/s]
      CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_Evaporation', atmos_fluxes%FrshFlux_Evaporation,     &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_FrshFlux_Evaporation', '[m/s]', 'atmos_fluxes_FrshFlux_Evaporation', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))  

      ! river runoff flux                                         [m/s]
      CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_Runoff', atmos_fluxes%FrshFlux_Runoff, &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('atmos_fluxes_FrshFlux_Runoff', '[m/s]', 'atmos_fluxes_FrshFlux_Runoff', datatype_flt),&
        &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))  

    ENDIF  !  coupled

    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_TotalSalt', atmos_fluxes%FrshFlux_TotalSalt, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_TotalSalt', '[m/s]', 'atmos_fluxes_FrshFlux_TotalSalt', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_TotalOcean', atmos_fluxes%FrshFlux_TotalOcean, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_TotalOcean', '[m/s]', 'atmos_fluxes_FrshFlux_TotalOcean', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_TotalIce', atmos_fluxes%FrshFlux_TotalIce, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_TotalIce', '[m/s]', 'atmos_fluxes_FrshFlux_TotalIce', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_VolumeIce', atmos_fluxes%FrshFlux_VolumeIce, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_VolumeIce', '[m/s]', 'atmos_fluxes_FrshFlux_VolumeIce', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
    !  atmos_fluxes%FrshFlux_VolumeTotal is zero due to disassociated pointer in ocean_bulk
    CALL add_var(ocean_default_list,'atmos_fluxes_FrshFlux_VolumeTotal', atmos_fluxes%FrshFlux_VolumeTotal, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('atmos_fluxes_FrshFlux_VolumeTotal', '[m/s]', 'atmos_fluxes_FrshFlux_VolumeTotal', datatype_flt),&
      &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))
   CALL add_var(ocean_default_list,'atmos_flux_cellThicknessUnderIce', atmos_fluxes%cellThicknessUnderIce, &
     &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
     &          t_cf_var('atmos_flux_cellThicknessUnderIce', 'm', 'atmos_flux_cellThicknessUnderIce', datatype_flt),&
     &          grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
     &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"))

    CALL message(TRIM(routine), 'end' )
  END SUBROUTINE construct_atmos_fluxes

  !-------------------------------------------------------------------------
  !
  !> ice_init
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE ice_init( p_patch_3D, p_os, ice, cellThicknessUnderIce)
    TYPE(t_patch_3D), TARGET, INTENT(in)  :: p_patch_3D
    TYPE(t_hydro_ocean_state)             :: p_os
    TYPE (t_sea_ice),      INTENT (INOUT) :: ice
    REAL(wp), DIMENSION(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks), &
      & INTENT(OUT)                       :: cellThicknessUnderIce

    !local variables
    REAL(wp), DIMENSION(nproma,ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) :: &
      & Tinterface, & ! temperature at snow-ice interface
      & Tfw           ! Ocean freezing temperature [C]

    TYPE(t_patch), POINTER                :: p_patch
    TYPE(t_patch_vert), POINTER           :: p_patch_vert

    !INTEGER i,j,k      ! counter for loops
    INTEGER k !, jb, jc, i_startidx_c, i_endidx_c! counter for loops
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:ice_init'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    p_patch => p_patch_3D%p_patch_2D(1)
    p_patch_vert => p_patch_3D%p_patch_1D(1)

    !Constructor basic init already done at this point
    !   CALL alloc_mem_commo_ice (ice, atmos_fluxes, atmos_fluxesAve)
    !   CALL ice_zero            (ice, atmos_fluxes, atmos_fluxesAve)

    ! FORALL(i=1:nproma, j=1:p_patch%alloc_cell_blocks, k=1:ice%kice)
    !    ice% hi    (i,j,k) = sictho (i,j)
    !    ice% hs    (i,j,k) = sicsno (i,j)
    ! END FORALL

    ! LL Note: this needs to be rewritten using subsets !
    IF ( no_tracer < 2 .OR. use_constant_tfreez ) THEN
      Tfw(:,:,:) = Tf
    ELSE
      DO k=1,ice%kice
        Tfw(:,k,:) = -mu*p_os%p_prog(nold(1))%tracer(:,1,:,2)
      ENDDO
    ENDIF

    ice% Tsurf(:,:,:)  = Tf
    ice% T1   (:,:,:)  = Tf
    ice% T2   (:,:,:)  = Tf
    ice% conc (:,:,:)  = 0.0_wp

    ! Stupid initialisation trick for Levitus initialisation
    IF (use_IceInitialization_fromTemperature) THEN
      WHERE (p_os%p_prog(nold(1))%tracer(:,1,:,1) <= init_analytic_temp_under_ice .and. v_base%lsm_c(:,1,:) <= sea_boundary )
        ice%hi(:,1,:)   = init_analytic_hi_param
        ice%hs(:,1,:)   = init_analytic_hs_param
        ice%conc(:,1,:) = init_analytic_conc_param
      ENDWHERE
    ! or constant initialization for ice, snow and concentration
    ELSE
      WHERE (v_base%lsm_c(:,1,:) <= sea_boundary )
        ice%hi(:,1,:)    = init_analytic_hi_param
        ice%hs(:,1,:)    = init_analytic_hs_param
        ice%conc(:,1,:)  = init_analytic_conc_param
      ENDWHERE
    ENDIF

    WHERE(ice% hi(:,:,:) > 0.0_wp)
      ice% Tsurf (:,:,:) = Tfw(:,:,:)
      ice% T1    (:,:,:) = Tfw(:,:,:)
      ice% T2    (:,:,:) = Tfw(:,:,:)
      Tinterface (:,:,:) = (Tfw(:,:,:) * (ki/ks * ice%hs(:,:,:)/ice%hi(:,:,:)) &
        &                + ice%Tsurf(:,:,:)) / (1.0_wp+ki/ks * ice%hs(:,:,:)/ice%hi(:,:,:))
      ice% conc  (:,:,:) = ice%conc(:,:,:)/REAL(ice%kice,wp)
      ice% T1    (:,:,:) = Tfw(:,:,:) + 2._wp/3._wp*(Tinterface(:,:,:)-Tfw(:,:,:))
      ice% T2    (:,:,:) = Tfw(:,:,:) + 1._wp/3._wp*(Tinterface(:,:,:)-Tfw(:,:,:))
      ice% draft (:,:,:) = (rhos * ice%hs(:,:,:) + rhoi * ice%hi(:,:,:)) / rho_ref
    END WHERE

    ice%concSum(:,:)   = SUM(ice%conc(:,:,:), 2)
    ice%draftave (:,:) = sum(ice%draft(:,:,:) * ice%conc(:,:,:),2)
    ice%zUnderIce(:,:) = p_patch_vert%prism_thick_flat_sfc_c(:,1,:) +  p_os%p_prog(nold(1))%h(:,:) - ice%draftave(:,:)

    cellThicknessUnderIce (:,:) = ice%zUnderIce(:,:)

    IF ( i_ice_dyn == 1 ) THEN ! AWI dynamics
      CALL ice_fem_grid_init(p_patch_3D)
      CALL ice_init_fem
      CALL ice_fem_grid_post(p_patch)
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('IceInit: hi       ' ,ice%hi       ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: hs       ' ,ice%hs       ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: conc     ' ,ice%conc     ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: draft    ' ,ice%draft    ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: draftave ' ,ice%draftave ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: zUnderIce' ,ice%zUnderIce,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: Tfw      ' ,Tfw          ,str_module, idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------
      
    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE ice_init
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! !  ice_fast: Ice routines for atmospheric time step. Sets air-ice fluxes and
  !!    calculates the development of the ice temperature field
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE ice_fast(i_startidx_c, i_endidx_c, nbdim, kice, pdtime, &
            &   Tsurf,          & ! Surface temperature [degC]
            &   T1,             & ! Temperature of upper layer [degC]
            &   T2,             & ! Temperature of lower layer [degC]
            &   hi,             & ! Ice thickness
            &   hs,             & ! Snow thickness
            &   Qtop,           & ! Energy flux available for surface melting [W/m2]
            &   Qbot,           & ! Energy flux available for bottom melting [W/m2]
            &   SWnet,          & ! Net shortwave flux [W/m^2]
            &   nonsolar,       & ! Latent and sensible heat flux and longwave radiation [W/m^2]
            &   dnonsolardT,    & ! Derivative of non-solar fluxes w.r.t. temperature [W/m^2/K]
            &   Tfw,            & ! Freezing temperature of the ocean
            &   albvisdir,      & ! Albedo VIS, direct/parallel
            &   albvisdif,      & ! Albedo VIS, diffuse
            &   albnirdir,      & ! Albedo NIR, direct/parallel
            &   albnirdif,      & ! Albedo NIR, diffuse
            &   doy)              ! Day of the year

    INTEGER, INTENT(IN)    :: i_startidx_c, i_endidx_c, nbdim, kice
    REAL(wp),INTENT(IN)    :: pdtime
    REAL(wp),INTENT(INOUT) :: Tsurf      (nbdim,kice)
    REAL(wp),INTENT(INOUT) :: T1         (nbdim,kice)
    REAL(wp),INTENT(INOUT) :: T2         (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hi         (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hs         (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qtop       (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qbot       (nbdim,kice)
    REAL(wp),INTENT(IN)    :: SWnet      (nbdim,kice)
    REAL(wp),INTENT(IN)    :: nonsolar   (nbdim,kice)
    REAL(wp),INTENT(IN)    :: dnonsolardT(nbdim,kice)
    REAL(wp),INTENT(IN)    :: Tfw        (nbdim)
    REAL(wp),INTENT(OUT)   :: albvisdir  (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: albvisdif  (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: albnirdir  (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: albnirdif  (nbdim,kice)

    INTEGER, OPTIONAL,INTENT(IN)  :: doy

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:ice_fast'

    !-------------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_ice_fast)

    ! #achim
    SELECT CASE (i_ice_therm)
    CASE (1)
      IF (use_no_flux_gradients) THEN
        CALL set_ice_temp_zero_nogradients(i_startidx_c, i_endidx_c, nbdim, kice, i_ice_therm, pdtime, &
              &   Tsurf,          &
              &   hi,             &
              &   hs,             &
              &   Qtop,           &
              &   Qbot,           &
              &   SWnet,          &
              &   nonsolar,       &
              &   dnonsolardT,    &
              &   Tfw)
      ELSE
        CALL set_ice_temp_zerolayer(i_startidx_c, i_endidx_c, nbdim, kice, i_ice_therm, pdtime, &
              &   Tsurf,          &
              &   hi,             &
              &   hs,             &
              &   Qtop,           &
              &   Qbot,           &
              &   SWnet,          &
              &   nonsolar,       &
              &   dnonsolardT,    &
              &   Tfw)
      ENDIF
    CASE (2)
      CALL set_ice_temp_winton(i_startidx_c, i_endidx_c, nbdim, kice, pdtime, &
            &   Tsurf,          &
            &   T1,             &
            &   T2,             &
            &   hi,             &
            &   hs,             &
            &   Qtop,           &
            &   Qbot,           &
            &   SWnet,          &
            &   nonsolar,       &
            &   dnonsolardT,    &
            &   Tfw)
    CASE (3)
      IF ( .NOT. PRESENT(doy) ) THEN
        CALL finish(TRIM(routine),'i_ice_therm = 3 not allowed in this context')
      ENDIF
      CALL set_ice_temp_zerolayer(i_startidx_c, i_endidx_c, nbdim, kice, i_ice_therm, pdtime, &
            &   Tsurf,          &
            &   hi,             &
            &   hs,             &
            &   Qtop,           &
            &   Qbot,           &
            &   SWnet,          &
            &   nonsolar,       &
            &   dnonsolardT,    &
            &   Tfw,            &
            &   doy=doy)
    CASE (4)
      WHERE ( hi(:,:) > 0._wp )
      Tsurf=min(0._wp, Tsurf + (SWnet+nonsolar + ki/hi*(Tf-Tsurf)) &
        &               / (ci*rhoi*hci_layer/pdtime-dnonsolardT+ki/hi))
      ELSEWHERE
        Tsurf(:,:) = Tf
      ENDWHERE
    END SELECT

    ! New albedo based on the new surface temperature
    CALL set_ice_albedo(i_startidx_c, i_endidx_c, nbdim, kice, Tsurf, hi, hs, &
      & albvisdir, albvisdif, albnirdir, albnirdif)

    IF (ltimer) CALL timer_stop(timer_ice_fast)

   END SUBROUTINE ice_fast
  !-------------------------------------------------------------------------------
  !
  !
  !>
  !! !  ice_slow: Ice routines for ocean time step. Calculates average of atmospheric
  ! !           time steps, ice velocity, ice growth rates and updates ice structure
  ! !           accordingly
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE ice_slow(p_patch_3D, p_os, ice, atmos_fluxes, p_op_coeff)
    TYPE(t_patch_3D), TARGET, INTENT(IN)    :: p_patch_3D
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os
    TYPE(t_sea_ice),          INTENT(INOUT) :: ice
    TYPE(t_atmos_fluxes),     INTENT(INOUT) :: atmos_fluxes
    TYPE(t_operator_coeff),   INTENT(IN)    :: p_op_coeff

    TYPE(t_patch),      POINTER :: p_patch
    TYPE(t_patch_vert), POINTER :: p_patch_vert
    TYPE(t_subset_range), POINTER :: all_cells
    !INTEGER :: jb, jc, i_startidx_c, i_endidx_c

    REAL(wp), DIMENSION(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) :: energyCheck, sst, zuipsnowf
    REAL(wp), POINTER             :: flat(:,:)
    INTEGER                       :: k, jb, jc, i_startidx_c, i_endidx_c

    !-------------------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_ice_slow)

    p_patch      => p_patch_3D%p_patch_2D(1)
    p_patch_vert => p_patch_3D%p_patch_1D(1)
    ! subset range pointer
    all_cells => p_patch%cells%all 
    flat      => p_patch_vert%prism_thick_flat_sfc_c(:,1,:)

    sst(:,:)        = 0.0_wp
    zuipsnowf(:,:)  = 0.0_wp

    CALL ice_zero       (ice)

    ice%hiold(:,:,:) = ice%hi(:,:,:)
    ice%hsold(:,:,:) = ice%hs(:,:,:)
    CALL dbg_print('IceSlow: hi before growth' ,ice%hi ,str_module,4, in_subset=p_patch%cells%owned)
    ! #achim
    IF      ( i_ice_therm == 2 ) THEN
      CALL ice_growth_winton    (p_patch, p_os, ice, atmos_fluxes%rpreci)!, atmos_fluxes%lat)
    ELSE IF ( i_ice_therm == 1 .OR. i_ice_therm == 3 ) THEN !2=zerolayer, 3=simple fluxes from dirk's thesis
      CALL ice_growth_zerolayer (p_patch, p_os, ice, atmos_fluxes%rpreci)
    END IF

    sst(:,:) = p_os%p_prog(nold(1))%tracer(:,1,:,1)
    !---  energy  -----------------------------------------------------------
    ! updates should be moved to respective routine ice_growth - done later in upperOceTS
    !  - 2015-01-19: ice growth/melt must not change zunderice yet
 !  ice%draftave (:,:) = (rhos * ice%hs(:,1,:) + rhoi * ice%hi(:,1,:)) * ice%conc(:,1,:) / rho_ref
 !  ice%zUnderIce(:,:) = flat(:,:) + p_os%p_prog(nold(1))%h(:,:) - ice%draftave(:,:)
    energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
      &             ice, sst(:,:), 0, info='AFT GROWTH')

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('IceSlow: hi after growth'   ,ice%hi       ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: Conc. after growth',ice%conc     ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: ice%u bef. dyn'    ,ice%u_prog   ,str_module, 4, in_subset=p_patch%verts%owned)
    CALL dbg_print('IceSlow: ice%v bef. dyn'    ,ice%v_prog   ,str_module, 4, in_subset=p_patch%verts%owned)
    CALL dbg_print('IceSlow: zUnderIce aft.gr',  ice%zUnderIce,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: heatOceI aftgrowth',ice%heatOceI ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: SST    aft. growth',sst          ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: energy aft. growth',energyCheck  ,str_module, 2, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    ! #slo# 2014-12: now use melt/growth energy heatOceI from below to change SST
    !  heatOceI=heatOceI*conc - was done later in upperOceTS - generates inconsistency
    !IF (i_Qio_type < 3) ice%heatOceI(:,1,:) = ice%heatOceI(:,1,:)*ice%conc(:,1,:)
    ! wrong - ice covered part only is used for heat exchange with ocean
    !IF (i_Qio_type == 3) ice%heatOceI(:,1,:) = ice%heatOceI(:,1,:)*ice%conc(:,1,:) 
    ice%heatOceI(:,1,:) = ice%heatOceI(:,1,:)*ice%conc(:,1,:)
    CALL dbg_print('IceSlow: heatOceI aftConCor',ice%heatOceI ,str_module, 4, in_subset=p_patch%cells%owned)

    ! #slo# 2015-01: now move sst-change back to surface module (mo_ocean_bulk)
  ! DO jb = all_cells%start_block, all_cells%end_block
  !   CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
  !   DO jc = i_startidx_c, i_endidx_c
  !     IF ( v_base%lsm_c(jc,1,jb) <= sea_boundary ) THEN
  !       sst(jc,jb) = p_os%p_prog(nold(1))%tracer(jc,1,jb,1)
  !       sst(jc,jb) = sst(jc,jb) + ice%heatOceI(jc,1,jb)*dtime/(clw*rho_ref*ice%zUnderIce(jc,jb))
  !       ! set new sst; heatOceI to zero
  !       p_os%p_prog(nold(1))%tracer(jc,1,jb,1) = sst(jc,jb)
  !       ice%heatOceI(jc,1,jb)                  = 0.0_wp
  !     ENDIF
  !   ENDDO
  ! ENDDO

    energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
      &             ice, sst(:,:), 0, info='AFT UPGROW')
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('IceSlow: heatOceI aftupGrow',ice%heatOceI ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: SST    aft. upGrow',sst          ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: energy aft. upGrow',energyCheck  ,str_module, 4, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    ! Ocean Heat Flux
    IF ( i_ice_therm >= 1 ) &
      & CALL upper_ocean_TS(p_patch,p_patch_vert,ice,p_os,atmos_fluxes)

    sst(:,:) = p_os%p_prog(nold(1))%tracer(:,1,:,1)  !  add corresponding flux or calculate zUnderIce
 !  ice%zUnderIce(:,:) = flat(:,:) + p_os%p_prog(nold(1))%h(:,:) &
 !    &                - (rhos * ice%hs(:,1,:) + rhoi * ice%hi(:,1,:)) * ice%conc(:,1,:) / rho_ref
    energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
      &             ice, sst(:,:), 0, info='AFT UOCETS')
    CALL dbg_print('IceSlow: zUnderIce aft.UPTS',ice%zUnderIce,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: energy aft.OceanTS',energyCheck  ,str_module, 2, in_subset=p_patch%cells%owned)

    ! #slo# 2014-12: now use energy atmos_fluxes%HeatFlux_Total from upperOceTS to change SST
    ! #slo# 2015-01: now move sst-change back to surface module (mo_ocean_bulk) using HeatFlux_Total
  ! DO jb = all_cells%start_block, all_cells%end_block
  !   CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
  !   DO jc = i_startidx_c, i_endidx_c
  !     IF ( v_base%lsm_c(jc,1,jb) <= sea_boundary ) THEN
  !       sst(jc,jb) = p_os%p_prog(nold(1))%tracer(jc,1,jb,1)
  !       sst(jc,jb) = sst(jc,jb) + atmos_fluxes%HeatFlux_Total(jc,jb)*dtime/(clw*rho_ref*ice%zUnderIce(jc,jb))
  !       ! set new sst; HeatFlux_Total to zero
  !       p_os%p_prog(nold(1))%tracer(jc,1,jb,1) = sst(jc,jb)
  !       atmos_fluxes%HeatFlux_Total(jc,jb)     = 0.0_wp
  !     ENDIF
  !   ENDDO
  ! ENDDO
    energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
      &             ice, sst(:,:), 0, info='AFT UPUPTS')
    CALL dbg_print('IceSlow: HeatTotal aft.UPTS',atmos_fluxes%HeatFlux_Total,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: SST    aft. upUPTS',sst          ,str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: energy aft. upUPTS',energyCheck  ,str_module, 2, in_subset=p_patch%cells%owned)

    ! Ice Concentration Change
    IF ( i_ice_therm >= 1 ) THEN
      CALL ice_conc_change(p_patch,ice,p_os)
    ELSE
      ! This is needed since vol and vols are calculated in ice_conc_change and cleanup will set all to zero!
      DO k=1,ice%kice
        ice%vol (:,k,:) = ice%hi(:,k,:)*ice%conc(:,k,:)*p_patch%cells%area(:,:)
        ice%vols(:,k,:) = ice%hs(:,k,:)*ice%conc(:,k,:)*p_patch%cells%area(:,:)
      ENDDO
    ENDIF

    ! #slo# 2015-01 - Test: update draft, zunderice now includes totalsnowfall as in ocets, inconsistent with h
    !               - update of draft/zunderice should be done whenever draft is changed
 !  ice%draft       (:,:,:) = (rhos * ice%hs(:,:,:) + rhoi * ice%hi(:,:,:)) / rho_ref
 !  ice%draftave    (:,:)   = sum(ice%draft(:,:,:) * ice%conc(:,:,:),2)
 !  ice%zUnderIce   (:,:)   = p_patch_vert%prism_thick_flat_sfc_c(:,1,:) &
 !    &                            + p_os%p_prog(nold(1))%h(:,:) &
 !    &                            - ice%draftave(:,:)& 
 !    &                            + ice%totalsnowfall(:,:)
    ! - Test: update zunderice only, includes totalsnowfall
    zuipsnowf    (:,:) = flat(:,:) + p_os%p_prog(nold(1))%h(:,:) &
      &                - (rhos * ice%hs(:,1,:) + rhoi * ice%hi(:,1,:)) * ice%conc(:,1,:) / rho_ref &
      &                            + ice%totalsnowfall(:,:)


    sst(:,:) = p_os%p_prog(nold(1))%tracer(:,1,:,1)  !  add corresponding flux or calculate zUnderIce
    ice%zUnderIce(:,:) = flat(:,:) + p_os%p_prog(nold(1))%h(:,:) &
      &                - (rhos * ice%hs(:,1,:) + rhoi * ice%hi(:,1,:)) * ice%conc(:,1,:) / rho_ref
    energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
      &             ice, sst(:,:), 0, info='AFT CONCCH')

    CALL dbg_print('IceSlow: energy aft. ConcCh',energyCheck  ,str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: zUnderIce a.ConcCh',ice%zUnderIce,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: zUI+snowf a.ConcCh',zuipsnowf,    str_module, 4, in_subset=p_patch%cells%owned)

    ! ocean stress calculated independent of ice dynamics
    CALL ice_ocean_stress( p_patch, atmos_fluxes, ice, p_os )

    CALL dbg_print('IceSlow: hi    bef.icedyn',ice%hi,       str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: hs    bef.icedyn',ice%hs,       str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: Conc. bef.icedyn',ice%conc     ,str_module, 3, in_subset=p_patch%cells%owned)

    IF ( i_ice_dyn >= 1 ) THEN
      ! AWI FEM model wrapper
      CALL fem_ice_wrap ( p_patch_3D, ice, p_os, atmos_fluxes, p_op_coeff, 1) ! not 1 but jstep should be passed to fem_ice_wrap as the last argument
                                                                              ! but it does not matter, usage is depreciated
      ! #vla# 2015-05 - debugging ice_advection routine
!      CALL dbg_print('IceSlow: p_ice%u'           ,ice%u_prog,             str_module, 3, in_subset=p_patch%verts%owned)
!      CALL dbg_print('IceSlow: p_ice%v'           ,ice%v_prog,             str_module, 3, in_subset=p_patch%verts%owned)
!      CALL dbg_print('IceSlow: ice%vol bef.adv'     ,ice%vol,                 str_module, 3, in_subset=p_patch%cells%owned)

      ! direct calculation of the sea-ice volume
!      ice_vol=0_wp
!      DO blockNo = p_patch%cells%owned%start_block, p_patch%cells%owned%end_block
!        CALL get_index_range(p_patch%cells%owned, blockNo, start_cell_index, end_cell_index)
!        !We are dealing with the surface layer first
!        DO jc =  start_cell_index, end_cell_index
!     ! northern hemisphere
!          IF (patch_2d%cells%center(jc,blockNo)%lat > equator) THEN
!            ice_vol  = ice_vol + p_patch%cells%area(jc,blockNo)*SUM(ice%vol(jc,:,blockNo))!*ice%conc(jc,:,blockNo))
!          ENDIF
!        END DO
!      END DO
!      IF (my_process_is_stdio()) then
!      write(0,*) global_sum_array(ice_vol)/1.0e9_wp
!      write(0,*) '-----------------------advection-------------------------'
!      ENDIF

!      CALL ice_advection( p_patch_3D, p_op_coeff, ice ) ! messy advection routine, bugs fixed; renamed as ice_advection_vla
      CALL ice_advection_vla( p_patch_3D, p_op_coeff, ice )

!      CALL dbg_print('IceSlow: ice%vol aft.adv'     ,ice%vol,                 str_module, 3, in_subset=p_patch%cells%owned)

    ELSE
      ice%u = 0._wp
      ice%v = 0._wp
    ENDIF

    CALL dbg_print('IceSlow: hi    bef.cleanup',ice%hi,       str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: hs    bef.cleanup',ice%hs,       str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: Conc. bef.cleanup',ice%conc     ,str_module, 3, in_subset=p_patch%cells%owned)

    ! the original routine ice_clean_up has been split into two:
    ! 1) ice_clean_up_dyn fixes possible overshoots in conc afther the advection step
    IF ( i_ice_dyn >= 1 ) THEN
        CALL ice_clean_up_dyn( p_patch_3D, ice )
    ENDIF
    ! ice_clean_up_thd
    ! 2) fixes undershoots in conc and limit sea ice thickness to seaice_limit of surface layer depth after changes
    !    due to the thermodynamic growth/melt;
    ! 3) Calculates the new freeboard.
    CALL ice_clean_up_thd( p_patch_3D, ice, atmos_fluxes, p_os )

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('IceSlow: hi endOf slow'     ,ice%hi,                 str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: hs endOf slow'     ,ice%hs,                 str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: ConcSumEndOf slow', ice%concSum,            str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: Conc.  EndOf slow', ice%conc,               str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: p_ice%u'           ,ice%u_prog,             str_module, 4, in_subset=p_patch%verts%owned)
    CALL dbg_print('IceSlow: p_ice%v'           ,ice%v_prog,             str_module, 4, in_subset=p_patch%verts%owned)
    CALL dbg_print('IceSlow: p_os%prog(nold)%vn',p_os%p_prog(nold(1))%vn,str_module, 5, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: p_os%prog(nnew)%vn',p_os%p_prog(nnew(1))%vn,str_module, 5, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: p_os%diag%u'       ,p_os%p_diag%u,          str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: p_os%diag%v'       ,p_os%p_diag%v,          str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: atmFlx%windStr-u' ,atmos_fluxes%topBoundCond_windStress_u,str_module,4,in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: atmFlx%windStr-v' ,atmos_fluxes%topBoundCond_windStress_v,str_module,4,in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    IF (ltimer) CALL timer_stop(timer_ice_slow)

  END SUBROUTINE ice_slow

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! !  ice_clean_up_dyn: Basic fix for overshoots in concentration that can appear due to ice convergence in areas with conc ~ 1
  !!
  SUBROUTINE ice_clean_up_dyn( p_patch_3D, p_ice)
    TYPE(t_patch_3D),TARGET,   INTENT(IN)    :: p_patch_3D
    TYPE(t_sea_ice),           INTENT(INOUT) :: p_ice

    ! Local variables
    ! patch
    TYPE(t_patch),POINTER                                                    :: p_patch

    ! subset range pointer
    p_patch      => p_patch_3D%p_patch_2D(1)

    ! Fix overshoots - ONLY for the one-ice-class case
    WHERE ( p_ice%conc(:,1,:) > 1._wp )
      p_ice%conc(:,1,:) = 1._wp

      ! New ice and snow thickness
      p_ice%hi   (:,1,:) = p_ice%vol (:,1,:)/( p_ice%conc(:,1,:)*p_patch%cells%area(:,:) )
      p_ice%hs   (:,1,:) = p_ice%vols(:,1,:)/( p_ice%conc(:,1,:)*p_patch%cells%area(:,:) )
    ENDWHERE

    ! Fix undershoots - ONLY for the one-ice-class case
    ! Quick fix, should be reformulated to occur at the advection stage
    WHERE ( ( p_ice%conc(:,1,:) < 0._wp ) .OR. ( p_ice%hi(:,1,:) < 0._wp ) )
      p_ice%conc(:,1,:) = 0._wp
      p_ice%hi(:,1,:)   = 0._wp
      p_ice%hs(:,1,:)   = 0._wp
    ENDWHERE

    p_ice%concSum                           = SUM(p_ice%conc, 2)

  END SUBROUTINE ice_clean_up_dyn

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! !  ice_clean_up_thd: Fix undershoots and beutify output after the thermodynamic growth/melt changes.
  !! !  Calculates the final freeboard
  !! @par Revision History
  !! Initial release by Einar Olason, MPI-M (2013-10).
  !! Modified by Vladimir Lapin,      MPI-M (2015-07).
  !!
  SUBROUTINE ice_clean_up_thd( p_patch_3D, p_ice, atmos_fluxes, p_os )
    TYPE(t_patch_3D),TARGET,   INTENT(IN)    :: p_patch_3D
    TYPE(t_sea_ice),           INTENT(INOUT) :: p_ice
    TYPE(t_atmos_fluxes),      INTENT(INOUT) :: atmos_fluxes
    TYPE(t_hydro_ocean_state), INTENT(IN)    :: p_os

    ! Local variables
    ! ranges
    TYPE(t_subset_range), POINTER                                            :: all_cells
    ! pathc
    TYPE(t_patch),POINTER                                                    :: p_patch
    TYPE(t_patch_vert),POINTER                                               :: p_patch_vert
    ! counters
    INTEGER                                                                  :: k, jb, jc, i_startidx_c, i_endidx_c
    ! Sea surface salinity
    REAL(wp), DIMENSION (nproma, p_patch_3d%p_patch_2D(1)%alloc_cell_blocks) :: sss
    REAL(wp)                                                                 :: z_smax

    ! subset range pointer
    p_patch      => p_patch_3D%p_patch_2D(1)
    p_patch_vert => p_patch_3D%p_patch_1D(1)
    all_cells    => p_patch%cells%all
    ! Sea surface salinity
    sss(:,:)  =  p_os%p_prog(nold(1))%tracer(:,1,:,2)

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c

        ! For prettier output we set speed to zero where there's no ice
        ! This does not affect the dynamics, since the EVP routine doesn't modify ice velocities where
        ! concentration is less than 0.01
        IF ( p_ice%hi(jc,1,jb) <= 0._wp ) THEN
          p_ice%u(jc,jb) = 0._wp
          p_ice%v(jc,jb) = 0._wp
        ENDIF

        DO k = 1, p_ice%kice

        ! Fix under shoots and remove ice where there's almost none left
        ! There should be no undershoots if the advection schene is monotonic and sign preserving
        ! Therefore we do not need to correct undershoots so far (slo 2015-06)
        ! IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary &
        !   &   .AND. ( p_ice%vol(jc,k,jb) <= 0._wp .OR. p_ice%conc(jc,k,jb) <= 1e-4_wp ) ) THEN
        !   ! Tracer flux due to removal
        !   atmos_fluxes%FrshFlux_TotalIce (jc,jb) = atmos_fluxes%FrshFlux_TotalIce (jc,jb)                      &
        !     & + (1._wp-sice/sss(jc,jb))*p_ice%hi(jc,k,jb)*p_ice%conc(jc,k,jb)*rhoi/(rho_ref*dtime)  & ! Ice
        !     & + p_ice%hs(jc,k,jb)*p_ice%conc(jc,k,jb)*rhos/(rho_ref*dtime)                           ! Snow
        !   ! Heat flux due to removal
        !   atmos_fluxes%HeatFlux_Total(jc,jb) = atmos_fluxes%HeatFlux_Total(jc,jb)   &
        !     & + p_ice%hi(jc,k,jb)*p_ice%conc(jc,k,jb)*alf*rhoi/dtime          & ! Ice
        !     & + p_ice%hs(jc,k,jb)*p_ice%conc(jc,k,jb)*alf*rhos/dtime            ! Snow
        !   p_ice%conc(jc,k,jb) = 0._wp
        !   p_ice%hi  (jc,k,jb) = 0._wp
        !   p_ice%vol (jc,k,jb) = 0._wp
        !   p_ice%hs  (jc,k,jb) = 0._wp
        !   p_ice%vols(jc,k,jb) = 0._wp
        ! ENDIF

        ! limit sea ice thickness to seaice_limit of surface layer depth, without elevation
          IF (limit_seaice) THEN
            z_smax = seaice_limit*p_patch_3D%p_patch_1D(1)%del_zlev_m(1)
                  IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary  .AND.  p_ice%hi(jc,k,jb) > z_smax ) THEN
                    ! Tracer flux due to removal
                    atmos_fluxes%FrshFlux_TotalIce (jc,jb) = atmos_fluxes%FrshFlux_TotalIce (jc,jb)                      &
                      & + (1._wp-sice/sss(jc,jb))*(p_ice%hi(jc,k,jb)-z_smax)*p_ice%conc(jc,k,jb)*rhoi/(rho_ref*dtime)  ! Ice
                    ! Heat flux due to removal
                    !  #slo# 2015-02: - this heat did not come from the ocean, but from atmosphere, heating ocean is wrong
                    !                 - check if conc must enter here as well, check energy for coupling
                    atmos_fluxes%HeatFlux_Total(jc,jb) = atmos_fluxes%HeatFlux_Total(jc,jb)   &
                      & + (p_ice%hi(jc,k,jb)-z_smax)*p_ice%conc(jc,k,jb)*alf*rhoi/dtime           ! Ice
                    p_ice%hi  (jc,k,jb) = z_smax
                    p_ice%vol (jc,k,jb) = p_ice%hi(jc,k,jb)*p_ice%conc(jc,k,jb)*p_patch%cells%area(jc,jb)
                  ENDIF
          END IF
          p_ice%draft(jc,k,jb)   = (rhos * p_ice%hs(jc,k,jb) + rhoi * p_ice%hi(jc,k,jb)) / rho_ref
        ENDDO
        p_ice%draftave (jc,jb) = sum(p_ice%draft(jc,:,jb) * p_ice%conc(jc,:,jb))
        p_ice%zUnderIce(jc,jb) = p_patch_vert%prism_thick_flat_sfc_c(jc,1,jb) + p_os%p_prog(nold(1))%h(jc,jb) &
          &                      - p_ice%draftave(jc,jb)  + p_ice%totalsnowfall(jc,jb)
        !  #slo# 2015-01: totalsnowfall is needed for correct salt update (in surface module)
        !                 since draft was increased by snowfall but water below ice is not effected by snowfall
        !                 snow to ice conversion does not effect draft
      ENDDO
    ENDDO

    p_ice%concSum                           = SUM(p_ice%conc, 2)
    atmos_fluxes%cellThicknessUnderIce(:,:) = p_ice%zUnderIce(:,:)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('iceClUp: hi aft. limiter'     ,p_ice%hi       ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('iceClUp: hs aft. limiter'     ,p_ice%hs       ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('iceClUp: Conc. aft. limiter'  ,p_ice%conc     ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('iceClUp: ConcSum aft. limit ' ,p_ice%concSum  ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceClUp: HeatTotal a. limit',  atmos_fluxes%HeatFlux_Total   ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceClUp: TotalIce  a. limit',  atmos_fluxes%FrshFlux_TotalIce,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceClUp: draft    '           ,p_ice%draft    ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceClUp: draftave '           ,p_ice%draftave ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceClUp: zUnderIce'           ,p_ice%zUnderIce,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('iceClUp: h-old'               ,p_os%p_prog(nold(1))%h,str_module,4, in_subset=p_patch%cells%owned)

  END SUBROUTINE ice_clean_up_thd


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! !  get_atmos_fluxes: Sets the atmospheric fluxes for the update of the ice
  ! !                 temperature
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE get_atmos_fluxes (p_patch, p_os,p_as,ice, atmos_fluxes)
    TYPE(t_patch),            INTENT(IN)    :: p_patch
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os
    TYPE(t_atmos_for_ocean),  INTENT(IN)    :: p_as
    TYPE (t_sea_ice),         INTENT(INOUT) :: ice
    TYPE (t_atmos_fluxes),    INTENT(INOUT) :: atmos_fluxes

!#ifdef coupled
    !atmos_fluxes% SWin   =
    !atmos_fluxes% LWin   =
    !atmos_fluxes% sens   =
    !atmos_fluxes% lat    =
    !atmos_fluxes% dsensdT =
    !atmos_fluxes% dlatdT  =
    !atmos_fluxes% dLWdT   =
!#elif defined CORE
    !CALL budget_core   (ice, atmos_fluxes)
!#else
    CALL calc_bulk_flux_oce(p_patch, p_as, p_os, atmos_fluxes)
    CALL calc_bulk_flux_ice(p_patch, p_as, ice , atmos_fluxes)
!#endif

  END SUBROUTINE get_atmos_fluxes
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! !   sum_fluxes: adds atmospheric fluxes for ocean time stepping. Necessary for
  !!      diagnosis, not for the ice model itself.
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE sum_fluxes        (atmos_fluxes, atmos_fluxesAve)
    TYPE (t_atmos_fluxes), INTENT (IN)    :: atmos_fluxes
    TYPE (t_atmos_fluxes), INTENT (INOUT) :: atmos_fluxesAve

    atmos_fluxesAve % sens   (:,:,:) = atmos_fluxesAve % sens   (:,:,:) + atmos_fluxes % sens   (:,:,:)
    atmos_fluxesAve % sensw  (:,:)   = atmos_fluxesAve % sensw  (:,:)   + atmos_fluxes % sensw  (:,:)
    atmos_fluxesAve % lat    (:,:,:) = atmos_fluxesAve % lat    (:,:,:) + atmos_fluxes % lat    (:,:,:)
    atmos_fluxesAve % latw   (:,:)   = atmos_fluxesAve % latw   (:,:)   + atmos_fluxes % latw   (:,:)
    atmos_fluxesAve % LWout  (:,:,:) = atmos_fluxesAve % LWout  (:,:,:) + atmos_fluxes % LWout  (:,:,:)
    atmos_fluxesAve % LWoutw (:,:)   = atmos_fluxesAve % LWoutw (:,:)   + atmos_fluxes % LWoutw (:,:)
    atmos_fluxesAve % LWnet  (:,:,:) = atmos_fluxesAve % LWnet  (:,:,:) + atmos_fluxes % LWnet  (:,:,:)
    atmos_fluxesAve % LWnetw (:,:)   = atmos_fluxesAve % LWnetw (:,:)   + atmos_fluxes % LWnetw (:,:)
    atmos_fluxesAve % SWnet  (:,:,:) = atmos_fluxesAve % SWnet  (:,:,:) + atmos_fluxes % SWnet  (:,:,:)
    atmos_fluxesAve % SWnetw (:,:)   = atmos_fluxesAve % SWnetw (:,:)   + atmos_fluxes % SWnetw (:,:)
    atmos_fluxesAve % LWin   (:,:)   = atmos_fluxesAve % LWin   (:,:)   + atmos_fluxes % LWin   (:,:)
    atmos_fluxesAve % rprecw (:,:)   = atmos_fluxesAve % rprecw (:,:)   + atmos_fluxes % rprecw (:,:)
    atmos_fluxesAve % rpreci (:,:)   = atmos_fluxesAve % rpreci (:,:)   + atmos_fluxes % rpreci (:,:)
    atmos_fluxesAve % counter        = atmos_fluxesAve % counter + 1

  END SUBROUTINE sum_fluxes
  !-------------------------------------------------------------------------------
  !
  !
  !>
  !! ! ave_fluxes: calculates the average of the atmospheric fluxes for ocean time
  !!   sum_fluxes: adds atmospheric fluxes for ocean time stepping. Necessary for
  !!   diagnosis, not for the ice model itself.
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE ave_fluxes (ice, atmos_fluxesAve)
    TYPE (t_sea_ice),      INTENT (INOUT) :: ice
    TYPE (t_atmos_fluxes), INTENT (INOUT) :: atmos_fluxesAve
    !
    !Local variables
    REAL(wp) :: ctr

    !-------------------------------------------------------------------------------

    ctr = REAL(atmos_fluxesAve% counter,wp)
    atmos_fluxesAve% sens   (:,:,:) = atmos_fluxesAve% sens  (:,:,:)  / ctr
    atmos_fluxesAve% sensw  (:,:)   = atmos_fluxesAve% sensw (:,:)    / ctr
    atmos_fluxesAve% lat    (:,:,:) = atmos_fluxesAve% lat   (:,:,:)  / ctr
    atmos_fluxesAve% latw   (:,:)   = atmos_fluxesAve% latw  (:,:)    / ctr
    atmos_fluxesAve% LWout  (:,:,:) = atmos_fluxesAve% LWout (:,:,:)  / ctr
    atmos_fluxesAve% LWoutw (:,:)   = atmos_fluxesAve% LWoutw(:,:)    / ctr
    atmos_fluxesAve% LWnet  (:,:,:) = atmos_fluxesAve% LWnet (:,:,:)  / ctr
    atmos_fluxesAve% LWnetw (:,:)   = atmos_fluxesAve% LWnetw(:,:)    / ctr
    atmos_fluxesAve% SWnet  (:,:,:) = atmos_fluxesAve% SWnet (:,:,:)  / ctr
    atmos_fluxesAve% SWnetw (:,:)   = atmos_fluxesAve% SWnetw(:,:)    / ctr
    atmos_fluxesAve% LWin   (:,:)   = atmos_fluxesAve% LWin  (:,:)    / ctr
    atmos_fluxesAve% rprecw (:,:)   = atmos_fluxesAve% rprecw(:,:)    / ctr
    atmos_fluxesAve% rpreci (:,:)   = atmos_fluxesAve% rpreci(:,:)    / ctr
    ice    % Qbot   (:,:,:) = ice    % Qbot  (:,:,:)  / ctr
    ice    % Qtop   (:,:,:) = ice    % Qtop  (:,:,:)  / ctr

  END SUBROUTINE ave_fluxes
  !-------------------------------------------------------------------------------
  !
  !
  !>
  !! ! ice_zero: set the avereged fluxes to zero
  !! @par Revision History
  !! Initial release by Einar Olason, MPI-M (2011-09). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE ice_zero (ice)
    TYPE (t_sea_ice),      INTENT (INOUT) :: ice
    !TYPE (t_atmos_fluxes), INTENT (INOUT) :: atmos_fluxes
    !TYPE (t_atmos_fluxes), INTENT (INOUT) :: atmos_fluxesAve

    !atmos_fluxes    % sens        (:,:,:) = 0._wp
    !atmos_fluxes    % sensw       (:,:)   = 0._wp
    !atmos_fluxes    % lat         (:,:,:) = 0._wp
    !atmos_fluxes    % latw        (:,:)   = 0._wp
    !atmos_fluxes    % LWout       (:,:,:) = 0._wp
    !atmos_fluxes    % LWoutw      (:,:)   = 0._wp
    !atmos_fluxes    % LWnet       (:,:,:) = 0._wp
    !atmos_fluxes    % LWnetw      (:,:)   = 0._wp
    !atmos_fluxes    % SWin        (:,:)   = 0._wp
    !atmos_fluxes    % LWin        (:,:)   = 0._wp
    !atmos_fluxes    % rprecw      (:,:)   = 0._wp
    !atmos_fluxes    % rpreci      (:,:)   = 0._wp

!    atmos_fluxesAve % sens        (:,:,:) = 0._wp
!    atmos_fluxesAve % sensw       (:,:)   = 0._wp
!    atmos_fluxesAve % lat         (:,:,:) = 0._wp
!    atmos_fluxesAve % latw        (:,:)   = 0._wp
!    atmos_fluxesAve % LWout       (:,:,:) = 0._wp
!    atmos_fluxesAve % LWoutw      (:,:)   = 0._wp
!    atmos_fluxesAve % LWnet       (:,:,:) = 0._wp
!    atmos_fluxesAve % LWnetw      (:,:)   = 0._wp
!    atmos_fluxesAve % SWin        (:,:)   = 0._wp
!    atmos_fluxesAve % LWin        (:,:)   = 0._wp
!    atmos_fluxesAve % rprecw      (:,:)   = 0._wp
!    atmos_fluxesAve % rpreci      (:,:)   = 0._wp
!    atmos_fluxesAve % counter             = 0

!    ice     % Qbot        (:,:,:) = 0._wp
!    ice     % Qtop        (:,:,:) = 0._wp
!ICON_OMP_PARALLEL
!ICON_OMP_WORKSHARE
    ice     % surfmelt    (:,:,:) = 0._wp
    ice     % surfmeltT   (:,:,:) = 0._wp
    ice     % evapwi      (:,:,:) = 0._wp
    ice     % hiold       (:,:,:) = 0._wp
    ice     % hsold       (:,:,:) = 0._wp
    ice     % snow_to_ice (:,:,:) = 0._wp
    ice     % heatOceI    (:,:,:) = 0._wp
!ICON_OMP_END_WORKSHARE
!ICON_OMP_END_PARALLEL

  END SUBROUTINE ice_zero

  !-------------------------------------------------------------------------------
  !
  !
  !>
  !! ! ice_albedo: set ice albedo
  !-------------------------------------------------------------------------------
  !
  !
  !>
  !! ! ice_albedo: set ice albedo
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE set_ice_albedo(i_startidx_c, i_endidx_c, nbdim, kice, Tsurf, hi, hs, &
      & albvisdir, albvisdif, albnirdir, albnirdif)
    INTEGER, INTENT(IN)  :: i_startidx_c, i_endidx_c, nbdim, kice
    REAL(wp),INTENT(IN)  :: Tsurf(nbdim,kice)
    REAL(wp),INTENT(IN)  :: hi   (nbdim,kice)
    REAL(wp),INTENT(IN)  :: hs   (nbdim,kice)
    REAL(wp),INTENT(OUT) :: albvisdir  (nbdim,kice)
    REAL(wp),INTENT(OUT) :: albvisdif  (nbdim,kice)
    REAL(wp),INTENT(OUT) :: albnirdir  (nbdim,kice)
    REAL(wp),INTENT(OUT) :: albnirdif  (nbdim,kice)


    !Local variables
    REAL(wp), PARAMETER :: albtrans   = 0.5_wp
    REAL(wp)            :: albflag, frac_snow
    INTEGER             :: jc,k
    !-------------------------------------------------------------------------------

    SELECT CASE (i_ice_albedo)
    CASE (1)
      ! This is Uwe's albedo expression from the old budget function
      DO k=1,kice
        DO jc = i_startidx_c,i_endidx_c

          albflag =  1.0_wp/ ( 1.0_wp+albtrans * (Tsurf(jc,k))**2 )

          IF ( hi(jc,k) > 0._wp ) THEN
            IF ( hs(jc,k) > 1.e-2_wp ) THEN
              albvisdir(jc,k) =  albflag * albsm + (1.0_wp-albflag) * albs
            ELSE
              albvisdir(jc,k) =  albflag * albim + (1.0_wp-albflag) * albi
            ENDIF
          ELSE
            albvisdir(jc,k) = 0._wp
          ENDIF

        ENDDO
      ENDDO

      ! all albedos are the same
      albvisdif = albvisdir
      albnirdir = albvisdir
      albnirdif = albvisdir

    CASE (2)
      ! This is the CCSM 3 albedo scheme
!PREVENT_INCONSISTENT_IFORT_FMA
      DO k=1,kice
        DO jc = i_startidx_c,i_endidx_c
          frac_snow = hs(jc,k)/( hs(jc,k)+0.02_wp )
          IF ( Tsurf(jc,k) > -1._wp ) THEN
            albvisdir(jc,k) = frac_snow*( alb_sno_vis - 0.100_wp*(Tsurf(jc,k)+1._wp) ) &
              &     + (1._wp-frac_snow)*( alb_ice_vis - 0.075_wp*(Tsurf(jc,k)+1._wp) )
            albnirdir(jc,k) = frac_snow*( alb_sno_nir - 0.150_wp*(Tsurf(jc,k)+1._wp) ) &
              &     + (1._wp-frac_snow)*( alb_ice_nir - 0.075_wp*(Tsurf(jc,k)+1._wp) )
          ELSE
            albvisdir(jc,k) = frac_snow*alb_sno_vis + (1._wp-frac_snow)*alb_ice_vis
            albnirdir(jc,k) = frac_snow*alb_sno_nir + (1._wp-frac_snow)*alb_ice_nir
          ENDIF
        ENDDO
      ENDDO

      ! diffuse and direct albedos are the same
      albvisdif = albvisdir
      albnirdif = albnirdir

    END SELECT

  END SUBROUTINE set_ice_albedo



  !-------------------------------------------------------------------------------
  !
  !
  !>
  !! ! upper_ocean_TS: Adjusts the temperature and salinity of the upper ocean grid
  !!                 cell according to atmospheric heat and fresh-water fluxes,
  !!                 surface melting, ice growth, etc. The upper ocean temperature
  !!                 is also changed in subroutine ice_conc_change and at the
  !!                beginning of subroutine ice_growth
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE upper_ocean_TS(p_patch, p_patch_vert, ice, p_os, atmos_fluxes)
    TYPE(t_patch),TARGET,             INTENT(IN)    :: p_patch
    TYPE(t_patch_vert),        INTENT(IN)    :: p_patch_vert
    TYPE(t_sea_ice),           INTENT(INOUT) :: ice
    TYPE(t_hydro_ocean_state), INTENT(IN)    :: p_os
    TYPE(t_atmos_fluxes),      INTENT(INOUT) :: atmos_fluxes

    !Local Variables
    REAL(wp), DIMENSION (nproma, p_patch%alloc_cell_blocks) ::   &
      & zUnderIceOld,  &! water in upper ocean grid cell below ice (prev. time)   [m]
      & heatOceI,      &! heat flux into ocean through formerly ice covered areas [W/m^2]
      & heatOceW,      &! heat flux into ocean through open water areas           [W/m^2]
      & Delhice,       &! average change in ice thickness within a grid cell      [m]
      & Delhsnow,      &! average change in snow thickness within a grid cell     [m]
      & snowiceave,    &! average snow to ice conversion within a grid cell       [m]
      & Tfw,           &! sea surface freezing temperature                        [C]
      & csst,          &! sea surface temperature - approx. after cooling (input) [C]
      & sss,           &! sea surface salinity (input only)                       [psu]
      & preci,         &! solid precipitation rate                                [m/s]
      & precw           ! liquid precipitation rate                               [m/s]
      !& evap,          &! evaporated water                                       [psu]

    REAL(wp), DIMENSION (nproma, p_patch%alloc_cell_blocks) :: h1,h2,h3,snowmelted,xheat,xnewice

    TYPE(t_subset_range), POINTER :: subset
    INTEGER :: block, cell, cellStart,cellEnd

    zUnderIceOld(:,:) = 0.0_wp
    heatOceI    (:,:) = 0.0_wp
    heatOceW    (:,:) = 0.0_wp
    Delhice     (:,:) = 0.0_wp
    Delhsnow    (:,:) = 0.0_wp
    snowiceave  (:,:) = 0.0_wp
    Tfw         (:,:) = 0.0_wp
    csst        (:,:) = 0.0_wp
    sss         (:,:) = 0.0_wp
    preci       (:,:) = 0.0_wp
    precw       (:,:) = 0.0_wp
    h1          (:,:) = 0.0_wp
    h2          (:,:) = 0.0_wp
    h3          (:,:) = 0.0_wp
    snowmelted  (:,:) = 0.0_wp
    xheat       (:,:) = 0.0_wp
    xnewice     (:,:) = 0.0_wp

    !-------------------------------------------------------------------------------
    CALL dbg_print('UpperOcTS beg: hi',   ice%hi,                              str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOcTS beg: hs',   ice%hs,                              str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOcTS beg: conc', ice%conc,                            str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOcTS beg: zUnderIce', ice%zUnderIce,                  str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOcTS beg: h',    p_os%p_prog(nold(1))%h,              str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOcTS beg: SST',  p_os%p_prog(nold(1))%tracer(:,1,:,1),str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOcTS beg: SSS',  p_os%p_prog(nold(1))%tracer(:,1,:,2),str_module, 4, in_subset=p_patch%cells%owned)
    !-------------------------------------------------------------------------------

    !TODOram: openmp
    subset => p_patch%cells%all
    DO block = subset%start_block, subset%end_block
      CALL get_index_range(subset, block, cellStart, cellEnd)
      DO cell = cellStart, cellEnd
        IF (subset%vertical_levels(cell,block) < 1) CYCLE
        ! #eoo# What is swsum?
        ! swsum = 0.0_wp
        !sao_top =>p_os%p_prog(nold(1))%tracer(:,1,:,2)

        ! Ocean points only
        ! Calculate change in water level 'zo' from liquid and solid precipitation and
        ! evaporation
        sss             (cell,block)   = p_os%p_prog(nold(1))%tracer(cell,1,block,2)
        precw           (cell,block)   = atmos_fluxes% rprecw (cell,block)
        preci           (cell,block)   = atmos_fluxes% rpreci (cell,block)
        !evap            (cell,block)   = (atmos_fluxes% latw(cell,block)/ alv * dtime * &
        !  &                       sum(ice%conc(cell,1,block), 2) +          &
        !  &                       sum(ice%evapwi(cell,1,block) * ice% conc(cell,1,block), 2)) /rho_ref

        ! Calculate the sea surface freezing temperature                        [C]
        IF ( no_tracer < 2 .OR. use_constant_tfreez ) THEN
          Tfw(cell,block) = Tf
        ELSE
          Tfw(cell,block) = -mu*sss(cell,block)
        ENDIF

        ! TODO: No temperature change due to precip yet. Should not be done here

        ! Calculate average draft and thickness of water underneath ice in upper ocean grid box
        !  - changes in hi and hs from sbrt ice_growth are included here
        zUnderIceOld    (cell,block)   = ice%zUnderIce(cell,block)   !  not needed anymore
      ! ice%draft       (cell,:,block) = (rhos * ice%hs(cell,:,block) + rhoi * ice%hi(cell,:,block)) / rho_ref
      ! ice%draftave    (cell,block)   = sum(ice%draft(cell,:,block) * ice%conc(cell,:,block))
        ice%totalsnowfall(cell,block) =  atmos_fluxes%rpreci(cell,block)*dtime*ice%conc(cell,1,block)
      ! ice%zUnderIce   (cell,block)   = p_patch_vert%prism_thick_flat_sfc_c(cell,1,block) &
      !   &                            + p_os%p_prog(nold(1))%h(cell,block)                &
      !   &                            - ice%draftave(cell,block)
        !  #slo# 2015-01: bugfix - totalsnowfall was added to zunderice - inconsistent to zui=h-draftave
        ! &                            + ice%totalsnowfall(cell,block)

        ! Calculate average change in ice thickness and the snow-to-ice conversion with respect to whole grid area
        ! #slo# 2015-01: could be expressed as water equivalent here
        Delhice   (cell,block) = SUM( ( ice%hi(cell,:,block) - ice%hiold(cell,:,block) )*ice%conc(cell,:,block))
        Delhsnow  (cell,block) = SUM( ( ice%hs(cell,:,block) - ice%hsold(cell,:,block) )*ice%conc(cell,:,block))

        ! snowiceave in averaged ice thickness (m)
        !snowiceave(cell,block) = SUM( ice%snow_to_ice(cell,:,block)*ice% conc(cell,:,block))
        ! snowiceave now in water equivalent, snow_to_ice has snow density
        snowiceave(cell,block) = SUM( ice%snow_to_ice(cell,:,block)*ice% conc(cell,:,block))*rhos/rho_ref

        ! Adjust change in snow and ice thickness for snow-ice formation, which is dealt with separately in the fresh-water balance
        !!Delhsnow  (cell,block) = Delhsnow(cell, block) + snowiceave (cell, block)  

        ! Calculate heat input through formerly ice covered and through open water areas
        ! #slo# 2014-12: bugfix for qio3: heatOceI is for the whole grid area
        !                heat needs scaling with concentration to be valid for whole area in all cases
      ! IF (3 == i_Qio_type) THEN
      !   ! If the whole heat from the ice is used to melt the ice, the scaling
      !   ! with the ice concentration has to be skipped here
      !   heatOceI(cell,block)   = sum(ice% heatOceI(cell,:,block))
      ! ELSE
      !   heatOceI(cell,block)   = sum(ice% heatOceI(cell,:,block) * ice% conc(cell,:,block))
      ! END IF
        ! #slo# 2014-12: rescaling of melting/growing energy with concentration is done before
        !                in ice_slow to use energy for SST change
        heatOceI(cell,block)   = sum(ice% heatOceI(cell,:,block))
        heatOceW(cell,block) = ( atmos_fluxes%SWnetw(cell,block) &
          &                    + atmos_fluxes%LWnetw(cell,block) &
          &                    + atmos_fluxes%sensw(cell,block)  &
          &                    + atmos_fluxes%latw(cell,block) ) &
                               *(1.0_wp-sum(ice%conc(cell,:,block)))

        ! Calculate possible super-cooling of the surface layer by heatOceW
        !  #slo# 2014-12-11: heatOceW is part of heat flux over open water, supercooled sst is calculated over whole grid area
        xheat(cell,block)= heatocew(cell,block)
        csst(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,1) &
          &                + dtime*heatOceW(cell,block)/( clw*rho_ref*ice%zUnderIce(cell,block) )

        ! #slo# 2015-01: control calculation of newice directly from heat flux heatOceW
        !  - check effect of supercooling is required
        !  - newice is in ice thickness
      ! IF (p_os%p_prog(nold(1))%tracer(cell,1,block,1)< Tfw(cell,block) .AND. v_base%lsm_c(cell,1,block) <= sea_boundary ) THEN
      !   xnewice(cell,block) = - heatOceW(cell,block)*dtime/(alf*rhoi)
      !   heatOceW(cell,block) = 0.0_wp
      ! ELSE
      !   xnewice(cell,block) = 0.0_wp
      ! ENDIF

        ! #slo# 2014-11: wrong comment, energy not yet added to SST, but must be added to heatOceW or heatOceI
        !                part of heatOceW used for forming newice is added here
        ! Add energy for new-ice formation due to supercooled ocean to ocean temperature, form new ice
        IF ( csst(cell,block) < Tfw(cell,block) .AND. v_base%lsm_c(cell,1,block) <= sea_boundary ) THEN
          ! New ice forming over open water due to super cooling
          ! Fixed 2. April (2014) - newice is now the volume of ice formed over open water
          ! #slo# 2014-12: bugfix: heatOceW is already reduced to (1-conc)*heat,
          !                        therefore newice is value for whole grid area, not multiplied by 1-conc
          !                        results in larger sea ice growth
          !  Attention: in ice_concentration_change as well, newice must not be multiplied by 1-conc
          !ice%newice(cell,block) = (1._wp-ice%concSum(cell,block))*( Tfw(cell,block) - csst(cell,block) ) &
          ice%newice(cell,block) = ( Tfw(cell,block) - csst(cell,block) ) &
             &                     * ice%zUnderIce(cell,block)*clw*rho_ref/( alf*rhoi )
          ! Flux required to warm the ocean to the freezing point
          ! #slo# 2015-01: bugfix for supercooled ocean: whole grid area is cooled below freezing point,
          !                no multiplication by 1-conc is required to calculate the part of heatOceW due to supercooling
          !                but check: doubled supercooling, since it is already taken into account in oce_ice_heatflx?
          heatOceW(cell,block)   = ( Tfw(cell,block) - p_os%p_prog(nold(1))%tracer(cell,1,block,1) )     &
            &                      *ice%zUnderIce(cell,block)*clw*rho_ref/dtime
         !  &                      *ice%zUnderIce(cell,block)*(1.0_wp-ice%concSum(cell,block))*clw*rho_ref/dtime

          ! #slo# 2014-11: bugfix: value of newice was kept from last timestep!
        ELSE
          ice%newice(cell,block) = 0.0_wp
        ENDIF

        ! test
      ! ice%newice(cell,block) = xnewice(cell,block)

        ! Diagnosis: collect the 4 parts of heat fluxes into the atmos_fluxes variables - no flux under ice:
        atmos_fluxes%HeatFlux_ShortWave(cell,block) = atmos_fluxes%SWnetw(cell,block)*(1.0_wp-sum(ice%conc(cell,:,block)))
        atmos_fluxes%HeatFlux_LongWave (cell,block) = atmos_fluxes%LWnetw(cell,block)*(1.0_wp-sum(ice%conc(cell,:,block)))
        atmos_fluxes%HeatFlux_Sensible (cell,block) = atmos_fluxes%sensw (cell,block)*(1.0_wp-sum(ice%conc(cell,:,block)))
        atmos_fluxes%HeatFlux_Latent   (cell,block) = atmos_fluxes%latw  (cell,block)*(1.0_wp-sum(ice%conc(cell,:,block)))

        ! #slo# 2013-06
        ! Change of upper ocean temperature according to heat fluxes is done in vertical diffusion equation
        !  - HeatFlux_Total is calculated here, provided to tracer eq. using topBoundCond_Temp_vdiff, calculated in update_sfcflx
        !  - topBoundCond_Temp_vdiff is 
        !p_os%p_prog(nold(1))%tracer(cell,1,block,1) = p_os%p_prog(nold(1))%tracer(cell,1,block,1)&
        !  &                                    + dtime*(heatOceI + heatOceW) /               &
        !  &                                    (clw*rho_ref * ice%zUnderIce)
        ! TODO: should we also divide with ice%zUnderIce / ( v_base%del_zlev_m(1) +  p_os%p_prog(nold(1))%h(cell,block) ) ?
        !atmos_fluxes%topBoundCond_Temp_vdiff = (heatOceI + heatOceW) / (clw*rho_ref)
        atmos_fluxes%HeatFlux_Total(cell,block) = heatOceI(cell,block) + heatOceW(cell,block)

        ! TODO:
        ! Temperature change of upper ocean grid cell due  to melt-water inflow and
        ! precipitation
        !p_os%p_prog(nold(1))%tracer(cell,1,block,1) = (p_os%p_prog(nold(1))%tracer(cell,1,block,1) &
        !  &                      *zUnderIceOld                                       &
        !  &                      + precw*p_as%tafo + preci*0.0_wp + &                             !!!!!!!!!Dirk: times 0.0 ????
        !  &                        sum(ice%surfmeltT(cell,1,block) * ice%surfmelt * ice%conc(cell,1,block),2)) / &
        !  &                        (zUnderIceOld + sum(ice%surfmelt*ice%conc(cell,1,block),2) +    &
        !  &                        precw + preci)
        !
        ! Change salinity of upper ocean grid box from ice growth/melt, snowice
        ! formation and precipitation
        !p_os%p_prog(nold(1))%tracer(cell,1,block,2) = p_os%p_prog(nold(1))%tracer(cell,1,block,2)  &
        !  &                                    + (Delhice(cell,block)*rhoi - snowiceave(cell,block)*rhos)/rho_ref *  &
        !  &                                    MIN(Sice, sao_top(cell,block)) / ice%zUnderIce(cell,block)

        ! #slo# 2013-06
        ! Change in salinity is calculated according to resulting freshwater flux due to sea ice change:
        !  - fw_ice_impl is flux in m/s >0 for Delhice<0, i.e. positive input of water = decrease of sea ice depth
        !atmos_fluxes%forc_fwsice(cell,block) = -Delhice(cell,block)*rhoi - snowiceave(cell,block)*rhos)/(rho_ref*dtime)

        ! Volmue flux
        ! Fixed 27. March
        ! unit: m/s
        ! TODO: VolumeFlux_UnderIce only from rain
        atmos_fluxes%FrshFlux_VolumeIce(cell,block) = precw(cell,block)*ice%concSum(cell,block)   !& ! Rain goes through

        ! Tracer flux
        ! Fixed 27. March
        ! -->> snow growth (Delhsnow > 0) does NOT change tracers, but snow melt does
        ! -->> rain is only included in the volumeice, totalice now only contains
        !      fluxes that do not change the total fresh-water volume in the grid cell
        !  #slo# 2015-01: better reformulate with ice-change, snowmelt, snoiceave - assign variables in sbrt icegrowth!
        IF (v_base%lsm_c(cell,1,block) <= sea_boundary ) THEN
          atmos_fluxes%FrshFlux_TotalIce(cell,block) = &
        &        - (1._wp-sice/sss(cell,block))*Delhice(cell,block)*rhoi/(rho_ref*dtime)&   ! Ice melt and growth
        &            - MERGE(&
        &         (Delhsnow(cell,block)-ice%totalsnowfall(cell,block)/rhos*rho_ref)*rhos/(rho_ref*dtime), &
        &         0.0_wp, &
        &         Delhsnow(cell,block)-ice%totalsnowfall(cell,block)/rhos*rho_ref < 0.0_wp) &  ! snow melt ONLY
        &        - (1._wp-sice/sss(cell,block))*ice%newice(cell,block)*rhoi/(rho_ref*dtime)  ! new ice formation over open water

          ! divide three processes:
          ! ice melt and growth including snow_to_ice conversion proportional to salt difference of water and ice
          h1(cell,block)= - (1._wp-sice/sss(cell,block))*Delhice(cell,block)*rhoi/(rho_ref*dtime)
          ! snow melt only, without snow_to_ice conversion, in m water equivalent
          snowmelted(cell,block) = (Delhsnow(cell,block)-ice%totalsnowfall(cell,block)/rhos*rho_ref)*rhos/rho_ref &
            &                      + snowiceave(cell,block)
          ! snow melt only, without snow_to_ice conversion, in m water equivalent
          IF (Delhsnow(cell,block)-ice%totalsnowfall(cell,block)/rhos*rho_ref < 0.0_wp) THEN
            h2(cell,block) = -(Delhsnow(cell,block)-ice%totalsnowfall(cell,block)/rhos*rho_ref)*rhos/(rho_ref*dtime)
          !IF (snowmelted(cell,block) < 0.0_wp) THEN
          !  h2(cell,block) = -snowmelted(cell,block)/dtime
          ELSE
            h2(cell,block) = 0.0_wp
          ENDIF
          ! new ice formation
          h3(cell,block) =-(1._wp-sice/sss(cell,block))*ice%newice(cell,block)*rhoi/(rho_ref*dtime)
          
          ! test the new total ice flux:
          !atmos_fluxes%FrshFlux_TotalIce(cell,block) = h1(cell,block) + h2(cell,block) + h3(cell,block)

        ENDIF

        !heatabs         (cell,block)   = swsum * atmos_fluxes% SWin(cell,block) * (1 - ice%concsum)

        ! set to zero on land points
        IF (v_base%lsm_c(cell,1,block) > sea_boundary ) THEN
          atmos_fluxes%HeatFlux_Total    (cell,block) = 0.0_wp
          atmos_fluxes%HeatFlux_ShortWave(cell,block) = 0.0_wp
          atmos_fluxes%HeatFlux_LongWave (cell,block) = 0.0_wp
          atmos_fluxes%HeatFlux_Sensible (cell,block) = 0.0_wp
          atmos_fluxes%HeatFlux_Latent   (cell,block) = 0.0_wp
        END IF

      END DO
    END DO

    CALL dbg_print('UpperOceTS: hi       ', ice%hi       ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: hiold    ', ice%hiold    ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: hs       ', ice%hs       ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: hsold    ', ice%hsold    ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: Delhice  ', Delhice      ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: Delhsnow ', Delhsnow     ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: draft    ', ice%draft    ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: draftave ', ice%draftave ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: zUndIceOld',zUnderIceOld, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: zUnderIce', ice%zUnderIce,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: csst-cooled', csst       ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: SSS      ', sss          ,str_module, 5, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: newice   ', ice%newice   ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: xnewice  ', xnewice      ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: heatOceW ', heatOceW     ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: OceW beg ', xheat        ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: heatOceI ', heatOceI     ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: SnIceAveWE', snowiceave  ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: totalsnowfall',ice%totalsnowfall           , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: snowmelted'   ,snowmelted                  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: SWnetw   ', atmos_fluxes%SWnetw            , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: LWnetw   ', atmos_fluxes%LWnetw            , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: sensw    ', atmos_fluxes%sensw             , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: latw     ', atmos_fluxes%latw              , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: ShortWave', atmos_fluxes%HeatFlux_ShortWave, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: LongWave ', atmos_fluxes%HeatFlux_LongWave , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: Sensible ', atmos_fluxes%HeatFlux_Sensible , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: Latent   ', atmos_fluxes%HeatFlux_Latent   , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: TotalHeat', atmos_fluxes%HeatFlux_Total    , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: TotalIce ', atmos_fluxes%FrshFlux_TotalIce , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: VolumeIce', atmos_fluxes%FrshFlux_VolumeIce, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: h1       ', h1                             , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: h2       ', h2                             , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: h3       ', h3                             , str_module, 4, in_subset=p_patch%cells%owned)

  END SUBROUTINE upper_ocean_TS
  !-------------------------------------------------------------------------------
  !
  !
  !>
  !! !! ice_conc_change: Calculates the changes in concentration as well as the grid-cell average
  !                     thickness of new ice forming in open-water areas
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !! Einar Olason, renamed and added support for changing concentration
  !!
  SUBROUTINE ice_conc_change(p_patch,ice, p_os)

    TYPE(t_patch),             INTENT(IN)    :: p_patch
    TYPE (t_sea_ice),          INTENT(INOUT) :: ice
    TYPE(t_hydro_ocean_state), INTENT(IN)    :: p_os

    INTEGER  :: k
 !  REAL(wp) :: sst(nproma,p_patch%alloc_cell_blocks)
 !  REAL(wp) :: sss(nproma,p_patch%alloc_cell_blocks)
    REAL(wp) :: Tfw(nproma,p_patch%alloc_cell_blocks) ! Ocean freezing temperature [C]

 !  REAL(wp) :: leadclose_2n

    CALL dbg_print('IceConcCh: IceConc beg' ,ice%conc, str_module, 4, in_subset=p_patch%cells%owned)

    ! Calculate the sea surface freezing temperature                        [C]
    IF ( no_tracer < 2 .OR. use_constant_tfreez ) THEN
      Tfw(:,:) = Tf
    ELSE
      Tfw(:,:) = -mu * p_os%p_prog(nold(1))%tracer(:,1,:,2)
    ENDIF

    ! This should not be needed
    ! TODO ram - remove all instances of p_patch%cells%area(:,:) and test
    ! See also dynamics_fem/mo_ice_fem_utils.f90
    DO k=1,ice%kice
      ice%vol (:,k,:) = ice%hi(:,k,:)*ice%conc(:,k,:)*p_patch%cells%area(:,:)
      ice%vols(:,k,:) = ice%hs(:,k,:)*ice%conc(:,k,:)*p_patch%cells%area(:,:)
    ENDDO

    CALL dbg_print('IceConcCh: vol  at beg' ,ice%vol , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: vols at beg' ,ice%vols, str_module, 4, in_subset=p_patch%cells%owned)

!ICON_OMP_PARALLEL
!ICON_OMP_WORKSHARE 
    ! Concentration change due to new ice formation
    WHERE ( ice%newice(:,:) > 0._wp .AND. v_base%lsm_c(:,1,:) <= sea_boundary )
      ! New volume - we just preserve volume:
      ! #slo# 2014-12-11: newice is grown over open ocean but already averaged over whole grid area
      !                   newice must not be multiplied by 1-conc
      ice%vol  (:,1,:) = ice%vol(:,1,:) + ice%newice(:,:)*p_patch%cells%area(:,:)

      ! Hibler's way to change the concentration 
      !  - the formulation here uses the default values of leadclose parameters 2 and 3 in MPIOM:
      !    1 and 0 respectively, which recovers the Hibler model: conc=conc+newice/hnull
      ! Fixed 2. April (2014) - we don't need to multiply with 1-A here, like Hibler does, because it's
      ! already included in newice (we use volume, but Hibler growth rate)
      !ice%conc (:,1,:) = min( 1._wp, ice%conc(:,1,:) + ice%newice(:,:)/hnull )

      ! New formulation of leadclose parameter leadclose_2n includes parameters 2 and 3 of MPIOM:
      ! leadclose_2n (=mpiom_leadclose(3)/mpiom_leadclose(2)
      ! standard value of mpiom is: mpiom_leadclose(3)=2. mpiom_leadclose(2)=mpiom_leadclose(3)+1.
      ! i.e. leadclose_2n=2./3. according to mpiom default
      ice%conc(:,1,:) = min( 1._wp, ice%conc(:,1,:) + &
        &                           ice%newice(:,:)/(hnull+leadclose_2n*(ice%hi(:,1,:)-hnull)) )

      ! New ice and snow thickness
      ice%hi   (:,1,:) = ice%vol (:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
      ice%hs   (:,1,:) = ice%vols(:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
      !TODO: Re-calculate temperatures to conserve energy when we change the ice thickness
    ENDWHERE
!ICON_OMP_END_WORKSHARE

#ifndef _OPENMP
    CALL dbg_print('IceConcCh: conc leadcl' ,ice%conc, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hi   leadcl' ,ice%hi  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hs   leadcl' ,ice%hs  , str_module, 4, in_subset=p_patch%cells%owned)
#endif

!ICON_OMP_WORKSHARE
    ! This is where concentration, and thickness change due to ice melt (we must conserve volume)
    ! A.k.a. lateral melt
    WHERE ( ice%hiold(:,1,:) > ice%hi(:,1,:) .AND. ice%hi(:,1,:) > 0._wp )
      ! Hibler's way to change the concentration due to lateral melting (leadclose parameter 1)
      ice%conc(:,1,:) = MAX( 0._wp, ice%conc(:,1,:) &
        &        - ( ice%hiold(:,1,:)-ice%hi(:,1,:) )*ice%conc(:,1,:)*leadclose_1/ice%hiold(:,1,:) )

      ! New ice and snow thickness
      ice%hi  (:,1,:) = ice%vol (:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
      ice%hs  (:,1,:) = ice%vols(:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
      !TODO: Re-calculate temperatures to conserve energy when we change the ice thickness
    ENDWHERE
!ICON_OMP_END_WORKSHARE

#ifndef _OPENMP
    CALL dbg_print('IceConcCh: conc latMlt' ,ice%conc, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hi   latMlt' ,ice%hi  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hs   latMlt' ,ice%hs  , str_module, 4, in_subset=p_patch%cells%owned)
#endif

    ! Ice cannot grow thinner than hmin
    ! Changed 27. March
!ICON_OMP_WORKSHARE
    WHERE ( ice%hi(:,1,:) < hmin .AND. ice%hi(:,1,:) > 0._wp )
      ice%hi  (:,1,:) = hmin
      ice%conc(:,1,:) = ice%vol(:,1,:) / ( ice%hi(:,1,:)*p_patch%cells%area(:,:) )
      ice%hs  (:,1,:) = ice%vols(:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
    ENDWHERE
!ICON_OMP_END_WORKSHARE

!ICON_OMP_WORKSHARE
    WHERE (ice%hi(:,1,:) <= 0._wp)
      ice%Tsurf(:,1,:) = Tfw(:,:)
      ice%T1   (:,1,:) = Tfw(:,:)
      ice%T2   (:,1,:) = Tfw(:,:)
      ice%conc (:,1,:) = 0.0_wp
      ice%hi   (:,1,:) = 0.0_wp
      ice%hs   (:,1,:) = 0.0_wp
      ice%E1   (:,1,:) = 0.0_wp
      ice%E2   (:,1,:) = 0.0_wp
      ice%vol  (:,1,:) = 0.0_wp
    ENDWHERE
!ICON_OMP_END_WORKSHARE
!ICON_OMP_END_PARALLEL

    ice%concSum(:,:)  = SUM(ice%conc(:,:,:),2)

    CALL dbg_print('IceConcCh: IceConc end' ,ice%conc, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hi   at end' ,ice%hi  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hs   at end' ,ice%hs  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: vol  at end' ,ice%vol , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: vols at end' ,ice%vols, str_module, 4, in_subset=p_patch%cells%owned)

  END SUBROUTINE ice_conc_change


  !-------------------------------------------------------------------------
  !
  !> Forcing_from_bulk equals sbr "Budget_omip" in MPIOM.
  !! Sets the atmospheric fluxes over *SEA ICE ONLY* for the update of the ice
  !! temperature and ice growth rates for OMIP forcing
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-07). Originally code written by
  !! Dirk Notz, following MPIOM. Code transfered to ICON.
  !! Einar Olason, split calc_atm_fluxes_from_bulk into calc_bulk_flux_ice and calc_bulk_flux_oce
  !! so that the ocean model can be run without the ice model, but with OMIP fluxes.
  !
  !  INPUT variables:
  !   p_as%pao
  !   p_as%ftdew
  !   tracer(:,1,:,1)  :  SST
  !
  !  OUTPUT variables:
  !   atmos_fluxes             :  heat fluxes, derivatives, wind stress
  !
  SUBROUTINE calc_bulk_flux_ice(p_patch, p_as, p_ice, atmos_fluxes)
    TYPE(t_patch),            INTENT(IN), TARGET    :: p_patch
    TYPE(t_atmos_for_ocean),  INTENT(IN)    :: p_as
    TYPE(t_sea_ice),          INTENT(IN)    :: p_ice
    TYPE(t_atmos_fluxes),     INTENT(INOUT) :: atmos_fluxes

    !Local variables
    REAL(wp), DIMENSION (nproma,p_patch%alloc_cell_blocks) ::           &
      & Tsurf,          &  ! Surface temperature                             [C]
      & tafoK,          &  ! Air temperature at 2 m in Kelvin                [K]
      & fu10lim,        &  ! wind speed at 10 m height in range 2.5...32     [m/s]
      & esta,           &  ! water vapor pressure at 2 m height              [Pa]
      & esti,           &  ! water vapor pressure at ice surface             [Pa]
      & sphumida,       &  ! Specific humididty at 2 m height
      & sphumidi,       &  ! Specific humididty at ice surface
      & ftdewC,         &  ! Dew point temperature in Celsius                [C]
      & rhoair,         &  ! air density                                     [kg/m^3]
      & dragl0,         &  ! part of dragl
      & dragl1,         &  ! part of dragl
      & dragl,          &  ! Drag coefficient for latent   heat flux
      & drags,          &  ! Drag coefficient for sensible heat flux (=0.95 dragl)
      & fakts,          &  ! Effect of cloudiness on LW radiation
      & humi,           &  ! Effect of air humidity on LW radiation
      & fa, fi,         &  ! Enhancment factor for vapor pressure
      & dsphumididesti, &  ! Derivative of sphumidi w.r.t. esti
      & destidT,        &  ! Derivative of esti w.r.t. T
      & dfdT,           &  ! Derivative of f w.r.t. T
      & wspeed             ! Wind speed                                      [m/s]

    INTEGER :: i, jb, jc, i_startidx_c, i_endidx_c
    REAL(wp) :: aw,bw,cw,dw,ai,bi,ci,di,AAw,BBw,CCw,AAi,BBi,CCi,alpha,beta
    REAL(wp) :: fvisdir, fvisdif, fnirdir, fnirdif
    ! For wind-stress ramping
    REAL(wp) :: ramp = 1.0_wp

    TYPE(t_subset_range), POINTER :: all_cells

    !CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:calc_bulk_flux_ice'
    !-------------------------------------------------------------------------
    !CALL message(TRIM(routine), 'start' )

    tafoK(:,:)  = p_as%tafo(:,:)  + tmelt  ! Change units of tafo  to Kelvin
    ftdewC(:,:) = p_as%ftdew(:,:) - tmelt  ! Change units of ftdew to C

    ! set to zero for NAG, for debug necessary only
    IF (idbg_mxmn > 3 .OR. idbg_val>3) THEN
      fi      (:,:) = 0.0_wp
      esti    (:,:) = 0.0_wp
      sphumidi(:,:) = 0.0_wp
      dragl   (:,:) = 0.0_wp
      drags   (:,:) = 0.0_wp
      dfdT    (:,:) = 0.0_wp
      destidT (:,:) = 0.0_wp
      dsphumididesti(:,:) = 0.0_wp
    ENDIF

    ! subset range pointer
    all_cells => p_patch%cells%all

    !-----------------------------------------------------------------------
    ! Compute water vapor pressure and specific humididty in 2m height (esta)
    ! and at water surface (estw) according to "Buck Research Manual (1996)
    ! (see manuals for instruments at http://www.buck-research.com/);
    ! updated from Buck, A. L., New equations for computing vapor pressure and
    ! enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981"
    !-----------------------------------------------------------------------
    ! #slo# 2015-03: above comment now valid - the values for open water below are from Buck (1981)
    !                the values for ice are not changed in Buck (1996) in comparison to Buck (1981)
  ! aw=611.21_wp; bw=18.729_wp; cw=257.87_wp; dw=227.3_wp
    ai=611.15_wp; bi=23.036_wp; ci=279.82_wp; di=333.7_wp
    ! here are the updated values for open water according to Buck (1996)
    aw=611.21_wp; bw=18.678_wp; cw=257.14_wp; dw=234.5_wp

    AAw=7.2e-4_wp; BBw=3.20e-6_wp; CCw=5.9e-10_wp
    AAi=2.2e-4_wp; BBi=3.83e-6_wp; CCi=6.4e-10_wp

    alpha=0.62197_wp; beta=0.37803_wp

    ! #slo# correction: pressure in enhancement formula is in mb (hPa) according to Buck 1981 and 1996
    fa(:,:)        = 1.0_wp+AAw+p_as%pao*0.01_wp*(BBw+CCw*ftdewC**2)
    esta(:,:)      = fa * aw*EXP((bw-ftdewC/dw)*ftdewC/(ftdewC+cw))

    sphumida(:,:)  = alpha * esta/(p_as%pao-beta*esta)
    !-----------------------------------------------------------------------
    !  Compute longwave radiation according to
    !         Berliand, M. E., and T. G. Berliand, 1952: Determining the net
    !         long-wave radiation of the Earth with consideration of the effect
    !         of cloudiness. Izv. Akad. Nauk SSSR, Ser. Geofiz., 1, 6478.
    !         cited by: Budyko, Climate and Life, 1974.
    !         Note that for humi, esta is given in [mmHg] in the original
    !         publication. Therefore, 0.05*sqrt(esta/100) is used rather than
    !         0.058*sqrt(esta)
    !  This is the formula used in MPI-OM when using the QLOBERL preprocessing option (currently
    !  the default usage).
    !-----------------------------------------------------------------------

    ! NB: Lwinw and LWoutw is a misleading nomenclature in this case, since
    ! Berliand & Berliand ('52) calculate only LWnet
    humi    = 0.39_wp - 0.05_wp*SQRT(esta/100._wp)
    fakts   =  1.0_wp - ( 0.5_wp + 0.4_wp/90._wp &
      &         *MIN(ABS(rad2deg*p_patch%cells%center(:,:)%lat),60._wp) ) * p_as%fclou**2
    atmos_fluxes%LWin(:,:) = 0._wp
    atmos_fluxes%LWout(:,:,:) = 0._wp
    !-----------------------------------------------------------------------
    !  Calculate bulk equations according to
    !      Kara, B. A., P. A. Rochford, and H. E. Hurlburt, 2002:
    !      Air-Sea Flux Estimates And The 19971998 Enso Event,  Bound.-Lay.
    !      Met., 103(3), 439-458, doi: 10.1023/A:1014945408605.
    !-----------------------------------------------------------------------

    rhoair(:,:) = 0._wp
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c,i_endidx_c

        rhoair(jc,jb) = p_as%pao(jc,jb)                &
          &            /(rd*tafoK(jc,jb)*(1.0_wp+0.61_wp*sphumida(jc,jb)) )

      END DO
    END DO

    fu10lim(:,:)    = MAX (2.5_wp, MIN(32.5_wp,p_as%fu10(:,:)) )
    dragl1(:,:)     = 1e-3_wp*(-0.0154_wp + 0.5698_wp/fu10lim(:,:) &
      &               - 0.6743_wp/(fu10lim(:,:) * fu10lim(:,:)))
    dragl0(:,:)     = 1e-3_wp*(0.8195_wp+0.0506_wp*fu10lim(:,:) &
      &               - 0.0009_wp*fu10lim(:,:)*fu10lim(:,:))

    ! Fractions of SWin in each band (from cice)
    fvisdir=0.28_wp; fvisdif=0.24_wp; fnirdir=0.31_wp; fnirdif=0.17_wp
    Tsurf(:,:) = 0._wp ! For debug output

    ! Over sea ice area only
    DO i = 1, p_ice%kice
      WHERE (p_ice%hi(:,i,:)>0._wp)
        atmos_fluxes%SWnet(:,i,:) = ( 1._wp-atmos_fluxes%albvisdir(:,i,:) )*fvisdir*p_as%fswr(:,:) +   &
          &                 ( 1._wp-atmos_fluxes%albvisdif(:,i,:) )*fvisdif*p_as%fswr(:,:) +   &
          &                 ( 1._wp-atmos_fluxes%albnirdir(:,i,:) )*fnirdir*p_as%fswr(:,:) +   &
          &                 ( 1._wp-atmos_fluxes%albnirdif(:,i,:) )*fnirdif*p_as%fswr(:,:)
        Tsurf(:,:)    = p_ice%Tsurf(:,i,:)
        ! #slo# correction: pressure in enhancement formula is in mb (hPa) according to Buck 1981 and 1996
        !fi(:,:)       = 1.0_wp+AAi+p_as%pao(:,:)*(BBi+CCi*Tsurf(:,:) **2)
        fi(:,:)       = 1.0_wp+AAi+p_as%pao(:,:)*0.01_wp*(BBi+CCi*Tsurf(:,:) **2)
        esti(:,:)     = fi(:,:)*ai*EXP((bi-Tsurf(:,:) /di)*Tsurf(:,:) /(Tsurf(:,:) +ci))
        sphumidi(:,:) = alpha*esti(:,:)/(p_as%pao(:,:)-beta*esti(:,:))
        ! This may not be the best drag parametrisation to use over ice
        dragl(:,:)    = dragl0(:,:) + dragl1(:,:) * (Tsurf(:,:)-p_as%tafo(:,:))
        ! A reasonable maximum and minimum is needed for dragl in case there's a large difference
        ! between the 2-m and surface temperatures.
        dragl(:,:)    = MAX(0.5e-3_wp, MIN(3.0e-3_wp,dragl(:,:)))
        drags(:,:)    = 0.95_wp * dragl(:,:)

        ! #eoo# 2012-12-14: another bugfix
        ! #slo# 2012-12-13: bugfix, corrected form
        atmos_fluxes%LWnet (:,i,:)  = - fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4 &
           &                  - 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:))
        ! same form as MPIOM:
        !atmos_fluxes%LWnet (:,i,:)  = - (fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4 &
        !  &         + 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:)))
        ! bug
        !atmos_fluxes%LWnet (:,i,:)  = fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4 &
        !  &     - 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:))
        atmos_fluxes%dLWdT (:,i,:)  = -4._wp*zemiss_def*stbo*tafoK(:,:)**3
        atmos_fluxes%sens  (:,i,:)  = drags(:,:) * rhoair(:,:)*cpd*p_as%fu10(:,:) * fr_fac &
          &                    * (p_as%tafo(:,:) -Tsurf(:,:))
        atmos_fluxes%lat   (:,i,:)  = dragl(:,:) * rhoair(:,:)* alf *p_as%fu10(:,:) * fr_fac &
          &                   * (sphumida(:,:)-sphumidi(:,:))

        atmos_fluxes%dsensdT(:,i,:) = 0.95_wp*cpd*rhoair(:,:)*p_as%fu10(:,:)&
          &                  *(dragl0(:,:) - 2.0_wp*dragl(:,:))
        dsphumididesti(:,:) = alpha/(p_as%pao(:,:)-beta*esti(:,:)) &
          &                   * (1.0_wp + beta*esti(:,:)/(p_as%pao(:,:)-beta*esti(:,:)))
        destidT(:,:)        = (bi*ci*di-Tsurf(:,:)*(2.0_wp*ci+Tsurf(:,:)))&
          &                   /(di*(ci+Tsurf(:,:))**2) * esti(:,:)
        dfdT(:,:)               = 2.0_wp*CCi*BBi*Tsurf(:,:)
        atmos_fluxes%dlatdT(:,i,:)  = alf*rhoair(:,:)*p_as%fu10(:,:)* &
          &                  ( (sphumida(:,:)-sphumidi(:,:))*dragl1(:,:) &
          &                    - dragl(:,:)*dsphumididesti(:,:)*(fi(:,:)*destidT(:,:) &
          &                    + esti(:,:)*dfdT(:,:)) )
      ENDWHERE
    ENDDO

    !Dirk: why zero ?
    atmos_fluxes%rpreci(:,:) = 0.0_wp
    atmos_fluxes%rprecw(:,:) = 0.0_wp

    IF (use_calculated_ocean_stress) THEN
      !-----------------------------------------------------------------------
      !  Calculate ice wind stress
      !-----------------------------------------------------------------------
      wspeed(:,:) = SQRT( p_as%u**2 + p_as%v**2 )
      atmos_fluxes%stress_x(:,:) = Cd_ia*rhoair(:,:)*wspeed(:,:)*p_as%u(:,:)
      atmos_fluxes%stress_y(:,:) = Cd_ia*rhoair(:,:)*wspeed(:,:)*p_as%v(:,:)
    ELSE
      ! use wind stress provided by OMIP data
      atmos_fluxes%stress_x(:,:) = p_as%topBoundCond_windStress_u(:,:)
      atmos_fluxes%stress_y(:,:) = p_as%topBoundCond_windStress_v(:,:)
    ENDIF

    ! Ramp for wind-stress - needed for ice-ocean momentum coupling during spinup
  ! IF ( PRESENT(datetime) ) THEN
  !   ramp = MIN(1._wp,(datetime%calday + datetime%caltime &
  !     - time_config%ini_datetime%calday - time_config%ini_datetime%caltime) / ramp_wind)
  !   IF (idbg_mxmn > 3 .OR. idbg_val>3) WRITE(0,*) ' RAMP = ',ramp
  !   atmos_fluxes%stress_x(:,:)  = ramp*atmos_fluxes%stress_x(:,:)
  !   atmos_fluxes%stress_y(:,:)  = ramp*atmos_fluxes%stress_y(:,:)
  ! ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('CalcBulkI:tafoK'           ,tafoK    ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:rhoair'          ,rhoair   ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:fa'              ,fa       ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:fi'              ,fi       ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:esta'            ,esta     ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:esti'            ,esti     ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:sphumida'        ,sphumida ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:sphumidi'        ,sphumidi ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:dragl'           ,dragl    ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:drags'           ,drags    ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:stress_x'        ,atmos_fluxes%stress_x,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:stress_y'        ,atmos_fluxes%stress_y,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:atmflx%lat ice'  ,atmos_fluxes%lat     ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:atmflx%dsensdT'  ,atmos_fluxes%dsensdT ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:atmflx%dlatdT'   ,atmos_fluxes%dlatdT  ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:dsphumididesti'  ,dsphumididesti       ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:destidT'         ,destidT              ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:dfdT'            ,dfdT                 ,str_module,idt_src, in_subset=p_patch%cells%owned)
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('CalcBulkI:Tsurf ice'       , Tsurf               ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:atmflx%LWnet ice', atmos_fluxes%LWnet  ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:atmflx%sens ice' , atmos_fluxes%sens   ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkI:atmflx%lat ice'  , atmos_fluxes%lat    ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE calc_bulk_flux_ice

  !-------------------------------------------------------------------------
  !
  !> Forcing_from_bulk equals sbr "Budget_omip" in MPIOM.
  !! Sets the atmospheric fluxes for the update of
  !! temperature of *OPEN WATER* for OMIP forcing.
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2012-08). Originally code written by
  !! Dirk Notz, following MPIOM. Code transfered to ICON.
  !
  !  INPUT variables:
  !   p_as%pao
  !   p_as%ftdew
  !   tracer(:,1,:,1)  :  SST
  !  OUTPUT variables:
  !   atmos_fluxes             :  heat fluxes, derivatives, wind stress
   
  SUBROUTINE calc_bulk_flux_oce(p_patch, p_as, p_os, atmos_fluxes)
    TYPE(t_patch),            INTENT(IN), TARGET    :: p_patch
    TYPE(t_atmos_for_ocean),  INTENT(IN)    :: p_as
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os
    TYPE(t_atmos_fluxes),     INTENT(INOUT) :: atmos_fluxes

    !Local variables
    REAL(wp), DIMENSION (nproma,p_patch%alloc_cell_blocks) ::           &
      & Tsurf,          &  ! Surface temperature                             [C]
      & tafoK,          &  ! Air temperature at 2 m in Kelvin                [K]
      & fu10lim,        &  ! wind speed at 10 m height in range 2.5...32     [m/s]
      & esta,           &  ! water vapor pressure at 2 m height              [Pa]
      & estw,           &  ! water vapor pressure at water surface           [Pa]
      & sphumida,       &  ! Specific humididty at 2 m height
      & sphumidw,       &  ! Specific humididty at water surface
      & ftdewC,         &  ! Dew point temperature in Celsius                [C]
      & rhoair,         &  ! air density                                     [kg/m^3]
      & dragl0,         &  ! part of dragl
      & dragl1,         &  ! part of dragl
      & dragl,          &  ! Drag coefficient for latent   heat flux
      & drags,          &  ! Drag coefficient for sensible heat flux (=0.95 dragl)
      & fakts,          &  ! Effect of cloudiness on LW radiation
      & humi,           &  ! Effect of air humidity on LW radiation
      & fa, fw,         &  ! Enhancment factor for vapor pressure
      & wspeed,         &  ! Wind speed                                      [m/s]
      & C_ao               ! Drag coefficient for atm-ocean stress           [m/s]

    INTEGER :: jb, jc, i_startidx_c, i_endidx_c
    REAL(wp) :: aw,bw,cw,dw,AAw,BBw,CCw,alpha,beta
    REAL(wp) :: fvisdir, fvisdif, fnirdir, fnirdif
    ! For wind-stress ramping
    REAL(wp) :: ramp = 1.0_wp

    TYPE(t_subset_range), POINTER :: all_cells

    !CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'sea_ice:calc_bulk_flux_oce
    !-------------------------------------------------------------------------
    !CALL message(TRIM(routine), 'start' )

    Tsurf(:,:)  = p_os%p_prog(nold(1))%tracer(:,1,:,1)  ! set surface temp = mixed layer temp
    tafoK(:,:)  = p_as%tafo(:,:)  + tmelt               ! Change units of tafo  to Kelvin
    ftdewC(:,:) = p_as%ftdew(:,:) - tmelt               ! Change units of ftdew to Celsius

    ! subset range pointer
    all_cells => p_patch%cells%all



    !-----------------------------------------------------------------------
    ! Compute water vapor pressure and specific humididty in 2m height (esta)
    ! and at water surface (estw) according to "Buck Research Manual (1996)
    ! (see manuals for instruments at http://www.buck-research.com/);
    ! updated from Buck, A. L., New equations for computing vapor pressure and
    ! enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981"
    !-----------------------------------------------------------------------
    ! #slo# 2015-03: above comment now valid - the values below are from Buck (1981)
  ! aw    = 611.21_wp; bw    = 18.729_wp;  cw  = 257.87_wp; dw = 227.3_wp
    AAw   = 7.2e-4_wp; BBw   = 3.20e-6_wp; CCw = 5.9e-10_wp
    alpha = 0.62197_wp; beta = 0.37803_wp
    ! here are the updated values according to Buck (1996)
    aw    = 611.21_wp; bw    = 18.678_wp;  cw  = 257.14_wp; dw = 234.5_wp

    ! #slo# correction: pressure in enhancement formula is in mb (hPa) according to Buck 1981 and 1996
   !fa(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*(BBw+CCw*ftdewC(:,:)**2)
    fa(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*0.01_wp*(BBw+CCw*ftdewC(:,:)**2)
    esta(:,:) = fa(:,:) * aw*EXP((bw-ftdewC(:,:)/dw)*ftdewC(:,:)/(ftdewC(:,:)+cw))
   !esta(:,:) =           aw*EXP((bw-ftdewC(:,:)/dw)*ftdewC(:,:)/(ftdewC(:,:)+cw))
   !fw(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*(BBw+CCw*Tsurf(:,:) **2)
    fw(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*0.01_wp*(BBw+CCw*Tsurf(:,:) **2)
   !estw(:,:) = fw(:,:) *aw*EXP((bw-Tsurf(:,:) /dw)*Tsurf(:,:) /(Tsurf(:,:) +cw))
    ! For a given surface salinity we should multiply estw with  1 - 0.000537*S
    ! #slo# correction according to MPIOM: lowering of saturation vapor pressure over saline water
    !       is taken constant to 0.9815
    estw(:,:) = 0.9815_wp*fw(:,:)*aw*EXP((bw-Tsurf(:,:) /dw)*Tsurf(:,:) /(Tsurf(:,:) +cw))

    sphumida(:,:)  = alpha * esta(:,:)/(p_as%pao(:,:)-beta*esta(:,:))
    sphumidw(:,:)  = alpha * estw(:,:)/(p_as%pao(:,:)-beta*estw(:,:))

    !-----------------------------------------------------------------------
    !  Compute longwave radiation according to
    !         Berliand, M. E., and T. G. Berliand, 1952: Determining the net
    !         long-wave radiation of the Earth with consideration of the effect
    !         of cloudiness. Izv. Akad. Nauk SSSR, Ser. Geofiz., 1, 6478.
    !         cited by: Budyko, Climate and Life, 1974.
    !         Note that for humi, esta is given in [mmHg] in the original
    !         publication. Therefore, 0.05*sqrt(esta/100) is used rather than
    !         0.058*sqrt(esta)
    !  This is the formula used in MPI-OM when using the QLOBERL preprocessing option (currently
    !  the default usage).
    !-----------------------------------------------------------------------

    humi(:,:)    = 0.39_wp - 0.05_wp*SQRT(esta(:,:)/100._wp)
    fakts(:,:)   =  1.0_wp - ( 0.5_wp + 0.4_wp/90._wp &
      &         *MIN(ABS(rad2deg*p_patch%cells%center(:,:)%lat),60._wp) ) * p_as%fclou(:,:)**2
    ! NB: Lwin and LWoutw is a misleading nomenclature in this case, since
    ! Berliand & Berliand ('52) calculate only LWnetw
    atmos_fluxes%LWin(:,:) = 0._wp
    atmos_fluxes%LWoutw(:,:) = 0._wp

    ! #eoo# 2012-12-14: another bugfix
    ! #slo# #hha# 2012-12-13: bugfix, corrected form
    atmos_fluxes%LWnetw(:,:) = - fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4  &
      &                - 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:))
    ! same form as MPIOM:
    !atmos_fluxes%LWnetw(:,:) = - (fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4  &
    !  &         + 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:)))
    ! bug
    !atmos_fluxes%LWnetw(:,:) = fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4  &
    !  &         - 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:))

    ! Fractions of SWin in each band (from cice)
    fvisdir=0.28_wp; fvisdif=0.24_wp; fnirdir=0.31_wp; fnirdif=0.17_wp
    atmos_fluxes%SWnetw(:,:) = ( 1._wp-atmos_fluxes%albvisdirw(:,:) )*fvisdir*p_as%fswr(:,:) +   &
      &                ( 1._wp-atmos_fluxes%albvisdifw(:,:) )*fvisdif*p_as%fswr(:,:) +   &
      &                ( 1._wp-atmos_fluxes%albnirdirw(:,:) )*fnirdir*p_as%fswr(:,:) +   &
      &                ( 1._wp-atmos_fluxes%albnirdifw(:,:) )*fnirdif*p_as%fswr(:,:)

    !-----------------------------------------------------------------------
    !  Calculate bulk equations according to
    !      Kara, B. A., P. A. Rochford, and H. E. Hurlburt, 2002:
    !      Air-Sea Flux Estimates And The 19971998 Enso Event,  Bound.-Lay.
    !      Met., 103(3), 439-458, doi: 10.1023/A:1014945408605.
    !-----------------------------------------------------------------------

    rhoair(:,:) = 0._wp
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c,i_endidx_c

        rhoair(jc,jb) = p_as%pao(jc,jb)                &
          &            /(rd*tafoK(jc,jb)*(1.0_wp+0.61_wp*sphumida(jc,jb)) )

      END DO
    END DO

    fu10lim(:,:)    = MAX (2.5_wp, MIN(32.5_wp,p_as%fu10(:,:)) )
    dragl1(:,:)     = 1e-3_wp*(-0.0154_wp + 0.5698_wp/fu10lim(:,:) &
      &               - 0.6743_wp/(fu10lim(:,:) * fu10lim(:,:)))
    dragl0(:,:)     = 1e-3_wp*(0.8195_wp+0.0506_wp*fu10lim(:,:) &
      &               - 0.0009_wp*fu10lim(:,:)*fu10lim(:,:))
    dragl(:,:)      = dragl0(:,:) + dragl1(:,:) * (Tsurf(:,:)-p_as%tafo(:,:))
    ! A reasonable maximum and minimum is needed for dragl in case there's a large difference
    ! between the 2-m and surface temperatures.
    dragl(:,:)      = MAX(0.5e-3_wp, MIN(3.0e-3_wp,dragl(:,:)))
    drags(:,:)      = 0.95_wp * dragl(:,:)
    atmos_fluxes%sensw(:,:) = drags(:,:)*rhoair(:,:)*cpd*p_as%fu10(:,:) * fr_fac &
      &               * (p_as%tafo(:,:) -Tsurf(:,:))
    atmos_fluxes%latw(:,:)  = dragl(:,:)*rhoair(:,:)*alv*p_as%fu10(:,:) * fr_fac &
      &               * (sphumida(:,:)-sphumidw(:,:))

    IF (use_calculated_ocean_stress) THEN
      !-----------------------------------------------------------------------
      !  Calculate oceanic wind stress according to:
      !   Gill (Atmosphere-Ocean Dynamics, 1982, Academic Press) (see also Smith, 1980, J. Phys
      !   Oceanogr., 10, 709-726)
      !-----------------------------------------------------------------------

      wspeed(:,:) = SQRT( p_as%u**2 + p_as%v**2 )
      C_ao(:,:)   = MIN( 2._wp, MAX(1.1_wp, 0.61_wp+0.063_wp*wspeed ) )*1e-3_wp
      atmos_fluxes%stress_xw(:,:) = C_ao(:,:)*rhoair*wspeed(:,:)*p_as%u(:,:)
      atmos_fluxes%stress_yw(:,:) = C_ao(:,:)*rhoair*wspeed(:,:)*p_as%v(:,:)
    ELSE
      ! use wind stress provided by OMIP data
      atmos_fluxes%stress_xw(:,:) = p_as%topBoundCond_windStress_u(:,:)
      atmos_fluxes%stress_yw(:,:) = p_as%topBoundCond_windStress_v(:,:)
    ENDIF

    ! Ramp for wind-stress - needed for ice-ocean momentum coupling during spinup
  ! IF ( PRESENT(datetime) ) THEN
  !   ramp = MIN(1._wp,(datetime%calday + datetime%caltime &
  !     - time_config%ini_datetime%calday - time_config%ini_datetime%caltime) / ramp_wind)
  !   IF (idbg_mxmn > 3 .OR. idbg_val>3) WRITE(0,*) ' RAMP = ',ramp
  !   atmos_fluxes%stress_xw(:,:) = ramp*atmos_fluxes%stress_xw(:,:)
  !   atmos_fluxes%stress_yw(:,:) = ramp*atmos_fluxes%stress_yw(:,:)
  ! ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5          , fix)
    CALL dbg_print('CalcBulkO:tafoK'              , tafoK                 , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:tafo'               , p_as%tafo             , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:ftdew'              , p_as%ftdew            , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:ftdewC'             , ftdewC                , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:pao'                , p_as%pao              , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:fa'                 , fa                    , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:fw'                 , fw                    , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:esta'               , esta                  , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:estw'               , estw                  , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:sphumida'           , sphumida              , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:sphumidw'           , sphumidw              , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:rhoair'             , rhoair                , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:dragl'              , dragl                 , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:drags'              , drags                 , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:fu10'               , p_as%fu10             , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:fu10lim'            , fu10lim               , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:stress_xw'          , atmos_fluxes%stress_xw, str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:stress_yw'          , atmos_fluxes%stress_yw, str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:p_as%windStr-u',p_as%topBoundCond_windStress_u,str_module,idt_src, in_subset=p_patch%cells%owned)
    idt_src=3  ! output print level (1-5          , fix)
    CALL dbg_print('CalcBulkO:Tsurf ocean'        , Tsurf                 , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:atmflx%LWnetw'      , atmos_fluxes%LWnetw   , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:atmflx%sensw'       , atmos_fluxes%sensw    , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulkO:atmflx%latw'        , atmos_fluxes%latw     , str_module, idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE calc_bulk_flux_oce

  SUBROUTINE update_ice_statistic(p_acc, p_ice, subset)
    TYPE(t_sea_ice_acc),  INTENT(INOUT) :: p_acc
    TYPE(t_sea_ice),      INTENT(IN)    :: p_ice
    TYPE(t_subset_range), INTENT(IN)    :: subset

    CALL add_fields(p_acc%hi  , p_ice%hi  , subset , levels=p_ice%kice)
    CALL add_fields(p_acc%hs  , p_ice%hs  , subset , levels=p_ice%kice)
    CALL add_fields(p_acc%conc, p_ice%conc, subset , levels=p_ice%kice)
    CALL add_fields(p_acc%u   , p_ice%u   , subset)
    CALL add_fields(p_acc%v   , p_ice%v   , subset)
  END SUBROUTINE update_ice_statistic
  SUBROUTINE compute_mean_ice_statistics(p_acc,nsteps_since_last_output)
    TYPE(t_sea_ice_acc), INTENT(INOUT) :: p_acc
    INTEGER,INTENT(IN)                 :: nsteps_since_last_output

    p_acc%hi                        = p_acc%hi  /REAL(nsteps_since_last_output,wp)
    p_acc%hs                        = p_acc%hs  /REAL(nsteps_since_last_output,wp)
    p_acc%u                         = p_acc%u   /REAL(nsteps_since_last_output,wp)
    p_acc%v                         = p_acc%v   /REAL(nsteps_since_last_output,wp)
    p_acc%conc                      = p_acc%conc/REAL(nsteps_since_last_output,wp)
  END SUBROUTINE compute_mean_ice_statistics
  SUBROUTINE reset_ice_statistics(p_acc)
    TYPE(t_sea_ice_acc), INTENT(INOUT) :: p_acc
    p_acc%hi                        = 0.0_wp
    p_acc%hs                        = 0.0_wp
    p_acc%u                         = 0.0_wp
    p_acc%v                         = 0.0_wp
    p_acc%conc                      = 0.0_wp
  END SUBROUTINE reset_ice_statistics
  SUBROUTINE finish_unless_allocate(ist, routine, tag)
    INTEGER :: ist
    CHARACTER(len=*) :: tag,routine
    IF (ist /= SUCCESS) CALL finish(TRIM(routine),'allocation of '//TRIM(tag)//' failed.')
  END SUBROUTINE finish_unless_allocate

  ! compute the salt content in the upper most layer based on the liquid water height from the ice model: zUnderIce
  FUNCTION salt_content_in_surface(p_patch, thickness, p_ice, p_os, surface_fluxes, zUnderIceOld,computation_type, info) &
      & RESULT(salt)
    TYPE(t_patch),POINTER                                 :: p_patch
    REAL(wp),DIMENSION(nproma,p_patch%alloc_cell_blocks), &
      & INTENT(IN)                                        :: thickness,zUnderIceOld
    TYPE (t_sea_ice),       INTENT(INOUT)                 :: p_ice
    TYPE(t_hydro_ocean_state)                             :: p_os
    TYPE(t_sfc_flx)                                       :: surface_fluxes
    INTEGER,INTENT(IN), OPTIONAL                          :: computation_type
    CHARACTER(len=*) , OPTIONAL                           :: info

    ! locals
    REAL(wp), DIMENSION(nproma,p_patch%alloc_cell_blocks) :: salt, salinityDiff
    REAL(wp), DIMENSION(nproma,p_patch%alloc_cell_blocks) :: saltInSeaice, saltInLiquidWater
    INTEGER                                               :: my_computation_type
    CHARACTER(len=20)                                     :: my_info
    TYPE(t_subset_range), POINTER                         :: subset
    INTEGER                                               :: block, cell, cellStart,cellEnd

    my_computation_type = 0
    my_info             = 'BEFORE'

    salinityDiff = 0.0_wp
    salt         = 0.0_wp

    ! exit if actual debug-level is < 2
    IF (idbg_mxmn < 2 .AND. idbg_val < 2) RETURN

    CALL assign_if_present(my_computation_type, computation_type)
    CALL assign_if_present(my_info, info)

    subset => p_patch%cells%owned
    DO block = subset%start_block, subset%end_block
      CALL get_index_range(subset, block, cellStart, cellEnd)
      DO cell = cellStart, cellEnd
        IF (subset%vertical_levels(cell,block) < 1) CYCLE
        SELECT CASE (my_computation_type)
        CASE (0)
          ! compute salt amount in the first layer
          saltInSeaice(cell,block)      = sice &
            &                    * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                    * p_patch%cells%area(cell,block)
          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * p_ice%zUnderIce(cell,block) &
            &                    * p_patch%cells%area(cell,block)
        CASE (1)
          ! compute salt amount in the first layer
          saltInSeaice(cell,block)      = sice &
            &                    * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                    * p_patch%cells%area(cell,block)
    !     saltInSeaice(cell,block)      = 0.0_wp
          salinityDiff(cell,block)      = surface_fluxes%FrshFlux_TotalSalt(cell,block)*(dtime/(thickness(cell,block) &
            &                                                                    + p_os%p_prog(nold(1))%h(cell,block)))
          p_os%p_prog(nold(1))%tracer(cell,1,block,2)  = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                                   + salinityDiff(cell,block)
          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * p_ice%zUnderIce(cell,block) &
            &                    * p_patch%cells%area(cell,block)
    !     saltInLiquidWater(cell,block) = 0.0_wp
        CASE (2)
          ! compute salt amount in the first layer
          saltInSeaice(cell,block)      = sice &
            &                    * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                    * p_patch%cells%area(cell,block)

          salinityDiff(cell,block)      = surface_fluxes%FrshFlux_TotalSalt(cell,block)*(dtime/p_ice%zUnderIce(cell,block))
          p_os%p_prog(nold(1))%tracer(cell,1,block,2)  = p_os%p_prog(nold(1))%tracer(cell,1,block,2) + salinityDiff(cell,block)

          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * p_ice%zUnderIce(cell,block) &
            &                    * p_patch%cells%area(cell,block)

        CASE (3) ! use zunderIce for volume in tracer change
          saltInSeaice(cell,block)      = sice*rhoi &
            &                    * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                    * p_patch%cells%area(cell,block)

          p_os%p_prog(nold(1))%tracer(cell,1,block,2) = (p_os%p_prog(nold(1))%tracer(cell,1,block,2)*zUnderIceOld(cell,block) &
            &                                            - dtime*surface_fluxes%FrshFlux_TotalSalt(cell,block)) &
            &                                           /p_ice%zUnderIce(cell,block)
          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * p_ice%zUnderIce(cell,block)*rho_ref &
            &                    * p_patch%cells%area(cell,block)
        CASE (4) ! use zunderIce for volume in tracer change, multiply flux with top layer salinity
          saltInSeaice(cell,block)      = sice*rhoi &
            &                    * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                    * p_patch%cells%area(cell,block)

          p_os%p_prog(nold(1))%tracer(cell,1,block,2) = (p_os%p_prog(nold(1))%tracer(cell,1,block,2)*zUnderIceOld(cell,block) &
            &                                            -   dtime &
            &                                              * surface_fluxes%FrshFlux_TotalSalt(cell,block) &
            &                                              * p_os%p_prog(nold(1))%tracer(cell,1,block,2)) &
            &                                           /p_ice%zUnderIce(cell,block)
          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * p_ice%zUnderIce(cell,block)*rho_ref &
            &                    * p_patch%cells%area(cell,block)
        CASE (5) ! use zunderIce for volume in tracer change, multiply flux with top layer salinity
          p_ice%zUnderIce(cell,block) = zUnderIceOld(cell,block)
          saltInSeaice(cell,block)      = sice*rhoi &
            &                    * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                    * p_patch%cells%area(cell,block)
        !!DN This is no longer needed since we now update surface salinity
        !directly
        !!DN   p_os%p_prog(nold(1))%tracer(cell,1,block,2) = (p_os%p_prog(nold(1))%tracer(cell,1,block,2)*zUnderIceOld(cell,block) &
        !!DN   &                                            -   dtime &
        !!DN   &                                              * surface_fluxes%FrshFlux_TotalSalt(cell,block) &
        !!DN   &                                              * p_os%p_prog(nold(1))%tracer(cell,1,block,2)) &
        !!DN   &                                           /p_ice%zUnderIce(cell,block)
          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * p_ice%zUnderIce(cell,block)*rho_ref &
            &                    * p_patch%cells%area(cell,block)
        END SELECT

        salt(cell,block) = saltInSeaice(cell,block) + saltInLiquidWater(cell,block)
      END DO
    END DO

     CALL dbg_print('IceBudget: saltinIce '//TRIM(info)  , &
      &            saltInSeaice , &
      &            str_module, 5, in_subset=p_patch%cells%owned)
     CALL dbg_print('IceBudget: saltinLiquid '//TRIM(info)  , &
      &            saltInLiquidWater , &
      &            str_module, 5, in_subset=p_patch%cells%owned)
     CALL dbg_print('IceBudget: salt '//TRIM(info)  , &
      &            salt , &
      &            str_module, 5, in_subset=p_patch%cells%owned)
     CALL dbg_print('IceBudget: salinityDiff '//TRIM(info)  , &
      &            salinityDiff , &
      &            str_module, 5, in_subset=p_patch%cells%owned)
     CALL dbg_print('IceBudget: zUnderice '//TRIM(info)  , &
      &            p_ice%zUnderIce , &
      &            str_module, 5, in_subset=p_patch%cells%owned)
    !
    ! compute liquid volume in the first layer incl. water prepresentative of sea ice
  END FUNCTION salt_content_in_surface

  ! compute the energy content in the upper most layer based on the liquid water height from the ice model: zUnderIce
  FUNCTION energy_content_in_surface(p_patch, thickness, hold, p_ice, sst, computation_type, info) &
    & RESULT(energy)

    TYPE(t_patch),POINTER                                 :: p_patch
    REAL(wp),DIMENSION(nproma,p_patch%alloc_cell_blocks), INTENT(IN) :: thickness, hold, sst
    TYPE (t_sea_ice), INTENT(IN)                          :: p_ice
    INTEGER,INTENT(IN), OPTIONAL                          :: computation_type
    CHARACTER(len=*) , OPTIONAL                           :: info

    ! locals
    REAL(wp), DIMENSION(nproma,p_patch%alloc_cell_blocks) :: energy
    INTEGER                                               :: my_computation_type
    CHARACTER(len=20)                                     :: my_info
    TYPE(t_subset_range), POINTER                         :: subset
    INTEGER                                               :: block, cell, cellStart,cellEnd
    REAL(wp)                                              :: t_base, zui

    my_computation_type = 0
    my_info             = 'BEFORE'
    energy              = 0.0_wp
  ! t_base              = Tf
  ! t_base              = -5.0_wp
    t_base              = t_heat_base  !  arbitrary temperature basis for calculation of surface heat content

    ! exit if actual debug-level is < 2
    IF (idbg_mxmn < 2 .AND. idbg_val < 2) RETURN

    CALL assign_if_present(my_computation_type, computation_type)
    CALL assign_if_present(my_info, info)

    subset => p_patch%cells%owned
    DO block = subset%start_block, subset%end_block
      CALL get_index_range(subset, block, cellStart, cellEnd)
      DO cell = cellStart, cellEnd
        IF (subset%vertical_levels(cell,block) < 1) CYCLE
        SELECT CASE (my_computation_type)
        CASE (0)
          ! compute energy content of surface layer plus melting energy of ice and snow water equivalent
          !  - relative to arbitrary temperature t_base (e.g. -5C for mostly positive values)
          !  - omit multiplication with area, calculation per unit area, units in Joule/m2
          !  - constant freezing temperature Tf, ice-temperature set to Tf
          !  - use zUnderIce+draftave
          !  = (sst-t_base)*zUnderIce*rhow*clw - draftave*rhow*alf + (Tf-t_base)*(hi*rhoi+hs*rhos)*rhow*clw
          energy(cell,block) = (sst(cell,block) - t_base) * p_ice%zUnderIce(cell,block)*rho_ref*clw &
            &                - (p_ice%draftave(cell,block)*rho_ref*alf) &
            &                + (Tf - t_base)*p_ice%draftave(cell,block)*rho_ref*clw
        CASE (1)
          !  compute energy content - use zUnderIce and hi, hs, conc
          !  = (sst-t_base)*zUnderIce*rhow*clw - (hi*rhoi+hs*rhos)*alf*conc + (Tf-t_base)*(hi*rhoi+hs*rhos)*conc*clw
          energy(cell,block) = (sst(cell,block) - t_base) * p_ice%zUnderIce(cell,block)*rho_ref*clw &
            &                - ((p_ice%hi(cell,1,block)*rhoi + p_ice%hs(cell,1,block)*rhos)*p_ice%conc(cell,1,block)*alf) &
            &                + (Tf - t_base)*(p_ice%hi(cell,1,block)*rhoi + p_ice%hs(cell,1,block)*rhos) &
            &                                *p_ice%conc(cell,1,block)*clw
        CASE (2)
          !  compute energy content - use hi, hs only, compute local zUnderIce
          !  = (sst-t_base)*zUnderIce*rhow*clw - (hi*rhoi+hs*rhos)*alf*conc + (Tf-t_base)*(hi*rhoi+hs*rhos)*conc*clw
          zui                = thickness(cell,block)+hold(cell,block) &
            &                - (rhos * p_ice%hs(cell,1,block) + rhoi * p_ice%hi(cell,1,block)) * p_ice%conc(cell,1,block) / rho_ref
          energy(cell,block) = (sst(cell,block) - t_base) *zui*rho_ref*clw &
            &                - ((p_ice%hi(cell,1,block)*rhoi + p_ice%hs(cell,1,block)*rhos)*p_ice%conc(cell,1,block)*alf) &
            &                + (Tf - t_base)*(p_ice%hi(cell,1,block)*rhoi + p_ice%hs(cell,1,block)*rhos) &
            &                               *p_ice%conc(cell,1,block)*clw
          CONTINUE
        CASE (3)
          write(0,*) " Nothing computed"
        CASE DEFAULT
          CALL finish ('mo_sea_ice:computation_type','option not supported')
        END SELECT
        !salt(cell,block) = saltInSeaice(cell,block) + saltInLiquidWater(cell,block)
      END DO
    END DO

    CALL dbg_print('enContSurf: energy '//TRIM(info),energy,str_module, 5, in_subset=p_patch%cells%owned)

  END FUNCTION energy_content_in_surface

END MODULE mo_sea_ice
