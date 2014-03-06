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
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
  USE mo_dynamics_config,     ONLY: nold, nnew
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_exception,           ONLY: finish, message
  USE mo_impl_constants,      ONLY: success, max_char_length, sea_boundary
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref, ki, ks, Tf, albi, albim, albsm, albs,   &
    &                               fr_fac, mu, alf, alv, albedoW_sim, clw, cpd, zemiss_def, rd, &
    &                               stbo, tmelt, ci, Cd_ia, sice, alb_sno_vis, alb_sno_nir,      &
    &                               alb_ice_vis, alb_ice_nir
  USE mo_math_constants,      ONLY: rad2deg
  USE mo_statistics,          ONLY: add_fields
  USE mo_ocean_nml,           ONLY: no_tracer, use_file_initialConditions, n_zlev
  USE mo_sea_ice_nml,         ONLY: i_ice_therm, i_ice_dyn, ramp_wind, hnull, hmin, hci_layer, &
    &                               i_ice_albedo, leadclose_1
  USE mo_oce_types,           ONLY: t_hydro_ocean_state
  USE mo_oce_state,           ONLY: v_base, &
    &                               ocean_restart_list, set_oce_tracer_info, ocean_default_list
  USE mo_var_list,            ONLY: add_var, add_ref
  USE mo_var_metadata,        ONLY: groups
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants,       ONLY: DATATYPE_FLT32, DATATYPE_PACK16,        &
    &                               GRID_UNSTRUCTURED_CELL, GRID_REFERENCE, &
    &                               GRID_CELL, ZA_GENERIC_ICE, ZA_SURFACE,  &
    &                               GRID_UNSTRUCTURED_VERT, GRID_VERTEX,    &
    &                               GRID_UNSTRUCTURED_EDGE, GRID_EDGE
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, &
    &                               t_atmos_for_ocean, t_sea_ice_acc
  USE mo_sea_ice_winton,      ONLY: ice_growth_winton, set_ice_temp_winton
  USE mo_sea_ice_zerolayer,   ONLY: ice_growth_zerolayer, set_ice_temp_zerolayer
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_ice_fem_utils,       ONLY: fem_ice_wrap, init_fem_wgts, destruct_fem_wgts,             &
    &                               ice_fem_grid_init, ice_fem_grid_post, ice_advection,        &
    &                               ice_ocean_stress
  USE mo_grid_config,         ONLY: n_dom
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_ice_fast, timer_ice_slow
  USE mo_datetime,            ONLY: t_datetime
  USE mo_time_config,         ONLY: time_config

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
  PUBLIC :: destruct_atmos_fluxes

  PUBLIC :: ice_init
!  PUBLIC :: set_ice_albedo
!  PUBLIC :: sum_fluxes
!  PUBLIC :: ave_fluxes
  PUBLIC :: ice_fast
  PUBLIC :: ice_slow
!  PUBLIC :: upper_ocean_TS
  PUBLIC :: calc_bulk_flux_ice
  PUBLIC :: calc_bulk_flux_oce
  PUBLIC :: update_ice_statistic, compute_mean_ice_statistics, reset_ice_statistics

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

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
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    p_patch => p_patch_3D%p_patch_2D(1)
    alloc_cell_blocks = p_patch%alloc_cell_blocks
    nblks_v = p_patch%nblks_v
    nblks_e = p_patch%nblks_e

    p_ice%kice = i_no_ice_thick_class

    CALL add_var(ocean_restart_list, 'alb', p_ice%alb ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('alb', '', 'albedo of snow-ice system', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'Tsurf', p_ice%Tsurf ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('Tsurf', '', 'surface temperature', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'T1', p_ice%T1 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('T1', 'C', 'Temperature upper layer', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'T2', p_ice%T2 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('T2', 'C', 'Temperature lower layer', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'E1', p_ice%E1 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('E1', 'Jm/kg', 'Energy content upper layer', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'E2', p_ice%E2 ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('E2', 'Jm/kg', 'Energy content lower layer', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'vol', p_ice%vol ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('vol', 'm^3', 'ice volume', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'vols', p_ice%vols ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('vols', 'm^3', 'snow volume', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'hi', p_ice%hi ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hi', 'm', 'ice thickness', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'hs', p_ice%hs ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hs', 'm', 'snow thickness', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'hiold', p_ice%hiold ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hiold', 'm', 'ice thickness (last timstep)', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'hsold', p_ice%hsold ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hsold', 'm', 'snow thickness (last timstep)', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'Qtop', p_ice%Qtop ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('Qtop', 'W/m^2', 'Energy flux available for surface melting', &
      &                   DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'Qbot', p_ice%Qbot ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('Qbot', 'W/m^2', 'Energy flux at ice-ocean interface', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'heatocei', p_ice%heatocei ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('heatocei', 'J', 'Energy to ocean when all ice is melted', &
      &                   DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'snow_to_ice', p_ice%snow_to_ice ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('snow_to_ice', 'm', 'amount of snow that is transformed to ice', &
      &                   DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'surfmelt', p_ice%surfmelt ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('surfmelt', 'm', 'surface melt water running into ocean', &
      &                   DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'surfmeltT', p_ice%surfmeltT ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('surfmeltT', 'C', 'Mean temperature of surface melt water', &
      &                   DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'evapwi', p_ice%evapwi ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('evapwi', 'kg/m^2', 'amount of evaporated water if no ice left', &
      &                   DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'conc', p_ice%conc ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('conc', '', 'ice concentration in each ice class', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'ice_u_prog', p_ice%u_prog ,&
      &          GRID_UNSTRUCTURED_VERT, ZA_SURFACE, &
      &          t_cf_var('ice_u_prog', 'm/s', 'zonal velocity', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_VERTEX),&
      &          ldims=(/nproma,nblks_v/), lrestart_cont=.TRUE.)
    CALL add_var(ocean_restart_list, 'ice_v_prog', p_ice%v_prog ,&
      &          GRID_UNSTRUCTURED_VERT, ZA_SURFACE, &
      &          t_cf_var('ice_v', 'm/s', 'meridional velocity', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_VERTEX),&
      &          ldims=(/nproma,nblks_v/), lrestart_cont=.TRUE.)
      
    CALL add_var(ocean_restart_list, 'ice_u', p_ice%u ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('ice_u', 'm/s', 'zonal velocity', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_diag"),&
      &          lrestart_cont=.FALSE.)
    CALL add_var(ocean_restart_list, 'ice_v', p_ice%v ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('ice_v', 'm/s', 'meridional velocity', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_diag"),&
      &          lrestart_cont=.FALSE.)

    CALL add_var(ocean_restart_list, 'ice_vn', p_ice%vn_e ,&
      &          GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, &
      &          t_cf_var('ice_vn', 'm/s', 'zonal velocity', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE),&
      &          ldims=(/nproma,nblks_e/),&
      &          lrestart_cont=.FALSE.)

    CALL add_var(ocean_restart_list, 'concSum', p_ice%concSum ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('concSum', '', 'total ice concentration', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'newice', p_ice%newice ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('newice', 'm', 'new ice groth in open water', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'zUnderIce', p_ice%zUnderIce ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('zUnderIce', 'm', 'water in upper ocean grid cell below ice', &
      &                   DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
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

    ! add accumulated fields
    CALL add_var(ocean_default_list, 'hi_acc', p_ice%acc%hi ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hi_acc', 'm', 'ice thickness', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_default"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_default_list, 'hs_acc', p_ice%acc%hs ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('hs_acc', 'm', 'snow thickness', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_default"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_default_list, 'conc_acc', p_ice%acc%conc ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('conc_acc', '', 'ice concentration in each ice class', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_default"),&
      &          lrestart_cont=.TRUE.)
    CALL add_var(ocean_default_list, 'ice_u_acc', p_ice%acc%u ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('ice_u_acc', 'm/s', 'zonal velocity', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_default"),&
      &          lrestart_cont=.FALSE.)
    CALL add_var(ocean_default_list, 'ice_v_acc', p_ice%acc%v ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('ice_v_acc', 'm/s', 'meridional velocity', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("ice_default"),&
      &          lrestart_cont=.FALSE.)

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE construct_sea_ice
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
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for tafo failed')
    END IF
    ALLOCATE(p_as%ftdew(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for ftdew failed')
    END IF
    ALLOCATE(p_as%fclou(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for fclou failed')
    END IF

    ALLOCATE(p_as%fu10(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for fu10 failed')
    END IF

    ALLOCATE(p_as%fswr(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for fswr failed')
    END IF

    ALLOCATE(p_as%pao(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for pao failed')
    END IF

    ALLOCATE(p_as%u(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for u failed')
    END IF
    ALLOCATE(p_as%v(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for v failed')
    END IF

    ALLOCATE(p_as%precip(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for precip failed')
    END IF

    ALLOCATE(p_as%evap(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for evap failed')
    END IF

    ALLOCATE(p_as%runoff(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for runoff failed')
    END IF


    p_as%tafo  (:,:) = 0.0_wp
    p_as%ftdew (:,:) = 0.0_wp
    p_as%fclou (:,:) = 0.0_wp
    p_as%fu10  (:,:) = 0.0_wp
    p_as%fswr  (:,:) = 0.0_wp
    p_as%pao   (:,:) = 0.0_wp
    p_as%u     (:,:) = 0.0_wp
    p_as%v     (:,:) = 0.0_wp
    p_as%precip(:,:) = 0.0_wp
    p_as%evap  (:,:) = 0.0_wp
    p_as%runoff(:,:) = 0.0_wp

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
  SUBROUTINE construct_atmos_fluxes(p_patch, p_atm_f, i_no_ice_thick_class)
    !
    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_atmos_fluxes ), INTENT(INOUT) :: p_atm_f
    INTEGER,               INTENT(IN)    :: i_no_ice_thick_class
    ! Local variables
    INTEGER :: ibits = DATATYPE_PACK16
    INTEGER :: alloc_cell_blocks, ist

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_atmos_fluxes'

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    alloc_cell_blocks = p_patch%alloc_cell_blocks

    ALLOCATE(p_atm_f%sens(nproma,i_no_ice_thick_class,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for sens failed')
    END IF

    ALLOCATE(p_atm_f%lat(nproma,i_no_ice_thick_class,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for lat failed')
    END IF

    ALLOCATE(p_atm_f%LWout(nproma,i_no_ice_thick_class,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for LWout failed')
    END IF

    ALLOCATE(p_atm_f%LWnet(nproma,i_no_ice_thick_class,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for LWnet failed')
    END IF

    ALLOCATE(p_atm_f%SWnet(nproma,i_no_ice_thick_class,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for SWnet failed')
    END IF

    ALLOCATE(p_atm_f%bot(nproma,i_no_ice_thick_class,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for sens failed')
    END IF

    ALLOCATE(p_atm_f%dsensdT(nproma,i_no_ice_thick_class,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for dsensdT failed')
    END IF

    ALLOCATE(p_atm_f%dlatdT(nproma,i_no_ice_thick_class,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for sens failed')
    END IF

    ALLOCATE(p_atm_f%dLWdT(nproma,i_no_ice_thick_class,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for dLWdT failed')
    END IF

    ALLOCATE(p_atm_f%stress_x(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for stress_x failed')
    END IF

    ALLOCATE(p_atm_f%stress_y(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for stress_y failed')
    END IF

    ALLOCATE(p_atm_f%rprecw(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for rprecw failed')
    END IF

    ALLOCATE(p_atm_f%rpreci(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for rpreci failed')
    END IF

    ALLOCATE(p_atm_f%sensw(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for sensw failed')
    END IF

    ALLOCATE(p_atm_f%latw(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for latw failed')
    END IF


    ALLOCATE(p_atm_f%LWoutw(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for LWoutw failed')
    END IF

    ALLOCATE(p_atm_f%LWnetw(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for LWnetw failed')
    END IF

    ALLOCATE(p_atm_f%SWnetw(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for SWnetw failed')
    END IF

    ALLOCATE(p_atm_f%LWin(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for LWin failed')
     END IF

    ALLOCATE(p_atm_f%stress_xw(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for stress_xw failed')
     END IF

    ALLOCATE(p_atm_f%stress_yw(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for stress_yw failed')
     END IF

    !albedos need to go into the restart
    CALL add_var(ocean_restart_list, 'albvisdirw', p_atm_f%albvisdirw ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('albvisdirw', '', 'albvisdirw', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albvisdifw', p_atm_f%albvisdifw ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('albvisdifw', '', 'albvisdifw', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albnirdirw', p_atm_f%albnirdirw ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('albnirdirw', '', 'albnirdirw', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albnirdifw', p_atm_f%albnirdifw ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('albnirdifw', '', 'albnirdifw', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albvisdir', p_atm_f%albvisdir ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('albvisdir', '', 'albvisdir', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albvisdif', p_atm_f%albvisdif ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('albvisdif', '', 'albvisdif', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albnirdir', p_atm_f%albnirdir ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('albnirdir', '', 'albnirdir', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)

    CALL add_var(ocean_restart_list, 'albnirdif', p_atm_f%albnirdif ,&
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, &
      &          t_cf_var('albnirdif', '', 'albnirdif', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,i_no_ice_thick_class,alloc_cell_blocks/),in_group=groups("ice_diag"),&
      &          lrestart_cont=.TRUE.)


    ! Initialise everything with zero
    p_atm_f%sens   (:,:,:) = 0.0_wp
    p_atm_f%lat    (:,:,:) = 0.0_wp
    p_atm_f%LWout  (:,:,:) = 0.0_wp
    p_atm_f%LWnet  (:,:,:) = 0.0_wp
    p_atm_f%bot    (:,:,:) = 0.0_wp
    p_atm_f%dsensdT(:,:,:) = 0.0_wp
    p_atm_f%dlatdT (:,:,:) = 0.0_wp
    p_atm_f%dLWdT  (:,:,:) = 0.0_wp
    p_atm_f%rprecw (:,:)   = 0.0_wp
    p_atm_f%rpreci (:,:)   = 0.0_wp
    p_atm_f%sensw  (:,:)   = 0.0_wp
    p_atm_f%latw   (:,:)   = 0.0_wp
    p_atm_f%LWoutw (:,:)   = 0.0_wp
    p_atm_f%LWnetw (:,:)   = 0.0_wp
    p_atm_f%SWnetw (:,:)   = 0.0_wp
    p_atm_f%SWnet  (:,:,:) = 0.0_wp
    p_atm_f%LWin   (:,:)   = 0.0_wp
    p_atm_f%counter        = 0
    ! Initialise the albedos sensibly
    p_atm_f%albvisdir (:,:,:) = albi
    p_atm_f%albvisdif (:,:,:) = albi
    p_atm_f%albnirdir (:,:,:) = albi
    p_atm_f%albnirdif (:,:,:) = albi
    p_atm_f%albvisdirw(:,:) = albedoW_sim
    p_atm_f%albvisdifw(:,:) = albedoW_sim
    p_atm_f%albnirdirw(:,:) = albedoW_sim
    p_atm_f%albnirdifw(:,:) = albedoW_sim

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE construct_atmos_fluxes
  !-------------------------------------------------------------------------
  !
  !>
  !! Destructor of atmos fluxes for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011)
  !
  SUBROUTINE destruct_atmos_fluxes(p_atm_f)
    !
    TYPE(t_atmos_fluxes )       :: p_atm_f
    ! Local variables
    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_atmos_fluxes'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )


    DEALLOCATE(p_atm_f%sens, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for sens failed')
    END IF

    DEALLOCATE(p_atm_f%lat, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for lat failed')
    END IF

    DEALLOCATE(p_atm_f%LWout, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for LWout failed')
    END IF

    DEALLOCATE(p_atm_f%LWnet, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for LWnet failed')
    END IF

    DEALLOCATE(p_atm_f%SWnet, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for SWnet failed')
    END IF

    DEALLOCATE(p_atm_f%bot, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for sens failed')
    END IF

    DEALLOCATE(p_atm_f%dsensdT, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for dsensdT failed')
    END IF

    DEALLOCATE(p_atm_f%dlatdT, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for sens failed')
    END IF

    DEALLOCATE(p_atm_f%dLWdT, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for dLWdT failed')
    END IF

    DEALLOCATE(p_atm_f%rprecw, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for rprecw failed')
    END IF

    DEALLOCATE(p_atm_f%rpreci, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for rpreci failed')
    END IF

    DEALLOCATE(p_atm_f%sensw, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for sensw failed')
    END IF

    DEALLOCATE(p_atm_f%latw, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for latw failed')
    END IF


    DEALLOCATE(p_atm_f%LWoutw, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for LWoutw failed')
    END IF

    DEALLOCATE(p_atm_f%LWnetw, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for LWnetw failed')
    END IF

    DEALLOCATE(p_atm_f%SWnetw, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for SWnetw failed')
    END IF

    DEALLOCATE(p_atm_f%LWin, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for LWin failed')
    END IF

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE destruct_atmos_fluxes

  !-------------------------------------------------------------------------
  !
  !> ice_init
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE ice_init( p_patch_3D, p_os, ice ) !, Qatm, QatmAve)
    TYPE(t_patch_3D), TARGET, INTENT(in)  :: p_patch_3D
    TYPE(t_hydro_ocean_state)             :: p_os
    TYPE (t_sea_ice),      INTENT (INOUT) :: ice
    !TYPE (t_atmos_fluxes), INTENT (INOUT) :: Qatm
    !TYPE (t_atmos_fluxes), INTENT (INOUT) :: QatmAve

    !local variables
    REAL(wp), DIMENSION(nproma,ice%kice, p_patch_3D%p_patch_2D(n_dom)%alloc_cell_blocks) :: &
      & Tinterface, & ! temperature at snow-ice interface
      & draft,      & ! position of ice-ocean interface below sea level
      & Tfw           ! Ocean freezing temperature [C]

    TYPE(t_patch), POINTER                :: p_patch

    !INTEGER i,j,k      ! counter for loops
    INTEGER k !, jb, jc, i_startidx_c, i_endidx_c! counter for loops
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:ice_init'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    p_patch => p_patch_3D%p_patch_2D(n_dom)

    !Constructor basic init already done at this point
    !   CALL alloc_mem_commo_ice (ice, Qatm, QatmAve)
    !   CALL ice_zero            (ice, Qatm, QatmAve)

    ! FORALL(i=1:nproma, j=1:p_patch%alloc_cell_blocks, k=1:ice%kice)
    !    ice% hi    (i,j,k) = sictho (i,j)
    !    ice% hs    (i,j,k) = sicsno (i,j)
    ! END FORALL

    ! LL Note: this needs to be rewritten using subsets !
    IF ( no_tracer >= 2 ) THEN
      DO k=1,ice%kice
        Tfw(:,k,:) = -mu*p_os%p_prog(nold(1))%tracer(:,1,:,2)
      ENDDO
    ELSE
      Tfw(:,:,:) = Tf
    ENDIF

    ice% Tsurf(:,:,:)  = Tf
    ice% T1   (:,:,:)  = Tf
    ice% T2   (:,:,:)  = Tf
    ice% conc (:,:,:)  = 0.0_wp
    draft     (:,:,:)  = 0.0_wp

    ! Stupid initialisation trick for Levitus initialisation
    IF (use_file_initialConditions) THEN
      WHERE (p_os%p_prog(nold(1))%tracer(:,1,:,1) <= -1.6_wp &
          &     .and. v_base%lsm_c(:,1,:) <= sea_boundary )
        ice%hi(:,1,:) = 2.0_wp
        ice%hs(:,1,:) = 0.2_wp
        ice%conc(:,1,:) = 0.95_wp
      ENDWHERE
!      IF ( no_tracer < 2 ) THEN
!        WHERE (p_os%p_prog(nold(1))%tracer(:,:,:,1) <= -1.0_wp    &
!          &     .and. v_base%lsm_c(:,:,:) <= sea_boundary )   &
!          &             p_os%p_prog(nold(1))%tracer(:,:,:,1) = Tf
!      ENDIF
    ENDIF

    WHERE(ice% hi(:,:,:) > 0.0_wp)
      ice% Tsurf (:,:,:) = Tfw(:,:,:)
      ice% T1    (:,:,:) = Tfw(:,:,:)
      ice% T2    (:,:,:) = Tfw(:,:,:)
      Tinterface (:,:,:) = (Tfw(:,:,:) * (ki/ks * ice%hs(:,:,:)/ice%hi(:,:,:))+&
        &                    ice%Tsurf(:,:,:)) / (1.0_wp+ki/ks * ice%hs(:,:,:)/ice%hi(:,:,:))
      ice% conc  (:,:,:) = 1.0_wp/REAL(ice%kice,wp)
      ice% T1    (:,:,:) = Tfw(:,:,:) + 2._wp/3._wp*(Tinterface(:,:,:)-Tfw(:,:,:))
      ice% T2    (:,:,:) = Tfw(:,:,:) + 1._wp/3._wp*(Tinterface(:,:,:)-Tfw(:,:,:))
      draft      (:,:,:) = (rhos * ice%hs(:,:,:) + rhoi * ice%hi(:,:,:)) / rho_ref
    END WHERE

    ! TODO: use prism_thick_flat_sfc_c instead of del_zlev_m
    ice%zUnderIce (:,:)   = v_base%del_zlev_m(1) +  p_os%p_prog(nold(1))%h(:,:) &
      &                      - sum(draft(:,:,:) * ice%conc(:,:,:),2)

    IF ( i_ice_dyn == 1 ) THEN ! AWI dynamics
      CALL ice_fem_grid_init(p_patch_3D)
      CALL ice_init_fem
      CALL ice_fem_grid_post(p_patch)
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('IceInit: hi       ' ,ice%hi       ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: conc     ' ,ice%conc     ,str_module, idt_src, in_subset=p_patch%cells%owned)
    idt_src=4  ! output print level (1-5, fix)        
    CALL dbg_print('IceInit: Tfw      ' ,Tfw          ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: draft    ' ,draft        ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceInit: zUnderIce' ,ice%zUnderIce,str_module, idt_src, in_subset=p_patch%cells%owned)
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
  SUBROUTINE ice_slow(p_patch_3D, p_os, p_as, ice, QatmAve, p_sfc_flx, p_op_coeff)
    TYPE(t_patch_3D), TARGET, INTENT(in) :: p_patch_3D
    !TYPE(t_patch),            INTENT(IN)     :: p_patch 
    TYPE(t_hydro_ocean_state),INTENT(INOUT)  :: p_os
    TYPE(t_atmos_for_ocean),  INTENT(IN)     :: p_as
    TYPE (t_sea_ice),         INTENT (INOUT) :: ice
    !TYPE (t_atmos_fluxes),    INTENT (INOUT) :: Qatm
    TYPE (t_atmos_fluxes),    INTENT (INOUT) :: QatmAve
    TYPE(t_sfc_flx),          INTENT (INOUT) :: p_sfc_flx
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff

    TYPE(t_patch), POINTER :: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: jb, jc, i_startidx_c, i_endidx_c

    !-------------------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_ice_slow)

    p_patch => p_patch_3D%p_patch_2D(n_dom)
    ! subset range pointer
    all_cells => p_patch%cells%all 

    !CALL ave_fluxes     (ice, QatmAve)
    CALL ice_zero       (ice)

    ice%hiold(:,:,:) = ice%hi(:,:,:)
    ice%hsold(:,:,:) = ice%hs(:,:,:)
    CALL dbg_print('IceSlow: hi before groth' ,ice%hi ,str_module,5, in_subset=p_patch%cells%owned)
    ! #achim
    IF      ( i_ice_therm == 2 ) THEN
      CALL ice_growth_winton    (p_patch, p_os, ice, QatmAve%rpreci)!, QatmAve%lat)
    ELSE IF ( i_ice_therm == 1 .OR. i_ice_therm == 3 ) THEN !2=zerolayer, 3=simple fluxes from dirk's thesis
      CALL ice_growth_zerolayer (p_patch, p_os, ice, QatmAve%rpreci)
    END IF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('IceSlow: hi after growth'       ,ice%hi   ,str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: Conc. after growth'    ,ice%conc ,str_module, idt_src, in_subset=p_patch%cells%owned)
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('IceSlow: p_ice%u bef. dyn'    ,ice%u_prog ,str_module, idt_src, in_subset=p_patch%verts%owned)
    CALL dbg_print('IceSlow: p_ice%v bef. dyn'    ,ice%v_prog ,str_module, idt_src, in_subset=p_patch%verts%owned)
    !---------------------------------------------------------------------

    CALL upper_ocean_TS (p_patch,p_os,ice, QatmAve, p_sfc_flx)
    CALL ice_conc_change(p_patch,ice, p_os,p_sfc_flx)

    CALL ice_ocean_stress( p_patch, QatmAve, p_sfc_flx, ice, p_os )

    IF ( i_ice_dyn >= 1 ) THEN
      ! AWI FEM model wrapper
      CALL fem_ice_wrap ( p_patch_3D, ice, p_os, QatmAve, p_op_coeff )
      CALL ice_advection( p_patch_3D, p_op_coeff, ice )
    ELSE
      ice%u = 0._wp
      ice%v = 0._wp
    ENDIF

    CALL ice_clean_up( p_patch_3D, ice, p_sfc_flx, p_os )

    !CALL ice_advection  (ice)
    !CALL write_ice      (ice,QatmAve,1,ie,je)
    !CALL ice_zero       (ice, QatmAve)
    !sictho = ice%hi   (:,:,1) * ice%conc (:,:,1)
    !sicomo = ice%conc (:,:,1)
    !sicsno = ice%hs   (:,:,1) * ice%conc (:,:,1)
    CALL dbg_print('IceSlow: p_ice%u'            ,ice%u_prog,             str_module,3, in_subset=p_patch%verts%owned)
    CALL dbg_print('IceSlow: p_ice%v'            ,ice%v_prog,             str_module,3, in_subset=p_patch%verts%owned)
    CALL dbg_print('IceSlow: hi endOf slow'      ,ice%hi,                 str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: Conc.  EndOf slow'  ,ice%conc,               str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: ConcSumEndOf slow',  ice%concSum,            str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: p_os%prog(nold)%vn' ,p_os%p_prog(nold(1))%vn,str_module,5, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: p_os%prog(nnew)%vn' ,p_os%p_prog(nnew(1))%vn,str_module,5, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: p_os%diag%u'        ,p_os%p_diag%u,          str_module,4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: p_os%diag%v'        ,p_os%p_diag%v,          str_module,4, in_subset=p_patch%cells%owned)

    IF (ltimer) CALL timer_stop(timer_ice_slow)

  END SUBROUTINE ice_slow

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! !  ice_clean_up: Fix over and under shoots and beutify output
  !! @par Revision History
  !! Initial release by Einar Olason, MPI-M (2013-10).
  !!
  SUBROUTINE ice_clean_up( p_patch_3D, p_ice, p_sfc_flx, p_os )
    TYPE(t_patch_3D),TARGET,   INTENT(IN)    :: p_patch_3D
    TYPE(t_sea_ice),           INTENT(INOUT) :: p_ice
    TYPE(t_sfc_flx),           INTENT(INOUT) :: p_sfc_flx
    TYPE(t_hydro_ocean_state), INTENT(INOUT) :: p_os

    ! Local variables
    ! ranges
    TYPE(t_subset_range), POINTER :: all_cells
    ! pathc
    TYPE(t_patch),POINTER    :: p_patch
    ! counters
    INTEGER :: k, jb, jc, i_startidx_c, i_endidx_c
    ! Sea surface salinity
    REAL(wp), DIMENSION (nproma, p_patch_3d%p_patch_2D(1)%alloc_cell_blocks) :: sss

    ! subset range pointer
    p_patch => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all 
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

        ! Fix over shoots - ONLY for the one-ice-class case
        IF ( p_ice%conc(jc,1,jb) > 1._wp ) p_ice%conc(jc,1,jb) = 1._wp

        ! Fix under shoots and remove ice where there's almost none left
        DO k = 1, p_ice%kice
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary &
            &   .AND. ( p_ice%vol(jc,k,jb) <= 0._wp .OR. p_ice%conc(jc,k,jb) <= 1e-4_wp ) ) THEN
            ! Volmue flux due to removal
            p_sfc_flx%forc_fw_ice_vol(jc,jb) = p_sfc_flx%forc_fw_ice_vol(jc,jb) &
              & + p_ice%hi(jc,k,jb)*p_ice%conc(jc,k,jb)*rhoi/(rho_ref*dtime)    & ! Ice
              & + p_ice%hs(jc,k,jb)*p_ice%conc(jc,k,jb)*rhos/(rho_ref*dtime)      ! Snow
            ! Tracer flux due to removal
            p_sfc_flx%forc_fw_bc_ice (jc,jb) = p_sfc_flx%forc_fw_bc_ice (jc,jb)                      &
              & + (1._wp-sice/sss(jc,jb))*p_ice%hi(jc,k,jb)*p_ice%conc(jc,k,jb)*rhoi/(rho_ref*dtime) & ! Ice
              & + p_ice%hs(jc,k,jb)*p_ice%conc(jc,k,jb)*rhos/(rho_ref*dtime)                           ! Snow
            ! Heat flux due to removal
            p_sfc_flx%forc_hflx      (jc,jb) = p_sfc_flx%forc_hflx(jc,jb)       &
              & + p_ice%hi(jc,k,jb)*p_ice%conc(jc,k,jb)*alf*rhoi/dtime          & ! Ice
              & + p_ice%hs(jc,k,jb)*p_ice%conc(jc,k,jb)*alf*rhos/dtime            ! Snow
            p_ice%conc(jc,k,jb) = 0._wp
            p_ice%hi  (jc,k,jb) = 0._wp
            p_ice%vol (jc,k,jb) = 0._wp
            p_ice%hs  (jc,k,jb) = 0._wp
            p_ice%vols(jc,k,jb) = 0._wp
          ENDIF
        ENDDO

      ENDDO
    ENDDO

    p_ice%concSum = SUM(p_ice%conc, 2)

  END SUBROUTINE ice_clean_up

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
  SUBROUTINE get_atmos_fluxes (p_patch, p_os,p_as,ice, Qatm)
    TYPE(t_patch),            INTENT(IN)    :: p_patch
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os
    TYPE(t_atmos_for_ocean),  INTENT(IN)    :: p_as
    TYPE (t_sea_ice),         INTENT(INOUT) :: ice
    TYPE (t_atmos_fluxes),    INTENT(INOUT) :: Qatm

!#ifdef coupled
    !Qatm% SWin   =
    !Qatm% LWin   =
    !Qatm% sens   =
    !Qatm% lat    =
    !Qatm% dsensdT =
    !Qatm% dlatdT  =
    !Qatm% dLWdT   =
!#elif defined CORE
    !CALL budget_core   (ice, Qatm)
!#else
    CALL calc_bulk_flux_oce(p_patch, p_as, p_os, Qatm)
    CALL calc_bulk_flux_ice(p_patch, p_as, ice , Qatm)
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
  SUBROUTINE sum_fluxes        (Qatm, QatmAve)
    TYPE (t_atmos_fluxes), INTENT (IN)    :: Qatm
    TYPE (t_atmos_fluxes), INTENT (INOUT) :: QatmAve

    QatmAve % sens   (:,:,:) = QatmAve % sens   (:,:,:) + Qatm % sens   (:,:,:)
    QatmAve % sensw  (:,:)   = QatmAve % sensw  (:,:)   + Qatm % sensw  (:,:)
    QatmAve % lat    (:,:,:) = QatmAve % lat    (:,:,:) + Qatm % lat    (:,:,:)
    QatmAve % latw   (:,:)   = QatmAve % latw   (:,:)   + Qatm % latw   (:,:)
    QatmAve % LWout  (:,:,:) = QatmAve % LWout  (:,:,:) + Qatm % LWout  (:,:,:)
    QatmAve % LWoutw (:,:)   = QatmAve % LWoutw (:,:)   + Qatm % LWoutw (:,:)
    QatmAve % LWnet  (:,:,:) = QatmAve % LWnet  (:,:,:) + Qatm % LWnet  (:,:,:)
    QatmAve % LWnetw (:,:)   = QatmAve % LWnetw (:,:)   + Qatm % LWnetw (:,:)
    QatmAve % SWnet  (:,:,:) = QatmAve % SWnet  (:,:,:) + Qatm % SWnet  (:,:,:)
    QatmAve % SWnetw (:,:)   = QatmAve % SWnetw (:,:)   + Qatm % SWnetw (:,:)
    QatmAve % LWin   (:,:)   = QatmAve % LWin   (:,:)   + Qatm % LWin   (:,:)
    QatmAve % rprecw (:,:)   = QatmAve % rprecw (:,:)   + Qatm % rprecw (:,:)
    QatmAve % rpreci (:,:)   = QatmAve % rpreci (:,:)   + Qatm % rpreci (:,:)
    QatmAve % counter        = QatmAve % counter + 1

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
  SUBROUTINE ave_fluxes (ice, QatmAve)
    TYPE (t_sea_ice),      INTENT (INOUT) :: ice
    TYPE (t_atmos_fluxes), INTENT (INOUT) :: QatmAve
    !
    !Local variables
    REAL(wp) :: ctr

    !-------------------------------------------------------------------------------

    ctr = REAL(QatmAve% counter,wp)
    QatmAve% sens   (:,:,:) = QatmAve% sens  (:,:,:)  / ctr
    QatmAve% sensw  (:,:)   = QatmAve% sensw (:,:)    / ctr
    QatmAve% lat    (:,:,:) = QatmAve% lat   (:,:,:)  / ctr
    QatmAve% latw   (:,:)   = QatmAve% latw  (:,:)    / ctr
    QatmAve% LWout  (:,:,:) = QatmAve% LWout (:,:,:)  / ctr
    QatmAve% LWoutw (:,:)   = QatmAve% LWoutw(:,:)    / ctr
    QatmAve% LWnet  (:,:,:) = QatmAve% LWnet (:,:,:)  / ctr
    QatmAve% LWnetw (:,:)   = QatmAve% LWnetw(:,:)    / ctr
    QatmAve% SWnet  (:,:,:) = QatmAve% SWnet (:,:,:)  / ctr
    QatmAve% SWnetw (:,:)   = QatmAve% SWnetw(:,:)    / ctr
    QatmAve% LWin   (:,:)   = QatmAve% LWin  (:,:)    / ctr
    QatmAve% rprecw (:,:)   = QatmAve% rprecw(:,:)    / ctr
    QatmAve% rpreci (:,:)   = QatmAve% rpreci(:,:)    / ctr
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
    !TYPE (t_atmos_fluxes), INTENT (INOUT) :: Qatm
    !TYPE (t_atmos_fluxes), INTENT (INOUT) :: QatmAve

    !Qatm    % sens        (:,:,:) = 0._wp
    !Qatm    % sensw       (:,:)   = 0._wp
    !Qatm    % lat         (:,:,:) = 0._wp
    !Qatm    % latw        (:,:)   = 0._wp
    !Qatm    % LWout       (:,:,:) = 0._wp
    !Qatm    % LWoutw      (:,:)   = 0._wp
    !Qatm    % LWnet       (:,:,:) = 0._wp
    !Qatm    % LWnetw      (:,:)   = 0._wp
    !Qatm    % SWin        (:,:)   = 0._wp
    !Qatm    % LWin        (:,:)   = 0._wp
    !Qatm    % rprecw      (:,:)   = 0._wp
    !Qatm    % rpreci      (:,:)   = 0._wp

!    QatmAve % sens        (:,:,:) = 0._wp
!    QatmAve % sensw       (:,:)   = 0._wp
!    QatmAve % lat         (:,:,:) = 0._wp
!    QatmAve % latw        (:,:)   = 0._wp
!    QatmAve % LWout       (:,:,:) = 0._wp
!    QatmAve % LWoutw      (:,:)   = 0._wp
!    QatmAve % LWnet       (:,:,:) = 0._wp
!    QatmAve % LWnetw      (:,:)   = 0._wp
!    QatmAve % SWin        (:,:)   = 0._wp
!    QatmAve % LWin        (:,:)   = 0._wp
!    QatmAve % rprecw      (:,:)   = 0._wp
!    QatmAve % rpreci      (:,:)   = 0._wp
!    QatmAve % counter             = 0

!    ice     % Qbot        (:,:,:) = 0._wp
!    ice     % Qtop        (:,:,:) = 0._wp
    ice     % surfmelt    (:,:,:) = 0._wp
    ice     % surfmeltT   (:,:,:) = 0._wp
    ice     % evapwi      (:,:,:) = 0._wp
    ice     % hiold       (:,:,:) = 0._wp
    ice     % hsold       (:,:,:) = 0._wp
    ice     % snow_to_ice (:,:,:) = 0._wp
    ice     % heatOceI    (:,:,:) = 0._wp

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
  SUBROUTINE upper_ocean_TS(p_patch, p_os,ice, QatmAve, p_sfc_flx)
    TYPE(t_patch),             INTENT(IN)    :: p_patch
    TYPE(t_hydro_ocean_state), INTENT(INOUT) :: p_os
    !TYPE(t_atmos_for_ocean),   INTENT(IN)    :: p_as
    TYPE(t_sea_ice),           INTENT(INOUT) :: ice
    TYPE(t_atmos_fluxes),      INTENT(INOUT) :: QatmAve
    TYPE(t_sfc_flx),           INTENT(INOUT) :: p_sfc_flx

    !Local Variables
    ! position of ice-ocean interface below sea level                       [m]
    REAL(wp) :: draft(nproma,ice%kice, p_patch%alloc_cell_blocks)

    REAL(wp), DIMENSION (nproma, p_patch%alloc_cell_blocks) ::   &
      & draftAve,      &! average draft of sea ice within a grid cell             [m]
      & zUnderIceOld,  &! water in upper ocean grid cell below ice (prev. time)   [m]
      & heatOceI,      &! heat flux into ocean through formerly ice covered areas [W/m^2]
      & heatOceW,      &! heat flux into ocean through open water areas           [W/m^2]
      & delHice,       &! average change in ice thickness within a grid cell      [m]
      & delHsnow,      &! average change in snow thickness within a grid cell     [m]
      & snowiceave,    &! average snow to ice conversion within a grid cell       [m]
      & Tfw,           &! sea surface freezing temperature                        [C]
      & sst,           &! sea surface temperature - approx. after cooling         [C]
      & sss,           &! sea surface salinity                                    [psu]
      & preci,         &! solid precipitation rate                                [m/s]
      & precw           ! liquid precipitation rate                               [m/s]
      !& evap,          &! evaporated water                                       [psu]

    ! Needs work with FB_BGC_OCE etc.
    !REAL(wp)         :: swsum
    !REAL(wp),POINTER :: sao_top(:,:)
    !-------------------------------------------------------------------------------

    ! #eoo# What is swsum?
    ! swsum = 0.0_wp
    !sao_top =>p_os%p_prog(nold(1))%tracer(:,1,:,2)

    ! Ocean points only
    ! Calculate change in water level 'zo' from liquid and solid precipitation and
    ! evaporation
    sss             (:,:)   = p_os%p_prog(nold(1))%tracer(:,1,:,2)
    precw           (:,:)   = QatmAve% rprecw (:,:)
    preci           (:,:)   = QatmAve% rpreci (:,:)
    !evap            (:,:)   = (QatmAve% latw(:,:)/ alv * dtime * &
    !  &                       sum(ice%conc(:,:,:), 2) +          &
    !  &                       sum(ice%evapwi(:,:,:) * ice% conc(:,:,:), 2)) /rho_ref
    
    ! Calculate the sea surface freezing temperature                        [C]
    if ( no_tracer >= 2 ) then
      Tfw(:,:) = -mu*sss(:,:)
    else
      Tfw(:,:) = Tf
    endif

    ! TODO: No temperature change due to precip yet. Should not be done here

    ! Calculate average draft and thickness of water underneath ice in upper ocean
    ! grid box
    zUnderIceOld    (:,:)   = ice%zUnderIce(:,:)
    draft           (:,:,:) = (rhos * ice%hs(:,:,:) + rhoi * ice%hi(:,:,:)) / rho_ref
    draftave        (:,:)   = sum(draft(:,:,:) * ice%conc(:,:,:),2)
    ice%zUnderIce   (:,:)   = v_base%del_zlev_m(1) + p_os%p_prog(nold(1))%h(:,:) - draftave(:,:)

    ! Calculate average change in ice thickness and the snow-to-ice conversion
    Delhice   (:,:) = SUM( ( ice%hi(:,:,:) - ice%hiold(:,:,:) )*ice%conc(:,:,:), 2 )
    Delhsnow  (:,:) = SUM( ( ice%hs(:,:,:) - ice%hsold(:,:,:) )*ice%conc(:,:,:), 2 )
    snowiceave(:,:) = SUM( ice%snow_to_ice(:,:,:)*ice% conc(:,:,:), 2 )

    ! Calculate heat input through formerly ice covered and through open water areas
    heatOceI(:,:)   = sum(ice% heatOceI(:,:,:) * ice% conc(:,:,:),2)
    heatOceW(:,:) = ( QatmAve%SWnetw(:,:)                       &
      &         + QatmAve%LWnetw(:,:) + QatmAve%sensw(:,:)+     &
      &                 QatmAve%latw(:,:) )*(1.0_wp-sum(ice%conc(:,:,:),2))

    ! Calculate possible super-cooling of the surface layer
    sst = p_os%p_prog(nold(1))%tracer(:,1,:,1) +        &
      &      dtime*heatOceW(:,:)/( clw*rho_ref*ice%zUnderIce(:,:) )

    ! Add energy for new-ice formation due to supercooled ocean to  ocean temperature, form new ice
    WHERE ( sst < Tfw(:,:) .AND. v_base%lsm_c(:,1,:) <= sea_boundary )
      ! New ice forming over open water due to super cooling
      ice%newice(:,:) = ( Tfw(:,:) - sst(:,:) )*ice%zUnderIce(:,:)*clw*rho_ref/( alf*rhoi )
      ! Flux required to cool the ocean to the freezing point
      heatOceW(:,:)   = ( Tfw(:,:) - p_os%p_prog(nold(1))%tracer(:,1,:,1) )     &
        &     *ice%zUnderIce(:,:)*(1.0_wp-ice%concSum(:,:))*clw*rho_ref/dtime
    ENDWHERE

    CALL dbg_print('UpperOceTS: Delhice  ', Delhice      ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: Delhsnow ', Delhsnow     ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: draft    ', draft        ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: zUnderIce', ice%zUnderIce,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: newice   ', ice%newice   ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpperOceTS: heatOceW ', heatOceW     ,str_module, 4, in_subset=p_patch%cells%owned)

    ! Diagnosis: collect the 4 parts of heat fluxes into the p_sfc_flx variables - no flux under ice:
    p_sfc_flx%forc_swflx(:,:) = QatmAve%SWnetw(:,:)*(1.0_wp-sum(ice%conc(:,:,:),2))
    p_sfc_flx%forc_lwflx(:,:) = QatmAve%LWnetw(:,:)*(1.0_wp-sum(ice%conc(:,:,:),2))
    p_sfc_flx%forc_ssflx(:,:) = QatmAve%sensw (:,:)*(1.0_wp-sum(ice%conc(:,:,:),2))
    p_sfc_flx%forc_slflx(:,:) = QatmAve%latw  (:,:)*(1.0_wp-sum(ice%conc(:,:,:),2))

    ! #slo# 2013-06
    ! Change of upper ocean temperature according to heat fluxes is done in vertical diffusion equation
    !  - forc_hflx is calculated here, provided to tracer eq. via forc_tracer(1), calculated in update_sfcflx
    !  - forc_tracer(1) is 
    !p_os%p_prog(nold(1))%tracer(:,1,:,1) = p_os%p_prog(nold(1))%tracer(:,1,:,1)&
    !  &                                    + dtime*(heatOceI + heatOceW) /               &
    !  &                                    (clw*rho_ref * ice%zUnderIce)
    ! TODO: should we also divide with ice%zUnderIce / ( v_base%del_zlev_m(1) +  p_os%p_prog(nold(1))%h(:,:) ) ?
    !p_sfc_flx%forc_tracer(:,:,1) = (heatOceI + heatOceW) / (clw*rho_ref)
    p_sfc_flx%forc_hflx(:,:) = heatOceI(:,:) + heatOceW(:,:)

    ! TODO:
    ! Temperature change of upper ocean grid cell due  to melt-water inflow and
    ! precipitation
    !p_os%p_prog(nold(1))%tracer(:,1,:,1) = (p_os%p_prog(nold(1))%tracer(:,1,:,1) &
    !  &                      *zUnderIceOld                                       &
    !  &                      + precw*p_as%tafo + preci*0.0_wp + &                             !!!!!!!!!Dirk: times 0.0 ????
    !  &                        sum(ice%surfmeltT(:,:,:) * ice%surfmelt * ice%conc(:,:,:),2)) / &
    !  &                        (zUnderIceOld + sum(ice%surfmelt*ice%conc(:,:,:),2) +    &
    !  &                        precw + preci)
    !
    ! Change salinity of upper ocean grid box from ice growth/melt, snowice
    ! formation and precipitation
    !p_os%p_prog(nold(1))%tracer(:,1,:,2) = p_os%p_prog(nold(1))%tracer(:,1,:,2)  &
    !  &                                    + (Delhice(:,:)*rhoi - snowiceave(:,:)*rhos)/rho_ref *  &
    !  &                                    MIN(Sice, sao_top(:,:)) / ice%zUnderIce(:,:)

    ! #slo# 2013-06
    ! Change in salinity is calculated according to resulting freshwater flux due to sea ice change:
    !  - fw_ice_impl is flux in m/s >0 for Delhice<0, i.e. positive input of water = decrease of sea ice depth
    !p_sfc_flx%forc_fwsice(:,:) = -Delhice(:,:)*rhoi - snowiceave(:,:)*rhos)/(rho_ref*dtime)

    ! Volmue flux
    p_sfc_flx%forc_fw_ice_vol (:,:) = -Delhice(:,:)* rhoi/(rho_ref*dtime)   & ! Ice melt
      &                               -Delhsnow(:,:)*rhos/(rho_ref*dtime)   & ! Snow melt
      &                              + precw(:,:)*ice%concSum(:,:)          & ! Rain goes through
      &       - (1._wp-ice%concSum(:,:))*ice%newice(:,:)*rhoi/(rho_ref*dtime) ! New-ice formation

    ! Tracer flux
    WHERE (v_base%lsm_c(:,1,:) <= sea_boundary )
      p_sfc_flx%forc_fw_bc_ice (:,:) = precw(:,:)*ice%concSum(:,:)           & ! Rain goes through
        &       - (1._wp-sice/sss(:,:))*Delhice(:,:)*rhoi/(rho_ref*dtime)    & ! Ice melt
        &       - Delhsnow(:,:)*rhos/(rho_ref*dtime)                         & ! Snow melt
        &       - (1._wp-sice/sss(:,:))*(1._wp-ice%concSum(:,:))             & ! New-ice formation
        &         *ice%newice(:,:)*rhoi/(rho_ref*dtime)                       
    ENDWHERE

    !heatabs         (:,:)   = swsum * QatmAve% SWin(:,:) * (1 - ice%concsum)

    ! set to zero on land points
    WHERE (v_base%lsm_c(:,1,:) > sea_boundary )
      p_sfc_flx%forc_hflx (:,:) = 0.0_wp
      p_sfc_flx%forc_swflx(:,:) = 0.0_wp
      p_sfc_flx%forc_lwflx(:,:) = 0.0_wp
      p_sfc_flx%forc_ssflx(:,:) = 0.0_wp
      p_sfc_flx%forc_slflx(:,:) = 0.0_wp
    END WHERE

    CALL dbg_print('UpperOceTS: FwBcIce  ', p_sfc_flx%forc_fw_bc_ice, str_module, 4, in_subset=p_patch%cells%owned)

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
  SUBROUTINE ice_conc_change(p_patch,ice, p_os,p_sfc_flx)

    TYPE(t_patch),             INTENT(IN)    :: p_patch
    TYPE (t_sea_ice),          INTENT(INOUT) :: ice
    TYPE(t_hydro_ocean_state), INTENT(INOUT) :: p_os
    !TYPE (t_atmos_fluxes),     INTENT(IN)    :: QatmAve
    TYPE(t_sfc_flx),           INTENT(INOUT) :: p_sfc_flx

    INTEGER  :: k
    REAL(wp) :: sst(nproma,p_patch%alloc_cell_blocks)
    REAL(wp) :: sss(nproma,p_patch%alloc_cell_blocks)
    REAL(wp) :: Tfw(nproma,p_patch%alloc_cell_blocks) ! Ocean freezing temperature [C]

    CALL dbg_print('IceConcCh: IceConc beg' ,ice%conc, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: vol  at beg' ,ice%vol , str_module, 4, in_subset=p_patch%cells%owned)

    if ( no_tracer >= 2 ) then
      Tfw(:,:) = -mu*p_os%p_prog(nold(1))%tracer(:,1,:,2)
    else
      Tfw(:,:) = Tf
    endif

    ! This should not be needed
    ! TODO ram - remove all instances of p_patch%cells%area(:,:) and test
    ! See also dynamics_fem/mo_ice_fem_utils.f90
    DO k=1,ice%kice
      ice%vol (:,k,:) = ice%hi(:,k,:)*ice%conc(:,k,:)*p_patch%cells%area(:,:)
      ice%vols(:,k,:) = ice%hs(:,k,:)*ice%conc(:,k,:)*p_patch%cells%area(:,:)
    ENDDO

    ! Concentration change due to new ice formation
    WHERE ( ice%newice(:,:) > 0._wp .AND. v_base%lsm_c(:,1,:) <= sea_boundary )
      ! New volume - we just preserve volume:
      ice%vol  (:,1,:) = ice%vol(:,1,:)         &
        &       + ( 1._wp-ice%conc(:,1,:) )*ice%newice(:,:)*p_patch%cells%area(:,:)

      ! Hibler's way to change the concentration 
      !  - the formulation here uses the default values of leadclose parameters 2 and 3 in MPIOM:
      !    1 and 0 respectively
      ice%conc (:,1,:) = min( 1._wp,    &
        &               ice%conc(:,1,:) + ice%newice(:,:)*( 1._wp-ice%conc(:,1,:) )/hnull )

      ! New ice and snow thickness
      ice%hi   (:,1,:) = ice%vol (:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
      ice%hs   (:,1,:) = ice%vols(:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
      !TODO: Re-calculate temperatures to conserve energy when we change the ice thickness
    ENDWHERE

    CALL dbg_print('IceConcCh: conc leadcl' ,ice%conc, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hi   leadcl' ,ice%hi  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hs   leadcl' ,ice%hs  , str_module, 4, in_subset=p_patch%cells%owned)

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

    CALL dbg_print('IceConcCh: conc latMlt' ,ice%conc, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hi   latMlt' ,ice%hi  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hs   latMlt' ,ice%hs  , str_module, 4, in_subset=p_patch%cells%owned)

    ! Ice cannot grow thinner than hmin
    WHERE ( ice%hi(:,1,:) > 0._wp )
      ice%hi  (:,1,:) = MAX( hmin, ice%hi(:,1,:) )
      ice%conc(:,1,:) = ice%vol(:,1,:) / ( ice%hi(:,1,:)*p_patch%cells%area(:,:) )
    ENDWHERE

    ice%concSum(:,:)  = SUM(ice%conc(:,:,:),2)

    WHERE (ice%hi(:,1,:) <= 0._wp)
      ice%Tsurf(:,1,:) = Tfw(:,:)
      ice%T1   (:,1,:) = Tfw(:,:)
      ice%T2   (:,1,:) = Tfw(:,:)
      ice%conc (:,1,:) = 0.0_wp
      ice%hi   (:,1,:) = 0.0_wp
      ice%E1   (:,1,:) = 0.0_wp
      ice%E2   (:,1,:) = 0.0_wp
      ice%vol  (:,1,:) = 0.0_wp
    ENDWHERE

    CALL dbg_print('IceConcCh: IceConc end' ,ice%conc, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hi   at end' ,ice%hi  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hs   at end' ,ice%hs  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: vol  at end' ,ice%vol , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: vols at end' ,ice%vols, str_module, 4, in_subset=p_patch%cells%owned)

  END SUBROUTINE ice_conc_change


  !-------------------------------------------------------------------------
  !
  !> Forcing_from_bulk equals sbr "Budget_omip" in MPIOM.
  !! Sets the atmospheric fluxes for the update of the ice
  !! temperature and ice growth rates for OMIP forcing
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-07). Originally code written by
  !! Dirk Notz, following MPIOM. Code transfered to ICON.
  !! Einar Olason, split calc_atm_fluxes_from_bulk into calc_bulk_flux_ice and calc_bulk_flux_oce
  !! so that the ocean model can be run without the ice model, but with OMIP fluxes.
  !
  SUBROUTINE calc_bulk_flux_ice(p_patch, p_as, p_ice, Qatm, datetime)
    TYPE(t_patch),            INTENT(IN), TARGET    :: p_patch
    TYPE(t_atmos_for_ocean),  INTENT(IN)    :: p_as
    TYPE(t_sea_ice),          INTENT(IN)    :: p_ice
    TYPE(t_atmos_fluxes),     INTENT(INOUT) :: Qatm

    TYPE(t_datetime), OPTIONAL, INTENT(IN)   :: datetime


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
    REAL(wp) :: ramp

    TYPE(t_subset_range), POINTER :: all_cells

    !CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:calc_bulk_flux_ice'
    !-------------------------------------------------------------------------
    !CALL message(TRIM(routine), 'start' )

    tafoK(:,:)  = p_as%tafo(:,:)  + tmelt               ! Change units of tafo  to Kelvin
    ftdewC(:,:) = p_as%ftdew(:,:) - tmelt                    ! Change units of ftdew to C

    ! subset range pointer
    all_cells => p_patch%cells%all

    !-----------------------------------------------------------------------
    ! Compute water vapor pressure and specific humididty in 2m height (esta)
    ! and at water surface (estw) according to "Buck Research Manual (1996)
    ! (see manuals for instruments at http://www.buck-research.com/);
    ! updated from Buck, A. L., New equations for computing vapor pressure and
    ! enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981"
    !-----------------------------------------------------------------------

    aw=611.21_wp; bw=18.729_wp; cw=257.87_wp; dw=227.3_wp
    ai=611.15_wp; bi=23.036_wp; ci=279.82_wp; di=333.7_wp

    AAw=7.2e-4_wp; BBw=3.20e-6_wp; CCw=5.9e-10_wp
    AAi=2.2e-4_wp; BBi=3.83e-6_wp; CCi=6.4e-10_wp

    alpha=0.62197_wp; beta=0.37803_wp

    fa   = 1.0_wp+AAw+p_as%pao*(BBw+CCw*ftdewC**2)
    esta = fa * aw*EXP((bw-ftdewC/dw)*ftdewC/(ftdewC+cw))

    sphumida  = alpha * esta/(p_as%pao-beta*esta)
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
    Qatm%LWin(:,:) = 0._wp
    Qatm%LWout(:,:,:) = 0._wp

    !-----------------------------------------------------------------------
    !  Calculate bulk equations according to
    !      Kara, B. A., P. A. Rochford, and H. E. Hurlburt, 2002:
    !      Air-Sea Flux Estimates And The 19971998 Enso Event,  Bound.-Lay.
    !      Met., 103(3), 439-458, doi: 10.1023/A:1014945408605.
    !-----------------------------------------------------------------------

    rhoair(:,:) = 0._wp
    DO jb = 1,p_patch%alloc_cell_blocks
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
    DO i = 1, p_ice%kice
      WHERE (p_ice%hi(:,i,:)>0._wp)
        Qatm%SWnet(:,i,:) = ( 1._wp-Qatm%albvisdir(:,i,:) )*fvisdir*p_as%fswr(:,:) +   &
          &                 ( 1._wp-Qatm%albvisdif(:,i,:) )*fvisdif*p_as%fswr(:,:) +   &
          &                 ( 1._wp-Qatm%albnirdir(:,i,:) )*fnirdir*p_as%fswr(:,:) +   &
          &                 ( 1._wp-Qatm%albnirdif(:,i,:) )*fnirdif*p_as%fswr(:,:)
        Tsurf(:,:)    = p_ice%Tsurf(:,i,:)
        fi(:,:)       = 1.0_wp+AAi+p_as%pao(:,:)*(BBi+CCi*Tsurf(:,:) **2)
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
        Qatm%LWnet (:,i,:)  = - fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4 &
           &                  - 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:))
        ! same form as MPIOM:
        !Qatm%LWnet (:,i,:)  = - (fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4 &
        !  &         + 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:)))
        ! bug
        !Qatm%LWnet (:,i,:)  = fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4 &
        !  &     - 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:))
        Qatm%dLWdT (:,i,:)  = -4._wp*zemiss_def*stbo*tafoK(:,:)**3
        Qatm%sens  (:,i,:)  = drags(:,:) * rhoair(:,:)*cpd*p_as%fu10(:,:) * fr_fac &
          &                    * (p_as%tafo(:,:) -Tsurf(:,:))
        Qatm%lat   (:,i,:)  = dragl(:,:) * rhoair(:,:)* alf *p_as%fu10(:,:) * fr_fac &
          &                   * (sphumida(:,:)-sphumidi(:,:))

        Qatm%dsensdT(:,i,:) = 0.95_wp*cpd*rhoair(:,:)*p_as%fu10(:,:)&
          &                  *(dragl0(:,:) - 2.0_wp*dragl(:,:))
        dsphumididesti(:,:) = alpha/(p_as%pao(:,:)-beta*esti(:,:)) &
          &                   * (1.0_wp + beta*esti(:,:)/(p_as%pao(:,:)-beta*esti(:,:)))
        destidT(:,:)        = (bi*ci*di-Tsurf(:,:)*(2.0_wp*ci+Tsurf(:,:)))&
          &                   /(di*(ci+Tsurf(:,:))**2) * esti(:,:)
        dfdT(:,:)               = 2.0_wp*CCi*BBi*Tsurf(:,:)
        Qatm%dlatdT(:,i,:)  = alf*rhoair(:,:)*p_as%fu10(:,:)* &
          &                  ( (sphumida(:,:)-sphumidi(:,:))*dragl1(:,:) &
          &                    - dragl(:,:)*dsphumididesti(:,:)*(fi(:,:)*destidT(:,:) &
          &                    + esti(:,:)*dfdT(:,:)) )
      ENDWHERE
    ENDDO

    !Dirk: why zero ?
    Qatm%rpreci(:,:) = 0.0_wp
    Qatm%rprecw(:,:) = 0.0_wp

    !-----------------------------------------------------------------------
    !  Calculate ice wind stress
    !-----------------------------------------------------------------------

    wspeed(:,:) = SQRT( p_as%u**2 + p_as%v**2 )
    Qatm%stress_x(:,:) = Cd_ia*rhoair(:,:)*wspeed(:,:)*p_as%u(:,:)
    Qatm%stress_y(:,:) = Cd_ia*rhoair(:,:)*wspeed(:,:)*p_as%v(:,:)

    ! Ramp for wind-stress - needed for ice-ocean momentum coupling during spinup
    IF ( PRESENT(datetime) ) THEN
      ramp = MIN(1._wp,(datetime%calday + datetime%caltime &
        - time_config%ini_datetime%calday - time_config%ini_datetime%caltime) / ramp_wind)
      Qatm%stress_x(:,:)  = ramp*Qatm%stress_x(:,:)
      Qatm%stress_y(:,:)  = ramp*Qatm%stress_y(:,:)
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('CalcBulk: stress_x'       ,Qatm%stress_x   ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: stress_y'       ,Qatm%stress_y   ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: tafoK'          ,tafoK           ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: sphumida'       ,sphumida        ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: rhoair'         ,rhoair          ,str_module,idt_src, in_subset=p_patch%cells%owned)
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('CalcBulk: Tsurf ice'      ,Tsurf           ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: Qatm%LWnet ice' ,Qatm%LWnet      ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: Qatm%sens ice'  ,Qatm%sens       ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: Qatm%lat ice'   ,Qatm%lat        ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE calc_bulk_flux_ice

  !-------------------------------------------------------------------------
  !
  !> Forcing_from_bulk equals sbr "Budget_omip" in MPIOM.
  !! Sets the atmospheric fluxes for the update of
  !! temperature of open water for OMIP forcing.
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2012-08). Originally code written by
  !! Dirk Notz, following MPIOM. Code transfered to ICON.
  !
  SUBROUTINE calc_bulk_flux_oce(p_patch, p_as, p_os, Qatm, datetime)
    TYPE(t_patch),            INTENT(IN), TARGET    :: p_patch
    TYPE(t_atmos_for_ocean),  INTENT(IN)    :: p_as
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os
    TYPE(t_atmos_fluxes),     INTENT(INOUT) :: Qatm

    TYPE(t_datetime), OPTIONAL, INTENT(IN)   :: datetime


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
    REAL(wp) :: ramp

    TYPE(t_subset_range), POINTER :: all_cells

    !CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_bulk:calc_bulk_flux_oce
    !-------------------------------------------------------------------------
    !CALL message(TRIM(routine), 'start' )

    Tsurf(:,:)  = p_os%p_prog(nold(1))%tracer(:,1,:,1)  ! set surface temp = mixed layer temp
    tafoK(:,:)  = p_as%tafo(:,:)  + tmelt               ! Change units of tafo  to Kelvin
    ftdewC(:,:) = p_as%ftdew(:,:) - tmelt                    ! Change units of ftdew to C

    ! subset range pointer
    all_cells => p_patch%cells%all



    !-----------------------------------------------------------------------
    ! Compute water vapor pressure and specific humididty in 2m height (esta)
    ! and at water surface (estw) according to "Buck Research Manual (1996)
    ! (see manuals for instruments at http://www.buck-research.com/);
    ! updated from Buck, A. L., New equations for computing vapor pressure and
    ! enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981"
    !-----------------------------------------------------------------------

    aw=611.21_wp; bw=18.729_wp; cw=257.87_wp; dw=227.3_wp
    AAw=7.2e-4_wp; BBw=3.20e-6_wp; CCw=5.9e-10_wp
    alpha=0.62197_wp; beta=0.37803_wp

    fa(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*(BBw+CCw*ftdewC(:,:)**2)
    esta(:,:) = fa(:,:) * aw*EXP((bw-ftdewC(:,:)/dw)*ftdewC(:,:)/(ftdewC(:,:)+cw))
    fw(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*(BBw+CCw*Tsurf(:,:) **2)
    estw(:,:) = fw(:,:) *aw*EXP((bw-Tsurf(:,:) /dw)*Tsurf(:,:) /(Tsurf(:,:) +cw))
    ! For a given surface salinity we should multiply estw with  1 - 0.000537*S

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
    Qatm%LWin(:,:) = 0._wp
    Qatm%LWoutw(:,:) = 0._wp

    ! #eoo# 2012-12-14: another bugfix
    ! #slo# #hha# 2012-12-13: bugfix, corrected form
    Qatm%LWnetw(:,:) = - fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4  &
      &                - 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:))
    ! same form as MPIOM:
    !Qatm%LWnetw(:,:) = - (fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4  &
    !  &         + 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:)))
    ! bug
    !Qatm%LWnetw(:,:) = fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4  &
    !  &         - 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:))

    ! Fractions of SWin in each band (from cice)
    fvisdir=0.28_wp; fvisdif=0.24_wp; fnirdir=0.31_wp; fnirdif=0.17_wp
    Qatm%SWnetw(:,:) = ( 1._wp-Qatm%albvisdirw(:,:) )*fvisdir*p_as%fswr(:,:) +   &
      &                ( 1._wp-Qatm%albvisdifw(:,:) )*fvisdif*p_as%fswr(:,:) +   &
      &                ( 1._wp-Qatm%albnirdirw(:,:) )*fnirdir*p_as%fswr(:,:) +   &
      &                ( 1._wp-Qatm%albnirdifw(:,:) )*fnirdif*p_as%fswr(:,:)

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
    Qatm%sensw(:,:) = drags(:,:)*rhoair(:,:)*cpd*p_as%fu10(:,:) * fr_fac &
      &               * (p_as%tafo(:,:) -Tsurf(:,:))
    Qatm%latw(:,:)  = dragl(:,:)*rhoair(:,:)*alv*p_as%fu10(:,:) * fr_fac &
      &               * (sphumida(:,:)-sphumidw(:,:))

    !-----------------------------------------------------------------------
    !  Calculate oceanic wind stress according to:
    !   Gill (Atmosphere-Ocean Dynamics, 1982, Academic Press) (see also Smith, 1980, J. Phys
    !   Oceanogr., 10, 709-726)
    !-----------------------------------------------------------------------

    wspeed(:,:) = SQRT( p_as%u**2 + p_as%v**2 )
    C_ao(:,:)   = MIN( 2._wp, MAX(1.1_wp, 0.61_wp+0.063_wp*wspeed ) )*1e-3_wp
    Qatm%stress_xw(:,:) = C_ao(:,:)*rhoair*wspeed(:,:)*p_as%u(:,:)
    Qatm%stress_yw(:,:) = C_ao(:,:)*rhoair*wspeed(:,:)*p_as%v(:,:)

    ! Ramp for wind-stress - needed for ice-ocean momentum coupling during spinup
    IF ( PRESENT(datetime) ) THEN
      ramp = MIN(1._wp,(datetime%calday + datetime%caltime &
        - time_config%ini_datetime%calday - time_config%ini_datetime%caltime) / ramp_wind)
      Qatm%stress_xw(:,:) = ramp*Qatm%stress_xw(:,:)
      Qatm%stress_yw(:,:) = ramp*Qatm%stress_yw(:,:)
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('CalcBulk: stress_xw'       ,Qatm%stress_xw  ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: stress_yw'       ,Qatm%stress_yw  ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: tafoK'           ,tafoK           ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: sphumida'        ,sphumida        ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: rhoair'          ,rhoair          ,str_module,idt_src, in_subset=p_patch%cells%owned)
    idt_src=3  ! output print level (1-5, fix)                                                                    
    CALL dbg_print('CalcBulk: Tsurf ocean'     ,Tsurf           ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: Tsurf ocean (nnew)', &
      &     p_os%p_prog(nnew(1))%tracer(:,1,:,1), str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: Qatm%LWnetw'     ,Qatm%LWnetw     ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: Qatm%sensw'      ,Qatm%sensw      ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('CalcBulk: Qatm%latw'       ,Qatm%latw       ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE calc_bulk_flux_oce

  SUBROUTINE update_ice_statistic(p_acc, p_ice, subset)
    TYPE(t_sea_ice_acc),  INTENT(INOUT) :: p_acc
    TYPE(t_sea_ice),      INTENT(IN)    :: p_ice
    TYPE(t_subset_range), INTENT(IN)    :: subset

    CALL add_fields(p_acc%hi  , p_ice%hi  , subset , p_ice%kice,force_level=.TRUE.)
    CALL add_fields(p_acc%hs  , p_ice%hs  , subset , p_ice%kice,force_level=.TRUE.)
    CALL add_fields(p_acc%conc, p_ice%conc, subset , p_ice%kice,force_level=.TRUE.)
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
END MODULE mo_sea_ice
