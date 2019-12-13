!>
!! Contains the implementation of the semi-implicit Adams-Bashforth timestepping
!! for the ICON ocean model using the z* vertical co-ordinate.
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
MODULE mo_ocean_ab_timestepping_zstar
  !-------------------------------------------------------------------------
 
  USE mo_kind,                      ONLY: wp
  USE mo_parallel_config,           ONLY: nproma
  USE mo_sync,                      ONLY: sync_e, sync_c, sync_c1, sync_patch_array, &
    & sync_patch_array_mult, global_sum_array

  USE mo_impl_constants,            ONLY: sea_boundary, max_char_length, min_dolic
  USE mo_dbg_nml,                   ONLY: idbg_mxmn
  USE mo_ocean_nml, ONLY: n_zlev, solver_tolerance,&
    & l_with_vert_tracer_advection, OceanReferenceDensity, OceanReferenceDensity_inv, &
    & para_surfRelax_Temp, para_surfRelax_Salt, l_relaxsal_ice, i_sea_ice, &
    & type_surfRelax_Temp, type_surfRelax_Salt, Analytical_Forcing, &
    & OMIP_FluxFromFile, lhamocc, lfb_bgc_oce, lswr_jerlov, &
    & limit_elevation, &
    & ab_const, ab_beta, ab_gam, iswm_oce, iforc_oce, &
    & no_tracer, l_rigid_lid, l_edge_based,               &
    & use_absolute_solver_tolerance, solver_max_restart_iterations, &
    & solver_max_iter_per_restart, dhdtw_abort, select_transfer, &
    & select_solver, select_gmres, select_gmres_r, select_mres, &
    & select_gmres_mp_r, select_cg, select_cgj, select_bcgs, &
    & select_legacy_gmres, use_continuity_correction, select_cg_mp, &
    & solver_max_iter_per_restart_sp, solver_tolerance_sp, No_Forcing, &
    & MASS_MATRIX_INVERSION_TYPE,            &
    & MASS_MATRIX_INVERSION_ADVECTION, solver_tolerance_comp, &
    & MASS_MATRIX_INVERSION_ALLTERMS, &
    & PPscheme_type, PPscheme_ICON_Edge_vnPredict_type, &
    & solver_FirstGuess, MassMatrix_solver_tolerance,     &
    & createSolverMatrix, l_solver_compare, solver_comp_nsteps
  USE mo_run_config,                ONLY: dtime, debug_check_level, nsteps, output_mode
  USE mo_timer, ONLY: timer_start, timer_stop, timers_level, timer_extra1, &
    & timer_extra2, timer_extra3, timer_extra4, timer_ab_expl, timer_ab_rhs4sfc, timer_total

  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_physical_constants,        ONLY: grav, clw, rho_ref, Tf
  USE mo_ocean_initialization,      ONLY: is_initial_timestep
  USE mo_ocean_types, ONLY: t_hydro_ocean_state
  USE mo_ocean_time_events,   ONLY: ocean_time_nextStep, isCheckpoint, isEndOfThisRun, newNullDatetime
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,            ONLY: t_external_data
  USE mo_exception,                 ONLY: message, finish, warning, message_text
  USE mo_util_dbg_prnt,             ONLY: dbg_print, debug_print_MaxMinMean
  USE mo_ocean_boundcond,           ONLY: VelocityBottomBoundaryCondition_onBlock, top_bound_cond_horz_veloc
  USE mo_ocean_thermodyn,           ONLY: calculate_density, calc_internal_press_grad_zstar
  USE mo_ocean_physics_types,       ONLY: t_ho_params
  USE mo_ocean_pp_scheme,           ONLY: ICON_PP_Edge_vnPredict_scheme
  USE mo_ocean_surface_types,       ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_scalar_product, ONLY: map_edges2edges_viacell_3d_const_z, map_edges2edges_viacell_2D_per_level, &
    & calc_scalar_product_veloc_3d, map_edges2edges_sc_zstar, map_edges2edges_3d_zstar
  USE mo_ocean_math_operators, ONLY: div_oce_3D_onTriangles_onBlock, smooth_onCells, div_oce_3d, &
    & grad_fd_norm_oce_2d_onblock, grad_fd_norm_oce_2d_3d, div_oce_3D_general_onBlock, &
    & div_oce_3D_onTriangles_onBlock, update_height_depdendent_variables
  USE mo_ocean_velocity_advection, ONLY: veloc_adv_horz_mimetic, veloc_adv_vert_mimetic
  USE mo_ocean_velocity_diffusion, ONLY: velocity_diffusion, velocity_diffusion_vertical_implicit_onBlock
  USE mo_ocean_types, ONLY: t_operator_coeff, t_solverCoeff_singlePrecision
  USE mo_grid_subset, ONLY: t_subset_range, get_index_range
  USE mo_grid_config, ONLY: n_dom
  USE mo_mpi, ONLY: work_mpi_barrier, my_process_is_stdio
  USE mo_statistics, ONLY: global_minmaxmean, print_value_location
  USE mo_ocean_solve, ONLY: t_ocean_solve, ocean_solve_ptr
  USE mo_ocean_solve_lhs_type, ONLY: t_lhs_agen, lhs_agen_ptr
  USE mo_ocean_solve_transfer, ONLY: t_transfer, transfer_ptr
  USE mo_ocean_solve_trivial_transfer, ONLY: t_trivial_transfer, trivial_transfer_ptr
  USE mo_ocean_solve_subset_transfer, ONLY: t_subset_transfer, subset_transfer_ptr
  USE mo_ocean_solve_aux, ONLY: t_destructible, t_ocean_solve_parm, solve_gmres, solve_cg, solve_mres, &
   & ocean_solve_clear, solve_precon_none, solve_precon_jac, solve_bcgs, solve_legacy_gmres, &
   & solve_trans_scatter, solve_trans_compact, solve_cell, solve_edge, solve_invalid
  USE mo_primal_flip_flop_lhs, ONLY: t_primal_flip_flop_lhs, lhs_primal_flip_flop_ptr
  USE mo_surface_height_lhs, ONLY: t_surface_height_lhs, lhs_surface_height_ptr
  USE mo_surface_height_lhs_zstar, ONLY: t_surface_height_lhs_zstar, lhs_surface_height_zstar_ptr
  USE mo_ocean_surface_types,    ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_sea_ice_types,          ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_name_list_output_init,  ONLY: isRegistered
  USE mo_hamocc_types,           ONLY: t_hamocc_state
  USE mtime,                     ONLY: datetime, timedelta, newTimedelta, getNoOfSecondsElapsedInDayDateTime, &
    & datetimeToString, deallocateTimedelta
  USE mo_ocean_tracer_transport_types,  ONLY: t_ocean_transport_state, t_ocean_tracer, t_tracer_collection
  USE mo_math_constants,         ONLY: pi, pi_2, rad2deg, deg2rad, dbl_eps
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff, no_primal_edges
  USE mo_ocean_tracer_transport_vert, ONLY: advect_flux_vertical 
  USE mo_memory_log,             ONLY: memory_log_add 
  USE mo_ocean_surface_refactor, ONLY: update_ocean_surface_refactor, update_atmos_fluxes
  USE mo_ocean_physics,         ONLY: update_ho_params
  USE mo_ocean_ab_timestepping_mimetic,  ONLY: calculate_explicit_term_ab, fill_rhs4surface_eq_ab
  USE mo_ice_interface,          ONLY: ice_fast_interface, ice_slow_interface
  USE mo_hydro_ocean_run, ONLY: update_time_g_n, update_time_indices
  USE mo_ocean_output, ONLY: output_ocean
  USE mo_ocean_ab_timestepping,  ONLY: calc_vert_velocity, calc_normal_velocity_ab
  USE mo_ocean_tracer,           ONLY: advect_ocean_tracers
  USE mo_ocean_diagnostics,      ONLY: diag_heat_salt_tendency
  USE mo_ocean_bulk_forcing,  ONLY: apply_surface_relaxation, update_ocean_surface_stress
  USE mo_swr_absorption,     ONLY: dynamic_swr_absorption
  USE mo_ocean_diagnostics,      ONLY: calc_fast_oce_diagnostics, calc_psi
  USE mo_derived_variable_handling, ONLY: update_statistics, reset_statistics
  USE mo_restart_attributes,     ONLY: t_RestartAttributeList, getAttributesForRestarting
  USE mo_master_config,          ONLY: isRestart
  USE mo_sea_ice_nml,            ONLY: i_ice_dyn
  USE mo_restart,                ONLY: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor
  USE mo_ice_fem_interface,      ONLY: ice_fem_init_vel_restart, ice_fem_update_vel_restart

 
  !-------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_vert_velocity_bottomup_zstar
  PUBLIC :: calc_normal_velocity_ab_zstar
  PUBLIC :: solve_free_surface_eq_zstar
  PUBLIC :: update_zstar_variables

! solver object (free ocean surface)
  CLASS(t_destructible), POINTER :: free_sfc_solver => NULL()
! solver object (free ocean surface)
  CLASS(t_destructible), POINTER :: free_sfc_solver_comp => NULL()
! left-hand-side object (free ocean surface)
  CLASS(t_destructible), POINTER :: free_sfc_solver_lhs => NULL()
! communication infrastructure object (free ocean surface)
  CLASS(t_destructible), POINTER :: free_sfc_solver_trans => NULL()
  CLASS(t_destructible), POINTER :: free_sfc_solver_comp_trans => NULL()
    
  TYPE(t_surface_height_lhs_zstar), POINTER :: lhs_sh => NULL()
!
 
  !-------------------------------------------------------------------------
CONTAINS

 

  !! Update stretching variables based on new surface height
  SUBROUTINE update_zstar_variables( patch_3d, eta_c, stretch_c, stretch_e)
    
    TYPE(t_patch_3d ), POINTER, INTENT(in)          :: patch_3d
    REAL(wp), INTENT(INOUT) :: eta_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(INOUT) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp), INTENT(INOUT) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    
    TYPE(t_patch), POINTER :: patch_2d
    
    REAL(wp) :: H_c  (nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! Column depth at cell 
    
    INTEGER  :: jb, jc, je, bt_lev 
    INTEGER  :: start_index, end_index 

    INTEGER, DIMENSION(:,:,:), POINTER :: idx, blk
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    INTEGER  :: id1, id2, bl1, bl2 
    REAL(wp) :: st1, st2 

    !------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    idx             => patch_3D%p_patch_2D(1)%edges%cell_idx
    blk             => patch_3D%p_patch_2D(1)%edges%cell_blk
    
    cells_in_domain  => patch_2D%cells%in_domain
    edges_in_domain  => patch_2D%edges%in_domain
    !------------------------------------------------------------------
 
!ICON_OMP_MASTER
      CALL sync_patch_array(sync_c, patch_2D, eta_c    )
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER


!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, bt_lev) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        
        bt_lev = patch_3d%p_patch_1d(1)%dolic_c(jc, jb)      
 
        H_c  (jc, jb)      = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, jb)
        if ( patch_3D%lsm_c(jc, 1, jb) <= sea_boundary ) THEN
          stretch_c(jc, jb)  = (H_c(jc, jb) + eta_c(jc, jb))/H_c(jc, jb) 
        else 
          stretch_c(jc, jb)  = 1.0_wp
        ENDIF
          
      END DO
    END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_MASTER
      CALL sync_patch_array(sync_c, patch_2D, stretch_c)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER


!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, je, jk, id1, id2, bl1, bl2, st1, st2) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, start_index, end_index)
        DO je = start_index, end_index
          id1 = idx(je, jb, 1)
          id2 = idx(je, jb, 2)
          bl1 = blk(je, jb, 1)
          bl2 = blk(je, jb, 2)
   
          st1 = stretch_c(id1, bl1) 
          st2 = stretch_c(id2, bl2) 
  
          !! FIXME: There seem to be edge cases where this does not work
          IF(patch_3D%lsm_e(je, 1, jb) <= sea_boundary)THEN
            stretch_e(je, jb) = 0.5_wp*(st1 + st2)
          ELSE
            stretch_e(je, jb) = 1.0_wp
          ENDIF
  
  
        ENDDO
      END DO
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_MASTER
      CALL sync_patch_array(sync_e, patch_2D, stretch_e)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER

  END SUBROUTINE update_zstar_variables
 
   
  !! Init variables related to surface height elliptic solver
  SUBROUTINE init_free_sfc(patch_3d, ocean_state, op_coeffs, solverCoeff_sp, str_e, free_sfc_solver, &
      & free_sfc_solver_comp, free_sfc_solver_lhs, free_sfc_solver_trans, free_sfc_solver_comp_trans&
      &, lhs_sh)
    TYPE(t_patch_3d ),POINTER, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(INOUT) :: ocean_state
    TYPE(t_operator_coeff), INTENT(IN), TARGET :: op_coeffs
    TYPE(t_solverCoeff_singlePrecision), INTENT(in), TARGET :: solverCoeff_sp
    REAL(wp), INTENT(IN), CONTIGUOUS :: str_e(:,:)
    CLASS(t_destructible), POINTER, INTENT(INOUT) :: free_sfc_solver 
    CLASS(t_destructible), POINTER, INTENT(INOUT) :: free_sfc_solver_comp 
    CLASS(t_destructible), POINTER, INTENT(INOUT) :: free_sfc_solver_lhs 
    CLASS(t_destructible), POINTER, INTENT(INOUT) :: free_sfc_solver_trans
    CLASS(t_destructible), POINTER, INTENT(INOUT) :: free_sfc_solver_comp_trans
    TYPE(t_surface_height_lhs_zstar), POINTER , INTENT(INOUT):: lhs_sh
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_trivial_transfer), POINTER :: trans_triv
    TYPE(t_subset_transfer), POINTER :: trans_subs
    TYPE(t_ocean_solve), POINTER :: solve, solve_comp
    CLASS(t_lhs_agen), POINTER :: lhs
    CLASS(t_transfer), POINTER :: trans
    TYPE(t_ocean_solve_parm) :: par, par_sp
    INTEGER :: trans_mode, sol_type
    CHARACTER(len=*), PARAMETER :: method_name='mo_ocean_ab_timestepping_mimetic:init_free_sfc_ab_mimetic'

!
    IF (ASSOCIATED(free_sfc_solver)) RETURN
    patch_2D => patch_3d%p_patch_2d(1)
    NULLIFY(trans_triv)
! allocate lhs object
    ALLOCATE(t_surface_height_lhs_zstar :: free_sfc_solver_lhs)
! init lhs object
    lhs_sh => lhs_surface_height_zstar_ptr(free_sfc_solver_lhs)

    CALL lhs_sh%construct(patch_3d, ocean_state%p_diag%thick_e, &
      & op_coeffs, solverCoeff_sp, str_e)
    lhs => lhs_agen_ptr(free_sfc_solver_lhs)
! allocate and init communication infrastructure object 
    SELECT CASE(select_transfer)
! all ocean workers are involved in solving (input is just copied to internal
! arrays)
    CASE(0)
      ALLOCATE(t_trivial_transfer :: free_sfc_solver_trans)
      trans_triv => trivial_transfer_ptr(free_sfc_solver_trans)
      CALL trans_triv%construct(solve_cell, patch_2D)
! solve only on a subset of workers
    CASE DEFAULT
      trans_mode = MERGE(solve_trans_compact, solve_trans_scatter, select_transfer .GT. 0)
      ALLOCATE(t_subset_transfer :: free_sfc_solver_trans)
      trans_subs => subset_transfer_ptr(free_sfc_solver_trans)
      CALL trans_subs%construct(solve_cell, patch_2D, ABS(select_transfer), trans_mode)
    END SELECT
    trans => transfer_ptr(free_sfc_solver_trans)
!prepare init of solver
    CALL par%init(solve_precon_none, 1, 800, patch_2d%cells%in_domain%end_block, &
      & patch_2D%alloc_cell_blocks, nproma, patch_2d%cells%in_domain%end_index, &
      & solver_tolerance, use_absolute_solver_tolerance)
    par_sp%nidx = -1 ! indicates not using sp-solver
    sol_type = solve_gmres
! decide which solver type to use
    SELECT CASE(select_solver)
    CASE(select_gmres) ! GMRES
    CASE(select_gmres_r) ! GMRES-restart
      par%m = solver_max_iter_per_restart
      par%nr = solver_max_restart_iterations
    CASE(select_gmres_mp_r) ! GMRES-restart(sp+wp)
      par%nr = solver_max_restart_iterations
      par_sp = par
      par_sp%m = solver_max_iter_per_restart_sp
      par_sp%tol = REAL(solver_tolerance_sp, wp)
      par%m = solver_max_iter_per_restart
    CASE(select_cg) ! CG (Fletcher-Reeves)
      sol_type = solve_cg
    CASE(select_cg_mp) ! CG (Fletcher-Reeves, sp+wp)
      sol_type = solve_cg
      par_sp = par
      par_sp%tol = REAL(solver_tolerance_sp, wp)
    CASE(select_cgj) ! CG + Jacobi preconditioner
      sol_type = solve_cg
      par%pt = solve_precon_jac
    CASE(select_bcgs) ! BiCG-Stab
      sol_type = solve_bcgs
    CASE(select_legacy_gmres)
      sol_type = solve_legacy_gmres
      par%m = solver_max_iter_per_restart
      par_sp%nr = solver_max_restart_iterations
    CASE(select_mres)
      sol_type = solve_mres
      par%nr = (par%m + 18) / 19
      par%m = 19
    CASE DEFAULT
      CALL finish(method_name, "Unknown solver")
    END SELECT ! solver
! alloc and init solver object
    ALLOCATE(t_ocean_solve :: free_sfc_solver)
    solve => ocean_solve_ptr(free_sfc_solver)
    CALL solve%construct(sol_type, par, par_sp, lhs, trans)
    IF (l_solver_compare) THEN
      IF (ASSOCIATED(trans_triv)) THEN
        NULLIFY(free_sfc_solver_comp_trans)
      ELSE
        ALLOCATE(t_trivial_transfer :: free_sfc_solver_comp_trans)
        trans_triv => trivial_transfer_ptr(free_sfc_solver_comp_trans)
        CALL trans_triv%construct(solve_cell, patch_2D)
        trans => transfer_ptr(free_sfc_solver_comp_trans)
      END IF
      par%tol = solver_tolerance_comp
      par_sp%nidx = -1
      ALLOCATE(t_ocean_solve :: free_sfc_solver_comp)
      solve_comp => ocean_solve_ptr(free_sfc_solver_comp)
      CALL solve_comp%construct(solve_legacy_gmres, par, par_sp, lhs, trans)
    END IF
  END SUBROUTINE init_free_sfc


  !-------------------------------------------------------------------------
  !>
  !! Computation of new velocity in Adams-Bashforth timestepping.
  !!
  !! Adapted to zstar
  !!
  SUBROUTINE calc_normal_velocity_ab_zstar(patch_3d,ocean_state, op_coeffs, eta)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs
    REAL(wp), INTENT(in)                 :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !  local variables
    INTEGER  :: start_edge_index, end_edge_index, je, jk, blockNo
    REAL(wp) :: gdt_x_ab_beta, one_minus_ab_gam
    REAL(wp) :: z_grad_h(nproma,patch_3d%p_patch_2d(1)%nblks_e), z_grad_h_block(nproma)
    TYPE(t_subset_range), POINTER :: edges_in_domain, owned_edges
    TYPE(t_patch), POINTER :: patch
    !----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    !-----------------------------------------------------------------------
    patch           => patch_3d%p_patch_2d(1)
    edges_in_domain => patch%edges%in_domain
    owned_edges     => patch%edges%owned
    one_minus_ab_gam = 1.0_wp - ab_gam
    gdt_x_ab_beta =  grav * dtime * ab_beta

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,start_edge_index,end_edge_index, je, jk, z_grad_h_block) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      ! Step 1) Compute normal derivative of new surface height
      CALL grad_fd_norm_oce_2d_onBlock(eta, patch, op_coeffs%grad_coeff(:,1, blockNo), &
        & z_grad_h_block(:), start_edge_index, end_edge_index, blockNo)
      ! Step 2) Calculate the new velocity from the predicted one and the new surface height
      DO je = start_edge_index, end_edge_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
           ocean_state%p_prog(nnew(1))%vn(je,jk,blockNo) = (ocean_state%p_diag%vn_pred(je,jk,blockNo) &
             & - gdt_x_ab_beta * z_grad_h_block(je))
        END DO          
      END DO
      DO je = start_edge_index, end_edge_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          ocean_state%p_diag%vn_time_weighted(je,jk,blockNo) = ab_gam * ocean_state%p_prog(nnew(1))%vn(je,jk,blockNo) &
            & + one_minus_ab_gam * ocean_state%p_prog(nold(1))%vn(je,jk,blockNo)
        END DO
      END DO
    END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO

    CALL sync_patch_array_mult(sync_e, patch, 2, ocean_state%p_prog(nnew(1))%vn, ocean_state%p_diag%vn_time_weighted)
  
  END SUBROUTINE calc_normal_velocity_ab_zstar
  !-------------------------------------------------------------------------
 
  
  !-------------------------------------------------------------------------
  !>
  !! Computation of new vertical velocity using continuity equation
  !! Calculate diagnostic vertical velocity from horizontal velocity using the
  !! incommpressibility condition in the continuity equation.
  !! vertical velocity is integrated from bottom to topLevel
  !! vertical velocity is negative for positive divergence
  !! of horizontal velocity
  !!
  !! Adapted to zstar
  !!
  SUBROUTINE calc_vert_velocity_bottomup_zstar( patch_3d, ocean_state, op_coeffs, stretch_c, stretch_e)
    TYPE(t_patch_3d), TARGET :: patch_3d       ! patch on which computation is performed
    TYPE(t_hydro_ocean_state) :: ocean_state
    TYPE(t_operator_coeff), INTENT(in) :: op_coeffs
    REAL(wp), INTENT(IN)               :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp), INTENT(IN)               :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    ! Local variables
    INTEGER :: jc, jk, blockNo, start_index, end_index
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells, cells_owned
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp), POINTER :: vertical_velocity(:,:,:)
    REAL(wp) :: div_m_c(nproma, n_zlev)
    REAL(wp) :: deta_dt, H_c 
    INTEGER  :: bt_lev
 
    !-----------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain  => patch_2D%cells%in_domain
    cells_owned      => patch_2D%cells%owned
    all_cells        => patch_2D%cells%all
    edges_in_domain  => patch_2D%edges%in_domain
    vertical_velocity=> ocean_state%p_diag%w
    ! due to nag -nan compiler-option:
    !------------------------------------------------------------------
    ! Step 1) Calculate divergence of horizontal velocity at all levels
    !------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    CALL map_edges2edges_3d_zstar( patch_3d, ocean_state%p_diag%vn_time_weighted, op_coeffs, &
      & stretch_e, ocean_state%p_diag%mass_flx_e)
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      CALL div_oce_3D_onTriangles_onBlock(ocean_state%p_diag%mass_flx_e, patch_3D, op_coeffs%div_coeff, &
        & ocean_state%p_diag%div_mass_flx_c(:,:,blockNo), blockNo=blockNo, start_index=start_index, &
        & end_index=end_index, start_level=1, end_level=n_zlev)
      div_m_c(:, :) = ocean_state%p_diag%div_mass_flx_c(:,:,blockNo)
      DO jc = start_index, end_index
        !use bottom boundary condition for vertical velocity at bottom of prism
        ! this should be awlays zero
        
        bt_lev = patch_3d%p_patch_1d(1)%dolic_c(jc, blockNo)      
        H_c    = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, blockNo)

        deta_dt = -SUM(div_m_c(jc, 1:patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)))
        vertical_velocity(jc, patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo) + 1, blockNo) = 0.0_wp
        DO jk = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), 1, -1
          vertical_velocity(jc,jk,blockNo) = vertical_velocity(jc,jk+1,blockNo) - &
            & ( ocean_state%p_diag%div_mass_flx_c(jc,jk,blockNo) + &
            & (1.0_wp/H_c)*                       &
            & deta_dt*patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, blockNo) &
            & )/stretch_c(jc, blockNo)
        END DO
      END DO
    END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO
    
    CALL sync_patch_array(sync_c,patch_2D,vertical_velocity)

  END SUBROUTINE calc_vert_velocity_bottomup_zstar
  !-------------------------------------------------------------------------
 
  
  !-------------------------------------------------------------------------
  !>
  !! Computation of velocity predictor in Adams-Bashforth timestepping.
  !!
  SUBROUTINE calculate_explicit_term_zstar( patch_3d, ocean_state, p_phys_param,&
    & is_first_timestep, op_coeffs, p_as, stretch_c, eta_c, stretch_e)
    TYPE(t_patch_3d ), POINTER, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE (t_ho_params)                   :: p_phys_param
    LOGICAL,INTENT(in)                   :: is_first_timestep
    TYPE(t_operator_coeff), INTENT(IN), TARGET :: op_coeffs
    TYPE(t_atmos_for_ocean), INTENT(inout) :: p_as
    REAL(wp), INTENT(IN) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! stretch factor 
    REAL(wp), INTENT(IN) :: eta_c    (nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! stretch factor 
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    TYPE(t_subset_range), POINTER :: owned_edges, owned_cells
    
    owned_edges     => patch_3d%p_patch_2d(n_dom)%edges%owned
    owned_cells     => patch_3d%p_patch_2d(n_dom)%cells%owned

    ! STEP 1: horizontal advection
    CALL veloc_adv_horz_mimetic( patch_3d,         &
      & ocean_state%p_prog(nold(1))%vn,    &
      & ocean_state%p_prog(nold(1))%vn,    &
      & ocean_state%p_diag,                &
      & ocean_state%p_diag%veloc_adv_horz, &
      & op_coeffs)
    
    ! STEP 2: compute 3D contributions: gradient of hydrostatic pressure and vertical velocity advection
      
    ! calculate density from EOS using temperature and salinity at timelevel n
    ! FIXME zstar: Uses depth to calculate pressure, will need fix 
    CALL calculate_density( patch_3d,                         &
     & ocean_state%p_prog(nold(1))%tracer(:,:,:,1:no_tracer),&
     & ocean_state%p_diag%rho(:,:,:) )

    CALL calc_internal_press_grad_zstar( patch_3d,&
       &                          ocean_state%p_diag%rho,&
       &                          ocean_state%p_diag%press_hyd,& 
       &                          ocean_state%p_aux%bc_total_top_potential, &
       &                          op_coeffs%grad_coeff,  &
       &                          stretch_c,             &
       &                          ocean_state%p_diag%press_grad)     
    ! calculate vertical velocity advection
    !! FIXME: All derivatives are calculated from level = 2
    !! Level=1 derivatives seem to be assumed to be 0
    !! With this assumption no changes would be required
    CALL veloc_adv_vert_mimetic(          &
      & patch_3d,                         &
      & ocean_state%p_diag,op_coeffs,     &
      & ocean_state%p_diag%veloc_adv_vert )

    ! STEP 3: compute harmonic or biharmoic laplacian diffusion of velocity.
    !         This term is discretized explicitly. Order and form of the laplacian
    !         are determined in mo_oce_diffusion according to namelist settings
    CALL velocity_diffusion(patch_3d,              &
      & ocean_state%p_prog(nold(1))%vn, &
      & p_phys_param,            &
      & ocean_state%p_diag,op_coeffs,  &
      & ocean_state%p_diag%laplacian_horz)
    
    CALL explicit_vn_pred_zstar( patch_3d, ocean_state, op_coeffs, p_phys_param, is_first_timestep, &
      & eta_c, stretch_e)
    
  END SUBROUTINE calculate_explicit_term_zstar

  !-------------------------------------------------------------------------
  SUBROUTINE calculate_explicit_term_g_n_onBlock_zstar( patch_3d, ocean_state, is_first_timestep, &
    & start_edge_index, end_edge_index, blockNo)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    LOGICAL,INTENT(in)                   :: is_first_timestep
    INTEGER,INTENT(in)                   :: start_edge_index, end_edge_index, blockNo
    INTEGER :: je, jk
    !-----------------------------------------------------------------------
    
    DO je = start_edge_index, end_edge_index
      DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
        ocean_state%p_aux%g_n(je, jk, blockNo) = &
          & - ocean_state%p_diag%press_grad    (je, jk, blockNo)  &
          & - ocean_state%p_diag%grad          (je, jk, blockNo)  &            
          & - ocean_state%p_diag%veloc_adv_horz(je, jk, blockNo)  &
          & - ocean_state%p_diag%veloc_adv_vert(je, jk, blockNo)  &
          & + ocean_state%p_diag%laplacian_horz(je, jk, blockNo)  
      END DO
    END DO
    
    IF(is_first_timestep)THEN
      ocean_state%p_aux%g_nimd(1:nproma,1:n_zlev, blockNo) = &
        & ocean_state%p_aux%g_n(1:nproma,1:n_zlev,blockNo)
    ELSE
      DO je = start_edge_index, end_edge_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          ocean_state%p_aux%g_nimd(je, jk,blockNo)                          &
            & = (1.5_wp+ab_const) * ocean_state%p_aux%g_n(je, jk,blockNo)   &
            & - (0.5_wp+ab_const) * ocean_state%p_aux%g_nm1(je, jk,blockNo)
        END DO
      END DO
    ENDIF

  END SUBROUTINE calculate_explicit_term_g_n_onBlock_zstar
 

  !-------------------------------------------------------------------------
  SUBROUTINE calculate_explicit_vn_pred_3D_onBlock_zstar( patch_3d, ocean_state, z_gradh_e, &
    & start_edge_index, end_edge_index, blockNo)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    REAL(wp) :: z_gradh_e(nproma)
    INTEGER, INTENT(in) :: start_edge_index, end_edge_index, blockNo
    INTEGER :: je, jk, bottom_level
    !-----------------------------------------------------------------------
    DO je = start_edge_index, end_edge_index
      DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
        ocean_state%p_diag%vn_pred(je,jk,blockNo) = ocean_state%p_prog(nold(1))%vn(je,jk,blockNo)  &
          & + dtime*(ocean_state%p_aux%g_nimd(je,jk,blockNo) &
          & - z_gradh_e(je))
      END DO
    END DO
    
    CALL VelocityBottomBoundaryCondition_onBlock(patch_3d, &
      & blockNo,start_edge_index, end_edge_index, &
      & ocean_state%p_prog(nold(1))%vn(:,:,blockNo), &
      & ocean_state%p_diag%vn_pred(:,:,blockNo),     &
      & ocean_state%p_aux%bc_bot_vn(:,blockNo))

    !IF surface forcing applied as topLevel boundary condition to vertical diffusion
    !The surface forcing is applied as volume forcing at rhs,
    !i.e. if it part of explicit term in momentum and tracer eqs.
    !in this case, topLevel boundary ondition of vertical Laplacians are homogeneous.
    !Below is the code that adds surface forcing to explicit term of momentum eq.
    !FIXME zstar: Do these boundary conditions make sense 
    DO je = start_edge_index, end_edge_index
      IF(patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)>=min_dolic) THEN
        ocean_state%p_diag%vn_pred(je,1,blockNo) =  ocean_state%p_diag%vn_pred(je,1,blockNo)      &
          & + dtime*ocean_state%p_aux%bc_top_vn(je,blockNo)                                  &
          & /patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(je,1,blockNo) ! Change to prism_thick_e ?
        bottom_level = patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
        ocean_state%p_diag%vn_pred(je,bottom_level,blockNo)                                  &
          & = ocean_state%p_diag%vn_pred(je,bottom_level,blockNo)                            &
          & - dtime*ocean_state%p_aux%bc_bot_vn(je,blockNo)                                  &
          & /patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(je,bottom_level,blockNo)
      ENDIF
    END DO
  END SUBROUTINE calculate_explicit_vn_pred_3D_onBlock_zstar

  
  !-------------------------------------------------------------------------
  !!Subroutine implements implicit vertical diffusion for horizontal velocity fields
  !!by inverting a scalar field..
  !------------------------------------------------------------------------
  SUBROUTINE velocity_diffusion_vertical_implicit_zstar( &
    & patch_3d,                            &
    & velocity,                            &
    & a_v, stretch_e,                      &
    & operators_coefficients,              &
    & start_index, end_index, edge_block)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)              :: velocity(:,:)   ! on edges, (nproma, levels)
    REAL(wp), INTENT(inout)              :: a_v(:,:)      ! on edges, (nproma, levels)
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    TYPE(t_operator_coeff),INTENT(IN) ,TARGET :: operators_coefficients
    INTEGER , INTENT(in):: start_index, end_index, edge_block
    !
    REAL(wp) :: dt_inv
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev), diagonal_product
    REAL(wp) :: column_velocity(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    REAL(wp) :: inv_str_e 

    INTEGER :: bottom_level
    INTEGER :: edge_index, level

    !-----------------------------------------------------------------------
    dt_inv = 1.0_wp/dtime

    DO edge_index = start_index, end_index
      bottom_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)

      inv_str_e    = 1._wp/stretch_e(edge_index, edge_block)

      IF (bottom_level < 2 ) CYCLE ! nothing to diffuse

      ! Note : the inv_prism_thick_e, inv_prism_center_dist_e should be updated in calculate_thickness
      DO level=1, bottom_level
        inv_prism_thickness(level)        = inv_str_e*patch_3d%p_patch_1d(1)%inv_prism_thick_e(edge_index,level,edge_block)
        inv_prisms_center_distance(level) = inv_str_e*patch_3d%p_patch_1d(1)%inv_prism_center_dist_e(edge_index,level,edge_block)
      ENDDO


      !------------------------------------
      ! Fill triangular matrix
      ! b is diagonal, a is the upper diagonal, c is the lower
      !   top level
      a(1) = 0.0_wp
      c(1) = -a_v(edge_index,2) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
      b(1) = dt_inv - c(1)
      DO level = 2, bottom_level-1
        a(level) = - a_v(edge_index,level)   * inv_prism_thickness(level) * inv_prisms_center_distance(level)
        c(level) = - a_v(edge_index,level+1) * inv_prism_thickness(level) * inv_prisms_center_distance(level+1)
        b(level) = dt_inv - a(level) - c(level)
      END DO
      ! bottom
      a(bottom_level) = -a_v(edge_index,bottom_level) * inv_prism_thickness(bottom_level) * &
        & inv_prisms_center_distance(bottom_level)
      b(bottom_level) = dt_inv - a(bottom_level)

      ! precondition: set diagonal equal to diagonal_product
      diagonal_product = PRODUCT(b(1:bottom_level))

      DO level = 1, bottom_level
        fact(level) = diagonal_product / b(level)
        a(level)  = a(level)  * fact(level)
        b(level)  = diagonal_product
        c(level)  = dt_inv * fact(level) - a(level) - b(level)
        column_velocity(level) = velocity(edge_index,level) * dt_inv * fact(level)
      ENDDO
      c(bottom_level) = 0.0_wp

      !------------------------------------
      ! solver from lapack
      !
      ! eliminate lower diagonal
      DO level = bottom_level-1, 1, -1
        fact(level+1)  = c( level ) / b( level+1 )
        b( level ) = b( level ) - fact(level+1) * a( level +1 )
        column_velocity( level ) = column_velocity( level ) - fact(level+1) * column_velocity( level+1 )
      ENDDO

      !     Back solve with the matrix U from the factorization.
      column_velocity( 1 ) = column_velocity( 1 ) / b( 1 )
      DO level = 2, bottom_level
        column_velocity( level ) = ( column_velocity( level ) - a( level ) * column_velocity( level-1 ) ) / b( level )
      ENDDO

      DO level = 1, bottom_level
        velocity(edge_index,level) = column_velocity(level)
      ENDDO

    END DO ! edge_index = start_index, end_index

  END SUBROUTINE velocity_diffusion_vertical_implicit_zstar
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE explicit_vn_pred_zstar( patch_3d, ocean_state, op_coeffs, p_phys_param, &
      & is_first_timestep, eta, stretch_e)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_operator_coeff), INTENT(IN) :: op_coeffs
    TYPE (t_ho_params) :: p_phys_param
    LOGICAL, INTENT(in)  :: is_first_timestep
    REAL(wp), INTENT(IN) :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    REAL(wp) :: z_gradh_e(nproma)
    TYPE(t_subset_range), POINTER :: edges_in_domain
    INTEGER :: start_edge_index, end_edge_index, blockNo
    TYPE(t_patch), POINTER :: patch_2D

    patch_2D        => patch_3d%p_patch_2d(n_dom)
    edges_in_domain => patch_3d%p_patch_2d(n_dom)%edges%in_domain
    ! STEP 4: calculate weighted gradient of surface height at previous timestep
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, z_gradh_e) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      z_gradh_e(:)  = 0.0_wp
      CALL grad_fd_norm_oce_2d_onBlock( eta, patch_2D, &
        & op_coeffs%grad_coeff(:,1,blockNo), z_gradh_e(:), start_edge_index, end_edge_index, blockNo)
      z_gradh_e(start_edge_index:end_edge_index) = &
         & (1.0_wp-ab_beta) * grav * z_gradh_e(start_edge_index:end_edge_index)
      ! STEP 5:
      CALL calculate_explicit_term_g_n_onBlock_zstar( patch_3d, ocean_state, is_first_timestep, &
        & start_edge_index, end_edge_index, blockNo)
        
      CALL calculate_explicit_vn_pred_3D_onBlock_zstar( patch_3d, ocean_state, z_gradh_e(:),    &
      & start_edge_index, end_edge_index, blockNo)
      ! calculate vertical friction, ie p_phys_param%a_veloc_v
      ! FIXME zstar: This is not modified for zstar
      ! requires density gradient, where is that calculated 
      IF (PPscheme_type == PPscheme_ICON_Edge_vnPredict_type) &
        CALL ICON_PP_Edge_vnPredict_scheme(patch_3d, blockNo, start_edge_index, end_edge_index, &
          & ocean_state, ocean_state%p_diag%vn_pred(:,:,blockNo))
      !In 3D case implicit vertical velocity diffusion is chosen
      CALL velocity_diffusion_vertical_implicit_zstar(patch_3d, &
        & ocean_state%p_diag%vn_pred(:,:,blockNo), p_phys_param%a_veloc_v(:,:,blockNo), stretch_e, &
        & op_coeffs, start_edge_index, end_edge_index, blockNo)

    END DO
!ICON_OMP_END_PARALLEL_DO
  END SUBROUTINE explicit_vn_pred_zstar
  
  
  !-------------------------------------------------------------------------
  !>
  !!  Calculation of right-hand side of elliptic surface equation.
  !!  This is used in semi implicit timelevel stepping.
  !!
  SUBROUTINE fill_rhs4surface_eq_zstar( patch_3d, ocean_state, op_coeffs, stretch_e, eta)
    ! Patch on which computation is performed
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_operator_coeff), INTENT(IN) :: op_coeffs
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    REAL(wp), INTENT(IN) :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: jc, blockNo, jk, je
    REAL(wp) ::inv_gdt2    !, delta_z
    REAL(wp) :: z_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: div_z_depth_int_c(nproma)
    REAL(wp) :: div_z_c(nproma,n_zlev)
    REAL(wp) :: z_vn_ab(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, owned_edges
    TYPE(t_patch), POINTER :: patch_2D

    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_3d%p_patch_2d(1)%cells%in_domain
    edges_in_domain => patch_3d%p_patch_2d(1)%edges%in_domain
    owned_edges     => patch_3d%p_patch_2d(1)%edges%owned
    inv_gdt2 = 1.0_wp / (grav*dtime*dtime)
    z_vn_ab(:,:,:edges_in_domain%start_block-1) = 0._wp

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, je, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      z_vn_ab(:,:,blockNo)  = 0.0_wp
      DO je = start_edge_index, end_edge_index
!DIR$ SIMD
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          z_vn_ab(je,jk,blockNo)=ab_gam * ocean_state%p_diag%vn_pred(je,jk,blockNo) &
            & + (1.0_wp -ab_gam) * ocean_state%p_prog(nold(1))%vn(je,jk,blockNo)
        END DO
      ENDDO
    END DO
!ICON_OMP_END_PARALLEL_DO

    z_vn_ab(:,:,edges_in_domain%end_block+1:) = 0._wp
    !
    ! calculate depth-integrated velocity z_e
    !  - edge-based and cell-based
    !  - 3d and 2d (surface)
    !-------------------------------------------------------------------------------
      
    CALL sync_patch_array(sync_e, patch_2D, z_vn_ab)
    CALL map_edges2edges_3d_zstar( patch_3d, z_vn_ab, op_coeffs, stretch_e, z_e )
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, jc, jk, div_z_depth_int_c, div_z_c) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      CALL div_oce_3D_onTriangles_onBlock( z_e, patch_3D, op_coeffs%div_coeff, div_z_c, &
        & blockNo=blockNo, start_index=start_cell_index, end_index=end_cell_index,      &
        & start_level=1, end_level=n_zlev)
      ! integrate div on columns
      div_z_depth_int_c(:) = 0.0_wp
      DO jc = start_cell_index, end_cell_index
        div_z_depth_int_c(jc) = SUM(div_z_c(jc, 1:patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)))
      END DO
      ocean_state%p_aux%p_rhs_sfc_eq(:,blockNo) = 0.0_wp
      DO jc = start_cell_index, end_cell_index
        IF (patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo) > 0) THEN
          ocean_state%p_aux%p_rhs_sfc_eq(jc,blockNo) = ( ( eta(jc,blockNo) &
            & - dtime * div_z_depth_int_c(jc)) * inv_gdt2)

        ENDIF
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO
  
  END SUBROUTINE fill_rhs4surface_eq_zstar

 
   

  SUBROUTINE solve_free_surface_eq_zstar( patch_3d, ocean_state, p_ext_data,  &
    & p_oce_sfc , p_as, p_phys_param, operators_coefficients, solvercoeff_sp, &
    & timestep, eta_c, stretch_c, stretch_e, eta_c_new, stretch_c_new)

      TYPE(t_patch_3d ), POINTER, INTENT(in)             :: patch_3d
      TYPE(t_hydro_ocean_state), TARGET, INTENT(inout)   :: ocean_state(n_dom)
      TYPE(t_external_data), TARGET, INTENT(in)          :: p_ext_data(n_dom)
      TYPE(t_ocean_surface)                              :: p_oce_sfc
      TYPE(t_atmos_for_ocean),  INTENT(inout)            :: p_as
      TYPE(t_ho_params)                                  :: p_phys_param 
      TYPE(t_operator_coeff),   INTENT(inout)            :: operators_coefficients
      TYPE(t_solvercoeff_singleprecision), INTENT(inout) :: solvercoeff_sp
   
      INTEGER , INTENT(IN   ) :: timestep 
      REAL(wp), INTENT(IN   ) :: eta_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
      REAL(wp), INTENT(IN   ) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
      REAL(wp), INTENT(INOUT) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
      REAL(wp), INTENT(INOUT) :: eta_c_new(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
      REAL(wp), INTENT(INOUT) :: stretch_c_new(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    
      TYPE(t_patch), POINTER :: patch_2d
    
      TYPE(t_subset_range), POINTER :: owned_cells, owned_edges
    
      TYPE(t_surface_height_lhs_zstar), POINTER :: lhs_sh => NULL()
      TYPE(t_ocean_solve), POINTER :: solve, solve_comp
    
      INTEGER :: n_it, n_it_sp, ret_status
    
      CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_timestepping_zstar:solve_free_sfc'
    
      CHARACTER(LEN=12)  :: str_module = 'zstar_surf'  ! Output of module for 1 line debug)
    
      REAL(wp) :: rn, minmaxmean(3)

      !------------------------------------------------------------------------
      patch_2d        => patch_3d%p_patch_2d(1)
      owned_cells => patch_2D%cells%owned
      owned_edges => patch_2D%edges%owned

      !------------------------------------------------------------------------
      ! solve for new free surface
      !------------------------------------------------------------------------
      
      !! FIXME: The RHS can be filled explicitly here, however, the LHS requires
      !! multiplying the Beta*Gamma term with (H+eta)/H
      !! The height goes in using map_edges2edges_viacell_2D 
      !! init sfc solver related objects, if necessary
      IF (.NOT.ASSOCIATED(free_sfc_solver)) &
        CALL init_free_sfc(patch_3d, ocean_state(1), operators_coefficients, solverCoeff_sp, stretch_e, &
        & free_sfc_solver, free_sfc_solver_comp, free_sfc_solver_lhs, free_sfc_solver_trans, &
        & free_sfc_solver_comp_trans, lhs_sh)
      solve => ocean_solve_ptr(free_sfc_solver)
      IF (l_solver_compare) solve_comp => ocean_solve_ptr(free_sfc_solver_comp)
      
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('test    : h var'   ,ocean_state(1)%p_prog(nold(1))%h ,str_module, 2, in_subset=owned_cells)
      CALL dbg_print('on entry: h-new'   ,eta_c_new                         ,str_module, 2, in_subset=owned_cells)
      CALL dbg_print('on entry: vn-new'  ,ocean_state(1)%p_prog(nnew(1))%vn,str_module, 2, in_subset=owned_edges)

      ! Apply windstress
      CALL top_bound_cond_horz_veloc(patch_3d, ocean_state(1), operators_coefficients, p_oce_sfc)
      
      CALL calculate_explicit_term_zstar(patch_3d, ocean_state(1), p_phys_param, &
        & is_initial_timestep(timestep), operators_coefficients, p_as, stretch_c, eta_c, stretch_e)
      
      ! Calculate RHS of surface equation
      CALL fill_rhs4surface_eq_zstar(patch_3d, ocean_state(1), operators_coefficients, stretch_e, eta_c)

      !! Update stretching co-efficient for LHS
      CALL lhs_sh%update(stretch_e)

      ! Solve surface equation with solver

      !!ICON_OMP PARALLEL WORKSHARE
      solve%x_loc_wp(:,:) = eta_c(:,:)
      !!ICON_OMP END PARALLEL WORKSHARE

      solve%b_loc_wp => ocean_state(1)%p_aux%p_rhs_sfc_eq
  
      ! call solver
      CALL solve%solve(n_it, n_it_sp)
      rn = MERGE(solve%res_loc_wp(1), 0._wp, n_it .NE. 0)
      
      IF (rn > solver_tolerance) THEN
        ret_status = 2
        CALL warning(routine, "NOT YET CONVERGED !!")
        RETURN
      ENDIF
      !!ICON_OMP PARALLEL WORKSHARE
      eta_c_new(:,:) = solve%x_loc_wp(:,:)
      !!ICON_OMP END PARALLEL WORKSHARE
 
      IF (createSolverMatrix) &
        CALL solve%dump_matrix(timestep)
  
      !-------- end of solver ---------------
      !---------------------------------------------------------------------
 
      CALL update_zstar_variables( patch_3d, eta_c_new, stretch_c_new, stretch_e)

      !------------------------------------------------------------------------
      ! end solve for free surface and update stretching
      !------------------------------------------------------------------------
     
  END SUBROUTINE solve_free_surface_eq_zstar


  
END MODULE mo_ocean_ab_timestepping_zstar
