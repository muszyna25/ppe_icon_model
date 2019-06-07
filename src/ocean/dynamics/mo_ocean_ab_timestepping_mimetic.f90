!>
!! Contains the implementation of the semi-implicit Adams-Bashforth timestepping
!! for the ICON ocean model based on the mimetic spatiaöl discretization approach.
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2010/04)
!!  Modified by Stephan Lorenz,     MPI-M (2010-06)
!!   - renaming and adjustment to ocean domain and patch
!!   - implementation of continuity equation for vertical velocities
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
#include "iconfor_dsl_definitions.inc"
#include "omp_definitions.inc"
#include "icon_definitions.inc"
#include "crayftn_ptr_fail.inc"
!----------------------------
MODULE mo_ocean_ab_timestepping_mimetic

  USE mo_kind,                      ONLY: wp
  USE mo_parallel_config,           ONLY: nproma
  USE mo_sync,                      ONLY: sync_e, sync_c, sync_patch_array, sync_patch_array_mult
  USE mo_impl_constants,            ONLY: sea_boundary, max_char_length, min_dolic
  USE mo_dbg_nml,                   ONLY: idbg_mxmn
  USE mo_ocean_nml, ONLY: n_zlev, solver_tolerance,&
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
  USE mo_run_config,                ONLY: dtime, debug_check_level
  USE mo_timer, ONLY: timer_start, timer_stop, timers_level, timer_extra1, &
    & timer_extra2, timer_extra3, timer_extra4, timer_ab_expl, timer_ab_rhs4sfc
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_physical_constants,        ONLY: grav
  USE mo_ocean_initialization,      ONLY: is_initial_timestep
  USE mo_ocean_types, ONLY: t_hydro_ocean_state
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,            ONLY: t_external_data
  USE mo_exception,                 ONLY: message, finish, warning, message_text
  USE mo_util_dbg_prnt,             ONLY: dbg_print, debug_print_MaxMinMean
  USE mo_ocean_boundcond,           ONLY: VelocityBottomBoundaryCondition_onBlock, top_bound_cond_horz_veloc
  USE mo_ocean_thermodyn,           ONLY: calculate_density, calc_internal_press_grad
  USE mo_ocean_physics_types,       ONLY: t_ho_params
  USE mo_ocean_pp_scheme,           ONLY: ICON_PP_Edge_vnPredict_scheme
  USE mo_ocean_surface_types,       ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_scalar_product, ONLY: map_edges2edges_viacell_3d_const_z, map_edges2edges_viacell_2D_per_level
  USE mo_ocean_math_operators, ONLY: div_oce_3D_onTriangles_onBlock, smooth_onCells, div_oce_3d, &
    & grad_fd_norm_oce_2d_onblock, grad_fd_norm_oce_2d_3d, div_oce_3D_general_onBlock, div_oce_3D_onTriangles_onBlock
  USE mo_ocean_velocity_advection, ONLY: veloc_adv_horz_mimetic, veloc_adv_vert_mimetic
  USE mo_ocean_velocity_diffusion, ONLY: velocity_diffusion, velocity_diffusion_vertical_implicit_onBlock
  USE mo_ocean_types, ONLY: t_operator_coeff, t_solverCoeff_singlePrecision
  USE mo_grid_subset, ONLY: t_subset_range, get_index_range
  USE mo_grid_config, ONLY: n_dom
  USE mo_mpi, ONLY: work_mpi_barrier
  USE mo_statistics, ONLY: global_minmaxmean, print_value_location
  USE mo_ocean_solve, ONLY: t_ocean_solve, ocean_solve_ptr
  USE mo_ocean_solve_lhs_type, ONLY: t_lhs_agen, lhs_agen_ptr
  USE mo_ocean_solve_transfer, ONLY: t_transfer, transfer_ptr
  USE mo_ocean_solve_trivial_transfer, ONLY: t_trivial_transfer, trivial_transfer_ptr
  USE mo_ocean_solve_subset_transfer, ONLY: t_subset_transfer, subset_transfer_ptr
  USE mo_ocean_solve_aux, ONLY: t_destructible, t_ocean_solve_parm, solve_gmres, solve_cg, solve_mres, &
   & ocean_solve_clear, solve_precon_none, solve_precon_jac, solve_bcgs, solve_legacy_gmres, &
   & solve_trans_scatter, solve_trans_compact, solve_cell, solve_edge
  USE mo_primal_flip_flop_lhs, ONLY: t_primal_flip_flop_lhs, lhs_primal_flip_flop_ptr
  USE mo_surface_height_lhs, ONLY: t_surface_height_lhs, lhs_surface_height_ptr

  IMPLICIT NONE
  
  PRIVATE  
  !
  PUBLIC :: solve_free_sfc_ab_mimetic
  PUBLIC :: calc_normal_velocity_ab_mimetic
  PUBLIC :: calc_vert_velocity_mim_bottomup
  PUBLIC :: invert_mass_matrix
  PUBLIC :: clear_ocean_ab_timestepping_mimetic
  !
! solver object (free ocean surface)
  CLASS(t_destructible), POINTER :: free_sfc_solver => NULL()
! solver object (free ocean surface)
  CLASS(t_destructible), POINTER :: free_sfc_solver_comp => NULL()
! left-hand-side object (free ocean surface)
  CLASS(t_destructible), POINTER :: free_sfc_solver_lhs => NULL()
! communication infrastructure object (free ocean surface)
  CLASS(t_destructible), POINTER :: free_sfc_solver_trans => NULL()
  CLASS(t_destructible), POINTER :: free_sfc_solver_comp_trans => NULL()
!  solver object (inversion of mass matrix)
  CLASS(t_destructible), POINTER :: inv_mm_solver => NULL()
! left-hand-side object (inversion of mass matrix)
  CLASS(t_destructible), POINTER :: inv_mm_solver_lhs => NULL()
! communication infrastructure object (inversion of mass matrix)
  CLASS(t_destructible), POINTER :: inv_mm_solver_trans => NULL()
  CHARACTER(LEN=12)  :: str_module = 'oceSTEPmimet' ! Output of module for 1 line debug
  INTEGER :: idt_src = 1                            ! Level of detail for 1 line debug
  REAL(wp), PARAMETER ::  min_top_height = 0.05_wp  ! we have to have at least 5cm water on topLevel of sea cells
  INTEGER, SAVE :: istep = 0

CONTAINS

! destruct and free all solver-related objects
  SUBROUTINE clear_ocean_ab_timestepping_mimetic()

    CALL ocean_solve_clear(free_sfc_solver_lhs)
    CALL ocean_solve_clear(free_sfc_solver_trans)
    CALL ocean_solve_clear(free_sfc_solver)
    CALL ocean_solve_clear(free_sfc_solver_comp)
    CALL ocean_solve_clear(free_sfc_solver_comp_trans)
    CALL ocean_solve_clear(inv_mm_solver_lhs)
    CALL ocean_solve_clear(inv_mm_solver_trans)
    CALL ocean_solve_clear(inv_mm_solver)
  END SUBROUTINE clear_ocean_ab_timestepping_mimetic

  SUBROUTINE init_free_sfc_ab_mimetic(patch_3d, ocean_state, op_coeffs, solverCoeff_sp)
    TYPE(t_patch_3d ),POINTER, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(INOUT) :: ocean_state
    TYPE(t_operator_coeff), INTENT(IN), TARGET :: op_coeffs
    TYPE(t_solverCoeff_singlePrecision), INTENT(in), TARGET :: solverCoeff_sp
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_surface_height_lhs), POINTER :: lhs_sh
    TYPE(t_trivial_transfer), POINTER :: trans_triv
    TYPE(t_subset_transfer), POINTER :: trans_subs
    TYPE(t_ocean_solve), POINTER :: solve, solve_comp
    CLASS(t_lhs_agen), POINTER :: lhs
    CLASS(t_transfer), POINTER :: trans
    TYPE(t_ocean_solve_parm) :: par, par_sp
    INTEGER :: trans_mode, sol_type
    CHARACTER(len=*), PARAMETER :: method_name='mo_ocean_ab_timestepping_mimetic:init_free_sfc_ab_mimetic'

    IF (ASSOCIATED(free_sfc_solver)) RETURN
    patch_2D => patch_3d%p_patch_2d(1)
    NULLIFY(trans_triv)
! allocate lhs object
    ALLOCATE(t_surface_height_lhs :: free_sfc_solver_lhs)
! init lhs object
    lhs_sh => lhs_surface_height_ptr(free_sfc_solver_lhs)
    CALL lhs_sh%construct(patch_3d, ocean_state%p_diag%thick_e, &
      & op_coeffs, solverCoeff_sp)
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
  END SUBROUTINE init_free_sfc_ab_mimetic

  !-------------------------------------------------------------------------
  !>
  !! !  Solves the free surface equation.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE solve_free_sfc_ab_mimetic(patch_3d, ocean_state, p_ext_data, p_as, p_oce_sfc, &
    & p_phys_param, timestep, op_coeffs, solverCoeff_sp, ret_status)
    TYPE(t_patch_3d ),POINTER, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(INOUT) :: ocean_state
    TYPE(t_external_data), TARGET, INTENT(in) :: p_ext_data
    TYPE(t_ocean_surface), INTENT(inout) :: p_oce_sfc
    TYPE(t_atmos_for_ocean), INTENT(inout) :: p_as
    TYPE (t_ho_params), INTENT(inout) :: p_phys_param
    INTEGER, INTENT(in) :: timestep
    TYPE(t_operator_coeff), INTENT(IN), TARGET :: op_coeffs
    TYPE(t_solverCoeff_singlePrecision), INTENT(in), TARGET :: solverCoeff_sp
    INTEGER :: n_it, n_it_sp, ret_status
    REAL(wp) :: rn, minmaxmean(3)
    LOGICAL :: l_is_compare_step
    CHARACTER(LEN=max_char_length) :: string
    TYPE(t_subset_range), POINTER :: owned_cells, owned_edges
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_ocean_solve), POINTER :: solve, solve_comp
    CHARACTER(len=*), PARAMETER :: method_name='mo_ocean_ab_timestepping_mimetic:solve_free_sfc_ab_mimetic'

    l_is_compare_step = .false.
    patch_2D => patch_3d%p_patch_2d(1)
    owned_cells => patch_2D%cells%owned
    owned_edges => patch_2D%edges%owned
! init sfc solver related objects, if necessary
    IF (.NOT.ASSOCIATED(free_sfc_solver)) &
      & CALL init_free_sfc_ab_mimetic(patch_3d, ocean_state, op_coeffs, solverCoeff_sp)
    solve => ocean_solve_ptr(free_sfc_solver)
    IF (l_solver_compare) solve_comp => ocean_solve_ptr(free_sfc_solver_comp)
! decide on how often solver is called (max), only relevant if GMRES-restart-alikes are chosen
    ret_status = 0
      !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('on entry: h-old'                ,ocean_state%p_prog(nold(1))%h ,str_module, 3, in_subset=owned_cells)
    CALL dbg_print('on entry: vn-old'               ,ocean_state%p_prog(nold(1))%vn,str_module, 3, in_subset=owned_edges)
    CALL dbg_print('on entry: h-new'                ,ocean_state%p_prog(nnew(1))%h ,str_module, 2, in_subset=owned_cells)
    CALL dbg_print('on entry: vn-new'               ,ocean_state%p_prog(nnew(1))%vn,str_module, 2, in_subset=owned_edges)
    ! Apply windstress
    CALL top_bound_cond_horz_veloc(patch_3d, ocean_state, op_coeffs, p_oce_sfc) ! ,     &
    start_timer(timer_ab_expl,3)
    CALL calculate_explicit_term_ab(patch_3d, ocean_state, p_phys_param, &
      & is_initial_timestep(timestep), op_coeffs, p_as)
    stop_timer(timer_ab_expl,3)
    IF(.NOT.l_rigid_lid)THEN
      ! Calculate RHS of surface equation
      start_detail_timer(timer_ab_rhs4sfc,5)
      CALL fill_rhs4surface_eq_ab(patch_3d, ocean_state, op_coeffs)
      stop_detail_timer(timer_ab_rhs4sfc,5)
      ! Solve surface equation with solver
! decide on guess to use in solver
      SELECT CASE (solver_FirstGuess)
      CASE (1)
        CALL smooth_onCells(patch_3d, ocean_state%p_prog(nold(1))%h, solve%x_loc_wp, &
          & (/ 0.5_wp, 0.5_wp /), .false., -999999.0_wp)
      CASE (2)
        !!ICON_OMP PARALLEL WORKSHARE
        solve%x_loc_wp(:,:) = ocean_state%p_prog(nold(1))%h(:,:)
        !!ICON_OMP END PARALLEL WORKSHARE
      CASE default
        !!ICON_OMP PARALLEL WORKSHARE
        solve%x_loc_wp(:,:) = 0.0_wp
        !!ICON_OMP END PARALLEL WORKSHARE
      END SELECT
      solve%b_loc_wp => ocean_state%p_aux%p_rhs_sfc_eq
      IF (l_solver_compare) THEN
        IF (istep .EQ. 0) l_is_compare_step = .true.
        istep = istep + 1
        IF (istep .GE. solver_comp_nsteps) istep = 0
      END IF
      IF (l_is_compare_step) THEN
        !!ICON_OMP PARALLEL WORKSHARE
        solve_comp%x_loc_wp(:,:) = solve%x_loc_wp(:,:)
        !!ICON_OMP END PARALLEL WORKSHARE
        solve_comp%b_loc_wp => solve%b_loc_wp
      END IF
      CALL dbg_print('bef ocean_solve('//TRIM(solve%sol_type_name)//'): h-old',ocean_state%p_prog(nold(1))%h(:,:) , &
        & str_module,idt_src,in_subset=owned_cells)
! call solver
      CALL solve%solve(n_it, n_it_sp)
      rn = MERGE(solve%res_loc_wp(1), 0._wp, n_it .NE. 0)
      ! output of sum of iterations every timestep
      IF (idbg_mxmn >= 0) THEN
        IF (n_it_sp .NE. -2) THEN
          WRITE(string,'(2(a,i4),2(a,e28.20),a)') &
            & 'SUM of ocean_solve iteration(sp,wp) = (', &
            & n_it_sp - 1, ', ', n_it - 1, ') , residual = (', &
            & solve%res_loc_wp(1), ', ', rn, ')'
        ELSE
          WRITE(string,'(a,i4,a,e28.20)') &
            &'SUM of ocean_solve iteration =', n_it - 1, ', residual =', rn
        END IF
        CALL message('ocean_solve('//TRIM(solve%sol_type_name)//'): surface height',TRIM(string))
      ENDIF
      IF (rn > solver_tolerance) THEN
        ret_status = 2
        CALL warning(method_name, "NOT YET CONVERGED !!")
        RETURN
      ENDIF
      !!ICON_OMP PARALLEL WORKSHARE
      ocean_state%p_prog(nnew(1))%h(:,:) = solve%x_loc_wp(:,:)
      !!ICON_OMP END PARALLEL WORKSHARE

      IF (l_is_compare_step) THEN
        CALL solve_comp%solve(n_it, n_it_sp)
        rn = MERGE(solve_comp%res_loc_wp(1), 0._wp, n_it .NE. 0)
        WRITE(string,'(a,i4,a,e28.20)') &
          & 'SUM of ocean_solve iteration =', n_it-1,', residual =', rn
        CALL message('ocean_solve('//TRIM(solve_comp%sol_type_name)//'): surface height',TRIM(string))
        solve_comp%x_loc_wp(:,:) = solve%x_loc_wp(:,:) - solve_comp%x_loc_wp(:,:)
        minmaxmean(:) = global_minmaxmean(values=solve_comp%x_loc_wp(:,:), in_subset=owned_cells)
        WRITE(string,"(a,3(e12.3,'  '))") "comparison of solutions: (min/max/mean)", minmaxmean(:)
        CALL message('ocean_solve('//TRIM(solve_comp%sol_type_name)//'): surface height',TRIM(string))
        solve_comp%x_loc_wp(:,:) = solve_comp%x_loc_wp(:,:) * solve_comp%x_loc_wp(:,:)
        minmaxmean(:) = global_minmaxmean(values=solve_comp%x_loc_wp(:,:), in_subset=owned_cells)
        WRITE(string,"(a,3(e12.3,'  '))") "comparison of solutions (squared): (min/max/mean)", SQRT(minmaxmean(:))
        CALL message('ocean_solve('//TRIM(solve_comp%sol_type_name)//'): surface height',TRIM(string))
      END IF

      IF (createSolverMatrix) &
        CALL solve%dump_matrix(timestep)
!        CALL createSolverMatrix_onTheFly(solve%lhs, ocean_state%p_aux%p_rhs_sfc_eq, timestep)

      !-------- end of solver ---------------
      CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_prog(nnew(1))%h)
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      ! idt_src=2  ! output print level (1-5, fix)
      !     z_h_c = lhs_surface_height( ocean_state%p_prog(nnew(1))%h, &
      !         & ocean_state%p_prog(nold(1))%h, &
      !         & patch_3d,             &
      !         & z_implcoeff,            &
      !         & ocean_state%p_diag%thick_e,    &
      !         & ocean_state%p_diag%thick_c,    &
      !         & op_coeffs)             &
      !         & -ocean_state%p_aux%p_rhs_sfc_eq
      !     CALL dbg_print('SolvSfc: residual h-res'    ,z_h_c                  ,str_module,idt_src)
!      vol_h(:,:) = patch_3d%p_patch_2d(n_dom)%cells%area(:,:) * ocean_state%p_prog(nnew(1))%h(:,:)
!      CALL dbg_print('after ocean_solve: vol_h(:,:)',vol_h ,str_module,idt_src, in_subset=owned_cells)
      !---------------------------------------------------------------------
      CALL dbg_print('vn-new',ocean_state%p_prog(nnew(1))%vn,str_module, 2,in_subset=owned_edges)
      CALL dbg_print('aft ocean_solve('//TRIM(solve%sol_type_name)//'): h-new', &
        & ocean_state%p_prog(nnew(1))%h(:,:) ,str_module,2,in_subset=owned_cells)
      minmaxmean(:) = global_minmaxmean(values=ocean_state%p_prog(nnew(1))%h(:,:), in_subset=owned_cells)
      CALL debug_print_MaxMinMean('after ocean_solve('//TRIM(solve%sol_type_name)//'): h-new', &
        & minmaxmean, str_module, 2)
      IF (minmaxmean(1) + patch_3D%p_patch_1D(1)%del_zlev_m(1) <= min_top_height) THEN
!          CALL finish(method_name, "height below min_top_height")
        CALL warning(method_name, "height below min_top_height")
        CALL print_value_location(ocean_state%p_prog(nnew(1))%h(:,:), minmaxmean(1), owned_cells)
        CALL print_value_location(ocean_state%p_prog(nnew(2))%h(:,:), minmaxmean(1), owned_cells)
        CALL work_mpi_barrier()
        ret_status = 1
        RETURN
      ENDIF
    ENDIF  ! l_rigid_lid
  END SUBROUTINE solve_free_sfc_ab_mimetic
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  !>
  !! Computation of velocity predictor in Adams-Bashforth timestepping.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE calculate_explicit_term_ab( patch_3d, ocean_state, p_phys_param,&
    & is_first_timestep, op_coeffs, p_as)
    TYPE(t_patch_3d ), POINTER, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE (t_ho_params)                   :: p_phys_param
    LOGICAL,INTENT(in)                   :: is_first_timestep
    TYPE(t_operator_coeff), INTENT(IN), TARGET :: op_coeffs
    TYPE(t_atmos_for_ocean), INTENT(inout) :: p_as
    TYPE(t_subset_range), POINTER :: owned_edges, owned_cells
    
    owned_edges     => patch_3d%p_patch_2d(n_dom)%edges%owned
    owned_cells     => patch_3d%p_patch_2d(n_dom)%cells%owned

    ! STEP 1: horizontal advection
    start_detail_timer(timer_extra1,4)
    IF(is_first_timestep)THEN
      CALL veloc_adv_horz_mimetic( patch_3d,         &
        & ocean_state%p_prog(nold(1))%vn,    &
        & ocean_state%p_prog(nold(1))%vn,    &
        & ocean_state%p_diag,                &
        & ocean_state%p_diag%veloc_adv_horz, &
        & op_coeffs)
    ELSE
      CALL veloc_adv_horz_mimetic( patch_3d,         &
        & ocean_state%p_prog(nold(1))%vn,    &
        & ocean_state%p_prog(nnew(1))%vn,    &
        & ocean_state%p_diag,                &
        & ocean_state%p_diag%veloc_adv_horz, &
        & op_coeffs)
    ENDIF
    stop_detail_timer(timer_extra1,4)
    ! STEP 2: compute 3D contributions: gradient of hydrostatic pressure and vertical velocity advection
    IF ( iswm_oce /= 1 ) THEN
      ! calculate density from EOS using temperature and salinity at timelevel n
      start_detail_timer(timer_extra2,4)
      CALL calculate_density( patch_3d,                         &
       & ocean_state%p_prog(nold(1))%tracer(:,:,:,1:no_tracer),&
       & ocean_state%p_diag%rho(:,:,:) )

!       IF(.NOT.l_partial_cells)THEN      
!         ! calculate hydrostatic pressure from density at timelevel nc
!         CALL calc_internal_press( patch_3d,                  &  ! in
!           & ocean_state%p_diag%rho,                          &  ! in
!           & patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c,   &  ! in
!           & ocean_state%p_prog(nold(1))%h,                   &  ! in
!           & ocean_state%p_diag%press_hyd)                       ! inout
!       
!         ! calculate gradient of hydrostatic pressure in 3D
!         CALL grad_fd_norm_oce_3d(               &
!           & ocean_state%p_diag%press_hyd,       &
!           & patch_3d,                           &
!           & op_coeffs%grad_coeff,               &
!           & ocean_state%p_diag%press_grad)
!       ELSE
!          CALL calc_internal_press_grad0( patch_3d,&
!          &                              ocean_state%p_diag%rho,&
!          &                              op_coeffs%grad_coeff,  &
!          &                              ocean_state%p_diag%press_grad)
!       ENDIF




     CALL calc_internal_press_grad( patch_3d,&
         &                          ocean_state%p_diag%rho,&
         &                          ocean_state%p_diag%press_hyd,& 
         &                          p_as%pao,&
         &                          op_coeffs%grad_coeff,  &
         &                          ocean_state%p_diag%press_grad)     
      ! calculate vertical velocity advection
      CALL veloc_adv_vert_mimetic(          &
        & patch_3d,                         &
        & ocean_state%p_diag,op_coeffs,     &
        & ocean_state%p_diag%veloc_adv_vert )
      stop_detail_timer(timer_extra2,4)
    ELSE  !  iswm_oce=1
      ocean_state%p_diag%veloc_adv_vert = 0.0_wp
      ocean_state%p_diag%laplacian_vert = 0.0_wp
    ENDIF
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src = 3  ! output print level (1-5, fix)
    CALL dbg_print('density'                   ,ocean_state%p_diag%rho           ,str_module,idt_src, &
      & in_subset = owned_cells)
    CALL dbg_print('internal pressure'         ,ocean_state%p_diag%press_hyd     ,str_module,idt_src, &
      in_subset = owned_cells)
    CALL dbg_print('internal press grad'       ,ocean_state%p_diag%press_grad    ,str_module,idt_src, &
      in_subset = owned_edges)
    idt_src = 4  ! output print level (1-5, fix)
    CALL dbg_print('kinetic energy'            ,ocean_state%p_diag%kin           ,str_module,idt_src, &
      in_subset = owned_cells)
    ! STEP 3: compute harmonic or biharmoic laplacian diffusion of velocity.
    !         This term is discretized explicitly. Order and form of the laplacian
    !         are determined in mo_oce_diffusion according to namelist settings
    start_detail_timer(timer_extra3,5)
    CALL velocity_diffusion(patch_3d,              &
      & ocean_state%p_prog(nold(1))%vn, &
      & p_phys_param,            &
      & ocean_state%p_diag,op_coeffs,  &
      & ocean_state%p_diag%laplacian_horz)
    stop_detail_timer(timer_extra3,5)
    start_detail_timer(timer_extra4,4)
    IF (   MASS_MATRIX_INVERSION_TYPE == MASS_MATRIX_INVERSION_ALLTERMS&
      &.OR.MASS_MATRIX_INVERSION_TYPE == MASS_MATRIX_INVERSION_ADVECTION) THEN
      CALL explicit_vn_pred_invert_mass_matrix( patch_3d, ocean_state, op_coeffs, p_phys_param, is_first_timestep)
    ELSE ! If no inversion of mass matrix
      CALL explicit_vn_pred( patch_3d, ocean_state, op_coeffs, p_phys_param, is_first_timestep)
    ENDIF!
    stop_detail_timer(timer_extra4,4)
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    idt_src = 4  ! output print level (1-5, fix)
    CALL dbg_print('bc_top_vn'   ,ocean_state%p_aux%bc_top_vn,str_module,idt_src, in_subset=owned_edges)
    CALL dbg_print('horizontal advection'      ,ocean_state%p_diag%veloc_adv_horz,str_module,idt_src, &
      & in_subset=owned_edges)
    CALL dbg_print('horizontal grad'          ,ocean_state%p_diag%grad,str_module,idt_src, &
      & in_subset=owned_edges)
    CALL dbg_print('vertical advection'        ,ocean_state%p_diag%veloc_adv_vert,str_module,idt_src, &
      & in_subset=owned_edges)
    CALL dbg_print('VelocDiff: LaPlacHorz'    ,ocean_state%p_diag%laplacian_horz  ,str_module,idt_src, in_subset=owned_edges)
    IF (iswm_oce /= 1) THEN
      CALL dbg_print('vn_pred'   ,ocean_state%p_diag%vn_pred       ,str_module,2, in_subset=owned_edges)
    ELSE
      CALL dbg_print('VelocDiff: LaPlacVert'    ,ocean_state%p_diag%laplacian_vert,str_module,idt_src, in_subset=owned_edges)
    ENDIF
    idt_src=5  ! output print level (1-5, fix)
    CALL dbg_print('vn(nold)'                 ,ocean_state%p_prog(nold(1))%vn ,str_module,idt_src, in_subset=owned_edges)
    CALL dbg_print('G_n+1/2 - g_nimd'         ,ocean_state%p_aux%g_nimd       ,str_module,idt_src, in_subset=owned_edges)
    CALL dbg_print('G_n'                      ,ocean_state%p_aux%g_n          ,str_module,idt_src, in_subset=owned_edges)
    CALL dbg_print('G_n-1'                    ,ocean_state%p_aux%g_nm1        ,str_module,idt_src, in_subset=owned_edges)
  END SUBROUTINE calculate_explicit_term_ab

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE explicit_vn_pred( patch_3d, ocean_state, op_coeffs, p_phys_param, is_first_timestep)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_operator_coeff), INTENT(IN) :: op_coeffs
    TYPE (t_ho_params) :: p_phys_param
    LOGICAL, INTENT(in) :: is_first_timestep
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
      CALL grad_fd_norm_oce_2d_onBlock( ocean_state%p_prog(nold(1))%h, patch_2D, &
        & op_coeffs%grad_coeff(:,1,blockNo), z_gradh_e(:), start_edge_index, end_edge_index, blockNo)
      z_gradh_e(start_edge_index:end_edge_index) = &
         & (1.0_wp-ab_beta) * grav * z_gradh_e(start_edge_index:end_edge_index)
      ! STEP 5:
      CALL calculate_explicit_term_g_n_onBlock( patch_3d, ocean_state, is_first_timestep, &
        & start_edge_index, end_edge_index, blockNo)
      IF ( iswm_oce /= 1) THEN
        CALL calculate_explicit_vn_pred_3D_onBlock( patch_3d, ocean_state, z_gradh_e(:),    &
        & start_edge_index, end_edge_index, blockNo)
        ! calculate vertical friction, ie p_phys_param%a_veloc_v
        IF (PPscheme_type == PPscheme_ICON_Edge_vnPredict_type) &
          CALL ICON_PP_Edge_vnPredict_scheme(patch_3d, blockNo, start_edge_index, end_edge_index, &
            & ocean_state, ocean_state%p_diag%vn_pred(:,:,blockNo))
        !In 3D case implicit vertical velocity diffusion is chosen
        CALL velocity_diffusion_vertical_implicit_onBlock(patch_3d, &
          & ocean_state%p_diag%vn_pred(:,:,blockNo), p_phys_param%a_veloc_v(:,:,blockNo), op_coeffs, &
          & start_edge_index, end_edge_index, blockNo)
      ELSE !( iswm_oce == 1) THEN
        CALL calculate_explicit_vn_pred_2D_onBlock( patch_3d, ocean_state, z_gradh_e(:), &
        & start_edge_index, end_edge_index, blockNo) 
      END IF
    END DO
!ICON_OMP_END_PARALLEL_DO
  END SUBROUTINE explicit_vn_pred
  
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE   explicit_vn_pred_invert_mass_matrix( patch_3d, ocean_state, op_coeffs, p_phys_param, is_first_timestep)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_operator_coeff), INTENT(IN) :: op_coeffs
    TYPE (t_ho_params) :: p_phys_param
    LOGICAL,INTENT(in) :: is_first_timestep
    REAL(wp) :: z_gradh_e(nproma)
    TYPE(t_subset_range), POINTER :: edges_in_domain
    INTEGER :: start_edge_index, end_edge_index, blockNo
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: z_e(nproma,n_zlev,patch_3d%p_patch_2d(n_dom)%nblks_e)
    !-----------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(n_dom)
    edges_in_domain => patch_3d%p_patch_2d(n_dom)%edges%in_domain
    IF (MASS_MATRIX_INVERSION_TYPE== MASS_MATRIX_INVERSION_ADVECTION)THEN
      !Here the inversion of the mass matrix is already carried out
!ICON_OMP PARALLEL WORKSHARE
      z_e(:,:,:) = ocean_state%p_diag%veloc_adv_horz(:,:,:) + ocean_state%p_diag%veloc_adv_vert(:,:,:)
!ICON_OMP END PARALLEL WORKSHARE
       WRITE(0,*)'ADV before:', maxval(z_e(:,1,:)),minval(z_e(:,1,:))
      ocean_state%p_diag%veloc_adv_horz = invert_mass_matrix(patch_3d, ocean_state, op_coeffs, z_e)
      CALL sync_patch_array(sync_e, patch_2D, z_e)
       WRITE(0,*)'ADV after:', maxval(ocean_state%p_diag%veloc_adv_horz(:,1,:)), &
         & minval(ocean_state%p_diag%veloc_adv_horz(:,1,:))
    END IF
    ! STEP 4: calculate weighted gradient of surface height at previous timestep
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, z_gradh_e) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      z_gradh_e(:)  = 0.0_wp
      CALL grad_fd_norm_oce_2d_onBlock( ocean_state%p_prog(nold(1))%h, patch_2D, &
        & op_coeffs%grad_coeff(:,1,blockNo), z_gradh_e(:), start_edge_index, end_edge_index, blockNo)
      z_gradh_e(start_edge_index:end_edge_index) = &
         & (1.0_wp-ab_beta) * grav * z_gradh_e(start_edge_index:end_edge_index)
      ! STEP 5:
      CALL calculate_explicit_term_g_n_onBlock( patch_3d, ocean_state, is_first_timestep, &
        & start_edge_index, end_edge_index, blockNo)
      IF ( iswm_oce /= 1) THEN
        CALL calculate_explicit_vn_pred_3D_onBlock( patch_3d, ocean_state, z_gradh_e(:),    &
        & start_edge_index, end_edge_index, blockNo)
        !In 3D case implicit vertical velocity diffusion is chosen
        CALL velocity_diffusion_vertical_implicit_onBlock( patch_3d, &
          & ocean_state%p_diag%vn_pred(:,:,blockNo), p_phys_param%a_veloc_v(:,:,blockNo), &
          & op_coeffs, start_edge_index, end_edge_index, blockNo)
      ELSE !( iswm_oce == 1)
        CALL calculate_explicit_vn_pred_2D_onBlock( patch_3d, ocean_state, z_gradh_e(:), &
        & start_edge_index, end_edge_index, blockNo) 
      END IF      
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    IF (MASS_MATRIX_INVERSION_TYPE == MASS_MATRIX_INVERSION_ALLTERMS )THEN
      !Here the inversion is just prepared
      WRITE(0,*)'vn_pred before:', maxval(ocean_state%p_diag%vn_pred(:,1,:)), &
        & minval(ocean_state%p_diag%vn_pred(:,1,:))
      CALL map_edges2edges_viacell_2D_per_level( patch_3d, &
        & ocean_state%p_diag%vn_pred(:,1,:), op_coeffs, &
        & ocean_state%p_diag%vn_pred_ptp(:,1,:),1 )
      WRITE(0,*)'vn_pred after:', maxval(ocean_state%p_diag%vn_pred_ptp(:,1,:)), &
        & minval(ocean_state%p_diag%vn_pred_ptp(:,1,:))
      ocean_state%p_diag%vn_pred=ocean_state%p_diag%vn_pred_ptp       
    END IF
  END SUBROUTINE   explicit_vn_pred_invert_mass_matrix

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE calculate_explicit_vn_pred_3D_onBlock( patch_3d, ocean_state, z_gradh_e, &
    & start_edge_index, end_edge_index, blockNo)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    REAL(wp) :: z_gradh_e(nproma)
    INTEGER, INTENT(in) :: start_edge_index, end_edge_index, blockNo
    INTEGER :: je, jk, bottom_level
    !-----------------------------------------------------------------------
    IF(.NOT.l_rigid_lid)THEN
        DO je = start_edge_index, end_edge_index
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
            ocean_state%p_diag%vn_pred(je,jk,blockNo) = ocean_state%p_prog(nold(1))%vn(je,jk,blockNo)  &
              & + dtime*(ocean_state%p_aux%g_nimd(je,jk,blockNo) &
              & - z_gradh_e(je))
          END DO
        END DO
    ELSE ! IF(l_rigid_lid)THEN
        DO je = start_edge_index, end_edge_index
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
            ocean_state%p_diag%vn_pred(je,jk,blockNo) = ocean_state%p_prog(nold(1))%vn(je,jk,blockNo)&
              & + dtime*ocean_state%p_aux%g_nimd(je,jk,blockNo)
        END DO
      END DO
    ENDIF!Rigid lid
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
  END SUBROUTINE calculate_explicit_vn_pred_3D_onBlock

  !-------------------------------------------------------------------------
  SUBROUTINE calculate_explicit_vn_pred_2D_onBlock( patch_3d, ocean_state, z_gradh_e, &
    &  start_edge_index, end_edge_index, blockNo)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    REAL(wp),INTENT(in)                  :: z_gradh_e(nproma)
    INTEGER, INTENT(in)                  :: start_edge_index, end_edge_index, blockNo
    INTEGER :: je, jk
    !-----------------------------------------------------------------------
    DO je = start_edge_index, end_edge_index
      DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
        ocean_state%p_diag%vn_pred(je,jk,blockNo) = (ocean_state%p_prog(nold(1))%vn(je,jk,blockNo)        &
          & + dtime*(ocean_state%p_aux%g_nimd(je,jk,blockNo)        &
          & - z_gradh_e(je)))
      END DO
    END DO
    !In case of Shallow-water with forcing and or damping
    IF ( iforc_oce /= No_Forcing) THEN
      DO je = start_edge_index, end_edge_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          ocean_state%p_diag%vn_pred(je,jk,blockNo) = (ocean_state%p_diag%vn_pred(je,jk,blockNo) &
            & + ocean_state%p_aux%bc_top_vn(je,blockNo)    &
            & - ocean_state%p_aux%bc_bot_vn(je,blockNo))
        END DO
      END DO
    ENDIF
  END SUBROUTINE calculate_explicit_vn_pred_2D_onBlock

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE calculate_explicit_term_g_n_onBlock( patch_3d, ocean_state, is_first_timestep, &
    & start_edge_index, end_edge_index, blockNo)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    LOGICAL,INTENT(in)                   :: is_first_timestep
    INTEGER,INTENT(in)                   :: start_edge_index, end_edge_index, blockNo
    INTEGER :: je, jk
    !-----------------------------------------------------------------------
    IF(MASS_MATRIX_INVERSION_TYPE/=MASS_MATRIX_INVERSION_ADVECTION)THEN
      DO je = start_edge_index, end_edge_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          ocean_state%p_aux%g_n(je, jk, blockNo) = &
            & - ocean_state%p_diag%press_grad    (je, jk, blockNo)  &
            & - ocean_state%p_diag%grad          (je, jk, blockNo)  &            
            & - ocean_state%p_diag%veloc_adv_horz(je, jk, blockNo)  &
            & - ocean_state%p_diag%veloc_adv_vert(je, jk, blockNo)  &
            & + ocean_state%p_diag%laplacian_horz(je, jk, blockNo)  !&
            !& + ocean_state%p_diag%laplacian_vert(je, jk, blockNo) !<- is done later implicitely 
        END DO
      END DO
    ELSEIF(MASS_MATRIX_INVERSION_TYPE==MASS_MATRIX_INVERSION_ADVECTION)THEN
      DO je = start_edge_index, end_edge_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          ocean_state%p_aux%g_n(je, jk, blockNo) = &
            & - ocean_state%p_diag%press_grad    (je, jk, blockNo)  &
            & - ocean_state%p_diag%grad          (je, jk, blockNo)  &            
            & - ocean_state%p_diag%veloc_adv_horz(je, jk, blockNo)  &
   !         & - ocean_state%p_diag%veloc_adv_vert(je, jk, blockNo)  &
            & + ocean_state%p_diag%laplacian_horz(je, jk, blockNo)  !&
            !& + ocean_state%p_diag%laplacian_vert(je, jk, blockNo) !<- is done later implicitely 
        END DO
      END DO
    ENDIF
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
  END SUBROUTINE calculate_explicit_term_g_n_onBlock
 
  !-------------------------------------------------------------------------
  !>
  !!  Calculation of right-hand side of elliptic surface equation.
  !!  This is used in semi implicit timelevel stepping.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE fill_rhs4surface_eq_ab( patch_3d, ocean_state, op_coeffs)
    ! Patch on which computation is performed
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_operator_coeff), INTENT(IN) :: op_coeffs
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
    IF(iswm_oce == 1)THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, je, jk) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
        z_vn_ab(:,:,blockNo)  = 0.0_wp
        DO je = start_edge_index, end_edge_index
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
            z_vn_ab(je,jk,blockNo)=ab_gam*ocean_state%p_diag%vn_pred(je,jk,blockNo)&
              & + (1.0_wp -ab_gam)* ocean_state%p_prog(nold(1))%vn(je,jk,blockNo)
          END DO
        ENDDO
      END DO
!ICON_OMP_END_PARALLEL_DO
    ELSE! IF(iswm_oce /= 1)THEN
!       IF(expl_vertical_velocity_diff==1)THEN
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
     ENDIF
     z_vn_ab(:,:,edges_in_domain%end_block+1:) = 0._wp
    !
    ! calculate depth-integrated velocity z_e
    !  - edge-based and cell-based
    !  - 3d and 2d (surface)
    !-------------------------------------------------------------------------------
    IF (l_edge_based) THEN
      !-------------------------------------------------------------------------------
      z_e(1:nproma,1:n_zlev,1:patch_2D%nblks_e) = 0.0_wp
      IF( iswm_oce /= 1 ) THEN !the 3D case
        DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
          DO je = start_edge_index, end_edge_index
!DIR$ SIMD
            DO jk=1,patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
              z_e(je,jk,blockNo)=  z_vn_ab(je,jk,blockNo) * patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,blockNo)
            END DO
          ENDDO
        END DO
      ELSEIF( iswm_oce == 1 ) THEN
        DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
          DO je = start_edge_index, end_edge_index
            DO jk=1,patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
              z_e(je,jk,blockNo)= z_vn_ab(je,jk,blockNo)*ocean_state%p_diag%thick_e(je,blockNo)
            END DO
          ENDDO
        END DO
      ENDIF
      ! !-------------------------------------------------------------------------------
    ELSE ! IF(.NOT. l_edge_based)THEN!NOT EDGE-BASED
      ! !-------------------------------------------------------------------------------
      CALL sync_patch_array(sync_e, patch_2D, z_vn_ab)
      IF( iswm_oce /= 1 ) THEN !the 3D case
        CALL map_edges2edges_viacell_3d_const_z( patch_3d, z_vn_ab, op_coeffs, z_e )
      ELSE ! IF( iswm_oce == 1 ) THEN
        CALL map_edges2edges_viacell_3d_const_z( patch_3d, z_vn_ab(:,1,:), op_coeffs, z_e(:,1,:) )
      ENDIF!( iswm_oce == 1 )
    ENDIF!EDGE-BASED
    IF( patch_2d%cells%max_connectivity == 3 )THEN
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
            ocean_state%p_aux%p_rhs_sfc_eq(jc,blockNo) = ((ocean_state%p_prog(nold(1))%h(jc,blockNo) &
              & - dtime * div_z_depth_int_c(jc)) * inv_gdt2)
          ENDIF
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO
    ELSE
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, jc, jk, div_z_depth_int_c, div_z_c) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
        CALL div_oce_3D_general_onBlock( z_e, patch_3D, op_coeffs%div_coeff, div_z_c, &
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
          ocean_state%p_aux%p_rhs_sfc_eq(jc,blockNo) = ((ocean_state%p_prog(nold(1))%h(jc,blockNo) &
          & - dtime * div_z_depth_int_c(jc)) * inv_gdt2)
        ENDIF
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO    
  ENDIF!patch_2d%cells%max_connectivity
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('RHS thick_e'               ,patch_3d%p_patch_1d(1)%prism_thick_e, str_module,idt_src, &
        & in_subset=patch_3d%p_patch_2d(1)%edges%owned )
    CALL dbg_print('RHS thick_c'               ,patch_3d%p_patch_1d(1)%prism_thick_c, str_module,idt_src, &
        & in_subset=patch_3d%p_patch_2d(1)%cells%owned)        
    CALL dbg_print('RHS z_vn_ab'               ,z_vn_ab                  ,str_module,idt_src, &
        & in_subset=owned_edges )
    CALL dbg_print('RHS z_e'                   ,z_e                      ,str_module,idt_src, &
        & in_subset=owned_edges )
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('RHS final'                 ,ocean_state%p_aux%p_rhs_sfc_eq  ,str_module,idt_src, &
      in_subset=patch_3d%p_patch_2d(1)%cells%owned)
  END SUBROUTINE fill_rhs4surface_eq_ab
  
  !-------------------------------------------------------------------------
  !>
  !! Computation of new velocity in Adams-Bashforth timestepping.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE calc_normal_velocity_ab_mimetic(patch_3d,ocean_state, op_coeffs)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs
    !  local variables
    INTEGER  :: start_edge_index, end_edge_index, je, jk, blockNo
    REAL(wp) :: gdt_x_ab_beta, one_minus_ab_gam
    REAL(wp) :: z_grad_h(nproma,patch_3d%p_patch_2d(1)%nblks_e), z_grad_h_block(nproma)
    TYPE(t_subset_range), POINTER :: edges_in_domain, owned_edges
    CHARACTER(LEN=*), PARAMETER ::     &
      & method_name='mo_ocean_ab_timestepping_mimetic: calc_normal_velocity_ab_mimetic'
    TYPE(t_patch), POINTER :: patch
    !----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    !-----------------------------------------------------------------------
    patch           => patch_3d%p_patch_2d(1)
    edges_in_domain => patch%edges%in_domain
    owned_edges     => patch%edges%owned
    one_minus_ab_gam = 1.0_wp - ab_gam
    gdt_x_ab_beta =  grav * dtime * ab_beta
    IF (iswm_oce == 1) THEN ! shallow water case
      z_grad_h(1:nproma,1:patch%nblks_e) = 0.0_wp
      ! Step 1) Compute normal derivative of new surface height
      CALL grad_fd_norm_oce_2d_3d( ocean_state%p_prog(nnew(1))%h, patch, &
        & op_coeffs%grad_coeff(:,1,:), z_grad_h(:,:))
      IF (MASS_MATRIX_INVERSION_TYPE .EQ. MASS_MATRIX_INVERSION_ALLTERMS) THEN
        ! Step 2) Calculate the new velocity from the predicted one and the new surface height
        DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
          DO jk = 1, n_zlev
            DO je = start_edge_index, end_edge_index
              IF(patch_3d%lsm_e(je,jk,blockNo) <= sea_boundary)THEN
                ocean_state%p_prog(nnew(1))%vn(je,jk,blockNo) = ocean_state%p_diag%vn_pred_ptp(je,jk,blockNo) &
                  & - gdt_x_ab_beta * z_grad_h(je,blockNo)
              END IF
            END DO
          END DO
        END DO
        WRITE(0,*)'New velocity before Inversion', maxval(ocean_state%p_prog(nnew(1))%vn(:,1,:)), &
          & minval(ocean_state%p_prog(nnew(1))%vn(:,1,:)), maxval(ocean_state%p_diag%vn_pred_ptp(:,1,:)), &
          & minval(ocean_state%p_diag%vn_pred_ptp(:,1,:))
         ocean_state%p_prog(nnew(1))%vn = invert_mass_matrix(patch_3d, ocean_state, op_coeffs, ocean_state%p_prog(nnew(1))%vn)
         WRITE(0,*)'New velocity after Inversion', maxval(ocean_state%p_prog(nnew(1))%vn(:,1,:)), &
           & minval(ocean_state%p_prog(nnew(1))%vn(:,1,:))
      ELSE
        ! Step 2) Calculate the new velocity from the predicted one and the new surface height
        DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
          DO jk = 1, n_zlev
            DO je = start_edge_index, end_edge_index
              IF (patch_3d%lsm_e(je,jk,blockNo) .LE. sea_boundary) THEN
                ocean_state%p_prog(nnew(1))%vn(je,jk,blockNo) = ocean_state%p_diag%vn_pred(je,jk,blockNo) &
                  & - gdt_x_ab_beta * z_grad_h(je,blockNo)
              END IF
            END DO
          END DO
        END DO
      END IF
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
        DO je = start_edge_index, end_edge_index
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
            ocean_state%p_diag%vn_time_weighted(je,jk,blockNo) = ab_gam*ocean_state%p_prog(nnew(1))%vn(je,jk,blockNo) &
              & + (1.0_wp -ab_gam)*ocean_state%p_prog(nold(1))%vn(je,jk,blockNo)          
          END DO
        END DO
      END DO
    ELSE !real 3d case
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,start_edge_index,end_edge_index, je, jk, z_grad_h_block) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
        ! Step 1) Compute normal derivative of new surface height
        CALL grad_fd_norm_oce_2d_onBlock(ocean_state%p_prog(nnew(1))%h, patch, op_coeffs%grad_coeff(:,1, blockNo), &
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
    ENDIF    
    CALL sync_patch_array_mult(sync_e, patch, 2, ocean_state%p_prog(nnew(1))%vn, ocean_state%p_diag%vn_time_weighted)
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('NormVel: vn_old'            ,ocean_state%p_prog(nold(1))%vn     ,str_module,idt_src, in_subset=owned_edges )
    CALL dbg_print('NormVel: vn_pred'           ,ocean_state%p_diag%vn_pred         ,str_module,idt_src, in_subset=owned_edges)
    CALL dbg_print('NormVel: vn_time_weighted'  ,ocean_state%p_diag%vn_time_weighted,str_module,idt_src, in_subset=owned_edges)
    CALL dbg_print('NormVel: vn_change'         ,ocean_state%p_prog(nnew(1))%vn - &
      & ocean_state%p_prog(nold(1))%vn     ,str_module,idt_src, in_subset=owned_edges)
    idt_src=2  ! outputm print level (1-5, fix)
    CALL dbg_print('NormVel: vn_new'            ,ocean_state%p_prog(nnew(1))%vn     ,str_module,idt_src, in_subset=owned_edges)
  END SUBROUTINE calc_normal_velocity_ab_mimetic
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
  !! @par Revision History
  !! Developed  by  Peter Korn,   MPI-M (2006).
  !!  Modified by Stephan Lorenz, MPI-M (2010-06)
  !TODO review
!<Optimize:inUse>
  SUBROUTINE calc_vert_velocity_mim_bottomup( patch_3d, ocean_state, op_coeffs )
    TYPE(t_patch_3d), TARGET :: patch_3d       ! patch on which computation is performed
    TYPE(t_hydro_ocean_state) :: ocean_state
    TYPE(t_operator_coeff), INTENT(in) :: op_coeffs
    ! Local variables
    INTEGER :: jc, jk, blockNo, je, z_dolic, start_index, end_index
    REAL(wp) :: z_c(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks), z_abort
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells, cells_owned
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp), POINTER :: vertical_velocity(:,:,:)
    CHARACTER(len=*), PARAMETER :: method_name='mo_ocean_ab_timestepping_mimetic:alc_vert_velocity_mim_bottomup'
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
    IF( l_edge_based )THEN
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
        DO jk = 1, n_zlev
          DO je = start_index, end_index
            ocean_state%p_diag%mass_flx_e(je,jk,blockNo) = ocean_state%p_diag%vn_time_weighted(je,jk,blockNo)&
              & * patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,blockNo)
          END DO
        END DO
      END DO
      CALL div_oce_3d( ocean_state%p_diag%mass_flx_e,    &
        & patch_3D,                   &
        & op_coeffs%div_coeff,      &
        & ocean_state%p_diag%div_mass_flx_c)
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
        DO jc = start_index, end_index
          z_dolic = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
          !use bottom boundary condition for vertical velocity at bottom of prism
          vertical_velocity(jc,z_dolic+1,blockNo)=0.0_wp
          DO jk = z_dolic, 1, -1
            vertical_velocity(jc,jk,blockNo) = vertical_velocity(jc,jk+1,blockNo) - ocean_state%p_diag%div_mass_flx_c(jc,jk,blockNo)
          END DO
        END DO
      END DO
    !-------------------------------------------------------------------------------
    ELSE ! IF(.NOT.l_edge_based)THEN
      CALL map_edges2edges_viacell_3d_const_z( patch_3d, ocean_state%p_diag%vn_time_weighted, op_coeffs, &
        & ocean_state%p_diag%mass_flx_e)!, patch_3D%p_patch_1D(1)%prism_thick_c)
      IF ( patch_2d%cells%max_connectivity == 3 )THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
          CALL div_oce_3D_onTriangles_onBlock(ocean_state%p_diag%mass_flx_e, patch_3D, op_coeffs%div_coeff, &
            & ocean_state%p_diag%div_mass_flx_c(:,:,blockNo), blockNo=blockNo, start_index=start_index, &
            & end_index=end_index, start_level=1, end_level=n_zlev)
          DO jc = start_index, end_index
            !use bottom boundary condition for vertical velocity at bottom of prism
            ! this should be awlays zero
            ! vertical_velocity(jc,z_dolic+1,blockNo)=0.0_wp
            DO jk = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), 1, -1
              vertical_velocity(jc,jk,blockNo) = vertical_velocity(jc,jk+1,blockNo) - &
                & ocean_state%p_diag%div_mass_flx_c(jc,jk,blockNo)
            END DO
          END DO
        END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO
      ELSE
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
          CALL div_oce_3D_general_onBlock(ocean_state%p_diag%mass_flx_e, patch_3D,op_coeffs%div_coeff, &
            & ocean_state%p_diag%div_mass_flx_c(:,:,blockNo), blockNo=blockNo, start_index=start_index, &
            & end_index=end_index, start_level=1, end_level=n_zlev)
          DO jc = start_index, end_index
            !use bottom boundary condition for vertical velocity at bottom of prism
            DO jk = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), 1, -1
              vertical_velocity(jc,jk,blockNo) = vertical_velocity(jc,jk+1,blockNo) - &
                & ocean_state%p_diag%div_mass_flx_c(jc,jk,blockNo)
            END DO
          END DO
        END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO
      ENDIF !patch_2d%cells%max_connectivity
    ENDIF  !  (l_EDGE_BASED)
    IF(l_rigid_lid)THEN
      vertical_velocity(:,1,:) = 0.0_wp
    ENDIF
    CALL sync_patch_array(sync_c,patch_2D,vertical_velocity)

    !-----------------------------------------------------
    IF (use_continuity_correction) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_index, end_index)
        DO jc = start_index, end_index
          ocean_state%p_prog(nnew(1))%h(jc,blockNo) = ocean_state%p_prog(nold(1))%h(jc,blockNo) + &
          &  vertical_velocity(jc,1,blockNo) * dtime
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO
      !---------------------------------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      ! slo - cells_owned for correct global mean
      CALL dbg_print('Vert veloc: w', vertical_velocity, str_module,idt_src, in_subset=cells_owned)
      CALL dbg_print('after cont-correct: h-new',ocean_state%p_prog(nnew(1))%h(:,:) ,str_module,idt_src, &
        & in_subset=cells_owned)
      CALL dbg_print('after cont-correct: vol_h', &
        & patch_3d%p_patch_2d(n_dom)%cells%area(:,:) * ocean_state%p_prog(nnew(1))%h(:,:), &
        & str_module,idt_src, in_subset=cells_owned)
    ENDIF
    !---------------------------------------------------------------------
    IF (debug_check_level > 8) THEN
      z_c          (1:nproma,1:patch_2D%alloc_cell_blocks)   = 0.0_wp
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
        DO jc = start_index, end_index
          z_c(jc,blockNo) = ((ocean_state%p_prog(nnew(1))%h(jc,blockNo)-ocean_state%p_prog(nold(1))%h(jc,blockNo))/dtime - &
            & vertical_velocity(jc,1,blockNo)) * patch_3d%wet_c(jc,1,blockNo)
        ENDDO
      ENDDO
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=1  ! output print level (1-5, fix)
      CALL dbg_print('CalcVertVelMimBU: d_h/dt-w',z_c, &
        & str_module,idt_src, in_subset=cells_in_domain)
      idt_src=2  ! output print level (1-5, fix)
      CALL dbg_print('CalcVertVelMimBU: vertical_velocity =W' ,vertical_velocity, &
        & str_module,idt_src, in_subset=cells_in_domain)
      idt_src=4  ! output print level (1-5, fix)
      CALL dbg_print('CalcVertVelMimBU: mass flx',ocean_state%p_diag%mass_flx_e,   &
        &  str_module,idt_src, in_subset=patch_2D%edges%owned)
      CALL dbg_print('CalcVertVelMimBU: div mass',ocean_state%p_diag%div_mass_flx_c,&
        & str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------
      ! Abort if largest mismatch in surface elevation due to solution of ocean_gmres-solver is > 1mm/year
      !   criterion is 1mm/year * dtime = 3.17e-11 m/s * dtime
      z_abort = dhdtw_abort*dtime
      IF (MAXVAL(ABS(z_c(:,:))) > z_abort) THEN
        CALL message('mo_ocean_ab_timestepping_mimetic:calc_vert_velocity_mim_bottomup', &
          & 'MISMATCH IN SURFACE EQUATION:')
        CALL message('mo_ocean_ab_timestepping_mimetic:calc_vert_velocity_mim_bottomup', &
          & 'Elevation change does not match vertical velocity')
        WRITE(message_text,'(2(a,e20.12))') ' (h_new-h_old)/dtime - w = ', MAXVAL(ABS(z_c(:,:))), &
          & ' z_abort=', z_abort
        CALL message ('mo_ocean_ab_timestepping_mimetic:calc_vert_velocity_mim_bottomup', message_text)
        CALL finish(TRIM('mo_ocean_ab_timestepping_mimetic:calc_vert_velocity_mim_bottomup'), &
          & 'MISMATCH in surface equation')
      ENDIF
    ENDIF ! (debug_check_level > 8)
    !---------------------------------------------------------------------
    IF (debug_check_level > 20) THEN
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
        DO jc = start_index, end_index
          IF(patch_3d%lsm_c(jc,1,blockNo) > sea_boundary) THEN
            IF (ocean_state%p_prog(nnew(1))%h(jc,blockNo) /= 0.0_wp) &
              & CALL finish("calc_vert_velocity_mim_bottomup", "h(jc,blockNo) /= 0 on land")
          ENDIF
        END DO
      END DO
    ENDIF
  END SUBROUTINE calc_vert_velocity_mim_bottomup
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !!  The result is NOT synced. Should be done in the calling method if required
  FUNCTION invert_mass_matrix(patch_3d, ocean_state, op_coeffs, rhs_e) result(inv_flip_flop_e)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_operator_coeff),INTENT(in), TARGET :: op_coeffs
    REAL(wp), INTENT(INOUT), TARGET :: rhs_e(:,:,:)!(nproma,n_zlev,patch_2D%nblks_e)
    REAL(wp) :: inv_flip_flop_e(SIZE(rhs_e,1),SIZE(rhs_e,2),SIZE(rhs_e,3))
    INTEGER :: jk, n_it, n_it_sp
    REAL(wp) :: rn
    CHARACTER(LEN=max_char_length) :: string
    TYPE(t_primal_flip_flop_lhs), POINTER :: lhs_pff
    TYPE(t_trivial_transfer), POINTER :: trans_triv
    CLASS(t_lhs_agen), POINTER :: lhs
    CLASS(t_transfer), POINTER :: trans
    TYPE(t_ocean_solve), POINTER :: solve
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_ocean_solve_parm) :: par, par_sp

! allocate + init solver related objects at first call of invert_mass_matrix() (GMRES is chosen)
    IF (.NOT. ASSOCIATED(inv_mm_solver)) THEN
      ALLOCATE(t_primal_flip_flop_lhs :: inv_mm_solver_lhs)
      lhs => lhs_agen_ptr(inv_mm_solver_lhs)
      ALLOCATE(t_trivial_transfer :: inv_mm_solver_trans)
      trans_triv => trivial_transfer_ptr(inv_mm_solver_trans)
      trans => transfer_ptr(inv_mm_solver_trans)
      patch_2d => patch_3d%p_patch_2d(1)
      CALL trans_triv%construct(solve_edge, patch_2d)
      CALL par%init(solve_precon_none, solver_max_restart_iterations, &
        & solver_max_iter_per_restart, patch_2d%cells%in_domain%end_block, &
        & SIZE(rhs_e,3), SIZE(rhs_e,1), patch_2d%edges%in_domain%end_index, &
        & MassMatrix_solver_tolerance, .true.)
      par_sp%nidx = -1
      ALLOCATE(t_ocean_solve :: inv_mm_solver)
      solve => ocean_solve_ptr(inv_mm_solver)
      CALL solve%construct(solve_gmres, par, par_sp, lhs, trans)
    END IF
    lhs_pff => lhs_primal_flip_flop_ptr(inv_mm_solver_lhs)
    solve => ocean_solve_ptr(inv_mm_solver)
    DO jk=1, n_zlev
! re-initialize lhs for each level...
      CALL lhs_pff%construct(patch_3d, op_coeffs, jk)
! fill rhs and guess->solution arrays
!!$OMP PARALLEL WORKSHARE
      solve%x_loc_wp(:,:) = 0._wp
!!$OMP END PARALLEL WORKSHARE
      solve%b_loc_wp => rhs_e(:,jk,:)
! call solver
      CALL solve%solve(n_it, n_it_sp)
      rn = MERGE((solve%res_loc_wp(1)), 0._wp, n_it .GT. 0)
! copy solution into result array
!!$OMP PARALLEL WORKSHARE
      inv_flip_flop_e(:,jk,:) = solve%x_loc_wp(:,:)
!!$OMP END PARALLEL WORKSHARE
      IF (idbg_mxmn >= 1) THEN
        WRITE(string,'(a,i4,a,e28.20)') &
          & 'ocean_restart_gmres iteration =', n_it - 1,', residual =', rn
        CALL message('invert_mass_matrix',TRIM(string))
      ENDIF
    END DO
  END FUNCTION invert_mass_matrix

END MODULE mo_ocean_ab_timestepping_mimetic
