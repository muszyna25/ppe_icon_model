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
!----------------------------
MODULE mo_ocean_ab_timestepping_mimetic

  USE mo_kind,                      ONLY: wp, sp
  USE mo_parallel_config,           ONLY: nproma, l_fast_sum
  USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  USE mo_sync,                      ONLY: sync_e, sync_c, sync_patch_array, sync_patch_array_mult
  USE mo_impl_constants,            ONLY: sea_boundary, &  !  sea,                          &
    & max_char_length, min_dolic
  USE mo_dbg_nml,                   ONLY: idbg_mxmn
  USE mo_ocean_nml,                 ONLY: n_zlev, solver_tolerance,&
    & ab_const, ab_beta, ab_gam, iswm_oce,                &
    & iforc_oce,                                          &
    & no_tracer, l_rigid_lid, l_edge_based,               &
    & use_absolute_solver_tolerance,                      &
    & solver_max_restart_iterations,                      &
    & solver_max_iter_per_restart, dhdtw_abort,           &
    & forcing_enable_freshwater, select_solver,           &
    & select_restart_gmres, select_gmres,                 &
    & use_continuity_correction,                          &
    & select_restart_mixedPrecision_gmres,                &
    & solver_max_iter_per_restart_sp,                     &
    & solver_tolerance_sp, l_partial_cells,               &
    & No_Forcing,                                         &
    & MASS_MATRIX_INVERSION_TYPE,NO_INVERSION,            &
    & MASS_MATRIX_INVERSION_ADVECTION,                    &
    & MASS_MATRIX_INVERSION_ALLTERMS,                     &
    & physics_parameters_type,                            &
    & physics_parameters_ICON_PP_Edge_vnPredict_type,     &
    & solver_FirstGuess, MassMatrix_solver_tolerance
    
  USE mo_run_config,                ONLY: dtime, ltimer, debug_check_level
  USE mo_timer  
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_physical_constants,        ONLY: grav,rho_inv
  USE mo_ocean_initialization,      ONLY: is_initial_timestep
  USE mo_ocean_types,               ONLY: t_hydro_ocean_state, t_hydro_ocean_diag
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,            ONLY: t_external_data
  USE mo_ocean_gmres,               ONLY: ocean_restart_gmres, gmres_oce_old, gmres_oce_e2e, &
    & ocean_restart_gmres_singlePrecesicion,ocean_restart_gmres_e2e
  USE mo_exception,                 ONLY: message, finish, warning, message_text
  USE mo_util_dbg_prnt,             ONLY: dbg_print, debug_print_MaxMinMean
  USE mo_ocean_boundcond,           ONLY: VelocityBottomBoundaryCondition_onBlock, top_bound_cond_horz_veloc
  USE mo_ocean_thermodyn,           ONLY: calculate_density, calc_internal_press, calc_internal_press_grad
  USE mo_ocean_physics_types,       ONLY: t_ho_params
  USE mo_ocean_pp_scheme,           ONLY: update_physics_parameters_ICON_PP_Edge_vnPredict_scheme
  USE mo_sea_ice_types,             ONLY: t_sfc_flx
  USE mo_scalar_product,            ONLY:   &
    & calc_scalar_product_veloc_3d,         &
    & map_edges2edges_viacell_3d_const_z,   &
    & map_edges2edges_viacell_2d_constZ_onTriangles_sp, &
    & map_edges2edges_viacell_2D_per_level,map_scalar_prismtop2center
  USE mo_ocean_math_operators,      ONLY: div_oce_3d, grad_fd_norm_oce_3d,        &
    & grad_fd_norm_oce_2d_3d, grad_fd_norm_oce_2d_3d_sp,                          &
    & grad_fd_norm_oce_2d_onBlock, div_oce_2D_onTriangles_onBlock, &
    & div_oce_3D_onTriangles_onBlock, div_oce_2D_onTriangles_onBlock_sp,          &
    & smooth_onCells, div_oce_2D_general_onBlock, div_oce_2D_general_onBlock_sp,  &
	& div_oce_3D_general_onBlock
  USE mo_ocean_veloc_advection,     ONLY: veloc_adv_horz_mimetic, veloc_adv_vert_mimetic
  
  USE mo_ocean_diffusion,           ONLY: velocity_diffusion,&
    & velocity_diffusion_vertical_implicit_onBlock
  USE mo_ocean_types,               ONLY: t_operator_coeff, t_solverCoeff_singlePrecision
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_grid_config,               ONLY: n_dom
  USE mo_parallel_config,           ONLY: p_test_run
  USE mo_mpi,                       ONLY: my_process_is_stdio, get_my_global_mpi_id, work_mpi_barrier ! my_process_is_mpi_parallel
  USE mo_statistics,                ONLY: global_minmaxmean, print_value_location
  IMPLICIT NONE
  
  PRIVATE  
  !
  PUBLIC :: solve_free_sfc_ab_mimetic
  PUBLIC :: calc_normal_velocity_ab_mimetic
  PUBLIC :: calc_vert_velocity_mim_bottomup
  PUBLIC :: construct_ho_lhs_fields_mimetic, destruct_ho_lhs_fields_mimetic
  PUBLIC :: invert_mass_matrix
  !
  
  INTEGER, PARAMETER :: topLevel=1
  CHARACTER(LEN=12)  :: str_module = 'oceSTEPmimet'  ! Output of module for 1 line debug
  INTEGER :: idt_src    = 1               ! Level of detail for 1 line debug
  
  ! these are allocated once for efficiency and used only by the lhs for the solver
  onCells_2D :: lhs_result
  onEdges_2D :: lhs_z_grad_h, lhs_z_e
  onCells_2D_RealPrecision(sp) :: lhs_result_sp
  onEdges_2D_RealPrecision(sp) :: lhs_z_grad_h_sp, lhs_z_e_sp
  onCells_2D :: z_h_c

  ! the same as above in single precision
!   REAL(wp), ALLOCATABLE, TARGET :: lhs_result(:,:)  ! (nproma,patch%alloc_cell_blocks)
!   REAL(wp), ALLOCATABLE :: lhs_z_grad_h(:,:)
!   REAL(wp), ALLOCATABLE :: lhs_z_e     (:,:)
!   ! the same as above in single precision
!   REAL(sp), ALLOCATABLE, TARGET :: lhs_result_sp(:,:)  ! (nproma,patch%alloc_cell_blocks)
!   REAL(sp), ALLOCATABLE :: lhs_z_grad_h_sp(:,:)
!   REAL(sp), ALLOCATABLE :: lhs_z_e_sp     (:,:)


  ! TYPE(t_cartesian_coordinates), ALLOCATABLE :: lhs_z_grad_h_cc(:,:)
  REAL(wp), PARAMETER ::  min_top_height = 0.05_wp ! we have to have at least 5cm water on topLevel of sea cells
  
CONTAINS
  
  
  !-------------------------------------------------------------------------
  !>
  !! !  Solves the free surface equation.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE solve_free_sfc_ab_mimetic(patch_3d, ocean_state, p_ext_data, p_sfc_flx, &
    & p_phys_param, timestep, op_coeffs, solverCoeff_sp, return_status)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)       :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET             :: ocean_state
    TYPE(t_external_data), TARGET, INTENT(in)     :: p_ext_data
    TYPE(t_sfc_flx), INTENT(inout)                :: p_sfc_flx
    TYPE (t_ho_params)                            :: p_phys_param
    INTEGER, INTENT(in)                           :: timestep
    TYPE(t_operator_coeff)                        :: op_coeffs
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp
    INTEGER :: return_status
    !
    !Local variables
    !
    INTEGER,PARAMETER :: nmax_iter   = 800      ! maximum number of iterations
    INTEGER :: n_iter                          ! actual number of iterations
    INTEGER :: iter_sum                        ! sum of iterations
    INTEGER :: jc,blockNo,je   ! ,jk,il_v1,il_v2,ib_v1,ib_v2
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: start_edge_index, end_edge_index
    ! REAL(wp) :: z_h_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    LOGICAL :: lprecon         = .FALSE.
    ! REAL(wp) :: z_implcoeff
    REAL(wp) :: zresidual(nmax_iter)    ! norms of the residual (convergence history);an argument of dimension at least m is required
    REAL(wp) :: residual_norm
    REAL(wp) :: relative_tolerance, absolute_tolerance               ! (relative or absolute) tolerance
    LOGICAL :: maxIterations_isReached     ! true if reached m iterations
    INTEGER :: gmres_restart_iterations
    CHARACTER(LEN=max_char_length) :: string
    TYPE(t_subset_range), POINTER :: all_cells, all_edges, owned_cells, owned_edges
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: vol_h(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) ::  minmaxmean(3)

    REAL(sp) :: zresidual_sp(nmax_iter)    ! norms of the residual (convergence history);an argument of dimension at least m is required
    REAL(sp) :: residual_norm_sp
    REAL(sp) :: z_h_c_sp(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(sp) :: h_old_sp(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(sp) :: rhs_sfc_eq_sp(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    CHARACTER(len=*), PARAMETER :: method_name='mo_ocean_ab_timestepping_mimetic:solve_free_sfc_ab_mimetic'
    !-------------------------------------------------------------------------------
    patch_2D     => patch_3d%p_patch_2d(1)
    all_cells    => patch_2D%cells%ALL
    all_edges    => patch_2D%edges%ALL
    owned_cells  => patch_2D%cells%owned
    owned_edges  => patch_2D%edges%owned
    return_status = 0
    !-------------------------------------------------------------------------------

      !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('on entry: h-old'                ,ocean_state%p_prog(nold(1))%h ,str_module, idt_src, in_subset=owned_cells)
    CALL dbg_print('on entry: vn-old'               ,ocean_state%p_prog(nold(1))%vn,str_module, idt_src, in_subset=owned_edges)

    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('on entry: h-new'                ,ocean_state%p_prog(nnew(1))%h ,str_module, idt_src, in_subset=owned_cells)
!    IF (MOD(timestep,12) == 0) THEN
!       idbg_mxmn=5
!       CALL dbg_print('on entry: vn-new'             ,ocean_state%p_prog(nnew(1))%vn,str_module, idt_src, in_subset=owned_edges)
!       CALL dbg_print('on entry: w'                  ,ocean_state%p_diag%w,str_module, idt_src,  in_subset=owned_cells)
!       idbg_mxmn=1
!    ELSE
    idt_src=2
       CALL dbg_print('on entry: vn-new'               ,ocean_state%p_prog(nnew(1))%vn,str_module, idt_src, in_subset=owned_edges)
!    ENDIF
    !---------------------------------------------------------------------
    
    ! abort condition for elevation and vn: - better use CFL criterion
    ! IF ( (MAXVAL(ocean_state%p_prog(nnew(1))%h)  >  1.e6_wp) .OR. &
    !   & (MINVAL(ocean_state%p_prog(nnew(1))%h)  < -1.e6_wp) .OR. &
    !   & (MAXVAL(ocean_state%p_prog(nold(1))%vn) >  1.e6_wp) .OR. &
    !   & (MINVAL(ocean_state%p_prog(nnew(1))%vn) < -1.e6_wp) ) THEN
    !   CALL message('Solve free surface AB mimetic: ',' INSTABLE VN or H - stop now ')
    !   CALL finish ('Solve free surface AB mimetic: ',' INSTABLE VN or H !!')
    ! END IF
    
    ! Apply windstress
    CALL top_bound_cond_horz_veloc(patch_3d, ocean_state, op_coeffs, p_sfc_flx) ! ,     &
 !    & ocean_state%p_aux%bc_top_u, ocean_state%p_aux%bc_top_v, &
 !    & ocean_state%p_aux%bc_top_veloc_cc)
    
    ! Apply bot boundary condition for horizontal velocity
    ! This is done when calculating the vn_pred
!     CALL bot_bound_cond_horz_veloc(patch_3d, ocean_state, p_phys_param, op_coeffs)

    ! CALL dbg_print('bc_top_vn', ocean_state%p_aux%bc_top_vn, str_module, idt_src,  in_subset=owned_edges)
    ! CALL dbg_print('bc_bot_vn', ocean_state%p_aux%bc_bot_vn, str_module, idt_src,  in_subset=owned_edges)
    !---------------------------------------------------------------------

    
    start_timer(timer_ab_expl,3)
    CALL calculate_explicit_term_ab(patch_3d, ocean_state, p_phys_param, &
      & is_initial_timestep(timestep), op_coeffs)
    stop_timer(timer_ab_expl,3)
    
    IF(.NOT.l_rigid_lid)THEN
      
      ! Calculate RHS of surface equation
      start_detail_timer(timer_ab_rhs4sfc,5)
      CALL fill_rhs4surface_eq_ab(patch_3d, ocean_state, p_sfc_flx, op_coeffs)
      stop_detail_timer(timer_ab_rhs4sfc,5)

      
      
      ! Solve surface equation with ocean_gmres solver
!       z_h_c = 0.0_wp
      SELECT CASE (solver_FirstGuess)
      CASE (1)
        CALL smooth_onCells(patch_3D=patch_3d,      &
          & in_value=ocean_state%p_prog(nold(1))%h, & 
          & out_value=z_h_c,                        &
          & smooth_weights=(/ 0.5_wp, 0.5_wp /),    &
          & has_missValue=.false., missValue=-999999.0_wp)
        CALL sync_patch_array(sync_c, patch_3d%p_patch_2d(1), z_h_c)

      CASE default
        z_h_c = 0.0_wp
      END SELECT

      CALL dbg_print('bef ocean_gmres: h-old',ocean_state%p_prog(nold(1))%h(:,:) ,str_module,idt_src,in_subset=owned_cells)
!       CALL dbg_print('p_rhs_sfc_eq',ocean_state%p_aux%p_rhs_sfc_eq, str_module,1, in_subset=owned_cells)

      
      SELECT CASE (select_solver)

      !-----------------------------------------------------------------------------------------
      CASE (select_gmres)
        IF(lprecon)THEN
          !ocean_state%p_aux%p_rhs_sfc_eq = ocean_state%p_aux%p_rhs_sfc_eq *patch%cells%area

          CALL gmres_oce_old( z_h_c(:,:),       &  ! arg 1 of lhs. x input is the first guess.
            & lhs_surface_height_ab_mim, &  ! function calculating l.h.s.
            & ocean_state%p_diag%thick_e,       &  ! edge thickness for LHS
            & ocean_state%p_diag%thick_c,       &  ! ocean_state%p_diag%thick_c, &
          ! arg 6 of lhs ocean_state%p_prog(nold(1))%h,
          ! ocean_state%p_diag%cons_thick_c(:,1,:),&
            & ocean_state%p_prog(nold(1))%h,    &  ! arg 2 of lhs !not used
            & patch_3d,                &  ! arg 3 of lhs
            ! & z_implcoeff,               &  ! arg 4 of lhs
            & op_coeffs,                &
            & ocean_state%p_aux%p_rhs_sfc_eq,   &  ! right hand side as input
            & solver_tolerance,                 &  ! relative tolerance
            & use_absolute_solver_tolerance,     &  ! NOT absolute tolerance
            & nmax_iter,                 &  ! max. # of iterations to do
            & maxIterations_isReached,                 &  ! out: .true. = not converged
            & n_iter,                    &  ! out: # of iterations done
            & zresidual,                 &  ! inout: the residual (array)
            & jacobi_precon )

        ELSEIF(.NOT.lprecon)THEN
          CALL gmres_oce_old( z_h_c(:,:),       &  ! arg 1 of lhs. x input is the first guess.
            & lhs_surface_height_ab_mim, &  ! function calculating l.h.s.
            & ocean_state%p_diag%thick_e,       &  ! edge thickness for LHS
            & ocean_state%p_diag%thick_c,       &  ! ocean_state%p_diag%thick_c, &
          ! arg 6 of lhs ocean_state%p_prog(nold(1))%h,
          ! ocean_state%p_diag%cons_thick_c(:,1,:),&
            & ocean_state%p_prog(nold(1))%h,    &  ! arg 2 of lhs !not used
            & patch_3d,                &  ! arg 3 of lhs
            ! & z_implcoeff,               &  ! arg 4 of lhs
            & op_coeffs,                &
            & ocean_state%p_aux%p_rhs_sfc_eq,   &  ! right hand side as input
            & solver_tolerance,                 &  ! relative tolerance
            & use_absolute_solver_tolerance,                 &  ! NOT absolute tolerance
!              & .FALSE.,                   &  ! NOT absolute tolerance
            & nmax_iter,                 &  ! max. # of iterations to do
            & maxIterations_isReached,                 &  ! out: .true. = not converged
            & n_iter,                    &  ! out: # of iterations done
            & zresidual )
        ENDIF

        IF (maxIterations_isReached) THEN
          return_status = 2
          CALL warning('GMRES_oce_old:', "NOT YET CONVERGED !!")
          RETURN
        ENDIF

        ! output print level idt_src used for ocean gmres output with call message:
        IF(n_iter==0)n_iter=1
        idt_src=1
        IF (idbg_mxmn >= idt_src) THEN
          WRITE(string,'(a,i4,a,e28.20)') &
            & 'iteration =', n_iter,', residual =', ABS(zresidual(n_iter))
          CALL message('GMRES_oce_old: surface height',TRIM(string))
        ENDIF

        IF(lprecon)THEN
          ocean_state%p_prog(nnew(1))%h = z_h_c!*patch_2D%cells%area
        ELSEIF(.NOT.lprecon)THEN
          ocean_state%p_prog(nnew(1))%h = z_h_c
        ENDIF

        iter_sum = n_iter*n_iter
        
      !-----------------------------------------------------------------------------------------
      CASE (select_restart_gmres)

        ! call the new gmre_oce, uses out of order global sum
        ! residual_norm = solver_tolerance + 1.0_wp
        gmres_restart_iterations = 0
        iter_sum                 = 0
        maxIterations_isReached  = .true.
        IF (use_absolute_solver_tolerance) THEN          
          absolute_tolerance      = solver_tolerance
        ELSE
          relative_tolerance      = solver_tolerance
          absolute_tolerance      = 0.0_wp
        ENDIF
        ! write(0,*) tolerance, solver_tolerance, residual_norm, gmres_restart_iterations, solver_max_restart_iterations
        DO WHILE(maxIterations_isReached .AND. gmres_restart_iterations < solver_max_restart_iterations)
          
          CALL ocean_restart_gmres( z_h_c(:,:),                   &  ! arg 1 of lhs. x input is the first guess.
            & lhs_surface_height_ab_mim,     &  ! function calculating l.h.s.
            & ocean_state%p_diag%thick_e,           &  ! edge thickness for LHS
            & ocean_state%p_diag%thick_c,           &  ! ocean_state%p_diag%thick_c, &
            & ocean_state%p_prog(nold(1))%h,        &  ! arg 2 of lhs !not used
            & patch_3d,                    &  ! arg 3 of lhs
            ! & z_implcoeff,                   &  ! arg 4 of lhs
            & op_coeffs,                    &
            & ocean_state%p_aux%p_rhs_sfc_eq,       &  ! right hand side as input
            & absolute_tolerance,                &     ! inout, if > 0 then is used as absolute_tolerance, out otherwise
            & relative_tolerance,                &
            & solver_max_iter_per_restart,       &  ! max. # of iterations to do
            & maxIterations_isReached,                     &  ! out: .true. = not converged
            & n_iter,                        &  ! out: # of iterations done
            & zresidual)

            
          IF(n_iter==0) THEN
            residual_norm = 0.0_wp
          ELSE
            residual_norm =  ABS(zresidual(n_iter))
          ENDIF
          
          iter_sum = iter_sum + n_iter
          ! output print level idt_src used for ocean_restart_gmres output with call message:
          idt_src=2
          IF (idbg_mxmn >= idt_src) THEN
            WRITE(string,'(a,i4,a,e28.20)') &
              & 'ocean_restart_gmres iteration =', n_iter,', residual =', residual_norm
            CALL message('GMRES_oce_new: surface height',TRIM(string))
          ENDIF
          
          gmres_restart_iterations = gmres_restart_iterations + 1
          
        END DO ! WHILE(tolerance >= solver_tolerance)
        
        ! output of sum of iterations every timestep
        idt_src=0
        IF (idbg_mxmn >= idt_src) THEN
          WRITE(string,'(a,i4,a,e28.20)') &
            & 'SUM of ocean_restart_gmres iteration =', iter_sum,', residual =', residual_norm
          CALL message('ocean_restart_gmres: surface height',TRIM(string))
        ENDIF
        
        IF (residual_norm > solver_tolerance) THEN
          return_status = 2
          CALL warning(method_name, "NOT YET CONVERGED !!")
          RETURN
        ENDIF
        
        ocean_state%p_prog(nnew(1))%h = z_h_c
        
      !-----------------------------------------------------------------------------------------
      CASE (select_restart_mixedPrecision_gmres)
        ! call the new gmre_oce, uses out of order global sum
        residual_norm_sp = solver_tolerance_sp + 1.0_sp
        gmres_restart_iterations = 0
        iter_sum                 = 0
        z_h_c_sp(:,:)            = 0.0_sp
        h_old_sp(:,:)            = REAL(ocean_state%p_prog(nold(1))%h(:,:), sp)
        rhs_sfc_eq_sp(:,:)       = REAL(ocean_state%p_aux%p_rhs_sfc_eq(:,:), sp)
        ! write(0,*) tolerance, solver_tolerance, residual_norm, gmres_restart_iterations, solver_max_restart_iterations
        DO WHILE(residual_norm_sp >= solver_tolerance_sp .AND. gmres_restart_iterations < solver_max_restart_iterations)

          CALL ocean_restart_gmres_singlePrecesicion( &
            & z_h_c_sp(:,:),                        &  ! arg 1 of lhs. x input is the first guess.
            & lhs_surface_height_ab_mim_sp,         &  ! function calculating l.h.s.
            & h_old_sp,                             &  ! arg 2 of lhs !not used
            & patch_3d,                             &
            & solverCoeff_sp,                       &
            & rhs_sfc_eq_sp,                        &  ! right hand side as input
            & solver_tolerance_sp,                  &  ! tolerance
            & use_absolute_solver_tolerance,        &  ! use absolute tolerance = true
            & solver_max_iter_per_restart_sp,       &  ! max. # of iterations to do
            & maxIterations_isReached,                            &  ! out: .true. = not converged
            & n_iter,                               &  ! out: # of iterations done
            & zresidual_sp )

          IF(n_iter==0) THEN
            residual_norm_sp = 0.0_sp
          ELSE
            residual_norm_sp =  ABS(zresidual_sp(n_iter))
          ENDIF

          iter_sum = iter_sum + n_iter
          ! output print level idt_src used for ocean_restart_gmres output with call message:
          idt_src=2
          IF (idbg_mxmn >= idt_src) THEN
            WRITE(string,'(a,i4,a,e28.20)') &
              & 'ocean_restart_gmres_sp iteration =', n_iter,', residual =', residual_norm_sp
            CALL message('GMRES_oce_new: surface height',TRIM(string))
          ENDIF

          gmres_restart_iterations = gmres_restart_iterations + 1

        END DO ! WHILE(tolerance >= solver_tolerance)

        ! output of sum of iterations every timestep
        idt_src=0
        IF (idbg_mxmn >= idt_src) THEN
          WRITE(string,'(a,i4,a,e28.20)') &
            & 'SUM of ocean_restart_gmres_sp iteration =', iter_sum,', residual =', residual_norm_sp
          CALL message('ocean_restart_gmres: surface height',TRIM(string))
        ENDIF


        ! call the double precision new gmre_oce, uses out of order global sum
        residual_norm = solver_tolerance + 1.0_wp
        gmres_restart_iterations = 0
        iter_sum                 = 0
        z_h_c(:,:)               = REAL(z_h_c_sp(:,:), wp)
        absolute_tolerance      = solver_tolerance
        ! write(0,*) tolerance, solver_tolerance, residual_norm, gmres_restart_iterations, solver_max_restart_iterations
        DO WHILE(residual_norm >= solver_tolerance .AND. gmres_restart_iterations < solver_max_restart_iterations)

          CALL ocean_restart_gmres( z_h_c(:,:),                   &  ! arg 1 of lhs. x input is the first guess.
            & lhs_surface_height_ab_mim,     &  ! function calculating l.h.s.
            & ocean_state%p_diag%thick_e,           &  ! edge thickness for LHS
            & ocean_state%p_diag%thick_c,           &  ! ocean_state%p_diag%thick_c, &
            & ocean_state%p_prog(nold(1))%h,        &  ! arg 2 of lhs !not used
            & patch_3d,                    &  ! arg 3 of lhs
            ! & z_implcoeff,                   &  ! arg 4 of lhs
            & op_coeffs,                    &
            & ocean_state%p_aux%p_rhs_sfc_eq,       &  ! right hand side as input
            & absolute_tolerance,                &     ! inout, if > 0 then is used as absolute_tolerance, out otherwise
            & relative_tolerance,                &
            & solver_max_iter_per_restart,   &  ! max. # of iterations to do
            & maxIterations_isReached,                     &  ! out: .true. = not converged
            & n_iter,                        &  ! out: # of iterations done
            & zresidual )

          IF(n_iter==0) THEN
            residual_norm = 0.0_wp
          ELSE
            residual_norm =  ABS(zresidual(n_iter))
          ENDIF

          iter_sum = iter_sum + n_iter
          ! output print level idt_src used for ocean_restart_gmres output with call message:
          idt_src=2
          IF (idbg_mxmn >= idt_src) THEN
            WRITE(string,'(a,i4,a,e28.20)') &
              & 'ocean_restart_gmres iteration =', n_iter,', residual =', residual_norm
            CALL message('GMRES_oce_new: surface height',TRIM(string))
          ENDIF

          gmres_restart_iterations = gmres_restart_iterations + 1

        END DO ! WHILE(tolerance >= solver_tolerance)

        ! output of sum of iterations every timestep
        idt_src=0
        IF (idbg_mxmn >= idt_src) THEN
          WRITE(string,'(a,i4,a,e28.20)') &
            & 'SUM of ocean_restart_gmres iteration =', iter_sum,', residual =', residual_norm
          CALL message('mixed ocean_restart_gmres: surface height',TRIM(string))
        ENDIF

        IF (residual_norm > solver_tolerance) THEN
          return_status = 2
          CALL warning(method_name, "NOT YET CONVERGED !!")
          RETURN
        ENDIF
        
        ocean_state%p_prog(nnew(1))%h = z_h_c

      !-----------------------------------------------------------------------------------------
      CASE default
        CALL finish(method_name, "Unknown solver")

      END SELECT ! solver

      !-------- end of solver ---------------
      CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_prog(nnew(1))%h)
      
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      ! idt_src=2  ! output print level (1-5, fix)
      !     z_h_c = lhs_surface_height_ab_mim( ocean_state%p_prog(nnew(1))%h, &
      !         & ocean_state%p_prog(nold(1))%h, &
      !         & patch_3d,             &
      !         & z_implcoeff,            &
      !         & ocean_state%p_diag%thick_e,    &
      !         & ocean_state%p_diag%thick_c,    &
      !         & op_coeffs)             &
      !         & -ocean_state%p_aux%p_rhs_sfc_eq
      !     CALL dbg_print('SolvSfc: residual h-res'    ,z_h_c                  ,str_module,idt_src)
!      vol_h(:,:) = patch_3d%p_patch_2d(n_dom)%cells%area(:,:) * ocean_state%p_prog(nnew(1))%h(:,:)
!      CALL dbg_print('after ocean_gmres: vol_h(:,:)',vol_h ,str_module,idt_src, in_subset=owned_cells)
      !---------------------------------------------------------------------
      minmaxmean(:) = global_minmaxmean(values=ocean_state%p_prog(nnew(1))%h(:,:), in_subset=owned_cells)
      idt_src=1  ! output print level (1-5, fix)
      CALL debug_print_MaxMinMean('after ocean_gmres: h-new', minmaxmean, str_module, idt_src)
      IF (minmaxmean(1) + patch_3D%p_patch_1D(1)%del_zlev_m(1) <= min_top_height) THEN
!          CALL finish(method_name, "height below min_top_height")
        CALL warning(method_name, "height below min_top_height")
        CALL print_value_location(ocean_state%p_prog(nnew(1))%h(:,:), minmaxmean(1), owned_cells)
        CALL print_value_location(ocean_state%p_prog(nnew(2))%h(:,:), minmaxmean(1), owned_cells)
        CALL work_mpi_barrier()
        return_status = 1
        RETURN
      ENDIF
      !---------------------------------------------------------------------
      
    ENDIF  ! l_rigid_lid
    
    ! write(0,*) "solve_free_sfc_ab_mimetic: sum(h)=", SUM(ocean_state%p_prog(nnew(1))%h(:,:))
    
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
    & is_first_timestep, op_coeffs)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE (t_ho_params)                   :: p_phys_param
    LOGICAL,INTENT(in)                   :: is_first_timestep
    TYPE(t_operator_coeff)               :: op_coeffs
    !
    TYPE(t_subset_range), POINTER :: owned_edges, owned_cells
    TYPE(t_patch), POINTER :: patch_2D
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !  &       routine = ('mo_ocean_ab_timestepping_mimetic:calculate_explicit_term_ab')
    !-----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    
    patch_2D        => patch_3d%p_patch_2d(n_dom)
    owned_edges     => patch_3d%p_patch_2d(n_dom)%edges%owned
    owned_cells     => patch_3d%p_patch_2d(n_dom)%cells%owned
    
    !---------------------------------------------------------------------
    ! STEP 1: horizontal advection
    !---------------------------------------------------------------------
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
    
    !---------------------------------------------------------------------
    ! STEP 2: compute 3D contributions: gradient of hydrostatic pressure and vertical velocity advection
    !---------------------------------------------------------------------
    
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
         &                          op_coeffs%grad_coeff,  &
         &                          ocean_state%p_diag%press_grad)     
      
      
      ! this is not needed
      ! CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_diag%press_grad)

      ! calculate vertical velocity advection
      CALL veloc_adv_vert_mimetic(          &
        & patch_3d,                         &
        & ocean_state%p_diag,op_coeffs,     &
        & ocean_state%p_diag%veloc_adv_vert )
          
      stop_detail_timer(timer_extra2,4)
      
      ! calculate vertical velocity diffusion
      !   For the alternative choice "expl_vertical_velocity_diff==1" see couples of
      !   lines below below
!       IF (expl_vertical_velocity_diff==0) THEN
!         CALL velocity_diffusion_vert_explicit( patch_3d,     &
!           & ocean_state%p_diag,            &
!           & ocean_state%p_aux,op_coeffs,  &
!           & p_phys_param,           &
!           & ocean_state%p_diag%laplacian_vert)
!       ENDIF
      
    ELSE  !  iswm_oce=1
      ocean_state%p_diag%veloc_adv_vert = 0.0_wp
      ocean_state%p_diag%laplacian_vert = 0.0_wp
    ENDIF
    
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src = 3  ! output print level (1-5, fix)
    CALL dbg_print('horizontal advection'      ,ocean_state%p_diag%veloc_adv_horz,str_module,idt_src, &
      & in_subset=owned_edges)
    CALL dbg_print('density'                   ,ocean_state%p_diag%rho           ,str_module,idt_src, &
      & in_subset = owned_cells)
    CALL dbg_print('internal pressure'         ,ocean_state%p_diag%press_hyd     ,str_module,idt_src, &
      in_subset = owned_cells)
    CALL dbg_print('internal press grad'       ,ocean_state%p_diag%press_grad    ,str_module,idt_src, &
      in_subset = owned_edges)
    idt_src = 4  ! output print level (1-5, fix)
    CALL dbg_print('kinetic energy'            ,ocean_state%p_diag%kin           ,str_module,idt_src, &
      in_subset = owned_cells)
    CALL dbg_print('vertical advection'        ,ocean_state%p_diag%veloc_adv_vert,str_module,idt_src, &
      & in_subset=owned_edges)
    !---------------------------------------------------------------------
    
    !---------------------------------------------------------------------
    ! STEP 3: compute harmonic or biharmoic laplacian diffusion of velocity.
    !         This term is discretized explicitly. Order and form of the laplacian
    !         are determined in mo_oce_diffusion according to namelist settings
    !---------------------------------------------------------------------
    
    start_detail_timer(timer_extra3,5)
    CALL velocity_diffusion(patch_3d,              &
      & ocean_state%p_prog(nold(1))%vn, &
      & p_phys_param,            &
      & ocean_state%p_diag,op_coeffs,  &
      & ocean_state%p_diag%laplacian_horz)
    stop_detail_timer(timer_extra3,5)
 
    ! CALL dbg_print('bc_top_vn'   ,ocean_state%p_aux%bc_top_vn       ,str_module,idt_src, in_subset=owned_edges)
    ! write(0,*) "MASS_MATRIX_INVERSION_TYPE=", MASS_MATRIX_INVERSION_TYPE
    ! CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_diag%laplacian_horz)
    !---------------------------------------------------------------------
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
    CALL dbg_print('VelocDiff: LaPlacHorz'    ,ocean_state%p_diag%laplacian_horz  ,str_module,idt_src, in_subset=owned_edges)
    IF (iswm_oce /= 1) THEN
      CALL dbg_print('ImplVelocDiff vertical'   ,ocean_state%p_diag%vn_pred       ,str_module,idt_src, in_subset=owned_edges)
    ELSE
      CALL dbg_print('VelocDiff: LaPlacVert'    ,ocean_state%p_diag%laplacian_vert,str_module,idt_src, in_subset=owned_edges)
    ENDIF

    idt_src=5  ! output print level (1-5, fix)
    CALL dbg_print('vn(nold)'                 ,ocean_state%p_prog(nold(1))%vn ,str_module,idt_src, in_subset=owned_edges)
    CALL dbg_print('G_n+1/2 - g_nimd'         ,ocean_state%p_aux%g_nimd       ,str_module,idt_src, in_subset=owned_edges)
    CALL dbg_print('G_n'                      ,ocean_state%p_aux%g_n          ,str_module,idt_src, in_subset=owned_edges)
    CALL dbg_print('G_n-1'                    ,ocean_state%p_aux%g_nm1        ,str_module,idt_src, in_subset=owned_edges)

    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('vn_pred'                   ,ocean_state%p_diag%vn_pred           ,str_module,idt_src, in_subset=owned_edges)
    !---------------------------------------------------------------------    
  END SUBROUTINE calculate_explicit_term_ab
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE explicit_vn_pred( patch_3d, ocean_state, op_coeffs, p_phys_param, is_first_timestep)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_operator_coeff)               :: op_coeffs
    TYPE (t_ho_params)                   :: p_phys_param
    LOGICAL,INTENT(in)                   :: is_first_timestep

    REAL(wp) :: z_gradh_e(nproma)
    INTEGER :: je, jk, blockNo
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_edges
    INTEGER :: start_edge_index, end_edge_index, dolic_e
    TYPE(t_patch), POINTER :: patch_2D
    !REAL(wp) :: z_e(nproma,n_zlev,patch_3d%p_patch_2d(n_dom)%nblks_e)
    !-----------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(n_dom)
    edges_in_domain => patch_3d%p_patch_2d(n_dom)%edges%in_domain
    
!     z_e=0.0_wp
!     CALL map_edges2edges_viacell_2D_per_level( patch_3d,                &
!                                             & ocean_state%p_diag%grad(:,1,:),&  
!                                             & op_coeffs,              &
!                                             &  z_e(:,1,:),1 )
! write(0,*)'before after',maxval(ocean_state%p_diag%grad),minval(ocean_state%p_diag%grad),maxval(z_e),minval(z_e)
!     ocean_state%p_diag%grad=z_e
    
    !---------------------------------------------------------------------
    ! STEP 4: calculate weighted gradient of surface height at previous timestep
    !---------------------------------------------------------------------
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, z_gradh_e) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      z_gradh_e(:)  = 0.0_wp
      CALL grad_fd_norm_oce_2d_onBlock(         &
        & ocean_state%p_prog(nold(1))%h,        &
        & patch_2D,                             &
        & op_coeffs%grad_coeff(:,1,blockNo),    &
        & z_gradh_e(:),                 &
        & start_edge_index, end_edge_index, blockNo)

      z_gradh_e(start_edge_index:end_edge_index) = &
         & (1.0_wp-ab_beta) * grav * z_gradh_e(start_edge_index:end_edge_index)

         
      !---------------------------------------------------------------------
      ! STEP 5:
      !---------------------------------------------------------------------      
      
      CALL calculate_explicit_term_g_n_onBlock( patch_3d, ocean_state, is_first_timestep, &
        & start_edge_index, end_edge_index, blockNo)
        
      IF ( iswm_oce /= 1) THEN
        CALL calculate_explicit_vn_pred_3D_onBlock( patch_3d, ocean_state, z_gradh_e(:),    &
        & start_edge_index, end_edge_index, blockNo)

        ! calculate vertical friction, ie p_phys_param%a_veloc_v
        IF (physics_parameters_type == physics_parameters_ICON_PP_Edge_vnPredict_type) &
          CALL update_physics_parameters_ICON_PP_Edge_vnPredict_scheme(patch_3d, &
               & blockNo, start_edge_index, end_edge_index, ocean_state, ocean_state%p_diag%vn_pred(:,:,blockNo))

        !In 3D case implicit vertical velocity diffusion is chosen
        CALL velocity_diffusion_vertical_implicit_onBlock( &
          & patch_3d,                                      &
          & ocean_state%p_diag%vn_pred(:,:,blockNo),       &
          & p_phys_param%a_veloc_v(:,:,blockNo),           &
          & op_coeffs,                                     &
          & start_edge_index, end_edge_index, blockNo)
        
      ELSE !( iswm_oce == 1) THEN
        CALL calculate_explicit_vn_pred_2D_onBlock( patch_3d, ocean_state, z_gradh_e(:), &
        & start_edge_index, end_edge_index, blockNo) 
      ENDIF
      
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    
  END SUBROUTINE explicit_vn_pred
  !-------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE   explicit_vn_pred_invert_mass_matrix( patch_3d, ocean_state, op_coeffs, p_phys_param, is_first_timestep)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_operator_coeff)               :: op_coeffs
    TYPE (t_ho_params)                   :: p_phys_param
    LOGICAL,INTENT(in)                   :: is_first_timestep

    REAL(wp) :: z_gradh_e(nproma)
    INTEGER :: je, jk, blockNo
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_edges
    INTEGER :: start_edge_index, end_edge_index, dolic_e
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: z_e(nproma,n_zlev,patch_3d%p_patch_2d(n_dom)%nblks_e)
    !-----------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(n_dom)
    edges_in_domain => patch_3d%p_patch_2d(n_dom)%edges%in_domain

    IF(MASS_MATRIX_INVERSION_TYPE== MASS_MATRIX_INVERSION_ADVECTION)THEN
       !Here the inversion of the mass matrix is already carried out        
       
       z_e=ocean_state%p_diag%veloc_adv_horz + ocean_state%p_diag%veloc_adv_vert
        Write(0,*)'ADV before:',&
        &maxval(z_e(:,1,:)),minval(z_e(:,1,:))
            
       ocean_state%p_diag%veloc_adv_horz = invert_mass_matrix(patch_3d, ocean_state, op_coeffs, z_e)

       CALL sync_patch_array(sync_e, patch_2D, z_e)
       
        Write(0,*)'ADV after:',&
        &maxval(ocean_state%p_diag%veloc_adv_horz(:,1,:)),& 
        &minval(ocean_state%p_diag%veloc_adv_horz(:,1,:))
     
    ENDIF
    
    !---------------------------------------------------------------------
    ! STEP 4: calculate weighted gradient of surface height at previous timestep
    !---------------------------------------------------------------------   
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, z_gradh_e) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      z_gradh_e(:)  = 0.0_wp
      CALL grad_fd_norm_oce_2d_onBlock(         &
        & ocean_state%p_prog(nold(1))%h,        &
        & patch_2D,                             &
        & op_coeffs%grad_coeff(:,1,blockNo),    &
        & z_gradh_e(:),                 &
        & start_edge_index, end_edge_index, blockNo)

      z_gradh_e(start_edge_index:end_edge_index) = &
         & (1.0_wp-ab_beta) * grav * z_gradh_e(start_edge_index:end_edge_index)
         
      !---------------------------------------------------------------------
      ! STEP 5:
      !---------------------------------------------------------------------      
      CALL calculate_explicit_term_g_n_onBlock( patch_3d, ocean_state, is_first_timestep, &
        & start_edge_index, end_edge_index, blockNo)
        
      IF ( iswm_oce /= 1) THEN
        CALL calculate_explicit_vn_pred_3D_onBlock( patch_3d, ocean_state, z_gradh_e(:),    &
        & start_edge_index, end_edge_index, blockNo)

        !In 3D case implicit vertical velocity diffusion is chosen
        CALL velocity_diffusion_vertical_implicit_onBlock( &
          & patch_3d,                                      &
          & ocean_state%p_diag%vn_pred(:,:,blockNo),       &
          & p_phys_param%a_veloc_v(:,:,blockNo),           &
          & op_coeffs,                                     &
          & start_edge_index, end_edge_index, blockNo)
        
      ELSE !( iswm_oce == 1)
        CALL calculate_explicit_vn_pred_2D_onBlock( patch_3d, ocean_state, z_gradh_e(:), &
        & start_edge_index, end_edge_index, blockNo) 
        
      ENDIF      
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF(MASS_MATRIX_INVERSION_TYPE == MASS_MATRIX_INVERSION_ALLTERMS )THEN
      !Here the inversion is just prepared
      Write(0,*)'vn_pred before:',&
      &maxval(ocean_state%p_diag%vn_pred(:,1,:)),& 
      &minval(ocean_state%p_diag%vn_pred(:,1,:))
     
!       CALL map_edges2edges_viacell_3d_const_z( patch_3d,                 &
!                                             & ocean_state%p_diag%vn_pred(:,1,:),&  
!                                             & op_coeffs,                 &
!                                             & ocean_state%p_diag%vn_pred_ptp(:,1,:) )
      CALL map_edges2edges_viacell_2D_per_level( patch_3d,                 &
                                            & ocean_state%p_diag%vn_pred(:,1,:),&  
                                            & op_coeffs,                 &
                                            & ocean_state%p_diag%vn_pred_ptp(:,1,:),1 )
                                            
      Write(0,*)'vn_pred after:',&
       &maxval(ocean_state%p_diag%vn_pred_ptp(:,1,:)),& 
       &minval(ocean_state%p_diag%vn_pred_ptp(:,1,:))
       
      ocean_state%p_diag%vn_pred=ocean_state%p_diag%vn_pred_ptp       
    END IF
  
    
  END SUBROUTINE   explicit_vn_pred_invert_mass_matrix
  !-------------------------------------------------------------------------
  

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

 
  !-------------------------------------------------------------------------
  !>
  !!  Calculation of right-hand side of elliptic surface equation.
  !!  This is used in semi implicit timelevel stepping.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE fill_rhs4surface_eq_ab( patch_3d, ocean_state, p_sfc_flx, op_coeffs)
    !
    ! Patch on which computation is performed
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_sfc_flx), INTENT(in)          :: p_sfc_flx
    TYPE(t_operator_coeff)               :: op_coeffs
    !
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
    !REAL(wp) :: thick
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !       & routine = ('mo_ocean_ab_timestepping_mimetic:fill_rhs4surface_eq_ab')
    !-------------------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_3d%p_patch_2d(1)%cells%in_domain
    edges_in_domain => patch_3d%p_patch_2d(1)%edges%in_domain
    owned_edges     => patch_3d%p_patch_2d(1)%edges%owned
    
    inv_gdt2 = 1.0_wp / (grav*dtime*dtime)
    
    ! LL: this should not be required
    ! CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_diag%vn_pred)
    ! CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_prog(nold(1))%vn)
    ! CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_diag%vn_impl_vert_diff)

!     IF (p_test_run) z_vn_ab(:,:,:)  = 0.0_wp

#ifdef NAGFOR
    z_vn_ab(:,:,:)  = 0.0_wp
#endif

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
        
!       ELSEIF(expl_vertical_velocity_diff==0)THEN
!         
!         DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!           CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
!           DO je = start_edge_index, end_edge_index
!             i_dolic_e =  patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
!             DO jk=1,i_dolic_e
!               z_vn_ab(je,jk,blockNo)=ab_gam*ocean_state%p_diag%vn_pred(je,jk,blockNo)&
!                 & + (1.0_wp -ab_gam)* ocean_state%p_prog(nold(1))%vn(je,jk,blockNo)
!             END DO
!           ENDDO
!         END DO
!         
!       ENDIF

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
    
    ! note: these two operators can be combined
    ! div_z_c(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)  = 0.0_wp
    
!	z_e=10.0_wp
!	!z_e(4,1,4)=10.0_wp

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
! write(123,*)'div',jc, div_z_c(jc, 1:1)		   		  
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
    ! write(456,*)'div',jc, div_z_c(jc, 1:1)
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
!     CALL dbg_print('RHS div_z_depth_int_c'     ,div_z_depth_int_c        ,str_module,idt_src, &
!       in_subset=patch_3d%p_patch_2d(1)%cells%owned)
!    CALL dbg_print('RHS div_z_c'     ,div_z_c        ,str_module,idt_src, &
!      in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('RHS final'                 ,ocean_state%p_aux%p_rhs_sfc_eq  ,str_module,idt_src, &
      in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    !---------------------------------------------------------------------
    
  END SUBROUTINE fill_rhs4surface_eq_ab
  !-------------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE construct_ho_lhs_fields_mimetic(patch_3d)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    
    TYPE(t_patch), POINTER :: patch    
    INTEGER :: return_status
    
    patch         => patch_3d%p_patch_2d(1)
    
    ALLOCATE( &
      & lhs_result(nproma,patch%alloc_cell_blocks), &
      & lhs_z_grad_h(nproma,patch%nblks_e),     &
      & lhs_z_e     (nproma,patch%nblks_e),     &
      & z_h_c(nproma,patch%alloc_cell_blocks), &
!      & lhs_z_e_top (nproma,patch%nblks_e),     &
!      & lhs_z_grad_h_cc(nproma,patch%alloc_cell_blocks),  &
      & stat = return_status)
    IF (return_status > 0) &
      & CALL finish("mo_ocean_ab_timestepping_mimetic:init_ho_lhs_fields", "Allocation failed")

    ALLOCATE(lhs_result_sp(nproma,patch%alloc_cell_blocks), &
      & lhs_z_grad_h_sp(nproma,patch%nblks_e),     &
      & lhs_z_e_sp     (nproma,patch%nblks_e),     &
      & stat = return_status)
    IF (return_status > 0) &
      & CALL finish("mo_ocean_ab_timestepping_mimetic:init_ho_lhs_fields", "sp Allocation failed")
    
    ! these are arrays used by the lhs routine
    lhs_result(:,:)   = 0.0_wp
    lhs_z_grad_h(:,:) = 0.0_wp
    lhs_z_e     (:,:) = 0.0_wp
    z_h_c       (:,:) = 0.0_wp
!    lhs_z_e_top (:,:) = 0.0_wp

    lhs_result_sp(:,:)   = 0.0_wp
    lhs_z_grad_h_sp(:,:) = 0.0_wp
    lhs_z_e_sp     (:,:) = 0.0_wp
    

!    lhs_z_grad_h_cc(:,:)%x(1) = 0.0_wp
!    lhs_z_grad_h_cc(:,:)%x(2) = 0.0_wp
!    lhs_z_grad_h_cc(:,:)%x(3) = 0.0_wp

  END SUBROUTINE construct_ho_lhs_fields_mimetic
  !-------------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE destruct_ho_lhs_fields_mimetic()

    DEALLOCATE(         &
      & lhs_result,     &
      & lhs_z_grad_h,   &
      & lhs_z_e     ,   &
      & z_h_c)

    DEALLOCATE(lhs_result_sp, &
      & lhs_z_grad_h_sp,     &
      & lhs_z_e_sp)

    NULLIFY (lhs_result)
    NULLIFY (lhs_z_grad_h)
    NULLIFY (lhs_z_e)
    NULLIFY (z_h_c)
    NULLIFY (lhs_result_sp)
    NULLIFY (lhs_z_grad_h_sp)
    NULLIFY (lhs_z_e_sp)

  END SUBROUTINE destruct_ho_lhs_fields_mimetic
  !-------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------
  !>
  !! Computation of left-hand side of the surface height equation
  !! This function calculates the left-hand side of the surface height equation in
  !! the Adams-Bashforth 2 timestepping. Function is called from iterative solver.
  !! Iteration of height as input
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  FUNCTION lhs_surface_height_ab_mim( x, h_old, patch_3d, thickness_e,&
    & thickness_c,op_coeffs) result(lhs)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp),    INTENT(inout)           :: x(:,:)    ! inout for sync, dimension: (nproma,patch%alloc_cell_blocks)
    REAL(wp),    INTENT(in)              :: h_old(:,:)
    ! REAL(wp),    INTENT(in)              :: coeff
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs
    REAL(wp),    INTENT(in)              :: thickness_e(:,:)
    REAL(wp),    INTENT(in)              :: thickness_c(:,:) !thickness of fluid column
    !  these are small (2D) arrays and allocated once for efficiency
    ! Left-hand side calculated from iterated height
    !
    REAL(wp) :: lhs(SIZE(x,1), SIZE(x,2))  ! (nproma,patch_2D%alloc_cell_blocks)
    
    ! local variables,
    REAL(wp) :: gdt2_inv, gam_times_beta
    INTEGER :: start_index, end_index
    INTEGER :: jc, blockNo, je
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D     ! patch_2D on which computation is performed
    !-----------------------------------------------------------------------
    start_detail_timer(timer_lhs,3)
    !-----------------------------------------------------------------------
    patch_2D           => patch_3d%p_patch_2d(1)
    cells_in_domain    => patch_2D%cells%in_domain
    edges_in_domain    => patch_2D%edges%in_domain
    
    
    gdt2_inv       = 1.0_wp / (grav*(dtime)**2)
    gam_times_beta = ab_gam * ab_beta
    
    lhs   (1:nproma,cells_in_domain%end_block:patch_2D%alloc_cell_blocks)  = 0.0_wp
    
!     CALL dbg_print('h_old', h_old, "lhs_surface_height_ab_mim", 1, &
!         & in_subset=patch_3d%p_patch_2d(1)%cells%owned)
!     CALL dbg_print('thickness_c', thickness_c, "lhs_surface_height_ab_mim", 1, &
!         & in_subset=patch_3d%p_patch_2d(1)%cells%owned)
!     CALL dbg_print('thickness_e', thickness_e, "lhs_surface_height_ab_mim", 1, &
!         & in_subset=patch_3d%p_patch_2d(1)%edges%owned)

    CALL sync_patch_array(sync_c, patch_2D, x(1:nproma,1:patch_2D%cells%all%end_block) )
    
    !---------------------------------------
    start_detail_timer(timer_extra31,6)
    IF(l_edge_based)THEN
    
      !Step 1) Calculate gradient of iterated height.
      CALL grad_fd_norm_oce_2d_3d( x, &
        & patch_2D,                       &
        & op_coeffs%grad_coeff(:,1,:),&
        & lhs_z_grad_h(:,:))
      
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
        DO je = start_index, end_index
          lhs_z_e(je,blockNo) = lhs_z_grad_h(je,blockNo) * thickness_e(je,blockNo)
          
!          IF(patch_3d%lsm_e(je,1,blockNo) > sea_boundary)THEN
!            IF (lhs_z_e(je,blockNo) /= 0.0_wp) &
!              CALL finish("lhs_surface_height_ab_mim", "lhs_z_e(je,blockNo) /= 0 on land")
!          ENDIF

        END DO
      END DO
      
    ELSE  !IF(.NOT.l_edge_based)THEN
      
      CALL grad_fd_norm_oce_2d_3d( x,     &
        & patch_2D,                       &
        & op_coeffs%grad_coeff(:,1,:),    &
        & lhs_z_grad_h(:,:),              &
        & subset_range=patch_2D%edges%gradIsCalculable)


      CALL map_edges2edges_viacell_3d_const_z( patch_3d, lhs_z_grad_h(:,:), op_coeffs, lhs_z_e(:,:))

      
    ENDIF ! l_edge_based
    !---------------------------------------
    stop_detail_timer(timer_extra31,6)
    
    start_detail_timer(timer_extra32,6)
    !Step 3) Calculate divergence
    ! store the div in lhs for reducing memory and improving performance

    IF( patch_2d%cells%max_connectivity == 3 )THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      
        CALL div_oce_2D_onTriangles_onBlock(lhs_z_e, patch_2D, op_coeffs%div_coeff, lhs(:,blockNo), &
          & level=topLevel, blockNo=blockNo, start_index=start_index, end_index=end_index)
    
        !Step 4) Finalize LHS calculations
        DO jc = start_index, end_index        
          !lhs(jc,blockNo) =(x(jc,blockNo) - gdt2 * ab_gam * ab_beta * lhs_div_z_c(jc,blockNo)) / gdt2 !rho_sfc(jc,blockNo)*rho_inv
          lhs(jc,blockNo) = x(jc,blockNo) * gdt2_inv - gam_times_beta * lhs(jc,blockNo)        
        END DO
      END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO
    ELSE

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
        CALL div_oce_2D_general_onBlock( lhs_z_e, patch_2D, op_coeffs%div_coeff, lhs(:,blockNo),&
                & level=topLevel,blockNo=blockNo, start_index=start_index, end_index=end_index)        
        !Step 4) Finalize LHS calculations
        DO jc = start_index, end_index
          !lhs(jc,blockNo) =(x(jc,blockNo) - gdt2 * ab_gam * ab_beta * lhs_div_z_c(jc,blockNo)) / gdt2 !rho_sfc(jc,blockNo)*rho_inv
          lhs(jc,blockNo) = x(jc,blockNo) * gdt2_inv - gam_times_beta * lhs(jc,blockNo)
        END DO
      END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO

    ENDIF

!     CALL dbg_print('x', x, "lhs_surface_height_ab_mim", 1, &
!         & in_subset=patch_3d%p_patch_2d(1)%cells%owned)
!     CALL dbg_print('lhs', lhs, "lhs_surface_height_ab_mim", 1, &
!         & in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    !---------------------------------------
    stop_detail_timer(timer_extra32,6)

    IF (debug_check_level > 20) THEN
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
        DO jc = start_index, end_index
          IF(patch_3d%lsm_c(jc,1,blockNo) > sea_boundary) THEN
            IF (lhs(jc,blockNo) /= 0.0_wp) &
              & CALL finish("lhs_surface_height_ab_mim", "lhs(jc,blockNo) /= 0 on land")
          ENDIF
        END DO
      END DO
    ENDIF
   
    stop_detail_timer(timer_lhs,3)
    
  END FUNCTION lhs_surface_height_ab_mim
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------
  ! as in lhs_surface_height_ab_mim in single precision
!<Optimize:inUse>
  FUNCTION lhs_surface_height_ab_mim_sp( x, h_old, patch_3d, solverCoeffs) result(lhs)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(sp),    INTENT(inout)           :: x(:,:)    ! inout for sync, dimension: (nproma,patch%alloc_cell_blocks)
    REAL(sp),    INTENT(in)              :: h_old(:,:)
    TYPE(t_solverCoeff_singlePrecision),INTENT(in)    :: solverCoeffs
    ! Left-hand side calculated from iterated height
    !
    REAL(sp) :: lhs(SIZE(x,1), SIZE(x,2))  ! (nproma,patch_2D%alloc_cell_blocks)

    ! local variables,
    REAL(sp) :: gdt2_inv, gam_times_beta
    INTEGER :: start_index, end_index
    INTEGER :: jc, blockNo, je
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D     ! patch_2D on which computation is performed
    REAL(wp) :: x_sync(SIZE(x,1), SIZE(x,2))    ! used to syn x, since we cannot synd single precision at the moment
    !-----------------------------------------------------------------------
    IF( patch_2d%cells%max_connectivity /= 3 )THEN
      CALL finish("lhs_surface_height_ab_mim_sp", "only works on triangles")
    ENDIF
    start_timer(timer_lhs_sp,2)
    !-----------------------------------------------------------------------
    patch_2D           => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain

    gdt2_inv       = REAL(1.0_wp / (grav*(dtime)**2),sp)
    gam_times_beta = REAL(ab_gam * ab_beta, sp)

    lhs  (1:nproma,cells_in_domain%end_block:patch_2D%alloc_cell_blocks)  = 0.0_sp

    x_sync(:,:) = REAL(x(:,:), wp)
    CALL sync_patch_array(sync_c, patch_2D, x_sync(1:nproma,1:patch_2D%cells%all%end_block) )
    x(:,:) = REAL(x_sync(:,:), sp)

    !---------------------------------------
    IF(l_edge_based)THEN

      !Step 1) Calculate gradient of iterated height.
      CALL grad_fd_norm_oce_2d_3d_sp( x,   &
        & patch_2D,                     &
        & solverCoeffs%grad_coeff(:,:), &
        & lhs_z_grad_h_sp(:,:))

      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
        DO je = start_index, end_index
          lhs_z_e_sp(je,blockNo) = lhs_z_grad_h_sp(je,blockNo) * solverCoeffs%edge_thickness(je,blockNo)
        END DO
      END DO

    ELSE  !IF(.NOT.l_edge_based)THEN

      ! the map_edges2edges_viacell_3d_const_z should be changed to calculate only
      ! on in_domain_edges. Still, edge values need to be synced
      !Step 1) Calculate gradient of iterated height.
      CALL grad_fd_norm_oce_2d_3d_sp( x,  &
        & patch_2D,                       &
        & solverCoeffs%grad_coeff(:,:),   &
        & lhs_z_grad_h_sp(:,:),           &
        & subset_range=patch_2D%edges%gradIsCalculable)
        
      CALL map_edges2edges_viacell_2d_constZ_onTriangles_sp( patch_3d, lhs_z_grad_h_sp(:,:), solverCoeffs, lhs_z_e_sp(:,:))

    ENDIF ! l_edge_based
    !---------------------------------------

    !Step 3) Calculate divergence
    ! store the div in lhs for reducing memory and improving performance

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
        CALL div_oce_2D_onTriangles_onBlock_sp(lhs_z_e_sp, patch_2D, solverCoeffs%div_coeff, lhs(:,blockNo), &
          & blockNo=blockNo, start_index=start_index, end_index=end_index)

        !Step 4) Finalize LHS calculations
        DO jc = start_index, end_index
          !lhs(jc,blockNo) =(x(jc,blockNo) - gdt2 * ab_gam * ab_beta * lhs_div_z_c(jc,blockNo)) / gdt2 !rho_sfc(jc,blockNo)*rho_inv
          lhs(jc,blockNo) = x(jc,blockNo) * gdt2_inv - gam_times_beta * lhs(jc,blockNo)
        END DO
      END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO	
    !---------------------------------------

    stop_timer(timer_lhs_sp,2)

  END FUNCTION lhs_surface_height_ab_mim_sp
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computation of new velocity in Adams-Bashforth timestepping.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE calc_normal_velocity_ab_mimetic(patch_3d,ocean_state, op_coeffs)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs
    !
    !  local variables
    INTEGER  :: start_edge_index, end_edge_index
    INTEGER  :: je, jk, blockNo
    REAL(wp) :: gdt_x_ab_beta, one_minus_ab_gam
    REAL(wp) :: z_grad_h(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_grad_h_block(nproma)
    TYPE(t_subset_range), POINTER :: edges_in_domain, owned_cells, owned_edges
    CHARACTER(LEN=*), PARAMETER ::     &
      & method_name='mo_ocean_ab_timestepping_mimetic: calc_normal_velocity_ab_mimetic'
    TYPE(t_patch), POINTER :: patch
    !----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    !-----------------------------------------------------------------------
    patch           => patch_3d%p_patch_2d(1)
    edges_in_domain => patch%edges%in_domain
    owned_cells     => patch%cells%owned
    owned_edges     => patch%edges%owned
    !----------------------------------------------------------------------    
    one_minus_ab_gam = 1.0_wp - ab_gam
    gdt_x_ab_beta =  grav * dtime * ab_beta
    
    IF (iswm_oce == 1) THEN ! shallow water case
    
      z_grad_h(1:nproma,1:patch%nblks_e) = 0.0_wp
      
      ! Step 1) Compute normal derivative of new surface height
      CALL grad_fd_norm_oce_2d_3d( ocean_state%p_prog(nnew(1))%h, &
        & patch,                                                  &
        & op_coeffs%grad_coeff(:,1,:),                            &
        & z_grad_h(:,:))         
        
      IF(MASS_MATRIX_INVERSION_TYPE==MASS_MATRIX_INVERSION_ALLTERMS)THEN
!        write(0,*)'New height grad before',&
!         &maxval(z_grad_h),minval(z_grad_h)
!         !CALL map_edges2edges_viacell_3d_const_z( patch_3D, z_grad_h, op_coeffs, z_grad_h)
!          CALL map_edges2edges_viacell_2D_per_level( patch_3D, z_grad_h, op_coeffs, z_grad_h,1 )
!        write(0,*)'New height grad after',&
!         &maxval(z_grad_h),minval(z_grad_h)

        
        ! Step 2) Calculate the new velocity from the predicted one and the new surface height
        DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
          DO jk = 1, n_zlev
            DO je = start_edge_index, end_edge_index
              IF(patch_3d%lsm_e(je,jk,blockNo) <= sea_boundary)THEN
                ocean_state%p_prog(nnew(1))%vn(je,jk,blockNo)   &
                  & = ocean_state%p_diag%vn_pred_ptp(je,jk,blockNo) &
                  & - gdt_x_ab_beta * z_grad_h(je,blockNo)
              ENDIF
            END DO
          END DO
        END DO
        write(0,*)'New velocity before Inversion',&
        &maxval(ocean_state%p_prog(nnew(1))%vn(:,1,:)),minval(ocean_state%p_prog(nnew(1))%vn(:,1,:)),&
        &maxval(ocean_state%p_diag%vn_pred_ptp(:,1,:)),minval(ocean_state%p_diag%vn_pred_ptp(:,1,:))
        
         ocean_state%p_prog(nnew(1))%vn    &
         & =invert_mass_matrix( patch_3d,   &
         &                     ocean_state,&
         &                     op_coeffs,  &
         &                     ocean_state%p_prog(nnew(1))%vn)              
         write(0,*)'New velocity after Inversion',&
         &maxval(ocean_state%p_prog(nnew(1))%vn(:,1,:)),minval(ocean_state%p_prog(nnew(1))%vn(:,1,:))
        
        
      ELSE
      
        ! Step 2) Calculate the new velocity from the predicted one and the new surface height
        DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
          DO jk = 1, n_zlev
            DO je = start_edge_index, end_edge_index
              IF(patch_3d%lsm_e(je,jk,blockNo) <= sea_boundary)THEN
                ocean_state%p_prog(nnew(1))%vn(je,jk,blockNo)   &
                  & = ocean_state%p_diag%vn_pred(je,jk,blockNo) &
                  & - gdt_x_ab_beta * z_grad_h(je,blockNo)
              ENDIF
            END DO
          END DO
        END DO
      ENDIF
      
      
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
        DO je = start_edge_index, end_edge_index
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
            ocean_state%p_diag%vn_time_weighted(je,jk,blockNo)         &
              & = ab_gam*ocean_state%p_prog(nnew(1))%vn(je,jk,blockNo) &
              & + (1.0_wp -ab_gam)*ocean_state%p_prog(nold(1))%vn(je,jk,blockNo)          
          END DO
        END DO
      END DO

    ELSE !real 3d case
      
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,start_edge_index,end_edge_index, je, jk, z_grad_h_block) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
        ! Step 1) Compute normal derivative of new surface height
        CALL grad_fd_norm_oce_2d_onBlock(            &
          & ocean_state%p_prog(nnew(1))%h,           &
          & patch,                                   &
          & op_coeffs%grad_coeff(:,1, blockNo),      &
          & z_grad_h_block(:),                       & ! we only need one block,
          & start_edge_index, end_edge_index, blockNo)

        ! Step 2) Calculate the new velocity from the predicted one and the new surface height
        DO je = start_edge_index, end_edge_index
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
             ocean_state%p_prog(nnew(1))%vn(je,jk,blockNo)  &
             & = (ocean_state%p_diag%vn_pred(je,jk,blockNo) &
             & - gdt_x_ab_beta * z_grad_h_block(je))

          END DO          
        END DO

        DO je = start_edge_index, end_edge_index
          DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
            ocean_state%p_diag%vn_time_weighted(je,jk,blockNo)      &
              & = ab_gam           * ocean_state%p_prog(nnew(1))%vn(je,jk,blockNo) &
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
!     IF (.NOT.l_rigid_lid) THEN
!       CALL dbg_print('NormVel: grad h-new'      ,z_grad_h                    ,str_module,idt_src, in_subset=owned_edges)
!     END IF
    CALL dbg_print('NormVel: vn_time_weighted'  ,ocean_state%p_diag%vn_time_weighted,str_module,idt_src, in_subset=owned_edges)
    CALL dbg_print('NormVel: vn_change'         ,ocean_state%p_prog(nnew(1))%vn - &
      & ocean_state%p_prog(nold(1))%vn     ,str_module,idt_src, in_subset=owned_edges)
    idt_src=2  ! outputm print level (1-5, fix)
    CALL dbg_print('NormVel: vn_new'            ,ocean_state%p_prog(nnew(1))%vn     ,str_module,idt_src, in_subset=owned_edges)
    !---------------------------------------------------------------------
    
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
    
    TYPE(t_patch_3d), TARGET                :: patch_3d       ! patch on which computation is performed
    TYPE(t_hydro_ocean_state)               :: ocean_state
    TYPE(t_operator_coeff),INTENT(in)       :: op_coeffs
    !
    !
    ! Local variables
    INTEGER :: jc, jk, blockNo, je
    INTEGER :: z_dolic
    INTEGER :: start_index, end_index
    REAL(wp) :: z_c(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_abort
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells
    REAL(wp) ::  minmaxmean(3)
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp),  POINTER  :: vertical_velocity(:,:,:)
    CHARACTER(len=*), PARAMETER :: method_name='mo_ocean_ab_timestepping_mimetic:alc_vert_velocity_mim_bottomup'
    !-----------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain  => patch_2D%cells%in_domain
    all_cells        => patch_2D%cells%all
    edges_in_domain  => patch_2D%edges%in_domain
    vertical_velocity=> ocean_state%p_diag%w

    ! due to nag -nan compiler-option:
    ! vertical_velocity(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks) = 0.0_wp
    !------------------------------------------------------------------
    ! Step 1) Calculate divergence of horizontal velocity at all levels
    !------------------------------------------------------------------
 !ocean_state%p_diag%vn_time_weighted(:,:,:)=0.0_wp   
    
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
          !use bottom boundary condition for vertical velocity at bottom
          !of prism
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

        IF( patch_2d%cells%max_connectivity == 3 )THEN

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
          DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
            CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
        
            CALL div_oce_3D_onTriangles_onBlock(                 &
              & ocean_state%p_diag%mass_flx_e,                   &
              & patch_3D,op_coeffs%div_coeff,                    &
              & ocean_state%p_diag%div_mass_flx_c(:,:,blockNo),  &
              & blockNo=blockNo, start_index=start_index, end_index=end_index,      &
              & start_level=1, end_level=n_zlev)
        
            DO jc = start_index, end_index          
              !use bottom boundary condition for vertical velocity at bottom of prism
              ! this should be awlays zero
              ! vertical_velocity(jc,z_dolic+1,blockNo)=0.0_wp
              DO jk = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), 1, -1
                vertical_velocity(jc,jk,blockNo) &
                  &= vertical_velocity(jc,jk+1,blockNo) - ocean_state%p_diag%div_mass_flx_c(jc,jk,blockNo)
              END DO
            END DO
        
          END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO 


        ELSE
		  
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
        
          CALL div_oce_3D_general_onBlock(                     &
            & ocean_state%p_diag%mass_flx_e,                   &
            & patch_3D,op_coeffs%div_coeff,                    &
            & ocean_state%p_diag%div_mass_flx_c(:,:,blockNo),  &
            & blockNo=blockNo, start_index=start_index, end_index=end_index,      &
            & start_level=1, end_level=n_zlev)
        
          DO jc = start_index, end_index

            !use bottom boundary condition for vertical velocity at bottom of prism
            ! this should be awlays zero
            ! vertical_velocity(jc,z_dolic+1,blockNo)=0.0_wp
            DO jk = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), 1, -1
              vertical_velocity(jc,jk,blockNo) &
              &= vertical_velocity(jc,jk+1,blockNo) - ocean_state%p_diag%div_mass_flx_c(jc,jk,blockNo)
            END DO
          END DO
        
        END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO 
		  		  
      ENDIF	!patch_2d%cells%max_connectivity  

    ENDIF  !  (l_EDGE_BASED)
    
    IF(l_rigid_lid)THEN
      vertical_velocity(:,1,:) = 0.0_wp
    ENDIF

    CALL sync_patch_array(sync_c,patch_2D,vertical_velocity)
    
    !CALL map_scalar_prismtop2center(patch_3d, vertical_velocity, op_coeffs, ocean_state%p_diag%w_prismcenter)
    
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
      CALL dbg_print('Vert veloc: w', &
        & vertical_velocity, str_module,idt_src, in_subset=cells_in_domain)
      
      CALL dbg_print('after cont-correct: h-new',ocean_state%p_prog(nnew(1))%h(:,:) ,str_module,idt_src, &
        & in_subset=cells_in_domain)
      CALL dbg_print('after cont-correct: vol_h', &
        & patch_3d%p_patch_2d(n_dom)%cells%area(:,:) * ocean_state%p_prog(nnew(1))%h(:,:), &
        & str_module,idt_src, in_subset=cells_in_domain)
!      minmaxmean(:) = global_minmaxmean(values=ocean_state%p_prog(nnew(1))%h(:,:), in_subset=cells_in_domain)
!      IF (my_process_is_stdio()) THEN
!        IF (minmaxmean(1) + patch_3D%p_patch_1D(1)%del_zlev_m(1) <= min_top_height) &
!          CALL finish(method_name, "height below min_top_height")
!      ENDIF
      !---------------------------------------------------------------------

    ENDIF
    !-----------------------------------------------------

    
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
              & CALL finish("lhs_surface_height_ab_mim", "lhs(jc,blockNo) /= 0 on land")
          ENDIF
        END DO
      END DO
    ENDIF
    !---------------------------------------------------------------------
    
  END SUBROUTINE calc_vert_velocity_mim_bottomup
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !!  The result is NOT synced. Should be done in the calling method if required
  FUNCTION invert_mass_matrix(patch_3d, ocean_state, op_coeffs, rhs_e) result(inv_flip_flop_e)

    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET      :: ocean_state
    TYPE(t_operator_coeff),INTENT(in)      :: op_coeffs
    REAL(wp)      :: rhs_e(:,:,:)!(nproma,n_zlev,patch_2D%nblks_e)
    REAL(wp)      :: inv_flip_flop_e(SIZE(rhs_e,1),SIZE(rhs_e,2),SIZE(rhs_e,3))
    !
    !LOCAL VARIABLES
    !TYPE(t_patch), TARGET :: patch_2D    
    INTEGER,PARAMETER :: nmax_iter= 800 ! maximum number of iterations
    REAL(wp) :: zimpl_coeff = 1.0_wp    !COEFF has to be set appropriately !!!!
    INTEGER  :: n_iter                  ! number of iterations
    REAL(wp) :: tolerance               ! (relative or absolute) tolerance
    REAL(wp) :: z_residual(nmax_iter)
    LOGICAL  :: lmax_iter               ! true if reached m iterations
    REAL(wp) :: rhstemp(nproma,patch_3d%p_patch_2D(1)%nblks_e)
    INTEGER  :: jk
    INTEGER  :: iter_sum    
    LOGICAL  :: lprecon         = .FALSE.
    REAL(wp) :: zresidual(nmax_iter)    ! norms of the residual (convergence history);an argument of dimension at least m is required
    REAL(wp) :: residual_norm
    REAL(wp) :: relative_tolerance, absolute_tolerance               ! (relative or absolute) tolerance
    LOGICAL  :: maxIterations_isReached, use_absolute_solver_tolerance     ! true if reached m iterations
    INTEGER  :: gmres_restart_iterations
    REAL(wp) :: vol_h(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) ::  minmaxmean(3)
    TYPE(t_subset_range), POINTER :: all_cells, all_edges, owned_cells, owned_edges
    CHARACTER(LEN=max_char_length) :: string

    !REAL(sp) :: zresidual_sp(nmax_iter)    ! norms of the residual (convergence history);an argument of dimension at least m is required
    !REAL(sp) :: residual_norm_sp
    !REAL(sp) :: z_h_c_sp(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(sp) :: h_old_sp(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(sp) :: rhs_sfc_eq_sp(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    !-----------------------------------------------------------------------
    
    tolerance                = MassMatrix_solver_tolerance
    inv_flip_flop_e(:,:,:)   = 0.0_wp
    use_absolute_solver_tolerance=.true.
	
	
    DO jk=1, n_zlev
    
        iter_sum=0
        gmres_restart_iterations = 0
        iter_sum                 = 0
        maxIterations_isReached  = .true.
        IF (use_absolute_solver_tolerance) THEN          
          absolute_tolerance      = tolerance!solver_tolerance
        ELSE
          relative_tolerance      = tolerance!solver_tolerance
          absolute_tolerance      = 0.0_wp
        ENDIF
        
!         write(0,*) "Tolerence:", absolute_tolerance, relative_tolerance
  
        SELECT CASE (select_solver)
        !-----------------------------------------------------------------------------------------
        CASE (select_gmres)
          CALL gmres_oce_e2e(                      &
            & x = inv_flip_flop_e(:,jk,:),         &
            & lhs = lhs_primal_flip_flop,          &
            & patch_3d = patch_3d,                 &
            & level = jk,                          &
            & p_op_coeff = op_coeffs,              &
            & b = rhs_e(:,jk,:),                   &
            & tolerance = tolerance,               &
            & abstol = use_absolute_solver_tolerance,&
            & m = nmax_iter,                       &
            & maxiterex = maxIterations_isReached, &
            & niter = n_iter,                      &
            & res = zresidual)

          IF (maxIterations_isReached) THEN
            WRITE(string,'(a,i4,a,e28.20)') &
              & 'gmres old iteration =', n_iter,', residual =', ABS(zresidual(n_iter))
            CALL message('invert_mass_matrix',TRIM(string))
            CALL finish('GMRES_oce_old: solver surface equation: ','NOT YET CONVERGED !!')
          ELSE
            ! output print level idt_src used for ocean gmres output with call message:
            IF(n_iter==0)n_iter=1
            idt_src=1
            IF (idbg_mxmn >= idt_src) THEN
              WRITE(string,'(a,i4,a,e28.20)') &
                & 'gmres old iteration =', n_iter,', residual =', ABS(zresidual(n_iter))
              CALL message('invert_mass_matrix',TRIM(string))
            ENDIF
          ENDIF
    
        CASE (select_restart_gmres)
          DO WHILE(maxIterations_isReached .AND. gmres_restart_iterations < solver_max_restart_iterations)

            CALL ocean_restart_gmres_e2e(     &
              & inv_flip_flop_e(:,jk,:),      &  ! arg 1 of lhs. x input is the first guess.
              & lhs_primal_flip_flop,         &  ! function calculating l.h.s.
              & patch_3d,                     &  ! arg 3 of lhs
              & op_coeffs,                    &
              & jk,                           &
              & rhs_e(:,jk,:),                &  ! right hand side as input
              & absolute_tolerance,                &     ! inout, if > 0 then is used as absolute_tolerance, out otherwise
              & relative_tolerance,                &
              & solver_max_iter_per_restart,       &  ! max. # of iterations to do
              & maxIterations_isReached,                     &  ! out: .true. = not converged
              & n_iter,                        &  ! out: # of iterations done
              & zresidual)

            IF(n_iter==0) THEN
              residual_norm = 0.0_wp
            ELSE
              residual_norm =  ABS(zresidual(n_iter))
            ENDIF

            iter_sum = iter_sum + n_iter
            ! output print level idt_src used for ocean_restart_gmres output with call message:
            idt_src=2
            IF (idbg_mxmn >= idt_src) THEN
              WRITE(string,'(a,i4,a,e28.20)') &
                & 'ocean_restart_gmres iteration =', n_iter,', residual =', residual_norm
              CALL message('invert_mass_matrix',TRIM(string))
            ENDIF

            gmres_restart_iterations = gmres_restart_iterations + 1

          END DO ! WHILE(tolerance >= solver_tolerance)
          
          idt_src=1
          IF (idbg_mxmn >= idt_src) THEN
            WRITE(string,'(a,i4,a,e28.20)') &
              & 'ocean_restart_gmres iteration =', iter_sum,', residual =', residual_norm
            CALL message('invert_mass_matrix',TRIM(string))
          ENDIF

            gmres_restart_iterations = gmres_restart_iterations + 1
        !-----------------------------------------------------------------------------------------
        END SELECT
       
      END DO    
  END FUNCTION invert_mass_matrix
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !!  results is valid only in in_domain edges
  FUNCTION lhs_primal_flip_flop( x, patch_3d, op_coeffs,jk) result(llhs)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp),INTENT(inout)               :: x(:,:)
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs
    INTEGER                              :: jk
    REAL(wp)                             :: llhs(SIZE(x,1), SIZE(x,2))!,SIZE(x,3))
    
    !local variables
    !REAL(wp) :: z_e(SIZE(x,1), SIZE(x,2))
    !-----------------------------------------------------------------------
    !edges_in_domain => patch_2D%edges%in_domain
    !CALL finish("lhs_primal_flip_flop", "not implemented")
    
    !llhs(:,:) = 0.5_wp*x(:,:)
    CALL sync_patch_array(sync_e, patch_3d%p_patch_2d(1), x )
    
    CALL map_edges2edges_viacell_2D_per_level( patch_3d, &
      & x(:,:), &
      & op_coeffs, llhs(:,:), jk)  
    !write(*,*)'max/min v:PTPv:', maxval(x(:,:)),minval(x(:,:)),maxval(llhs(:,:)),minval(llhs(:,:))
    
  END FUNCTION lhs_primal_flip_flop
  !--------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE jacobi_precon( p_jp, patch_3d, op_coeffs,thick_e) !RESULT(p_jp)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp),INTENT(inout)                        :: p_jp(:,:)    ! inout for sync, dimension: (nproma,patch%alloc_cell_blocks)
    TYPE(t_operator_coeff),INTENT(in)             :: op_coeffs
    REAL(wp),INTENT(in)                           :: thick_e(:,:)
    !
    ! Left-hand side calculated from iterated height
    !REAL(wp) :: p_jp(SIZE(x,1), SIZE(x,2))  ! (nproma,patch%alloc_cell_blocks)
    !
    ! local variables
    REAL(wp) :: gdt2
    INTEGER :: start_index, end_index
    INTEGER :: jc, blockNo  !, je
    !REAL(wp) :: z1,z2,z3
    INTEGER :: edge_1_idx, edge_1_blk, edge_2_idx, edge_2_blk, edge_3_idx, edge_3_blk
    REAL(wp) :: jacobi_diagnostic(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_cells  !, cells_in_domain, all_edges
    TYPE(t_patch), POINTER :: patch     ! patch on which computation is performed
    !-----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start - iteration by ocean_gmres')
    patch =>  patch_3d%p_patch_2d(1)
    !write(*,*)'inside jacobi'
    WRITE(*,*)'residual before',MAXVAL(p_jp),MINVAL(p_jp)!,maxvalp_jp),minval(p_jp)
    all_cells => patch%cells%ALL
    !cells_in_domain => patch%cells%in_domain
    !all_edges => patch%edges%all

    !p_jp(1:nproma,1:patch%alloc_cell_blocks)  = 0.0_wp

    gdt2 = grav*(dtime)**2


    CALL sync_patch_array(sync_c, patch, p_jp )


    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        !IF ( v_base%lsm_c(jc,1,blockNo) <= sea_boundary ) THEN
        !IF(patch_3d%lsm_c(jc,1,blockNo)>=MIN_DOLIC)THEN
        edge_1_idx = patch%cells%edge_idx(jc,blockNo,1)
        edge_1_blk = patch%cells%edge_blk(jc,blockNo,1)
        edge_2_idx = patch%cells%edge_idx(jc,blockNo,2)
        edge_2_blk = patch%cells%edge_blk(jc,blockNo,2)
        edge_3_idx = patch%cells%edge_idx(jc,blockNo,3)
        edge_3_blk = patch%cells%edge_blk(jc,blockNo,3)
        !
        !         jacobi_diagnostic(jc,blockNo) = patch%cells%area(jc,blockNo)-&
        !                &gdt2*ab_gam*ab_beta*((op_coeffs%div_coeff(jc,1,blockNo,1)&!*thick_e(edge_1_idx,edge_1_blk)&
        !                &*op_coeffs%grad_coeff(edge_1_idx,1,edge_1_blk)&
        !                &+&
        !                &op_coeffs%div_coeff(jc,1,blockNo,2)&!*thick_e(edge_2_idx,edge_2_blk)&
        !                &*op_coeffs%grad_coeff(edge_2_idx,1,edge_2_blk)&
        !                &+&
        !                &op_coeffs%div_coeff(jc,1,blockNo,3)&!*thick_e(edge_3_idx,edge_3_blk)&
        !                &*op_coeffs%grad_coeff(edge_3_idx,1,edge_3_blk)))&
        !                &*patch%cells%area(jc,blockNo)
        !
        ! z1=gdt2*ab_gam*ab_beta*(&
        !   & op_coeffs%div_coeff(jc,1,blockNo,1)*thick_e(edge_1_idx,edge_1_blk)&
        !   &*op_coeffs%grad_coeff(edge_1_idx,1,edge_1_blk))
        !
        !
        ! z2=gdt2*ab_gam*ab_beta*(op_coeffs%div_coeff(jc,1,blockNo,2)*thick_e(edge_2_idx,edge_2_blk)&
        !    &*op_coeffs%grad_coeff(edge_2_idx,1,edge_2_blk))
        !
        ! z3=gdt2*ab_gam*ab_beta*(op_coeffs%div_coeff(jc,1,blockNo,3)*thick_e(edge_3_idx,edge_3_blk)&
        !   &*op_coeffs%grad_coeff(edge_3_idx,1,edge_3_blk))
        !
        ! jacobi_diagnostic(jc,blockNo)=1.0_wp-sqrt(z1*z1+z2*z2+z3*z3)gdt2

        !         jacobi_diagnostic(jc,blockNo) = (1.0_wp-&
        !                &gdt2*ab_gam*ab_beta*(&
        !                & op_coeffs%div_coeff(jc,1,blockNo,1)*thick_e(edge_1_idx,edge_1_blk)&
        !                &*op_coeffs%grad_coeff(edge_1_idx,1,edge_1_blk)&
        !                &+&
        !                &op_coeffs%div_coeff(jc,1,blockNo,2)*thick_e(edge_2_idx,edge_2_blk)&
        !                &*op_coeffs%grad_coeff(edge_2_idx,1,edge_2_blk)&
        !                &+&
        !                &op_coeffs%div_coeff(jc,1,blockNo,3)*thick_e(edge_3_idx,edge_3_blk)&
        !                &*op_coeffs%grad_coeff(edge_3_idx,1,edge_3_blk)))*patch%cells%area(jc,blockNo)/gdt2
        !
        jacobi_diagnostic(jc,blockNo)=1.0_wp- &
          & gdt2*ab_gam*ab_beta*(patch%edges%primal_edge_length(edge_1_idx,edge_1_blk)&
          & *patch%edges%inv_dual_edge_length(edge_1_idx,edge_1_blk)&
          & +patch%edges%primal_edge_length(edge_2_idx,edge_2_blk)&
          & *patch%edges%inv_dual_edge_length(edge_2_idx,edge_2_blk)&
          & +patch%edges%primal_edge_length(edge_3_idx,edge_3_blk) &
          & *patch%edges%inv_dual_edge_length(edge_1_idx,edge_1_blk))&
          & /(patch%cells%area(jc,blockNo)*gdt2)
        !  jacobi_diagnostic(jc,blockNo) = (1.0_wp - gdt2*ab_gam*ab_beta*(patch%edges%primal_edge_length(edge_1_idx,edge_1_blk) &
        !     !          & *patch%edges%inv_dual_edge_length(edge_1_idx,edge_1_blk)&
        !               &+patch%edges%primal_edge_length(edge_2_idx,edge_2_blk) &
        !     !          & *patch%edges%inv_dual_edge_length(edge_2_idx,edge_2_blk)&
        !               &+ patch%edges%primal_edge_length(edge_3_idx,edge_3_blk))) &
        !     !          & *patch%edges%inv_dual_edge_length(edge_3_idx,edge_3_blk)))!&
        !               &/(patch%cells%area(jc,blockNo)*gdt2)
        !        IF(jacobi_diagnostic/=0.0_wp)THEN
        p_jp(jc,blockNo) = p_jp(jc,blockNo)*jacobi_diagnostic(jc,blockNo)

        !        ENDIF
        !ENDIF
      END DO
    END DO
    WRITE(*,*)'diag element',MAXVAL(jacobi_diagnostic),MINVAL(jacobi_diagnostic)
    WRITE(*,*)'residual after',MAXVAL(p_jp),MINVAL(p_jp)!,maxvalp_jp),minval(p_jp)
  END SUBROUTINE jacobi_precon
  !-------------------------------------------------------------------------

END MODULE mo_ocean_ab_timestepping_mimetic
