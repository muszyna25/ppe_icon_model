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
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_oce_ab_timestepping_mimetic

  USE mo_kind,                      ONLY: wp, sp
  USE mo_parallel_config,           ONLY: nproma, l_fast_sum
  USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  USE mo_sync,                      ONLY: sync_e, sync_c, sync_patch_array
  USE mo_impl_constants,            ONLY: sea_boundary, &  !  sea,                          &
    & max_char_length, min_dolic
  USE mo_dbg_nml,                   ONLY: idbg_mxmn
  USE mo_ocean_nml,                 ONLY: n_zlev, solver_tolerance, l_inverse_flip_flop,    &
    & ab_const, ab_beta, ab_gam, iswm_oce,                &
    & expl_vertical_velocity_diff, iforc_oce,             &
    & no_tracer, l_rigid_lid, l_edge_based,               &
    & use_absolute_solver_tolerance,                      &
    & solver_max_restart_iterations,                      &
    & solver_max_iter_per_restart, dhdtw_abort,           &
    & forcing_enable_freshwater, select_solver,           &
    & select_restart_gmres, select_gmres,                 &
    & use_continuity_correction,                          &
    & select_restart_mixedPrecision_gmres,                &
    & solver_max_iter_per_restart_sp,                     &
    & solver_tolerance_sp
  
  USE mo_run_config,                ONLY: dtime, ltimer, debug_check_level
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_ab_expl,           &
    & timer_ab_rhs4sfc, timer_lhs, timer_lhs_sp
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_physical_constants,        ONLY: grav,rho_inv
  USE mo_ocean_initialization,      ONLY: is_initial_timestep
  USE mo_oce_types,                 ONLY: t_hydro_ocean_state, t_hydro_ocean_diag
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,            ONLY: t_external_data
  USE mo_ocean_gmres,               ONLY: ocean_restart_gmres, gmres_oce_old, gmres_oce_e2e, &
    & ocean_restart_gmres_singlePrecesicion
  USE mo_exception,                 ONLY: message, finish, message_text
  USE mo_util_dbg_prnt,             ONLY: dbg_print, debug_print_MaxMinMean
  USE mo_oce_boundcond,             ONLY: bot_bound_cond_horz_veloc, top_bound_cond_horz_veloc
  USE mo_oce_thermodyn,             ONLY: calc_density, calc_internal_press
  USE mo_oce_physics,               ONLY: t_ho_params
  USE mo_sea_ice_types,             ONLY: t_sfc_flx
  USE mo_scalar_product,            ONLY: map_edges2edges_viacell_3d, & ! map_cell2edges_3D,&
    & calc_scalar_product_veloc_3d,&
    & map_edges2edges_viacell_3d_const_z, map_edges2edges_viacell_2d_constZ_sp
  USE mo_oce_math_operators,        ONLY: div_oce_3d, grad_fd_norm_oce_3d,&
    & grad_fd_norm_oce_2d_3d, grad_fd_norm_oce_2d_3d_sp, calculate_thickness, &
    & div_oce_2d_sp
  USE mo_oce_veloc_advection,       ONLY: veloc_adv_horz_mimetic, veloc_adv_vert_mimetic
  
  USE mo_oce_diffusion,             ONLY: velocity_diffusion,&
    & velocity_diffusion_vert_explicit,  &
    & velocity_diffusion_vertical_implicit
  USE mo_oce_types,                 ONLY: t_operator_coeff, t_solverCoeff_singlePrecision
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_grid_config,               ONLY: n_dom
!  USE mo_parallel_config,           ONLY: p_test_run
  USE mo_mpi,                       ONLY: my_process_is_stdio, get_my_global_mpi_id ! my_process_is_mpi_parallel
  USE mo_statistics,                ONLY: global_minmaxmean
  IMPLICIT NONE
  
  PRIVATE
  
  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  !
  ! PUBLIC INTERFACE
  !
  PUBLIC :: solve_free_sfc_ab_mimetic
  PUBLIC :: calc_normal_velocity_ab_mimetic
  PUBLIC :: calc_vert_velocity_mim_bottomup
  PUBLIC :: init_ho_lhs_fields_mimetic
  !
  ! Private implemenation
  !
  PRIVATE :: fill_rhs4surface_eq_ab
  PRIVATE :: calculate_explicit_term_ab   ! calc_velocity_predictor
  PRIVATE :: lhs_surface_height_ab_mim
  PRIVATE :: inverse_primal_flip_flop
  PRIVATE :: jacobi_precon
  
  INTEGER, PARAMETER :: top=1
  CHARACTER(LEN=12)  :: str_module = 'oceSTEPmimet'  ! Output of module for 1 line debug
  INTEGER :: idt_src    = 1               ! Level of detail for 1 line debug
  
  ! TRUE=staggering between thermodynamic and dynamic part, offset of half timestep
  ! between dynamic and thermodynamic variables thermodynamic and dnamic variables are colocated in time
  LOGICAL, PUBLIC,PARAMETER :: l_staggered_timestep = .FALSE.
  
  
  ! these are allocated once for efficiency and used only by the lhs for the solver
  REAL(wp), ALLOCATABLE, TARGET :: lhs_result(:,:)  ! (nproma,patch%alloc_cell_blocks)
  REAL(wp), ALLOCATABLE :: lhs_z_grad_h(:,:)
  REAL(wp), ALLOCATABLE :: lhs_z_e     (:,:)
  ! REAL(wp), ALLOCATABLE :: lhs_z_e_top (:,:)
  ! the same as above in single precision
  REAL(sp), ALLOCATABLE, TARGET :: lhs_result_sp(:,:)  ! (nproma,patch%alloc_cell_blocks)
  REAL(sp), ALLOCATABLE :: lhs_z_grad_h_sp(:,:)
  REAL(sp), ALLOCATABLE :: lhs_z_e_sp     (:,:)


  ! TYPE(t_cartesian_coordinates), ALLOCATABLE :: lhs_z_grad_h_cc(:,:)
  REAL(wp), PARAMETER ::  min_top_height = 0.05_wp ! we have to have at least 5cm water on top of sea cells
  
CONTAINS
  
  
  !-------------------------------------------------------------------------
  !>
  !! !  Solves the free surface equation.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize_Used>
  SUBROUTINE solve_free_sfc_ab_mimetic(patch_3d, ocean_state, p_ext_data, p_sfc_flx, &
    & p_phys_param, timestep, op_coeffs, solverCoeff_sp)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout)       :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET             :: ocean_state
    TYPE(t_external_data), TARGET, INTENT(in)     :: p_ext_data
    TYPE(t_sfc_flx), INTENT(inout)                :: p_sfc_flx
    TYPE (t_ho_params)                            :: p_phys_param
    INTEGER, INTENT(in)                           :: timestep
    TYPE(t_operator_coeff)                        :: op_coeffs
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp
    !
    !Local variables
    !
    INTEGER,PARAMETER :: nmax_iter   = 800      ! maximum number of iterations
    REAL(wp) :: tolerance =0.0_wp               ! (relative or absolute) tolerance
    INTEGER :: n_iter                          ! actual number of iterations
    INTEGER :: iter_sum                        ! sum of iterations
    INTEGER :: jc,jb,je   ! ,jk,il_v1,il_v2,ib_v1,ib_v2
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: edge_start_idx, edge_end_idx
    REAL(wp) :: z_h_c(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_h_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    LOGICAL :: lprecon         = .FALSE.
    ! REAL(wp) :: z_implcoeff
    REAL(wp) :: zresidual(nmax_iter)    ! norms of the residual (convergence history);an argument of dimension at least m is required
    REAL(wp) :: residual_norm
    LOGICAL :: l_maxiter     ! true if reached m iterations
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

    CHARACTER(len=*), PARAMETER :: method_name='mo_oce_ab_timestepping_mimetic:solve_free_sfc_ab_mimetic'
    !-------------------------------------------------------------------------------
    patch_2D   => patch_3d%p_patch_2d(1)
    all_cells    => patch_2D%cells%ALL
    all_edges    => patch_2D%edges%ALL
    owned_cells  => patch_2D%cells%owned
    owned_edges  => patch_2D%edges%owned
    !-------------------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    tolerance                            = solver_tolerance
    z_h_c(1:nproma,1:patch_2D%alloc_cell_blocks) = 0.0_wp
    z_h_e(1:nproma,1:patch_2D%nblks_e)           = 0.0_wp
    
    CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_prog(nold(1))%h)
    CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_prog(nold(1))%vn)
    
    IF (is_initial_timestep(timestep) ) THEN
      IF (l_staggered_timestep ) CALL calc_scalar_product_veloc_3D( patch_3d,&
                                                                  & ocean_state%p_prog(nold(1))%vn,&
                                                                  & ocean_state%p_prog(nold(1))%vn,&
                                                                  & ocean_state%p_diag,            &
                                                                  & op_coeffs)

    ENDIF

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
    idt_src=1
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
    CALL bot_bound_cond_horz_veloc(patch_3d, ocean_state, p_phys_param, op_coeffs)

    ! CALL dbg_print('bc_top_vn', ocean_state%p_aux%bc_top_vn, str_module, idt_src,  in_subset=owned_edges)
    ! CALL dbg_print('bc_bot_vn', ocean_state%p_aux%bc_bot_vn, str_module, idt_src,  in_subset=owned_edges)
    !---------------------------------------------------------------------

    
    IF (ltimer) CALL timer_start(timer_ab_expl)
    CALL calculate_explicit_term_ab(patch_3d, ocean_state, p_phys_param, &
      & is_initial_timestep(timestep), op_coeffs)
    IF (ltimer) CALL timer_stop(timer_ab_expl)
    
    IF(.NOT.l_rigid_lid)THEN
      
      ! Calculate RHS of surface equation
      IF (ltimer) CALL timer_start(timer_ab_rhs4sfc)
      CALL fill_rhs4surface_eq_ab(patch_3d, ocean_state, p_sfc_flx, op_coeffs)
      IF (ltimer) CALL timer_stop(timer_ab_rhs4sfc)
      
      
      ! Solve surface equation with ocean_gmres solver
      z_h_c =0.0_wp!ocean_state%p_prog(nold(1))%h! 0.0_wp !potentially better choice: ocean_state%p_prog(nold(1))%h
      
      !The lhs needs different thicknesses at edges for 3D and SWE (including bathymetry)
      !IF( iswm_oce == 1 ) THEN
      !  z_h_e = ocean_state%p_diag%h_e! ocean_state%p_diag%thick_e
      !  ELSEIF( iswm_oce /= 1 ) THEN
      !  z_h_e = ocean_state%p_diag%h_e         ! #slo# 2011-02-21 bugfix (for mimetic only)
      !ENDIF
      
      CALL sync_patch_array(sync_e, patch_2D, z_h_e)
      CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_diag%thick_c)
      CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_prog(nold(1))%h)
      
      CALL dbg_print('bef ocean_gmres: h-old',ocean_state%p_prog(nold(1))%h(:,:) ,str_module,idt_src, in_subset=owned_cells)

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
              & tolerance,                 &  ! relative tolerance
              & use_absolute_solver_tolerance,     &  ! NOT absolute tolerance
              & nmax_iter,                 &  ! max. # of iterations to do
              & l_maxiter,                 &  ! out: .true. = not converged
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
              & tolerance,                 &  ! relative tolerance
              & use_absolute_solver_tolerance,                 &  ! NOT absolute tolerance
!              & .FALSE.,                   &  ! NOT absolute tolerance
              & nmax_iter,                 &  ! max. # of iterations to do
              & l_maxiter,                 &  ! out: .true. = not converged
              & n_iter,                    &  ! out: # of iterations done
              & zresidual )
          ENDIF
          
          IF (l_maxiter) THEN
            CALL finish('GMRES_oce_old: solver surface equation: ','NOT YET CONVERGED !!')
          ELSE
            ! output print level idt_src used for ocean gmres output with call message:
            IF(n_iter==0)n_iter=1
            idt_src=1
            IF (idbg_mxmn >= idt_src) THEN
              WRITE(string,'(a,i4,a,e28.20)') &
                & 'iteration =', n_iter,', residual =', ABS(zresidual(n_iter))
              CALL message('GMRES_oce_old: surface height',TRIM(string))
            ENDIF
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
        residual_norm = solver_tolerance + 1.0_wp
        gmres_restart_iterations = 0
        iter_sum                 = 0
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
            & solver_tolerance,              &  ! tolerance
            & use_absolute_solver_tolerance, &  ! use absolute tolerance = true
            & solver_max_iter_per_restart,   &  ! max. # of iterations to do
            & l_maxiter,                     &  ! out: .true. = not converged
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
          CALL message('ocean_restart_gmres: surface height',TRIM(string))
        ENDIF
        
        IF (residual_norm > solver_tolerance) &
          & CALL finish(method_name,'NOT YET CONVERGED !!')
        
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
            & l_maxiter,                            &  ! out: .true. = not converged
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
            & solver_tolerance,              &  ! tolerance
            & use_absolute_solver_tolerance, &  ! use absolute tolerance = true
            & solver_max_iter_per_restart,   &  ! max. # of iterations to do
            & l_maxiter,                     &  ! out: .true. = not converged
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

        IF (residual_norm > solver_tolerance) &
          & CALL finish(method_name,'NOT YET CONVERGED !!')

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
      IF (my_process_is_stdio()) THEN
        IF (minmaxmean(1) + patch_3D%p_patch_1D(1)%del_zlev_m(1) <= min_top_height) &
          CALL finish(method_name, "height below min_top_height")
      ENDIF
      !---------------------------------------------------------------------
      
    ENDIF  ! l_rigid_lid
    
    ! write(0,*) "solve_free_sfc_ab_mimetic: sum(h)=", SUM(ocean_state%p_prog(nnew(1))%h(:,:))
    
  END SUBROUTINE solve_free_sfc_ab_mimetic
  !-------------------------------------------------------------------------
  
  
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
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc, jb  !, je
    !REAL(wp) :: z1,z2,z3
    INTEGER :: edge_1_idx, edge_1_blk, edge_2_idx, edge_2_blk, edge_3_idx, edge_3_blk
    REAL(wp) :: p_diag(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
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
    
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
      DO jc = i_startidx, i_endidx
        !IF ( v_base%lsm_c(jc,1,jb) <= sea_boundary ) THEN
        !IF(patch_3d%lsm_c(jc,1,jb)>=MIN_DOLIC)THEN
        edge_1_idx = patch%cells%edge_idx(jc,jb,1)
        edge_1_blk = patch%cells%edge_blk(jc,jb,1)
        edge_2_idx = patch%cells%edge_idx(jc,jb,2)
        edge_2_blk = patch%cells%edge_blk(jc,jb,2)
        edge_3_idx = patch%cells%edge_idx(jc,jb,3)
        edge_3_blk = patch%cells%edge_blk(jc,jb,3)
        !
        !         p_diag(jc,jb) = patch%cells%area(jc,jb)-&
        !                &gdt2*ab_gam*ab_beta*((op_coeffs%div_coeff(jc,1,jb,1)&!*thick_e(edge_1_idx,edge_1_blk)&
        !                &*op_coeffs%grad_coeff(edge_1_idx,1,edge_1_blk)&
        !                &+&
        !                &op_coeffs%div_coeff(jc,1,jb,2)&!*thick_e(edge_2_idx,edge_2_blk)&
        !                &*op_coeffs%grad_coeff(edge_2_idx,1,edge_2_blk)&
        !                &+&
        !                &op_coeffs%div_coeff(jc,1,jb,3)&!*thick_e(edge_3_idx,edge_3_blk)&
        !                &*op_coeffs%grad_coeff(edge_3_idx,1,edge_3_blk)))&
        !                &*patch%cells%area(jc,jb)
        !
        ! z1=gdt2*ab_gam*ab_beta*(&
        !   & op_coeffs%div_coeff(jc,1,jb,1)*thick_e(edge_1_idx,edge_1_blk)&
        !   &*op_coeffs%grad_coeff(edge_1_idx,1,edge_1_blk))
        !
        !
        ! z2=gdt2*ab_gam*ab_beta*(op_coeffs%div_coeff(jc,1,jb,2)*thick_e(edge_2_idx,edge_2_blk)&
        !    &*op_coeffs%grad_coeff(edge_2_idx,1,edge_2_blk))
        !
        ! z3=gdt2*ab_gam*ab_beta*(op_coeffs%div_coeff(jc,1,jb,3)*thick_e(edge_3_idx,edge_3_blk)&
        !   &*op_coeffs%grad_coeff(edge_3_idx,1,edge_3_blk))
        !
        ! p_diag(jc,jb)=1.0_wp-sqrt(z1*z1+z2*z2+z3*z3)gdt2
        
        !         p_diag(jc,jb) = (1.0_wp-&
        !                &gdt2*ab_gam*ab_beta*(&
        !                & op_coeffs%div_coeff(jc,1,jb,1)*thick_e(edge_1_idx,edge_1_blk)&
        !                &*op_coeffs%grad_coeff(edge_1_idx,1,edge_1_blk)&
        !                &+&
        !                &op_coeffs%div_coeff(jc,1,jb,2)*thick_e(edge_2_idx,edge_2_blk)&
        !                &*op_coeffs%grad_coeff(edge_2_idx,1,edge_2_blk)&
        !                &+&
        !                &op_coeffs%div_coeff(jc,1,jb,3)*thick_e(edge_3_idx,edge_3_blk)&
        !                &*op_coeffs%grad_coeff(edge_3_idx,1,edge_3_blk)))*patch%cells%area(jc,jb)/gdt2
        !
        p_diag(jc,jb)=1.0_wp- &
          & gdt2*ab_gam*ab_beta*(patch%edges%primal_edge_length(edge_1_idx,edge_1_blk)&
          & *patch%edges%inv_dual_edge_length(edge_1_idx,edge_1_blk)&
          & +patch%edges%primal_edge_length(edge_2_idx,edge_2_blk)&
          & *patch%edges%inv_dual_edge_length(edge_2_idx,edge_2_blk)&
          & +patch%edges%primal_edge_length(edge_3_idx,edge_3_blk) &
          & *patch%edges%inv_dual_edge_length(edge_1_idx,edge_1_blk))&
          & /(patch%cells%area(jc,jb)*gdt2)
        !  p_diag(jc,jb) = (1.0_wp - gdt2*ab_gam*ab_beta*(patch%edges%primal_edge_length(edge_1_idx,edge_1_blk) &
        !     !          & *patch%edges%inv_dual_edge_length(edge_1_idx,edge_1_blk)&
        !               &+patch%edges%primal_edge_length(edge_2_idx,edge_2_blk) &
        !     !          & *patch%edges%inv_dual_edge_length(edge_2_idx,edge_2_blk)&
        !               &+ patch%edges%primal_edge_length(edge_3_idx,edge_3_blk))) &
        !     !          & *patch%edges%inv_dual_edge_length(edge_3_idx,edge_3_blk)))!&
        !               &/(patch%cells%area(jc,jb)*gdt2)
        !        IF(p_diag/=0.0_wp)THEN
        p_jp(jc,jb) = p_jp(jc,jb)*p_diag(jc,jb)
        
        !        ENDIF
        !ENDIF
      END DO
    END DO
    WRITE(*,*)'diag element',MAXVAL(p_diag),MINVAL(p_diag)
    WRITE(*,*)'residual after',MAXVAL(p_jp),MINVAL(p_jp)!,maxvalp_jp),minval(p_jp)
  END SUBROUTINE jacobi_precon
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computation of velocity predictor in Adams-Bashforth timestepping.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize_Used>
  SUBROUTINE calculate_explicit_term_ab( patch_3d, ocean_state, p_phys_param,&
    & l_initial_timestep, op_coeffs)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE (t_ho_params)                   :: p_phys_param
    LOGICAL,INTENT(in)                   :: l_initial_timestep
    TYPE(t_operator_coeff)               :: op_coeffs
    !
    !local variables
    !
    INTEGER :: je, jk, jb
    REAL(wp) :: gdt
    REAL(wp) :: z_gradh_e(nproma, patch_3d%p_patch_2d(n_dom)%nblks_e)
    REAL(wp) :: z_e(nproma,n_zlev,patch_3d%p_patch_2d(n_dom)%nblks_e)
    !REAL(wp) :: z_h_old_rho(nproma, patch_3d%p_patch_2d(n_dom)%nblks_c)
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_edges, owned_edges, owned_cells
    INTEGER :: edge_start_idx, edge_end_idx, dolic_e
    TYPE(t_patch), POINTER :: patch_2D
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !  &       routine = ('mo_oce_ab_timestepping_mimetic:calculate_explicit_term_ab')
    !-----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    
    patch_2D      => patch_3d%p_patch_2d(n_dom)
    edges_in_domain => patch_3d%p_patch_2d(n_dom)%edges%in_domain
    all_edges       => patch_3d%p_patch_2d(n_dom)%edges%ALL
    owned_edges     => patch_3d%p_patch_2d(n_dom)%edges%owned
    owned_cells     => patch_3d%p_patch_2d(n_dom)%cells%owned
    
    gdt             = grav * dtime    
    !---------------------------------------------------------------------
    ! STEP 1: horizontal advection
    !---------------------------------------------------------------------
    
    IF(l_initial_timestep)THEN
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
    
    !---------------------------------------------------------------------
    ! STEP 2: compute 3D contributions: gradient of hydrostatic pressure and vertical velocity advection
    !---------------------------------------------------------------------
    
    IF ( iswm_oce /= 1 ) THEN
      ! calculate density from EOS using temperature and salinity at timelevel n
      CALL calc_density( patch_3d,                                    &
        & ocean_state%p_prog(nold(1))%tracer(:,:,:,1:no_tracer),&
        & ocean_state%p_diag%rho(:,:,:) )
      
      ! calculate hydrostatic pressure from density at timelevel nc
      CALL calc_internal_press( patch_3d,                  &  ! in
        & ocean_state%p_diag%rho,                          &  ! in
        & patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c,   &  ! in
        & ocean_state%p_prog(nold(1))%h,                   &  ! in
        & ocean_state%p_diag%press_hyd)                       ! inout
      
      ! calculate gradient of hydrostatic pressure in 3D
      CALL grad_fd_norm_oce_3d( ocean_state%p_diag%press_hyd,  &
        & patch_3d,             &
        & op_coeffs%grad_coeff,  &
        & ocean_state%p_diag%press_grad)
      CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_diag%press_grad)
      
      ! calculate vertical velocity advection
      CALL veloc_adv_vert_mimetic( patch_3d,                 &
        & ocean_state%p_diag,op_coeffs,     &
        & ocean_state%p_diag%veloc_adv_vert )
      
      ! calculate vertical velocity diffusion
      !   For the alternative choice "expl_vertical_velocity_diff==1" see couples of
      !   lines below below
      IF (expl_vertical_velocity_diff==0) THEN
        CALL velocity_diffusion_vert_explicit( patch_3d,     &
          & ocean_state%p_diag,            &
          & ocean_state%p_aux,op_coeffs,  &
          & p_phys_param,           &
          & ocean_state%p_diag%laplacian_vert)
      ENDIF
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
    
    CALL velocity_diffusion(patch_3d,              &
      & ocean_state%p_prog(nold(1))%vn, &
      & p_phys_param,            &
      & ocean_state%p_diag,op_coeffs,  &
      & ocean_state%p_diag%laplacian_horz)
    
    CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_diag%laplacian_horz)

    !---------------------------------------------------------------------
    ! STEP 4: calculate weighted gradient of surface height at previous timestep
    !---------------------------------------------------------------------
    ! zero all for the nag compiler
    z_gradh_e(:,:)  = 0.0_wp
        
    !ocean_state%p_prog(nold(1))%h(1:nproma, 1:patch_3d%p_patch_2d(n_dom)%nblks_c)&
    !&=rho_inv* ocean_state%p_diag%rho(1:nproma,1, 1:patch_3d%p_patch_2d(n_dom)%nblks_c)&  
    !&*ocean_state%p_prog(nold(1))%h(1:nproma, 1:patch_3d%p_patch_2d(n_dom)%nblks_c) 
    
    CALL grad_fd_norm_oce_2d_3d(ocean_state%p_prog(nold(1))%h, &
      & patch_2D,                  &
      & op_coeffs%grad_coeff(:,1,:),  &
      & z_gradh_e(:,:))
    CALL dbg_print('old height gradient'  ,z_gradh_e, str_module,idt_src, in_subset=owned_edges)

    DO jb = owned_edges%start_block, owned_edges%end_block
       CALL get_index_range(edges_in_domain, jb, edge_start_idx, edge_end_idx)
       z_gradh_e(edge_start_idx:edge_end_idx, jb) = &
         & (1.0_wp-ab_beta) * grav * z_gradh_e(edge_start_idx:edge_end_idx, jb)
    ENDDO
    CALL sync_patch_array(sync_e, patch_2D, z_gradh_e(:,:))

    !---------------------------------------------------------------------
    ! STEP 5:
    !---------------------------------------------------------------------
    
    IF (l_inverse_flip_flop) THEN
      
      IF ( iswm_oce /= 1 ) THEN
        z_e = inverse_primal_flip_flop(patch_2D,patch_3d,op_coeffs,&
          & ocean_state%p_diag%veloc_adv_horz, ocean_state%p_diag%h_e)
      ELSE
        z_e = inverse_primal_flip_flop(patch_2D,patch_3d, op_coeffs, &
          & ocean_state%p_diag%veloc_adv_horz, ocean_state%p_diag%thick_e)
      ENDIF
      
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src = 5  ! output print level (1-5, fix)
      CALL dbg_print('bef.dual-flip-fl: LaPlaHorz',ocean_state%p_diag%laplacian_horz,str_module,idt_src)
      !---------------------------------------------------------------------
      
      IF(l_staggered_timestep)THEN
        
        DO jb = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, jb, edge_start_idx, edge_end_idx)
          DO jk = 1, n_zlev
            
            ocean_state%p_aux%g_n(edge_start_idx:edge_end_idx, jk, jb) = &
              & - z_e(edge_start_idx:edge_end_idx,jk,jb)  &
              & - ocean_state%p_diag%veloc_adv_vert(edge_start_idx:edge_end_idx,jk,jb)  &
              & + ocean_state%p_diag%laplacian_horz(edge_start_idx:edge_end_idx,jk,jb)  &
              & + ocean_state%p_diag%laplacian_vert(edge_start_idx:edge_end_idx,jk,jb)
          END DO
        END DO
        
      ELSEIF(.NOT.l_staggered_timestep)THEN
        DO jb = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, jb, edge_start_idx, edge_end_idx)
          DO jk = 1, n_zlev
            ocean_state%p_aux%g_n(edge_start_idx:edge_end_idx, jk, jb) = &
              & -ocean_state%p_diag%press_grad(edge_start_idx:edge_end_idx, jk, jb)       &
              & - z_e(edge_start_idx:edge_end_idx, jk, jb)  &
              & - ocean_state%p_diag%veloc_adv_vert(edge_start_idx:edge_end_idx, jk, jb)  &
              & + ocean_state%p_diag%laplacian_horz(edge_start_idx:edge_end_idx, jk, jb)  &
              & + ocean_state%p_diag%laplacian_vert(edge_start_idx:edge_end_idx, jk, jb)
          END DO
        END DO
      ENDIF
      
      CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_aux%g_n)
      
    ELSE ! IF(.NOT.(l_inverse_flip_flop))THEN
      
      IF(l_staggered_timestep)THEN
        DO jb = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, jb, edge_start_idx, edge_end_idx)
          DO jk = 1, n_zlev
            ocean_state%p_aux%g_n(edge_start_idx:edge_end_idx, jk, jb) =&!-ocean_state%p_diag%press_grad(:,jk,:)      &
              & - ocean_state%p_diag%veloc_adv_horz(edge_start_idx:edge_end_idx, jk, jb)  &
              & - ocean_state%p_diag%veloc_adv_vert(edge_start_idx:edge_end_idx, jk, jb)  &
              & + ocean_state%p_diag%laplacian_horz(edge_start_idx:edge_end_idx, jk, jb)  &
              & + ocean_state%p_diag%laplacian_vert(edge_start_idx:edge_end_idx, jk, jb)
          END DO
        END DO
      ELSE ! IF(.NOT.l_staggered_timestep)THEN
        DO jb = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, jb, edge_start_idx, edge_end_idx)
          DO jk = 1, n_zlev
            ocean_state%p_aux%g_n(edge_start_idx:edge_end_idx, jk, jb) = &
              & - ocean_state%p_diag%press_grad(edge_start_idx:edge_end_idx, jk, jb)      &
              & - ocean_state%p_diag%veloc_adv_horz(edge_start_idx:edge_end_idx, jk, jb)  &
              & - ocean_state%p_diag%veloc_adv_vert(edge_start_idx:edge_end_idx, jk, jb)  &
              & + ocean_state%p_diag%laplacian_horz(edge_start_idx:edge_end_idx, jk, jb)  &
              & + ocean_state%p_diag%laplacian_vert(edge_start_idx:edge_end_idx, jk, jb)
          END DO
        END DO
      ENDIF
      
    ENDIF!(L_INVERSE_FLIP_FLOP)
    
    IF(l_initial_timestep)THEN
      ocean_state%p_aux%g_nimd(1:nproma,1:n_zlev,1:patch_2D%nblks_e) = &
        & ocean_state%p_aux%g_n(1:nproma,1:n_zlev,1:patch_2D%nblks_e)
    ELSE
      ocean_state%p_aux%g_nimd(1:nproma,1:n_zlev,1:patch_2D%nblks_e) &
        & = (1.5_wp+ab_const)* ocean_state%p_aux%g_n(1:nproma,1:n_zlev,1:patch_2D%nblks_e)&
        & - (0.5_wp+ab_const)* ocean_state%p_aux%g_nm1(1:nproma,1:n_zlev,1:patch_2D%nblks_e)
    ENDIF
    
    IF ( iswm_oce /= 1) THEN
      
      IF(.NOT.l_rigid_lid)THEN
        
        IF(l_staggered_timestep)THEN
          DO jb = all_edges%start_block, all_edges%end_block
            CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
            DO jk = 1, n_zlev
              DO je = edge_start_idx, edge_end_idx
                
                IF(patch_3d%p_patch_1d(1)%dolic_e(je,jb)>=min_dolic)THEN
                  ocean_state%p_diag%vn_pred(je,jk,jb) = ocean_state%p_prog(nold(1))%vn(je,jk,jb)    &
                    & + dtime*(ocean_state%p_aux%g_nimd(je,jk,jb)     &
                    & -ocean_state%p_diag%press_grad(je,jk,jb)        &
                    & - z_gradh_e(je ,jb))
                ENDIF

              END DO
            END DO
          END DO
          
        ELSE   ! IF(.NOT.l_staggered_timestep)THEN
          
          !---------DEBUG DIAGNOSTICS-------------------------------------------
          idt_src = 5  ! output print level (1-5, fix)
          CALL dbg_print('Bef fin term: vn_pred'      ,ocean_state%p_diag%vn_pred       ,str_module,idt_src, &
             & in_subset=owned_edges)
          CALL dbg_print('Bef fin term: g_nimd'       ,ocean_state%p_aux%g_nimd         ,str_module,idt_src, &
             & in_subset=owned_edges)
          !---------------------------------------------------------------------
          
          DO jb = all_edges%start_block, all_edges%end_block
            CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
            DO jk = 1, n_zlev
              DO je = edge_start_idx, edge_end_idx
                
                IF(patch_3d%p_patch_1d(1)%dolic_e(je,jb)>=min_dolic)THEN
                  ocean_state%p_diag%vn_pred(je,jk,jb) = ocean_state%p_prog(nold(1))%vn(je,jk,jb)  &
                    & + dtime*(ocean_state%p_aux%g_nimd(je,jk,jb) &
                    & - z_gradh_e(je, jb))
                ENDIF

              END DO
            END DO
          END DO
        ENDIF!Staggered
        
        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src = 5  ! output print level (1-5, fix)
        CALL dbg_print('Aft fin term: vn_pred'      ,ocean_state%p_diag%vn_pred       ,str_module,idt_src, in_subset=owned_edges)
        !---------------------------------------------------------------------
        
      ELSEIF(l_rigid_lid)THEN
        
        IF(l_staggered_timestep)THEN
          
          DO jb = all_edges%start_block, all_edges%end_block
            CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
            DO jk = 1, n_zlev
              DO je = edge_start_idx, edge_end_idx
                
                IF(patch_3d%p_patch_1d(1)%dolic_e(je,jb)>=min_dolic)THEN
                  ocean_state%p_diag%vn_pred(je,jk,jb) = ocean_state%p_prog(nold(1))%vn(je,jk,jb)  &
                    & + dtime*(ocean_state%p_aux%g_nimd(je,jk,jb) &
                    & - ocean_state%p_diag%press_grad(je,jk,jb))
                ENDIF
                
              END DO
            END DO
          END DO
          
        ELSEIF(.NOT.l_staggered_timestep)THEN
          
          DO jb = all_edges%start_block, all_edges%end_block
            CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
            DO jk = 1, n_zlev
              DO je = edge_start_idx, edge_end_idx
                
                IF(patch_3d%p_patch_1d(1)%dolic_e(je,jb)>=min_dolic)THEN
                  ocean_state%p_diag%vn_pred(je,jk,jb) = ocean_state%p_prog(nold(1))%vn(je,jk,jb)&
                    & + dtime*ocean_state%p_aux%g_nimd(je,jk,jb)
                ENDIF
                
              END DO
            END DO
          END DO
          
        ENDIF!Staggered
      ENDIF!Rigid lid
      
      !IF surface forcing applied as top boundary condition to vertical diffusion
      !The surface forcing is applied as volume forcing at rhs,
      !i.e. if it part of explicit term in momentum and tracer eqs.
      !in this case, top boundary ondition of vertical Laplacians are homogeneous.
      !Below is the code that adds surface forcing to explicit term of momentum eq.
      IF(expl_vertical_velocity_diff==1)THEN
        
         DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
          DO je = edge_start_idx, edge_end_idx
            
            IF(patch_3d%p_patch_1d(1)%dolic_e(je,jb)>=min_dolic)THEN
              ocean_state%p_diag%vn_pred(je,1,jb) =  ocean_state%p_diag%vn_pred(je,1,jb)      &
                & + dtime*ocean_state%p_aux%bc_top_vn(je,jb)&
                & /patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(je,1,jb) ! Change to prism_thick_e ?
              
              dolic_e=patch_3d%p_patch_1d(1)%dolic_e(je,jb)
              ocean_state%p_diag%vn_pred(je,dolic_e,jb)    &
                & = ocean_state%p_diag%vn_pred(je,dolic_e,jb)&
                & - dtime*ocean_state%p_aux%bc_bot_vn(je,jb)               &
                & /patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(je,dolic_e,jb)
            ENDIF
            
          END DO
        END DO
      ELSE
      
        CALL finish(TRIM('mo_oce_ab_timestepping_mimetic:calculate_explicit_term_ab'), &
              &  'explicit vertical diffusion no longer supported')
      ENDIF
      
      !In the SW-case the external forcing is applied as volume force.
      !This force is stored in data type top-boundary-condition.
    ELSEIF ( iswm_oce == 1)THEN! .AND. iforc_oce==11) THEN
      
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
        DO jk = 1, n_zlev
          DO je = edge_start_idx, edge_end_idx
            !IF(patch_3d%lsm_c(je,jk,jb) <= sea_boundary)THEN
            ocean_state%p_diag%vn_pred(je,jk,jb) = (ocean_state%p_prog(nold(1))%vn(je,jk,jb)        &
              & + dtime*(ocean_state%p_aux%g_nimd(je,jk,jb)        &
              & - z_gradh_e(je, jb)))
          END DO
        END DO
      END DO
      
      !In case of Shallow-water with forcing and or damping
      IF ( iforc_oce/=10) THEN
        
        DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
          DO jk = 1, n_zlev
            DO je = edge_start_idx, edge_end_idx
              ocean_state%p_diag%vn_pred(je,jk,jb) = (ocean_state%p_diag%vn_pred(je,jk,jb) &
                & + ocean_state%p_aux%bc_top_vn(je,jb)    &
                & - ocean_state%p_aux%bc_bot_vn(je,jb))
              !ENDIF
            END DO
          END DO
        END DO
      ENDIF
    ENDIF

    CALL dbg_print('Bef veloc_diff_vert: vn_pred', ocean_state%p_diag%vn_pred, str_module,idt_src, &
      & in_subset=owned_edges)
    !In 3D case and if implicit vertical velocity diffusion is chosen
    IF(iswm_oce /= 1.AND.expl_vertical_velocity_diff==1)THEN
      
      !Surface forcing is implemented as volume forcing in top layer.
      !In this case homogeneous boundary conditions for vertical Laplacian
      
      CALL velocity_diffusion_vertical_implicit( patch_3d,             &
        & ocean_state%p_diag%vn_pred,      &
        & p_phys_param%a_veloc_v,   &
        & op_coeffs) !,               &
!        & ocean_state%p_diag%vn_impl_vert_diff)
!      IF(l_rigid_lid)THEN
!        ocean_state%p_diag%vn_pred(1:nproma,1:n_zlev,1:patch_2D%nblks_e) &
!          & = ocean_state%p_diag%vn_impl_vert_diff(1:nproma,1:n_zlev,1:patch_2D%nblks_e)
!      ENDIF
      
    ENDIF
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('vn(nold)'                  ,ocean_state%p_prog(nold(1))%vn       ,str_module,idt_src, &
      & in_subset=owned_edges)
    CALL dbg_print('VelocAdvHorizontal'        ,ocean_state%p_diag%veloc_adv_horz    ,str_module,idt_src, &
      & in_subset=owned_edges)
    CALL dbg_print('VelocLaPlac horizontal'    ,ocean_state%p_diag%laplacian_horz    ,str_module,idt_src, &
      & in_subset=owned_edges)
    
    IF (expl_vertical_velocity_diff == 1 .AND. iswm_oce /= 1) THEN
      CALL dbg_print('ImplVelocDiff vertical'    ,ocean_state%p_diag%vn_pred ,str_module,idt_src, &
        & in_subset=owned_edges )
    ELSE
      CALL dbg_print('VelocLaPlac vertical'      ,ocean_state%p_diag%laplacian_vert    ,str_module,idt_src, &
        & in_subset=owned_edges )
    ENDIF
    IF (l_inverse_flip_flop) &
      & CALL dbg_print('dual-flip-flop Adv. horz'  ,z_e                           ,str_module,idt_src, &
        & in_subset=owned_edges )
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('G_n+1/2 - g_nimd'          ,ocean_state%p_aux%g_nimd             ,str_module,idt_src, &
        & in_subset=owned_edges )
    CALL dbg_print('G_n'                       ,ocean_state%p_aux%g_n                ,str_module,idt_src, &
        & in_subset=owned_edges )
    CALL dbg_print('G_n-1'                     ,ocean_state%p_aux%g_nm1              ,str_module,idt_src, &
        & in_subset=owned_edges )
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('vn_pred'                   ,ocean_state%p_diag%vn_pred           ,str_module,idt_src, &
        & in_subset=owned_edges )
    !---------------------------------------------------------------------
    
    !-------------------------------------------
    !CALL message (TRIM(routine), 'end')
    
  END SUBROUTINE calculate_explicit_term_ab
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !!  Calculation of right-hand side of elliptic surface equation.
  !!  This is used in semi implicit timelevel stepping.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize_Used>
  SUBROUTINE fill_rhs4surface_eq_ab( patch_3d, ocean_state, p_sfc_flx, op_coeffs)
    !
    ! Patch on which computation is performed
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_sfc_flx), INTENT(in)          :: p_sfc_flx
    TYPE(t_operator_coeff)               :: op_coeffs
    !
    !  local variables
    !
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: edge_start_idx, edge_end_idx
    INTEGER :: jc, jb, jk, je
    INTEGER :: i_dolic_c,i_dolic_e
    REAL(wp) :: gdt2    !, delta_z
    REAL(wp) :: z_e2d(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: div_z_depth_int_c(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_z_c(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_vn_ab(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain, all_edges, owned_edges
    TYPE(t_patch), POINTER :: patch_2D
    !REAL(wp) :: thick
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !       & routine = ('mo_oce_ab_timestepping_mimetic:fill_rhs4surface_eq_ab')
    !-------------------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    
    patch_2D        => patch_3d%p_patch_2d(1)
    all_cells       => patch_3d%p_patch_2d(1)%cells%ALL
    cells_in_domain => patch_3d%p_patch_2d(1)%cells%in_domain
    all_edges       => patch_3d%p_patch_2d(1)%edges%ALL
    owned_edges     => patch_3d%p_patch_2d(1)%edges%owned
    
    gdt2 = grav*dtime*dtime
    
    z_vn_ab(1:nproma,1:n_zlev,1:patch_2D%nblks_e)  = 0.0_wp
    z_e2d(1:nproma,1:patch_2D%nblks_e)             = 0.0_wp
    z_e(1:nproma,1:n_zlev,1:patch_2D%nblks_e)      = 0.0_wp
    div_z_depth_int_c(1:nproma,1:patch_2D%alloc_cell_blocks) = 0.0_wp
    div_z_c(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)  = 0.0_wp
    
    
    ! LL: this should not be required
    CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_diag%vn_pred)
    CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_prog(nold(1))%vn)
    ! CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_diag%vn_impl_vert_diff)
    
    IF(iswm_oce == 1)THEN
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
        DO je = edge_start_idx, edge_end_idx
          DO jk=1,n_zlev
            z_vn_ab(je,jk,jb)=ab_gam*ocean_state%p_diag%vn_pred(je,jk,jb)&
              & + (1.0_wp -ab_gam)* ocean_state%p_prog(nold(1))%vn(je,jk,jb)
          END DO
        ENDDO
      END DO
    ELSEIF(iswm_oce /= 1)THEN
      
      IF(expl_vertical_velocity_diff==1)THEN
        DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
          DO je = edge_start_idx, edge_end_idx
            i_dolic_e =  patch_3d%p_patch_1d(1)%dolic_e(je,jb)
            DO jk=1,i_dolic_e
              z_vn_ab(je,jk,jb)=ab_gam * ocean_state%p_diag%vn_pred(je,jk,jb) &
                & + (1.0_wp -ab_gam) * ocean_state%p_prog(nold(1))%vn(je,jk,jb)
            END DO
          ENDDO
        END DO
        
      ELSEIF(expl_vertical_velocity_diff==0)THEN
        
        DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
          DO je = edge_start_idx, edge_end_idx
            i_dolic_e =  patch_3d%p_patch_1d(1)%dolic_e(je,jb)
            DO jk=1,i_dolic_e
              z_vn_ab(je,jk,jb)=ab_gam*ocean_state%p_diag%vn_pred(je,jk,jb)&
                & + (1.0_wp -ab_gam)* ocean_state%p_prog(nold(1))%vn(je,jk,jb)
            END DO
          ENDDO
        END DO
        
      ENDIF
    ENDIF
    
    !
    ! calculate depth-integrated velocity z_e
    !  - edge-based and cell-based
    !  - 3d and 2d (surface)
    
    ! !-------------------------------------------------------------------------------
    IF (l_edge_based) THEN
      ! !-------------------------------------------------------------------------------
      
      IF( iswm_oce /= 1 ) THEN !the 3D case
        
        DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
          DO je = edge_start_idx, edge_end_idx
            i_dolic_e =  patch_3d%p_patch_1d(1)%dolic_e(je,jb)
            DO jk=1,i_dolic_e
              z_e(je,jk,jb)= z_vn_ab(je,jk,jb)*patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,jb)
            END DO
          ENDDO
        END DO
        
      ELSEIF( iswm_oce == 1 ) THEN
        
        DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
          DO je = edge_start_idx, edge_end_idx
            i_dolic_e =  patch_3d%p_patch_1d(1)%dolic_e(je,jb)
            DO jk=1,i_dolic_e
              z_e(je,jk,jb)= z_vn_ab(je,jk,jb)*ocean_state%p_diag%thick_e(je,jb)
            END DO
          ENDDO
        END DO
        
      ENDIF

      ! !-------------------------------------------------------------------------------
    ELSEIF(.NOT. l_edge_based)THEN!NOT EDGE-BASED
      ! !-------------------------------------------------------------------------------
      
      IF( iswm_oce /= 1 ) THEN !the 3D case
        CALL map_edges2edges_viacell_3d_const_z( patch_3d, z_vn_ab, op_coeffs, z_e )
      ELSEIF( iswm_oce == 1 ) THEN
        !    CALL map_edges2edges_viacell_3D( patch_3d,    &
        !                                    & z_vn_ab(:,1,:),&
        !                                    & op_coeffs,    &
        !                                    & z_e(:,1,:),    &
        !                                    & ocean_state%p_diag%thick_c, level=1)
        CALL map_edges2edges_viacell_3d_const_z( patch_3d, z_vn_ab(:,1,:), op_coeffs, z_e(:,1,:) )
        
      ENDIF!( iswm_oce == 1 )
      
    ENDIF!EDGE-BASED
    
    CALL div_oce_3d( z_e, patch_2D,op_coeffs%div_coeff, div_z_c, subset_range=cells_in_domain )
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        i_dolic_c = patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
        DO jk=1,i_dolic_c
          IF(patch_3d%lsm_c(jc,jk,jb) <= sea_boundary)THEN
            div_z_depth_int_c(jc,jb)=div_z_depth_int_c(jc,jb) +div_z_c(jc,jk,jb)
          ENDIF
        END DO
      END DO
    END DO
    
    !-------------------------------------------------------------------------
    ! Apply net surface freshwater flux to elevation - incorrect?
    !IF(forcing_enable_freshwater)THEN
    !  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    !    CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
    !    DO jc = i_startidx_c, i_endidx_c
    !      IF(patch_3d%lsm_c(jc,1,jb) <= sea_boundary)THEN
    !        ocean_state%p_aux%p_rhs_sfc_eq(jc,jb) = ((ocean_state%p_prog(nold(1))%h(jc,jb)     &
    !                                       & - dtime*(div_z_depth_int_c(jc,jb) + &
    !                                       &          p_sfc_flx%forc_fwfx(jc,jb)) )/gdt2)
    !      ENDIF
    !    ENDDO
    !  END DO
    
    !ELSEIF(.NOT.forcing_enable_freshwater)THEN
    
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(patch_3d%lsm_c(jc,1,jb) <= sea_boundary)THEN
          ocean_state%p_aux%p_rhs_sfc_eq(jc,jb) = ((ocean_state%p_prog(nold(1))%h(jc,jb)&
            & - dtime*div_z_depth_int_c(jc,jb))/gdt2)
        ENDIF
      ENDDO
    END DO
    !ENDIF
    
    CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_aux%p_rhs_sfc_eq )
    
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('RHS thick_e'               ,patch_3d%p_patch_1d(1)%prism_thick_e, str_module,idt_src, &
        & in_subset=owned_edges )
    CALL dbg_print('RHS thick_c'               ,patch_3d%p_patch_1d(1)%prism_thick_c, str_module,idt_src, &
        & in_subset=patch_3d%p_patch_2d(1)%cells%owned)        
    CALL dbg_print('RHS z_vn_ab'               ,z_vn_ab                  ,str_module,idt_src, &
        & in_subset=owned_edges )
    CALL dbg_print('RHS z_e'                   ,z_e                      ,str_module,idt_src, &
        & in_subset=owned_edges )
    CALL dbg_print('RHS div_z_depth_int_c'     ,div_z_depth_int_c        ,str_module,idt_src, &
      in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    CALL dbg_print('RHS div_z_c'     ,div_z_c        ,str_module,idt_src, &
      in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('RHS final'                 ,ocean_state%p_aux%p_rhs_sfc_eq  ,str_module,idt_src, &
      in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    !---------------------------------------------------------------------
    
  END SUBROUTINE fill_rhs4surface_eq_ab
  !-------------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------------
!<Optimize_Used>
  SUBROUTINE init_ho_lhs_fields_mimetic(patch_3d)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    
    TYPE(t_patch), POINTER :: patch    
    INTEGER :: return_status
    
    patch         => patch_3d%p_patch_2d(1)
    
    ALLOCATE(lhs_result(nproma,patch%alloc_cell_blocks), &
      & lhs_z_grad_h(nproma,patch%nblks_e),     &
      & lhs_z_e     (nproma,patch%nblks_e),     &
!      & lhs_z_e_top (nproma,patch%nblks_e),     &
!      & lhs_z_grad_h_cc(nproma,patch%alloc_cell_blocks),  &
      & stat = return_status)
    IF (return_status > 0) &
      & CALL finish("mo_oce_ab_timestepping_mimetic:init_ho_lhs_fields", "Allocation failed")

    ALLOCATE(lhs_result_sp(nproma,patch%alloc_cell_blocks), &
      & lhs_z_grad_h_sp(nproma,patch%nblks_e),     &
      & lhs_z_e_sp     (nproma,patch%nblks_e),     &
      & stat = return_status)
    IF (return_status > 0) &
      & CALL finish("mo_oce_ab_timestepping_mimetic:init_ho_lhs_fields", "sp Allocation failed")
    
    ! these are arrays used by the lhs routine
    lhs_result(:,:)   = 0.0_wp
    lhs_z_grad_h(:,:) = 0.0_wp
    lhs_z_e     (:,:) = 0.0_wp
!    lhs_z_e_top (:,:) = 0.0_wp

    lhs_result_sp(:,:)   = 0.0_wp
    lhs_z_grad_h_sp(:,:) = 0.0_wp
    lhs_z_e_sp     (:,:) = 0.0_wp
    

!    lhs_z_grad_h_cc(:,:)%x(1) = 0.0_wp
!    lhs_z_grad_h_cc(:,:)%x(2) = 0.0_wp
!    lhs_z_grad_h_cc(:,:)%x(3) = 0.0_wp

  END SUBROUTINE init_ho_lhs_fields_mimetic
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
  !!  The result is NOT synced. Should be done in the calling method if required
  !-------------------------------------------------------------------------
!<Optimize_Used>
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
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc, jb, je
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D     ! patch_2D on which computation is performed
    !-----------------------------------------------------------------------
    IF (ltimer) CALL timer_start(timer_lhs)
    !-----------------------------------------------------------------------
    patch_2D           => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain
    
    
    gdt2_inv       = 1.0_wp / (grav*(dtime)**2)
    gam_times_beta = ab_gam * ab_beta
    
    lhs   (1:nproma,cells_in_domain%end_block:patch_2D%alloc_cell_blocks)  = 0.0_wp
    
    CALL sync_patch_array(sync_c, patch_2D, x(1:nproma,1:patch_2D%cells%all%end_block) )
    
    !Step 1) Calculate gradient of iterated height.
    CALL grad_fd_norm_oce_2d_3d( x, &
      & patch_2D,                       &
      & op_coeffs%grad_coeff(:,1,:),&
      & lhs_z_grad_h(:,:))
    
    ! the result lhs_z_grad_h is computed on in_domain edges
    ! CALL sync_patch_array(sync_e, patch_2D, lhs_z_grad_h(:,:) )
    
    
    !---------------------------------------
    IF(l_edge_based)THEN
      
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
        DO je = i_startidx, i_endidx
          lhs_z_e(je,jb) = lhs_z_grad_h(je,jb) * thickness_e(je,jb)
          
!          IF(patch_3d%lsm_e(je,1,jb) > sea_boundary)THEN
!            IF (lhs_z_e(je,jb) /= 0.0_wp) &
!              CALL finish("lhs_surface_height_ab_mim", "lhs_z_e(je,jb) /= 0 on land")
!          ENDIF

        END DO
      END DO
      ! no need to sync since we will compute only cells in domain
      !CALL sync_patch_array(SYNC_E, patch_2D, lhs_z_e(:,:) )
      
    ELSE  !IF(.NOT.l_edge_based)THEN
      
      ! the map_edges2edges_viacell_3d_const_z should be changes to calculate only
      ! on in_domai_edges. Still, edge valkues need to be synced
      CALL sync_patch_array(sync_e, patch_2D, lhs_z_grad_h(:,:) )
      
      
      IF( iswm_oce /= 1 ) THEN
        
        CALL map_edges2edges_viacell_3d_const_z( patch_3d, lhs_z_grad_h(:,:), op_coeffs, lhs_z_e(:,:))
        
      ELSEIF( iswm_oce == 1 ) THEN
        
        !CALL map_edges2edges_viacell_3D( patch_3d, lhs_z_grad_h, op_coeffs, lhs_z_e,thickness_c, level=top)
        CALL map_edges2edges_viacell_3d_const_z( patch_3d, lhs_z_grad_h(:,:), op_coeffs, lhs_z_e(:,:))
      ENDIF!( iswm_oce == 1 )
      
    ENDIF ! l_edge_based
    !---------------------------------------
    
    !Step 3) Calculate divergence
    ! store the div in lhs for reducing memory and improving performance
    CALL div_oce_3d( lhs_z_e, patch_2D, op_coeffs%div_coeff, lhs, &
      & level=top, subset_range=cells_in_domain  )
    
    !Step 4) Finalize LHS calculations
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
      DO jc = i_startidx, i_endidx
        
        !lhs(jc,jb) =(x(jc,jb) - gdt2 * ab_gam * ab_beta * lhs_div_z_c(jc,jb)) / gdt2 !rho_sfc(jc,jb)*rho_inv
        lhs(jc,jb) = x(jc,jb) * gdt2_inv - gam_times_beta * lhs(jc,jb)
        
      END DO
    END DO
    !---------------------------------------

    IF (debug_check_level > 20) THEN
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
        DO jc = i_startidx, i_endidx
          IF(patch_3d%lsm_c(jc,1,jb) > sea_boundary) THEN
            IF (lhs(jc,jb) /= 0.0_wp) &
              & CALL finish("lhs_surface_height_ab_mim", "lhs(jc,jb) /= 0 on land")
          ENDIF
        END DO
      END DO
    ENDIF

    
    IF (ltimer) CALL timer_stop(timer_lhs)
    
  END FUNCTION lhs_surface_height_ab_mim
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------
  ! as in lhs_surface_height_ab_mim in single precision
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
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc, jb, je
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D     ! patch_2D on which computation is performed
    REAL(wp) :: x_sync(SIZE(x,1), SIZE(x,2))    ! used to syn x, since we cannot synd single precision at the moment
    !-----------------------------------------------------------------------
    IF (ltimer) CALL timer_start(timer_lhs_sp)
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

    !Step 1) Calculate gradient of iterated height.
    CALL grad_fd_norm_oce_2d_3d_sp( x,   &
      & patch_2D,                     &
      & solverCoeffs%grad_coeff(:,:), &
      & lhs_z_grad_h_sp(:,:))

    !---------------------------------------
    IF(l_edge_based)THEN

      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
        DO je = i_startidx, i_endidx
          lhs_z_e_sp(je,jb) = lhs_z_grad_h_sp(je,jb) * solverCoeffs%edge_thickness(je,jb)
        END DO
      END DO

    ELSE  !IF(.NOT.l_edge_based)THEN

      ! the map_edges2edges_viacell_3d_const_z should be changes to calculate only
      ! on in_domain_edges. Still, edge valkues need to be synced
      lhs_z_grad_h = REAL(lhs_z_grad_h_sp, wp)
      CALL sync_patch_array(sync_e, patch_2D, lhs_z_grad_h(:,:) )
      lhs_z_grad_h_sp = REAL(lhs_z_grad_h, sp)

      CALL map_edges2edges_viacell_2d_constZ_sp( patch_3d, lhs_z_grad_h_sp(:,:), solverCoeffs, lhs_z_e_sp(:,:))

    ENDIF ! l_edge_based
    !---------------------------------------

    !Step 3) Calculate divergence
    ! store the div in lhs for reducing memory and improving performance
    CALL div_oce_2d_sp( lhs_z_e_sp, patch_2D, solverCoeffs%div_coeff, lhs, &
      & subset_range=cells_in_domain  )

    !Step 4) Finalize LHS calculations
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
      DO jc = i_startidx, i_endidx

        lhs(jc,jb) = x(jc,jb) * gdt2_inv - gam_times_beta * lhs(jc,jb)

      END DO
    END DO
    !---------------------------------------

    IF (ltimer) CALL timer_stop(timer_lhs_sp)

  END FUNCTION lhs_surface_height_ab_mim_sp
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computation of new velocity in Adams-Bashforth timestepping.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize_Used>
  SUBROUTINE calc_normal_velocity_ab_mimetic(patch_3d,ocean_state, op_coeffs, solverCoeff_sp, p_ext_data)
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp
    TYPE(t_external_data), TARGET        :: p_ext_data
    !
    !  local variables
    INTEGER :: edge_start_idx, edge_end_idx
    INTEGER :: je, jk, jb
    REAL(wp) :: gdt
    REAL(wp) :: z_grad_h(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    !REAL(wp) :: z_h_new_rho(nproma, patch_3d%p_patch_2d(n_dom)%nblks_c)
    TYPE(t_subset_range), POINTER :: edges_in_domain, owned_cells, owned_edges
    CHARACTER(LEN=*), PARAMETER ::     &
      & method_name='mo_oce_ab_timestepping_mimetic: calc_normal_velocity_ab_mimetic'
    TYPE(t_patch), POINTER :: patch
    !----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    !-----------------------------------------------------------------------
    patch           => patch_3d%p_patch_2d(1)
    edges_in_domain => patch%edges%in_domain
    owned_cells     => patch%cells%owned
    owned_edges     => patch%edges%owned
    !----------------------------------------------------------------------
    
    z_grad_h(1:nproma,1:patch%nblks_e) = 0.0_wp
    
    gdt=grav*dtime
    
    ! Step 1) Compute normal derivative of new surface height
    IF(.NOT.l_rigid_lid.OR.iswm_oce == 1) THEN
      CALL grad_fd_norm_oce_2d_3d( ocean_state%p_prog(nnew(1))%h, &
        & patch,                                                  &
        & op_coeffs%grad_coeff(:,1,:),                           &
        & z_grad_h(:,:))
    ENDIF
    
    
    ! Step 2) Calculate the new velocity from the predicted one and the new surface height
    IF (iswm_oce == 1) THEN ! shallow water case
      
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, edge_start_idx, edge_end_idx)
        DO jk = 1, n_zlev
          DO je = edge_start_idx, edge_end_idx
            IF(patch_3d%lsm_e(je,jk,jb) <= sea_boundary)THEN
              ocean_state%p_prog(nnew(1))%vn(je,jk,jb) = (ocean_state%p_diag%vn_pred(je,jk,jb) &
                & - gdt*ab_beta*z_grad_h(je,jb))
            ENDIF
          END DO
        END DO
      END DO
      
    ELSE !real 3d case
      
      IF (.NOT.l_rigid_lid) THEN
        
        DO jb = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, jb, edge_start_idx, edge_end_idx)
          DO je = edge_start_idx, edge_end_idx
            DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
                ocean_state%p_prog(nnew(1))%vn(je,jk,jb) = (ocean_state%p_diag%vn_pred(je,jk,jb) &
                  & - gdt*ab_beta*z_grad_h(je,jb))
            END DO
          END DO
        END DO
        
      ELSE
        CALL finish(method_name,"l_rigid_lid case has a bug")
        !ocean_state%p_prog(nnew(1))%vn(:,:,:) = ocean_state%p_diag%vn_pred*v_base%wet_e(:,:,:)
      ENDIF
    ENDIF
    
    
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, edge_start_idx, edge_end_idx)
      DO je = edge_start_idx, edge_end_idx
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
          ocean_state%p_diag%vn_time_weighted(je,jk,jb)      &
            & = ab_gam*ocean_state%p_prog(nnew(1))%vn(je,jk,jb) &
            & + (1.0_wp -ab_gam)*ocean_state%p_prog(nold(1))%vn(je,jk,jb)
        END DO
      END DO
    END DO
    
    CALL sync_patch_array(sync_e, patch, ocean_state%p_prog(nnew(1))%vn)
    CALL sync_patch_array(sync_e, patch, ocean_state%p_diag%vn_time_weighted)
    
    ! slo: z_grad_h not out of sync, but sync error with global_max in dbg_print?
    !CALL sync_patch_array(SYNC_E, patch, z_grad_h)
    
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('NormVel: vn_old'            ,ocean_state%p_prog(nold(1))%vn     ,str_module,idt_src, in_subset=owned_edges )
    CALL dbg_print('NormVel: vn_pred'           ,ocean_state%p_diag%vn_pred         ,str_module,idt_src, in_subset=owned_edges)
    IF (.NOT.l_rigid_lid) THEN
      CALL dbg_print('NormVel: grad h-new'      ,z_grad_h                    ,str_module,idt_src, in_subset=owned_edges)
    END IF
    CALL dbg_print('NormVel: vn_time_weighted'  ,ocean_state%p_diag%vn_time_weighted,str_module,idt_src, in_subset=owned_edges)
    CALL dbg_print('NormVel: vn_change'         ,ocean_state%p_prog(nnew(1))%vn - &
      & ocean_state%p_prog(nold(1))%vn     ,str_module,idt_src, in_subset=owned_edges)
    idt_src=2  ! outputm print level (1-5, fix)
    CALL dbg_print('NormVel: vn_new'            ,ocean_state%p_prog(nnew(1))%vn     ,str_module,idt_src, in_subset=owned_edges)
    !---------------------------------------------------------------------
    
    ! Update of scalar product quantities
    IF(l_staggered_timestep)THEN
      !CALL height_related_quantities(patch_3d, ocean_state, p_ext_data)
      CALL calculate_thickness(patch_3d, ocean_state, p_ext_data, op_coeffs, solverCoeff_sp)
      !   CALL calc_scalar_product_veloc_3D( patch,                &
      !                                    & ocean_state%p_prog(nnew(1))%vn,&
      !                                    & ocean_state%p_prog(nnew(1))%vn,&
      !                                    & ocean_state%p_diag,            &
      !                                    & op_coeffs)
      
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=4  ! output print level (1-5, fix)
      CALL dbg_print('NorVel: Staggered, kin'    ,ocean_state%p_diag%kin        ,str_module,idt_src, &
        in_subset = owned_cells)
      CALL dbg_print('NorVel: Staggered, ptp_vn' ,ocean_state%p_diag%ptp_vn     ,str_module,idt_src)
      !---------------------------------------------------------------------
    ENDIF
    !CALL message (TRIM(routine), 'end')
    
  END SUBROUTINE calc_normal_velocity_ab_mimetic
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computation of new vertical velocity using continuity equation
  
  !! Calculate diagnostic vertical velocity from horizontal velocity using the
  !! incommpressibility condition in the continuity equation.
  !! vertical velocity is integrated from bottom to top
  !! vertical velocity is negative for positive divergence
  !! of horizontal velocity
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn,   MPI-M (2006).
  !!  Modified by Stephan Lorenz, MPI-M (2010-06)
  !TODO review
!<Optimize_Used>
  SUBROUTINE calc_vert_velocity_mim_bottomup( patch_3d, ocean_state, op_coeffs )
    
    TYPE(t_patch_3d), TARGET, INTENT(in) :: patch_3d       ! patch on which computation is performed
    TYPE(t_hydro_ocean_state)            :: ocean_state
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs
    !
    !
    ! Local variables
    INTEGER :: jc, jk, jb, je
    INTEGER :: z_dolic
    INTEGER :: i_startidx, i_endidx
    REAL(wp) :: div_depth_int(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: z_vn_e(nproma,n_zlev,patch_3d%p_patch_2D(1)%nblks_e)
    REAL(wp) :: z_c(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_abort
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells
    REAL(wp) ::  minmaxmean(3)
    TYPE(t_patch), POINTER :: patch_2D
    ! TYPE(t_hydro_ocean_diag), POINTER   :: ocean_state%p_diag
    REAL(wp),  POINTER  :: vertical_velocity(:,:,:)

    CHARACTER(len=*), PARAMETER :: method_name='mo_oce_ab_timestepping_mimetic:alc_vert_velocity_mim_bottomup'
    
    !-----------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain  => patch_2D%cells%in_domain
    all_cells        => patch_2D%cells%all
    edges_in_domain  => patch_2D%edges%in_domain
    ! ocean_state%p_diag         = ocean_state%p_diag
    vertical_velocity => ocean_state%p_diag%w

    ! due to nag -nan compiler-option:
    ! vertical_velocity(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks) = 0.0_wp
    div_depth_int(1:nproma,1:patch_2D%alloc_cell_blocks)   = 0.0_wp
    z_c          (1:nproma,1:patch_2D%alloc_cell_blocks)   = 0.0_wp
    !------------------------------------------------------------------
    ! Step 1) Calculate divergence of horizontal velocity at all levels
    !------------------------------------------------------------------
    
    !-------------------------------------------------------------------------------
    IF( l_edge_based )THEN
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
        DO jk = 1, n_zlev
          DO je = i_startidx, i_endidx
            ocean_state%p_diag%mass_flx_e(je,jk,jb) = ocean_state%p_diag%vn_time_weighted(je,jk,jb)&
              & * patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,jb)
          END DO
        END DO
      END DO
     
      CALL div_oce_3d( ocean_state%p_diag%mass_flx_e,    &
        & patch_2D,                   &
        & op_coeffs%div_coeff,      &
        & ocean_state%p_diag%div_mass_flx_c,&
        & subset_range=cells_in_domain)

      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
        DO jc = i_startidx, i_endidx
          
          z_dolic = patch_3d%p_patch_1d(1)%dolic_c(jc,jb)         
          !use bottom boundary condition for vertical velocity at bottom
          !of prism
          vertical_velocity(jc,z_dolic+1,jb)=0.0_wp
          DO jk = z_dolic, 1, -1

            vertical_velocity(jc,jk,jb) = vertical_velocity(jc,jk+1,jb) - ocean_state%p_diag%div_mass_flx_c(jc,jk,jb)
          END DO
        END DO
      END DO
      
    !-------------------------------------------------------------------------------
    ELSE ! IF(.NOT.l_edge_based)THEN
      
      CALL map_edges2edges_viacell_3d_const_z( patch_3d, ocean_state%p_diag%vn_time_weighted, op_coeffs, &
        & ocean_state%p_diag%mass_flx_e)
      
      CALL div_oce_3d( ocean_state%p_diag%mass_flx_e,      &
        & patch_2D,op_coeffs%div_coeff,&
        & ocean_state%p_diag%div_mass_flx_c,  &
        & subset_range=cells_in_domain)
      
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
        DO jc = i_startidx, i_endidx
          
          z_dolic = patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          
          !use bottom boundary condition for vertical velocity at bottom of prism
          ! this should be awlays zero
          ! vertical_velocity(jc,z_dolic+1,jb)=0.0_wp

          DO jk = z_dolic, 1, -1
            vertical_velocity(jc,jk,jb) = vertical_velocity(jc,jk+1,jb) - ocean_state%p_diag%div_mass_flx_c(jc,jk,jb)
          END DO
        END DO
      END DO
      
      ! CALL sync_patch_array(sync_e,patch_2D,ocean_state%p_diag%mass_flx_e)
      ! CALL sync_patch_array(sync_c,patch_2D,ocean_state%p_diag%div_mass_flx_c)

    ENDIF  !  (l_EDGE_BASED)
    
    IF(l_rigid_lid)THEN
      vertical_velocity(:,1,:) = 0.0_wp
    ENDIF

    CALL sync_patch_array(sync_c,patch_2D,vertical_velocity)
    
    !-----------------------------------------------------
    IF (use_continuity_correction) THEN

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
        DO jc = i_startidx, i_endidx

          IF ( patch_3d%p_patch_1d(1)%dolic_c(jc,jb) > 0) &
            ocean_state%p_prog(nnew(1))%h(jc,jb) = ocean_state%p_prog(nold(1))%h(jc,jb) + &
              &  vertical_velocity(jc,1,jb) * dtime

        END DO
      END DO

      !---------------------------------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
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
    IF (debug_check_level > 3) THEN

      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
        DO jc = i_startidx, i_endidx

          z_c(jc,jb) = ((ocean_state%p_prog(nnew(1))%h(jc,jb)-ocean_state%p_prog(nold(1))%h(jc,jb))/dtime - &
            & vertical_velocity(jc,1,jb)) * patch_3d%wet_c(jc,1,jb)

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
        CALL message('mo_oce_ab_timestepping_mimetic:calc_vert_velocity_mim_bottomup', &
          & 'MISMATCH IN SURFACE EQUATION:')
        CALL message('mo_oce_ab_timestepping_mimetic:calc_vert_velocity_mim_bottomup', &
          & 'Elevation change does not match vertical velocity')
        WRITE(message_text,'(2(a,e20.12))') ' (h_new-h_old)/dtime - w = ', MAXVAL(ABS(z_c(:,:))), &
          & ' z_abort=', z_abort
        CALL message ('mo_oce_ab_timestepping_mimetic:calc_vert_velocity_mim_bottomup', message_text)
        CALL finish(TRIM('mo_oce_ab_timestepping_mimetic:calc_vert_velocity_mim_bottomup'), &
          & 'MISMATCH in surface equation')
      ENDIF
    ENDIF ! (debug_check_level > 3)
    !---------------------------------------------------------------------

    
  END SUBROUTINE calc_vert_velocity_mim_bottomup
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !!  The result is NOT synced. Should be done in the calling method if required
  FUNCTION inverse_primal_flip_flop(patch_2D, patch_3d, op_coeffs, rhs_e, h_e) result(inv_flip_flop_e)
    !
    TYPE(t_patch), TARGET                  :: patch_2D
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    TYPE(t_operator_coeff),INTENT(in)      :: op_coeffs
    REAL(wp)      :: rhs_e(:,:,:)!(nproma,n_zlev,patch_2D%nblks_e)
    REAL(wp)      :: h_e(:,:)  !(nproma,patch_2D%nblks_e)
    REAL(wp)      :: inv_flip_flop_e(SIZE(rhs_e,1),SIZE(rhs_e,2),SIZE(rhs_e,3))
    !
    !LOCAL VARIABLES
    INTEGER,PARAMETER :: nmax_iter= 800 ! maximum number of iterations
    REAL(wp) :: zimpl_coeff = 1.0_wp    !COEFF has to be set appropriately !!!!
    REAL(wp) :: zimpl_prime_coeff
    INTEGER :: n_iter                  ! number of iterations
    REAL(wp) :: tolerance               ! (relative or absolute) tolerance
    REAL(wp) :: z_residual(nmax_iter)
    LOGICAL :: lmax_iter               ! true if reached m iterations
    !LOGICAL  :: lverbose = .TRUE.
    !CHARACTER(len=MAX_CHAR_LENGTH) :: string
    REAL(wp) :: rhstemp(nproma,patch_2D%nblks_e)
    !REAL(wp), ALLOCATABLE :: inv_flip_flop_e2(:,:)!(nproma,patch_2D%nblks_e)
    REAL(wp) :: z_e(nproma,n_zlev,patch_2D%nblks_e)
    INTEGER :: jk
    !INTEGER :: i_startblk_e, i_endblk_e, edge_start_idx, edge_end_idx
    !-----------------------------------------------------------------------
    
    tolerance                = 1.0e-12_wp  ! solver_tolerance
    inv_flip_flop_e(:,:,:)   = 0.0_wp
    zimpl_prime_coeff = (1.0_wp-zimpl_coeff)
    
    rhstemp(:,:)          = 0.0_wp
    
    DO jk=1, n_zlev
      rhstemp(:,:) = rhs_e(:,jk,:)&
        & -zimpl_coeff*lhs_primal_flip_flop(inv_flip_flop_e(:,jk,:), patch_2D, patch_3d, op_coeffs,jk,zimpl_coeff, h_e)
      
      IF (MAXVAL (ABS (rhstemp (:,:))) <= tolerance) THEN
        inv_flip_flop_e(:,jk,:) = lhs_primal_flip_flop(inv_flip_flop_e(:,jk,:), patch_2D, patch_3d, op_coeffs,jk,zimpl_coeff, h_e)
        PRINT*, "Inv_flipflop gmres_oce_e2e solved by initial guess!",&
          & jk,MAXVAL(rhstemp(:,:)), MINVAL(rhstemp(:,:)),MAXVAL(rhs_e(:,jk,:)), MINVAL(rhs_e(:,jk,:))
      ELSE
        inv_flip_flop_e(:,jk,:)= 0.0_wp!rhs_e(:,jk,:)
        !write(*,*)'RHS', maxvaL(rhs_e(:,jk,:)),minvaL(rhs_e(:,jk,:))
        
        CALL gmres_oce_e2e( inv_flip_flop_e(:,jk,:), &  ! arg 1 of lhs. x input is the first guess.
          & lhs_primal_flip_flop,      &  ! function calculating l.h.s.
          & h_e,                       &  ! edge thickness for LHS
          & jk,                        &
          & patch_2D, patch_3d,       &  !arg 3 of lhs
          & zimpl_coeff,               &  !arg 4 of lhs
          & op_coeffs,                &
          & rhs_e(:,jk,:),             &  ! right hand side as input
          & tolerance,                 &  ! relative tolerance
          & .FALSE.,                   &  ! NOT absolute tolerance
          & nmax_iter,                 &  ! max. # of iterations to do
          & lmax_iter,                 &  ! out: .true. = not converged
          & n_iter,                    &  ! out: # of iterations done
          & z_residual)                  ! inout: the residual (array)
        
        rhstemp(:,:) = rhs_e(:,jk,:)-lhs_primal_flip_flop(inv_flip_flop_e(:,jk,:),patch_2D, patch_3d,op_coeffs,&
          & jk,zimpl_coeff,h_e)
        !WRITE(*,*)'max/min residual of inverse primal-flip-flop:',&
        !  &        jk, maxval(rhstemp),minval(rhstemp)
        idt_src=2  ! output print level (1-5, fix)
        CALL dbg_print('residual of inv_flip_flop'   ,rhstemp ,str_module,idt_src)
        !write(*,*)'sol', maxvaL(inv_flip_flop_e(:,jk,:)),minvaL(inv_flip_flop_e(:,jk,:))
        z_e(:,jk,:)=rhstemp(:,:)
        IF (MAXVAL (ABS (rhstemp (:,:))) >= tolerance) lmax_iter = .TRUE.
        idt_src=1
        IF (idbg_mxmn >= idt_src) THEN
          IF (lmax_iter) THEN
            WRITE (0, '(1x,a, I4.2, 1x, a,E8.2,1x, a,E8.2,1x, E8.2, 1x, a)') &
              & 'Inv_flipflop gmres_oce_e2e #Iter', n_iter, 'Tol ',tolerance, 'Res ',&
              & ABS(z_residual(n_iter)),MAXVAL (ABS(rhstemp(:,:))), 'ocean_gmres PROBLEM!!!!!!!!!!!!'
          ELSE
            WRITE (0, '(1x,a, I4.2, 1x, a,E8.2,1x, a,E8.2,1x, E8.2)') &
              & 'Inv_flipflop gmres_oce_e2e #Iter', n_iter, 'Tol ',tolerance, 'Res ',&
              & ABS(z_residual(n_iter)),MAXVAL (ABS(rhstemp(:,:)))
          ENDIF
        ENDIF
      END IF
    END DO
    
  END FUNCTION inverse_primal_flip_flop
  !--------------------------------------------------------------------
  
  !--------------------------------------------------------------------
  !!  results is valid only in in_domain edges
  FUNCTION lhs_primal_flip_flop( x, patch_2D, patch_3d, op_coeffs,jk,coeff, h_e) result(llhs)
    !
    TYPE(t_patch), TARGET, INTENT(in)             :: patch_2D
    TYPE(t_patch_3d ),TARGET, INTENT(in)          :: patch_3d
    REAL(wp),INTENT(inout)                        :: x(:,:)
    TYPE(t_operator_coeff),INTENT(in)             :: op_coeffs
    INTEGER ,INTENT(in)                           :: jk
    REAL(wp),INTENT(in)                           :: coeff
    REAL(wp),OPTIONAL,INTENT(in)                  :: h_e(SIZE(x,1), SIZE(x,2))!(:,:)
    REAL(wp)                                      :: llhs(SIZE(x,1), SIZE(x,2))
    
    !local variables
    REAL(wp) :: z_x_out(SIZE(x,1), 1,SIZE(x,2))!(nproma,patch_2D%nblks_e)
    REAL(wp) :: z_e(SIZE(x,1), 1,SIZE(x,2))!(nproma,patch_2D%nblks_e)
    !-----------------------------------------------------------------------
    !edges_in_domain => patch_2D%edges%in_domain
    
    z_x_out(:,1,:) = x
    
    CALL sync_patch_array(sync_e, patch_2D, x)
    
    CALL map_edges2edges_viacell_3d( patch_3d,    &
      & z_x_out(:,1,:),&
      & op_coeffs,    &
      & z_e(:,1,:),    &
    !& patch_3d%p_patch_1D(n_dom)%prism_thick_c(:,1,:),&
      & level=1)
    
    
    llhs(:,:) = coeff * z_e(:,1,:)
    !write(*,*)'max/min in', maxval(x(:,:)),minval(x(:,:))
    !rite(*,*)'max/min LHS', maxval(llhs(:,:)),minval(llhs(:,:))
    
  END FUNCTION lhs_primal_flip_flop
  !--------------------------------------------------------------------

END MODULE mo_oce_ab_timestepping_mimetic
