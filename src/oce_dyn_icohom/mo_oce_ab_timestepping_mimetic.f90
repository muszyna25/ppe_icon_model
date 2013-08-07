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
!! mpi parallelized LL
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
!---------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2  
!    modified for ICON project, DWD/MPI-M 2006
!
!---------------------------------------------------------------------------
!
!
!
USE mo_kind,                      ONLY: wp
USE mo_parallel_config,           ONLY: nproma, l_fast_sum
USE mo_math_utilities,            ONLY: t_cartesian_coordinates
USE mo_sync,                      ONLY: sync_e, sync_c, sync_patch_array
!USE mo_mpi,                       ONLY: my_process_is_mpi_parallel
USE mo_impl_constants,            ONLY: sea_boundary, &  !  sea,                          &
  &                                     max_char_length, MIN_DOLIC
USE mo_dbg_nml,                   ONLY: idbg_mxmn
USE mo_ocean_nml,                 ONLY: n_zlev, solver_tolerance, l_inverse_flip_flop,    &
  &                                     ab_const, ab_beta, ab_gam, iswm_oce,              &
  &                                     expl_vertical_velocity_diff, iforc_oce,           &
  &                                     no_tracer, l_RIGID_LID, l_edge_based,             &
  &                                     use_absolute_solver_tolerance, solver_start_tolerance, &
  &                                     solver_tolerance_decrease_ratio,                  &
  &                                     solver_max_restart_iterations,                    &
  &                                     solver_max_iter_per_restart, dhdtw_abort,         &
  &                                     l_forc_freshw

USE mo_run_config,                ONLY: dtime, ltimer
USE mo_timer,                     ONLY: timer_start, timer_stop, timer_ab_expl,           &
  &                                     timer_ab_rhs4sfc, timer_lhs
USE mo_dynamics_config,           ONLY: nold, nnew
USE mo_physical_constants,        ONLY: grav
USE mo_oce_state,                 ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, is_initial_timestep !,&
! &                                     set_lateral_boundary_values
USE mo_model_domain,              ONLY: t_patch, t_patch_3D
USE mo_ext_data_types,            ONLY: t_external_data
USE mo_ocean_gmres,               ONLY: ocean_gmres, gmres_oce_old
USE mo_exception,                 ONLY: message, finish, message_text
USE mo_util_dbg_prnt,             ONLY: dbg_print
USE mo_oce_boundcond,             ONLY: bot_bound_cond_horz_veloc, top_bound_cond_horz_veloc
USE mo_oce_thermodyn,             ONLY: calc_density, calc_internal_press
USE mo_oce_physics,               ONLY: t_ho_params
USE mo_sea_ice_types,             ONLY: t_sfc_flx
USE mo_scalar_product,            ONLY: map_edges2edges_viacell_3d, & ! map_cell2edges_3D,&
  &                                     calc_scalar_product_veloc_3D,&
!  &                                     nonlinear_coriolis_3d, nonlinear_coriolis_3d_old,&
  &                                     map_edges2edges_viacell_3d_const_z
USE mo_oce_math_operators,        ONLY: div_oce_3D, grad_fd_norm_oce_3D,&
  &                                     grad_fd_norm_oce_2d_3D, calc_thickness! , height_related_quantities
USE mo_oce_veloc_advection,       ONLY: veloc_adv_horz_mimetic, veloc_adv_vert_mimetic
 
USE mo_oce_diffusion,             ONLY: velocity_diffusion,& 
  &                                     velocity_diffusion_vert_mimetic,  &
  &                                     veloc_diffusion_vert_impl_hom
USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
USE mo_grid_config,               ONLY: n_dom
USE mo_parallel_config,           ONLY: p_test_run
IMPLICIT NONE

PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'
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
PRIVATE :: Jacobi_precon

INTEGER, PARAMETER :: top=1
CHARACTER(len=12)  :: str_module = 'oceSTEPmimet'  ! Output of module for 1 line debug
INTEGER            :: idt_src    = 1               ! Level of detail for 1 line debug

! TRUE=staggering between thermodynamic and dynamic part, offset of half timestep
! between dynamic and thermodynamic variables thermodynamic and dnamic variables are colocated in time
LOGICAL, PUBLIC,PARAMETER :: l_STAGGERED_TIMESTEP = .FALSE. 


! these are allocated once for efficeincy and used only by the lhs for the solver
REAL(wp), ALLOCATABLE, TARGET :: lhs_result(:,:)  ! (nproma,patch%nblks_c)
REAL(wp), ALLOCATABLE :: lhs_z_grad_h(:,:)
REAL(wp), ALLOCATABLE :: lhs_z_e     (:,:)
REAL(wp), ALLOCATABLE :: lhs_z_e_top (:,:)
REAL(wp), ALLOCATABLE :: lhs_div_z_c(:,:)  ! (nproma,1,patch%nblks_c)
TYPE(t_cartesian_coordinates), ALLOCATABLE :: lhs_z_grad_h_cc(:,:)

CONTAINS


!-------------------------------------------------------------------------
!>
!! !  Solves the free surface equation.
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
!!  mpi parallelized LL (2012)
!!
!
SUBROUTINE solve_free_sfc_ab_mimetic(p_patch_3D, p_os, p_ext_data, p_sfc_flx, &
    &                                  p_phys_param, timestep, p_op_coeff)
  
  TYPE(t_patch_3D ),TARGET, INTENT(INOUT)   :: p_patch_3D
  TYPE(t_hydro_ocean_state), TARGET             :: p_os
  TYPE(t_external_data), TARGET, INTENT(in)     :: p_ext_data
  TYPE(t_sfc_flx), INTENT(INOUT)                :: p_sfc_flx
  TYPE (t_ho_params)                            :: p_phys_param
  INTEGER                                       :: timestep
  !TYPE(t_int_state),TARGET,INTENT(IN)           :: p_int  
  TYPE(t_operator_coeff)                        :: p_op_coeff
  !
  !Local variables
  !
  INTEGER,PARAMETER :: nmax_iter   = 200      ! maximum number of iterations
  REAL(wp) :: tolerance =0.0_wp               ! (relative or absolute) tolerance
  INTEGER  :: n_iter                          ! actual number of iterations 
  INTEGER  :: iter_sum                        ! sum of iterations 
  INTEGER  :: jc,jb,je   ! ,jk,il_v1,il_v2,ib_v1,ib_v2
  INTEGER  :: i_startidx_c, i_endidx_c
  INTEGER  :: i_startidx_e, i_endidx_e
  REAL(wp) :: z_h_c(nproma,p_patch_3D%p_patch_2D(1)%nblks_c)
  REAL(wp) :: z_h_e(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
  LOGICAL  :: lprecon         = .FALSE.
  REAL(wp) :: z_implcoeff
  REAL(wp) :: zresidual(nmax_iter)    ! norms of the residual (convergence history);an argument of dimension at least m is required
  REAL(wp) :: residual_norm
  LOGICAL  :: l_maxiter     ! true if reached m iterations
  INTEGER  :: gmres_restart_iterations
  CHARACTER(len=max_char_length) :: string
  TYPE(t_subset_range), POINTER :: all_cells, all_edges
  TYPE(t_patch), POINTER :: patch_horz

  !CHARACTER(len=max_char_length), PARAMETER :: &
  !       & routine = ('mo_oce_ab_timestepping_mimetic:solve_free_sfc_ab_mimetic')
  !-------------------------------------------------------------------------------
  patch_horz => p_patch_3D%p_patch_2D(1)
  all_cells  => p_patch_3D%p_patch_2D(1)%cells%all
  all_edges  => p_patch_3D%p_patch_2D(1)%edges%all
  !-------------------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')
  tolerance                            = solver_tolerance
  z_h_c(1:nproma,1:patch_horz%nblks_c) = 0.0_wp
  z_h_e(1:nproma,1:patch_horz%nblks_e) = 0.0_wp

  CALL sync_patch_array(sync_c, patch_horz, p_os%p_prog(nold(1))%h)
  CALL sync_patch_array(sync_e, patch_horz, p_os%p_prog(nold(1))%vn)

  IF (is_initial_timestep(timestep) ) THEN

    IF (l_STAGGERED_TIMESTEP ) &
     & CALL calc_scalar_product_veloc_3D( p_patch_3D,&
                                      & p_os%p_prog(nold(1))%vn,&
                                      & p_os%p_prog(nold(1))%vn,&
                                      & p_os%p_diag,            &
                                      & p_op_coeff)
  ENDIF


  !Update prism thickness. The prism-thickness below the surface is
  !not updated it is initialized in construct_hydro_ocean_diag
  !with z-coordinate-thickness.
    !1) Thickness at cells
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
          p_patch_3D%p_patch_1D(n_dom)%prism_thick_c(jc,1,jb) &
          &= p_patch_3D%p_patch_1D(n_dom)%prism_thick_flat_sfc_c(jc,1,jb) +p_os%p_prog(nold(1))%h(jc,jb)
        ELSE
          !Surfacethickness over land remains zero
          !p_os%p_diag%prism_thick_c(jc,1,jb) = 0.0_wp
          p_patch_3D%p_patch_1D(n_dom)%prism_thick_c(jc,1,jb)= 0.0_wp
        ENDIF
      END DO
    END DO
    !2) Thickness at edges
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      DO je = i_startidx_e, i_endidx_e
        IF(p_patch_3D%lsm_e(je,1,jb) <= sea_boundary)THEN
          p_patch_3D%p_patch_1D(n_dom)%prism_thick_e(je,1,jb)&
          & = p_patch_3D%p_patch_1D(n_dom)%prism_thick_flat_sfc_e(je,1,jb) +p_os%p_diag%h_e(je,jb)
        ELSE
          !Surfacethickness over land remains zero
          p_patch_3D%p_patch_1D(n_dom)%prism_thick_e(je,1,jb)= 0.0_wp
        ENDIF
      END DO
    END DO 

  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=2  ! output print level (1-5, fix)
  CALL dbg_print('on entry: h-old'                ,p_os%p_prog(nold(1))%h ,str_module,idt_src)
  CALL dbg_print('on entry: h-new'                ,p_os%p_prog(nnew(1))%h ,str_module,idt_src)
  CALL dbg_print('on entry: vn-old'               ,p_os%p_prog(nold(1))%vn,str_module,idt_src)
  idt_src=1  ! output print level (1-5, fix)
  CALL dbg_print('on entry: vn-new'               ,p_os%p_prog(nnew(1))%vn,str_module,idt_src)
  !---------------------------------------------------------------------

  ! abort condition for elevation and vn:
  IF ( (maxval(p_os%p_prog(nnew(1))%h)  >  1.e6_wp) .or. &
    &  (minval(p_os%p_prog(nnew(1))%h)  < -1.e6_wp) .or. &
    &  (maxval(p_os%p_prog(nold(1))%vn) >  1.e6_wp) .or. &
    &  (minval(p_os%p_prog(nnew(1))%vn) < -1.e6_wp) ) THEN
    CALL message('Solve free surface AB mimetic: ',' INSTABLE VN or H - stop now ')
    CALL finish ('Solve free surface AB mimetic: ',' INSTABLE VN or H !!')
  END IF

  ! Apply windstress 
  CALL top_bound_cond_horz_veloc(p_patch_3D, p_os, p_op_coeff, p_sfc_flx,     &
    &                            p_os%p_aux%bc_top_u, p_os%p_aux%bc_top_v, &
    &                            p_os%p_aux%bc_top_veloc_cc)


  ! Apply bot boundary condition for horizontal velocity
  CALL bot_bound_cond_horz_veloc(p_patch_3D, p_os, p_phys_param, p_op_coeff)

  ! Calculate explicit terms of Adams-Bashforth timestepping
  IF (ltimer) CALL timer_start(timer_ab_expl)
  CALL calculate_explicit_term_ab(p_patch_3D, p_os, p_phys_param, &
    &                             is_initial_timestep(timestep), p_op_coeff)
  IF (ltimer) CALL timer_stop(timer_ab_expl)

  IF(.NOT.l_RIGID_LID)THEN

    ! Calculate RHS of surface equation
    IF (ltimer) CALL timer_start(timer_ab_rhs4sfc)
    CALL fill_rhs4surface_eq_ab(p_patch_3D, p_os, p_sfc_flx, p_op_coeff)
    IF (ltimer) CALL timer_stop(timer_ab_rhs4sfc)


    ! Solve surface equation with ocean_gmres solver
    z_h_c =0.0_wp!p_os%p_prog(nold(1))%h! 0.0_wp !potentially better choice: p_os%p_prog(nold(1))%h

    !The lhs needs different thicknesses at edges for 3D and SWE (including bathymetry)
    !IF( iswm_oce == 1 ) THEN  
    !  z_h_e = p_os%p_diag%h_e! p_os%p_diag%thick_e
    !  ELSEIF( iswm_oce /= 1 ) THEN 
    !  z_h_e = p_os%p_diag%h_e         ! #slo# 2011-02-21 bugfix (for mimetic only)
    !ENDIF 

    CALL sync_patch_array(SYNC_E, patch_horz, z_h_e)
    CALL sync_patch_array(SYNC_C, patch_horz, p_os%p_diag%thick_c)
    CALL sync_patch_array(SYNC_C, patch_horz, p_os%p_prog(nold(1))%h)

    IF (p_test_run .or. .not. l_fast_sum ) THEN
      ! the new gmres_oce uses out of order global sum for efficiency,
      ! when running in p_test_run mode use the old one
      IF(lprecon)THEN
        !p_os%p_aux%p_rhs_sfc_eq = p_os%p_aux%p_rhs_sfc_eq *patch%cells%area

        CALL gmres_oce_old( z_h_c(:,:),       &  ! arg 1 of lhs. x input is the first guess.
          &        lhs_surface_height_ab_mim, &  ! function calculating l.h.s.
          &        p_os%p_diag%thick_e,       &  ! edge thickness for LHS
          &        p_os%p_diag%thick_c,       &  ! p_os%p_diag%thick_c, & 
                                                 ! arg 6 of lhs p_os%p_prog(nold(1))%h,
                                                 ! p_os%p_diag%cons_thick_c(:,1,:),&
          &        p_os%p_prog(nold(1))%h,    &  ! arg 2 of lhs !not used
          &        p_patch_3D,                &  ! arg 3 of lhs
          &        z_implcoeff,               &  ! arg 4 of lhs
          &        p_op_coeff,                &
          &        p_os%p_aux%p_rhs_sfc_eq,   &  ! right hand side as input
          &        tolerance,                 &  ! relative tolerance
          &        .FALSE.,                   &  ! NOT absolute tolerance
          &        nmax_iter,                 &  ! max. # of iterations to do
          &        l_maxiter,                 &  ! out: .true. = not converged
          &        n_iter,                    &  ! out: # of iterations done
          &        zresidual,                 &  ! inout: the residual (array)  
          &        Jacobi_precon )

      ELSEIF(.NOT.lprecon)THEN
        CALL gmres_oce_old( z_h_c(:,:),       &  ! arg 1 of lhs. x input is the first guess.
          &        lhs_surface_height_ab_mim, &  ! function calculating l.h.s.
          &        p_os%p_diag%thick_e,       &  ! edge thickness for LHS
          &        p_os%p_diag%thick_c,       &  ! p_os%p_diag%thick_c, & 
                                                 ! arg 6 of lhs p_os%p_prog(nold(1))%h,
                                                 ! p_os%p_diag%cons_thick_c(:,1,:),&
          &        p_os%p_prog(nold(1))%h,    &  ! arg 2 of lhs !not used
          &        p_patch_3D,                &  ! arg 3 of lhs
          &        z_implcoeff,               &  ! arg 4 of lhs
          &        p_op_coeff,                &
          &        p_os%p_aux%p_rhs_sfc_eq,   &  ! right hand side as input
          &        tolerance,                 &  ! relative tolerance
          &        .FALSE.,                   &  ! NOT absolute tolerance
          &        nmax_iter,                 &  ! max. # of iterations to do
          &        l_maxiter,                 &  ! out: .true. = not converged
          &        n_iter,                    &  ! out: # of iterations done
          &        zresidual )
      ENDIF
      
      IF (l_maxiter) THEN
        CALL finish('GMRES_oce_old: solver surface equation: ','NOT YET CONVERGED !!')
      ELSE
        ! output print level idt_src used for ocean_gmres output with call message:
        IF(n_iter==0)n_iter=1
        idt_src=1
        IF (idbg_mxmn >= idt_src) THEN
          WRITE(string,'(a,i4,a,e28.20)') &
            'iteration =', n_iter,', residual =', ABS(zresidual(n_iter))
          CALL message('GMRES_oce_old: surface height',TRIM(string))
        ENDIF
      ENDIF
      
      IF(lprecon)THEN
        p_os%p_prog(nnew(1))%h = z_h_c!*patch_horz%cells%area
      ELSEIF(.NOT.lprecon)THEN
        p_os%p_prog(nnew(1))%h = z_h_c
      ENDIF

      iter_sum = n_iter*n_iter

    ELSE
    
      ! call the new gmre_oce, uses out of order global sum
      tolerance = solver_start_tolerance
      residual_norm = solver_tolerance + 1.0_wp
      gmres_restart_iterations = 0
      iter_sum                 = 0
      ! write(0,*) tolerance, solver_tolerance, residual_norm, gmres_restart_iterations, solver_max_restart_iterations
      DO  WHILE(residual_norm >= solver_tolerance .AND. gmres_restart_iterations < solver_max_restart_iterations)

        CALL ocean_gmres( z_h_c(:,:),                   &  ! arg 1 of lhs. x input is the first guess.
          &        lhs_surface_height_ab_mim,     &  ! function calculating l.h.s.
          &        p_os%p_diag%thick_e,           &  ! edge thickness for LHS
          &        p_os%p_diag%thick_c,           &  ! p_os%p_diag%thick_c, & 
          &        p_os%p_prog(nold(1))%h,        &  ! arg 2 of lhs !not used
          &        p_patch_3D,                    &  ! arg 3 of lhs
          &        z_implcoeff,                   &  ! arg 4 of lhs
          &        p_op_coeff,                    & 
          &        p_os%p_aux%p_rhs_sfc_eq,       &  ! right hand side as input
          &        tolerance,                     &  ! tolerance
          &        use_absolute_solver_tolerance, &  ! absolute/relative tolerance
          &        solver_max_iter_per_restart,   &  ! max. # of iterations to do
          &        l_maxiter,                     &  ! out: .true. = not converged
          &        n_iter,                        &  ! out: # of iterations done
          &        zresidual )


        IF(n_iter==0) THEN
          residual_norm = 0.0_wp
        ELSE
          residual_norm =  ABS(zresidual(n_iter))
        ENDIF

        iter_sum = iter_sum + n_iter
        ! output print level idt_src used for ocean_gmres output with call message:
        idt_src=2
        IF (idbg_mxmn >= idt_src) THEN
          WRITE(string,'(a,i4,a,e28.20)') &
            'ocean_gmres iteration =', n_iter,', residual =', residual_norm
          CALL message('GMRES_oce_new: surface height',TRIM(string))
        ENDIF

        IF (tolerance > solver_tolerance) &
          tolerance = MAX(tolerance * solver_tolerance_decrease_ratio, solver_tolerance)

        gmres_restart_iterations = gmres_restart_iterations + 1

      END DO ! WHILE(tolerance >= solver_tolerance)

      ! output of sum of iterations every timestep
      idt_src=0
      IF (idbg_mxmn >= idt_src) THEN
        WRITE(string,'(a,i4,a,e28.20)') &
          'SUM of ocean_gmres iteration =', iter_sum,', residual =', residual_norm
        CALL message('GMRES_oce_new: surface height',TRIM(string))
      ENDIF
 
      IF (residual_norm > solver_tolerance) &
          CALL finish('GMRES_oce_new: solver surface equation: ','NOT YET CONVERGED !!')

      p_os%p_prog(nnew(1))%h = z_h_c

    ENDIF
    !-------- end of solver ---------------

    CALL sync_patch_array(SYNC_C, patch_horz, p_os%p_prog(nnew(1))%h)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
!     z_h_c = lhs_surface_height_ab_mim( p_os%p_prog(nnew(1))%h, &
!         & p_os%p_prog(nold(1))%h, &
!         & p_patch_3D,             &
!         & z_implcoeff,            &
!         & p_os%p_diag%thick_e,    &
!         & p_os%p_diag%thick_c,    &
!         & p_op_coeff)             &
!         & -p_os%p_aux%p_rhs_sfc_eq
!     CALL dbg_print('SolvSfc: residual h-res'    ,z_h_c                  ,str_module,idt_src)
    idt_src=1  ! output print level (1-5, fix)
    CALL dbg_print('SolvSfc: after ocean_gmres: h-new',p_os%p_prog(nnew(1))%h ,str_module,idt_src)
    !---------------------------------------------------------------------

  ENDIF  ! l_rigid_lid

  ! write(0,*) "solve_free_sfc_ab_mimetic: sum(h)=", SUM(p_os%p_prog(nnew(1))%h(:,:))

END SUBROUTINE solve_free_sfc_ab_mimetic
!-------------------------------------------------------------------------  


!-------------------------------------------------------------------------
SUBROUTINE Jacobi_precon( p_jp, p_patch_3D, p_op_coeff,thick_e) !RESULT(p_jp)
!
TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
REAL(wp),INTENT(INOUT)                        :: p_jp(:,:)    ! inout for sync, dimension: (nproma,patch%nblks_c)
TYPE(t_operator_coeff),INTENT(in)             :: p_op_coeff
REAL(wp),INTENT(in)                           :: thick_e(:,:)   
!
! Left-hand side calculated from iterated height
!REAL(wp) :: p_jp(SIZE(p_x,1), SIZE(p_x,2))  ! (nproma,patch%nblks_c)
!
! local variables
REAL(wp) :: gdt2
INTEGER :: i_startidx, i_endidx
INTEGER :: jc, jb  !, je
!REAL(wp) :: z1,z2,z3
INTEGER :: edge_1_idx, edge_1_blk, edge_2_idx, edge_2_blk, edge_3_idx, edge_3_blk
REAL(wp) :: p_diag(nproma,p_patch_3D%p_patch_2D(1)%nblks_c)
TYPE(t_subset_range), POINTER :: all_cells  !, cells_in_domain, all_edges
TYPE(t_patch), POINTER :: patch     ! patch on which computation is performed
!-----------------------------------------------------------------------  
!CALL message (TRIM(routine), 'start - iteration by ocean_gmres')
patch =>  p_patch_3D%p_patch_2D(1)
!write(*,*)'inside jacobi'
write(*,*)'residual before',maxval(p_jp),minval(p_jp)!,maxvalp_jp),minval(p_jp)
all_cells => patch%cells%all
!cells_in_domain => patch%cells%in_domain
!all_edges => patch%edges%all

!p_jp(1:nproma,1:patch%nblks_c)  = 0.0_wp

gdt2 = grav*(dtime)**2


  CALL sync_patch_array(SYNC_C, patch, p_jp )


  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
    DO jc = i_startidx, i_endidx
      !IF ( v_base%lsm_c(jc,1,jb) <= sea_boundary ) THEN
      !IF(p_patch_3D%lsm_c(jc,1,jb)>=MIN_DOLIC)THEN
        edge_1_idx = patch%cells%edge_idx(jc,jb,1)
        edge_1_blk = patch%cells%edge_blk(jc,jb,1)
        edge_2_idx = patch%cells%edge_idx(jc,jb,2)
        edge_2_blk = patch%cells%edge_blk(jc,jb,2)
        edge_3_idx = patch%cells%edge_idx(jc,jb,3)
        edge_3_blk = patch%cells%edge_blk(jc,jb,3)
! 
!         p_diag(jc,jb) = patch%cells%area(jc,jb)-&
!                &gdt2*ab_gam*ab_beta*((p_op_coeff%div_coeff(jc,1,jb,1)&!*thick_e(edge_1_idx,edge_1_blk)&
!                &*p_op_coeff%grad_coeff(edge_1_idx,1,edge_1_blk)&
!                &+&
!                &p_op_coeff%div_coeff(jc,1,jb,2)&!*thick_e(edge_2_idx,edge_2_blk)&
!                &*p_op_coeff%grad_coeff(edge_2_idx,1,edge_2_blk)&
!                &+&
!                &p_op_coeff%div_coeff(jc,1,jb,3)&!*thick_e(edge_3_idx,edge_3_blk)&
!                &*p_op_coeff%grad_coeff(edge_3_idx,1,edge_3_blk)))&
!                &*patch%cells%area(jc,jb)
!
! z1=gdt2*ab_gam*ab_beta*(&
!   & p_op_coeff%div_coeff(jc,1,jb,1)*thick_e(edge_1_idx,edge_1_blk)&
!   &*p_op_coeff%grad_coeff(edge_1_idx,1,edge_1_blk))
! 
! 
! z2=gdt2*ab_gam*ab_beta*(p_op_coeff%div_coeff(jc,1,jb,2)*thick_e(edge_2_idx,edge_2_blk)&
!    &*p_op_coeff%grad_coeff(edge_2_idx,1,edge_2_blk))
! 
! z3=gdt2*ab_gam*ab_beta*(p_op_coeff%div_coeff(jc,1,jb,3)*thick_e(edge_3_idx,edge_3_blk)&
!   &*p_op_coeff%grad_coeff(edge_3_idx,1,edge_3_blk))
! 
! p_diag(jc,jb)=1.0_wp-sqrt(z1*z1+z2*z2+z3*z3)gdt2

!         p_diag(jc,jb) = (1.0_wp-&
!                &gdt2*ab_gam*ab_beta*(&
!                & p_op_coeff%div_coeff(jc,1,jb,1)*thick_e(edge_1_idx,edge_1_blk)&
!                &*p_op_coeff%grad_coeff(edge_1_idx,1,edge_1_blk)&
!                &+&
!                &p_op_coeff%div_coeff(jc,1,jb,2)*thick_e(edge_2_idx,edge_2_blk)&
!                &*p_op_coeff%grad_coeff(edge_2_idx,1,edge_2_blk)&
!                &+&
!                &p_op_coeff%div_coeff(jc,1,jb,3)*thick_e(edge_3_idx,edge_3_blk)&
!                &*p_op_coeff%grad_coeff(edge_3_idx,1,edge_3_blk)))*patch%cells%area(jc,jb)/gdt2
!   
 p_diag(jc,jb)=1.0_wp- &
 &gdt2*ab_gam*ab_beta*(patch%edges%primal_edge_length(edge_1_idx,edge_1_blk)&
 & *patch%edges%inv_dual_edge_length(edge_1_idx,edge_1_blk)&
 &+patch%edges%primal_edge_length(edge_2_idx,edge_2_blk)&
& *patch%edges%inv_dual_edge_length(edge_2_idx,edge_2_blk)&
 &+patch%edges%primal_edge_length(edge_3_idx,edge_3_blk) &
&*patch%edges%inv_dual_edge_length(edge_1_idx,edge_1_blk))&
&/(patch%cells%area(jc,jb)*gdt2)
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
write(*,*)'diag element',maxval(p_diag),minval(p_diag)
write(*,*)'residual after',maxval(p_jp),minval(p_jp)!,maxvalp_jp),minval(p_jp)
END SUBROUTINE Jacobi_precon
!-------------------------------------------------------------------------  
!
!  
!>
!! Computation of velocity predictor in Adams-Bashforth timestepping.
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
!!  mpi parallelized LL
SUBROUTINE calculate_explicit_term_ab( p_patch_3D, p_os, p_phys_param,&
                                     & l_initial_timestep, p_op_coeff)

  !TYPE(t_patch), TARGET, INTENT(in)             :: patch
  TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
  TYPE(t_hydro_ocean_state), TARGET             :: p_os
  TYPE (t_ho_params)                            :: p_phys_param
  !TYPE(t_int_state),TARGET,INTENT(IN), OPTIONAL :: p_int
  LOGICAL,INTENT(IN)                            :: l_initial_timestep
  TYPE(t_operator_coeff)                        :: p_op_coeff
  !
  !local variables
  !
  INTEGER  :: je, jk, jb
  REAL(wp) :: gdt
  REAL(wp) :: z_gradh_e(nproma,1,p_patch_3D%p_patch_2D(n_dom)%nblks_e)
  REAL(wp) :: z_e(nproma,n_zlev,p_patch_3D%p_patch_2D(n_dom)%nblks_e)
  TYPE(t_subset_range), POINTER :: edges_in_domain, all_edges
  INTEGER  :: i_startidx_e, i_endidx_e, dolic_e
  TYPE(t_patch), POINTER :: patch_horz
  !CHARACTER(len=max_char_length), PARAMETER :: &
  !  &       routine = ('mo_oce_ab_timestepping_mimetic:calculate_explicit_term_ab')
  !-----------------------------------------------------------------------  
  !CALL message (TRIM(routine), 'start')        

  patch_horz    => p_patch_3D%p_patch_2D(n_dom)
  edges_in_domain => p_patch_3D%p_patch_2D(n_dom)%edges%in_domain
  all_edges       => p_patch_3D%p_patch_2D(n_dom)%edges%all

  z_gradh_e(:,:,:) = 0.0_wp
  gdt              = grav*dtime

  !---------------------------------------------------------------------
  ! STEP 1: calculate gradient of surface height at previous timestep
  !---------------------------------------------------------------------

  CALL grad_fd_norm_oce_2d_3D( p_os%p_prog(nold(1))%h,     &
         &                  patch_horz,                  &
         &                  p_op_coeff%grad_coeff(:,1,:),  &
         &                  z_gradh_e(:,1,:))
  CALL sync_patch_array(sync_e, patch_horz, z_gradh_e(:,1,:))

  !---------------------------------------------------------------------
  ! STEP 2: horizontal advection
  !---------------------------------------------------------------------

  IF(l_initial_timestep)THEN
    CALL veloc_adv_horz_mimetic( p_patch_3D,         &
           &             p_os%p_prog(nold(1))%vn,    &
           &             p_os%p_prog(nold(1))%vn,    &
           &             p_os%p_diag,                &
           &             p_os%p_diag%veloc_adv_horz, &
           &             p_op_coeff)
  ELSE
    CALL veloc_adv_horz_mimetic( p_patch_3D,         &
           &             p_os%p_prog(nold(1))%vn,    &
           &             p_os%p_prog(nnew(1))%vn,    &
           &             p_os%p_diag,                &
           &             p_os%p_diag%veloc_adv_horz, &
           &             p_op_coeff)
  ENDIF

  !---------------------------------------------------------------------
  ! STEP 3: compute 3D contributions: gradient of hydrostatic pressure and vertical velocity advection
  !---------------------------------------------------------------------

  IF ( iswm_oce /= 1 ) THEN
    ! calculate density from EOS using temperature and salinity at timelevel n
     CALL calc_density( p_patch_3D,                                    &
            &           p_os%p_prog(nold(1))%tracer(:,:,:,1:no_tracer),&
            &           p_os%p_diag%rho(:,:,:) )

    ! calculate hydrostatic pressure from density at timelevel nc
     CALL calc_internal_press( p_patch_3D,                        &
            &                  p_os%p_diag%rho,                   &
            &                  p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c,&
            &                  p_os%p_prog(nold(1))%h,            &
            &                  p_os%p_diag%press_hyd)

    ! calculate gradient of hydrostatic pressure in 3D
    CALL grad_fd_norm_oce_3D( p_os%p_diag%press_hyd,  &
           &                  p_patch_3D,             &
           &                  p_op_coeff%grad_coeff,  &
           &                  p_os%p_diag%press_grad)
    CALL sync_patch_array(SYNC_E, patch_horz, p_os%p_diag%press_grad)

    ! calculate vertical velocity advection
    CALL veloc_adv_vert_mimetic( p_patch_3D,                 &
          &                      p_os%p_diag,p_op_coeff,     &
          &                      p_os%p_diag%veloc_adv_vert )

    ! calculate vertical velocity diffusion
    !   For the alternative choice "expl_vertical_velocity_diff==1" see couples of
    !   lines below below
    IF (expl_vertical_velocity_diff==0) THEN
        CALL velocity_diffusion_vert_mimetic( p_patch_3D,     &
        &                             p_os%p_diag,            &
        &                             p_os%p_aux,p_op_coeff,  &
        &                             p_phys_param,           &
        &                             p_os%p_diag%laplacian_vert)
    ENDIF
  ELSE  !  iswm_oce=1
    p_os%p_diag%veloc_adv_vert = 0.0_wp
    p_os%p_diag%laplacian_vert = 0.0_wp
  ENDIF

  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src = 3  ! output print level (1-5, fix)
  CALL dbg_print('old height gradient'       ,z_gradh_e                 ,str_module,idt_src)
  CALL dbg_print('horizontal advection'      ,p_os%p_diag%veloc_adv_horz,str_module,idt_src)
  CALL dbg_print('density'                   ,p_os%p_diag%rho           ,str_module,idt_src, &
    & in_subset=p_patch_3D%p_patch_2D(n_dom)%cells%owned)
  CALL dbg_print('internal pressure'         ,p_os%p_diag%press_hyd     ,str_module,idt_src)
  CALL dbg_print('internal press grad'       ,p_os%p_diag%press_grad    ,str_module,idt_src)
  idt_src = 4  ! output print level (1-5, fix)
  CALL dbg_print('kinetic energy'            ,p_os%p_diag%kin           ,str_module,idt_src)
  CALL dbg_print('vertical advection'        ,p_os%p_diag%veloc_adv_vert,str_module,idt_src)
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! STEP 4: compute harmonic or biharmoic laplacian diffusion of velocity.
  !         This term is discretized explicitly. Order and form of the laplacian
  !         are determined in mo_oce_diffusion according to namelist settings
  !---------------------------------------------------------------------

  CALL velocity_diffusion(p_patch_3D,              &
                        & p_os%p_prog(nold(1))%vn, &
                        & p_phys_param,            &
                        & p_os%p_diag,p_op_coeff,  &
                        & p_os%p_diag%laplacian_horz)

  CALL sync_patch_array(SYNC_E, patch_horz, p_os%p_diag%laplacian_horz)

  IF (L_INVERSE_FLIP_FLOP) THEN

    IF ( iswm_oce /= 1 ) THEN
      z_e = inverse_primal_flip_flop(patch_horz,p_patch_3D,p_op_coeff, p_os%p_diag%veloc_adv_horz, p_os%p_diag%h_e)
    ELSE
      z_e = inverse_primal_flip_flop(patch_horz,p_patch_3D, p_op_coeff,p_os%p_diag%veloc_adv_horz, p_os%p_diag%thick_e)
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src = 5  ! output print level (1-5, fix)
    CALL dbg_print('bef.dual-flip-fl: LaPlaHorz',p_os%p_diag%laplacian_horz,str_module,idt_src)
    !---------------------------------------------------------------------

    IF(l_STAGGERED_TIMESTEP)THEN

      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        DO jk = 1, n_zlev

          p_os%p_aux%g_n(i_startidx_e:i_endidx_e, jk, jb) = &      
            &  - z_e(i_startidx_e:i_endidx_e,jk,jb)  &
            &  - p_os%p_diag%veloc_adv_vert(i_startidx_e:i_endidx_e,jk,jb)  &
            &  + p_os%p_diag%laplacian_horz(i_startidx_e:i_endidx_e,jk,jb)  &
            &  + p_os%p_diag%laplacian_vert(i_startidx_e:i_endidx_e,jk,jb)
        END DO
      END DO

    ELSEIF(.NOT.l_STAGGERED_TIMESTEP)THEN
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        DO jk = 1, n_zlev
          p_os%p_aux%g_n(i_startidx_e:i_endidx_e, jk, jb) = &
            &  -p_os%p_diag%press_grad(i_startidx_e:i_endidx_e, jk, jb)       &
            &  - z_e(i_startidx_e:i_endidx_e, jk, jb)  &
            &  - p_os%p_diag%veloc_adv_vert(i_startidx_e:i_endidx_e, jk, jb)  &
            &  + p_os%p_diag%laplacian_horz(i_startidx_e:i_endidx_e, jk, jb)  &
            &  + p_os%p_diag%laplacian_vert(i_startidx_e:i_endidx_e, jk, jb)
        END DO
      END DO
    ENDIF

    CALL sync_patch_array(SYNC_E, patch_horz, p_os%p_aux%g_n)

  ELSEIF(.NOT.(L_INVERSE_FLIP_FLOP))THEN

    IF(l_STAGGERED_TIMESTEP)THEN
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        DO jk = 1, n_zlev
          p_os%p_aux%g_n(i_startidx_e:i_endidx_e, jk, jb) =&!-p_os%p_diag%press_grad(:,jk,:)      &
            &  - p_os%p_diag%veloc_adv_horz(i_startidx_e:i_endidx_e, jk, jb)  &
            &  - p_os%p_diag%veloc_adv_vert(i_startidx_e:i_endidx_e, jk, jb)  &
            &  + p_os%p_diag%laplacian_horz(i_startidx_e:i_endidx_e, jk, jb)  &
            &  + p_os%p_diag%laplacian_vert(i_startidx_e:i_endidx_e, jk, jb)
        END DO
      END DO
    ELSEIF(.NOT.l_STAGGERED_TIMESTEP)THEN
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        DO jk = 1, n_zlev
          p_os%p_aux%g_n(i_startidx_e:i_endidx_e, jk, jb) = &
            & - p_os%p_diag%press_grad(i_startidx_e:i_endidx_e, jk, jb)      &
            & - p_os%p_diag%veloc_adv_horz(i_startidx_e:i_endidx_e, jk, jb)  &
            & - p_os%p_diag%veloc_adv_vert(i_startidx_e:i_endidx_e, jk, jb)  &
            & + p_os%p_diag%laplacian_horz(i_startidx_e:i_endidx_e, jk, jb)  &
            & + p_os%p_diag%laplacian_vert(i_startidx_e:i_endidx_e, jk, jb)
        END DO
      END DO
    ENDIF

  ENDIF!(L_INVERSE_FLIP_FLOP)

  IF(l_initial_timestep)THEN
    p_os%p_aux%g_nimd(1:nproma,1:n_zlev,1:patch_horz%nblks_e) = p_os%p_aux%g_n(1:nproma,1:n_zlev,1:patch_horz%nblks_e)
  ELSE
    p_os%p_aux%g_nimd(1:nproma,1:n_zlev,1:patch_horz%nblks_e) &
    &= (1.5_wp+AB_const)* p_os%p_aux%g_n(1:nproma,1:n_zlev,1:patch_horz%nblks_e)&
    &- (0.5_wp+AB_const)* p_os%p_aux%g_nm1(1:nproma,1:n_zlev,1:patch_horz%nblks_e)
  ENDIF

  IF ( iswm_oce /= 1) THEN

    IF(.NOT.l_RIGID_LID)THEN

      IF(l_STAGGERED_TIMESTEP)THEN
        DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
          DO jk = 1, n_zlev
            DO je = i_startidx_e, i_endidx_e

              IF(p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)>=MIN_DOLIC)THEN
              !IF(v_base%lsm_e(je,jk,jb) <= sea_boundary ) THEN
                p_os%p_diag%vn_pred(je,jk,jb) = p_os%p_prog(nold(1))%vn(je,jk,jb)    &
                &                           + dtime*(p_os%p_aux%g_nimd(je,jk,jb)     &
                &                                   -p_os%p_diag%press_grad(je,jk,jb)  &
                &                           - (1.0_wp-ab_beta) * grav*z_gradh_e(je,1,jb))
              ENDIF
            END DO
          END DO
        END DO

      ELSEIF(.NOT.l_STAGGERED_TIMESTEP)THEN

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src = 5  ! output print level (1-5, fix)
      CALL dbg_print('Bef fin term: vn_pred'      ,p_os%p_diag%vn_pred       ,str_module,idt_src)
      CALL dbg_print('Bef fin term: g_nimd'       ,p_os%p_aux%g_nimd         ,str_module,idt_src)
      CALL dbg_print('Bef fin term: gradh_e'      ,z_gradh_e                 ,str_module,idt_src)
      !---------------------------------------------------------------------

        DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
          DO jk = 1, n_zlev
            DO je = i_startidx_e, i_endidx_e

              IF(p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)>=MIN_DOLIC)THEN
              !IF(p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary)THEN
                p_os%p_diag%vn_pred(je,jk,jb) = p_os%p_prog(nold(1))%vn(je,jk,jb)  &
                &                             + dtime*(p_os%p_aux%g_nimd(je,jk,jb) &
                &                             - (1.0_wp-ab_beta) * grav*z_gradh_e(je,1,jb))
              ENDIF
            END DO
          END DO
        END DO
      ENDIF!Staggered

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src = 5  ! output print level (1-5, fix)
      CALL dbg_print('Aft fin term: vn_pred'      ,p_os%p_diag%vn_pred       ,str_module,idt_src)
      !---------------------------------------------------------------------

    ELSEIF(l_RIGID_LID)THEN

      IF(l_STAGGERED_TIMESTEP)THEN

        DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
          DO jk = 1, n_zlev
            DO je = i_startidx_e, i_endidx_e

              IF(p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)>=MIN_DOLIC)THEN
                p_os%p_diag%vn_pred(je,jk,jb) = p_os%p_prog(nold(1))%vn(je,jk,jb)  &
                &                             + dtime*(p_os%p_aux%g_nimd(je,jk,jb) &
                &                             - p_os%p_diag%press_grad(je,jk,jb))
              ENDIF

            END DO
          END DO
        END DO

      ELSEIF(.NOT.l_STAGGERED_TIMESTEP)THEN

        DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
          DO jk = 1, n_zlev
            DO je = i_startidx_e, i_endidx_e

              IF(p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)>=MIN_DOLIC)THEN
                p_os%p_diag%vn_pred(je,jk,jb) = p_os%p_prog(nold(1))%vn(je,jk,jb)&
                &                             + dtime*p_os%p_aux%g_nimd(je,jk,jb)
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
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
        DO je = i_startidx_e, i_endidx_e

          IF(p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)>=MIN_DOLIC)THEN
            p_os%p_diag%vn_pred(je,1,jb) =  p_os%p_diag%vn_pred(je,1,jb)      &
            &                            + dtime*p_os%p_aux%bc_top_vn(je,jb)&
            &/p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_e(je,1,jb)

            dolic_e=p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)
            p_os%p_diag%vn_pred(je,dolic_e,jb)    &
            & = p_os%p_diag%vn_pred(je,dolic_e,jb)&
            & - dtime*p_os%p_aux%bc_bot_vn(je,jb)               &
            &/p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_e(je,dolic_e,jb)
          ENDIF

        END DO
      END DO
    ENDIF 

  !In the SW-case the external forcing is applied as volume force.
  !This force is stored in data type top-boundary-condition. 
  ELSEIF ( iswm_oce == 1)THEN! .AND. iforc_oce==11) THEN

    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      DO jk = 1, n_zlev
        DO je = i_startidx_e, i_endidx_e
          !IF(p_patch_3D%lsm_c(je,jk,jb) <= sea_boundary)THEN
            p_os%p_diag%vn_pred(je,jk,jb) = (p_os%p_prog(nold(1))%vn(je,jk,jb)        &
             &                            + dtime*(p_os%p_aux%g_nimd(je,jk,jb)        &
             &                            - (1.0_wp-ab_beta)*grav*z_gradh_e(je,1,jb)))
        END DO
      END DO
    END DO

    !In case of Shallow-water with forcing and or damping
    IF ( iforc_oce/=10) THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
        DO jk = 1, n_zlev
          DO je = i_startidx_e, i_endidx_e
            !IF(p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary)THEN
              p_os%p_diag%vn_pred(je,jk,jb) = (p_os%p_diag%vn_pred(je,jk,jb) &
               &                            + p_os%p_aux%bc_top_vn(je,jb)    &
               &                            - p_os%p_aux%bc_bot_vn(je,jb))   
            !ENDIF
          END DO
        END DO
      END DO
    ENDIF 
  ENDIF

  !In 3D case and if implicit vertical velocity diffusion is chosen
  IF(iswm_oce /= 1.AND.expl_vertical_velocity_diff==1)THEN

    !Surface forcing is implemented as volume forcing in top layer.
    !In this case homogeneous boundary conditions for vertical Laplacian

     CALL veloc_diffusion_vert_impl_hom( p_patch_3D,             &
                                     & p_os%p_diag%vn_pred,      &
                                     & p_os%p_diag%h_e,          &
                                     & p_phys_param%A_veloc_v,   &
                                     & p_op_coeff,               &
                                     & p_os%p_diag%vn_impl_vert_diff)
    IF(l_RIGID_LID)THEN
      p_os%p_diag%vn_pred(1:nproma,1:n_zlev,1:patch_horz%nblks_e) &
      &= p_os%p_diag%vn_impl_vert_diff(1:nproma,1:n_zlev,1:patch_horz%nblks_e)
    ENDIF

  ENDIF
  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=3  ! output print level (1-5, fix)
  CALL dbg_print('vn(nold)'                  ,p_os%p_prog(nold(1))%vn       ,str_module,idt_src)
  CALL dbg_print('VelocAdvHorizontal'        ,p_os%p_diag%veloc_adv_horz    ,str_module,idt_src)
  CALL dbg_print('VelocLaPlac horizontal'    ,p_os%p_diag%laplacian_horz    ,str_module,idt_src)

  IF (expl_vertical_velocity_diff == 1 .AND. iswm_oce /= 1) THEN
   CALL dbg_print('ImplVelocDiff vertical'    ,p_os%p_diag%vn_impl_vert_diff ,str_module,idt_src)
  ELSE
  CALL dbg_print('VelocLaPlac vertical'      ,p_os%p_diag%laplacian_vert    ,str_module,idt_src)
  ENDIF
  IF (l_inverse_flip_flop) &
& CALL dbg_print('dual-flip-flop Adv. horz'  ,z_e                           ,str_module,idt_src)
  idt_src=4  ! output print level (1-5, fix)
  CALL dbg_print('G_n+1/2 - g_nimd'          ,p_os%p_aux%g_nimd             ,str_module,idt_src)
  CALL dbg_print('G_n'                       ,p_os%p_aux%g_n                ,str_module,idt_src)
  CALL dbg_print('G_n-1'                     ,p_os%p_aux%g_nm1              ,str_module,idt_src)
  idt_src=2  ! output print level (1-5, fix)
  CALL dbg_print('vn_pred'                   ,p_os%p_diag%vn_pred           ,str_module,idt_src)
  !---------------------------------------------------------------------

  !-------------------------------------------
  !CALL message (TRIM(routine), 'end')        

END SUBROUTINE calculate_explicit_term_ab
!-------------------------------------------------------------------------  
!>
!! 
!!  Calculation of right-hand side of elliptic surface equation.
!!  This is used in semi implicit timelevel stepping. 
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
  !!  mpi parallelized LL
SUBROUTINE fill_rhs4surface_eq_ab( p_patch_3D, p_os, p_sfc_flx, p_op_coeff)
!
! Patch on which computation is performed
TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
!
! Type containing ocean state
TYPE(t_hydro_ocean_state), TARGET :: p_os
TYPE(t_sfc_flx), INTENT(IN)       :: p_sfc_flx
TYPE(t_operator_coeff)            :: p_op_coeff
!
!  local variables
!
INTEGER :: i_startidx_c, i_endidx_c
INTEGER :: i_startidx_e, i_endidx_e
INTEGER :: jc, jb, jk, je
INTEGER :: i_dolic_c,i_dolic_e
REAL(wp) :: gdt2    !, delta_z
REAL(wp) :: z_e2D(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
REAL(wp) :: z_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
REAL(wp) :: div_z_depth_int_c(nproma,p_patch_3D%p_patch_2D(1)%nblks_c)
REAL(wp) :: div_z_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)
REAL(wp) :: z_vn_ab(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)

!TYPE(t_cartesian_coordinates) :: z_u_pred_depth_int_cc(nproma,patch%nblks_c)
TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain, all_edges
TYPE(t_patch), POINTER :: patch_horz
!REAL(wp) :: thick
!CHARACTER(len=max_char_length), PARAMETER :: &
!       & routine = ('mo_oce_ab_timestepping_mimetic:fill_rhs4surface_eq_ab')
!-------------------------------------------------------------------------------
!CALL message (TRIM(routine), 'start') 

  patch_horz      => p_patch_3D%p_patch_2D(1)
  all_cells       => p_patch_3D%p_patch_2D(1)%cells%all
  cells_in_domain => p_patch_3D%p_patch_2D(1)%cells%in_domain
  all_edges       => p_patch_3D%p_patch_2D(1)%edges%all

  gdt2 = grav*dtime*dtime

  z_vn_ab(1:nproma,1:n_zlev,1:patch_horz%nblks_e)  = 0.0_wp
  z_e2D(1:nproma,1:patch_horz%nblks_e)             = 0.0_wp
  z_e(1:nproma,1:n_zlev,1:patch_horz%nblks_e)      = 0.0_wp
  div_z_depth_int_c(1:nproma,1:patch_horz%nblks_c) = 0.0_wp
  div_z_c(1:nproma,1:n_zlev,1:patch_horz%nblks_c)  = 0.0_wp
  

  ! LL: this should not be required
  CALL sync_patch_array(SYNC_E, patch_horz, p_os%p_diag%vn_pred)
  CALL sync_patch_array(SYNC_E, patch_horz, p_os%p_prog(nold(1))%vn)
  CALL sync_patch_array(SYNC_E, patch_horz, p_os%p_diag%vn_impl_vert_diff)

  IF(iswm_oce == 1)THEN
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      DO je = i_startidx_e, i_endidx_e
        DO jk=1,n_zlev 
          !IF(p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary)THEN
          z_vn_ab(je,jk,jb)=ab_gam*p_os%p_diag%vn_pred(je,jk,jb)&
          &+ (1.0_wp -ab_gam)* p_os%p_prog(nold(1))%vn(je,jk,jb)
          !ENDIF
        END DO
      ENDDO
    END DO
  ELSEIF(iswm_oce /= 1)THEN

    IF(expl_vertical_velocity_diff==1)THEN
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
        DO je = i_startidx_e, i_endidx_e
          i_dolic_e =  p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)! v_base%dolic_e(je,jb)
          DO jk=1,i_dolic_e 
            !IF(p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary)THEN
            z_vn_ab(je,jk,jb)=ab_gam*p_os%p_diag%vn_impl_vert_diff(je,jk,jb)&
            &+ (1.0_wp -ab_gam)* p_os%p_prog(nold(1))%vn(je,jk,jb)
            !ENDIF
          END DO
        ENDDO
      END DO

    ELSEIF(expl_vertical_velocity_diff==0)THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
        DO je = i_startidx_e, i_endidx_e
          i_dolic_e =  p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)! v_base%dolic_e(je,jb)
          DO jk=1,i_dolic_e 
            !IF(p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary)THEN
            z_vn_ab(je,jk,jb)=ab_gam*p_os%p_diag%vn_pred(je,jk,jb)&
            &+ (1.0_wp -ab_gam)* p_os%p_prog(nold(1))%vn(je,jk,jb)
            !ENDIF
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
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
        DO je = i_startidx_e, i_endidx_e
          i_dolic_e =  p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)
          DO jk=1,i_dolic_e 
            z_e(je,jk,jb)= z_vn_ab(je,jk,jb)*p_patch_3D%p_patch_1D(1)%prism_thick_e(je,jk,jb)
          END DO
        ENDDO
      END DO

    ELSEIF( iswm_oce == 1 ) THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
        DO je = i_startidx_e, i_endidx_e
          i_dolic_e =  p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)
          DO jk=1,i_dolic_e 
            z_e(je,jk,jb)= z_vn_ab(je,jk,jb)*p_os%p_diag%thick_e(je,jb)
          END DO
        ENDDO
      END DO

    ENDIF
  
  ! !-------------------------------------------------------------------------------
  ELSEIF(.NOT. l_edge_based)THEN!NOT EDGE-BASED
  ! !-------------------------------------------------------------------------------  

    IF( iswm_oce /= 1 ) THEN !the 3D case

      CALL map_edges2edges_viacell_3d_const_z( p_patch_3D, z_vn_ab, p_op_coeff, z_e )
  
    ELSEIF( iswm_oce == 1 ) THEN
  !    CALL map_edges2edges_viacell_3D( p_patch_3D,    &
  !                                    & z_vn_ab(:,1,:),&
  !                                    & p_op_coeff,    &
  !                                    & z_e(:,1,:),    &
  !                                    & p_os%p_diag%thick_c, level=1)
      CALL map_edges2edges_viacell_3d_const_z( p_patch_3D, z_vn_ab(:,1,:), p_op_coeff, z_e(:,1,:) )

    ENDIF!( iswm_oce == 1 )

  ENDIF!EDGE-BASED

  CALL div_oce_3d( z_e, patch_horz,p_op_coeff%div_coeff, div_z_c, subset_range=cells_in_domain )
  
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
    DO jc = i_startidx_c, i_endidx_c
      i_dolic_c = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
      DO jk=1,i_dolic_c 
        IF(p_patch_3D%lsm_c(jc,jk,jb) <= sea_boundary)THEN
          div_z_depth_int_c(jc,jb)=div_z_depth_int_c(jc,jb) +div_z_c(jc,jk,jb)
        ENDIF
      END DO
    END DO
  END DO

  !-------------------------------------------------------------------------
  ! Apply net surface freshwater flux to elevation - incorrect?
 !IF(l_forc_freshw)THEN
 !  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
 !    CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
 !    DO jc = i_startidx_c, i_endidx_c
 !      IF(p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
 !        p_os%p_aux%p_rhs_sfc_eq(jc,jb) = ((p_os%p_prog(nold(1))%h(jc,jb)     &
 !                                       & - dtime*(div_z_depth_int_c(jc,jb) + &
 !                                       &          p_sfc_flx%forc_fwfx(jc,jb)) )/gdt2)
 !      ENDIF
 !    ENDDO
 !  END DO

 !ELSEIF(.NOT.l_forc_freshw)THEN

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
          p_os%p_aux%p_rhs_sfc_eq(jc,jb) = ((p_os%p_prog(nold(1))%h(jc,jb)&
                                         & - dtime*div_z_depth_int_c(jc,jb))/gdt2)
        ENDIF
      ENDDO
    END DO
 !ENDIF

  CALL sync_patch_array(SYNC_C, patch_horz, p_os%p_aux%p_rhs_sfc_eq )

  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=3  ! output print level (1-5, fix)
  CALL dbg_print('RHS thick_e'               ,p_os%p_diag%thick_e      ,str_module,idt_src)
  CALL dbg_print('RHS z_vn_ab'               ,z_vn_ab                  ,str_module,idt_src) 
  CALL dbg_print('RHS z_e'                   ,z_e                      ,str_module,idt_src)
  CALL dbg_print('RHS div_z_depth_int_c'     ,div_z_depth_int_c        ,str_module,idt_src)
  idt_src=2  ! output print level (1-5, fix)
  CALL dbg_print('RHS final'                 ,p_os%p_aux%p_rhs_sfc_eq  ,str_module,idt_src)
  !---------------------------------------------------------------------

END SUBROUTINE fill_rhs4surface_eq_ab
!-------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------
SUBROUTINE init_ho_lhs_fields_mimetic(p_patch_3D)

  TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D

  TYPE(t_patch), POINTER :: patch     ! patch on which computation is performed
  INTEGER :: return_status
  
  patch         => p_patch_3D%p_patch_2D(1)

  ALLOCATE(lhs_result(nproma,patch%nblks_c), &
    lhs_z_grad_h(nproma,patch%nblks_e),     &
    lhs_z_e     (nproma,patch%nblks_e),     &
    lhs_z_e_top (nproma,patch%nblks_e),     &
    lhs_div_z_c (nproma,patch%nblks_c),     & 
    lhs_z_grad_h_cc(nproma,patch%nblks_c),  &
    stat = return_status)

    IF (return_status > 0) &
      CALL finish("mo_oce_ab_timestepping_mimetic:init_ho_lhs_fields", "Allocation failed")
    
    lhs_result(:,:)   = 0.0_wp
    lhs_z_grad_h(:,:) = 0.0_wp
    lhs_z_e     (:,:) = 0.0_wp
    lhs_z_e_top (:,:) = 0.0_wp
    lhs_div_z_c (:,:) = 0.0_wp
    lhs_z_grad_h_cc(:,:)%x(1) = 0.0_wp 
    lhs_z_grad_h_cc(:,:)%x(2) = 0.0_wp
    lhs_z_grad_h_cc(:,:)%x(3) = 0.0_wp
    
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
!!  mpi parallelized, the result is NOT synced. Should be done in the calling method if required
!-------------------------------------------------------------------------
FUNCTION lhs_surface_height_ab_mim( p_x, h_old, p_patch_3D,coeff, thickness_e,&
                                  & thickness_c,p_op_coeff) RESULT(p_lhs)

TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
REAL(wp),    INTENT(INOUT)       :: p_x(:,:)    ! inout for sync, dimension: (nproma,patch%nblks_c)
REAL(wp),    INTENT(IN)          :: h_old(:,:)
REAL(wp),    INTENT(in)          :: coeff
TYPE(t_operator_coeff),INTENT(in):: p_op_coeff
REAL(wp),    INTENT(in)          :: thickness_e(:,:) 
REAL(wp),    INTENT(in)          :: thickness_c(:,:) !thickness of fluid column
!  these are small (2D) arrays and allocated once for efficiency
! Left-hand side calculated from iterated height
!
REAL(wp) :: p_lhs(SIZE(p_x,1), SIZE(p_x,2))  ! (nproma,p_patch%nblks_c)

! local variables,
REAL(wp) :: gdt2_inv, gam_times_beta
INTEGER :: i_startidx, i_endidx
INTEGER :: jc, jb, je
TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
TYPE(t_patch), POINTER :: patch     ! patch on which computation is performed
  !-----------------------------------------------------------------------  
  IF (ltimer) CALL timer_start(timer_lhs)
  !-----------------------------------------------------------------------  
  patch           => p_patch_3D%p_patch_2D(1)
  cells_in_domain => patch%cells%in_domain
  edges_in_domain => patch%edges%in_domain


  gdt2_inv       = 1.0_wp / (grav*(dtime)**2)
  gam_times_beta = ab_gam * ab_beta

  p_lhs   (1:nproma,patch%nblks_c)  = 0.0_wp

  CALL sync_patch_array(SYNC_C, patch, p_x )

  !Step 1) Calculate gradient of iterated height.
  CALL grad_fd_norm_oce_2d_3D( p_x, &
           &                   patch,                       &
           &                   p_op_coeff%grad_coeff(:,1,:),&
           &                   lhs_z_grad_h(:,:))

  ! the result lhs_z_grad_h is computed on in_domain edges
  CALL sync_patch_array(SYNC_E, patch, lhs_z_grad_h(:,:) )


  !TODO check
  IF(l_edge_based)THEN

    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
      DO je = i_startidx, i_endidx
          lhs_z_e(je,jb) = lhs_z_grad_h(je,jb) * thickness_e(je,jb)
!         IF(p_patch_3D%lsm_e(je,1,jb) > sea_boundary)THEN
!           IF (lhs_z_e(je,jb) /= 0.0_wp) &
!             CALL finish("lhs_surface_height_ab_mim", "lhs_z_e(je,jb) /= 0 on land")
!         ENDIF

      END DO
    END DO
    ! no need to sync since we will compute only cells in domain
    !CALL sync_patch_array(SYNC_E, patch, lhs_z_e(:,:) )

  ELSE  !IF(.NOT.l_edge_based)THEN

   ! the map_edges2edges_viacell_3d_const_z should be changes to calculate only
   ! on in_domai_edges. Still, edge valkues need to be synced
   CALL sync_patch_array(SYNC_E, patch, lhs_z_grad_h(:,:) )


    IF( iswm_oce /= 1 ) THEN

      CALL map_edges2edges_viacell_3d_const_z( p_patch_3D, lhs_z_grad_h(:,:), p_op_coeff, lhs_z_e(:,:))

    ELSEIF( iswm_oce == 1 ) THEN

      !CALL map_edges2edges_viacell_3D( p_patch_3D, lhs_z_grad_h, p_op_coeff, lhs_z_e,thickness_c, level=top)
      CALL map_edges2edges_viacell_3d_const_z( p_patch_3D, lhs_z_grad_h(:,:), p_op_coeff, lhs_z_e(:,:))
    ENDIF!( iswm_oce == 1 )

  ENDIF ! l_edge_based

  !Step 3) Calculate divergence
  CALL div_oce_3D( lhs_z_e, patch, p_op_coeff%div_coeff, lhs_div_z_c, &
    & level=top, subset_range=cells_in_domain  )

  !Step 4) Finalize LHS calculations
  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
    DO jc = i_startidx, i_endidx

      !p_lhs(jc,jb) =(p_x(jc,jb) - gdt2 * ab_gam * ab_beta * lhs_div_z_c(jc,jb)) / gdt2
      p_lhs(jc,jb) = p_x(jc,jb) * gdt2_inv - gam_times_beta * lhs_div_z_c(jc,jb)

      IF(p_patch_3D%lsm_c(jc,1,jb) > sea_boundary) THEN
        IF (p_lhs(jc,jb) /= 0.0_wp) &
          CALL finish("lhs_surface_height_ab_mim", "p_lhs(jc,jb) /= 0 on land")
      ENDIF
    END DO
  END DO

  IF (ltimer) CALL timer_stop(timer_lhs)
  !-----------------------------------------------------------------------  

END FUNCTION lhs_surface_height_ab_mim
!-------------------------------------------------------------------------  

!-------------------------------------------------------------------------  
!>
!! Computation of new velocity in Adams-Bashforth timestepping.
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
  !!  mpi parallelized LL
SUBROUTINE calc_normal_velocity_ab_mimetic(p_patch_3D,p_os, p_op_coeff, p_ext_data)
!
  TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
  TYPE(t_hydro_ocean_state), TARGET      :: p_os
  TYPE(t_operator_coeff),INTENT(IN)      :: p_op_coeff
  TYPE(t_external_data), TARGET          :: p_ext_data
  !
  !  local variables
  INTEGER  :: i_startidx_e, i_endidx_e
  INTEGER  :: je, jk, jb
  REAL(wp) :: gdt
  REAL(wp) :: z_grad_h(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
  TYPE(t_subset_range), POINTER :: edges_in_domain!, cells_in_domain
  CHARACTER(len=*), PARAMETER ::     &
    &      method_name='mo_oce_ab_timestepping_mimetic: calc_normal_velocity_ab_mimetic'
  TYPE(t_patch), POINTER :: patch
  !----------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')
  !-----------------------------------------------------------------------
  patch         => p_patch_3D%p_patch_2D(1)
  edges_in_domain => patch%edges%in_domain
  !cells_in_domain => patch%cells%in_domain
  !----------------------------------------------------------------------

  z_grad_h(1:nproma,1:patch%nblks_e) = 0.0_wp

  gdt=grav*dtime

  ! Step 1) Compute normal derivative of new surface height
  IF(.NOT.l_RIGID_LID.OR.iswm_oce == 1) THEN
    CALL grad_fd_norm_oce_2d_3D( p_os%p_prog(nnew(1))%h,&
      &                  patch,                     &
      &                  p_op_coeff%grad_coeff(:,1,:),&
      &                  z_grad_h(:,:))
  ENDIF


  ! Step 2) Calculate the new velocity from the predicted one and the new surface height
  IF (iswm_oce == 1) THEN ! shallow water case

    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
#ifdef __SX__
  !CDIR UNROLL=6
#endif
      DO jk = 1, n_zlev
        DO je = i_startidx_e, i_endidx_e
          IF(p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary)THEN
          p_os%p_prog(nnew(1))%vn(je,jk,jb) = (p_os%p_diag%vn_pred(je,jk,jb) &
                                            &- gdt*ab_beta*z_grad_h(je,jb))
          ENDIF
        END DO
      END DO
    END DO
    
  ELSE !real 3d case

    IF (.NOT.l_RIGID_LID) THEN

      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
#ifdef __SX__
  !CDIR UNROLL=6
#endif
        DO jk = 1, n_zlev
          DO je = i_startidx_e, i_endidx_e
            IF(p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary)THEN
            p_os%p_prog(nnew(1))%vn(je,jk,jb) = (p_os%p_diag%vn_pred(je,jk,jb) &
                                            &- gdt*ab_beta*z_grad_h(je,jb))
            ENDIF
          END DO
        END DO
      END DO

    ELSE
      CALL finish(method_name,"l_RIGID_LID case has a bug")
      !p_os%p_prog(nnew(1))%vn(:,:,:) = p_os%p_diag%vn_pred*v_base%wet_e(:,:,:)
    ENDIF
  ENDIF


      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
#ifdef __SX__
  !CDIR UNROLL=6
#endif
        DO jk = 1, n_zlev
          DO je = i_startidx_e, i_endidx_e
            IF(p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary)THEN
              p_os%p_diag%vn_time_weighted(je,jk,jb)      &
              &= ab_gam*p_os%p_prog(nnew(1))%vn(je,jk,jb) &
              &+ (1.0_wp -ab_gam)*p_os%p_prog(nold(1))%vn(je,jk,jb)
            ENDIF  
          END DO
        END DO
      END DO

  CALL sync_patch_array(SYNC_E, patch, p_os%p_prog(nnew(1))%vn)
  CALL sync_patch_array(SYNC_E, patch, p_os%p_diag%vn_time_weighted)

  ! slo: z_grad_h not out of sync, but sync error with global_max in dbg_print?
  !CALL sync_patch_array(SYNC_E, patch, z_grad_h)

  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=3  ! output print level (1-5, fix)
  CALL dbg_print('NormVel: vn_old'            ,p_os%p_prog(nold(1))%vn     ,str_module,idt_src)
  CALL dbg_print('NormVel: vn_pred'           ,p_os%p_diag%vn_pred         ,str_module,idt_src)
  IF (.NOT.l_rigid_lid) THEN                
    CALL dbg_print('NormVel: grad h-new'      ,z_grad_h                    ,str_module,idt_src)
  END IF
  CALL dbg_print('NormVel: vn_time_weighted'  ,p_os%p_diag%vn_time_weighted,str_module,idt_src)
  CALL dbg_print('NormVel: vn_change'         ,p_os%p_prog(nnew(1))%vn - &
    &                                          p_os%p_prog(nold(1))%vn     ,str_module,idt_src)
  idt_src=2  ! outputm print level (1-5, fix)
  CALL dbg_print('NormVel: vn_new'            ,p_os%p_prog(nnew(1))%vn     ,str_module,idt_src)
  !---------------------------------------------------------------------

  ! Update of scalar product quantities
  IF(l_STAGGERED_TIMESTEP)THEN
    !CALL height_related_quantities(p_patch_3D, p_os, p_ext_data)
    CALL calc_thickness(p_patch_3D, p_os, p_ext_data)
!   CALL calc_scalar_product_veloc_3D( patch,                &
!                                    & p_os%p_prog(nnew(1))%vn,&
!                                    & p_os%p_prog(nnew(1))%vn,&
!                                    & p_os%p_diag,            &
!                                    & p_op_coeff)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('NorVel: Staggered, kin'    ,p_os%p_diag%kin        ,str_module,idt_src)
    CALL dbg_print('NorVel: Staggered, ptp_vn' ,p_os%p_diag%ptp_vn     ,str_module,idt_src)
    !---------------------------------------------------------------------
  ENDIF
  !CALL message (TRIM(routine), 'end')

END SUBROUTINE calc_normal_velocity_ab_mimetic
  !-------------------------------------------------------------------------
  !
  !  
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
  SUBROUTINE calc_vert_velocity_mim_bottomup( p_patch_3D, p_os, p_diag,p_op_coeff,pw_c )

    TYPE(t_patch_3D), TARGET, INTENT(IN) :: p_patch_3D       ! patch on which computation is performed
    TYPE(t_hydro_ocean_state)         :: p_os
    TYPE(t_hydro_ocean_diag)          :: p_diag
    TYPE(t_operator_coeff),INTENT(IN) :: p_op_coeff
    REAL(wp),         INTENT(INOUT)   :: pw_c (nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%nblks_c) ! vertical velocity on cells
    !
    !
    ! Local variables
    INTEGER :: jc, jk, jb, je
    INTEGER :: z_dolic
    INTEGER :: i_startidx, i_endidx
    REAL(wp) :: div_depth_int(nproma,p_patch_3D%p_patch_2D(1)%nblks_c)
    !REAL(wp) :: z_vn_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: z_c(nproma,p_patch_3D%p_patch_2D(1)%nblks_c)
    REAL(wp) :: z_abort
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: p_patch
    
    !-----------------------------------------------------------------------  
    p_patch         => p_patch_3D%p_patch_2D(1)
    cells_in_domain => p_patch%cells%in_domain
    edges_in_domain => p_patch%edges%in_domain
    
    ! due to nag -nan compiler-option:
    pw_c(1:nproma,1:n_zlev+1,1:p_patch%nblks_c) = 0.0_wp
    div_depth_int(1:nproma,1:p_patch%nblks_c)   = 0.0_wp
    z_c          (1:nproma,1:p_patch%nblks_c)   = 0.0_wp
    !------------------------------------------------------------------
    ! Step 1) Calculate divergence of horizontal velocity at all levels
    !------------------------------------------------------------------
    
    ! !-------------------------------------------------------------------------------
    IF(l_EDGE_BASED )THEN
    ! !-------------------------------------------------------------------------------
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
        DO jk = 1, n_zlev
          DO je = i_startidx, i_endidx
            p_os%p_diag%mass_flx_e(je,jk,jb) = p_diag%vn_time_weighted(je,jk,jb)&
            &* p_patch_3D%p_patch_1D(1)%prism_thick_e(je,jk,jb)
          END DO
        END DO
      END DO
    
      CALL div_oce_3D( p_os%p_diag%mass_flx_e,    &
                     & p_patch,                   &
                     & p_op_coeff%div_coeff,      &
                     & p_os%p_diag%div_mass_flx_c,&
                     & subset_range=cells_in_domain)


      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
        DO jc = i_startidx, i_endidx
    
          z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
    
          !use bottom boundary condition for vertical velocity at bottom
          !of prism
          pw_c(jc,z_dolic+1,jb)=0.0_wp
          DO jk = z_dolic, 1, -1
            !IF(p_patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            pw_c(jc,jk,jb) = pw_c(jc,jk+1,jb) - p_os%p_diag%div_mass_flx_c(jc,jk,jb)
            !ENDIF
          END DO
        END DO
      END DO
    
    ! !-------------------------------------------------------------------------------
    ELSEIF(.NOT.l_EDGE_BASED)THEN
    ! !-------------------------------------------------------------------------------
    
      CALL map_edges2edges_viacell_3d_const_z( p_patch_3D, p_diag%vn_time_weighted, p_op_coeff, p_os%p_diag%mass_flx_e)
    
    
      CALL sync_patch_array(SYNC_E,p_patch,p_os%p_diag%mass_flx_e)
    
      CALL div_oce_3D( p_os%p_diag%mass_flx_e,      &
                     & p_patch,p_op_coeff%div_coeff,&
                     & p_os%p_diag%div_mass_flx_c,  &
                     & subset_range=cells_in_domain)
      CALL sync_patch_array(SYNC_C,p_patch,p_os%p_diag%div_mass_flx_c)
    
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
        DO jc = i_startidx, i_endidx
    
          z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
    
          !use bottom boundary condition for vertical velocity at bottom
          !of prism
          pw_c(jc,z_dolic+1,jb)=0.0_wp
          DO jk = z_dolic, 1, -1
            !IF(p_patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            pw_c(jc,jk,jb) = pw_c(jc,jk+1,jb) - p_os%p_diag%div_mass_flx_c(jc,jk,jb)
            !ENDIF
          END DO
        END DO
      END DO

    ENDIF  !  (l_EDGE_BASED)
    
    IF(l_RIGID_LID)THEN
      pw_c(:,1,:) = 0.0_wp
    ENDIF
    CALL sync_patch_array(SYNC_C,p_patch,pw_c)
    
    !  DO jk=1,n_zlev
    !  write(*,*)'vert veloc',jk,minval(pw_c(:,jk,:)),maxval(pw_c(:,jk,:))
    !  END DO
    !  DO jk=1,n_zlev
    !  write(*,*)'div-mass-flux',jk,minval(p_os%p_diag%div_mass_flx_c(:,jk,:)),&
    !  &maxval(p_os%p_diag%div_mass_flx_c(:,jk,:))
    !  END DO
    
    z_c(:,:) = ((p_os%p_prog(nnew(1))%h(:,:)-p_os%p_prog(nold(1))%h(:,:))/dtime - pw_c(:,1,:))*p_patch_3D%wet_c(:,1,:)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=1  ! output print level (1-5, fix)
    CALL dbg_print('CalcVertVelMimBU: d_h/dt-w',z_c                       ,str_module,idt_src)
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('CalcVertVelMimBU: pw_c =W' ,pw_c                      ,str_module,idt_src)
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('CalcVertVelMimBU: mass flx',p_os%p_diag%mass_flx_e,    str_module,idt_src)
    CALL dbg_print('CalcVertVelMimBU: div mass',p_os%p_diag%div_mass_flx_c,str_module,idt_src)
    !---------------------------------------------------------------------

    ! Abort if largest mismatch in surface elevation due to solution of ocean_gmres-solver is > 1mm/year
    !   criterion is 1mm/year * dtime = 3.17e-11 m/s * dtime
    z_abort = dhdtw_abort*dtime
    IF (MAXVAL(ABS(z_c(:,:))) > z_abort) THEN
      CALL message('mo_oce_ab_timestepping_mimetic:calc_vert_velocity_mim_bottomup', &
        &          'MISMATCH IN SURFACE EQUATION:')
      CALL message('mo_oce_ab_timestepping_mimetic:calc_vert_velocity_mim_bottomup', &
        &          'Elevation change does not match vertical velocity')
      WRITE(message_text,'(2(a,e20.12))') ' (h_new-h_old)/dtime - w = ', MAXVAL(ABS(z_c(:,:))), &
        &           ' z_abort=', z_abort
      CALL message ('mo_oce_ab_timestepping_mimetic:calc_vert_velocity_mim_bottomup', message_text)
      CALL finish(TRIM('mo_oce_ab_timestepping_mimetic:calc_vert_velocity_mim_bottomup'), &
        &            'MISMATCH in surface equation')
    ENDIF

  END SUBROUTINE calc_vert_velocity_mim_bottomup
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!  mpi parallelized, the result is NOT synced. Should be done in the calling method if required
  FUNCTION inverse_primal_flip_flop(p_patch, p_patch_3D, p_op_coeff, rhs_e, h_e) result(inv_flip_flop_e)
    !
    TYPE(t_patch), TARGET :: p_patch 
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    TYPE(t_operator_coeff),INTENT(IN)             :: p_op_coeff
    REAL(wp)      :: rhs_e(:,:,:)!(nproma,n_zlev,p_patch%nblks_e)
    REAL(wp)      :: h_e(:,:)  !(nproma,p_patch%nblks_e)
    REAL(wp)      :: inv_flip_flop_e(SIZE(rhs_e,1),SIZE(rhs_e,2),SIZE(rhs_e,3))
    !
    !LOCAL VARIABLES
    INTEGER,PARAMETER :: nmax_iter= 800 ! maximum number of iterations
    REAL(wp) :: zimpl_coeff = 1.0_wp    !COEFF has to be set appropriately !!!!
    REAL(wp) :: zimpl_prime_coeff
    INTEGER  :: n_iter                  ! number of iterations
    REAL(wp) :: tolerance               ! (relative or absolute) tolerance
    REAL(wp) :: z_residual(nmax_iter)
    LOGICAL  :: lmax_iter               ! true if reached m iterations
    !LOGICAL  :: lverbose = .TRUE.
    !CHARACTER(len=MAX_CHAR_LENGTH) :: string
    REAL(wp) :: rhstemp(nproma,p_patch%nblks_e)
    !REAL(wp), ALLOCATABLE :: inv_flip_flop_e2(:,:)!(nproma,p_patch%nblks_e)
    REAL(wp) :: z_e(nproma,n_zlev,p_patch%nblks_e)
    INTEGER :: jk
    !INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    !-----------------------------------------------------------------------

    tolerance                = 1.0e-12_wp  ! solver_tolerance
    inv_flip_flop_e(:,:,:)   = 0.0_wp
    zimpl_prime_coeff = (1.0_wp-zimpl_coeff)

    rhstemp(:,:)          = 0.0_wp

    DO jk=1, n_zlev
      rhstemp(:,:) = rhs_e(:,jk,:)&
      & -zimpl_coeff*lhs_primal_flip_flop(inv_flip_flop_e(:,jk,:), p_patch, p_patch_3D, p_op_coeff,jk,zimpl_coeff, h_e)

      If (maxval (ABS (rhstemp (:,:))) <= tolerance) THEN
        inv_flip_flop_e(:,jk,:) = lhs_primal_flip_flop(inv_flip_flop_e(:,jk,:), p_patch, p_patch_3D, p_op_coeff,jk,zimpl_coeff, h_e)
        print*, "Inv_flipflop ocean_gmres solved by initial guess!",&
          & jk,MAXVAL(rhstemp(:,:)), MINVAL(rhstemp(:,:)),MAXVAL(rhs_e(:,jk,:)), MINVAL(rhs_e(:,jk,:))
      ELSE
        inv_flip_flop_e(:,jk,:)= 0.0_wp!rhs_e(:,jk,:)
        !write(*,*)'RHS', maxvaL(rhs_e(:,jk,:)),minvaL(rhs_e(:,jk,:))

        CALL ocean_gmres( inv_flip_flop_e(:,jk,:), &  ! arg 1 of lhs. x input is the first guess.
        &        lhs_primal_flip_flop,      &  ! function calculating l.h.s.
        &        h_e,                       &  ! edge thickness for LHS
        &        jk,                        &
        &        p_patch, p_patch_3D,       &  !arg 3 of lhs 
        &        zimpl_coeff,               &  !arg 4 of lhs
        &        p_op_coeff,                &
        &        rhs_e(:,jk,:),             &  ! right hand side as input
        &        tolerance,                 &  ! relative tolerance
        &        .FALSE.,                   &  ! NOT absolute tolerance
        &        nmax_iter,                 &  ! max. # of iterations to do
        &        lmax_iter,                 &  ! out: .true. = not converged
        &        n_iter,                    &  ! out: # of iterations done
        &        z_residual)                  ! inout: the residual (array)  

        rhstemp(:,:) = rhs_e(:,jk,:)-lhs_primal_flip_flop(inv_flip_flop_e(:,jk,:),p_patch, p_patch_3D,p_op_coeff,&
          &            jk,zimpl_coeff,h_e)
        !WRITE(*,*)'max/min residual of inverse primal-flip-flop:',&
        !  &        jk, maxval(rhstemp),minval(rhstemp) 
        idt_src=2  ! output print level (1-5, fix)
        CALL dbg_print('residual of inv_flip_flop'   ,rhstemp ,str_module,idt_src)
        !write(*,*)'sol', maxvaL(inv_flip_flop_e(:,jk,:)),minvaL(inv_flip_flop_e(:,jk,:))
        z_e(:,jk,:)=rhstemp(:,:)
        If (maxval (ABS (rhstemp (:,:))) >= tolerance) lmax_iter = .true.
          idt_src=1
          IF (idbg_mxmn >= idt_src) THEN
            IF (lmax_iter) THEN
              WRITE (0, '(1x,a, I4.2, 1x, a,E8.2,1x, a,E8.2,1x, E8.2, 1x, a)') &
              &'Inv_flipflop ocean_gmres #Iter', n_iter, 'Tol ',tolerance, 'Res ',&
              &  ABS(z_residual(n_iter)),MAXVAL (ABS(rhstemp(:,:))), 'ocean_gmres PROBLEM!!!!!!!!!!!!'
            ELSE
              WRITE (0, '(1x,a, I4.2, 1x, a,E8.2,1x, a,E8.2,1x, E8.2)') &
              &'Inv_flipflop ocean_gmres #Iter', n_iter, 'Tol ',tolerance, 'Res ',&
              &  ABS(z_residual(n_iter)),MAXVAL (ABS(rhstemp(:,:)))
            ENDIF
          ENDIF
       END IF
    END DO

  END FUNCTION inverse_primal_flip_flop
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !!  mpi parallelized LL, results is valid only in in_domain edges
  FUNCTION lhs_primal_flip_flop( x, p_patch, p_patch_3D, p_op_coeff,jk,coeff, h_e) RESULT(llhs)
    !
    TYPE(t_patch), TARGET, INTENT(in)             :: p_patch
    TYPE(t_patch_3D ),TARGET, INTENT(IN)          :: p_patch_3D
    REAL(wp),INTENT(inout)                        :: x(:,:)
    TYPE(t_operator_coeff),INTENT(IN)             :: p_op_coeff
    INTEGER ,INTENT(in)                           :: jk
    REAL(wp),INTENT(in)                           :: coeff
    REAL(wp),OPTIONAL,INTENT(in)                  :: h_e(SIZE(x,1), SIZE(x,2))!(:,:)
    REAL(wp)                                      :: llhs(SIZE(x,1), SIZE(x,2))

    !local variables
    REAL(wp) :: z_x_out(SIZE(x,1), 1,SIZE(x,2))!(nproma,p_patch%nblks_e)
    REAL(wp) :: z_e(SIZE(x,1), 1,SIZE(x,2))!(nproma,p_patch%nblks_e)
    !TYPE(t_cartesian_coordinates) :: z_vn_cc(nproma,p_patch%nblks_c)
    !TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    !edges_in_domain => p_patch%edges%in_domain

    z_x_out(:,1,:) = x 

    CALL sync_patch_array(SYNC_E, p_patch, x)

    CALL map_edges2edges_viacell_3D( p_patch_3D,    &
                                   & z_x_out(:,1,:),&
                                   & p_op_coeff,    &
                                   & z_e(:,1,:),    &
                                   !& p_patch_3D%p_patch_1D(n_dom)%prism_thick_c(:,1,:),&
                                   & level=1)


    llhs(:,:) = coeff * z_e(:,1,:)
    !write(*,*)'max/min in', maxval(x(:,:)),minval(x(:,:)) 
    !rite(*,*)'max/min LHS', maxval(llhs(:,:)),minval(llhs(:,:)) 

  END FUNCTION lhs_primal_flip_flop
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
END MODULE mo_oce_ab_timestepping_mimetic
