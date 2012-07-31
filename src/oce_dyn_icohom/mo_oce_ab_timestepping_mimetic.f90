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
USE mo_parallel_config,           ONLY: nproma
USE mo_math_utilities,            ONLY: t_cartesian_coordinates
USE mo_sync,                      ONLY: sync_e, sync_c, sync_patch_array
USE mo_mpi,                       ONLY: my_process_is_mpi_parallel
USE mo_impl_constants,            ONLY: sea_boundary,                                     &
! &                                     min_rlcell, min_rledge, min_rlcell,               &
  &                                     max_char_length, MIN_DOLIC
USE mo_dbg_nml,                   ONLY: idbg_mxmn
USE mo_ocean_nml,                 ONLY: n_zlev, solver_tolerance, l_inverse_flip_flop,    &
  &                                     ab_const, ab_beta, ab_gam, iswm_oce,              &
  &                                     expl_vertical_velocity_diff, iforc_oce, EOS_TYPE, &
  &                                     no_tracer, l_RIGID_LID, l_edge_based,             &
  &                                     FLUX_CALCULATION_HORZ, MIMETIC,CENTRAL,UPWIND
USE mo_run_config,                ONLY: dtime, ltimer
USE mo_timer,                     ONLY: timer_start, timer_stop, timer_ab_expl,           &
  &                                     timer_ab_rhs4sfc
USE mo_dynamics_config,           ONLY: nold, nnew
USE mo_physical_constants,        ONLY: grav
USE mo_oce_state,                 ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, v_base,  &
  &                                     set_lateral_boundary_values, is_initial_timestep
USE mo_model_domain,              ONLY: t_patch
USE mo_ext_data_types,            ONLY: t_external_data
USE mo_gmres,                     ONLY: gmres
USE mo_exception,                 ONLY: message, finish!, message_text
USE mo_util_dbg_prnt,             ONLY: dbg_print
USE mo_oce_boundcond,             ONLY: bot_bound_cond_horz_veloc, top_bound_cond_horz_veloc
USE mo_oce_thermodyn,             ONLY: calc_density, calc_internal_press
USE mo_oce_physics,               ONLY: t_ho_params
USE mo_sea_ice_types,             ONLY: t_sfc_flx
USE mo_scalar_product,            ONLY: map_cell2edges_3D,&
  &                                     map_edges2cell_3D, calc_scalar_product_veloc_3D
USE mo_oce_math_operators,        ONLY: grad_fd_norm_oce_2d, &
  &                                     div_oce_3D, grad_fd_norm_oce_3D,&
  &                                     grad_fd_norm_oce_2d_3D,  &
  &                                     height_related_quantities
USE mo_oce_veloc_advection,       ONLY: veloc_adv_horz_mimetic, veloc_adv_vert_mimetic
USE mo_intp_data_strc,            ONLY: t_int_state 
USE mo_oce_diffusion,             ONLY: velocity_diffusion,& 
  &                                     velocity_diffusion_vert_mimetic,  &
  &                                     veloc_diffusion_vert_impl_hom
USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
USE mo_grid_subset,               ONLY: t_subset_range, get_index_range

IMPLICIT NONE

PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'
!
! PUBLIC INTERFACE
!
PUBLIC :: solve_free_sfc_ab_mimetic
PUBLIC :: calc_normal_velocity_ab_mimetic
PUBLIC :: calc_vert_velocity_mimetic
PUBLIC :: calc_vert_velocity_mim_topdown
!
! Private implemenation
!
PRIVATE :: fill_rhs4surface_eq_ab
PRIVATE :: calculate_explicit_term_ab   ! calc_velocity_predictor
PRIVATE :: lhs_surface_height_ab_mim
PRIVATE :: inverse_primal_flip_flop
PRIVATE :: update_column_thickness
INTEGER, PARAMETER :: top=1
LOGICAL, PARAMETER :: l_forc_freshw = .FALSE.
CHARACTER(len=12)  :: str_module = 'oceSTEPmimet'  ! Output of module for 1 line debug
INTEGER            :: idt_src    = 1               ! Level of detail for 1 line debug

! TRUE=staggering between thermodynamic and dynamic part, offset of half timestep
! between dynamic and thermodynamic variables thermodynamic and dnamic variables are colocated in time
LOGICAL, PUBLIC,PARAMETER :: l_STAGGERED_TIMESTEP = .FALSE. 

CONTAINS
!-------------------------------------------------------------------------


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
SUBROUTINE solve_free_sfc_ab_mimetic(p_patch, p_os, p_ext_data, p_sfc_flx, &
    &                                  p_phys_param, timestep, p_op_coeff,p_int)
  !
  TYPE(t_patch), TARGET, INTENT(in)             :: p_patch
  TYPE(t_hydro_ocean_state), TARGET             :: p_os
  TYPE(t_external_data), TARGET, INTENT(in)     :: p_ext_data
  TYPE(t_sfc_flx), INTENT(INOUT)                :: p_sfc_flx
  TYPE (t_ho_params)                            :: p_phys_param
  INTEGER                                       :: timestep
  TYPE(t_int_state),TARGET,INTENT(IN)           :: p_int  
  TYPE(t_operator_coeff)                        :: p_op_coeff
  !
  !Local variables
  !
  !REAL(wp),ALLOCATABLE :: rhs(:,:)
  !REAL(wp) :: aveh, totarea
  !REAL(wp), POINTER :: p_height_area(:,:)
  !GMRS
  INTEGER,PARAMETER :: nmax_iter   = 200      ! maximum number of iterations
  REAL(wp) :: tolerance =0.0_wp  ! 1.0e-6_wp  ! (relative or absolute) tolerance
  INTEGER  :: n_iter                          ! number of iterations 

  REAL(wp) :: z_h_c(nproma,p_patch%nblks_c)
  REAL(wp) :: z_h_e(nproma,p_patch%nblks_e)

  REAL(wp) :: z_implcoeff
  REAL(wp) :: zresidual(nmax_iter)    ! norms of the residual (convergence history); 
  ! an argument of dimension at least m is required
  LOGICAL  :: l_maxiter                 ! true if reached m iterations
  !LOGICAL  :: lverbose         = .TRUE.
  !LOGICAL  :: l_first_timestep = .FALSE.
  CHARACTER(len=max_char_length) :: string
  !INTEGER  :: rl_start_c, rl_end_c, i_startblk_c, i_endblk_c, jc,jb, i_startidx_c, i_endidx_c
  !INTEGER  :: rl_start_e, rl_end_e, i_startblk_e, i_endblk_e, je,jb, i_startidx_e, i_endidx_e
  !CHARACTER(len=max_char_length), PARAMETER :: &
  !       & routine = ('mo_oce_ab_timestepping_mimetic:solve_free_sfc_ab_mimetic')
  !-------------------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')
  tolerance = solver_tolerance
  z_h_c = 0.0_wp
  z_h_e = 0.0_wp

  !CALL update_column_thickness( p_patch, p_os, p_os%p_diag, p_op_coeff, timestep )

  CALL sync_patch_array(sync_c, p_patch, p_os%p_prog(nold(1))%h)
  CALL sync_patch_array(sync_c, p_patch, p_os%p_diag%thick_c)
  CALL sync_patch_array(sync_e, p_patch, p_os%p_prog(nold(1))%vn)

  IF (is_initial_timestep(timestep) ) THEN

    CALL height_related_quantities(p_patch, p_os, p_ext_data)

    !This is required in top boundary condition for
    !vertical velocity: the time derivative of the surface height
    !is used there and needs special treatment in the first timestep.
    !see sbr top_bound_cond_vert_veloc in mo_ho_boundcond
    p_os%p_prog(nnew(1))%h=p_os%p_prog(nold(1))%h

    !Call set_lateral_boundary_values(p_patch, p_os%p_prog(nold(1))%vn)


    IF (l_STAGGERED_TIMESTEP ) &
     & CALL calc_scalar_product_veloc_3D( p_patch,&
                                      & p_os%p_prog(nold(1))%vn,&
                                      & p_os%p_prog(nold(1))%vn,&
                                      & p_os%p_diag%h_e,        &
                                      & p_os%p_diag,            &
                                      & p_op_coeff)
  ENDIF


  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=2  ! output print level (1-5, fix)
  CALL dbg_print('on entry: h-old'                ,p_os%p_prog(nold(1))%h ,str_module,idt_src)
  CALL dbg_print('on entry: h-new'                ,p_os%p_prog(nnew(1))%h ,str_module,idt_src)
  CALL dbg_print('on entry: vn-old'               ,p_os%p_prog(nold(1))%vn,str_module,idt_src)
  idt_src=1  ! output print level (1-5, fix)
  CALL dbg_print('on entry: vn-new'               ,p_os%p_prog(nnew(1))%vn,str_module,idt_src)
  !---------------------------------------------------------------------

  ! abort condition for elevation and vn:
  IF ( (maxval(p_os%p_prog(nnew(1))%h)  >  1.e20_wp) .or. &
    &  (minval(p_os%p_prog(nnew(1))%h)  < -1.e20_wp) .or. &
    &  (maxval(p_os%p_prog(nold(1))%vn) >  1.e20_wp) .or. &
    &  (minval(p_os%p_prog(nnew(1))%vn) < -1.e20_wp) ) THEN
    CALL message('Solve free surface AB mimetic: ',' INSTABLE VN or H - stop now ')
    CALL finish ('Solve free surface AB mimetic: ',' INSTABLE VN or H !!')
  END IF

  ! Apply windstress 
  CALL top_bound_cond_horz_veloc(p_patch, p_os, p_op_coeff, p_sfc_flx,     &
    &                            p_os%p_aux%bc_top_u, p_os%p_aux%bc_top_v, &
    &                            p_os%p_aux%bc_top_veloc_cc)


  ! Apply bot boundary condition for horizontal velocity
  CALL bot_bound_cond_horz_veloc(p_patch, p_os, p_phys_param, p_op_coeff%div_coeff)

  ! Calculate explicit terms of Adams-Bashforth timestepping
  IF (ltimer) CALL timer_start(timer_ab_expl)
  CALL calculate_explicit_term_ab(p_patch, p_os, p_phys_param, p_int,&
    &                             is_initial_timestep(timestep), p_op_coeff)
  IF (ltimer) CALL timer_stop(timer_ab_expl)

  IF(.NOT.l_RIGID_LID)THEN

    ! Calculate RHS of surface equation
    IF (ltimer) CALL timer_start(timer_ab_rhs4sfc)
    CALL fill_rhs4surface_eq_ab( p_patch, p_os, p_sfc_flx, p_op_coeff)
    IF (ltimer) CALL timer_stop(timer_ab_rhs4sfc)

    ! Solve surface equation with GMRES solver
    z_h_c = 0.0_wp !potentially better choice: p_os%p_prog(nold(1))%h

    !The lhs needs different thicknesses at edges for 3D and SWE (including bathymetry)
    IF( iswm_oce == 1 ) THEN  
      z_h_e = p_os%p_diag%h_e! p_os%p_diag%thick_e
      ELSEIF( iswm_oce /= 1 ) THEN 
      z_h_e = p_os%p_diag%h_e         ! #slo# 2011-02-21 bugfix (for mimetic only)
    ENDIF 

  CALL sync_patch_array(SYNC_E, p_patch, z_h_e)
  CALL sync_patch_array(SYNC_C, p_patch, p_os%p_diag%thick_c)
  CALL sync_patch_array(SYNC_C, p_patch, p_os%p_prog(nold(1))%h)


  CALL gmres( z_h_c(:,:),             &  ! arg 1 of lhs. x input is the first guess.
      &        lhs_surface_height_ab_mim, &  ! function calculating l.h.s.
      &        p_os%p_diag%thick_e,       &  ! z_h_e, &  !arg 5 of lhs  !not used
      &        p_os%p_diag%thick_c,       &  ! p_os%p_diag%thick_c, & 
                                             ! arg 6 of lhs p_os%p_prog(nold(1))%h,
                                             ! p_os%p_diag%cons_thick_c(:,1,:),&
      &        p_os%p_prog(nold(1))%h,    &  !arg 2 of lhs !not used
      &        p_patch,                   &  !arg 3 of lhs 
      &        z_implcoeff,               &  !arg 4 of lhs
      &        p_op_coeff,                &
      &        p_os%p_aux%p_rhs_sfc_eq,   &  ! right hand side as input
      &        tolerance,                 &  ! relative tolerance
      &        .FALSE.,                   &  ! NOT absolute tolerance
      &        nmax_iter,                 &  ! max. # of iterations to do
      &        l_maxiter,                 &  ! out: .true. = not converged
      &        n_iter,                    &  ! out: # of iterations done
      &        zresidual )                   ! inout: the residual (array)  

    IF (l_maxiter) THEN
      CALL finish('GMRES solver surface equation: ','NOT YET CONVERGED !!')
    ELSE
      ! output print level idt_src used for GMRES output with call message:
      IF(n_iter==0)n_iter=1
      idt_src=0
      IF (idbg_mxmn >= idt_src) THEN
        WRITE(string,'(a,i4,a,e28.20)') &
          'iteration =', n_iter,', residual =', ABS(zresidual(n_iter))
        CALL message('GMRES surface height',TRIM(string))
      ENDIF
    ENDIF 

    p_os%p_prog(nnew(1))%h = z_h_c
    z_h_c = lhs_surface_height_ab_mim( p_os%p_prog(nnew(1))%h, &
      & p_os%p_prog(nold(1))%h, &
      & p_patch,                &
      & z_implcoeff,            &
      & z_h_e,                  &
      & p_os%p_diag%thick_c,    &
      & p_op_coeff)             &
      & -p_os%p_aux%p_rhs_sfc_eq

    CALL sync_patch_array(SYNC_C, p_patch, p_os%p_prog(nnew(1))%h)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('SolvSfc: residual h-res'    ,z_h_c                  ,str_module,idt_src)
    idt_src=1  ! output print level (1-5, fix)
    CALL dbg_print('SolvSfc: after GMRES: h-new',p_os%p_prog(nnew(1))%h ,str_module,idt_src)
    !---------------------------------------------------------------------

  ENDIF  ! l_rigid_lid

END SUBROUTINE solve_free_sfc_ab_mimetic
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
SUBROUTINE calculate_explicit_term_ab( p_patch, p_os, p_phys_param,&
                                     & p_int, l_initial_timestep, p_op_coeff)

  TYPE(t_patch), TARGET, INTENT(in)             :: p_patch
  TYPE(t_hydro_ocean_state), TARGET             :: p_os
  TYPE (t_ho_params)                            :: p_phys_param
  TYPE(t_int_state),TARGET,INTENT(IN), OPTIONAL :: p_int
  LOGICAL,INTENT(IN)                            :: l_initial_timestep
  TYPE(t_operator_coeff)                        :: p_op_coeff
  !
  !local variables
  !
  INTEGER  :: je, jk, jb
  REAL(wp) :: gdt
  REAL(wp) :: z_gradh_e(nproma,1,p_patch%nblks_e)
  REAL(wp) :: z_e(nproma,n_zlev,p_patch%nblks_e)
  !INTEGER  :: rl_start_c, rl_end_c, i_startblk_c, i_endblk_c,i_startidx_c, i_endidx_c
  TYPE(t_subset_range), POINTER :: edges_in_domain, all_edges
  INTEGER  :: i_startidx_e, i_endidx_e

  CHARACTER(len=max_char_length), PARAMETER :: &
    &       routine = ('mo_oce_ab_timestepping_mimetic:calculate_explicit_term_ab')
  !-----------------------------------------------------------------------  
  !CALL message (TRIM(routine), 'start')        
  edges_in_domain => p_patch%edges%in_domain
  all_edges => p_patch%edges%all

  z_gradh_e(:,:,:) = 0.0_wp
  gdt              = grav*dtime

  !---------------------------------------------------------------------
  ! STEP 1: calculate gradient of surface height at previous timestep
  !---------------------------------------------------------------------

  !IF ( iswm_oce == 1 ) THEN
  ! LL: already synced
  !CALL sync_patch_array(sync_c, p_patch, p_os%p_prog(nold(1))%h)

  CALL grad_fd_norm_oce_2d_3D( p_os%p_prog(nold(1))%h, &
         &                  p_patch,                &
         &                  p_op_coeff%grad_coeff,  &
         &                  z_gradh_e(:,1,:))
  CALL sync_patch_array(sync_e, p_patch, z_gradh_e(:,1,:))

  !---------------------------------------------------------------------
  ! STEP 2: horizontal advection
  !---------------------------------------------------------------------

  IF(l_initial_timestep)THEN
    CALL veloc_adv_horz_mimetic( p_patch,            &
           &             p_os%p_prog(nold(1))%vn,    &
           &             p_os%p_prog(nold(1))%vn,    &
           &             p_os%p_diag,                &
           &             p_os%p_diag%veloc_adv_horz, &
           &             p_op_coeff, p_int )
  ELSE
    CALL veloc_adv_horz_mimetic( p_patch,            &
           &             p_os%p_prog(nold(1))%vn,    &
           &             p_os%p_prog(nnew(1))%vn,    &
           &             p_os%p_diag,                &
           &             p_os%p_diag%veloc_adv_horz, &
           &             p_op_coeff, p_int )
  ENDIF

  !---------------------------------------------------------------------
  ! STEP 3: compute 3D contributions: gradient of hydrostatic pressure and vertical velocity advection
  !---------------------------------------------------------------------

  IF ( iswm_oce /= 1 ) THEN
    ! calculate density from EOS using temperature and salinity at timelevel n
     CALL calc_density( p_patch,                                 &
            &           p_os%p_prog(nold(1))%tracer(:,:,:,1:no_tracer),  &
            &           p_os%p_diag%rho(:,:,:) )

    ! calculate hydrostatic pressure from density at timelevel nc
     CALL calc_internal_press( p_patch,               &
            &                  p_os%p_diag%rho,       &
            &                  p_os%p_prog(nold(1))%h,&
            &                  p_os%p_diag%press_hyd)

    ! calculate gradient of hydrostatic pressure in 3D
    CALL grad_fd_norm_oce_3D( p_os%p_diag%press_hyd,  &
           &                  p_patch,                &
           &                  p_op_coeff%grad_coeff,  &
           &                  p_os%p_diag%press_grad)
    CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%press_grad)

    ! calculate vertical velocity advection
    CALL veloc_adv_vert_mimetic( p_patch,          &
         &             p_os%p_diag,p_op_coeff,     &
         &             p_os%p_diag%veloc_adv_vert )

    ! calculate vertical velocity diffusion
    !   For the alternative choice "expl_vertical_velocity_diff==1" see couples of
    !   lines below below
    IF (expl_vertical_velocity_diff==0) THEN
        CALL velocity_diffusion_vert_mimetic( p_patch,        &
        &                             p_os%p_diag,            &
        &                             p_os%p_aux,p_op_coeff,  &
        &                             p_os%p_prog(nold(1))%h, &
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
  CALL dbg_print('density'                   ,p_os%p_diag%rho           ,str_module,idt_src)
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

  CALL velocity_diffusion(p_patch, &
                        & p_os%p_prog(nold(1))%vn, &
                        & p_phys_param,            &
                        & p_os%p_diag,p_op_coeff,  &
                        & p_os%p_diag%laplacian_horz)

  CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%laplacian_horz)

  IF (L_INVERSE_FLIP_FLOP) THEN

    IF ( iswm_oce /= 1 ) THEN
      z_e = inverse_primal_flip_flop(p_patch, p_os%p_diag%veloc_adv_horz, p_os%p_diag%h_e)
    ELSE
      z_e = inverse_primal_flip_flop(p_patch, p_os%p_diag%veloc_adv_horz, p_os%p_diag%thick_e)
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

    CALL sync_patch_array(SYNC_E, p_patch, p_os%p_aux%g_n)

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
    p_os%p_aux%g_nimd(1:nproma,1:n_zlev,1:p_patch%nblks_e) = p_os%p_aux%g_n(1:nproma,1:n_zlev,1:p_patch%nblks_e)
  ELSE
    p_os%p_aux%g_nimd(1:nproma,1:n_zlev,1:p_patch%nblks_e) = (1.5_wp+AB_const)* p_os%p_aux%g_n(1:nproma,1:n_zlev,1:p_patch%nblks_e)&
    &                        - (0.5_wp+AB_const)* p_os%p_aux%g_nm1(1:nproma,1:n_zlev,1:p_patch%nblks_e)
  ENDIF

  IF ( iswm_oce /= 1) THEN

    IF(.NOT.l_RIGID_LID)THEN

      IF(l_STAGGERED_TIMESTEP)THEN
        DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
          DO jk = 1, n_zlev
            DO je = i_startidx_e, i_endidx_e

              IF(v_base%dolic_e(je,jb)>=MIN_DOLIC)THEN
              !IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN

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
              IF(v_base%dolic_e(je,jb)>=MIN_DOLIC)THEN
              !IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
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

              IF(v_base%dolic_e(je,jb)>=MIN_DOLIC)THEN
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

              IF(v_base%dolic_e(je,jb)>=MIN_DOLIC)THEN
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

          IF(v_base%dolic_e(je,jb)>=MIN_DOLIC)THEN
            p_os%p_diag%vn_pred(je,1,jb) =  p_os%p_diag%vn_pred(je,1,jb)      &
            &                            + dtime*p_os%p_aux%bc_top_vn(je,jb)/v_base%del_zlev_m(1)

            p_os%p_diag%vn_pred(je,v_base%dolic_e(je,jb),jb)    &
            & = p_os%p_diag%vn_pred(je,v_base%dolic_e(je,jb),jb)&
            & - dtime*p_os%p_aux%bc_bot_vn(je,jb)               &
            &/v_base%del_zlev_m(v_base%dolic_e(je,jb))
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
          !IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
            p_os%p_diag%vn_pred(je,jk,jb) = (p_os%p_prog(nold(1))%vn(je,jk,jb)        &
             &                            + dtime*(p_os%p_aux%g_nimd(je,jk,jb)        &
             &                            - (1.0_wp-ab_beta)*grav*z_gradh_e(je,1,jb)))&
             &                            *v_base%wet_e(je,jk,jb)
        END DO
      END DO
    END DO

    !In case of Shallow-water with forcing and or damping
    IF ( iforc_oce/=10) THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
        DO jk = 1, n_zlev
          DO je = i_startidx_e, i_endidx_e
            !IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
              p_os%p_diag%vn_pred(je,jk,jb) = (p_os%p_diag%vn_pred(je,jk,jb) &
               &                            + p_os%p_aux%bc_top_vn(je,jb)    &
               &                            - p_os%p_aux%bc_bot_vn(je,jb))   &
               &                            *v_base%wet_e(je,jk,jb)
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
     CALL veloc_diffusion_vert_impl_hom( p_patch,                &
                                     & p_os%p_diag%vn_pred,      &
                                     & p_os%p_diag%h_e,          &
                                     & p_phys_param%A_veloc_v,   &
                                     & p_os%p_diag%vn_impl_vert_diff)

    IF(l_RIGID_LID)THEN
      p_os%p_diag%vn_pred(1:nproma,1:n_zlev,1:p_patch%nblks_e) = p_os%p_diag%vn_impl_vert_diff(1:nproma,1:n_zlev,1:p_patch%nblks_e)
    ENDIF

  ENDIF

  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=3  ! output print level (1-5, fix)
  CALL dbg_print('vn(nold)'                  ,p_os%p_prog(nold(1))%vn       ,str_module,idt_src)
  CALL dbg_print('VelocAdvHorizontal'        ,p_os%p_diag%veloc_adv_horz    ,str_module,idt_src)
  CALL dbg_print('VelocLaPlac horizontal'    ,p_os%p_diag%laplacian_horz    ,str_module,idt_src)
  CALL dbg_print('VelocLaPlac vertical'      ,p_os%p_diag%laplacian_vert    ,str_module,idt_src)
  IF (expl_vertical_velocity_diff == 1 .AND. iswm_oce /= 1) &
& CALL dbg_print('ImplVelocDiff vertical'    ,p_os%p_diag%vn_impl_vert_diff ,str_module,idt_src)
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
SUBROUTINE fill_rhs4surface_eq_ab( p_patch, p_os, p_sfc_flx, p_op_coeff)
!
! Patch on which computation is performed
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
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
REAL(wp) :: z_e(nproma,p_patch%nblks_e)
REAL(wp) :: z_e2(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: div_z_depth_int_c(nproma,p_patch%nblks_e)
REAL(wp) :: gdt2, delta_z
REAL(wp) :: div_z_c_2D(nproma,p_patch%nblks_c)
REAL(wp) :: div_z_c(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_vn_ab(nproma,n_zlev,p_patch%nblks_e)
TYPE(t_cartesian_coordinates) :: z_u_pred_cc(nproma,n_zlev,p_patch%nblks_c)
TYPE(t_cartesian_coordinates) :: z_u_cc1(nproma,n_zlev,p_patch%nblks_c)
TYPE(t_cartesian_coordinates) :: z_u_cc2(nproma,n_zlev,p_patch%nblks_c)
TYPE(t_cartesian_coordinates) :: z_u_pred_depth_int_cc(nproma,p_patch%nblks_c)
TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain, all_edges
!REAL(wp) :: thick
!CHARACTER(len=max_char_length), PARAMETER :: &
!       & routine = ('mo_oce_ab_timestepping_mimetic:fill_rhs4surface_eq_ab')
!-------------------------------------------------------------------------------
!CALL message (TRIM(routine), 'start')        
  all_cells       => p_patch%cells%all
  cells_in_domain => p_patch%cells%in_domain
  all_edges       => p_patch%edges%all

  gdt2 = grav*(dtime)**2

  z_vn_ab(1:nproma,1:n_zlev,1:p_patch%nblks_e)  = 0.0_wp
  z_e(1:nproma,1:p_patch%nblks_e)               = 0.0_wp
  z_e2(1:nproma,1:n_zlev,1:p_patch%nblks_e)     = 0.0_wp
  div_z_depth_int_c(1:nproma,1:p_patch%nblks_c) = 0.0_wp
  div_z_c_2D(1:nproma,1:p_patch%nblks_c)        = 0.0_wp
  div_z_c(1:nproma,1:n_zlev,1:p_patch%nblks_c)  = 0.0_wp
  
  z_u_pred_depth_int_cc(1:nproma,1:p_patch%nblks_c)%x(1) = 0.0_wp
  z_u_pred_depth_int_cc(1:nproma,1:p_patch%nblks_c)%x(2) = 0.0_wp
  z_u_pred_depth_int_cc(1:nproma,1:p_patch%nblks_c)%x(3) = 0.0_wp
  z_u_pred_cc(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(1) = 0.0_wp
  z_u_pred_cc(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(2) = 0.0_wp
  z_u_pred_cc(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(3) = 0.0_wp

  ! LL: this should not be required
  CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%vn_pred)
  CALL sync_patch_array(SYNC_E, p_patch, p_os%p_prog(nold(1))%vn)
  CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%vn_impl_vert_diff)

  IF(iswm_oce == 1)THEN
    z_vn_ab(1:nproma,1:n_zlev,1:p_patch%nblks_e) = ab_gam*p_os%p_diag%vn_pred(1:nproma,1:n_zlev,1:p_patch%nblks_e) &
                                          &+ (1.0_wp -ab_gam)* p_os%p_prog(nold(1))%vn(1:nproma,1:n_zlev,1:p_patch%nblks_e)
  ELSEIF(iswm_oce /= 1)THEN
    IF(expl_vertical_velocity_diff==1)THEN
      z_vn_ab(1:nproma,1:n_zlev,1:p_patch%nblks_e) = ab_gam*p_os%p_diag%vn_impl_vert_diff(1:nproma,1:n_zlev,1:p_patch%nblks_e) &
                                            &+ (1.0_wp -ab_gam)* p_os%p_prog(nold(1))%vn(1:nproma,1:n_zlev,1:p_patch%nblks_e)
    ELSEIF(expl_vertical_velocity_diff==0)THEN
      z_vn_ab(1:nproma,1:n_zlev,1:p_patch%nblks_e) = ab_gam*p_os%p_diag%vn_pred(1:nproma,1:n_zlev,1:p_patch%nblks_e) &
                                            &+ (1.0_wp -ab_gam)* p_os%p_prog(nold(1))%vn(1:nproma,1:n_zlev,1:p_patch%nblks_e)
    ENDIF
  ENDIF

IF (l_edge_based) THEN

!  IF( iswm_oce /= 1 ) THEN !the 3D case
!   !calculate depth-integrated velocity 
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      DO je = i_startidx_e, i_endidx_e
        i_dolic_e = v_base%dolic_e(je,jb)
        DO jk=1,i_dolic_e !1,i_dolic_c
          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
            z_e2(je,jk,jb)= z_vn_ab(je,jk,jb)*v_base%del_zlev_m(jk)
            IF(jk==1)THEN
              z_e2(je,jk,jb)=z_vn_ab(je,jk,jb)&
              &*(v_base%del_zlev_m(jk)+p_os%p_diag%h_e(je,jb))
            ENDIF
          ENDIF
        END DO
      ENDDO
    END DO
   !ENDIF
  CALL div_oce_3d(z_e2, p_patch,p_op_coeff%div_coeff, div_z_c,&
                & subset_range=cells_in_domain )
   DO jb = all_cells%start_block, all_cells%end_block
     CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
     DO jc = i_startidx_c, i_endidx_c
       i_dolic_c = v_base%dolic_c(jc,jb)
       DO jk=1,i_dolic_c !1,i_dolic_c
         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
            div_z_depth_int_c(jc,jb)&
            &=div_z_depth_int_c(jc,jb) +div_z_c(jc,jk,jb)
           ENDIF
        END DO
    END DO
  END DO
!    CALL div_oce_3d(z_vn_ab, p_patch,p_op_coeff%div_coeff, div_z_c,&
!                  & subset_range=cells_in_domain )
!    DO jb = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
!      DO jc = i_startidx_c, i_endidx_c
!        i_dolic_c = v_base%dolic_c(jc,jb)
!        DO jk=1,i_dolic_c !1,i_dolic_c
!          IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
!            delta_z=v_base%del_zlev_m(jk)
!            IF(jk==1)delta_z=v_base%del_zlev_m(jk)+p_os%p_prog(nold(1))%h(jc,jb)
!             div_z_depth_int_c(jc,jb)&
!             &=div_z_depth_int_c(jc,jb) +div_z_c(jc,jk,jb)*delta_z
!            ENDIF
!         END DO
!     END DO
!   END DO

  IF(l_forc_freshw)THEN

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
!         p_os%p_aux%p_rhs_sfc_eq(jc,jb) = ((p_os%p_prog(nold(1))%h(jc,jb)&
!                                        & - dtime*(div_z_c(jc,jb)        &
!                                        & + p_sfc_flx%forc_tracer(jc,jb,2) ))/gdt2)&
!                                        &  *v_base%wet_c(jc,1,jb) !last idx=2 for freshwater
          p_os%p_aux%p_rhs_sfc_eq(jc,jb) = ((p_os%p_prog(nold(1))%h(jc,jb)&
                                         & - dtime*(div_z_depth_int_c(jc,jb)        &
                                         & + p_sfc_flx%forc_tracer(jc,jb,2) ))/gdt2)&
                                         &  *v_base%wet_c(jc,1,jb) !last idx=2 for freshwater
        ENDIF
       !write(*,*)'RHS:',jc,jb,p_os%p_aux%p_rhs_sfc_eq(jc,jb), div_z_c(jc,1,jb)
      ENDDO
    END DO
  ELSEIF(.NOT.l_forc_freshw)THEN

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN 
!             p_os%p_aux%p_rhs_sfc_eq(jc,jb) = ((p_os%p_prog(nold(1))%h(jc,jb)&
!                                            & - dtime*div_z_c(jc,jb))/gdt2)!&
!                                            !& *v_base%wet_c(jc,1,jb)
          p_os%p_aux%p_rhs_sfc_eq(jc,jb) = ((p_os%p_prog(nold(1))%h(jc,jb)&
                                         & - dtime*div_z_depth_int_c(jc,jb))/gdt2)!&
                                         !& *v_base%wet_c(jc,1,jb)
        ENDIF
      ENDDO
    END DO
  ENDIF

ELSE!NOT EDGE-BASED
!  CALL map_edges2cell_3D( p_patch,                 &
!                        & p_os%p_prog(nold(1))%vn, &
!                        & p_op_coeff, z_u_pred_cc)
! 
! ! CALL map_edges2cell_3D( p_patch,                 &
! !                       & p_os%p_prog(nold(1))%vn, &
! !                       & p_op_coeff, z_u_cc2)
! ! IF(iswm_oce == 1)THEN
! !   CALL map_edges2cell_3D( p_patch,           &
! !                         & p_os%p_diag%vn_pred ,&
! !                         & p_op_coeff, z_u_cc1)
! ! 
! ! ELSEIF(iswm_oce /= 1)THEN
! ! 
! !   IF(expl_vertical_velocity_diff==0)THEN
! !     CALL map_edges2cell_3D( p_patch,           &
! !                         & p_os%p_diag%vn_pred ,&
! !                         & p_op_coeff, z_u_cc1)
! ! 
! !   ELSEIF(expl_vertical_velocity_diff==1)THEN
! !     CALL map_edges2cell_3D( p_patch,                    &
! !                         & p_os%p_diag%vn_impl_vert_diff,&
! !                         & p_op_coeff, z_u_cc1)
! !   ENDIF
! ! 
! ! ENDIF
! ! 
! !  DO jb = all_cells%start_block, all_cells%end_block
! !     CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
! !     DO jc = i_startidx_c, i_endidx_c
! !       i_dolic_c = v_base%dolic_c(jc,jb)
! !       DO jk=1,i_dolic_c
! !         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! !           z_u_pred_cc(jc,jk,jb)%x = ab_gam*z_u_cc1(jc,jk,jb)%x&
! !           & + (1.0_wp -ab_gam)*z_u_cc2(jc,jk,jb)%x
! !         ENDIF
! !       END DO
! !     ENDDO
! !   END DO


! ! Step 1) Do within each layer a edge to cell mapping of vn_pred
! ! For below-surface cells, no height has to be provided, and the reconstructions
! ! are normalized by cell area (see primal flip-flop).
! ! For surface cells, the working hypothesis is to do it as in the SW-code, i.e.
! ! we take also the height at the edges into account:

  CALL map_edges2cell_3D( p_patch,              &
                         & z_vn_ab ,            &
                         & p_op_coeff, z_u_pred_cc)

  IF( iswm_oce /= 1 ) THEN !the 3D case

    !calculate depth-integrated velocity 
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF (v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
        !Surface elevation
        delta_z = v_base%del_zlev_m(1) + p_os%p_prog(nold(1))%h(jc,jb)
        z_u_pred_depth_int_cc(jc,jb)%x = z_u_pred_depth_int_cc(jc,jb)%x &
                                  &+ z_u_pred_cc(jc,1,jb)%x        &
                                  &* delta_z
        ENDIF
        i_dolic_c = v_base%dolic_c(jc,jb)
        DO jk=2,i_dolic_c
          IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
            delta_z = v_base%del_zlev_m(jk)
            z_u_pred_depth_int_cc(jc,jb)%x = z_u_pred_depth_int_cc(jc,jb)%x &
                                       &+z_u_pred_cc(jc,jk,jb)%x      &
                                       &* delta_z
          ENDIF
        END DO
      ENDDO
    END DO

  ELSEIF( iswm_oce == 1 ) THEN !the shallow-water case

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        z_u_pred_depth_int_cc(jc,jb)%x = z_u_pred_cc(jc,1,jb)%x&
                                       &*p_os%p_diag%thick_c(jc,jb)
      ENDDO
     END DO
  ENDIF!( iswm_oce /= 1 ) THEN


  CALL map_cell2edges_3D(p_patch,z_u_pred_depth_int_cc, z_e,p_op_coeff, level=1) 
  CALL div_oce_3d( z_e, p_patch,p_op_coeff%div_coeff, div_z_c_2D,&
               & level=1, subset_range=cells_in_domain )



  IF(l_forc_freshw)THEN

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
          p_os%p_aux%p_rhs_sfc_eq(jc,jb) = ((p_os%p_prog(nold(1))%h(jc,jb)  &
                                         & - dtime*(div_z_c_2D(jc,jb)        &
                                         & + p_sfc_flx%forc_tracer(jc,jb,2)))/gdt2)&
                                         &  *v_base%wet_c(jc,1,jb) !last idx=2 for freshwater
        ENDIF
      ENDDO
    END DO

  ELSEIF(.NOT.l_forc_freshw)THEN

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN 
             p_os%p_aux%p_rhs_sfc_eq(jc,jb) = ((p_os%p_prog(nold(1))%h(jc,jb)&
                                            & - dtime*div_z_c_2D(jc,jb))/gdt2)!&
                                            !& *v_base%wet_c(jc,1,jb)
        ENDIF
      ENDDO
    END DO
  ENDIF

ENDIF!EDGE-BASED

 CALL sync_patch_array(SYNC_C, p_patch, p_os%p_aux%p_rhs_sfc_eq )


  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=3  ! output print level (1-5, fix)
  CALL dbg_print('RHS thick_e'               ,p_os%p_diag%thick_e      ,str_module,idt_src)
  CALL dbg_print('RHS z_vn_ab'               ,z_vn_ab                  ,str_module,idt_src)
  IF (l_edge_based) &
    & CALL dbg_print('RHS div_z_depth_int_c' ,div_z_depth_int_c        ,str_module,idt_src)
  IF (.NOT. l_edge_based) &
    & CALL dbg_print('RHS div_z_c_2d'        ,div_z_c_2d               ,str_module,idt_src)
 !CALL dbg_print('RHS div_z_c'               ,div_z_c                  ,str_module,idt_src)
  idt_src=2  ! output print level (1-5, fix)
  CALL dbg_print('RHS final'                 ,p_os%p_aux%p_rhs_sfc_eq  ,str_module,idt_src)
  !---------------------------------------------------------------------


END SUBROUTINE fill_rhs4surface_eq_ab

!-------------------------------------------------------------------------------------
!
!  
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

FUNCTION lhs_surface_height_ab_mim( p_x, h_old, p_patch, coeff, h_e,&
                                  & thickness_c,p_op_coeff) RESULT(p_lhs)
!
REAL(wp),    INTENT(INOUT)       :: p_x(:,:) ! inout for sync, dimension: (nproma,p_patch%nblks_c)
REAL(wp),    INTENT(IN)          :: h_old(:,:)
TYPE(t_patch), TARGET, INTENT(in):: p_patch          ! patch on which computation is performed
REAL(wp),    INTENT(in)          :: coeff
TYPE(t_operator_coeff),INTENT(in):: p_op_coeff
REAL(wp),    INTENT(in)          :: h_e(:,:)         !SW-case: thickness at edges
                                              !3D-case: surface height above zero at edges
REAL(wp),    INTENT(in)   :: thickness_c(:,:) !thickness of fluid column    
!
! Left-hand side calculated from iterated height
!
REAL(wp) :: p_lhs(SIZE(p_x,1), SIZE(p_x,2))  ! (nproma,p_patch%nblks_c)
!
! local variables
REAL(wp) :: gdt2
REAL(wp) :: z_grad_h(nproma,p_patch%nblks_e)
REAL(wp) :: z_e(nproma,p_patch%nblks_e)
REAL(wp) :: div_z_c(SIZE(p_x,1), SIZE(p_x,2))  ! (nproma,1,p_patch%nblks_c)
TYPE(t_cartesian_coordinates) :: z_grad_h_cc(nproma,p_patch%nblks_c)
INTEGER :: i_startidx, i_endidx
INTEGER :: jc, jb
TYPE(t_subset_range), POINTER :: cells_in_domain, all_cells

! CHARACTER(len=max_char_length), PARAMETER ::     &
!   &      routine = ('mo_oce_ab_timestepping_mimetic: lhs_surface_height_ab_mim')
!-----------------------------------------------------------------------  
!CALL message (TRIM(routine), 'start - iteration by GMRES')        

cells_in_domain => p_patch%cells%in_domain
all_cells => p_patch%cells%all

p_lhs(:,:)    = 0.0_wp
div_z_c(:,:)  = 0.0_wp
z_e(:,:)      = 0.0_wp
z_grad_h(:,:) = 0.0_wp

gdt2 = grav*(dtime)**2
!z_u_c     = 0.0_wp
!z_v_c     = 0.0_wp

CALL sync_patch_array(SYNC_C, p_patch, p_x )

!Step 1) Calculate gradient of iterated height.
CALL grad_fd_norm_oce_2d_3D( p_x, &
         &                   p_patch,                &
         &                   p_op_coeff%grad_coeff,  &
         &                   z_grad_h(:,:))

CALL sync_patch_array(SYNC_E, p_patch, z_grad_h(:,:) )

IF(l_edge_based)THEN
  z_e(:,:)=z_grad_h(:,:)*h_e(:,:)
ELSE
!Step 2) map the gradient to the cell center, multiply it
!by fluid thickness and map the result back to edges 
IF( iswm_oce /= 1 ) THEN !the 3D case

  CALL map_edges2cell_3D( p_patch,       &
                   & z_grad_h,           &
                   & p_op_coeff,         &
                   & z_grad_h_cc,        &
                   & level=top        )


ELSEIF( iswm_oce == 1 ) THEN !the shallow-water case

 CALL map_edges2cell_3D( p_patch,        &
                   & z_grad_h,           &
                   & p_op_coeff,         &
                   & z_grad_h_cc,        &
                   & level=top  )

ENDIF

DO jb = all_cells%start_block, all_cells%end_block
  CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
  DO jc = i_startidx, i_endidx
    IF ( v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
    z_grad_h_cc(jc,jb)%x = z_grad_h_cc(jc,jb)%x *thickness_c(jc,jb)
    ENDIF 
  END DO
END DO

 CALL map_cell2edges_3D( p_patch,    &
                    & z_grad_h_cc,   &
                    & z_e,p_op_coeff,&
                    & level=top)
ENDIF
!write(*,*)'LHS:max/min :', maxval(z_e(:,1,:)), minval(z_e(:,1,:))

!Step 3) Calculate divergence
CALL div_oce_3D( z_e, p_patch, p_op_coeff%div_coeff, div_z_c, &
  & level=top, subset_range=cells_in_domain  )

!Step 4) Finalize LHS calculations
DO jb = cells_in_domain%start_block, cells_in_domain%end_block
  CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
  DO jc = i_startidx, i_endidx
    IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
      p_lhs(jc,jb) =(p_x(jc,jb)- gdt2*ab_gam*ab_beta*div_z_c(jc,jb))/gdt2
   ELSE
      p_lhs(jc,jb) =0.0_wp
   ENDIF
  END DO
END DO

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
SUBROUTINE calc_normal_velocity_ab_mimetic(p_patch, p_os, p_op_coeff, p_ext_data, p_phys_param)
!
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch
  TYPE(t_hydro_ocean_state), TARGET :: p_os
  TYPE(t_operator_coeff),INTENT(IN) :: p_op_coeff
  TYPE(t_external_data), TARGET     :: p_ext_data
  TYPE (t_ho_params)                :: p_phys_param
  !
  !  local variables
  INTEGER :: i_startidx_e, i_endidx_e
  INTEGER :: i_startidx_c, i_endidx_c
  INTEGER :: je, jk, jb, jc
  REAL(wp) :: z_grad_h(nproma,p_patch%nblks_e)
  REAL(wp) :: gdt, delta_z
  TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
  !TYPE(t_cartesian_coordinates) :: z_u_cc(nproma,n_zlev,p_patch%nblks_c)
  !TYPE(t_cartesian_coordinates) :: z_u_cc2(nproma,n_zlev,p_patch%nblks_c)
  REAL(wp) :: z_vn_ab(nproma,n_zlev,p_patch%nblks_e)
  CHARACTER(len=*), PARAMETER ::     &
    &      method_name='mo_oce_ab_timestepping_mimetic: calc_normal_velocity_ab_mimetic'

  !----------------------------------------------------------------------
  !CALL message (TRIM(routine), 'start')
  !-----------------------------------------------------------------------
  edges_in_domain => p_patch%edges%in_domain
  cells_in_domain => p_patch%cells%in_domain

  z_grad_h(1:nproma,1:p_patch%nblks_e) = 0.0_wp

  gdt=grav*dtime

  ! Step 1) Compute normal derivative of new surface height
  IF(.NOT.l_RIGID_LID.OR.iswm_oce == 1) THEN
    CALL grad_fd_norm_oce_2d_3D( p_os%p_prog(nnew(1))%h,&
      &                  p_patch,                  &
      &                  p_op_coeff%grad_coeff,    &
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
          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
            p_os%p_prog(nnew(1))%vn(je,jk,jb) = (p_os%p_diag%vn_pred(je,jk,jb) &
                                              &- gdt*ab_beta*z_grad_h(je,jb))&
                                              &*v_base%wet_e(je,jk,jb)
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

          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
            p_os%p_prog(nnew(1))%vn(je,jk,jb) = (p_os%p_diag%vn_pred(je,jk,jb) &
                                              &- gdt*ab_beta*z_grad_h(je,jb))&
                                              &*v_base%wet_e(je,jk,jb)
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
            IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
              p_os%p_diag%vn_time_weighted(je,jk,jb) &
              &= ab_gam*p_os%p_prog(nnew(1))%vn(je,jk,jb) &
              &+ (1.0_wp -ab_gam)*p_os%p_prog(nold(1))%vn(je,jk,jb)
            ENDIF  
          END DO
        END DO
      END DO

  CALL sync_patch_array(SYNC_E, p_patch, p_os%p_prog(nnew(1))%vn)
  CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%vn_time_weighted)

  ! slo: z_grad_h not out of sync, but sync error with global_max in dbg_print?
  !CALL sync_patch_array(SYNC_E, p_patch, z_grad_h)

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

  !CALL height_related_quantities(p_patch, p_os, p_ext_data)
  ! Update of scalar product quantities
  IF(l_STAGGERED_TIMESTEP)THEN
    CALL height_related_quantities(p_patch, p_os, p_ext_data)

    ! #slo# vn-new is already multiplied with wet_c - not necessary?
    Call set_lateral_boundary_values( p_patch, &
                                    & p_os%p_prog(nnew(1))%vn)

    CALL calc_scalar_product_veloc_3D( p_patch,                &
                                     & p_os%p_prog(nnew(1))%vn,&
                                     & p_os%p_prog(nnew(1))%vn,&
                                     & p_os%p_diag%h_e,        &
                                     & p_os%p_diag,            &
                                     & p_op_coeff)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('NorVel: Staggered, kin'    ,p_os%p_diag%kin        ,str_module,idt_src)
    CALL dbg_print('NorVel: Staggered, ptp_vn' ,p_os%p_diag%ptp_vn     ,str_module,idt_src)
    !---------------------------------------------------------------------

  ENDIF

  !CALL message (TRIM(routine), 'end')

END SUBROUTINE calc_normal_velocity_ab_mimetic
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!>
!! 
!! @par Revision History
!! Developed  by  Peter Korn,   MPI-M (2012).
!! 
  !!  mpi parallelized LL
SUBROUTINE update_column_thickness( p_patch, p_os, p_diag, p_op_coeff, timestep)
!
TYPE(t_patch), TARGET, INTENT(IN) :: p_patch       ! patch on which computation is performed
TYPE(t_hydro_ocean_state)         :: p_os
TYPE(t_hydro_ocean_diag)          :: p_diag
TYPE(t_operator_coeff),INTENT(IN) :: p_op_coeff
INTEGER                           :: timestep

! Local variables
INTEGER  :: jc, jk, jb
INTEGER  :: i_startidx_c, i_endidx_c
INTEGER  :: i_dolic_c
!REAL(wp) :: delta_z
REAL(wp) :: z_thick_c_old(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_div_c(nproma,p_patch%nblks_c)
REAL(wp) :: z_e(nproma,p_patch%nblks_e)

TYPE(t_cartesian_coordinates) :: z_vn_c(nproma,n_zlev,p_patch%nblks_c)
TYPE(t_cartesian_coordinates) ::z_u_depth_int_cc(nproma,p_patch%nblks_c)
TYPE(t_subset_range), POINTER :: cells_in_domain
CHARACTER(len=*), PARAMETER :: &
   & method_name = 'mo_oce_ab_timestepping_mimetic:update_column_thickness'
!-----------------------------------------------------------------------  
  cells_in_domain => p_patch%cells%in_domain


  IF(timestep==1)THEN
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      !DO jk = 1, n_zlev
        !delta_z = v_base%del_zlev_m(jk)
        DO jc = i_startidx_c, i_endidx_c
        !IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
        !  p_os%p_diag%depth_c(jc,jk,jb) = delta_z
          p_os%p_diag%cons_thick_c(jc,1,jb)=p_os%p_diag%thick_c(jc,jb)
        !ENDIF
        END DO
      !END DO
    END DO
  ENDIF
  z_thick_c_old = p_os%p_diag%cons_thick_c


  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
    DO jc = i_startidx_c, i_endidx_c
      i_dolic_c = v_base%dolic_c(jc,jb)
      z_u_depth_int_cc(jc,jb)%x(1:3) = 0.0_wp
      DO jk=1,i_dolic_c
         z_u_depth_int_cc(jc,jb)%x = z_u_depth_int_cc(jc,jb)%x&
                                  &+ p_os%p_diag%p_vn(jc,jk,jb)%x&
                                  &* v_base%del_zlev_m(jk)
      END DO

     !Add surface elevation
     !For a linear surface this can be eleminated
      z_u_depth_int_cc(jc,jb)%x = z_u_depth_int_cc(jc,jb)%x&
                               &+ p_os%p_diag%p_vn(jc,1,jb)%x&
                               &* p_os%p_prog(nold(1))%h(jc,jb)
    ENDDO
  END DO

  CALL map_cell2edges_3D( p_patch,          &
                      & z_u_depth_int_cc,   &
                      & z_e, p_op_coeff,    &
                      & level=1)

  CALL div_oce_3D( z_e, p_patch,p_op_coeff%div_coeff, z_div_c,&
               & level=1, subset_range=cells_in_domain )

  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
    DO jc = i_startidx_c, i_endidx_c

      IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
        p_os%p_diag%cons_thick_c(jc,1,jb) = (z_thick_c_old(jc,1,jb)& !p_os%p_prog(nold(1))%h(jc,jb)&
                                       & - dtime*z_div_c(jc,jb))&
                                       &  *v_base%wet_c(jc,1,jb)            
       ENDIF
    ENDDO
  END DO
!write(*,*)'update column thickness', maxval(p_os%p_diag%cons_thick_c),&
!&minval(p_os%p_diag%cons_thick_c)



END SUBROUTINE update_column_thickness
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!>
!! Computation of new vertical velocity using continuity equation

!! Calculate diagnostic vertical velocity from horizontal velocity using the
!! incommpressibility condition in the continuity equation.
!! For the case of the semi-implicit-AB scheme the land-sea-mask may be applied
!! at least after collecting the whole explicit term.
!! 
!! @par Revision History
!! Developed  by  Peter Korn,   MPI-M (2006).
!!  Modified by Stephan Lorenz, MPI-M (2010-06)
!!  Modified by Stephan Lorenz, MPI-M (2010-08)
!!   - velocities and sea level are passed through the interface
!!   - no calculation of new sea level here
!! 
  !!  mpi parallelized LL
SUBROUTINE calc_vert_velocity_mimetic( p_patch, p_os, p_diag, p_op_coeff,&
                                     & ph_c, ph_e, bot_bc_w, pw_c )
!
TYPE(t_patch), TARGET, INTENT(IN) :: p_patch       ! patch on which computation is performed
TYPE(t_hydro_ocean_state)         :: p_os
TYPE(t_hydro_ocean_diag)          :: p_diag
TYPE(t_operator_coeff),INTENT(IN) :: p_op_coeff
REAL(wp),         INTENT(INOUT)   :: ph_c(:,:)
REAL(wp),         INTENT(INOUT)   :: ph_e(:,:)
!REAL(wp),            INTENT(IN)   :: top_bc_w(:,:) !bottom boundary condition for vertical velocity
REAL(wp),            INTENT(IN)   :: bot_bc_w(:,:) !bottom boundary condition for vertical velocity
REAL(wp),         INTENT(INOUT)   :: pw_c (:,:,:)  ! vertical velocity on cells
!
!
! Local variables
INTEGER :: jc, jk, jb, je
INTEGER :: z_dolic
INTEGER :: i_startidx, i_endidx
REAL(wp) :: delta_z
REAL(wp) :: z_div_c(nproma,n_zlev+1,p_patch%nblks_c)
REAL(wp) :: z_div_int_c(nproma,p_patch%nblks_c)
REAL(wp) :: z_vn(nproma,n_zlev,p_patch%nblks_e)
TYPE(t_cartesian_coordinates):: z_vn_c(nproma,n_zlev,p_patch%nblks_c)
TYPE(t_subset_range), POINTER :: cells_in_domain
TYPE(t_subset_range), POINTER :: edges_in_domain
!INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
CHARACTER(len=*), PARAMETER :: &
   & method_name = 'mo_oce_ab_timestepping_mimetic:calc_vert_velocity_mimetic'
!-----------------------------------------------------------------------  
cells_in_domain => p_patch%cells%in_domain
edges_in_domain => p_patch%edges%in_domain

! due to nag -nan compiler-option:
z_div_int_c(:,:) = 0.0_wp
z_div_c  (:,:,:) = 0.0_wp
pw_c     (:,:,:) = 0.0_wp
z_vn     (:,:,:) = 0.0_wp
!------------------------------------------------------------------
! Step 1) Calculate divergence of horizontal velocity at all levels
!------------------------------------------------------------------
IF(l_edge_based)THEN
  DO jb = edges_in_domain%start_block, edges_in_domain%end_block
    CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
    DO jk = 1, n_zlev
      delta_z = v_base%del_zlev_m(jk)
      DO je = i_startidx, i_endidx
        IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
          IF(jk==1)THEN
             delta_z=v_base%del_zlev_m(jk) + p_os%p_diag%h_e(je,jb)!&
          ENDIF
          z_vn(je,jk,jb) = p_diag%vn_time_weighted(je,jk,jb)*delta_z
        ENDIF
      END DO
    END DO
  END DO
 CALL div_oce_3D(z_vn, p_patch,p_op_coeff%div_coeff, z_div_c, subset_range=cells_in_domain)
 ! !Note we are summing from bottom up to one layer below top.
  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
    DO jc = i_startidx, i_endidx

      z_dolic = v_base%dolic_c(jc,jb)
      IF ( z_dolic>=MIN_DOLIC)THEN !

        !use bottom boundary condition for vertical velocity at bottom
        !of prism
        !delta_z = v_base%del_zlev_m(z_dolic)
        pw_c(jc,z_dolic,jb) = bot_bc_w(jc,jb) - z_div_c(jc,z_dolic,jb)!*delta_z

         DO jk = z_dolic-1, 1, -1
          IF(v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          !delta_z = v_base%del_zlev_m(jk)
          !IF (jk == 1)THEN 
          ! !at surface level (toplev) add surface elevation h_c
          !delta_z = v_base%del_zlev_m(jk) + ph_c(jc,jb)
          !ENDIF

          ! vertical velocity is integrated from bottom to top
          ! vertical velocity is negative for positive divergence
          ! of horizontal velocity
          pw_c(jc,jk,jb) = pw_c(jc,jk+1,jb) - z_div_c(jc,jk,jb)!*delta_z
          ENDIF
        END DO
      END IF
    END DO
  END DO

ELSE
  CALL map_edges2cell_3D(  p_patch,                      &
                         & p_os%p_diag%vn_time_weighted, &
                         & p_op_coeff, z_vn_c)
  !CALL map_edges2cell( p_patch, p_diag%vn_time_weighted, z_vn_c)!, ph_e)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
        DO jc = i_startidx, i_endidx
          z_dolic = v_base%dolic_c(jc,jb)
          IF ( z_dolic>=MIN_DOLIC)THEN !v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          DO jk = 1,z_dolic
            delta_z = v_base%del_zlev_m(jk)
            IF (jk == 1) delta_z = v_base%del_zlev_m(jk) + p_os%p_prog(nold(1))%h(jc,jb)
            z_vn_c(jc,jk,jb)%x=delta_z*z_vn_c(jc,jk,jb)%x
          END DO
        END IF
      END DO
    END DO
  CALL map_cell2edges_3D( p_patch, z_vn_c, z_vn,p_op_coeff)
  CALL div_oce_3D(z_vn, p_patch,p_op_coeff%div_coeff, z_div_c, subset_range=cells_in_domain)

  !CALL div_oce_3D(p_diag%vn_time_weighted, p_patch,p_op_coeff%div_coeff, z_div_c,&
  !                & subset_range=cells_in_domain)
  ! !Note we are summing from bottom up to one layer below top.
  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
    DO jc = i_startidx, i_endidx

      z_dolic = v_base%dolic_c(jc,jb)
      IF ( z_dolic>=MIN_DOLIC)THEN
        !use bottom boundary condition for vertical velocity at bottom
        !of prism
        !delta_z = v_base%del_zlev_m(z_dolic)
        pw_c(jc,z_dolic,jb) = bot_bc_w(jc,jb) - z_div_c(jc,z_dolic,jb)*delta_z

         DO jk = z_dolic-1, 2, -1
          IF(v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          delta_z = v_base%del_zlev_m(jk)
          IF (jk == 1)THEN 
          ! !at surface level (toplev) add surface elevation h_c
          delta_z = v_base%del_zlev_m(jk) + ph_c(jc,jb)
          ENDIF

          ! vertical velocity is integrated from bottom to top
          ! vertical velocity is negative for positive divergence
          ! of horizontal velocity
          pw_c(jc,jk,jb) = pw_c(jc,jk+1,jb) - z_div_c(jc,jk,jb)*delta_z
          IF(jk>=1)&
          &z_div_int_c(jc,jb)=z_div_int_c(jc,jb)+z_div_c(jc,jk,jb)

          ENDIF
        END DO
      END IF
    END DO
  END DO
ENDIF ! l_edge_based


IF(l_RIGID_LID)THEN
  pw_c(:,2,:) = 0.0_wp
ENDIF

CALL sync_patch_array(SYNC_C, p_patch, pw_c)

!---------DEBUG DIAGNOSTICS-------------------------------------------
idt_src=4  ! output print level (1-5, fix)
CALL dbg_print('CalcVertVelMim: z_vn'        ,z_vn                     ,str_module,idt_src)
CALL dbg_print('CalcVertVelMim: div(v)_c'    ,z_div_c                  ,str_module,idt_src)
CALL dbg_print('CalcVertVelMim: div_int_c'   ,z_div_int_c              ,str_module,idt_src)
idt_src=3  ! output print level (1-5, fix)
CALL dbg_print('CalcVertVelMim: pw_c'        ,pw_c                     ,str_module,idt_src)
!---------------------------------------------------------------------

END SUBROUTINE calc_vert_velocity_mimetic
!-------------------------------------------------------------------------
!
!  
!>
!! Computation of new vertical velocity using continuity equation

!! Calculate diagnostic vertical velocity from horizontal velocity using the
!! incommpressibility condition in the continuity equation.
!! For the case of the semi-implicit-AB scheme the land-sea-mask may be applied
!! at least after collecting the whole explicit term.
!! 
!! @par Revision History
!! Developed  by  Peter Korn,   MPI-M (2006).
!!  Modified by Stephan Lorenz, MPI-M (2010-06)
!!  Modified by Stephan Lorenz, MPI-M (2010-08)
!!   - velocities and sea level are passed through the interface
!!   - no calculation of new sea level here
!! mpi not used
SUBROUTINE calc_vert_velocity_mim_topdown( p_patch, p_os, p_diag,p_op_coeff, &
                                          &ph_e, top_bc_w, bot_bc_w, pw_c )
!
TYPE(t_patch), TARGET, INTENT(IN) :: p_patch       ! patch on which computation is performed
TYPE(t_hydro_ocean_state)         :: p_os
TYPE(t_hydro_ocean_diag)          :: p_diag
TYPE(t_operator_coeff),INTENT(IN) :: p_op_coeff
REAL(wp),         INTENT(INOUT)   :: ph_e(:,:)  ! 
REAL(wp),            INTENT(IN)   :: top_bc_w(:,:) !bottom boundary condition for vertical velocity
REAL(wp),            INTENT(IN)   :: bot_bc_w(:,:) !bottom boundary condition for vertical velocity
REAL(wp),         INTENT(INOUT)   :: pw_c (:,:,:)  ! vertical velocity on cells
!
!
! Local variables
INTEGER :: jc, jk, jb, je
INTEGER :: z_dolic!jic, jib
INTEGER :: i_startidx, i_endidx, i_startblk, i_endblk
INTEGER :: i_startidx_c, i_endidx_c, i_startblk_c, i_endblk_c
!INTEGER :: rl_start, rl_end
REAL(wp) :: delta_z
REAL(wp) :: z_div_c(nproma,n_zlev+1,p_patch%nblks_c)
REAL(wp) :: z_pw_bottom_up(nproma,n_zlev+1,p_patch%nblks_c)
REAL(wp) :: z_vn(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_div_int_c(nproma,p_patch%nblks_c)
REAL(wp) :: div_z_c_2D(nproma,p_patch%nblks_c)
REAL(wp) :: z_e(nproma,p_patch%nblks_e)
INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc 
!REAL(wp) :: div_z_c_3D_accum(nproma,n_zlev,p_patch%nblks_c)
TYPE(t_cartesian_coordinates):: flux_sum
TYPE(t_cartesian_coordinates):: z_vn_c(nproma,n_zlev,p_patch%nblks_c)
TYPE(t_cartesian_coordinates):: z_u_depth_int_cc(nproma,p_patch%nblks_c)
TYPE(t_subset_range), POINTER:: all_cells, cells_in_domain, edges_in_domain
CHARACTER(len=*), PARAMETER :: &
  & method_name = ('mo_oce_ab_timestepping_mimetic:calc_vert_velocity_mimetic')
!-----------------------------------------------------------------------  

  all_cells       => p_patch%cells%all
  cells_in_domain => p_patch%cells%in_domain
  edges_in_domain => p_patch%edges%in_domain

! due to nag -nan compiler-option:
  pw_c           = 0.0_wp
  z_div_int_c    = 0.0_wp
  z_pw_bottom_up = 0.0_wp
  z_div_c(:,:,:) = 0.0_wp
  z_vn(1:nproma,1:n_zlev,1:p_patch%nblks_e)=0.0_wp
  !div_z_c_3D_accum=0.0_wp

  z_vn_c(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(1) = 0.0_wp
  z_vn_c(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(2) = 0.0_wp
  z_vn_c(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(3) = 0.0_wp

!------------------------------------------------------------------
! Step 1) Calculate divergence of horizontal velocity at all levels
!------------------------------------------------------------------
IF(l_EDGE_BASED)THEN

  DO jb = edges_in_domain%start_block, edges_in_domain%end_block
    CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
    DO jk = 1, n_zlev
      delta_z = v_base%del_zlev_m(jk)
      DO je = i_startidx, i_endidx
        IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
          IF(jk==1)THEN
             delta_z=v_base%del_zlev_m(jk) + p_os%p_diag%h_e(je,jb)!&
          ENDIF
          z_vn(je,jk,jb) = p_diag%vn_time_weighted(je,jk,jb)*delta_z
        ENDIF
      END DO
    END DO
  END DO
 CALL div_oce_3D(z_vn, p_patch,p_op_coeff%div_coeff, z_div_c, subset_range=cells_in_domain)
 ! !Note we are summing from bottom up to one layer below top.
  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
    DO jc = i_startidx, i_endidx
      z_dolic = v_base%dolic_c(jc,jb)
      IF ( z_dolic>=MIN_DOLIC)THEN

        !Top layer
        pw_c(jc,1,jb)=-sum(z_div_c(jc,:,jb))

        DO jk = 2,z_dolic
          IF(v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          pw_c(jc,jk,jb)     = pw_c(jc,jk-1,jb) + z_div_c(jc,jk-1,jb)
          ENDIF
        END DO
        IF(v_base%lsm_oce_c(jc,z_dolic,jb) <= sea_boundary ) THEN
          pw_c(jc,z_dolic+1,jb) = pw_c(jc,z_dolic,jb) + z_div_c(jc,z_dolic,jb)
        ENDIF

      END IF
    END DO
  END DO
  !write(*,*)'max/min difference',&
  ! &maxval((p_os%p_prog(nnew(1))%h-p_os%p_prog(nold(1))%h)/dtime - pw_c(:,1,:))


ELSEIF(.NOT.l_EDGE_BASED)THEN
!CALL sync_patch_array(SYNC_C,p_patch,pw_c)
!----------------------------------------------------------------------------------
!  !z_u_depth_int_cc(1:nproma,1:p_patch%nblks_c)%x(1) = 0.0_wp
!  !z_u_depth_int_cc(1:nproma,1:p_patch%nblks_c)%x(2) = 0.0_wp
!  !z_u_depth_int_cc(1:nproma,1:p_patch%nblks_c)%x(3) = 0.0_wp!
!   CALL map_edges2cell_3D( p_patch, p_diag%vn_time_weighted,p_op_coeff, z_vn_c)
!   !calculate depth-integrated velocity 
!     DO jk=n_zlev,2,-1
!       DO jb = all_cells%start_block, all_cells%end_block
!         CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
!         DO jc = i_startidx_c, i_endidx_c
!           IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
!             delta_z = v_base%del_zlev_m(jk)
!             z_u_depth_int_cc(jc,jb)%x = z_u_depth_int_cc(jc,jb)%x &
!                                        &+ z_vn_c(jc,jk,jb)%x      &
!                                        &* delta_z
!           ENDIF
!         END DO
!       ENDDO
!       CALL map_cell2edges_3D(p_patch,z_u_depth_int_cc, z_e,p_op_coeff, level=jk) 
!       CALL div_oce_3d( z_e, p_patch,p_op_coeff%div_coeff, div_z_c_2D,&
!                      & level=jk, subset_range=cells_in_domain )
!        div_z_c_3D_accum(:,jk,:)=div_z_c_3D_accum(:,jk,:)-div_z_c_2D(:,:)
!     END DO
!     DO jb = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
!       DO jc = i_startidx_c, i_endidx_c
!        IF (v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
!        !Surface elevation
!        delta_z = v_base%del_zlev_m(1) + p_os%p_prog(nold(1))%h(jc,jb)
!        z_u_depth_int_cc(jc,jb)%x = z_u_depth_int_cc(jc,jb)%x &
!                                  &+ z_vn_c(jc,1,jb)%x        &
!                                  &* delta_z
!        ENDIF
!       END DO
!     END DO
!        CALL map_cell2edges_3D(p_patch,z_u_depth_int_cc, z_e,p_op_coeff, level=1) 
!        CALL div_oce_3d( z_e, p_patch,p_op_coeff%div_coeff, div_z_c_2D,&
!                      & level=1, subset_range=cells_in_domain )
!        div_z_c_3D_accum(:,1,:)=div_z_c_3D_accum(:,1,:)-div_z_c_2D(:,:)
!    CALL map_cell2edges_3D(p_patch,z_u_depth_int_cc,z_e,p_op_coeff, level=1) 
!    CALL div_oce_3d( z_e, p_patch,p_op_coeff%div_coeff, div_z_c_2D,&
!                 & level=1, subset_range=cells_in_domain )

   CALL map_edges2cell_3D( p_patch, p_diag%vn_time_weighted, p_op_coeff, z_vn_c)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
      DO jc = i_startidx, i_endidx
        z_dolic = v_base%dolic_c(jc,jb)
        DO jk = 1,z_dolic
          IF (v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
            delta_z = v_base%del_zlev_m(jk)
            IF(jk==1)THEN
              delta_z=v_base%del_zlev_m(jk) + p_os%p_prog(nold(1))%h(jc,jb)
            ENDIF
            z_vn_c(jc,jk,jb)%x=z_vn_c(jc,jk,jb)%x *delta_z
          END IF
        END DO
      END DO
    END DO
  CALL sync_patch_array(SYNC_C,p_patch,z_vn_c(:,:,:)%x(1))
  CALL sync_patch_array(SYNC_C,p_patch,z_vn_c(:,:,:)%x(2))
  CALL sync_patch_array(SYNC_C,p_patch,z_vn_c(:,:,:)%x(3))

  SELECT CASE(FLUX_CALCULATION_HORZ)

  CASE(MIMETIC) 

    CALL map_cell2edges_3D( p_patch, z_vn_c, z_vn,p_op_coeff)

  CASE(CENTRAL)

    iilc => p_patch%edges%cell_idx
    iibc => p_patch%edges%cell_blk

    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
      DO jk = 1, n_zlev
        delta_z = v_base%del_zlev_m(jk)
        DO je = i_startidx, i_endidx
          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
            flux_sum%x=0.5_wp*(z_vn_c(iilc(je,jb,1),jk,iibc(je,jb,1))%x&
                             &+z_vn_c(iilc(je,jb,2),jk,iibc(je,jb,2))%x)
            z_vn(je,jk,jb)=DOT_PRODUCT(flux_sum%x,p_patch%edges%primal_cart_normal(je,jb)%x)
          ENDIF
        END DO
      END DO
    END DO

  CASE(UPWIND)

    iilc => p_patch%edges%cell_idx
    iibc => p_patch%edges%cell_blk

    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
      DO jk = 1, n_zlev
        delta_z = v_base%del_zlev_m(jk)
        DO je = i_startidx, i_endidx
          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
            IF ( p_os%p_diag%vn_time_weighted(je,jk,jb) > 0.0_wp ) THEN

              z_vn(je,jk,jb)=DOT_PRODUCT(z_vn_c(iilc(je,jb,1),jk,iibc(je,jb,1))%x,&
                            &p_patch%edges%primal_cart_normal(je,jb)%x)

            ELSEIF(p_os%p_diag%vn_time_weighted(je,jk,jb) <= 0.0_wp ) THEN
              z_vn(je,jk,jb)=DOT_PRODUCT(z_vn_c(iilc(je,jb,2),jk,iibc(je,jb,2))%x,&
                            &p_patch%edges%primal_cart_normal(je,jb)%x)
            ENDIF
          ENDIF
        END DO
      END DO
    END DO

  END SELECT
  CALL sync_patch_array(SYNC_E,p_patch,z_vn)
  CALL div_oce_3D(z_vn, p_patch,p_op_coeff%div_coeff, z_div_c, subset_range=cells_in_domain)
  CALL sync_patch_array(SYNC_C,p_patch,z_div_c)


  ! !Note we are summing from top to bottom.
  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
    DO jc = i_startidx, i_endidx

      z_dolic = v_base%dolic_c(jc,jb)
      IF(z_dolic>0)THEN
!            !WRITE(1234,*)'----------------------',z_dolic 
!           DO jk = z_dolic,1,-1
!             IF(v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
!             z_pw_bottom_up(jc,jk,jb) = z_pw_bottom_up(jc,jk+1,jb) - z_div_c(jc,jk,jb)
!           ! z_pw_bottom_up(jc,jk,jb) = div_z_c_3D_accum(jc,jk,jb)!z_div_c(jc,jk,jb)
!            !z_div_int_c(jc,jb) = z_div_int_c(jc,jb)-z_div_c(jc,jk,jb)
!            !pw_c(jc,jk,jb) = pw_c(jc,jk+1,jb) - z_div_c(jc,jk,jb)
!            ENDIF
!          END DO

         !Top
         pw_c(jc,1,jb)=-sum(z_div_c(jc,:,jb))

         DO jk = 2,z_dolic
          IF(v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          pw_c(jc,jk,jb)     = pw_c(jc,jk-1,jb) + z_div_c(jc,jk-1,jb)
          ENDIF
        END DO
!          DO jk = 1,z_dolic
!           IF(v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! write(1234,*)'details',jc,jk,jb, z_pw_bottom_up(jc,jk,jb),pw_c(jc,jk,jb),&
! !& pw_c(jc,jk-1,jb),!&z_div_c(jc,jk,jb),!&div_z_c_2D(jc,jb),z_div_int_c(jc,jb),&
! !&div_z_c_3D_accum(jc,jk,jb), &
! !& z_div_int_c(jc,jb)- pw_c(jc,1,jb),&
! !& z_div_int_c(jc,jb)- z_pw_bottom_up(jc,1,jb),&
! &(p_os%p_prog(nnew(1))%h(jc,jb)-p_os%p_prog(nold(1))%h(jc,jb))/dtime - z_pw_bottom_up(jc,1,jb)!,&
! !&(p_os%p_prog(nnew(1))%h(jc,jb)-p_os%p_prog(nold(1))%h(jc,jb))/dtime -div_z_c_3D_accum(jc,1,jb)!div_z_c_2D(jc,jb)!,&!
!           ENDIF
!          END DO
! WRITE(1234,*)'----------------------'
      END IF
    END DO
  END DO

ENDIF

IF(l_RIGID_LID)THEN
  pw_c(:,1,:) = 0.0_wp
ENDIF
CALL sync_patch_array(SYNC_C,p_patch,pw_c)

!---------DEBUG DIAGNOSTICS-------------------------------------------
idt_src=4  ! output print level (1-5, fix)
CALL dbg_print('CalcVertVelMimTD: z_vn'      ,z_vn                     ,str_module,idt_src)
CALL dbg_print('CalcVertVelMimTD: div(v)_c'  ,z_div_c                  ,str_module,idt_src)
idt_src=3  ! output print level (1-5, fix)
CALL dbg_print('CalcVertVelMimTD: pw_c'      ,pw_c                     ,str_module,idt_src)
!---------------------------------------------------------------------

END SUBROUTINE calc_vert_velocity_mim_topdown
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!!  mpi parallelized, the result is NOT synced. Should be done in the calling method if required
FUNCTION inverse_primal_flip_flop(p_patch, rhs_e, h_e) result(inv_flip_flop_e)
   !
   TYPE(t_patch), TARGET :: p_patch 
   REAL(wp)      :: rhs_e(:,:,:)!(nproma,n_zlev,p_patch%nblks_e)
   REAL(wp)      :: h_e(:,:)  !(nproma,p_patch%nblks_e)
   REAL(wp)      :: inv_flip_flop_e(SIZE(rhs_e,1),SIZE(rhs_e,2),SIZE(rhs_e,3))
   !
   !LOCAL VARIABLES
   INTEGER,PARAMETER :: nmax_iter= 200 ! maximum number of iterations
   REAL(wp) :: zimpl_coeff = 1.0_wp    !COEFF has to be set appropriately !!!!
   REAL(wp) :: zimpl_prime_coeff
   INTEGER  :: n_iter                  ! number of iterations
   REAL(wp) :: tolerance               ! (relative or absolute) tolerance
   REAL(wp) :: z_residual(nmax_iter)
   LOGICAL  :: lmax_iter               ! true if reached m iterations
   !LOGICAL  :: lverbose = .TRUE.
   !CHARACTER(len=MAX_CHAR_LENGTH) :: string
   REAL(wp) :: rhstemp(nproma,p_patch%nblks_e)
   REAL(wp), ALLOCATABLE :: inv_flip_flop_e2(:,:)!(nproma,p_patch%nblks_e)
   REAL(wp) :: z_e(nproma,n_zlev,p_patch%nblks_e)
   INTEGER :: jk
   !INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
   !INTEGER :: rl_start_e, rl_end_e, je,jb

   !-----------------------------------------------------------------------

   tolerance                = 1.0e-10_wp  ! solver_tolerance
   inv_flip_flop_e(:,:,:)   = 0.0_wp
   zimpl_prime_coeff = (1.0_wp-zimpl_coeff)

   ALLOCATE (inv_flip_flop_e2(nproma,p_patch%nblks_e))
   inv_flip_flop_e2(:,:) = 0.0_wp
   rhstemp(:,:)          = 0.0_wp


   DO jk=1, n_zlev
     inv_flip_flop_e2(:,:) = rhs_e(:,jk,:)
     rhstemp(:,:) = rhs_e(:,jk,:)&
       & -zimpl_prime_coeff*lhs_primal_flip_flop(inv_flip_flop_e2,p_patch,jk,zimpl_coeff,h_e) 

     If (maxval (ABS (rhstemp (:,:))) <= tolerance) THEN
       inv_flip_flop_e(:,jk,:) = inv_flip_flop_e2(:,:)
       print*, "Inv_flipflop GMRES solved by initial guess!",&
         & jk,MAXVAL(ABS(rhstemp(:,:))), MAXVAL(ABS(rhs_e(:,jk,:)))
     ELSE
       inv_flip_flop_e2(:,:) = 0.0_wp!rhs_e(:,jk,:)
       CALL gmres( inv_flip_flop_e2(:,:), &! Input is the first guess
                    & lhs_primal_flip_flop,   &! sbr. calculating l.h.s.
                    & h_e,                    &
                    & p_patch,                &! used for calculating l.h.s.
                    & jk,                     &!idx of vertical level
                    & p_patch%edges%in_domain,&
                    & zimpl_coeff,            &! used for calculating l.h.s.
                    & rhs_e(:,jk,:),          &! right hand side as input
                    & tolerance,              &! relative tolerance
                    & .true.,                 & !absolute tolerance
                    & nmax_iter,              &! max. # of iterations to do
                    & lmax_iter,              &! out: .true. = not converged
                    & n_iter,                 &! out: # of iterations done
                    & z_residual              &! out: the residual (array)
                     )

       rhstemp(:,:) = rhs_e(:,jk,:)-lhs_primal_flip_flop(inv_flip_flop_e2(:,:),p_patch, &
         &            jk,zimpl_coeff,h_e)
    !  WRITE(*,*)'max/min residual of inverse primal-flip-flop:',&
    !    &        jk, maxval(rhstemp),minval(rhstemp) 
       idt_src=2  ! output print level (1-5, fix)
       CALL dbg_print('resitual of inv_flip_flop'   ,rhstemp                  ,str_module,idt_src)

       z_e(:,jk,:)=rhstemp(:,:)
       If (maxval (ABS (rhstemp (:,:))) >= tolerance) lmax_iter = .true.
         idt_src=1
         IF (idbg_mxmn >= idt_src) THEN
           IF (lmax_iter) THEN
             WRITE (0, '(1x,a, I4.2, 1x, a,E8.2,1x, a,E8.2,1x, E8.2, 1x, a)') &
             &'Inv_flipflop GMRES #Iter', n_iter, 'Tol ',tolerance, 'Res ',&
             &  ABS(z_residual(n_iter)),MAXVAL (ABS(rhstemp(:,:))), 'GMRES PROBLEM!!!!!!!!!!!!'
           ELSE
             WRITE (0, '(1x,a, I4.2, 1x, a,E8.2,1x, a,E8.2,1x, E8.2)') &
             &'Inv_flipflop GMRES #Iter', n_iter, 'Tol ',tolerance, 'Res ',&
             &  ABS(z_residual(n_iter)),MAXVAL (ABS(rhstemp(:,:)))
           ENDIF
         ENDIF
         inv_flip_flop_e(:,jk,:) = inv_flip_flop_e2(:,:)
      END IF
   END DO

   DEALLOCATE (inv_flip_flop_e2)

   END FUNCTION inverse_primal_flip_flop
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !!  mpi parallelized LL, results is valid only in in_domain edges
   FUNCTION lhs_primal_flip_flop( x, p_patch, jk, p_coeff, h_e) RESULT(llhs)
    !
    TYPE(t_patch), TARGET, INTENT(in) :: p_patch
    REAL(wp),INTENT(inout)       :: x(:,:)
    INTEGER ,INTENT(in)          :: jk
    REAL(wp),INTENT(in)          :: p_coeff
    REAL(wp),OPTIONAL,INTENT(in) :: h_e(SIZE(x,1), SIZE(x,2))!(:,:)
    REAL(wp)                     :: llhs(SIZE(x,1), SIZE(x,2))

    !locl variables
    !INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
    !INTEGER :: rl_start_c, rl_end_c
    !INTEGER :: jc,jb
    !REAL(wp) :: z_x_e(SIZE(x,1),SIZE(x,2))!(nproma,p_patch%nblks_e)
    REAL(wp) :: z_x_out(SIZE(x,1), SIZE(x,2))!(nproma,p_patch%nblks_e)
    !TYPE(t_cartesian_coordinates) :: z_vn_cc(nproma,p_patch%nblks_c)

    !INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    !INTEGER :: il_e, ib_e
    !INTEGER :: ie,je  
    !INTEGER, PARAMETER :: no_cell_edges = 3
    INTEGER :: jb, i_startidx_e, i_endidx_e
    !INTEGER :: rl_start_e, rl_end_e 
    !REAL(wp) :: z_weight
    !REAL(wp) :: z_thick_e
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    edges_in_domain => p_patch%edges%in_domain

    z_x_out(:,:) = 0.0_wp

    CALL sync_patch_array(SYNC_E, p_patch, x)

!     CALL map_edges2edges( p_patch,   &
!                         & x,    &
!                         & z_x_out,   &
!                         & h_e,      &
!                         & level=1 )
  !    rl_start_c = 1
  !    rl_end_c  = min_rlcell
  !    rl_start_e = 1
  !    rl_end_e  = min_rledge
  !    i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
  !    i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
  !    i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
  !    i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
  ! 
  !     DO jb = i_startblk_c, i_endblk_c
  !      CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
  !                        i_startidx_c, i_endidx_c, rl_start_c, rl_end_c) 
  !      DO jc = i_startidx_c, i_endidx_c
  !        z_vn_cc(jc,jb)%x(1:3) = 0.0_wp
  !      END DO
  !    END DO
  !    z_x_e(:,:) = 0.0_wp
  !    z_x_e(:,:) = x(:,:)
  ! !write(*,*)'call with level',jk, x(16,1280), SIZE(x,1), SIZE(x,2)
  ! !write(*,*)'max/min input', maxval(z_x_e(:,:)),minval(z_x_e(:,:)), maxval(x(:,:)),minval(x(:,:)) 
  ! !--------------------------------------------------------------------
  !  CELL_BLK_LOOP: DO jb = i_startblk_c, i_endblk_c
  !     CALL get_indices_c( p_patch, jb,&
  !                       & i_startblk_c, i_endblk_c,&
  !                       & i_startidx_c, i_endidx_c,&
  !                       & rl_start_c, rl_end_c)
  ! 
  !     IF(jk==1)THEN
  !     !We are dealing with the surface layer first
  !     CELL_IDX_LOOP_TOP: DO jc =  i_startidx_c, i_endidx_c
  !       z_weight         = 0.0_wp
  !       z_vn_cc(jc,jb)%x = 0.0_wp
  !       z_thick_e        = 0.0_wp
  !       IF(v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary)THEN
  !         DO ie=1, no_cell_edges
  !           il_e = p_patch%cells%edge_idx(jc,jb,ie)
  !           ib_e = p_patch%cells%edge_blk(jc,jb,ie)
  ! 
  !           z_thick_e = v_base%del_zlev_m(jk) + h_e(il_e,ib_e) 
  !           z_weight = z_weight + v_base%variable_vol_norm(jc,jb,ie)!*z_thick_e
  ! 
  !            z_vn_cc(jc,jb)%x = z_vn_cc(jc,jb)%x&
  !                           & + v_base%edge2cell_coeff_cc(jc,jb,ie)%x&
  !                           & * z_x_e(il_e,ib_e)! * z_thick_e
  !         END DO
  ! 
  !         z_vn_cc(jc,jb)%x = z_vn_cc(jc,jb)%x / z_weight
  !       ELSE
  !        z_vn_cc(jc,jb)%x=0.0_wp 
  !       ENDIF
  !     END DO CELL_IDX_LOOP_TOP
  !     ELSEIF(jk>1)THEN
  !     !Now we calculate at the levels below the surface
  !       CELL_IDX_LOOP: DO jc =  i_startidx_c, i_endidx_c
  !       IF(v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary)THEN
  !         z_vn_cc(jc,jb)%x = 0.0_wp
  !         DO ie=1, no_cell_edges
  !           il_e = p_patch%cells%edge_idx(jc,jb,ie)
  !           ib_e = p_patch%cells%edge_blk(jc,jb,ie)
  !           z_vn_cc(jc,jb)%x = z_vn_cc(jc,jb)%x&
  !                             & + v_base%edge2cell_coeff_cc(jc,jb,ie)%x&
  !                             & * z_x_e(il_e,ib_e)
  !       END DO
  ! 
  !         z_vn_cc(jc,jb)%x = z_vn_cc(jc,jb)%x/v_base%fixed_vol_norm(jc,jb)
  !       ELSE
  !         z_vn_cc(jc,jb)%x = 0.0_wp
  !       ENDIF
  !     END DO CELL_IDX_LOOP
  !     ENDIF
  ! END DO CELL_BLK_LOOP
  ! 
  ! EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e
  ! 
  !   CALL get_indices_e(p_patch, jb,&
  !                    & i_startblk_e, i_endblk_e,&
  !                    & i_startidx_e, i_endidx_e,&
  !                    & rl_start_e, rl_end_e)
  ! 
  !     EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e
  !       IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
  ! 
  !         !Get indices of two adjacent triangles
  !         il_c1 = p_patch%edges%cell_idx(je,jb,1)
  !         ib_c1 = p_patch%edges%cell_blk(je,jb,1)
  !         il_c2 = p_patch%edges%cell_idx(je,jb,2)
  !         ib_c2 = p_patch%edges%cell_blk(je,jb,2)
  ! !write(*,*)'input', je,jb,z_x_e(je,jb), x(je,jb)
  !         z_x_e(je,jb) =&
  !       & DOT_PRODUCT(z_vn_cc(il_c1,ib_c1)%x,v_base%edge2cell_coeff_cc_t(je,jb,1)%x)&
  !       &+DOT_PRODUCT(z_vn_cc(il_c2,ib_c2)%x,v_base%edge2cell_coeff_cc_t(je,jb,2)%x)
  !        ELSE
  !          z_x_e(je,jb) = 0.0_wp
  !        ENDIF
  !     END DO EDGE_IDX_LOOP
  ! END DO EDGE_BLK_LOOP
  !    llhs(1:nproma,1:p_patch%nblks_e) = p_coeff*z_x_e(1:nproma,1:p_patch%nblks_e)

      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        llhs(i_startidx_e:i_endidx_e, jb) = p_coeff * z_x_out(i_startidx_e:i_endidx_e, jb)
      ENDDO
  !write(*,*)'max/min LHS', maxval(llhs(:,:)),minval(llhs(:,:)) 

  END FUNCTION lhs_primal_flip_flop
  !--------------------------------------------------------------------



! ! !------------------------------------------------------------
! ! !
! ! !  
! ! !>
! ! !! Computation of new vertical velocity using continuity equation
! ! 
! ! !! Calculate diagnostic vertical velocity from horizontal velocity using the
! ! !! incommpressibility condition in the continuity equation.
! ! !! For the case of the semi-implicit-AB scheme the land-sea-mask may be applied
! ! !! at least after collecting the whole explicit term.
! ! !! 
! ! !! @par Revision History
! ! !! Developed  by  Peter Korn,   MPI-M (2006).
! ! !!  Modified by Stephan Lorenz, MPI-M (2010-06)
! ! !!  Modified by Stephan Lorenz, MPI-M (2010-08)
! ! !!   - velocities and sea level are passed through the interface
! ! !!   - no calculation of new sea level here
! ! !! 
! ! SUBROUTINE calc_vert_velocity_mimetic_old( p_patch, pvn_e, ph_c, ph_e, bot_bc_w, pw_c )
! ! !
! ! TYPE(t_patch), TARGET, INTENT(IN) :: p_patch       ! patch on which computation is performed
! ! ! #slo# pvn_e is inout since div_oce sets boundary values to zero
! ! REAL(wp),         INTENT(INOUT) :: pvn_e(:,:,:)  ! horizontal velocity on edges
! ! REAL(wp),            INTENT(IN) :: ph_c (:,:)    ! surface elevation on cells 
! ! REAL(wp),            INTENT(IN) :: ph_e (:,:)    ! surface elevation on edges 
! ! REAL(wp),            INTENT(IN) :: bot_bc_w(:,:) !bottom boundary condition for vertical velocity
! ! REAL(wp),         INTENT(INOUT) :: pw_c (:,:,:)  ! vertical velocity on cells
! ! !
! ! !
! ! ! Local variables
! ! !
! ! INTEGER :: jc, jk, jb
! ! !INTEGER :: jic, jib
! ! INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
! ! INTEGER :: rl_start, rl_end
! ! 
! ! REAL(wp) :: delta_z
! ! REAL(wp) :: z_div_c(nproma,n_zlev,p_patch%nblks_c) ! div_oce needs 3 dimensions
! ! ! REAL(wp) :: z_u_c(nproma,n_zlev,p_patch%nblks_c)
! ! ! REAL(wp) :: z_v_c(nproma,n_zlev,p_patch%nblks_c)
! ! ! REAL(wp) :: z_vn_e(nproma,n_zlev,p_patch%nblks_e)
! ! !CHARACTER(len=max_char_length), PARAMETER :: &
! ! !       & routine = ('mo_oce_ab_timestepping_mimetic:calc_vert_velocity_mimetic')
! ! !-----------------------------------------------------------------------  
! ! !CALL message (TRIM(routine), 'start')        
! ! 
! ! ! #slo# due to nag -nan compiler-option:
! ! z_div_c(:,:,:) = 0.0_wp
! ! 
! ! 
! !   !------------------------------------------------------------------
! !   ! Step 1) Calculate divergence of horizontal velocity at level jk
! !   !------------------------------------------------------------------
! !   CALL div_oce( pvn_e, p_patch, z_div_c)
! ! 
! ! !Note we are summing from bottom to top.
! ! LEVEL_LOOP: DO jk = n_zlev, 1, -1
! ! 
! !   !------------------------------------------------------------------
! !   ! Step 2) Use the divergence and the vertical velocity at the previous deeper
! !   !         layer to calculate the new vertical velocity at cell centers
! !   !------------------------------------------------------------------
! ! 
! !   ! Loop over cells
! !   rl_start   = 1
! !   rl_end     = min_rlcell
! !   i_startblk = p_patch%cells%start_blk(rl_start,1)
! !   i_endblk   = p_patch%cells%end_blk(rl_end,1)
! ! 
! !   DO jb = i_startblk, i_endblk
! ! 
! !   CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
! ! 
! !     DO jc = i_startidx, i_endidx
! ! 
! !       IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! ! 
! !         IF (jk == 1)THEN 
! !           ! at surface level (toplev) add surface elevation h_c
! !           ! #slo# bugfix 11-25
! !           !delta_z = delta_z + ph_c(jc,jb)
! !           delta_z = v_base%del_zlev_m(jk) + ph_c(jc,jb)
! !         ELSE
! !           ! depth of cell
! !           delta_z = v_base%del_zlev_m(jk)
! !         ENDIF
! ! 
! !         IF ( jk < v_base%dolic_c(jc,jb) ) THEN
! ! 
! !           ! vertical velocity is integrated from bottom to top and multiplied by
! !           ! depth of the layer
! !           ! #slo# 2010-10-18: vertical velocity is negative for positive divergence
! !           !                   of horizontal velocity
! !           pw_c(jc,jk,jb) = pw_c(jc,jk+1,jb) - z_div_c(jc,jk,jb) * delta_z
! !           !write(*,*)'vert veloc A:',jc,jb, jk,v_base%dolic_c(jc,jb), pw_c(jc,jk,jb)
! !         ELSEIF ( jk == v_base%dolic_c(jc,jb) ) THEN
! !           !use bottom boundary condition for vertical velocity at bottom
! !           !of prism
! !           pw_c(jc,jk,jb) = bot_bc_w(jc,jb) - z_div_c(jc,jk,jb) * delta_z
! ! 
! !         ELSEIF ( jk > v_base%dolic_c(jc,jb) ) THEN
! !           ! Set vertical velocity at center of cell bottom 
! !           ! to zero if vertical layer is below dolic_c
! !           pw_c(jc,jk,jb) = 0.0_wp
! !           !write(*,*)'vert veloc C:',jc,jb, jk,v_base%dolic_c(jc,jb), pw_c(jc,jk,jb)
! !         END IF 
! !       END IF
! !     END DO
! !   END DO
! ! 
! ! END DO LEVEL_LOOP
! ! DO jk = 1,n_zlev
! ! write(*,*)'max/min vert veloc',jk, maxval(pw_c(:,jk,:)), minval(pw_c(:,jk,:))
! ! END DO
! ! 
! !    !jib = i_oct_blk
! !    !jic = i_oct_idx
! !    !98 format(3(a,g25.9))
! !    !  WRITE(message_text,98) &
! !    !    &     'div(jic,lv1) =', z_div_c   (jic,1,jib), &
! !    !    &    '     div(lv2) =', z_div_c   (jic,2,jib), &
! !    !    &    '     div(lv3) =', z_div_c   (jic,3,jib)
! !    !  CALL message (' ', message_text)
! !    !  WRITE(message_text,98) &
! !    !    &     'pwc(jic,lv1) =',    pw_c   (jic,1,jib), &
! !    !    &    '     pwc(lv2) =',    pw_c   (jic,2,jib), &
! !    !    &    '     pwc(lv4) =',    pw_c   (jic,4,jib)
! !    !  CALL message (' ', message_text)
! ! 
! ! !CALL message (TRIM(routine), 'end')        
! ! 
! ! END SUBROUTINE calc_vert_velocity_mimetic_old

END MODULE mo_oce_ab_timestepping_mimetic
