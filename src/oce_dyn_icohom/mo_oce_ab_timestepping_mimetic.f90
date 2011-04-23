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
!-------------------------------------------------------------------------  
!
!    ProTeX FORTRAN source: Style 2  
!    modified for ICON project, DWD/MPI-M 2006
!  
!-------------------------------------------------------------------------  
!  
!   
! 
USE mo_kind,                      ONLY: wp
!USE mo_mpi,                       ONLY: p_pe, p_io
USE mo_run_nml,                   ONLY: nproma
USE mo_math_utilities,            ONLY: t_cartesian_coordinates!, gc2cc
USE mo_impl_constants,            ONLY: sea_boundary,                                     &
  &                                     min_rlcell, min_rledge, min_rlcell,               &
  &                                     max_char_length
USE mo_ocean_nml,                 ONLY: n_zlev, solver_tolerance, L_INVERSE_FLIP_FLOP,    &
                                    &   ab_const, ab_beta, ab_gam, iswm_oce, idisc_scheme,&
                                    &   expl_vertical_velocity_diff
USE mo_run_nml,                   ONLY: dtime
USE mo_dynamics_nml,              ONLY: nold,nnew
USE mo_physical_constants,        ONLY: grav!, re
USE mo_oce_state,                 ONLY: t_hydro_ocean_state, t_hydro_ocean_diag,          &
  &                                     set_lateral_boundary_values
USE mo_model_domain,              ONLY: t_patch
USE mo_oce_linear_solver,         ONLY: gmres_oce, gmres_e2e
USE mo_exception,                 ONLY: message, finish!, message_text
USE mo_loopindices,               ONLY: get_indices_c, get_indices_e !, get_indices_v
USE mo_oce_boundcond,             ONLY: bot_bound_cond_horz_veloc, top_bound_cond_horz_veloc,&
  &                                     bot_bound_cond_vert_veloc, top_bound_cond_vert_veloc
USE mo_oce_thermodyn,             ONLY: calc_density, calc_internal_press
USE mo_oce_physics,               ONLY: t_ho_params
USE mo_oce_forcing,               ONLY: t_ho_sfc_flx, update_ho_sfcflx
USE mo_scalar_product,            ONLY: map_cell2edges, map_edges2cell, map_edges2edges, &
  &                                     calc_scalar_product_for_veloc
USE mo_oce_math_operators,        ONLY: div_oce, grad_fd_norm_oce, grad_fd_norm_oce_2d,  &
  &                                     height_related_quantities

!USE mo_math_gradients,            ONLY: grad_fd_norm
USE mo_oce_veloc_advection,       ONLY: veloc_adv_horz_mimetic, veloc_adv_vert_mimetic

USE mo_interpolation,             ONLY: t_int_state
USE mo_oce_index,                 ONLY: print_mxmn_3d, print_mxmn_2d
USE mo_oce_diffusion,             ONLY: velocity_diffusion_horz_mimetic,&
  &                                     velocity_diffusion_vert_mimetic, &
                                        veloc_diffusion_vert_impl

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
! Private implemenation
!
PRIVATE :: fill_rhs4surface_eq_ab
PRIVATE :: calculate_explicit_term_ab   ! calc_velocity_predictor
PRIVATE :: lhs_surface_height_ab
PRIVATE :: inverse_primal_flip_flop

INTEGER, PARAMETER  :: top=1

LOGICAL, PARAMETER :: l_forc_freshw = .FALSE.
CONTAINS
!-------------------------------------------------------------------------  
!
!  
!>
!! !  Solves the free surface equation.
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
SUBROUTINE solve_free_sfc_ab_mimetic(p_patch, p_os, p_sfc_flx, p_phys_param, timestep, p_int)
!
TYPE(t_patch), TARGET, INTENT(in)             :: p_patch
TYPE(t_hydro_ocean_state), TARGET             :: p_os 
TYPE(t_ho_sfc_flx), INTENT(INOUT)             :: p_sfc_flx
TYPE (t_ho_params)                            :: p_phys_param
INTEGER                                       :: timestep
TYPE(t_int_state),TARGET,INTENT(IN), OPTIONAL :: p_int
!
!Local variables
!
!REAL(wp),ALLOCATABLE :: rhs(:,:)
!REAL(wp) :: aveh, totarea
!REAL(wp), POINTER :: p_height_area(:,:)
!GMRS
INTEGER,PARAMETER :: nmax_iter   = 100      ! maximum number of iterations
REAL(wp) :: tolerance =0.0_wp!1.0E-6            ! (relative or absolute) tolerance
INTEGER  :: n_iter,jk                        ! number of iterations 

REAL(wp) :: z_h_c(nproma,p_patch%nblks_c)
REAL(wp) :: z_h_e(nproma,p_patch%nblks_e)

REAL(wp) :: z_implcoeff
REAL(wp) :: zresidual(nmax_iter)    ! norms of the residual (convergence history); 
                                    ! an argument of dimension at least m is required
LOGICAL  :: l_maxiter                 ! true if reached m iterations
LOGICAL  :: lverbose = .TRUE. 
LOGICAL  :: l_first_timestep 
REAL(wp) :: z_l2
CHARACTER(len=max_char_length) :: string
!INTEGER  :: rl_start_c, rl_end_c, i_startblk_c, i_endblk_c, jc,jb, i_startidx_c, i_endidx_c
!TYPE(t_cartesian_coordinates)  :: z_temp(nproma,p_patch%nblks_c)
!CHARACTER(len=max_char_length), PARAMETER :: &
!       & routine = ('mo_oce_ab_timestepping_mimetic:solve_free_sfc_ab_mimetic')
!-------------------------------------------------------------------------------
!CALL message (TRIM(routine), 'start')
tolerance = solver_tolerance
z_h_c = 0.0_wp
z_h_e = 0.0_wp

CALL height_related_quantities(p_patch, p_os)

IF(timestep==1)THEN
  l_first_timestep = .TRUE.
CALL map_edges2edges( p_patch, p_os%p_prog(nold(1))%vn, p_os%p_prog(nnew(1))%vn,&
                             & p_os%p_diag%h_e,opt_slev=1, opt_elev=1 )
  !This is required in top boundary condition for
  !vertical velocity: the time derivative of the surface height
  !is used there and needs special treatment in the first timestep.
  !see sbr top_bound_cond_vert_veloc in mo_ho_boundcond
  p_os%p_prog(nnew(1))%h=p_os%p_prog(nold(1))%h

  Call set_lateral_boundary_values(p_patch, p_os%p_prog(nold(1))%vn)

  CALL calc_scalar_product_for_veloc( p_patch,                &
                                    & p_os%p_prog(nold(1))%vn,&
                                    & p_os%p_prog(nold(1))%vn,&
                                    & p_os%p_diag%h_e,        &
                                    & p_os%p_diag)

! p_os%p_prog(nold(1))%vn= p_os%p_diag%ptp_vn
!   CALL calc_scalar_product_for_veloc( p_patch,                &
!                                     & p_os%p_prog(nold(1))%vn,&
!                                     & p_os%p_prog(nold(1))%vn,&
!                                     & p_os%p_diag%h_e,        &
!                                     & p_os%p_diag)
ELSE
  l_first_timestep = .FALSE.
ENDIF
write(*,*)'on entry: height:',&
& maxval(p_os%p_prog(nnew(1))%h),minval(p_os%p_prog(nnew(1))%h),&
& maxval(p_os%p_prog(nold(1))%h),minval(p_os%p_prog(nold(1))%h),&
& timestep
write(*,*)'on entry: vn:',&
& maxval(p_os%p_prog(nnew(1))%vn),minval(p_os%p_prog(nnew(1))%vn),&
& maxval(p_os%p_prog(nold(1))%vn),minval(p_os%p_prog(nold(1))%vn)


!P.K.: Currently (Nov 2010) shallow-water means ATMOSPHERIC shallow-water.
!This is sufficient for the purpose of testing. For an oceanic
!shallow-water model including wind-forcing part of the
!boundary conditions below should be included. 
IF ( iswm_oce /= 1 ) THEN

  CALL update_ho_sfcflx(p_patch, p_sfc_flx)

  ! Apply windstress 
  CALL top_bound_cond_horz_veloc(p_patch, p_os, p_phys_param,p_sfc_flx, &
    &                            p_os%p_aux%bc_top_u, p_os%p_aux%bc_top_v,&
    &                            p_os%p_aux%bc_top_veloc_cc)!z_temp  )
!   rl_start_c = 1
!   rl_end_c = min_rlcell
!   i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
!   i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
!   DO jb = i_startblk_c, i_endblk_c
!     CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
!       &                             rl_start_c, rl_end_c)
!    DO jc = i_startidx_c, i_endidx_c
!      !p_os%p_aux%bc_top_veloc_cc(jc,jb)%x = z_temp(jc,jb)%x
!      IF(jb==900)THEN
!        write(*,*)'top bc A:',jc,jb,p_os%p_aux%bc_top_u(jc,jb), p_os%p_aux%bc_top_veloc_cc(jc,jb)%x!z_temp(jc,jb)%x,
!      ENDIF
!    END DO
!   END DO
  ! Apply bot boundary condition for horizontal velocity
  CALL bot_bound_cond_horz_veloc(p_patch, p_os, p_phys_param)

  ! Apply top boundary condition for vertical velocity, this
  !takes height field into account
  CALL top_bound_cond_vert_veloc( p_patch, p_os, p_os%p_aux%bc_top_w,timestep)

  !Apply bottom boundary condition for vertical velocity
  !This boundary condition is time-invariant, to be computed only once 
   IF(l_first_timestep)THEN
     CALL bot_bound_cond_vert_veloc( p_patch, p_os, p_os%p_aux%bc_bot_w )
  ENDIF
ENDIF


! Calculate explicit terms of Adams-Bashforth timestepping
!CALL message (TRIM(routine), 'call calculate_explicit_term_ab')        
CALL calculate_explicit_term_ab(p_patch, p_os, p_phys_param, p_int, l_first_timestep )


! Calculate RHS of surface equation
!CALL message (TRIM(routine), 'call fill_rhs4surface_eq_ab')        
CALL fill_rhs4surface_eq_ab( p_patch, p_os, p_sfc_flx)


! Solve surface equation with GMRES solver
!CALL message (TRIM(routine), 'call GMRES solver')

!use old height as first guess
!z_h_c = p_os%p_prog(nold(1))%h  ! test #slo# 2011-02-21
z_h_c = 0.0_wp                   ! no worse performance

!The lhs needs different thicknesses at edges for 3D and SWE (including bathymetry)
IF( iswm_oce == 1 ) THEN  
  z_h_e = p_os%p_diag%thick_e
ELSEIF( iswm_oce /= 1 ) THEN 
  z_h_e = p_os%p_diag%h_e         ! #slo# 2011-02-21 bugfix (for mimetic only)
ENDIF

CALL gmres_oce( z_h_c,                  &  ! x. Input is the first guess
      &        lhs_surface_height_ab,   &  ! function calculating l.h.s.
      &        z_h_e,                   &
      &        p_os%p_diag%thick_c,     &
      &        p_os%p_prog(nold(1))%h,  &
      &        p_patch,                 & 
      &        p_patch%nblks_c,         &
      &        p_patch%npromz_c,        &
      &        z_implcoeff,             &
      &        p_os%p_aux%p_rhs_sfc_eq, &  ! right hand side as input
      &        tolerance,               &  ! relative tolerance
      &        .FALSE.,                 &  ! NOT absolute tolerance
      &        nmax_iter,               &  ! max. # of iterations to do
      &        l_maxiter,               &  ! out: .true. = not converged
      &        n_iter,                  &  ! out: # of iterations done
      &        zresidual )                 ! inout: the residual (array)  

IF (l_maxiter) THEN
CALL finish('GMRES solver surface equation: ','NOT YET CONVERGED !!')
ELSE
  IF (lverbose) THEN
     WRITE(string,'(a,i4,a,e20.10)') &
     'iteration ', n_iter,', residual = ', ABS(zresidual(n_iter))
     CALL message('GMRES surface height: ',TRIM(string))
     !write(*,*)'residual:',    REAL(zresidual,wp)!write(*,'(1x,e20.10)') zresidual
  ENDIF
ENDIF 

p_os%p_prog(nnew(1))%h = z_h_c!- p_patch%patch_oce%bathymetry_c
z_h_c=0.0_wp
z_h_c = lhs_surface_height_ab( p_os%p_prog(nnew(1))%h,&
                             & p_os%p_prog(nold(1))%h, &
                             & p_patch, &
                             & z_implcoeff,&
                             & p_os%p_diag%thick_e,&
                             & p_os%p_diag%thick_c)&
     & -p_os%p_aux%p_rhs_sfc_eq
 write(*,*)'MIN/MAX:h-residual:',&
 & minval(z_h_c(1:nproma,1:p_patch%nblks_c)),&
 & maxval(z_h_c(1:nproma,1:p_patch%nblks_c)) 


!   rl_start_c = 1
!   rl_end_c = min_rlcell
!   i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
!   i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
!   DO jb = i_startblk_c, i_endblk_c
!     CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
!       &                             rl_start_c, rl_end_c)
!    DO jc = i_startidx_c, i_endidx_c
!      IF(jb==900)THEN
!        write(*,*)'top bc A:',jc,jb,p_os%p_aux%bc_top_veloc_cc(jc,jb)%x
!      ENDIF
! !      write(*,*)'height:old:new:', jc,jb,&
! !      & p_os%p_prog(nold(1))%h(jc,jb),&
! !      & p_os%p_prog(nnew(1))%h(jc,jb)!,&
! !   !   & p_os%p_diag%thick_c(jc,jb)
!    END DO
!   END DO

 write(*,*)'MIN/MAX:h-new:',&
 & minval(p_os%p_prog(nnew(1))%h(1:nproma,1:p_patch%nblks_c)),&
 & maxval(p_os%p_prog(nnew(1))%h(1:nproma,1:p_patch%nblks_c))


!CALL message (TRIM(routine), 'end')
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
SUBROUTINE calculate_explicit_term_ab( p_patch, p_os, p_phys_param, p_int, l_first_timestep)
!
TYPE(t_patch), TARGET, INTENT(in)             :: p_patch
TYPE(t_hydro_ocean_state), TARGET             :: p_os
TYPE (t_ho_params)                            :: p_phys_param
TYPE(t_int_state),TARGET,INTENT(IN), OPTIONAL :: p_int
LOGICAL                                       :: l_first_timestep 
!
!local variables
!
!TYPE(t_ho_params)     :: p_param
INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
INTEGER  :: rl_start, rl_end
INTEGER  :: je, jk, jb, jc
REAL(wp) :: gdt
REAL(wp) :: z_flip_flop_e(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_gradh_e(nproma,1,p_patch%nblks_e)
REAL(wp) :: z_ptp_gradh(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_e(nproma,n_zlev,p_patch%nblks_e)
!TYPE(t_cartesian_coordinates) :: p_vn_c(nproma,n_zlev,p_patch%nblks_c)

INTEGER  :: rl_start_c, rl_end_c, i_startblk_c, i_endblk_c,i_startidx_c, i_endidx_c
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_ab_timestepping_mimetic:calculate_explicit_term_ab')
!-----------------------------------------------------------------------  
!CALL message (TRIM(routine), 'start')        
! #slo# 2011-02-17: move invert of gdt=grav*dt into module for constants (tbd)
gdt=grav*dtime

rl_start   = 1
rl_end     = min_rledge
i_startblk = p_patch%edges%start_blk(rl_start,1)
i_endblk   = p_patch%edges%end_blk(rl_end,1)


!STEP 1: calculate gradient of surface height at previous timestep
!CALL message (TRIM(routine), 'call grad_fd_norm_oce_2d')
CALL grad_fd_norm_oce_2d( p_os%p_prog(nold(1))%h, &
       &                  p_patch,                &
       &                  z_gradh_e(:,1,:))

! STEP 2: calculate nonlinear Coriolis term.
! horiz. advcetion= gradient kinetic energy + nonlinear Coriolis term
!(nonlinear Coriolis term corresponds to (Coriolis param + curl of veloc)*tangential veloc).
!
!STEP 2: horizontal advection
!CALL message (TRIM(routine), 'call veloc_adv_horz')
CALL veloc_adv_horz_mimetic( p_patch,                    &
       &             p_os%p_prog(nold(1))%vn,    &
       &             p_os%p_diag,                &
       &             p_os%p_diag%veloc_adv_horz, & !contains nonlinear Coriolis term, gradient kinetic energy is calculated seperately
       &             p_int                      )


! STEP 3: compute 3D contributions: gradient of hydrostatic pressure and vertical velocity advection
IF ( iswm_oce /= 1 ) THEN
  ! calculate density from EOS using temperature and salinity at timelevel n
  !CALL message (TRIM(routine), 'call calc_density')
  CALL calc_density( p_patch,                                 &
         &           p_os%p_prog(nold(1))%tracer(:,:,:,1:2),  &
         &           p_os%p_diag%rho(:,:,:) )

  ! calculate hydrostatic pressure from density at timelevel nc
  !CALL message (TRIM(routine), 'call calc_internal_pressure')
  CALL calc_internal_press( p_patch,  p_phys_param, &
         &                  p_os%p_diag%rho,       &
         &                  p_os%p_prog(nold(1))%h,&
         &                  p_os%p_diag%press_hyd)

!    rl_start_c = 1
!    rl_end_c = min_rlcell
!    i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
!    i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
!    DO jk=1,n_zlev
!      DO jb = i_startblk_c, i_endblk_c
!        CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
!          &                i_startidx_c, i_endidx_c, &
!          &                rl_start_c, rl_end_c)
!         DO jc = i_startidx_c, i_endidx_c
!           p_vn_c(jc,jk,jb)%x = 0.0_wp
!         END DO
!      END DO
!    ENDDO
  ! calculate gradient of hydrostatic pressure in 3D
  !CALL message (TRIM(routine), 'call grad_fd_norm_oce')
  CALL grad_fd_norm_oce( p_os%p_diag%press_hyd,  &
         &               p_patch,                &
         &               p_os%p_diag%press_grad)

  !CALL message (TRIM(routine), 'call veloc_adv_vert')
  CALL veloc_adv_vert_mimetic( p_patch,          &
       &             p_os%p_aux,                 &
       &             p_os%p_diag,                &
       &             p_os%p_aux%bc_top_w,        &
       &             p_os%p_aux%bc_bot_w,        &
       &             p_os%p_diag%veloc_adv_vert )

  IF(expl_vertical_velocity_diff==0)THEN
  !For alternative "expl_vertical_velocity_diff==1" see couples of
  !lines below below
      CALL velocity_diffusion_vert_mimetic( p_patch,        &
      &                             p_os%p_diag,            &
      &                             p_os%p_aux,             &
      &                             p_os%p_prog(nold(1))%h, &
      &                             p_phys_param,&
      &                             p_os%p_diag%laplacian_vert)
  ENDIF

  WRITE(*,*)'MAX/MIN old height gradient:', MAXVAL(z_gradh_e(:,1,:)),&
  & MINVAL(z_gradh_e(:,1,:)) 
  CALL print_mxmn_2d (1,z_gradh_e(:,1,:),p_patch%nblks_e,'old height gradient')
  DO jk=1, n_zlev
    WRITE(*,*)'MAX/MIN density:',       jk, &
      &        MAXVAL(p_os%p_diag%rho(:,jk,:)),&
      &        MINVAL(p_os%p_diag%rho(:,jk,:)) 
  END DO
  DO jk=1, n_zlev
    WRITE(*,*)'MAX/MIN internal press:',jk, &
      &        MAXVAL(p_os%p_diag%press_hyd(:,jk,:)),&
      &        MINVAL(p_os%p_diag%press_hyd(:,jk,:)) 
  END DO
  DO jk=1, n_zlev
    WRITE(*,*)'MAX/MIN internal press grad:',    jk, &
      &        MAXVAL(p_os%p_diag%press_grad(:,jk,:)),&
      &        MINVAL(p_os%p_diag%press_grad(:,jk,:)) 
  END DO

  DO jk=1, n_zlev
    write(*,*)'MAX/MIN kin energy:', jk,MAXVAL(p_os%p_diag%kin(:,jk,:)),&
                                       &MINVAL(p_os%p_diag%kin(:,jk,:)) 
  END DO

 DO jk=1, n_zlev
    write(*,*)'MAX/MIN vert advection:', jk,MAXVAL(p_os%p_diag%veloc_adv_vert(:,jk,:)),&
                                       &MINVAL(p_os%p_diag%veloc_adv_vert(:,jk,:)) 
  END DO
ELSEIF( iswm_oce == 1 ) THEN
  p_os%p_diag%press_grad     = 0.0_wp
  p_os%p_diag%veloc_adv_vert = 0.0_wp
  p_os%p_diag%laplacian_vert = 0.0_wp
ENDIF

! STEP 3: compute laplacian diffusion of velocity
! calculate horizontal laplacian of horizontal velocity, provided this term is discretized explicitly
! ! CALL message (TRIM(routine), 'call laplacian diffusion')        
! ! IF(horizontal_diffusion_veloc==EXPLICIT)THEN
    CALL velocity_diffusion_horz_mimetic(p_patch,&
                                       & p_os%p_prog(nold(1))%vn,&
                                       & p_phys_param,&
                                       & p_os%p_diag,&
                                       & p_os%p_diag%laplacian_horz)
DO jk=1, n_zlev
  write(*,*)'MAX/MIN horz diffusion:', jk,MAXVAL(p_os%p_diag%laplacian_horz(:,jk,:)),&
                                         &MINVAL(p_os%p_diag%laplacian_horz(:,jk,:)) 
END DO

! Put it all together to obtain explicit term at intermediate timelevel
IF(L_INVERSE_FLIP_FLOP)THEN

  z_e(:,:,:)            = 0.0_wp
  z_ptp_gradh(:,1,:)    = 0.0_wp
  z_flip_flop_e(:,:,:)  = 0.0_wp

  ! Apply inverse of flip-flop to g_nimd. The results is a edge variable
  !Commented out to gain speed, with GMRS-inversion at each vertical level code 
  !is too slow, but actually no time for thoughts about speeding up this part. P.K.
  !CALL message (TRIM(routine), 'call inverse_flip_flop')
   z_e(:,:,:)= p_os%p_diag%press_grad(:,:,:)!&
            !& +p_os%p_diag%grad

  CALL map_edges2edges( p_patch,                &
                     & z_e,                     &
                     & p_os%p_diag%press_grad,  &
                     & p_os%p_diag%h_e )

 p_os%p_aux%g_n(:,:,:) =- p_os%p_diag%press_grad(:,:,:)      &
     &                  - p_os%p_diag%veloc_adv_horz(:,:,:)  &
     &                  - p_os%p_diag%veloc_adv_vert(:,:,:)  &
     &                   + p_os%p_diag%laplacian_horz(:,:,:) &
     &                   + p_os%p_diag%laplacian_vert(:,:,:)

  z_ptp_gradh= 0.0_wp
  CALL map_edges2edges( p_patch,       &
                     & z_gradh_e,      &
                     & z_ptp_gradh,    &
                     & p_os%p_diag%h_e,&
                     & opt_slev=1, opt_elev=1)


  IF(l_first_timestep)THEN
    p_os%p_aux%g_nimd(:,:,:) = p_os%p_aux%g_n(:,:,:)
  ELSE
    p_os%p_aux%g_nimd(:,:,:) = (1.5_wp+AB_const)* p_os%p_aux%g_n(:,:,:)   &
    &                        - (0.5_wp+AB_const)* p_os%p_aux%g_nm1(:,:,:)
  ENDIF


  p_os%p_diag%vn_pred(:,:,:) =  p_os%p_diag%ptp_vn(:,:,:)&
                          & -dtime*(grav* (1.0_wp-ab_beta)*z_ptp_gradh(:,:,:)&
                          &       + p_os%p_aux%g_nimd(:,:,:))
!   !p_os%p_diag%vn_pred = p_os%p_prog(nold(1))%vn !p_os%p_diag%ptp_vn
!   DO jk=1, n_zlev
!     write(*,*)'MAX/MIN before INVERSE FLIPFLOP:', jk,MAXVAL(p_os%p_diag%vn_pred(:,jk,:)),&
!                                                     &MINVAL(p_os%p_diag%vn_pred(:,jk,:)) 
!   END DO
!   
! 
!   p_os%p_aux%g_nm1 = inverse_primal_flip_flop(p_patch, p_os%p_diag%vn_pred, p_os%p_diag%h_e)
! 
!   DO jk=1, n_zlev
!     write(*,*)'MAX/MIN AFTER INVERSE FLIPFLOP:', jk,MAXVAL(p_os%p_aux%g_nm1(:,jk,:)),&
!                                                  &  MINVAL(p_os%p_aux%g_nm1(:,jk,:)) 
!   END DO
ELSEIF(.NOT.(L_INVERSE_FLIP_FLOP))THEN
! DO jk = 1, n_zlev
!     write(*,*)'MIN/MAX before dual-flip-flop:',jk,minval(p_os%p_diag%veloc_adv_horz(:,jk,:) ),&
!                                    maxval(p_os%p_diag%veloc_adv_horz(:,jk,:) ) 
!   END DO
!   z_e = inverse_primal_flip_flop(p_patch, p_os%p_diag%veloc_adv_horz, p_os%p_diag%thick_e)
!   p_os%p_diag%veloc_adv_horz = z_e
! 
! DO jk = 1, n_zlev
!     write(*,*)'MIN/MAX after dual-flip-flop:',jk,minval(p_os%p_diag%veloc_adv_horz(:,jk,:) ),&
!                                    maxval(p_os%p_diag%veloc_adv_horz(:,jk,:) ) 
!   END DO

  p_os%p_aux%g_n(:,:,:) =-p_os%p_diag%press_grad(:,:,:)      &
    &                   - p_os%p_diag%veloc_adv_horz(:,:,:)  &
    &                   - p_os%p_diag%veloc_adv_vert(:,:,:)  &
    &                   + p_os%p_diag%laplacian_horz(:,:,:)  &
    &                   + p_os%p_diag%laplacian_vert(:,:,:)


  IF(l_first_timestep)THEN
    p_os%p_aux%g_nimd(:,:,:) = p_os%p_aux%g_n(:,:,:)
  ELSE
    p_os%p_aux%g_nimd(:,:,:) = (1.5_wp+AB_const)* p_os%p_aux%g_n(:,:,:)   &
    &                        - (0.5_wp+AB_const)* p_os%p_aux%g_nm1(:,:,:)
  ENDIF

  DO jb = i_startblk, i_endblk
    CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = 1, n_zlev
      DO je = i_startidx, i_endidx
        IF(p_patch%patch_oce%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN

         p_os%p_diag%vn_pred(je,jk,jb) = p_os%p_prog(nold(1))%vn(je,jk,jb)       &
           &                           + dtime*(p_os%p_aux%g_nimd(je,jk,jb)     &
           &                           - (1.0_wp-ab_beta) * grav*z_gradh_e(je,1,jb))
        ELSE
          p_os%p_diag%vn_pred(je,jk,jb) = 0.0_wp
        ENDIF
      END DO
    END DO
  END DO
  IF(expl_vertical_velocity_diff==1)THEN

    CALL veloc_diffusion_vert_impl( p_patch,                  &
                                  & p_os%p_diag%vn_pred,      &
                                  & p_os%p_aux,               &
                                  & p_phys_param%A_veloc_v,   &
                                  & p_os%p_diag%vn_impl_vert_diff)
  ENDIF

  DO jk = 1, n_zlev
    write(*,*)'min/max vn_pred:',jk,     minval(p_os%p_diag%vn_pred(:,jk,:)), &
      &                                  maxval(p_os%p_diag%vn_pred(:,jk,:))

    write(*,*)'min/max vn_pred Term 1:',jk,     minval(p_os%p_aux%g_nimd(:,jk,:)), &
      &                                         maxval(p_os%p_aux%g_nimd(:,jk,:))

    IF(ab_beta/=1.0_wp)THEN
      write(*,*)'min/max vn_pred Term 2:',jk,     minval((1.0_wp-ab_beta) * grav*z_gradh_e), &
        &                                         maxval((1.0_wp-ab_beta) * grav*z_gradh_e)
    ENDIF
    write(*,*)'min/max vn_old :',jk,  minval(p_os%p_prog(nold(1))%vn(:,jk,:)), &
      &                               maxval(p_os%p_prog(nold(1))%vn(:,jk,:))

    write(*,*)'min/max G_n+1/2:',jk,  minval(p_os%p_aux%g_nimd(:,jk,:)),&
                                    & maxval(p_os%p_aux%g_nimd(:,jk,:))

    write(*,*)'min/max G_n:',    jk,  minval(p_os%p_aux%g_n(:,jk,:)),&
                                    & maxval(p_os%p_aux%g_n(:,jk,:))

    write(*,*)'min/max G_n-1:',  jk,  minval(p_os%p_aux%g_nm1(:,jk,:)),&
                                    & maxval(p_os%p_aux%g_nm1(:,jk,:))

    write(*,*)'min/max adv:',    jk,  minval(p_os%p_diag%veloc_adv_horz(:,jk,:)), &
      &                               maxval(p_os%p_diag%veloc_adv_horz(:,jk,:))
    IF(expl_vertical_velocity_diff==1)THEN
      write(*,*)'min/max adv:',    jk,  minval(p_os%p_diag%vn_impl_vert_diff(:,jk,:)), &
      &                                 maxval(p_os%p_diag%vn_impl_vert_diff(:,jk,:))
    ENDIF
    write(*,*)
  END DO
  write(*,*)'min/max -dt*g*grad_h:',minval(dtime*grav*z_gradh_e), maxval(dtime* grav*z_gradh_e)
ENDIF!(L_INVERSE_FLIP_FLOP)

!CALL message (TRIM(routine), 'end')        
END SUBROUTINE calculate_explicit_term_ab
!-------------------------------------------------------------------------  
!
!  
!>
!! 
!!  Calculation of right-hand side of elliptic surface equation.
!!  This is used in semi implicit timelevel stepping. 
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
SUBROUTINE fill_rhs4surface_eq_ab( p_patch, p_os, p_sfc_flx)
!
! Patch on which computation is performed
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
!
! Type containing ocean state
TYPE(t_hydro_ocean_state), TARGET :: p_os
TYPE(t_ho_sfc_flx), INTENT(IN)    :: p_sfc_flx
!
!  local variables
!
INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
INTEGER :: i_startblk_e, i_endblk_e!, i_startidx_e, i_endidx_e
INTEGER :: rl_start_c, rl_end_c
INTEGER :: rl_start_e, rl_end_e
INTEGER :: jc, jb, jk!, je
REAL(wp) :: z_e(nproma,1,p_patch%nblks_e)
REAL(wp) :: gdt2
REAL(wp) :: div_z_c(nproma,1,p_patch%nblks_c) 
REAL(wp) :: z_vn_ab(nproma,n_zlev,p_patch%nblks_e)
TYPE(t_cartesian_coordinates) :: z_u_pred_cc(nproma,n_zlev,p_patch%nblks_c)
TYPE(t_cartesian_coordinates) :: z_u_pred_depth_int_cc(nproma,1,p_patch%nblks_c)
!REAL(wp) :: thick
!CHARACTER(len=max_char_length), PARAMETER :: &
!       & routine = ('mo_oce_ab_timestepping_mimetic:fill_rhs4surface_eq_ab')
!-------------------------------------------------------------------------------
!CALL message (TRIM(routine), 'start')        

! #slo# 2011-02-17: move invert of gdt=grav*dt into module for constants (tbd)
gdt2 = grav*(dtime)**2

rl_start_c = 1
rl_end_c = min_rlcell

rl_start_e = 1
rl_end_e   = min_rledge

i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c( p_patch, jb,             &
                   &  i_startblk_c, i_endblk_c,&
                   &  i_startidx_c, i_endidx_c,&
                   &  rl_start_c, rl_end_c)
  DO jc = i_startidx_c, i_endidx_c
    z_u_pred_depth_int_cc(jc,1,jb)%x(:) = 0.0_wp
    z_u_pred_cc(jc,1,jb)%x(:)           =0.0_wp
  END DO
END DO
div_z_c(:,:,:)=0.0_wp



IF(expl_vertical_velocity_diff==0)THEN
  z_vn_ab = ab_gam*p_os%p_diag%vn_pred + (1.0_wp -ab_gam)* p_os%p_prog(nold(1))%vn
ELSEIF(expl_vertical_velocity_diff==1)THEN
  z_vn_ab = ab_gam*p_os%p_diag%vn_impl_vert_diff + (1.0_wp -ab_gam)* p_os%p_prog(nold(1))%vn
ENDIF

! DO jk = 1, n_zlev
!   write(*,*)'MAX/MIN z_vn_ab:', jk,maxval(z_vn_ab(:,jk,:)),minval(z_vn_ab(:,jk,:))
! END DO
!Step 1) Do within each layer a edge to cell mapping of vn_pred
!For below-surface cells, no height has to be provided, and the reconstructions
!are normalized by cell area (see primal flip-flop).
!For surface cells, the working hypothesis is to do it as in the SW-code, i.e.
!we take also the height at the edges into account:
!
IF( iswm_oce /= 1 ) THEN !the 3D case
  CALL map_edges2cell( p_patch,   &
                    & z_vn_ab,    &
                    & z_u_pred_cc,&
                    & p_os%p_diag%h_e )

!calculate depth-integarted velocity 
  DO jb = i_startblk_c, i_endblk_c

    CALL get_indices_c( p_patch, jb,             &
                     &  i_startblk_c, i_endblk_c,&
                     &  i_startidx_c, i_endidx_c,&
                     &  rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      DO jk=1,n_zlev
        z_u_pred_depth_int_cc(jc,1,jb)%x = z_u_pred_depth_int_cc(jc,1,jb)%x&
                                        &+ z_u_pred_cc(jc,jk,jb)%x&
                                        &* p_patch%patch_oce%del_zlev_m(jk)
      END DO

     !Add surface elevation
     !For a linear surface this can be eleminated
     z_u_pred_depth_int_cc(jc,1,jb)%x = z_u_pred_depth_int_cc(jc,1,jb)%x&
                                    &+ z_u_pred_cc(jc,1,jb)%x&
                                    &* p_os%p_prog(nold(1))%h(jc,jb)
    ENDDO
  END DO

!write(*,*)'MAX/MIN z_e 1:', maxval(p_os%p_diag%u_pred(:,1,:)),minval(p_os%p_diag%u_pred(:,1,:)), &
!& maxval(p_os%p_diag%v_pred(:,1,:)),minval(p_os%p_diag%v_pred(:,1,:))


ELSEIF( iswm_oce == 1 ) THEN !the shallow-water case
 CALL map_edges2cell( p_patch,   &
                   & z_vn_ab,    &
                   & z_u_pred_cc,&
                   & p_os%p_diag%thick_e )

DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c( p_patch, jb,&
                    & i_startblk_c, i_endblk_c,&
                    & i_startidx_c, i_endidx_c,&
                    & rl_start_c, rl_end_c)
  DO jc = i_startidx_c, i_endidx_c

    z_u_pred_depth_int_cc(jc,1,jb)%x = z_u_pred_cc(jc,1,jb)%x&
                                    &* p_os%p_diag%thick_c(jc,jb)
  ENDDO
END DO

!z_e(:,1,:) = z_vn_ab(:,1,:)*p_os%p_diag%thick_e(:,:)
ENDIF

 CALL map_cell2edges( p_patch,              &
                    & z_u_pred_depth_int_cc,&
                    & z_e,                  &
                    & opt_slev=1, opt_elev=1)



write(*,*)'MAX/MIN depth_int_cc:',1,maxval(z_u_pred_depth_int_cc(:,1,:)%x(1)), &
  &  minval(z_u_pred_depth_int_cc(:,1,:)%x(1))
! DO jk = 1, n_zlev
! write(*,*)'MAX/MIN z_e:', jk,maxval(p_os%p_diag%u_pred(:,jk,:)),minval(p_os%p_diag%u_pred(:,jk,:)), &
! & maxval(p_os%p_diag%v_pred(:,jk,:)),minval(p_os%p_diag%v_pred(:,jk,:))
! END DO
! write(*,*)'MAX/MIN depth integrated mass vector:', maxval(z_u_pred_depth_int(:,1,:)),minval(z_u_pred_depth_int(:,1,:)), &
! & maxval(z_v_pred_depth_int(:,1,:)),minval(z_v_pred_depth_int(:,1,:))
!write(*,*)'MAX/MIN depth integrated mass flux:', maxval(z_e(:,1,:)),minval(z_e(:,1,:)) 

! Step 4) Calculate divergence of predicted horizontal velocity
! #slo# due to nag -nan compiler-option:
!  - set divergence to zero otherwise last block contains undef values
!div_z_c(:,1,:) = 0.0_wp

CALL div_oce( z_e, p_patch, div_z_c, 1,1 ) ! to be included surface forcing* +dtime*(P_E)

IF(l_forc_freshw)THEN

  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                     & i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      IF(p_patch%patch_oce%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
        p_os%p_aux%p_rhs_sfc_eq(jc,jb) = (p_os%p_prog(nold(1))%h(jc,jb)&
                                       & - dtime*(div_z_c(jc,1,jb)     &
                                       & + p_sfc_flx%forc_freshw(jc,jb) ))/gdt2
       ELSE
         p_os%p_aux%p_rhs_sfc_eq(jc,jb) = 0.0_wp
       ENDIF
       !write(*,*)'RHS:',jc,jb,p_os%p_aux%p_rhs_sfc_eq(jc,jb), div_z_c(jc,1,jb)
    ENDDO
  END DO


ELSEIF(.NOT.l_forc_freshw)THEN

  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                     & i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      IF(p_patch%patch_oce%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
        p_os%p_aux%p_rhs_sfc_eq(jc,jb) = (p_os%p_prog(nold(1))%h(jc,jb)&
                                       & - dtime*div_z_c(jc,1,jb))/gdt2
       ELSE
         p_os%p_aux%p_rhs_sfc_eq(jc,jb) = 0.0_wp
       ENDIF
       !write(*,*)'RHS:',jc,jb,p_os%p_aux%p_rhs_sfc_eq(jc,jb), div_z_c(jc,1,jb)
    ENDDO
  END DO

ENDIF

write(*,*)'MAX/MIN thick_e:', maxval(p_os%p_diag%thick_e),&
&minval(p_os%p_diag%thick_e) 
write(*,*)'MAX/MIN RHS z_e:', maxval(z_e(:,1,:)),&
&minval(z_e(:,1,:)) 
write(*,*)'MAX/MIN div_c:', maxval(div_z_c(:,1,:)),&
&minval(div_z_c(:,1,:)) 
write(*,*)'MAX/MIN RHS:', maxval(p_os%p_aux%p_rhs_sfc_eq(:,:)),&
&minval(p_os%p_aux%p_rhs_sfc_eq(:,:)) 
!stop

!CALL message (TRIM(routine), 'end')        
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
!-------------------------------------------------------------------------

FUNCTION lhs_surface_height_ab( p_x, h_old, p_patch, coeff, h_e, thickness_c) RESULT(p_lhs)
!
REAL(wp),    INTENT(IN)   :: p_x(:,:)         ! dimension: (nproma,p_patch%nblks_c)
REAL(wp),    INTENT(IN)   :: h_old(:,:)
TYPE(t_patch), INTENT(in) :: p_patch          ! patch on which computation is performed
REAL(wp),    INTENT(in)   :: coeff
REAL(wp),    INTENT(in)   :: h_e(:,:)         !SW-case: thickness at edges
                                              !3D-case: surface height above zero at edges
REAL(wp),    INTENT(in)   :: thickness_c(:,:) !thickness of fluid column    
!
! Left-hand side calculated from iterated height
!
REAL(wp) :: p_lhs(SIZE(p_x,1), SIZE(p_x,2))  ! (nproma,p_patch%nblks_c)
!
! local variables
REAL(wp) :: gdt2
REAL(wp) :: z_grad_h(nproma,1,p_patch%nblks_e)
REAL(wp) :: z_e(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_e2(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_u_c(nproma,1,p_patch%nblks_c)
REAL(wp) :: z_v_c(nproma,1,p_patch%nblks_c)
REAL(wp) :: div_z_c(SIZE(p_x,1), 1, SIZE(p_x,2))  ! (nproma,1,p_patch%nblks_c)
TYPE(t_cartesian_coordinates) :: z_grad_h_cc(nproma,1,p_patch%nblks_c)
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
INTEGER :: rl_start, rl_end
INTEGER :: jc, jb
! CHARACTER(len=max_char_length), PARAMETER ::     &
!   &      routine = ('mo_oce_ab_timestepping_mimetic: lhs_surface_height_ab')
!-----------------------------------------------------------------------  
!CALL message (TRIM(routine), 'start - iteration by GMRES')        

gdt2 = grav*(dtime)**2
z_u_c = 0.0_wp
z_v_c = 0.0_wp
div_z_c = 0.0_wp
z_e2(:,:,:)=0.0_wp

rl_start = 1
rl_end   = min_rlcell
i_startblk = p_patch%cells%start_blk(rl_start,1)
i_endblk   = p_patch%cells%end_blk(rl_end,1)

!Step 1) Calculate gradient of iterated height.
CALL grad_fd_norm_oce_2D( p_x, p_patch, z_grad_h(:,1,:))

!Step 2) map the gradient to the cell center, multiply it
!by fluid thickness and map the result back to edges 
IF( iswm_oce /= 1 ) THEN !the 3D case

 CALL map_edges2cell( p_patch,           &
                   & z_grad_h,           &
                   & z_grad_h_cc,        &
                   & h_e,                &
                   & opt_slev=top,opt_elev=top )

ELSEIF( iswm_oce == 1 ) THEN !the shallow-water case

 CALL map_edges2cell( p_patch,           &
                   & z_grad_h,           &
                   & z_grad_h_cc,        &
                   & h_e,                &
                   & opt_slev=top,opt_elev=top )

ENDIF
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
  DO jc = i_startidx, i_endidx
    z_grad_h_cc(jc,1,jb)%x = z_grad_h_cc(jc,1,jb)%x *thickness_c(jc,jb)
  END DO
END DO

CALL map_cell2edges( p_patch,    &
                   & z_grad_h_cc,&
                   & z_e,        &
                   & opt_slev=1, opt_elev=1)

!z_e(:,1,:)=z_grad_h(:,1,:)*h_e

!write(*,*)'LHS:max/min :', maxval(z_e(:,1,:)), minval(z_e(:,1,:))
!  rl_start = 1
!  rl_end = min_rledge
!  i_startblk = p_patch%edges%start_blk(rl_start,1)
!  i_endblk   = p_patch%edges%end_blk(rl_end,1)
!  DO jb = i_startblk, i_endblk
!     CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
!   DO jc = i_startidx, i_endidx
!   write(*,*)'z_e:', jc,jb,z_e(jc,1,jb), z_e2(jc,1,jb)
!     END DO
!   END DO
! stop
!Step 3) Calculate divergence
CALL div_oce( z_e, p_patch, div_z_c, top,top )

!write(*,*)'div', div_z_c(:,1,:)

!Step 4) Finalize LHS calculations
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
  DO jc = i_startidx, i_endidx
    IF(p_patch%patch_oce%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
      p_lhs(jc,jb) = (p_x(jc,jb)- gdt2*ab_gam*ab_beta*div_z_c(jc,1,jb))/gdt2
   ELSE
      p_lhs(jc,jb) =0.0_wp
   ENDIF
  END DO
END DO
! p_lhs(:,:) = (p_x(:,:)- gdt2*ab_gam*ab_beta*div_z_c(:,1,:))/gdt2

!  rl_start = 1
!  rl_end = min_rlcell
!  i_startblk = p_patch%cells%start_blk(rl_start,1)
!  i_endblk   = p_patch%cells%end_blk(rl_end,1)
!  DO jb = i_startblk, i_endblk
!     CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
!   DO jc = i_startidx, i_endidx
! write(*,*)'P-lhs',jc,jb,p_lhs(jc,jb)
! ! ! write(*,*)'thickness', jc,jb,thickness_c(jc,jb)
! !  write(*,*)'vectors', jc,jb,z_u_c(jc,1,jb), z_v_c(jc,1,jb)
!     END DO
!   END DO
!  write(*,*)'max/min LHS z_e', maxval(z_e),minval(z_e)  
!  write(*,*)'max/min LHS h_e', maxval(h_e),minval(h_e) 
!  write(*,*)'max/min LHS',     maxval(p_lhs),minval(p_lhs) 


END FUNCTION lhs_surface_height_ab
!-------------------------------------------------------------------------  
!
!  
!>
!! Computation of new velocity in Adams-Bashforth timestepping.
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
SUBROUTINE calc_normal_velocity_ab_mimetic(p_patch, p_os)
!
! Patch on which computation is performed
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
!
! Type containing ocean state
TYPE(t_hydro_ocean_state), TARGET :: p_os
!
!
!  local variables
!
INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
INTEGER :: rl_start_e, rl_end_e
INTEGER :: je, jk, jb
REAL(wp) :: z_grad_h(nproma,1,p_patch%nblks_e)
REAL(wp) :: gdt
TYPE(t_cartesian_coordinates) :: z_p_cc(nproma,1,p_patch%nblks_c)
!REAL(wp) :: z_tmp_1_h_c(nproma,1,p_patch%nblks_c)
REAL(wp) :: z_tmp_h_e(nproma,1,p_patch%nblks_e)
!CHARACTER(len=max_char_length), PARAMETER ::     &
!  &      routine = ('mo_oce_ab_timestepping_mimetic: calc_normal_velocity_ab_mimetic')
!-----------------------------------------------------------------------  
!CALL message (TRIM(routine), 'start')        
gdt=grav*dtime

rl_start_e = 1
rl_end_e = min_rledge
i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

! Step 1) Compute normal derivative of new surface height
 CALL grad_fd_norm_oce_2D(p_os%p_prog(nnew(1))%h, &
   &                      p_patch,                &
   &                      z_grad_h(:,1,:))

IF(L_INVERSE_FLIP_FLOP)THEN


  CALL map_edges2edges( p_patch,       &
                     & z_grad_h,       &
                     & z_tmp_h_e,      &
                     & p_os%p_diag%h_e,&
                     & opt_slev=1, opt_elev=1 )


  DO jb = i_startblk_e, i_endblk_e
    CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, &
                      &i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
#ifdef __SX__
!CDIR UNROLL=6
#endif

    DO jk = 1, n_zlev
      DO je = i_startidx_e, i_endidx_e
          p_os%p_diag%vn_pred(je,jk,jb) = p_os%p_diag%vn_pred(je,jk,jb)&
                                        & - gdt*ab_beta*z_tmp_h_e(je,1,jb)
      END DO 
    END DO
  END DO
  DO jk = 1, n_zlev
    write(*,*)'MIN/MAX before dual-flip-flop:',jk,minval(p_os%p_diag%vn_pred(:,jk,:) ),&
                                   maxval(p_os%p_diag%vn_pred(:,jk,:) ) 
  END DO

  p_os%p_prog(nnew(1))%vn = inverse_primal_flip_flop(p_patch, p_os%p_diag%vn_pred, p_os%p_diag%h_e)

  DO jk = 1, n_zlev
    write(*,*)'MIN/MAX after dual-flip-flop:',jk,minval(p_os%p_prog(nnew(1))%vn(:,jk,:) ),&
                                   maxval(p_os%p_prog(nnew(1))%vn(:,jk,:) ) 
  END DO


ELSEIF(.NOT.L_INVERSE_FLIP_FLOP)THEN
! Step 2) Calculate the new velocity from the predicted one and the new surface height

  DO jb = i_startblk_e, i_endblk_e
    CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, &
                      &i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
#ifdef __SX__
!CDIR UNROLL=6
#endif

    DO jk = 1, n_zlev
      DO je = i_startidx_e, i_endidx_e
          p_os%p_prog(nnew(1))%vn(je,jk,jb) = p_os%p_diag%vn_pred(je,jk,jb)  &
            &                               - gdt*ab_beta*z_grad_h(je,1,jb)

      END DO 
    END DO
  END DO
ENDIF




write(*,*)'MIN/MAX new height gradient:', minval(z_grad_h),&
                               & maxval(z_grad_h) 
write(*,*)'MIN/MAX h-contrib to veloc:',minval(- ab_beta*gdt*z_grad_h),&
                                maxval(- ab_beta*gdt*z_grad_h) 

DO jk = 1, n_zlev
  write(*,*)'MIN/MAX vn old:',jk,minval(p_os%p_prog(nold(1))%vn(:,jk,:) ),&
                                 maxval(p_os%p_prog(nold(1))%vn(:,jk,:) ) 
END DO
DO jk = 1, n_zlev
  write(*,*)'MIN/MAX vn new:',jk,minval(p_os%p_prog(nnew(1))%vn(:,jk,:) ),&
                                 maxval(p_os%p_prog(nnew(1))%vn(:,jk,:) )!, p_os%p_prog(nnew(1))%vn(1:10,jk,1) 
END DO

DO jk = 1, n_zlev
  write(*,*)'MIN/MAX vn change:',jk,&
  &minval(p_os%p_prog(nnew(1))%vn(:,jk,:)-p_os%p_prog(nold(1))%vn(:,jk,:) ),&
  &maxval(p_os%p_prog(nnew(1))%vn(:,jk,:)-p_os%p_prog(nold(1))%vn(:,jk,:))
END DO


! Step 3) Update explicit term for Adams-Bashforth stepping 
 p_os%p_aux%g_nm1 = p_os%p_aux%g_n
 p_os%p_aux%g_n   = 0.0_wp


  !Update of scalar product quantities
  !CALL height_related_quantities(p_patch, p_os)
  CALL calc_scalar_product_for_veloc( p_patch,                &
                                    & p_os%p_prog(nold(1))%vn,&
                                    & p_os%p_prog(nnew(1))%vn,&
                                    & p_os%p_diag%h_e,        &
                                    & p_os%p_diag)
  Call set_lateral_boundary_values(p_patch, &
                                  &p_os%p_prog(nnew(1))%vn)

!CALL message (TRIM(routine), 'end')
END SUBROUTINE calc_normal_velocity_ab_mimetic
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
!! 
SUBROUTINE calc_vert_velocity_mimetic( p_patch, p_diag, top_bc_w, bot_bc_w, pw_c )
!
TYPE(t_patch), TARGET, INTENT(IN) :: p_patch       ! patch on which computation is performed
TYPE(t_hydro_ocean_diag)          :: p_diag
REAL(wp),            INTENT(IN)   :: top_bc_w(:,:) !bottom boundary condition for vertical velocity
REAL(wp),            INTENT(IN)   :: bot_bc_w(:,:) !bottom boundary condition for vertical velocity
REAL(wp),         INTENT(INOUT)   :: pw_c (:,:,:)  ! vertical velocity on cells
!
!
! Local variables
INTEGER :: jc, jk, jb
!INTEGER :: jic, jib
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
INTEGER :: rl_start, rl_end

INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
INTEGER :: rl_start_e, rl_end_e, je
REAL(wp) :: delta_z
REAL(wp) :: z_div_c(nproma,n_zlev,p_patch%nblks_c)
!REAL(wp) :: z_vn(nproma,n_zlev,p_patch%nblks_e)
TYPE(t_cartesian_coordinates):: z_vn_c(nproma,n_zlev,p_patch%nblks_c)
!CHARACTER(len=max_char_length), PARAMETER :: &
!       & routine = ('mo_oce_ab_timestepping_mimetic:calc_vert_velocity_mimetic')
!-----------------------------------------------------------------------  

! #slo# due to nag -nan compiler-option:
z_div_c(:,:,:) = 0.0_wp

rl_start   = 1
rl_end     = min_rlcell
i_startblk = p_patch%cells%start_blk(rl_start,1)
i_endblk   = p_patch%cells%end_blk(rl_end,1)

!------------------------------------------------------------------
! Step 1) Determine vertical velocity in top layer via boundary condition
!------------------------------------------------------------------
pw_c(:,1,:) = top_bc_w(:,:)

!------------------------------------------------------------------
! Step 2) Calculate divergence of horizontal velocity at all levels
!------------------------------------------------------------------
CALL div_oce(  p_diag%ptp_vn, p_patch, z_div_c)


!------------------------------------------------------------------
! Step 3) Use the divergence and the vertical velocity at the previous deeper
!         layer to calculate the new vertical velocity at cell centers
!------------------------------------------------------------------

! !Note we are summing from bottom up to one layer below top.
! !In top layer vertical velocity is given by boundary condition
LEVEL_LOOP2: DO jk = n_zlev, 2, -1

  delta_z = p_patch%patch_oce%del_zlev_m(jk)

  DO jb = i_startblk, i_endblk
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,&
                     & i_startidx, i_endidx, rl_start, rl_end)
    DO jc = i_startidx, i_endidx

      IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

        IF ( jk < p_patch%patch_oce%dolic_c(jc,jb) ) THEN
          ! vertical velocity is integrated from bottom to top and multiplied by
          ! depth of the layer
          ! vertical velocity is negative for positive divergence
          ! of horizontal velocity
          pw_c(jc,jk,jb) = pw_c(jc,jk+1,jb) - z_div_c(jc,jk,jb)*delta_z
          !write(*,*)'vert veloc A:',jc,jb, jk,pw_c(jc,jk,jb), z_div_c(jc,jk,jb)
        ELSEIF ( jk == p_patch%patch_oce%dolic_c(jc,jb) ) THEN
          !use bottom boundary condition for vertical velocity at bottom
          !of prism
          pw_c(jc,jk,jb) = bot_bc_w(jc,jb) - z_div_c(jc,jk,jb)*delta_z

        ELSEIF ( jk > p_patch%patch_oce%dolic_c(jc,jb) ) THEN
          ! Set vertical velocity at center of cell bottom 
          ! to zero if vertical layer is below dolic_c
          pw_c(jc,jk,jb) = 0.0_wp
          !write(*,*)'vert veloc C:',jc,jb, jk,p_patch%patch_oce%dolic_c(jc,jb), pw_c(jc,jk,jb)
        END IF 
      END IF
! IF(pw_c(jc,jk,jb)/=0.0_wp)THEN
! write(92,*)'MIMETIC: vert veloc', jk,jc,jb, pw_c(jc,jk,jb)
! ENDIF
    END DO
  END DO
END DO LEVEL_LOOP2

DO jk = 1,n_zlev
write(*,*)'max/min vert veloc',jk, maxval(pw_c(:,jk,:)), minval(pw_c(:,jk,:))
END DO 

   !jib = i_oct_blk
   !jic = i_oct_idx
   !98 format(3(a,g25.9))
   !  WRITE(message_text,98) &
   !    &     'div(jic,lv1) =', z_div_c   (jic,1,jib), &
   !    &    '     div(lv2) =', z_div_c   (jic,2,jib), &
   !    &    '     div(lv3) =', z_div_c   (jic,3,jib)
   !  CALL message (' ', message_text)
   !  WRITE(message_text,98) &
   !    &     'pwc(jic,lv1) =',    pw_c   (jic,1,jib), &
   !    &    '     pwc(lv2) =',    pw_c   (jic,2,jib), &
   !    &    '     pwc(lv4) =',    pw_c   (jic,4,jib)
   !  CALL message (' ', message_text)

!CALL message (TRIM(routine), 'end')        

END SUBROUTINE calc_vert_velocity_mimetic
!-------------------------------------------------------------------------
FUNCTION inverse_primal_flip_flop(p_patch, rhs_e, h_e) result(inv_flip_flop_e)
   !
   TYPE(t_patch) :: p_patch 
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
   LOGICAL  :: lverbose = .TRUE.
   !CHARACTER(len=MAX_CHAR_LENGTH) :: string
   REAL(wp) :: rhstemp(nproma,p_patch%nblks_e)
   REAL(wp), ALLOCATABLE :: inv_flip_flop_e2(:,:)!(nproma,p_patch%nblks_e)
   INTEGER :: jk
   INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
   INTEGER :: rl_start_e, rl_end_e, je,jb

   !-----------------------------------------------------------------------
   rl_start_e = 1
   rl_end_e  = min_rledge
   i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
   i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

   tolerance                = 1.0E-5!solver_tolerance
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
      CALL gmres_e2e( inv_flip_flop_e2(:,:),&! Input is the first guess
                    & lhs_primal_flip_flop,   &! sbr. calculating l.h.s.
                    & h_e,                    &
                    & p_patch,                &! used for calculating l.h.s.
                    & jk,                     &!idx of vertical level
                    & p_patch%nblks_e,        &
                    & p_patch%npromz_e,       &
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
       WRITE(*,*)'max/min residual of inverse primal-flip-flop:',&
      &jk, maxval(rhstemp),minval(rhstemp) 

        If (maxval (ABS (rhstemp (:,:))) >= tolerance) lmax_iter = .true.
          IF (lverbose) THEN
            IF (lmax_iter) THEN
              WRITE (6, '(1x,a, I4.2, 1x, a,E8.2,1x, a,E8.2,1x, E8.2, 1x, a)') &
              &'Inv_flipflop GMRES #Iter', n_iter, 'Tol ',tolerance, 'Res ',&
              &  ABS(z_residual(n_iter)),MAXVAL (ABS(rhstemp(:,:))), 'GMRES PROBLEM!!!!!!!!!!!!'
            ELSE
              WRITE (6, '(1x,a, I4.2, 1x, a,E8.2,1x, a,E8.2,1x, E8.2)') &
              &'Inv_flipflop GMRES #Iter', n_iter, 'Tol ',tolerance, 'Res ',&
              &  ABS(z_residual(n_iter)),MAXVAL (ABS(rhstemp(:,:)))
            ENDIF
        ENDIF
       inv_flip_flop_e(:,jk,:) = inv_flip_flop_e2(:,:)
      END IF
   END DO

! DO jk=1, 1
!   DO jb = i_startblk_e, i_endblk_e
!     CALL get_indices_e(p_patch, jb,&
!                      & i_startblk_e, i_endblk_e,&
!                      & i_startidx_e, i_endidx_e,&
!                      & rl_start_e, rl_end_e)
!     DO je =  i_startidx_e, i_endidx_e
!       IF(rhs_e(je,jk,jb)/=0.0_wp)THEN
!       write(*,*)'RHS:solution:', jk,je,jb,rhs_e(je,jk,jb), inv_flip_flop_e(je,jk,jb) 
!       ENDIF
!     END DO
!   END DO
! END DO


   DEALLOCATE (inv_flip_flop_e2)

   END FUNCTION inverse_primal_flip_flop
   !--------------------------------------------------------------------
   FUNCTION lhs_primal_flip_flop( x, p_patch, jk, p_coeff, h_e) RESULT(llhs) 
   !
   TYPE(t_patch),INTENT(in)     :: p_patch
   REAL(wp),INTENT(inout)       :: x(:,:)
   INTEGER ,INTENT(in)          :: jk
   REAL(wp),INTENT(in)          :: p_coeff
   REAL(wp),OPTIONAL,INTENT(in) :: h_e(SIZE(x,1), SIZE(x,2))!(:,:)
   REAL(wp)                     :: llhs(SIZE(x,1), SIZE(x,2))

   !locl variables
   INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
   INTEGER :: rl_start_c, rl_end_c
   INTEGER :: jc,jb
   REAL(wp) :: z_x_e(SIZE(x,1), SIZE(x,2))!(nproma,p_patch%nblks_e)
   TYPE(t_cartesian_coordinates) :: z_vn_cc(nproma,p_patch%nblks_c)

  INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  INTEGER :: il_e, ib_e
  INTEGER :: ie,je  
  INTEGER, PARAMETER :: no_cell_edges = 3
  INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  INTEGER :: rl_start_e, rl_end_e 
  REAL(wp) :: z_weight
  REAL(wp) :: z_thick_e
   !-----------------------------------------------------------------------
   rl_start_c = 1
   rl_end_c  = min_rlcell
   rl_start_e = 1
   rl_end_e  = min_rledge

   i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
   i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
   i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
   i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

    DO jb = i_startblk_c, i_endblk_c
     CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
                       i_startidx_c, i_endidx_c, rl_start_c, rl_end_c) 
     DO jc = i_startidx_c, i_endidx_c
       z_vn_cc(jc,jb)%x(1:3) = 0.0_wp
     END DO
   END DO
   z_x_e(:,:) = 0.0_wp
   z_x_e(:,:) = x(:,:)
!write(*,*)'call with level',jk, x(16,1280), SIZE(x,1), SIZE(x,2)
!write(*,*)'max/min input', maxval(z_x_e(:,:)),minval(z_x_e(:,:)), maxval(x(:,:)),minval(x(:,:)) 
!--------------------------------------------------------------------
 CELL_BLK_LOOP: DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb,&
                      & i_startblk_c, i_endblk_c,&
                      & i_startidx_c, i_endidx_c,&
                      & rl_start_c, rl_end_c)

    IF(jk==1)THEN
    !We are dealing with the surface layer first
    CELL_IDX_LOOP_TOP: DO jc =  i_startidx_c, i_endidx_c
      z_weight         = 0.0_wp
      z_vn_cc(jc,jb)%x = 0.0_wp
      z_thick_e        = 0.0_wp
      IF(p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary)THEN
        DO ie=1, no_cell_edges
          il_e = p_patch%cells%edge_idx(jc,jb,ie)
          ib_e = p_patch%cells%edge_blk(jc,jb,ie)

          z_thick_e = p_patch%patch_oce%del_zlev_m(jk) + h_e(il_e,ib_e) 
          z_weight = z_weight + p_patch%patch_oce%variable_vol_norm(jc,jb,ie)!*z_thick_e

           z_vn_cc(jc,jb)%x = z_vn_cc(jc,jb)%x&
                          & + p_patch%patch_oce%edge2cell_coeff_cc(jc,jb,ie)%x&
                          & * z_x_e(il_e,ib_e)! * z_thick_e
        END DO

        z_vn_cc(jc,jb)%x = z_vn_cc(jc,jb)%x / z_weight
      ELSE
       z_vn_cc(jc,jb)%x=0.0_wp 
      ENDIF
    END DO CELL_IDX_LOOP_TOP
    ELSEIF(jk>1)THEN
    !Now we calculate at the levels below the surface
      CELL_IDX_LOOP: DO jc =  i_startidx_c, i_endidx_c
      IF(p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary)THEN
        z_vn_cc(jc,jb)%x = 0.0_wp
        DO ie=1, no_cell_edges
          il_e = p_patch%cells%edge_idx(jc,jb,ie)
          ib_e = p_patch%cells%edge_blk(jc,jb,ie)
          z_vn_cc(jc,jb)%x = z_vn_cc(jc,jb)%x&
                            & + p_patch%patch_oce%edge2cell_coeff_cc(jc,jb,ie)%x&
                            & * z_x_e(il_e,ib_e)
      END DO

        z_vn_cc(jc,jb)%x = z_vn_cc(jc,jb)%x/p_patch%patch_oce%fixed_vol_norm(jc,jb)
      ELSE
        z_vn_cc(jc,jb)%x = 0.0_wp
      ENDIF
    END DO CELL_IDX_LOOP
    ENDIF
END DO CELL_BLK_LOOP

EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e

  CALL get_indices_e(p_patch, jb,&
                   & i_startblk_e, i_endblk_e,&
                   & i_startidx_e, i_endidx_e,&
                   & rl_start_e, rl_end_e)

    EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e
      IF(p_patch%patch_oce%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN

        !Get indices of two adjacent triangles
        il_c1 = p_patch%edges%cell_idx(je,jb,1)
        ib_c1 = p_patch%edges%cell_blk(je,jb,1)
        il_c2 = p_patch%edges%cell_idx(je,jb,2)
        ib_c2 = p_patch%edges%cell_blk(je,jb,2)
!write(*,*)'input', je,jb,z_x_e(je,jb), x(je,jb)
        z_x_e(je,jb) =&
      & DOT_PRODUCT(z_vn_cc(il_c1,ib_c1)%x,p_patch%patch_oce%edge2cell_coeff_cc_t(je,jb,1)%x)&
      &+DOT_PRODUCT(z_vn_cc(il_c2,ib_c2)%x,p_patch%patch_oce%edge2cell_coeff_cc_t(je,jb,2)%x)
       ELSE
         z_x_e(je,jb) = 0.0_wp
       ENDIF
    END DO EDGE_IDX_LOOP
END DO EDGE_BLK_LOOP

   llhs(1:nproma,1:p_patch%nblks_e) = p_coeff*z_x_e(1:nproma,1:p_patch%nblks_e)

!write(*,*)'max/min LHS', maxval(llhs(:,:)),minval(llhs(:,:)) 

END FUNCTION lhs_primal_flip_flop
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
! !       IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! ! 
! !         IF (jk == 1)THEN 
! !           ! at surface level (toplev) add surface elevation h_c
! !           ! #slo# bugfix 11-25
! !           !delta_z = delta_z + ph_c(jc,jb)
! !           delta_z = p_patch%patch_oce%del_zlev_m(jk) + ph_c(jc,jb)
! !         ELSE
! !           ! depth of cell
! !           delta_z = p_patch%patch_oce%del_zlev_m(jk)
! !         ENDIF
! ! 
! !         IF ( jk < p_patch%patch_oce%dolic_c(jc,jb) ) THEN
! ! 
! !           ! vertical velocity is integrated from bottom to top and multiplied by
! !           ! depth of the layer
! !           ! #slo# 2010-10-18: vertical velocity is negative for positive divergence
! !           !                   of horizontal velocity
! !           pw_c(jc,jk,jb) = pw_c(jc,jk+1,jb) - z_div_c(jc,jk,jb) * delta_z
! !           !write(*,*)'vert veloc A:',jc,jb, jk,p_patch%patch_oce%dolic_c(jc,jb), pw_c(jc,jk,jb)
! !         ELSEIF ( jk == p_patch%patch_oce%dolic_c(jc,jb) ) THEN
! !           !use bottom boundary condition for vertical velocity at bottom
! !           !of prism
! !           pw_c(jc,jk,jb) = bot_bc_w(jc,jb) - z_div_c(jc,jk,jb) * delta_z
! ! 
! !         ELSEIF ( jk > p_patch%patch_oce%dolic_c(jc,jb) ) THEN
! !           ! Set vertical velocity at center of cell bottom 
! !           ! to zero if vertical layer is below dolic_c
! !           pw_c(jc,jk,jb) = 0.0_wp
! !           !write(*,*)'vert veloc C:',jc,jb, jk,p_patch%patch_oce%dolic_c(jc,jb), pw_c(jc,jk,jb)
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
