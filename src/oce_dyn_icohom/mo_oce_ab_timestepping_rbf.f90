!>                           
!! Contains the implementation of the semi-implicit Adams-Bashforth timestepping
!! for the ICON ocean model using the RBF based spatial discetization.
!! 
!! 
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/01)
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
MODULE mo_oce_ab_timestepping_rbf
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
!USE mo_math_utilities,            ONLY: t_cartesian_coordinates, gc2cc
USE mo_parallel_config,           ONLY: nproma
USE mo_impl_constants,            ONLY: sea_boundary, &
  &                                     min_rlcell, min_rledge, min_rlcell, &
  &                                     max_char_length, MIN_DOLIC
USE mo_ocean_nml,                 ONLY: n_zlev, solver_tolerance, iforc_oce,&! toplev, &
  &                                     ab_const, ab_beta, ab_gam, iswm_oce,&
                                    &   expl_vertical_velocity_diff, itestcase_oce, l_RIGID_LID
USE mo_run_config,                ONLY: dtime
USE mo_dynamics_config,           ONLY: nold, nnew
USE mo_physical_constants,        ONLY: grav!, re
!USE mo_math_constants,            ONLY: pi, deg2rad
USE mo_oce_state,                 ONLY: t_hydro_ocean_state, & ! t_hydro_ocean_diag
  &                                     v_base, set_lateral_boundary_values
USE mo_model_domain,              ONLY: t_patch
USE mo_ext_data_types,            ONLY: t_external_data
USE mo_oce_linear_solver,         ONLY: gmres_oce
USE mo_exception,                 ONLY: message, finish!, message_text
USE mo_loopindices,               ONLY: get_indices_c, get_indices_e !, get_indices_v
USE mo_oce_boundcond,             ONLY: bot_bound_cond_horz_veloc, top_bound_cond_horz_veloc
USE mo_oce_thermodyn,             ONLY: calc_density, calc_internal_press
USE mo_oce_physics,               ONLY: t_ho_params
!USE mo_oce_forcing,               ONLY: t_sfc_flx!, update_sfcflx
USE mo_sea_ice_types,             ONLY: t_sfc_flx
USE mo_oce_math_operators,        ONLY: div_oce, grad_fd_norm_oce, grad_fd_norm_oce_2d,&
  &                                     height_related_quantities
USE mo_oce_diffusion,             ONLY: velocity_diffusion_horz_RBF, velocity_diffusion_vert_RBF, &
  &                                     veloc_diffusion_vert_impl_hom
USE mo_oce_veloc_advection,       ONLY: veloc_adv_horz_RBF, veloc_adv_vert_RBF
USE mo_intp_data_strc,            ONLY: t_int_state
USE mo_intp_rbf,                  ONLY: rbf_vec_interpol_edge, &
  &                                     rbf_vec_interpol_cell !, verts2edges_scalar
USE mo_math_utilities,            ONLY: gvec2cvec
!USE mo_scalar_product,            ONLY: map_edges2edges
USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
IMPLICIT NONE

PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'

!
! PUBLIC INTERFACE
!
PUBLIC :: solve_free_sfc_ab_RBF
PUBLIC :: calc_normal_velocity_ab_RBF
PUBLIC :: calc_vert_velocity_RBF
! Private implemenation
!
PRIVATE :: fill_rhs4surface_eq_ab_RBF
PRIVATE :: calculate_explicit_term_ab_RBF   ! calc_velocity_predictor
PRIVATE :: lhs_surface_height_ab_RBF

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
SUBROUTINE solve_free_sfc_ab_RBF(p_patch, p_os, p_ext_data, p_sfc_flx, &
  &                              p_phys_param, timestep, p_int)
!
TYPE(t_patch), TARGET, INTENT(in)         :: p_patch
TYPE(t_hydro_ocean_state), TARGET         :: p_os 
TYPE(t_external_data), TARGET, INTENT(in) :: p_ext_data
TYPE(t_sfc_flx), INTENT(INOUT)            :: p_sfc_flx
TYPE (t_ho_params)                        :: p_phys_param
INTEGER                                   :: timestep
TYPE(t_int_state),TARGET,INTENT(IN)       :: p_int
!
!Local variables
!
!REAL(wp),ALLOCATABLE :: rhs(:,:)
!REAL(wp) :: aveh, totarea
!REAL(wp), POINTER :: p_height_area(:,:)
!GMRS
TYPE(t_operator_coeff)                        :: p_op_coeff
INTEGER,PARAMETER :: nmax_iter   = 100      ! maximum number of iterations
REAL(wp):: tolerance                      ! (relative or absolute) tolerance
INTEGER :: n_iter                           ! number of iterations 

REAL(wp):: z_h_c(nproma,p_patch%nblks_c)
REAL(wp):: z_h_e(nproma,p_patch%nblks_e)

REAL(wp):: z_implcoeff
REAL(wp):: zresidual(nmax_iter)    ! norms of the residual (convergence history); 
                                    ! an argument of dimension at least m is required
LOGICAL :: l_maxiter                ! true if reached m iterations
LOGICAL :: lverbose = .TRUE. 
LOGICAL :: l_first_timestep 
CHARACTER(len=max_char_length) :: string
INTEGER :: rl_start_c, rl_end_c, i_startblk_c, i_endblk_c, jc,jb, jk,i_startidx_c, i_endidx_c
!CHARACTER(len=max_char_length), PARAMETER :: &
!       & routine = ('mo_oce_ab_timestepping_rbf:solve_free_sfc_ab_RBF')
!-------------------------------------------------------------------------------
tolerance = solver_tolerance
write(876,*)'-----------------------------',timestep
write(123,*)'---------------------', timestep
z_h_c = 0.0_wp
z_h_e = 0.0_wp
!CALL message (TRIM(routine), 'start')

IF(timestep==1)THEN
  l_first_timestep = .TRUE.

  Call set_lateral_boundary_values(p_patch, p_os%p_prog(nold(1))%vn)
  !Calculation of surface height at edges and total fluid thicknes 
  !of each column
  CALL height_related_quantities(p_patch, p_os, p_ext_data)

ELSE
  l_first_timestep = .FALSE.
  !CALL height_related_quantities(p_patch, p_os, p_ext_data)
ENDIF

write(*,*)'on entry: height:',&
& maxval(p_os%p_prog(nnew(1))%h),minval(p_os%p_prog(nnew(1))%h),&
& maxval(p_os%p_prog(nold(1))%h),minval(p_os%p_prog(nold(1))%h),&
& timestep
write(*,*)'on entry: vn:',&
& maxval(p_os%p_prog(nnew(1))%vn),minval(p_os%p_prog(nnew(1))%vn),&
& maxval(p_os%p_prog(nold(1))%vn),minval(p_os%p_prog(nold(1))%vn)

! #slo# 2011-05-23 - new abort condition for elevation and vn:
IF ( (maxval(p_os%p_prog(nnew(1))%h)  >  1.e20_wp) .or. &
  &  (minval(p_os%p_prog(nnew(1))%h)  < -1.e20_wp) .or. &
  &  (maxval(p_os%p_prog(nold(1))%vn) >  1.e20_wp) .or. &
  &  (minval(p_os%p_prog(nnew(1))%vn) < -1.e20_wp) ) THEN
     CALL message('Solve free surface AB mimetic: ',' INSTABIL VN or H - stop now ')
     CALL finish ('Solve free surface AB mimetic: ',' INSTABIL VN or H !!')
END IF

!calculte velocity reconstructions to be used in the sbr below. Results
!of thus reconstruction processa re stored in diag%vt and in diag%u and diag%v
CALL rbf_vec_interpol_edge( p_os%p_prog(nold(1))%vn,&
                          & p_patch,                &
                          & p_int,                  &
                          & p_os%p_diag%vt,         &
                          & opt_slev=1, opt_elev=n_zlev)
CALL rbf_vec_interpol_cell( p_os%p_prog(nold(1))%vn,&
                          & p_patch,&
                          & p_int,&
                          & p_os%p_diag%u,  &
                          & p_os%p_diag%v, &
                          & opt_slev=1, opt_elev=n_zlev)
rl_start_c = 1
rl_end_c = min_rlcell
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
DO jk=1,n_zlev
write(*,*)'main RBF max/min u, v:',&
& maxval(p_os%p_diag%u),minval(p_os%p_diag%u),&
& maxval(p_os%p_diag%v),minval(p_os%p_diag%v) 

write(876,*)'main RBF max/min u, v:',&
& maxval(p_os%p_diag%u),minval(p_os%p_diag%u),&
& maxval(p_os%p_diag%v),minval(p_os%p_diag%v) 

  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
      &                             rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c

      CALL gvec2cvec(p_os%p_diag%u(jc,jk,jb),        &
                   & p_os%p_diag%v(jc,jk,jb),        &
                   & p_patch%cells%center(jc,jb)%lon,&
                   & p_patch%cells%center(jc,jb)%lat,&
                   & p_os%p_diag%p_vn(jc,jk,jb)%x(1),&
                   & p_os%p_diag%p_vn(jc,jk,jb)%x(2),&
                   & p_os%p_diag%p_vn(jc,jk,jb)%x(3) )
   END DO
  END DO
END DO


IF ( iswm_oce /= 1 ) THEN
  ! Apply windstress 
  CALL top_bound_cond_horz_veloc(p_patch, p_os, p_sfc_flx,  &
                                & p_os%p_aux%bc_top_u, p_os%p_aux%bc_top_v,&
                                & p_os%p_aux%bc_top_veloc_cc) 
  CALL bot_bound_cond_horz_veloc(p_patch, p_os, p_phys_param)
ENDIF

! Calculate explicit terms of Adams-Bashforth timestepping
CALL calculate_explicit_term_ab_RBF(p_patch, p_os, p_phys_param, p_int, l_first_timestep )

IF(.NOT.l_RIGID_LID)THEN
  ! Calculate RHS of surface equation
  CALL fill_rhs4surface_eq_ab_RBF( p_patch, p_os, p_sfc_flx)


  ! Solve surface equation with GMRES solver
  z_h_c = 0.0_wp!p_os%p_prog(nold(1))%h

  !The lhs needs complete thicknesses at edges for 3D and SWE (including bathymetry)
  !IF( iswm_oce == 1 ) THEN  
  z_h_e = p_os%p_diag%thick_e   ! fluid thickness
  !ELSEIF( iswm_oce /= 1 ) THEN 
  !  z_h_e = p_os%p_diag%h_e      ! #slo# 2011-02-21 bugfix (for mimetic only)
  !ENDIF

!   CALL gmres_oce( z_h_c,                   &  ! arg 1 of lhs. x input is the first guess
!         &        lhs_surface_height_ab_RBF,&  ! function calculating l.h.s.
!         &        z_h_e,                   &   !arg 5 of lhs
!         &        p_os%p_diag%thick_c,     &   !arg 6 of lhs, not used, just for compatibility with interface gmres_oce (also used for mimetic)
!         &        p_os%p_prog(nold(1))%h,  &   !arg 2 of lhs
!         &        p_patch,                 &   !arg 3 of lhs
!         &        p_patch%nblks_c,         &
!         &        p_patch%npromz_c,        &
!         &        z_implcoeff,             &
!         &        p_op_coeff,              &
!         &        p_os%p_aux%p_rhs_sfc_eq, &  ! right hand side as input
!         &        tolerance,               &  ! relative tolerance
!         &       .FALSE.,                 &  ! NOT absolute tolerance
!         &        nmax_iter,               &  ! max. # of iterations to do
!         &        l_maxiter,               &  ! out: .true. = not converged
!         &        n_iter,                  &  ! out: # of iterations done
!         &        zresidual )                 ! inout: the residual (array)  

  IF (l_maxiter) THEN
  ! #slo# caution - NAG: variables n_iter and l_maxiter referenced but never set!
  ! rbf-interpolation incorrect (currently not used - 2012/04/25)
    CALL finish('GMRES solver surface equation: ','NOT YET CONVERGED !!')
  ELSE
    IF (lverbose) THEN
       WRITE(string,'(a,i4,a,e20.10)') &
       'iteration ', n_iter,', residual = ', ABS(zresidual(n_iter))
       CALL message('GMRES surface height: ',TRIM(string))
       !write(*,*)'residual:',    REAL(zresidual,wp)!write(*,'(1x,e20.10)') zresidual
    ENDIF
  ENDIF 

  p_os%p_prog(nnew(1))%h = z_h_c

  z_h_c = lhs_surface_height_ab_RBF( p_os%p_prog(nnew(1))%h,&
                                   & p_os%p_prog(nold(1))%h,&
                                   & p_patch,&
                                   & z_implcoeff,&
                                   & p_os%p_diag%thick_e,&
                                   & p_os%p_diag%thick_c,p_op_coeff)&
                                   & -p_os%p_aux%p_rhs_sfc_eq
   write(*,*)'MIN/MAX:h-residual:',&
   & minval(z_h_c(1:nproma,1:p_patch%nblks_c)),&
   & maxval(z_h_c(1:nproma,1:p_patch%nblks_c)) 


 write(*,*)'MIN/MAX:h-new:',&
 & minval(p_os%p_prog(nnew(1))%h(1:nproma,1:p_patch%nblks_c)),&
 & maxval(p_os%p_prog(nnew(1))%h(1:nproma,1:p_patch%nblks_c))
write(876,*)'MIN/MAX:h-new:',&
 & minval(p_os%p_prog(nnew(1))%h(1:nproma,1:p_patch%nblks_c)),&
 & maxval(p_os%p_prog(nnew(1))%h(1:nproma,1:p_patch%nblks_c))
ENDIF
END SUBROUTINE solve_free_sfc_ab_RBF
!-------------------------------------------------------------------------  
!
!  
!>
!! Computation of velocity predictor in Adams-Bashforth timestepping.
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
SUBROUTINE calculate_explicit_term_ab_RBF( p_patch, p_os, p_phys_param, p_int, l_first_timestep)
!
TYPE(t_patch), TARGET, INTENT(in)             :: p_patch
TYPE(t_hydro_ocean_state), TARGET             :: p_os
TYPE (t_ho_params)                            :: p_phys_param
TYPE(t_int_state),TARGET,INTENT(IN), OPTIONAL :: p_int
LOGICAL                                       :: l_first_timestep 
!
!local variables
INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
INTEGER  :: rl_start, rl_end
!INTEGER :: il_c1, ib_c1,il_c2, ib_c2
INTEGER  :: je, jk, jb
REAL(wp) :: gdt
REAL(wp) :: z_gradh_e(nproma,p_patch%nblks_e)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_ab_timestepping_rbf:calculate_explicit_term_ab_RBF')
!-----------------------------------------------------------------------  
gdt=grav*dtime

rl_start = 1
rl_end = min_rledge
i_startblk = p_patch%edges%start_blk(rl_start,1)
i_endblk   = p_patch%edges%end_blk(rl_end,1)


!STEP 1: calculate gradient of surface height at previous timestep
CALL grad_fd_norm_oce_2d( p_os%p_prog(nold(1))%h, &
       &                  p_patch,                &
       &                  z_gradh_e)
! WRITE(*,*)'MAX/MIN old height gradient:', MAXVAL(z_gradh_e),&
! & MINVAL(z_gradh_e) 
! STEP 2: calculate horiz. advcetion
! horiz. advcetion= gradient kinetic energy + nonlinear Coriolis term
!(nonlinear Coriolis term corresponds to (Coriolis param + curl of veloc)*tangential veloc).
CALL veloc_adv_horz_RBF( p_patch,                    &
       &             p_os%p_prog(nold(1))%vn,    &
       &             p_os%p_diag,                &
       &             p_os%p_diag%veloc_adv_horz, & !contains nonlinear Coriolis term, gradient kinetic energy is calculated seperately
       &             p_int                      )

! STEP 3: compute 3D contributions: gradient of hydrostatic pressure and vertical velocity advection
IF ( iswm_oce /= 1 ) THEN
  ! calculate density from EOS using temperature and salinity at timelevel n
  CALL calc_density( p_patch,                                 &
         &           p_os%p_prog(nold(1))%tracer(:,:,:,1:2),  &
         &           p_os%p_diag%rho(:,:,:) )

  ! calculate hydrostatic pressure from density at timelevel nc
   CALL calc_internal_press( p_patch,               &
          &                  p_os%p_diag%rho,       &
          &                  p_os%p_prog(nold(1))%h,&
          &                  p_os%p_diag%press_hyd)

  ! calculate gradient of hydrostatic pressure in 3D
  CALL grad_fd_norm_oce( p_os%p_diag%press_hyd,  &
         &               p_patch,                &
         &               p_os%p_diag%press_grad)


  CALL veloc_adv_vert_RBF( p_patch,                  &
       &             p_os%p_diag%u,              &
       &             p_os%p_diag%v,              &
       &             p_os%p_diag%w,              &
       &             p_os%p_aux%bc_top_u,        &
       &             p_os%p_aux%bc_top_v,        &
       &             p_os%p_aux%bc_bot_u,        &
       &             p_os%p_aux%bc_bot_v,        &
       &             p_os%p_aux%bc_top_w,        &
       &             p_os%p_aux%bc_bot_w,        &
       &             p_os%p_diag%veloc_adv_vert )
  IF(expl_vertical_velocity_diff==0)THEN
  !For alternative "expl_vertical_velocity_diff==1" see couples of
  !lines below below
    CALL velocity_diffusion_vert_RBF( p_patch,            &
    &                             p_os%p_diag%u,          &
    &                             p_os%p_diag%v,          &
    &                             p_os%p_prog(nold(1))%h, &
    &                             p_os%p_aux%bc_top_u,    &
    &                             p_os%p_aux%bc_top_v,    &
    &                             p_os%p_aux%bc_bot_u,    &
    &                             p_os%p_aux%bc_bot_v,    &
    &                             p_phys_param,           &
    &                             p_os%p_diag%laplacian_vert)
  ENDIF
  DO jk=1, n_zlev
    WRITE(*,*)'MAX/MIN density:',       jk, &
      &        MAXVAL(p_os%p_diag%rho(:,jk,:)),&
      & MINVAL(p_os%p_diag%rho(:,jk,:)) 
  END DO
  DO jk=1, n_zlev
    WRITE(*,*)'MAX/MIN internal press:',jk, &
      &        MAXVAL(p_os%p_diag%press_hyd(:,jk,:)),&
      & MINVAL(p_os%p_diag%press_hyd(:,jk,:)) 
  END DO
  DO jk=1, n_zlev
    WRITE(*,*)'MAX/MIN press grad:',    jk, &
      &        MAXVAL(p_os%p_diag%press_grad(:,jk,:)),&
      & MINVAL(p_os%p_diag%press_grad(:,jk,:)) 
  END DO

  DO jk=1, n_zlev
    write(*,*)'MAX/MIN kin energy:', jk,MAXVAL(p_os%p_diag%kin(:,jk,:)),&
   &MINVAL(p_os%p_diag%kin(:,jk,:)) 
  END DO

  DO jk=1, n_zlev
    write(*,*)'MAX/MIN vert advection:', jk,MAXVAL(p_os%p_diag%veloc_adv_vert(:,jk,:)),&
                                       &MINVAL(p_os%p_diag%veloc_adv_vert(:,jk,:)) 
  END DO
 DO jk=1, n_zlev
    write(*,*)'MAX/MIN vert diffusion:', jk,MAXVAL(p_os%p_diag%laplacian_vert(:,jk,:)),&
                                       &MINVAL(p_os%p_diag%laplacian_vert(:,jk,:)) 
  END DO
ELSEIF( iswm_oce == 1 ) THEN
  p_os%p_diag%press_grad     = 0.0_wp
  p_os%p_diag%veloc_adv_vert = 0.0_wp
  p_os%p_diag%laplacian_vert = 0.0_wp
ENDIF

! ! STEP 3: compute laplacian diffusion of velocity
! !  
! ! calculate horizontal laplacian of horizontal velocity, provided this term is discretized explicitly
! ! IF(horizontal_diffusion_veloc==EXPLICIT)THEN
     CALL velocity_diffusion_horz_RBF( p_patch,             &
     &                             p_os%p_prog(nold(1))%vn, &
     &                             p_phys_param,            & 
     &                             p_os%p_diag,             &
     &                             p_os%p_diag%laplacian_horz)
!p_os%p_diag%laplacian_horz=0.0_wp
DO jk=1, n_zlev
    write(*,*)'MAX/MIN horz diffusion:', jk,MAXVAL(p_os%p_diag%laplacian_horz(:,jk,:)),&
                                       &MINVAL(p_os%p_diag%laplacian_horz(:,jk,:)) 
END DO
! ! ELSE
! !  p_os%p_diag%laplacian_horz(:,:,:) = 0.0_wp
! ! ENDIF   

  DO jk = 1, n_zlev

  p_os%p_aux%g_n(:,jk,:) =-p_os%p_diag%press_grad(:,jk,:)       &
     &                   - p_os%p_diag%veloc_adv_horz(:,jk,:)  &
     &                   - p_os%p_diag%veloc_adv_vert(:,jk,:)  &
     &                   + p_os%p_diag%laplacian_horz(:,jk,:)  &
     &                   + p_os%p_diag%laplacian_vert(:,jk,:)

 END DO



IF(l_first_timestep)THEN
  p_os%p_aux%g_nimd(:,:,:) = p_os%p_aux%g_n(:,:,:)
ELSE
  p_os%p_aux%g_nimd(:,:,:) = (1.5_wp+AB_const)* p_os%p_aux%g_n(:,:,:)   &
    &                      - (0.5_wp+AB_const)* p_os%p_aux%g_nm1(:,:,:)
ENDIF
WRITE(*,*)'MAX/MIN old height gradient:', MAXVAL(z_gradh_e),&
& MINVAL(z_gradh_e) 

IF ( iswm_oce /= 1) THEN
  IF(.NOT.l_RIGID_LID)THEN
    DO jb = i_startblk, i_endblk
    CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
    DO jk = 1, n_zlev
      DO je = i_startidx, i_endidx
      IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN

          p_os%p_diag%vn_pred(je,jk,jb) = p_os%p_prog(nold(1))%vn(je,jk,jb)       &
           &                           + dtime*(p_os%p_aux%g_nimd(je,jk,jb)   &
           &                           - (1.0_wp-ab_beta) * grav*z_gradh_e(je,jb))
        ELSE
          p_os%p_diag%vn_pred(je,jk,jb) = 0.0_wp
        ENDIF
      END DO
    END DO
    END DO
  ELSEIF(l_RIGID_LID)THEN

    DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1, n_zlev
        DO je = i_startidx, i_endidx
          IF(v_base%dolic_e(je,jb)>=3)THEN

            p_os%p_diag%vn_pred(je,jk,jb) = p_os%p_prog(nold(1))%vn(je,jk,jb)       &
            &                           + dtime*p_os%p_aux%g_nimd(je,jk,jb)
          ELSE
            p_os%p_diag%vn_pred(je,jk,jb) = 0.0_wp
          ENDIF
        END DO
      END DO
    END DO
  ENDIF

!In the SW-case the external forcing is applied as volume force.
!This force is stored in data type top-boundary-condition. 
ELSEIF ( iswm_oce == 1)THEN! .AND. iforc_oce==11) THEN



  DO jb = i_startblk, i_endblk
    CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
    DO jk = 1, n_zlev
      DO je = i_startidx, i_endidx
        IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN

        p_os%p_diag%vn_pred(je,jk,jb) = p_os%p_prog(nold(1))%vn(je,jk,jb)        &
           &                           + dtime*(p_os%p_aux%g_nimd(je,jk,jb)      &
           &                           - (1.0_wp-ab_beta) * grav*z_gradh_e(je,jb))&
                                      & + p_os%p_aux%bc_top_vn(je,jb)            &
                                      & - p_os%p_aux%bc_bot_vn(je,jb)
        ELSE
          p_os%p_diag%vn_pred(je,jk,jb) = 0.0_wp
        ENDIF
      END DO
    END DO
  END DO
ENDIF
    !IF surface forcing applied as top boundary condition to vertical diffusion
    !The surface forcing is applied as volume forcing at rhs, 
    !i.e. if it part of explicit term in momentum and tracer eqs.
    !in this case, top boundary ondition of vertical Laplacians are homogeneous.
    !Below is the code that adds surface forcing to explicit term of momentum eq. 
    IF(expl_vertical_velocity_diff==1)THEN
      DO jb = i_startblk, i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
          &                rl_start, rl_end)
        DO je = i_startidx, i_endidx
          IF(v_base%dolic_e(je,jb)>=MIN_DOLIC)THEN
            p_os%p_diag%vn_pred(je,1,jb) =  p_os%p_diag%vn_pred(je,1,jb)      &
            &                            + dtime*p_os%p_aux%bc_top_vn(je,jb)/v_base%del_zlev_m(1)

            p_os%p_diag%vn_pred(je,v_base%dolic_e(je,jb),jb)&
            & = p_os%p_diag%vn_pred(je,v_base%dolic_e(je,jb),jb)       &
            & - dtime*p_os%p_aux%bc_bot_vn(je,jb)&
            &/v_base%del_zlev_m(v_base%dolic_e(je,jb))


          ENDIF
        END DO
      END DO
    ENDIF 



IF(expl_vertical_velocity_diff==1.AND.iswm_oce /= 1)THEN
  CALL veloc_diffusion_vert_impl_hom( p_patch,                  &
                                & p_os%p_diag%vn_pred,      &
                                 & p_os%p_diag%h_e,         &
                                & p_phys_param%A_veloc_v,   &
                                & p_os%p_diag%vn_impl_vert_diff)
  IF(l_RIGID_LID)THEN
    p_os%p_diag%vn_pred = p_os%p_diag%vn_impl_vert_diff
  ENDIF
ENDIF



DO jk = 1, n_zlev
write(*,*)'min/max vn_pred:',jk,     minval(p_os%p_diag%vn_pred(:,jk,:)), &
  &                               maxval(p_os%p_diag%vn_pred(:,jk,:))

write(*,*)'min/max vn_pred Term 1:',jk,     minval(p_os%p_aux%g_nimd(:,jk,:)), &
  &                                         maxval(p_os%p_aux%g_nimd(:,jk,:))
write(*,*)'min/max vn_pred Term 2:',jk,     minval((1.0_wp-ab_beta) * grav*z_gradh_e), &
  &                                         maxval((1.0_wp-ab_beta) * grav*z_gradh_e)

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
  &                               maxval(p_os%p_diag%vn_impl_vert_diff(:,jk,:))
ENDIF
write(*,*)
END DO
write(*,*)'min/max -dt*g*grad_h:',minval(dtime*grav*z_gradh_e),&
& maxval(dtime* grav*z_gradh_e)
!-------------------------------------------
  DO jk = 1, n_zlev
    write(876,*)'min/max vn_pred:',jk,     minval(p_os%p_diag%vn_pred(:,jk,:)), &
      &                                  maxval(p_os%p_diag%vn_pred(:,jk,:))

    write(876,*)'min/max vn_pred Term 1:',jk,     minval(p_os%p_aux%g_nimd(:,jk,:)), &
      &                                         maxval(p_os%p_aux%g_nimd(:,jk,:))

    IF(ab_beta/=1.0_wp)THEN
      write(876,*)'min/max vn_pred Term 2:',jk,     minval((1.0_wp-ab_beta) * grav*z_gradh_e), &
        &                                         maxval((1.0_wp-ab_beta) * grav*z_gradh_e)
    ENDIF
    write(876,*)'min/max vn_old :',jk,  minval(p_os%p_prog(nold(1))%vn(:,jk,:)), &
      &                               maxval(p_os%p_prog(nold(1))%vn(:,jk,:))

    write(876,*)'min/max G_n+1/2:',jk,  minval(p_os%p_aux%g_nimd(:,jk,:)),&
                                    & maxval(p_os%p_aux%g_nimd(:,jk,:))

    write(876,*)'min/max G_n:',    jk,  minval(p_os%p_aux%g_n(:,jk,:)),&
                                    & maxval(p_os%p_aux%g_n(:,jk,:))

    write(876,*)'min/max G_n-1:',  jk,  minval(p_os%p_aux%g_nm1(:,jk,:)),&
                                    & maxval(p_os%p_aux%g_nm1(:,jk,:))

      write(876,*)'min/max adv:',    jk,minval(p_os%p_diag%veloc_adv_horz(:,jk,:)), &
      &                               maxval(p_os%p_diag%veloc_adv_horz(:,jk,:))

  END DO
  write(876,*)'min/max -dt*g*grad_h:',minval(dtime*grav*z_gradh_e), maxval(dtime* grav*z_gradh_e)
!-------------------------------------------

!CALL message (TRIM(routine), 'end')        
END SUBROUTINE calculate_explicit_term_ab_RBF
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
SUBROUTINE fill_rhs4surface_eq_ab_RBF( p_patch, p_os, p_sfc_flx)
!
! Patch on which computation is performed
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
!
! Type containing ocean state
TYPE(t_hydro_ocean_state), TARGET :: p_os
TYPE(t_sfc_flx), INTENT(INOUT)    :: p_sfc_flx
!
!  local variables
!
INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
INTEGER :: rl_start_e, rl_end_e
INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
INTEGER :: rl_start_c, rl_end_c
INTEGER :: je, jk, jb, jc
INTEGER :: i_dolic 
REAL(wp) :: gdt2
REAL(wp) :: z_e(nproma,1,p_patch%nblks_e)
REAL(wp) :: div_z_c(nproma,1,p_patch%nblks_c)
REAL(wp) :: z_vn_ab(nproma,n_zlev,p_patch%nblks_e)
!CHARACTER(len=max_char_length), PARAMETER :: &
!       & routine = ('mo_oce_ab_timestepping_rbf:fill_rhs4surface_eq_ab')
!-------------------------------------------------------------------------------
gdt2 = grav*(dtime)**2

rl_start_e = 1
rl_end_e = min_rledge
rl_start_c = 1
rl_end_c = min_rlcell


i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

IF(expl_vertical_velocity_diff==0.OR. iswm_oce == 1)THEN
  z_vn_ab = ab_gam*p_os%p_diag%vn_pred + (1.0_wp -ab_gam)* p_os%p_prog(nold(1))%vn
ELSEIF(expl_vertical_velocity_diff==1)THEN
  z_vn_ab = ab_gam*p_os%p_diag%vn_impl_vert_diff + (1.0_wp -ab_gam)* p_os%p_prog(nold(1))%vn
ENDIF
  DO jk = 1, n_zlev
  write(*,*)'MAX/MIN z_vn_ab:', jk,maxval(z_vn_ab(:,jk,:)),minval(z_vn_ab(:,jk,:))
  write(876,*)'MAX/MIN z_vn_ab:', jk,maxval(z_vn_ab(:,jk,:)),minval(z_vn_ab(:,jk,:))
  END DO
IF( iswm_oce /= 1 ) THEN !the 3D case

  z_e(:,1,:) = 0.0_wp

  DO jb = i_startblk_e, i_endblk_e
    CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
      &                rl_start_e, rl_end_e)
    DO je = i_startidx_e, i_endidx_e
      i_dolic = v_base%dolic_e(je,jb)
      DO jk = 1, i_dolic !n_zlev
        z_e(je,1,jb)  = z_e(je,1,jb) + z_vn_ab(je,jk,jb)*v_base%del_zlev_m(jk)
      END DO
      !Take surface elevation into account
      z_e(je,1,jb)  = z_e(je,1,jb) + z_vn_ab(je,1,jb)*p_os%p_diag%h_e(je,jb)
    ENDDO
  END DO


ELSEIF( iswm_oce == 1 ) THEN !the shallow-water case
  z_e(:,1,:) = z_vn_ab(:,1,:)*p_os%p_diag%thick_e(:,:)
ENDIF

CALL div_oce( z_e, p_patch, div_z_c, 1,1 ) ! to be included surface forcing* +dtime*(P_E)
IF(l_forc_freshw)THEN

  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
      &                rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
        p_os%p_aux%p_rhs_sfc_eq(jc,jb) = (p_os%p_prog(nold(1))%h(jc,jb) &
                                       &- dtime*(div_z_c(jc,1,jb)        &
                                       & + p_sfc_flx%forc_tracer(jc,jb,2) ))/gdt2 !last idx=2 for freshwater
      ELSE
        p_os%p_aux%p_rhs_sfc_eq(jc,jb) = 0.0_wp
      ENDIF
    ENDDO
  END DO

ELSEIF(.NOT.l_forc_freshw)THEN

  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
      &                rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
        p_os%p_aux%p_rhs_sfc_eq(jc,jb) = (p_os%p_prog(nold(1))%h(jc,jb) &
                                       &- dtime*div_z_c(jc,1,jb))/gdt2
      ELSE
        p_os%p_aux%p_rhs_sfc_eq(jc,jb) = 0.0_wp
      ENDIF
    ENDDO
  END DO
ENDIF
   write(*,*)'MAX/MIN thick_e:', maxval(p_os%p_diag%thick_e),minval(p_os%p_diag%thick_e) 
   write(*,*)'MAX/MIN thick_c:', maxval(p_os%p_diag%thick_c),minval(p_os%p_diag%thick_c) 
   write(*,*)'MAX/MIN RHS z_e:', maxval(z_e(:,1,:)),minval(z_e(:,1,:)) 
   write(*,*)'MAX/MIN div_c:', maxval(div_z_c(:,1,:)),minval(div_z_c(:,1,:)) 
   write(*,*)'MAX/MIN RHS:', maxval(p_os%p_aux%p_rhs_sfc_eq(:,:)),&
 &minval(p_os%p_aux%p_rhs_sfc_eq(:,:)) 
!stop
 write(876,*)'MAX/MIN thick_e:', maxval(p_os%p_diag%thick_e),&
 &minval(p_os%p_diag%thick_e) 
   write(876,*)'MAX/MIN thick_c:', maxval(p_os%p_diag%thick_c),&
 &minval(p_os%p_diag%thick_c) 
 write(876,*)'MAX/MIN RHS z_e:', maxval(z_e(:,1,:)),&
 &minval(z_e(:,1,:)) 
 write(876,*)'MAX/MIN div_c:', maxval(div_z_c(:,1,:)),&
 &minval(div_z_c(:,1,:)) 
 write(876,*)'MAX/MIN RHS:', maxval(p_os%p_aux%p_rhs_sfc_eq(:,:)),&
 &minval(p_os%p_aux%p_rhs_sfc_eq(:,:)) 

!CALL message (TRIM(routine), 'end')        
END SUBROUTINE fill_rhs4surface_eq_ab_RBF
!-------------------------------------------------------------------------
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

FUNCTION lhs_surface_height_ab_RBF( p_x, h_old,&
                                 & p_patch,    &
                                 & coeff,      &
                                 & thickness_e,&
                                 & thickness_c,p_op_coeff) RESULT(p_lhs)
!
REAL(wp),    INTENT(IN)   :: p_x(:,:)         ! actual iteration of height, dimension: (nproma,p_patch%nblks_c)
REAL(wp),    INTENT(IN)   :: h_old(:,:)       ! dimension: (nproma,p_patch%nblks_c)
TYPE(t_patch), INTENT(in) :: p_patch          ! patch on which computation is performed
REAL(wp),    INTENT(in)   :: coeff
REAL(wp),    INTENT(in)   :: thickness_e(:,:) !thickness of fluid column at edges
REAL(wp),    INTENT(in)   :: thickness_c(:,:) !not used, just for compatibility with interface gmres_oce (also used for mimetic)
TYPE(t_operator_coeff)    :: p_op_coeff
!
! Left-hand side calculated from iterated height
REAL(wp) :: p_lhs(SIZE(p_x,1), SIZE(p_x,2))  ! (nproma,p_patch%nblks_c)
!
! local variables
REAL(wp) :: gdt2
REAL(wp) :: z_grad_h(nproma,1,p_patch%nblks_e)
REAL(wp) :: z_e(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: div_z_c(SIZE(p_x,1), 1, SIZE(p_x,2))  ! (nproma,1,p_patch%nblks_c)
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
INTEGER :: rl_start, rl_end
INTEGER :: jc, jb
!-----------------------------------------------------------------------  
gdt2 = grav*(dtime)**2
div_z_c = 0.0_wp
rl_start = 1
rl_end   = min_rlcell
i_startblk = p_patch%cells%start_blk(rl_start,1)
i_endblk   = p_patch%cells%end_blk(rl_end,1)

!Step 1) Calculate gradient of iterated height.
CALL grad_fd_norm_oce_2D( p_x, p_patch, z_grad_h(:,1,:))

z_e(:,1,:)=z_grad_h(:,1,:)*thickness_e

!write(*,*)'LHS:max/min :', maxval(z_e(:,1,:)), minval(z_e(:,1,:))
! rl_start = 1
! rl_end = min_rledge
! i_startblk = p_patch%edges%start_blk(rl_start,1)
! i_endblk   = p_patch%edges%end_blk(rl_end,1)
! DO jb = 1,1!i_startblk, i_endblk
!    CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
!  DO jc = i_startidx, i_endidx
!  write(*,*)'z_e', jc,jb,z_e(jc,1,jb)
!    END DO
!  END DO

!Step 3) Calculate divergence
CALL div_oce( z_e, p_patch, div_z_c, top,top )


!Step 4) Finalize LHS calculations
DO jb = i_startblk, i_endblk
  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
  DO jc = i_startidx, i_endidx
    IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
      p_lhs(jc,jb) = (p_x(jc,jb)- gdt2*ab_gam*ab_beta*div_z_c(jc,1,jb))/gdt2
   ELSE
      p_lhs(jc,jb) =0.0_wp
   ENDIF
  END DO
END DO

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
! write(*,*)'max/min LHS h_e', maxval(h_e),minval(h_e) 
!  write(*,*)'max/min LHS', maxval(p_lhs),minval(p_lhs) 

END FUNCTION lhs_surface_height_ab_RBF
!-------------------------------------------------------------------------  
!
!  
!>
!! Computation of new velocity in Adams-Bashforth timestepping.
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
SUBROUTINE calc_normal_velocity_ab_RBF(p_patch, p_os, p_ext_data)
!
! Patch on which computation is performed
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
!
! Type containing ocean state
TYPE(t_hydro_ocean_state), TARGET :: p_os
!
! Type containing external data
TYPE(t_external_data), TARGET, INTENT(in) :: p_ext_data
!
!
!  local variables
!
INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
INTEGER :: rl_start_e, rl_end_e
INTEGER :: je, jk, jb
REAL(wp) :: z_grad_h(nproma,p_patch%nblks_e)
REAL(wp) :: gdt
!CHARACTER(len=max_char_length), PARAMETER ::     &
!  &      routine = ('mo_oce_ab_timestepping_rbf: calc_normal_velocity_ab_RBF')
!-----------------------------------------------------------------------  
!CALL message (TRIM(routine), 'start')        
gdt=grav*dtime
! Step 1) Compute normal derivative of new surface height
 CALL grad_fd_norm_oce_2D(p_os%p_prog(nnew(1))%h, &
   &                      p_patch,                &
   &                      z_grad_h)

! write(*,*)'height:',&
! &p_patch%edges%cell_idx(9,28,1),p_patch%edges%cell_blk(9,28,1),&
! &p_patch%edges%cell_idx(9,28,2),p_patch%edges%cell_blk(9,28,2),&
! & p_os%p_prog(nnew(1))%h(p_patch%edges%cell_idx(9,28,1),p_patch%edges%cell_blk(9,28,1)),&
! & p_os%p_prog(nnew(1))%h(p_patch%edges%cell_idx(9,28,2),p_patch%edges%cell_blk(9,28,2)),&
! &p_os%p_prog(nold(1))%h(p_patch%edges%cell_idx(9,28,1),p_patch%edges%cell_blk(9,28,1)),&
! &p_os%p_prog(nold(1))%h(p_patch%edges%cell_idx(9,28,2),p_patch%edges%cell_blk(9,28,2))

! Step 2) Calculate the new velocity from the predicted one and the new surface height
rl_start_e = 1
rl_end_e = min_rledge
i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

IF(.NOT.l_RIGID_LID)THEN
  DO jb = i_startblk_e, i_endblk_e
    CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
      &                rl_start_e, rl_end_e)
#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = 1, n_zlev
      DO je = i_startidx_e, i_endidx_e

          p_os%p_prog(nnew(1))%vn(je,jk,jb) = p_os%p_diag%vn_pred(je,jk,jb)     &
            &                               - gdt*ab_beta*z_grad_h(je,jb)
      END DO 
    END DO
  END DO

ELSEIF(l_RIGID_LID)THEN

  p_os%p_prog(nnew(1))%vn = p_os%p_diag%vn_pred

ENDIF


write(*,*)
IF(.NOT.l_RIGID_LID)THEN
 write(*,*)'MIN/MAX new height gradient:',minval(z_grad_h),&
                               & maxval(z_grad_h) 
 write(*,*)'MIN/MAX h-contrib to veloc:',minval(- ab_beta*gdt*z_grad_h),&
                                maxval(- ab_beta*gdt*z_grad_h) 
ENDIF
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
&maxval(p_os%p_prog(nnew(1))%vn(:,jk,:) -p_os%p_prog(nold(1))%vn(:,jk,:))!, p_os%p_prog(nnew(1))%vn(1:10,jk,1) 
END DO
IF(.NOT.l_RIGID_LID)THEN
write(876,*)'MIN/MAX new height gradient:', minval(z_grad_h),&
                               & maxval(z_grad_h) 
write(876,*)'MIN/MAX h-contrib to veloc:',minval(- ab_beta*gdt*z_grad_h),&
                                maxval(- ab_beta*gdt*z_grad_h) 
ENDIF
DO jk = 1, n_zlev
  write(876,*)'MIN/MAX vn old:',jk,minval(p_os%p_prog(nold(1))%vn(:,jk,:) ),&
                                 maxval(p_os%p_prog(nold(1))%vn(:,jk,:) ) 
END DO
DO jk = 1, n_zlev
  write(876,*)'MIN/MAX vn new:',jk,minval(p_os%p_prog(nnew(1))%vn(:,jk,:) ),&
                                 maxval(p_os%p_prog(nnew(1))%vn(:,jk,:) )!, p_os%p_prog(nnew(1))%vn(1:10,jk,1) 
END DO

DO jk = 1, n_zlev
  write(876,*)'MIN/MAX vn change:',jk,&
  &minval(p_os%p_prog(nnew(1))%vn(:,jk,:)-p_os%p_prog(nold(1))%vn(:,jk,:) ),&
  &maxval(p_os%p_prog(nnew(1))%vn(:,jk,:)-p_os%p_prog(nold(1))%vn(:,jk,:))
END DO


!stop
! Step 3) Update explicit term for Adams-Bashforth stepping 
 p_os%p_aux%g_nm1 = p_os%p_aux%g_n
 p_os%p_aux%g_n   = 0.0_wp


  CALL height_related_quantities(p_patch, p_os, p_ext_data)
  Call set_lateral_boundary_values(p_patch, &
                                  &p_os%p_prog(nnew(1))%vn)
!CALL message (TRIM(routine), 'end')
!stop
END SUBROUTINE calc_normal_velocity_ab_RBF
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
SUBROUTINE calc_vert_velocity_RBF( p_patch, pvn_new, pvn_old, ph_c, top_bc_w, bot_bc_w, pw_c )
!
TYPE(t_patch), TARGET, INTENT(IN) :: p_patch       ! patch on which computation is performed
! #slo# pvn_e is inout since div_oce sets boundary values to zero
REAL(wp),         INTENT(INOUT) :: pvn_new(:,:,:)  ! new horizontal velocity on edges timestep n+1
REAL(wp),         INTENT(INOUT) :: pvn_old(:,:,:)  ! old horizontal velocity on edges timestep n
REAL(wp),         INTENT(INOUT) :: ph_c(:,:)  ! 
REAL(wp),         INTENT(IN)    :: top_bc_w(:,:)   !top boundary condition for vertical velocity
REAL(wp),         INTENT(IN)    :: bot_bc_w(:,:)   !bottom boundary condition for vertical velocity
REAL(wp),         INTENT(INOUT) :: pw_c (:,:,:)    ! vertical velocity on cells
!
!
! Local variables
!
INTEGER :: jc, jk, jb
!INTEGER :: jic, jib
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
INTEGER :: rl_start, rl_end

!INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
!INTEGER :: rl_start_e, rl_end_e, je
REAL(wp) :: delta_z
REAL(wp) :: z_div_c(nproma,n_zlev,p_patch%nblks_c) ! div_oce needs 3 dimensions
REAL(wp) :: z_vn_ab(nproma,n_zlev,p_patch%nblks_e)
!REAL(wp) :: z_vn_ab2(nproma,n_zlev,p_patch%nblks_e)
!CHARACTER(len=max_char_length), PARAMETER :: &
!       & routine = ('mo_oce_ab_timestepping_rbf:calc_vert_velocity_RBF')
!-----------------------------------------------------------------------  
rl_start = 1
rl_end = min_rlcell
i_startblk = p_patch%cells%start_blk(rl_start,1)
i_endblk   = p_patch%cells%end_blk(rl_end,1)

! #slo# due to nag -nan compiler-option:
z_div_c(:,:,:) = 0.0_wp
z_vn_ab(:,:,:) = 0.0_wp
!z_vn_ab2(:,:,:) = 0.0_wp


!------------------------------------------------------------------
! Step 1) Calculate divergence of horizontal velocity at level jk
!------------------------------------------------------------------
z_vn_ab = ab_gam*pvn_new + (1.0_wp -ab_gam)* pvn_old
! z_vn_ab2 = ab_gam*pvn_new + (1.0_wp -ab_gam)* pvn_old
! CALL map_edges2edges( p_patch, z_vn_ab2, z_vn_ab)
CALL div_oce( z_vn_ab, p_patch, z_div_c)

!Note we are summing from bottom up to one layer below top.
!In top layer vertical velocity is given by boundary condition
LEVEL_LOOP: DO jk = n_zlev, 1, -1

  delta_z = v_base%del_zlev_m(jk)

  !------------------------------------------------------------------
  ! Step 2) Use the divergence and the vertical velocity at the previous deeper
  !         layer to calculate the new vertical velocity at cell centers
  !------------------------------------------------------------------
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,&
                     & i_startidx, i_endidx, rl_start, rl_end)
    DO jc = i_startidx, i_endidx

      IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

!        IF (jk == 1)THEN 
          ! at surface level (toplev) add surface elevation h_c
!          delta_z = v_base%del_zlev_m(jk) + ph_c(jc,jb)
!        ELSE
!          ! depth of cell
!          delta_z = v_base%del_zlev_m(jk)
!        ENDIF

        IF ( jk < v_base%dolic_c(jc,jb) ) THEN

          ! vertical velocity is integrated from bottom to top and multiplied by
          ! depth of the layer
          ! vertical velocity is negative for positive divergence
          ! of horizontal velocity
          pw_c(jc,jk,jb) = pw_c(jc,jk+1,jb) - z_div_c(jc,jk,jb) * delta_z
          !write(*,*)'vert veloc A:',jc,jb, jk,v_base%dolic_c(jc,jb), pw_c(jc,jk,jb)
        ELSEIF ( jk == v_base%dolic_c(jc,jb) ) THEN
          !use bottom boundary condition for vertical velocity at bottom of prism
          pw_c(jc,jk,jb) = bot_bc_w(jc,jb) - z_div_c(jc,jk,jb) * delta_z
          !write(*,*)'vert veloc B:',jc,jb, jk,  pw_c(jc,jk,jb),bot_bc_w(jc,jb), z_div_c(jc,jk,jb)
        ELSEIF ( jk > v_base%dolic_c(jc,jb) ) THEN
          ! Set vertical velocity at center of cell bottom 
          ! to zero if vertical layer is below dolic_c
          pw_c(jc,jk,jb) = 0.0_wp
          !write(*,*)'vert veloc C:',jc,jb, jk,v_base%dolic_c(jc,jb), pw_c(jc,jk,jb)
        END IF 
      END IF
! IF(pw_c(jc,jk,jb)/=0.0_wp)THEN
! write(94,*)'RBF: vert veloc', jk,jc,jb, pw_c(jc,jk,jb)
! ENDIF
    END DO
  END DO
END DO LEVEL_LOOP
write(*,*)
DO jk = 1,n_zlev
write(*,*)'max/min vert veloc',jk, maxval(pw_c(:,jk,:)), minval(pw_c(:,jk,:)),&
&maxval(z_div_c(:,jk,:)), minval(z_div_c(:,jk,:))
write(876,*)'max/min vert veloc',jk, maxval(pw_c(:,jk,:)), minval(pw_c(:,jk,:)),&
&maxval(z_div_c(:,jk,:)), minval(z_div_c(:,jk,:))

END DO

END SUBROUTINE calc_vert_velocity_RBF
!-------------------------------------------------------------------------  


END MODULE mo_oce_ab_timestepping_rbf
