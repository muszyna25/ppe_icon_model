!>
!! Contains the implementation of the tracer transport routines for the ICON ocean model.
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
MODULE mo_oce_tracer_transport
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
USE mo_math_utilities,            ONLY: t_cartesian_coordinates!, gc2cc
USE mo_impl_constants,            ONLY: sea_boundary, &
  &                                     min_rlcell, min_rledge, min_rlcell !, &
!  &                                     max_char_length
USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer, idisc_scheme,    &
                                    &   ab_const, ab_gam, expl_vertical_tracer_diff,&
                                    &   iswm_oce
USE mo_parallel_config,  ONLY: nproma
USE mo_dynamics_config,           ONLY: nold, nnew 
USE mo_run_config,                ONLY: dtime
USE mo_oce_state,                 ONLY: t_hydro_ocean_state!, t_hydro_ocean_diag
USE mo_model_domain,              ONLY: t_patch
!USE mo_exception,                 ONLY:  finish, message_text,message
USE mo_loopindices,               ONLY: get_indices_c, get_indices_e !, get_indices_v
USE mo_oce_boundcond,             ONLY: top_bound_cond_tracer!,&
USE mo_oce_physics
USE mo_oce_forcing,               ONLY: t_ho_sfc_flx!, update_ho_sfcflx
 USE mo_scalar_product,           ONLY:  map_cell2edges

USE mo_oce_math_operators,        ONLY: div_oce!, grad_fd_norm_oce, grad_fd_norm_oce_2d
!USE mo_interpolation,             ONLY: t_int_state
!USE mo_oce_index,                 ONLY: c_i, c_b, c_k, ne_b, ne_i, nc_b, nc_i, form4ar
USE mo_advection_utils,           ONLY: laxfr_upflux, laxfr_upflux_v
USE mo_oce_diffusion,             ONLY: tracer_diffusion_horz, tracer_diffusion_vert_expl,&
                                        & tracer_diffusion_vert_impl
IMPLICIT NONE

PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'

!
! PUBLIC INTERFACE
!
PUBLIC :: advect_tracer_ab
! Private implemenation
!
PRIVATE :: upwind_vflux_oce
PRIVATE :: upwind_hflux_oce
PRIVATE :: advect_individual_tracer_ab

INTEGER, PARAMETER  :: top=1

INTEGER            :: FLUX_CALCULATION 
INTEGER, PARAMETER :: UPWIND = 1
INTEGER, PARAMETER :: CENTRAL= 2
INTEGER, PARAMETER :: MIMETIC= 3

CONTAINS
!-------------------------------------------------------------------------  
!
!  
!>
!! !  SUBROUTINE advects the tracers present in the ocean model.
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
SUBROUTINE advect_tracer_ab(p_patch, p_os, p_param, p_sfc_flx, timestep)
!
!
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
TYPE(t_hydro_ocean_state), TARGET :: p_os
TYPE(t_ho_params), INTENT(inout)  :: p_param
TYPE(t_ho_sfc_flx), INTENT(INOUT) :: p_sfc_flx
INTEGER                           :: timestep! Actual timestep (to distinghuish initial step from others)

!
!Local variables
INTEGER :: jk, jc,jb, i_no_t
REAL(wp) :: max_val, min_val
INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, rl_start_c, rl_end_c
!INTEGER :: ntracers
!CHARACTER(len=max_char_length), PARAMETER :: &
!       & routine = ('mo_tracer_advection:advect_tracer_ab')
!-------------------------------------------------------------------------------
!

IF(idisc_scheme==1)THEN!Mimetic discretization
  FLUX_CALCULATION = MIMETIC
ELSEIF(idisc_scheme==2)THEN!RBF discretization
  FLUX_CALCULATION = CENTRAL
ENDIF

CALL update_ho_params(p_patch, p_os, p_param)

DO i_no_t = 1,no_tracer
  !First tracer is temperature
  !Second tracer is salinity
  IF( iswm_oce /= 1) THEN
    CALL top_bound_cond_tracer( p_patch,  &
                              & p_os,     &
                              & i_no_t,   &
                              & p_sfc_flx,&
                              & p_os%p_aux%bc_top_tracer)
  ENDIF

  CALL advect_individual_tracer_ab( p_patch,                                &
                               & p_os%p_prog(nold(1))%tracer(:,:,:,i_no_t), &
                               & p_os,                                      &
                               & p_os%p_aux%g_n_c_h(:,:,:,i_no_t),          &
                               & p_os%p_aux%g_nm1_c_h(:,:,:,i_no_t),        &
                               & p_os%p_aux%g_nimd_c_h(:,:,:,i_no_t),       &
                               & p_os%p_aux%g_n_c_v(:,:,:,i_no_t),          &
                               & p_os%p_aux%g_nm1_c_v(:,:,:,i_no_t),        &
                               & p_os%p_aux%g_nimd_c_v(:,:,:,i_no_t),       &
                               & p_os%p_aux%bc_top_tracer(:,:,i_no_t),      &
                               & p_os%p_aux%bc_bot_tracer(:,:,i_no_t),      &
                               & p_param%K_tracer_h(:,:,:,i_no_t ),         &
                               & p_param%A_tracer_v(:,:,:, i_no_t),         &
                               & p_os%p_prog(nnew(1))%tracer(:,:,:,i_no_t), &
                               & timestep )

	!  CALL advect_individual_tracer_ab_old( p_patch,                                &
	!                                & p_os%p_prog(nold(1))%tracer(:,:,:,i_no_t), &
	!                                & p_os,                                      &
	!                                & p_os%p_aux%g_n_c_h(:,:,:,i_no_t),          &
	!                                & p_os%p_aux%g_nm1_c_h(:,:,:,i_no_t),        &
	!                                & p_os%p_aux%g_nimd_c_h(:,:,:,i_no_t),       &
	!                                & p_os%p_aux%g_n_c_v(:,:,:,i_no_t),          &
	!                                & p_os%p_aux%g_nm1_c_v(:,:,:,i_no_t),        &
	!                                & p_os%p_aux%g_nimd_c_v(:,:,:,i_no_t),       &
	!                                & p_os%p_aux%bc_top_tracer(:,:,i_no_t),      &
	!                                & p_os%p_aux%bc_bot_tracer(:,:,i_no_t),      &
	!                                & p_param%K_tracer_h(:,:,:,i_no_t ),         &
	!                                & p_param%A_tracer_v(:,:,:, i_no_t),         &
	!                                & p_os%p_prog(nnew(1))%tracer(:,:,:,i_no_t), &
	!                                & timestep )
! 
END DO

rl_start_c   = 1
rl_end_c     = min_rlcell
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

  DO jk = 1, n_zlev
  max_val=0.0_wp
  min_val=20.0_wp
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
        IF(p_os%p_prog(nnew(1))%tracer(jc,jk,jb,1)>max_val)&
          &max_val=p_os%p_prog(nnew(1))%tracer(jc,jk,jb,1)
        IF(p_os%p_prog(nnew(1))%tracer(jc,jk,jb,1)<min_val)&
          &min_val=p_os%p_prog(nnew(1))%tracer(jc,jk,jb,1)

 !       IF(p_os%p_prog(nnew(1))%tracer(jc,jk,jb,1)==0.0_wp)&
 !       & write(*,*)'tracer =0:',jc,jk,jb 

      ENDIF
     END DO
   END DO
   WRITE(*,*)'MAX/MIN trac new:      ', jk, max_val, min_val
    WRITE(*,*)'MAX/MIN tracold1:      ', jk, maxval(p_os%p_prog(nold(1))%tracer(:,jk,:,1)), &
      &                                      minval(p_os%p_prog(nold(1))%tracer(:,jk,:,1))
    WRITE(*,*)'MAX/MIN tracnew1:      ',jk,  maxval(p_os%p_prog(nnew(1))%tracer(:,jk,:,1)), &
      &                                      minval(p_os%p_prog(nnew(1))%tracer(:,jk,:,1))
!     WRITE(*,*)'MAX/MIN tracer 2:      ', jk, maxval(p_os%p_prog(nnew(1))%tracer(:,jk,:,2)), &
!       &                                      minval(p_os%p_prog(nnew(1))%tracer(:,jk,:,2))

END DO
END SUBROUTINE advect_tracer_ab
!-------------------------------------------------------------------------  
!
!
!>
!! !  SUBROUTINE advects the tracers present in the ocean model.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
SUBROUTINE advect_individual_tracer_ab(p_patch, trac_old,           &
                                     & p_os, G_n_c_h,G_nm1_c_h,G_nimd_c_h, &
                                     & G_n_c_v, G_nm1_c_v, G_nimd_c_v, &
                                     & bc_top_tracer, bc_bot_tracer,&
                                     & K_h, A_v,                    &
                                     & trac_new, timestep)
!
!
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp)                          :: trac_old(:,:,:)
TYPE(t_hydro_ocean_state), TARGET :: p_os
REAL(wp) :: G_n_c_h   (nproma, n_zlev,   p_patch%nblks_c)  !G^n
REAL(wp) :: G_nm1_c_h (nproma, n_zlev,   p_patch%nblks_c)  !G^(n-1)
REAL(wp) :: G_nimd_c_h(nproma, n_zlev,   p_patch%nblks_c)  !G^(n+1/2)
REAL(wp) :: G_n_c_v   (nproma, n_zlev,   p_patch%nblks_c)  !G^n
REAL(wp) :: G_nm1_c_v (nproma, n_zlev,   p_patch%nblks_c)  !G^(n-1)
REAL(wp) :: G_nimd_c_v(nproma, n_zlev,   p_patch%nblks_c)  !G^(n+1/2)
REAL(wp) :: bc_top_tracer(nproma, p_patch%nblks_c)
REAL(wp) :: bc_bot_tracer(nproma, p_patch%nblks_c)
REAL(wp) :: K_h(:,:,:)                                  !horizontal mixing coeff
REAL(wp) :: A_v(:,:,:)                                   !vertical mixing coeff
REAL(wp) :: trac_new(:,:,:)                              !new tracer 
INTEGER  :: timestep                                     ! Actual timestep (to distinghuish initial step from others)
!
!Local variables
REAL(wp) :: delta_t
REAL(wp) :: h_tmp(nproma,n_zlev, p_patch%nblks_c)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_tracer_advection:advect_individual_tracer')
!-------------------------------------------------------------------------------
delta_t= dtime
h_tmp = 0.0_wp
CALL advect_horizontal(p_patch, trac_old,           &
                     & p_os, G_n_c_h,G_nm1_c_h,G_nimd_c_h, &
                     & K_h,                    &
                     & trac_new, timestep, delta_t, h_tmp)

IF( iswm_oce /= 1) THEN
   CALL advect_vertical(p_patch, trac_new,              &
                      & p_os,                           &
                      & G_n_c_v, G_nm1_c_v, G_nimd_c_v, &
                      & bc_top_tracer, bc_bot_tracer,   &
                      & A_v,                            &
                      & trac_new, timestep, delta_t, h_tmp)
ENDIF

END SUBROUTINE advect_individual_tracer_ab
  !-----------------------------------------------------------------------
!
!
!>
!! !  SUBROUTINE advects horizontally the tracers present in the ocean model.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!! 
SUBROUTINE advect_horizontal(p_patch, trac_old,           &
                           & p_os, G_n_c_h,G_nm1_c_h,G_nimd_c_h, &
                           & K_h,                    &
                           & trac_new, timestep, delta_t, h_tmp)
!
!
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp) :: trac_old(:,:,:)
TYPE(t_hydro_ocean_state), TARGET :: p_os
!3 arrays for explicit part for tracer in Adams-Bashford  stepping,
!stores information across different timelevels 
REAL(wp) :: G_n_c_h   (nproma, n_zlev,   p_patch%nblks_c)  !G^n
REAL(wp) :: G_nm1_c_h (nproma, n_zlev,   p_patch%nblks_c)  !G^(n-1)
REAL(wp) :: G_nimd_c_h(nproma, n_zlev,   p_patch%nblks_c)  !G^(n+1/2)
REAL(wp) :: K_h(:,:,:)                                   !horizontal mixing coeff
REAL(wp) :: trac_new(:,:,:)                              !new tracer 
INTEGER  :: timestep                                     ! Actual timestep (to distinghuish initial step from others)
REAL(wp) :: delta_t
REAL(wp), OPTIONAL :: h_tmp(nproma,n_zlev, p_patch%nblks_c)
!
!Local variables
REAL(wp) :: delta_z!, delta_z2
!INTEGER  :: ctr, ctr_total
INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, rl_start_c, rl_end_c
INTEGER  :: jc, jk, jb!, jkp1        !< index of edge, vert level, block
!INTEGER  :: z_dolic
REAL(wp) :: z_adv_flux_h(nproma,n_zlev,p_patch%nblks_e)  ! horizontal advective tracer flux
REAL(wp) :: z_div_adv_h(nproma,n_zlev,p_patch%nblks_c)   ! horizontal tracer divergence
REAL(wp) :: z_div_diff_h(nproma,n_zlev,p_patch%nblks_c)  ! horizontal tracer divergence
REAL(wp) :: z_diff_flux_h(nproma,n_zlev,p_patch%nblks_e) ! horizontal diffusive tracer flux
REAL(wp) :: z_transport_vn(nproma,n_zlev,p_patch%nblks_e)! horizontal transport velocity

REAL(wp) :: z_mass_flux_h(nproma,n_zlev, p_patch%nblks_e)
REAL(wp) :: z_trac_c(nproma,n_zlev, p_patch%nblks_c)
REAL(wp) :: z_div_mass_flux_h(nproma,n_zlev, p_patch%nblks_c)
REAL(wp) :: z_h(nproma,n_zlev, p_patch%nblks_c)
REAL(wp) :: z_h2(nproma,n_zlev, p_patch%nblks_c)
REAL(wp) :: z_h_tmp(nproma,n_zlev, p_patch%nblks_c)

!REAL(wp) :: max_val, min_val, dtime2
REAL(wp)           :: z_tol
LOGICAL :: ldbg = .true.
TYPE(t_cartesian_coordinates):: z_vn_c(nproma,n_zlev,p_patch%nblks_c)
TYPE(t_cartesian_coordinates):: z_vn_c2(nproma,n_zlev,p_patch%nblks_c)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_tracer_advection:advect_individual_tracer')
!-------------------------------------------------------------------------------
z_tol= 0.0_wp!1.0E-13

rl_start_c   = 1
rl_end_c     = min_rlcell
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
! rl_start_e   = 1
! rl_end_e     = min_rledge
! i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
! i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

z_adv_flux_h  = 0.0_wp
z_div_adv_h   = 0.0_wp
z_div_diff_h  = 0.0_wp
z_diff_flux_h = 0.0_wp
z_mass_flux_h = 0.0_wp
z_h           = 0.0_wp
z_h2          = 0.0_wp
z_h_tmp       =0.0_wp
z_trac_c      =0.0_wp
z_div_mass_flux_h=0.0_wp
h_tmp  = 0.0_wp

!Generate 3D array of cell heights and of cell height times tracer content
IF ( iswm_oce /= 1) THEN
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jk = 1, n_zlev
      delta_z = p_patch%patch_oce%del_zlev_m(jk)
      DO jc = i_startidx_c, i_endidx_c
        IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          IF (jk == 1) delta_z = p_patch%patch_oce%del_zlev_m(jk)&
                             & + p_os%p_prog(nold(1))%h(jc,jb)
            z_h_tmp(jc,jk,jb) = delta_z 
            z_trac_c(jc,jk,jb)=trac_old(jc,jk,jb)*delta_z
            !write(*,*)'tracer',jc,jb,jk,z_trac_c(jc,jk,jb),trac_old(jc,jk,jb),delta_z
        ENDIF
      END DO
    END DO
  END DO
ELSEIF ( iswm_oce == 1) THEN

  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
      DO jc = i_startidx_c, i_endidx_c
        IF ( p_patch%patch_oce%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
         delta_z           = p_os%p_diag%thick_c(jc,jb)
         z_h_tmp(jc,1,jb) = delta_z 
         z_trac_c(jc,1,jb)= trac_old(jc,1,jb)*delta_z
         !write(123,*)'tracer',jc,jb,z_trac_c(jc,1,jb),trac_old(jc,1,jb),delta_z
       ELSE
         z_h_tmp(jc,1,jb) = 0.0_wp 
         z_trac_c(jc,1,jb)= 0.0_wp
      ENDIF
    END DO
  END DO
ENDIF
!Step 1) Horizontal advection and diffusion
SELECT CASE(FLUX_CALCULATION)

CASE(UPWIND)
  !produce weighted transport velocity
  z_transport_vn = ab_gam*p_os%p_prog(nnew(1))%vn + (1.0_wp-ab_gam)*p_os%p_prog(nold(1))%vn

  !upwind estimate of mass flux
  CALL upwind_hflux_oce( p_patch,        &
                       & z_h_tmp,        &
                       & z_transport_vn, &
                       & z_mass_flux_h )
  !upwind estimate of tracer flux
  CALL upwind_hflux_oce( p_patch,        &
                       & z_trac_c,       &
                       & z_transport_vn, &
                       & z_adv_flux_h )
CASE(CENTRAL)
  !produce weighted transport velocity
  z_transport_vn = ab_gam*p_os%p_prog(nnew(1))%vn + (1.0_wp-ab_gam)*p_os%p_prog(nold(1))%vn
  !central estimate of mass flux
  CALL central_hflux_oce( p_patch,       &
                       & z_h_tmp,        &
                       & z_transport_vn, &
                       & z_mass_flux_h )
  !central estimate of tracer flux
  CALL central_hflux_oce( p_patch,       &
                       & z_trac_c,       &
                       & z_transport_vn, &
                       & z_adv_flux_h )
CASE(MIMETIC)

  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jk = 1, n_zlev
      DO jc = i_startidx_c, i_endidx_c

        z_vn_c(jc,jk,jb)%x(:) =0.0_wp
        z_vn_c2(jc,jk,jb)%x(:)=0.0_wp
        IF (jk == 1)THEN
          delta_z = p_patch%patch_oce%del_zlev_m(jk)+p_os%p_prog(nold(1))%h(jc,jb)

          z_vn_c(jc,jk,jb)%x  = trac_old(jc,jk,jb)*delta_z*p_os%p_diag%p_vn(jc,jk,jb)%x
          z_vn_c2(jc,jk,jb)%x = delta_z*p_os%p_diag%p_vn(jc,jk,jb)%x
         ELSE
          z_vn_c(jc,jk,jb)%x = trac_old(jc,jk,jb)*p_os%p_diag%p_vn(jc,jk,jb)%x
          z_vn_c2(jc,jk,jb)%x= p_os%p_diag%p_vn(jc,jk,jb)%x
         ENDIF
      END DO
    END DO
  END DO

  CALL map_cell2edges( p_patch, z_vn_c, z_adv_flux_h )
  CALL map_cell2edges( p_patch, z_vn_c2, z_mass_flux_h )

END SELECT
 
!Step 4) calculate horizontal divergence of mass and tracer flux
CALL div_oce( z_adv_flux_h, p_patch, z_div_adv_h)
CALL div_oce( z_mass_flux_h, p_patch, z_div_mass_flux_h)

 DO jk = 1, n_zlev
 write(*,*)'max/min adv-flux:',jk,&
 &maxval(z_adv_flux_h(:,jk,:)),&
 &minval(z_adv_flux_h(:,jk,:))
 END DO
 DO jk = 1, n_zlev
 write(*,*)'max/min mass-flux:',jk,&
 &maxval(z_mass_flux_h(:,jk,:)),&
 &minval(z_mass_flux_h(:,jk,:))
 END DO
!calculate (dummy) height consistent with divergence of mass flux
DO jk=1,n_zlev
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
      z_h(jc,jk,jb) = z_h_tmp(jc,jk,jb) - delta_t*z_div_mass_flux_h(jc,jk,jb)
      ELSE
      z_h(jc,jk,jb) = 0.0_wp
      ENDIF
    END DO
  END DO!write(*,*)'max/min h_dummy:',maxval(z_h(:,jk,:)),minval(z_h(:,jk,:))
END DO 
 DO jk = 1, n_zlev
  write(*,*)'max/min dummy height:',jk,&
  &maxval(z_h(:,jk,:)),&
  &minval(z_h(:,jk,:))
  END DO
IF(present(h_tmp))THEN
  h_tmp = z_h 
ENDIF
!The diffusion part: calculate horizontal diffusive flux
CALL tracer_diffusion_horz( p_patch, &
 &                          trac_old,&
 &                          p_os,    &
 &                          K_h,     & 
 &                          z_diff_flux_h)

!Calculate divergence of diffusive flux
CALL div_oce( z_diff_flux_h, p_patch, z_div_diff_h)

!Calculate explicit part, stored in G_n, later
!AB-extrapolation us computed
DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                   & i_startidx_c, i_endidx_c,&
                   & rl_start_c, rl_end_c)
  DO jk = 1, n_zlev
    DO jc = i_startidx_c, i_endidx_c
      IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
        G_n_c_h(jc,jk,jb) = z_div_diff_h(jc,jk,jb)
        !write(*,*)'div',jc,jk,jb,G_n_c_h(jc,jk,jb),z_div_diff_h(jc,jk,jb)
      ENDIF
    END DO
  END DO
END DO

IF(timestep/=1)THEN
    G_nimd_c_h(:,:,:) = (1.5_wp+AB_const)* G_n_c_h(:,:,:)   &
        &             - (0.5_wp+AB_const)*G_nm1_c_h(:,:,:)
  ELSE
    G_nimd_c_h(:,:,:) = G_n_c_h(:,:,:)
ENDIF

IF(ldbg)THEN
  DO jk = 1, n_zlev
    write(*,*)'max/min G_T:',jk,&
    &maxval(G_n_c_h(:,jk,:)), minval(G_n_c_h(:,jk,:))
  END DO
  DO jk = 1, n_zlev
    write(*,*)'max/min G_T+1/2:',jk,maxval(G_nimd_c_h(:,jk,:)),&
                                   &minval(G_nimd_c_h(:,jk,:))
  END DO
  DO jk = 1, n_zlev
    write(*,*)'max/min adv flux & div h:',jk,maxval(z_adv_flux_h(:,jk,:)),&
                                           & minval(z_adv_flux_h(:,jk,:)),&
                                           & maxval(z_div_adv_h(:,jk,:)),&
                                           & minval(z_div_adv_h(:,jk,:))
  END DO
  DO jk = 1, n_zlev
    write(*,*)'max/min mass flux & div h:',jk, maxval(z_mass_flux_h(:,jk,:)),&
                                             & minval(z_mass_flux_h(:,jk,:)),&
                                             & maxval(z_div_mass_flux_h(:,jk,:)),&
                                             & minval(z_div_mass_flux_h(:,jk,:))
  END DO
ENDIF

!Final step: calculate new tracer values
!write(123,*)'timestep----------------------',timestep
DO jk = 1, n_zlev
  CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                   & i_startidx_c, i_endidx_c,&
                   & rl_start_c, rl_end_c)
    DO jb = i_startblk_c, i_endblk_c
      DO jc = i_startidx_c, i_endidx_c
        IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

!           IF(jk==1)THEN
!             delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+p_patch%patch_oce%del_zlev_m(top)
!             !delta_z2  = p_os%p_prog(nnew(1))%h(jc,jb)+p_patch%patch_oce%del_zlev_m(top)
!           ELSE
!             delta_z  = p_patch%patch_oce%del_zlev_m(jk)
!            !delta_z2  = p_patch%patch_oce%del_zlev_m(jk) 
!           ENDIF

          IF(z_h(jc,jk,jb)/=0.0_wp)THEN
            trac_new(jc,jk,jb) = (z_trac_c(jc,jk,jb)             &
                               & -delta_t*(z_div_adv_h(jc,jk,jb)  &
                               &           -G_nimd_c_h(jc,jk,jb)))&
                               &/z_h(jc,jk,jb)
! write(123,*)'trac calc:',jc,jb,jk,&
! &trac_new(jc,jk,jb),&
! &trac_old(jc,jk,jb),&
! &z_trac_c(jc,jk,jb),&
! &z_div_adv_h(jc,jk,jb),&
! &G_nimd_c_h(jc,jk,jb),&
! &z_h(jc,jk,jb),&
! &z_trac_c(jc,jk,jb)/z_h(jc,jk,jb),&
! &delta_t*(z_div_adv_h(jc,jk,jb)-G_nimd_c_h(jc,jk,jb))/z_h(jc,jk,jb)
          ELSE
            trac_new(jc,jk,jb) = 0.0_wp
          ENDIF

        ENDIF
      END DO
    END DO
END DO

!swap the G_n-values: n -> (n-1)
G_nm1_c_h(:,:,:)  = G_n_c_h(:,:,:)
G_nimd_c_h(:,:,:) = 0.0_wp
G_n_c_h(:,:,:)    = 0.0_wp

!coarse check
!  DO jk = 1, n_zlev
  !max_val=0.0_wp
  !min_val=100.0_wp
!  DO jb = i_startblk_c, i_endblk_c
!    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
!    &                   rl_start_c, rl_end_c)
!    DO jc = i_startidx_c, i_endidx_c
!       IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
!         IF(trac_new(jc,jk,jb)>max_val)&
!           &max_val=trac_new(jc,jk,jb)
!         IF(trac_new(jc,jk,jb)<min_val)&
!           &min_val=trac_new(jc,jk,jb)
! 
!         IF(trac_new(jc,jk,jb)==0.0_wp)&
!         & write(*,*)'tracer =0:',jc,jk,jb 
!       ENDIF
!     END DO
!   END DO
!   !WRITE(*,*)'MAX/MIN tracer after horizontal transport:', jk, max_val, min_val
!END DO
END SUBROUTINE advect_horizontal
!-------------------------------------------------------------------------  
!
!
!>
!! !  SUBROUTINE advects vertically the tracers present in the ocean model.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
SUBROUTINE advect_vertical(p_patch, trac_in,               &
                         & p_os,                           &
                         & G_n_c_v, G_nm1_c_v, G_nimd_c_v, &
                         & bc_top_tracer, bc_bot_tracer,   &
                         & A_v,                            &
                         & trac_out, timestep, delta_t, h_tmp)
!
!
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp)                          :: trac_in(:,:,:)
TYPE(t_hydro_ocean_state), TARGET :: p_os
REAL(wp) :: G_n_c_v   (nproma, n_zlev,   p_patch%nblks_c)  !G^n
REAL(wp) :: G_nm1_c_v (nproma, n_zlev,   p_patch%nblks_c)  !G^(n-1)
REAL(wp) :: G_nimd_c_v(nproma, n_zlev,   p_patch%nblks_c)  !G^(n+1/2)
REAL(wp) :: bc_top_tracer(nproma, p_patch%nblks_c)
REAL(wp) :: bc_bot_tracer(nproma, p_patch%nblks_c)
REAL(wp) :: A_v(:,:,:)                                   !vertical mixing coeff
REAL(wp) :: trac_out(:,:,:)                              !new tracer 
INTEGER  :: timestep                                     ! Actual timestep (to distinghuish initial step from others)
REAL(wp) :: delta_t
REAL(wp), OPTIONAL :: h_tmp(nproma,n_zlev, p_patch%nblks_c)
!
!Local variables
REAL(wp) :: delta_z!, delta_z2
!INTEGER  :: ctr, ctr_total
INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, rl_start_c, rl_end_c
INTEGER  :: jc, jk, jb, jkp1        !< index of edge, vert level, block
INTEGER  :: z_dolic
REAL(wp) :: z_adv_flux_v (nproma, n_zlev+1, p_patch%nblks_c)  ! vertical advective tracer flux
REAL(wp) :: z_div_adv_v (nproma, n_zlev,p_patch%nblks_c)        ! vertical tracer divergence
REAL(wp) :: z_div_diff_v (nproma, n_zlev,p_patch%nblks_c)        ! vertical tracer divergence
REAL(wp) :: z_diff_flux_v(nproma, n_zlev+1,p_patch%nblks_c)   ! vertical diffusive tracer flux
REAL(wp) :: z_transport_w(nproma,n_zlev+1,p_patch%nblks_c)  ! vertical transport velocity
REAL(wp) :: z_trac_c(nproma,n_zlev, p_patch%nblks_c)
REAL(wp) :: z_h(nproma,n_zlev, p_patch%nblks_c)
REAL(wp) :: z_G_n_c_v   (nproma, n_zlev, p_patch%nblks_c)
REAL(wp) :: z_G_nm1_c_v (nproma, n_zlev, p_patch%nblks_c)
REAL(wp) :: z_G_nimd_c_v(nproma, n_zlev, p_patch%nblks_c)
!REAL(wp) :: max_val, min_val!, dtime2
REAL(wp) :: z_tol
!LOGICAL  :: ldbg = .false.
!TYPE(t_cartesian_coordinates):: z_vn_c(nproma,n_zlev,p_patch%nblks_c)
!TYPE(t_cartesian_coordinates):: z_vn_c2(nproma,n_zlev,p_patch%nblks_c)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_tracer_advection:advect_individual_tracer')
!-------------------------------------------------------------------------------
z_tol= 0.0_wp!1.0E-13

rl_start_c   = 1
rl_end_c     = min_rlcell
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

! rl_start_e   = 1
! rl_end_e     = min_rledge
! i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
! i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

z_adv_flux_v  = 0.0_wp
z_div_adv_v   = 0.0_wp
z_div_diff_v  = 0.0_wp
z_diff_flux_v = 0.0_wp
z_h           = 0.0_wp
z_transport_w =0.0_wp
z_trac_c      =0.0_wp
z_G_n_c_v        = 0.0_wp
z_G_nm1_c_v      = 0.0_wp
z_G_nimd_c_v     = 0.0_wp

!tracer times heigth; this includes free surface
!height in top cell
DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
  &                   rl_start_c, rl_end_c)
  DO jk = 1, n_zlev
    delta_z = p_patch%patch_oce%del_zlev_m(jk)
    DO jc = i_startidx_c, i_endidx_c
      IF (jk == 1) THEN
        delta_z = p_patch%patch_oce%del_zlev_m(jk)&
              & + p_os%p_prog(nold(1))%h(jc,jb)
      ENDIF
      z_trac_c(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z
    END DO
  END DO
END DO

!Produce time-weighted vertical transport velocity
z_transport_w  = ab_gam*p_os%p_diag%w + (1.0_wp-ab_gam)*p_os%p_diag%w_old

SELECT CASE(FLUX_CALCULATION)

CASE(UPWIND)
  CALL upwind_vflux_oce( p_patch,      &
                       & trac_in,      &
                       & z_transport_w,&
                       & z_adv_flux_v )
CASE(CENTRAL)
  CALL central_vflux_oce( p_patch,     &
                       & trac_in,      &
                       & z_transport_w,&
                       & z_adv_flux_v )
CASE(MIMETIC)
  CALL mimetic_vflux_oce( p_patch,      &
                       & trac_in,       &
                       & z_transport_w, &
                       & z_adv_flux_v )
END SELECT
 
!calculate mass flux
DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
  &                   rl_start_c, rl_end_c)
  DO jc = i_startidx_c, i_endidx_c
    z_dolic = p_patch%patch_oce%dolic_c(jc,jb)
    IF(z_dolic>2)THEN
       z_h(jc,1,jb) = (p_patch%patch_oce%del_zlev_m(1) + p_os%p_prog(nold(1))%h(jc,jb))&
                   &- (z_transport_w(jc,1,jb)-z_transport_w(jc,2,jb))
       DO jk = 2, z_dolic-1
         z_h(jc,jk,jb) = p_patch%patch_oce%del_zlev_m(jk)&
                     & - (z_transport_w(jc,jk,jb)-z_transport_w(jc,jk+1,jb))
       END DO
       z_h(jc,z_dolic,jb) = p_patch%patch_oce%del_zlev_m(z_dolic)
   ENDIF
  END DO
END DO

!The diffusion part
!vertical diffusion explicit = 0 or implicit =1
IF(expl_vertical_tracer_diff==1)THEN

  !divergence is calculated for advective fluxes
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      !interior: from one below surface to the ground
      z_dolic = p_patch%patch_oce%dolic_c(jc,jb)
      IF(z_dolic>0)THEN
        DO jk = 1, z_dolic
          jkp1 = jk + 1
          ! positive vertical divergence in direction of w (upward positive)
          z_div_adv_v(jc,jk,jb) = &
          & (z_adv_flux_v(jc,jk,jb)-z_adv_flux_v(jc,jkp1,jb))

        END DO!level-loop
      ENDIF!(z_dolic>0)
    END DO!idx-loop
  END DO !blk-loop

  !Add advective part
  DO jk = 1, n_zlev 
    CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                     & i_startidx_c, i_endidx_c,&
                     & rl_start_c, rl_end_c)
    DO jb = i_startblk_c, i_endblk_c
      DO jc = i_startidx_c, i_endidx_c
        IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          IF(jk==1)THEN
            delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+p_patch%patch_oce%del_zlev_m(top)
          ELSE
            delta_z  = p_patch%patch_oce%del_zlev_m(jk)
          ENDIF
            z_trac_c(jc,jk,jb) = (trac_in(jc,jk,jb)*delta_z     &
                               & -delta_t*z_div_adv_v(jc,jk,jb))  &
                               &/(z_h(jc,jk,jb)+z_tol)
        ENDIF
      END DO
    END DO
  END DO
  !calculate vert diffusion impicit: result is stored in trac_out
  CALL tracer_diffusion_vert_impl( p_patch,     &
                              & z_trac_c(:,:,:),&
                              & bc_top_tracer,  & 
                              & bc_bot_tracer,  &
                              & A_v,            &
                              & trac_out(:,:,:))

!vertival diffusion is calculated explicitely
ELSEIF(expl_vertical_tracer_diff==0)THEN

  IF(present(h_tmp))THEN

  !calculation of diffusive fluxes
    CALL tracer_diffusion_vert_expl( p_patch,       &
                                  &  trac_in,       &
                                  & h_tmp,          &
                                  &  bc_top_tracer, &
                                  &  bc_bot_tracer, &
                                  &  A_v,           &
                                  &  z_diff_flux_v)
  ELSE
    !calculation of vertical diffusive fluxes
    CALL tracer_diffusion_vert_expl( p_patch,       &
                                  &  trac_in,       &
                                  &  z_h,           &
                                  &  bc_top_tracer, &
                                  &  bc_bot_tracer, & 
                                  &  A_v,           &
                                  &  z_diff_flux_v)

  ENDIF
  DO jk = 1, n_zlev+1
    write(*,*)'max/min adv & diff tracer flux v:',jk,&
    &maxval(z_adv_flux_v(:,jk,:)), minval(z_adv_flux_v(:,jk,:)),&
    &maxval(z_diff_flux_v(:,jk,:)),minval(z_diff_flux_v(:,jk,:))
  END DO
  DO jk = 1, n_zlev
    write(*,*)'max/min vertical mass flux:',jk,&
    &maxval(z_h(:,jk,:)),minval(z_h(:,jk,:))
  END DO

  ! divergence is calculated for advective & diffusive fluxes
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      !interior: from one below surface to the ground
      z_dolic = p_patch%patch_oce%dolic_c(jc,jb)
      IF(z_dolic>0)THEN
        DO jk = 1, z_dolic
          jkp1 = jk + 1 

          IF(jk==1)THEN
            delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+p_patch%patch_oce%del_zlev_m(top)
          ELSE
            delta_z  = p_patch%patch_oce%del_zlev_m(jk)
          ENDIF
          ! positive vertical divergence in direction of w (upward positive)
          z_div_adv_v(jc,jk,jb) = z_adv_flux_v(jc,jk,jb)-z_adv_flux_v(jc,jkp1,jb)

          G_n_c_v(jc,jk,jb) = z_diff_flux_v(jc,jk,jb)-z_diff_flux_v(jc,jkp1,jb)

        END DO!level-loop
      ENDIF!(z_dolic>0)
    END DO!idx-loop
  END DO !blk-loop

  IF(timestep/=1)THEN
    G_nimd_c_v(:,:,:) = (1.5_wp+AB_const)* G_n_c_v(:,:,:)   &
      &               - (0.5_wp+AB_const)*G_nm1_c_v(:,:,:)
  ELSE
    G_nimd_c_v(:,:,:) = G_n_c_v(:,:,:)
  ENDIF

  DO jk = 1, n_zlev
    write(*,*)'max/min div adv & diffusive  flux v:',&
    &jk,maxval(z_div_adv_v(:,jk,:)), minval(z_div_adv_v(:,jk,:)),&
    &maxval(G_n_c_v(:,jk,:)),minval(G_n_c_v(:,jk,:))
  END DO

  !Final step: calculate new tracer values
  DO jk = 1, n_zlev 
    CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                     & i_startidx_c, i_endidx_c,&
                     & rl_start_c, rl_end_c)
    DO jb = i_startblk_c, i_endblk_c
      DO jc = i_startidx_c, i_endidx_c
        IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

          IF(jk==1)THEN
            delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+p_patch%patch_oce%del_zlev_m(jk)
          ELSE
            delta_z  = p_patch%patch_oce%del_zlev_m(jk)
          ENDIF
            trac_out(jc,jk,jb) = (trac_in(jc,jk,jb)*delta_z     &
                               & -delta_t*(z_div_adv_v(jc,jk,jb)  &
                               & -G_nimd_c_v(jc,jk,jb)))         &
                               &/(z_h(jc,jk,jb)+z_tol)
        ENDIF
      END DO
    END DO
  END DO

  G_nm1_c_v(:,:,:)  = G_n_c_v(:,:,:)
  G_nimd_c_v(:,:,:) = 0.0_wp
  G_n_c_v(:,:,:)    = 0.0_wp

ENDIF!(lvertical_diff_implicit)THEN

! rl_start_c   = 1
! rl_end_c     = min_rlcell
! i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
! i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
! 
!   DO jk = 1, n_zlev
!   max_val=0.0_wp
!   min_val=100.0_wp
!   DO jb = i_startblk_c, i_endblk_c
!     CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
!     &                   rl_start_c, rl_end_c)
!     DO jc = i_startidx_c, i_endidx_c
!       IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
!         IF(trac_new(jc,jk,jb)>max_val)&
!           &max_val=trac_new(jc,jk,jb)
!         IF(trac_new(jc,jk,jb)<min_val)&
!           &min_val=trac_new(jc,jk,jb)
! 
!         IF(trac_new(jc,jk,jb)==0.0_wp)&
!         & write(*,*)'tracer =0:',jc,jk,jb 
!       ENDIF
!      END DO
!    END DO
!    WRITE(*,*)'MAX/MIN tracer after vertical transport:', jk, max_val, min_val
! END DO
!    DO jk = 1, n_zlev
!     ctr=0
!     ctr_total=0
!     DO jb = i_startblk_c, i_endblk_c
!       CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
!                        & i_startidx_c, i_endidx_c,&
!                        & rl_start_c, rl_end_c)
!    DO jc = i_startidx_c, i_endidx_c
!      IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary )THEN
!      ctr_total=ctr_total+1
!      ENDIF
!  
!     IF(trac_new(jc,jk,jb)>t_prof(jk)) THEN
!     ctr=ctr+1
! !  write(*,*)'indices',jc,jb,jk, p_patch%patch_oce%lsm_oce_c(jc,top,jb),&
! ! ! !&trac_old(jc,jk,jb),
! ! ! &(delta_z/delta_z2)*trac_old(jc,top,jb),&
! !  & dtime*G_nimd_c(jc,top,jb),G_n_c_h(jc,top,jb),z_div_adv_h(jc,top,jb),&
! !  trac_new(jc,top,jb), ctr
!   ENDIF
!        END DO
!    END DO
!    write(*,*)'counter > initial:',jk, ctr, ctr_total
!  END DO


END SUBROUTINE advect_vertical
!-----------------------------------------------------------------------
  !>
  !! First order upwind scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !! Developed by L.Bonaventura  (2004).
  !! Adapted to new grid structure by L. Bonaventura, MPI-M, August 2005.
  !! Modification by Daniel Reinert, DWD (2010-02-09)
  !! - transferred to separate subroutine
  !! Modification by Stephan Lorenz, MPI (2010-09-06)
  !! - adapted to hydrostatic ocean core
  !!
  SUBROUTINE upwind_hflux_oce( ppatch, pvar_c, pvn_e, pupflux_e, opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch      !< patch on which computation is performed
    REAL(wp), INTENT(INOUT)  :: pvar_c(:,:,:)      !< advected cell centered variable
    REAL(wp), INTENT(INOUT)  :: pvn_e(:,:,:)       !< normal velocity on edges
    !EAL(wp), INTENT(INOUT), OPTIONAL :: ph_e (:,:)         !< surface elevation on edges
    REAL(wp), INTENT(INOUT)  :: pupflux_e(:,:,:)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_slev    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_elev    ! optional vertical end level
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER  :: slev, elev
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF

    rl_start = 1
    rl_end   = min_rledge

    i_startblk = ppatch%edges%start_blk(rl_start,1)
    i_endblk   = ppatch%edges%end_blk(rl_end,1)
    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    !
    ! for no-slip boundary conditions, boundary treatment for tracer (zero at leteral walls) 
    !is implicit done via velocity boundary conditions
    !
    ! line and block indices of two neighboring cells
    iilc => ppatch%edges%cell_idx
    iibc => ppatch%edges%cell_blk

    ! loop through all patch edges (and blocks)
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(ppatch, jb, i_startblk, i_endblk,   &
        &                i_startidx, i_endidx, rl_start,rl_end)

#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator 
          IF ( ppatch%patch_oce%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
            pupflux_e(je,jk,jb) =  &
            &  laxfr_upflux( pvn_e(je,jk,jb), pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), &
            &                                 pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )
!write(*,*)'upwind flux -h',je,jk,jb,pupflux_e(je,jk,jb), pvn_e(je,jk,jb),&
!&pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2))
          ELSE
            pupflux_e(je,jk,jb) = 0.0_wp
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
!    IF (p_pe == p_io .AND. ldbg) THEN
!      WRITE(*,form4ar) &
!       &  '  zhorflx  E1 =', zhorflx_e(ne_i(1),c_k,ne_b(1)), &
!       & '         ph_e  =', ph_e(ne_i(1),ne_b(1))
!      WRITE(*,form4ar) &
!       &  '  zhorflx  E2 =', zhorflx_e(ne_i(2),c_k,ne_b(2)), &
!       & '         ph_e  =', ph_e(ne_i(2),ne_b(2))
!      WRITE(*,form4ar) &
!       &  '  zhorflx  E3 =', zhorflx_e(ne_i(3),c_k,ne_b(3)), &
!       & '         ph_e  =', ph_e(ne_i(3),ne_b(3))
!       WRITE(*,form4ar) &
!         &  '  pupflx_h E1 =', pupflux_e(ne_i(1),c_k,ne_b(1)), &
!         &  '    pvar_c(1) =', pvar_c   (nc_i(1),c_k,nc_b(1)), &
!         &  '    pvar_c(2) =', pvar_c   (nc_i(2),c_k,nc_b(2))!     END IF
  END SUBROUTINE upwind_hflux_oce
  !-----------------------------------------------------------------------
  !>
  !! Central scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using central fluxes
  !!
  !! @par Revision History
  !! Peter korn, MPI-M, 2011
  !!
  SUBROUTINE central_hflux_oce( ppatch, pvar_c, pvn_e, pupflux_e )

    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch      !< patch on which computation is performed
    REAL(wp), INTENT(INOUT)  :: pvar_c(:,:,:)      !< advected cell centered variable
    REAL(wp), INTENT(INOUT)  :: pvn_e(:,:,:)       !< normal velocity on edges
    !REAL(wp), INTENT(INOUT)  :: ph_e (:,:)         !< surface elevation on edges
    REAL(wp), INTENT(INOUT)  :: pupflux_e(:,:,:)   !< variable in which the upwind flux is stored

    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block
    !-----------------------------------------------------------------------
    rl_start   = 1
    rl_end     = min_rledge
    i_startblk = ppatch%edges%start_blk(rl_start,1)
    i_endblk   = ppatch%edges%end_blk(rl_end,1)

    ! line and block indices of two neighboring cells
    iilc => ppatch%edges%cell_idx
    iibc => ppatch%edges%cell_blk

    ! loop through all patch edges (and blocks)
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(ppatch, jb, i_startblk, i_endblk,&
                   & i_startidx, i_endidx, rl_start, rl_end)
#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO jk = 1, n_zlev
        DO je = i_startidx, i_endidx
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          IF ( ppatch%patch_oce%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
            pupflux_e(je,jk,jb) =  0.5_wp*pvn_e(je,jk,jb)        &
              &        *( pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)) &
              &          +pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
  END SUBROUTINE central_hflux_oce
  !-------------------------------------------------------------------------
  !>
  !! First order upwind scheme for vertical tracer advection
  !!
  !! Calculation of vertical tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !! Initial revision by Jochen Foerstner, DWD (2008-05-15)
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized to height based vertical coordinate systems
  !! Modification by Stephan Lorenz, MPI (2010-09-07)
  !! - adapted to hydrostatic ocean core
  !!
  SUBROUTINE upwind_vflux_oce( ppatch, pvar_c, pw_c, pupflux_i )

    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch      !< patch on which computation is performed
    REAL(wp), INTENT(INOUT)  :: pvar_c(:,:,:)      !< advected cell centered variable
    REAL(wp), INTENT(INOUT)  :: pw_c(:,:,:)        !< vertical velocity on cells
    REAL(wp), INTENT(INOUT)  :: pupflux_i(:,:,:)   !< variable in which the upwind flux is stored
                                                   !< dim: (nproma,n_zlev+1,nblks_c)
    ! local variables
    ! height based but reversed (downward increasing depth) coordinate system,
    ! grid coefficient is negative (same as pressure based atmospheric coordinate system
    REAL(wp), PARAMETER :: zcoeff_grid = -1.0_wp
    INTEGER :: z_dolic

    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end
    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: jkm1                     !< jk - 1
    !-------------------------------------------------------------------------
    rl_start = 1
    rl_end   = min_rlcell
    i_startblk = ppatch%cells%start_blk(rl_start,1)
    i_endblk   = ppatch%cells%end_blk(rl_end,1)

    ! no fluxes at top boundary
    pupflux_i(:,1,:) = pvar_c(:,1,:)*pw_c(:,1,:) !0.0_wp
    ! no fluxes at bottom boundary
    pupflux_i(:,n_zlev+1,:) = 0.0_wp
    !
    DO jb =  i_startblk, i_endblk
      CALL get_indices_c(ppatch, jb, i_startblk, i_endblk,&
                   & i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
        z_dolic = ppatch%patch_oce%dolic_c(jc,jb)
        IF(z_dolic>=2)THEN
          DO jk = 2, n_zlev
            jkm1 = jk - 1
              ! calculate vertical tracer flux using upwind method
              pupflux_i(jc,jk,jb) =                 &
                &  laxfr_upflux_v( pw_c(jc,jk,jb),  &
                &                  pvar_c(jc,jkm1,jb), pvar_c(jc,jk,jb), zcoeff_grid )
            END DO
        ENDIF
      END DO
    END DO 
!    IF (p_pe == p_io .AND. ldbg) THEN
!       jk=c_k
!       jkm1=jk+1
!       WRITE(*,form4ar) &
!         &  '  pupflx_v E1 =', pupflux_i(c_i,jk  ,c_b), &
!         &  '    pvar_c(1) =', pvar_c   (c_i,jk  ,c_b), &
!         &  '    pvar_c(2) =', pvar_c   (c_i,jkm1,c_b)
!     END IF
  END SUBROUTINE upwind_vflux_oce
  !-------------------------------------------------------------------------
  !>
  !! First order upwind scheme for vertical tracer advection
  !!
  !! Calculation of vertical tracer fluxes 
  !!
  !! @par Revision History
  !!
  SUBROUTINE mimetic_vflux_oce( ppatch, pvar_c, pw_c, pupflux_i )

    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch      !< patch on which computation is performed
    REAL(wp), INTENT(INOUT)  :: pvar_c(:,:,:)      !< advected cell centered variable
    REAL(wp), INTENT(INOUT)  :: pw_c(:,:,:)        !< vertical velocity on cells
    REAL(wp), INTENT(INOUT)  :: pupflux_i(:,:,:)   !< variable in which the upwind flux is stored
                                                   !< dim: (nproma,n_zlev+1,nblks_c)
    ! local variables
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end
    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: jkm1, jkp1, z_dolic                    !< jk - 1
    REAL(wp) :: w_ave(n_zlev)
    REAL(wp) :: w_avep1(n_zlev)
    !-------------------------------------------------------------------------
    rl_start = 1
    rl_end   = min_rlcell
    i_startblk = ppatch%cells%start_blk(rl_start,1)
    i_endblk   = ppatch%cells%end_blk(rl_end,1)

    w_ave(:)  = 0.0_wp
    w_avep1(:)= 0.0_wp
    ! no fluxes at top boundary
    pupflux_i(:,1,:) = pvar_c(:,1,:)*pw_c(:,1,:) !0.0_wp
    ! no fluxes at bottom boundary
    pupflux_i(:,n_zlev+1,:) = 0.0_wp

    DO jb =  i_startblk, i_endblk
      CALL get_indices_c(ppatch, jb, i_startblk, i_endblk,&
                   & i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          z_dolic = ppatch%patch_oce%dolic_c(jc,jb)
          IF(z_dolic>=2)THEN
            DO jk = 2, z_dolic-1
              jkm1 = jk - 1
              jkp1 = jk + 1

              w_ave(jk)=0.5_wp*(pw_c(jc,jk,jb)+pw_c(jc,jkm1,jb))&
                    &*pvar_c(jc,jk,jb)

             w_avep1(jk)=0.5_wp*(pw_c(jc,jkp1,jb)+pw_c(jc,jk,jb))&
                    &*pvar_c(jc,jkp1,jb)

            pupflux_i(jc,jk,jb) = 0.5_wp*(w_ave(jk)  *ppatch%patch_oce%del_zlev_m(jk) &
                                &        +w_avep1(jk)*ppatch%patch_oce%del_zlev_m(jkp1))&
                                &        /ppatch%patch_oce%del_zlev_i(jk)
            END DO
        ENDIF 
      END DO 
    END DO
!     DO jb =  i_startblk, i_endblk
!       CALL get_indices_c(ppatch, jb, i_startblk, i_endblk,&
!                    & i_startidx, i_endidx, rl_start, rl_end)
!         DO jc = i_startidx, i_endidx
!           DO jk = 2, n_zlev
!            ! index of top & bottom half level
!            jkm1 = jk - 1
!            jkp1 = jk + 1
!              w_ave(jk)=0.5_wp*(pw_c(jc,jk,jb)+pw_c(jc,jkm1,jb))&
!                    &*pvar_c(jc,jk,jb)
!            END DO ! end cell loop
!        END DO ! end level loop
!      END DO ! end block loop
!     DO jb =  i_startblk, i_endblk
!       CALL get_indices_c(ppatch, jb, i_startblk, i_endblk,&
!                    & i_startidx, i_endidx, rl_start, rl_end)
!         DO jc = i_startidx, i_endidx
!          DO jk = 1, n_zlev-1
!            ! index of top & bottom half level
!            jkm1 = jk - 1
!            jkp1 = jk + 1
!             w_avep1(jk)=0.5_wp*(pw_c(jc,jkp1,jb)+pw_c(jc,jk,jb))&
!                    &*pvar_c(jc,jkp1,jb)
!            END DO ! end cell loop
!        END DO ! end level loop
!      END DO ! end block loop
!     DO jb =  i_startblk, i_endblk
!       CALL get_indices_c(ppatch, jb, i_startblk, i_endblk,&
!                    & i_startidx, i_endidx, rl_start, rl_end)
!         DO jc = i_startidx, i_endidx 
!          DO jk = 1, n_zlev-1
!            ! index of top & bottom half level
!            jkm1 = jk - 1
!            jkp1 = jk + 1
!            pupflux_i(jc,jk,jb) = 0.5_wp*(w_ave(jk)  *ppatch%patch_oce%del_zlev_m(jk) &
!                                &        +w_avep1(jk)*ppatch%patch_oce%del_zlev_m(jkp1))/ppatch%patch_oce%del_zlev_i(jk)
!            END DO ! end cell loop
!        END DO ! end level loop
!      END DO ! end block loop
  END SUBROUTINE mimetic_vflux_oce
 !-------------------------------------------------------------------------
  !>
  !!
  !! Calculation of central vertical tracer fluxes
  !!
  !! @par Revision History
  !! Petter Korn, MPI-M
  !!
  SUBROUTINE central_vflux_oce( ppatch, pvar_c, pw_c, pupflux_i )

    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch      !< patch on which computation is performed
    REAL(wp), INTENT(INOUT)  :: pvar_c(:,:,:)      !< advected cell centered variable
    REAL(wp), INTENT(INOUT)  :: pw_c(:,:,:)        !< vertical velocity on cells
    REAL(wp), INTENT(INOUT)  :: pupflux_i(:,:,:)   !< variable in which the upwind flux is stored
                                                   !< dim: (nproma,n_zlev+1,nblks_c)
    ! local variables
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end
    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: jkm1                     !< jk - 1
    INTEGER  :: z_dolic
    !-------------------------------------------------------------------------
    rl_start = 1
    rl_end   = min_rlcell
    i_startblk = ppatch%cells%start_blk(rl_start,1)
    i_endblk   = ppatch%cells%end_blk(rl_end,1)

    !fluxes at top boundary
    pupflux_i(:,1,:) = pvar_c(:,1,:)* pw_c(:,1,:)    !0.0_wp
    ! no fluxes at bottom boundary
    pupflux_i(:,n_zlev+1,:) = 0.0_wp

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(ppatch, jb, i_startblk, i_endblk,&
                      & i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          z_dolic = ppatch%patch_oce%dolic_c(jc,jb)
          IF(z_dolic>=2)THEN
            DO jk = 2, z_dolic
              ! index of top half level
              jkm1 = jk - 1

              ! calculate vertical tracer flux using upwind method
              pupflux_i(jc,jk,jb) = 0.5_wp*pw_c(jc,jk,jb) &
                &  * (pvar_c(jc,jkm1,jb)+ pvar_c(jc,jk,jb) )
! IF(jb==959.and.jc==10)THEN 
! IF(ppatch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary)THEN
! write(*,*)'vert fluxes:',jc,jk,jb,pupflux_i(jc,jk,jb),pw_c(jc,jk,jb),&
! & pvar_c(jc,jkm1,jb), pvar_c(jc,jk,jb) 
! ENDIF
! ENDIF 
          END DO
        ENDIF
      END DO
    END DO
  END SUBROUTINE central_vflux_oce
!------------------------------------------------------------------------


!!>
! !  SUBROUTINE advects the tracers present in the ocean model.
!
! @par Revision History
! Developed  by  Peter Korn, MPI-M (2010).
! 
SUBROUTINE advect_individual_tracer_ab_old(p_patch, trac_old,           &
                                     & p_os, G_n_c_h,G_nm1_c_h,G_nimd_c_h, &
                                     & G_n_c_v, G_nm1_c_v, G_nimd_c_v, &
                                     & bc_top_tracer, bc_bot_tracer,&
                                     & K_h, A_v,                    &
                                     & trac_new, timestep)


TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp)                          :: trac_old(:,:,:)
TYPE(t_hydro_ocean_state), TARGET :: p_os
REAL(wp) :: G_n_c_h   (nproma, n_zlev,   p_patch%nblks_c)  !G^n
REAL(wp) :: G_nm1_c_h (nproma, n_zlev,   p_patch%nblks_c)  !G^(n-1)
REAL(wp) :: G_nimd_c_h(nproma, n_zlev,   p_patch%nblks_c)  !G^(n+1/2)
REAL(wp) :: G_n_c_v   (nproma, n_zlev,   p_patch%nblks_c)  !G^n
REAL(wp) :: G_nm1_c_v (nproma, n_zlev,   p_patch%nblks_c)  !G^(n-1)
REAL(wp) :: G_nimd_c_v(nproma, n_zlev,   p_patch%nblks_c)  !G^(n+1/2)
REAL(wp) :: bc_top_tracer(nproma, p_patch%nblks_c)
REAL(wp) :: bc_bot_tracer(nproma, p_patch%nblks_c)
REAL(wp) :: K_h(:,:,:)                                  !horizontal mixing coeff
REAL(wp) :: A_v(:,:,:)                                   !vertical mixing coeff
REAL(wp) :: trac_new(:,:,:)                              !new tracer 
INTEGER  :: timestep                                     ! Actual timestep (to distinghuish initial step from others)

!Local variables
REAL(wp) :: delta_z!, delta_z2
!INTEGER  :: ctr, ctr_total
INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, rl_start_c, rl_end_c
INTEGER  :: jc, jk, jb, jkp1        !< index of edge, vert level, block
INTEGER  :: z_dolic
REAL(wp) :: z_adv_flux_h (nproma, n_zlev,   p_patch%nblks_e)  ! horizontal advective tracer flux
REAL(wp) :: z_adv_flux_v (nproma, n_zlev+1, p_patch%nblks_c)  ! vertical advective tracer flux
REAL(wp) :: z_div_adv_h (nproma, n_zlev,  p_patch%nblks_c)        ! horizontal tracer divergence
REAL(wp) :: z_div_adv_v (nproma, n_zlev,p_patch%nblks_c)        ! vertical tracer divergence
REAL(wp) :: z_div_diff_h (nproma, n_zlev,  p_patch%nblks_c)        ! horizontal tracer divergence
REAL(wp) :: z_div_diff_v (nproma, n_zlev,p_patch%nblks_c)        ! vertical tracer divergence

REAL(wp) :: z_diff_flux_h(nproma, n_zlev,  p_patch%nblks_e)   ! horizontal diffusive tracer flux
REAL(wp) :: z_diff_flux_v(nproma, n_zlev+1,p_patch%nblks_c)   ! vertical diffusive tracer flux


REAL(wp) :: z_transport_vn(nproma,n_zlev, p_patch%nblks_e)  ! horizontal transport velocity
REAL(wp) :: z_transport_w(nproma,n_zlev+1,p_patch%nblks_c)  ! vertical transport velocity

REAL(wp) :: z_mass_flux_h(nproma,n_zlev, p_patch%nblks_e)
REAL(wp) :: z_trac_c(nproma,n_zlev, p_patch%nblks_c)
REAL(wp) :: z_div_mass_flux_h(nproma,n_zlev, p_patch%nblks_c)
REAL(wp) :: z_h(nproma,n_zlev, p_patch%nblks_c)
REAL(wp) :: z_h2(nproma,n_zlev, p_patch%nblks_c)
REAL(wp) :: z_h_tmp(nproma,n_zlev, p_patch%nblks_c)

REAL(wp) :: z_G_n_c_v   (nproma, n_zlev, p_patch%nblks_c)
REAL(wp) :: z_G_nm1_c_v (nproma, n_zlev, p_patch%nblks_c)
REAL(wp) :: z_G_nimd_c_v(nproma, n_zlev, p_patch%nblks_c)
REAL(wp) :: max_val, min_val, dtime2

INTEGER            :: FLUX_CALCULATION = 3 
INTEGER, PARAMETER :: UPWIND = 1
INTEGER, PARAMETER :: CENTRAL= 2
INTEGER, PARAMETER :: MIMETIC= 3
REAL(wp)           :: z_tol
LOGICAL :: ldbg = .false.
TYPE(t_cartesian_coordinates):: z_vn_c(nproma,n_zlev,p_patch%nblks_c)
TYPE(t_cartesian_coordinates):: z_vn_c2(nproma,n_zlev,p_patch%nblks_c)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_tracer_advection:advect_individual_tracer')
!-------------------------------------------------------------------------------
z_tol= 0.0_wp!1.0E-13

dtime2= 1.0_wp*dtime

rl_start_c   = 1
rl_end_c     = min_rlcell
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

! ! rl_start_e   = 1
! ! rl_end_e     = min_rledge
! ! i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
! ! i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

z_adv_flux_h  = 0.0_wp
z_div_adv_h   = 0.0_wp
z_div_diff_h  = 0.0_wp
z_diff_flux_h = 0.0_wp
z_adv_flux_v  = 0.0_wp
z_div_adv_v   = 0.0_wp
z_div_diff_v  = 0.0_wp
z_diff_flux_v = 0.0_wp
z_mass_flux_h = 0.0_wp
z_h           = 0.0_wp
z_h2          = 0.0_wp
z_h_tmp       =0.0_wp
z_transport_w =0.0_wp
z_trac_c      =0.0_wp
z_div_mass_flux_h=0.0_wp
z_G_n_c_v        = 0.0_wp
z_G_nm1_c_v      = 0.0_wp
z_G_nimd_c_v     = 0.0_wp

DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
  &                   rl_start_c, rl_end_c)
  DO jk = 1, n_zlev
    delta_z = p_patch%patch_oce%del_zlev_m(jk)
    DO jc = i_startidx_c, i_endidx_c
      IF (jk == 1) delta_z = p_patch%patch_oce%del_zlev_m(jk)&
                         & + p_os%p_prog(nold(1))%h(jc,jb)
      z_h_tmp(jc,jk,jb) = delta_z 
      z_trac_c(jc,jk,jb)=trac_old(jc,jk,jb)*delta_z
    END DO
  END DO
END DO

!Step 1) Horizontal advection and diffusion
SELECT CASE(FLUX_CALCULATION)

CASE(UPWIND)
  z_transport_vn = ab_gam*p_os%p_prog(nnew(1))%vn + (1.0_wp-ab_gam)*p_os%p_prog(nold(1))%vn

  CALL upwind_hflux_oce( p_patch,        &
                       & z_h_tmp,        &
                       & z_transport_vn, &
                       & z_mass_flux_h )
  CALL upwind_hflux_oce( p_patch,        &
                       & z_trac_c,       &
                       & z_transport_vn, &
                       & z_adv_flux_h )
CASE(CENTRAL)
  z_transport_vn = ab_gam*p_os%p_prog(nnew(1))%vn + (1.0_wp-ab_gam)*p_os%p_prog(nold(1))%vn

  CALL central_hflux_oce( p_patch,       &
                       & z_h_tmp,        &
                       & z_transport_vn, &
                       & z_mass_flux_h )
  CALL central_hflux_oce( p_patch,       &
                       & z_trac_c,       &
                       & z_transport_vn, &
                       & z_adv_flux_h )
CASE(MIMETIC)

  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jk = 1, n_zlev
      DO jc = i_startidx_c, i_endidx_c

        z_vn_c(jc,jk,jb)%x(:) =0.0_wp
        z_vn_c2(jc,jk,jb)%x(:)=0.0_wp
        IF (jk == 1)THEN
          delta_z = p_os%p_prog(nold(1))%h(jc,jb) +p_patch%patch_oce%del_zlev_m(jk)!+p_os%p_prog(nold(1))%h(jc,jb)

          z_vn_c(jc,jk,jb)%x  = trac_old(jc,jk,jb)*delta_z*p_os%p_diag%p_vn(jc,jk,jb)%x&
          &*p_patch%cells%area(jc,jb)
          z_vn_c2(jc,jk,jb)%x = delta_z*p_os%p_diag%p_vn(jc,jk,jb)%x*p_patch%cells%area(jc,jb)
         ELSE
          z_vn_c(jc,jk,jb)%x = trac_old(jc,jk,jb)*p_os%p_diag%p_vn(jc,jk,jb)%x
          z_vn_c2(jc,jk,jb)%x= p_os%p_diag%p_vn(jc,jk,jb)%x
         ENDIF
      END DO
    END DO
  END DO

   CALL map_cell2edges( p_patch, z_vn_c, z_adv_flux_h )
   CALL map_cell2edges( p_patch, z_vn_c2, z_mass_flux_h )
END SELECT
 
!Step 4) calculate horizontal divergence
CALL div_oce( z_adv_flux_h, p_patch, z_div_adv_h)
CALL div_oce( z_mass_flux_h, p_patch, z_div_mass_flux_h)

DO jk = 1, n_zlev
write(*,*)'max/min adv-flux:',jk,&
&maxval(z_adv_flux_h(:,jk,:)),&
&minval(z_adv_flux_h(:,jk,:))
END DO

DO jk=1,n_zlev
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      z_h(jc,jk,jb) = z_h_tmp(jc,jk,jb) - dtime2*z_div_mass_flux_h(jc,jk,jb)
    END DO
  END DO!write(*,*)'max/min h_dummy:',maxval(z_h(:,jk,:)),minval(z_h(:,jk,:))
END DO

!The diffusion part
!Calculate horizontal diffusive flux
CALL tracer_diffusion_horz( p_patch, &
 &                          trac_old,&
 &                          p_os,    &
 &                          K_h,     & 
 &                          z_diff_flux_h)


CALL div_oce( z_diff_flux_h, p_patch, z_div_diff_h)


DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                   & i_startidx_c, i_endidx_c,&
                   & rl_start_c, rl_end_c)
  DO jk = 1, n_zlev
    DO jc = i_startidx_c, i_endidx_c
      G_n_c_h(jc,jk,jb) = z_div_diff_h(jc,jk,jb)
    END DO
  END DO
END DO

IF(timestep/=1)THEN
    G_nimd_c_h(:,:,:) = (1.5_wp+AB_const)* G_n_c_h(:,:,:)   &
        &             - (0.5_wp+AB_const)*G_nm1_c_h(:,:,:)
  ELSE
    G_nimd_c_h(:,:,:) = G_n_c_h(:,:,:)
ENDIF

IF(ldbg)THEN
  DO jk = 1, n_zlev
    write(*,*)'max/min G_T:',jk,&
    &maxval(G_n_c_h(:,jk,:)), minval(G_n_c_h(:,jk,:))
  END DO

  DO jk = 1, n_zlev
    write(*,*)'max/min G_T+1/2:',jk,maxval(G_nimd_c_h(:,jk,:)),&
                                   &minval(G_nimd_c_h(:,jk,:))
  END DO
  DO jk = 1, n_zlev
    write(*,*)'max/min adv flux & div h:',jk,maxval(z_adv_flux_h(:,jk,:)),&
                                           & minval(z_adv_flux_h(:,jk,:)),&
                                           & maxval(z_div_adv_h(:,jk,:)),&
                                           & minval(z_div_adv_h(:,jk,:))
  END DO
  DO jk = 1, n_zlev
    write(*,*)'max/min mass flux & div h:',jk, maxval(z_mass_flux_h(:,jk,:)),&
                                              & minval(z_mass_flux_h(:,jk,:)),&
                                             & maxval(z_div_mass_flux_h(:,jk,:)),&
                                             & minval(z_div_mass_flux_h(:,jk,:))
  END DO
ENDIF

!Final step: calculate new tracer values
DO jk = 1, n_zlev
  CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                   & i_startidx_c, i_endidx_c,&
                   & rl_start_c, rl_end_c)
    DO jb = i_startblk_c, i_endblk_c
      DO jc = i_startidx_c, i_endidx_c
        IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

          IF(jk==1)THEN
            delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+p_patch%patch_oce%del_zlev_m(top)
            !delta_z2  = p_os%p_prog(nnew(1))%h(jc,jb)+p_patch%patch_oce%del_zlev_m(top)
          ELSE
            delta_z  = p_patch%patch_oce%del_zlev_m(jk)
           !delta_z2  = p_patch%patch_oce%del_zlev_m(jk) 
          ENDIF
          trac_new(jc,jk,jb) = (z_trac_c(jc,jk,jb)             &
                             & -dtime2*(z_div_adv_h(jc,jk,jb)  &
                             &            -G_nimd_c_h(jc,jk,jb)))&
                             &/(z_h(jc,jk,jb)+z_tol)
        ENDIF
      END DO
    END DO
END DO

  G_nm1_c_h(:,:,:)  = G_n_c_h(:,:,:)
  G_nimd_c_h(:,:,:) = 0.0_wp
  G_n_c_h(:,:,:)    = 0.0_wp

  DO jk = 1, n_zlev
  max_val=0.0_wp
  min_val=100.0_wp
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
        IF(trac_new(jc,jk,jb)>max_val)&
          &max_val=trac_new(jc,jk,jb)
        IF(trac_new(jc,jk,jb)<min_val)&
          &min_val=trac_new(jc,jk,jb)

        IF(trac_new(jc,jk,jb)==0.0_wp)&
        & write(*,*)'tracer =0:',jc,jk,jb 
      ENDIF
     END DO
   END DO
   WRITE(*,*)'MAX/MIN tracer after horizontal transport:', jk, max_val, min_val
END DO

!-----------------------------------------------------------------------------------------------------------------------------
!Step 2) Vertical advection and diffusion


z_transport_w  = ab_gam*p_os%p_diag%w + (1.0_wp-ab_gam)*p_os%p_diag%w_old

SELECT CASE(FLUX_CALCULATION)

CASE(UPWIND)
  CALL upwind_vflux_oce( p_patch,      &
                       & trac_new,     &
                       & z_transport_w,&
                       & z_adv_flux_v )
CASE(CENTRAL)
  CALL central_vflux_oce( p_patch,     &
                       & trac_new,     &
                       & z_transport_w,&
                       & z_adv_flux_v )
CASE(MIMETIC)
  CALL mimetic_vflux_oce( p_patch,      &
                       & trac_new,     &
                       & z_transport_w,&
                       & z_adv_flux_v )
END SELECT

!calculate mass flux
DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
  &                   rl_start_c, rl_end_c)
  DO jc = i_startidx_c, i_endidx_c
    z_dolic = p_patch%patch_oce%dolic_c(jc,jb)
    IF(z_dolic/=0)THEN
      DO jk = 1, z_dolic
        z_h2(jc,jk,jb) = z_h_tmp(jc,jk,jb)&
        & - (z_transport_w(jc,jk,jb) -z_transport_w(jc,jk+1,jb))
      END DO
      z_h2(jc,z_dolic,jb) = z_h_tmp(jc,z_dolic,jb)
   ENDIF
  END DO
END DO

!The diffusion part
IF(expl_vertical_tracer_diff==1)THEN

  !vertical divergence is calculated for vertical advective fluxes
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      !interior: from one below surface to the ground
      z_dolic = p_patch%patch_oce%dolic_c(jc,jb)
      IF(z_dolic>0)THEN
        DO jk = 1, z_dolic
          jkp1 = jk + 1
          !positive vertical divergence in direction of w (upward positive)
          z_div_adv_v(jc,jk,jb) = &
          & (z_adv_flux_v(jc,jk,jb)-z_adv_flux_v(jc,jkp1,jb))

        END DO!level-loop
      ENDIF!(z_dolic>0)
    END DO!idx-loop
  END DO !blk-loop

!Final step: calculate new tracer values
DO jk = 1, n_zlev 
  CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                   & i_startidx_c, i_endidx_c,&
                   & rl_start_c, rl_end_c)
    DO jb = i_startblk_c, i_endblk_c
      DO jc = i_startidx_c, i_endidx_c
        IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          IF(jk==1)THEN
            delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+p_patch%patch_oce%del_zlev_m(top)
            !delta_z2  = p_os%p_prog(nnew(1))%h(jc,jb)+p_patch%patch_oce%del_zlev_m(top)
          ELSE
            delta_z  = p_patch%patch_oce%del_zlev_m(jk)  !delta_z2  = p_patch%patch_oce%del_zlev_m(jk) 
          ENDIF
            z_trac_c(jc,jk,jb) = (trac_new(jc,jk,jb)*delta_z     &
                               & -dtime2*z_div_adv_v(jc,jk,jb))  &
                               &/(z_h2(jc,jk,jb)+z_tol)
        ENDIF
      END DO
    END DO
END DO

  CALL tracer_diffusion_vert_impl( p_patch,        &
                              & z_trac_c(:,:,:),&
                              & bc_top_tracer,  & 
                              & bc_bot_tracer,  &
                              & A_v,            &
                              & trac_new(:,:,:))

ELSEIF(expl_vertical_tracer_diff==0)THEN

  !calculation of vertical diffusive fluxes
  CALL tracer_diffusion_vert_expl( p_patch,      &
                                &  trac_old, z_h, &
                                &  bc_top_tracer, &
                                &  bc_bot_tracer,  & 
                                &  A_v,           &
                                &  z_diff_flux_v)

  DO jk = 1, n_zlev+1
    write(*,*)'max/min adv & diff tracer flux v:',jk,&
    &maxval(z_adv_flux_v(:,jk,:)), minval(z_adv_flux_v(:,jk,:)),&
    &maxval(z_diff_flux_v(:,jk,:)),minval(z_diff_flux_v(:,jk,:))
  END DO
  DO jk = 1, n_zlev
    write(*,*)'max/min vertical mass flux:',jk,&
    &maxval(z_h2(:,jk,:)),minval(z_h2(:,jk,:))
  END DO

  !vertical divergence is calculated for vertical advective & diffusive fluxes
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      !interior: from one below surface to the ground
      z_dolic = p_patch%patch_oce%dolic_c(jc,jb)
      IF(z_dolic>0)THEN
        DO jk = 1, z_dolic
          jkp1 = jk + 1 

          IF(jk==1)THEN
            delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+p_patch%patch_oce%del_zlev_m(top)
          ELSE
            delta_z  = p_patch%patch_oce%del_zlev_m(jk)
          ENDIF
          !positive vertical divergence in direction of w (upward positive)
          z_div_adv_v(jc,jk,jb) = &
          & (z_adv_flux_v(jc,jk,jb)-z_adv_flux_v(jc,jkp1,jb))

          G_n_c_v(jc,jk,jb) = &   !z_div_diff_v(jc,jk,jb) = &
          & (z_diff_flux_v(jc,jk,jb)-z_diff_flux_v(jc,jkp1,jb))
        END DO!level-loop
      ENDIF!(z_dolic>0)
    END DO!idx-loop
  END DO !blk-loop

  IF(timestep/=1)THEN
    G_nimd_c_v(:,:,:) = (1.5_wp+AB_const)* G_n_c_v(:,:,:)   &
      &               - (0.5_wp+AB_const)*G_nm1_c_v(:,:,:)
  ELSE
    G_nimd_c_v(:,:,:) = G_n_c_v(:,:,:)
  ENDIF



  DO jk = 1, n_zlev
    write(*,*)'max/min div adv & diffusive  flux v:',&
    &jk,maxval(z_div_adv_v(:,jk,:)), minval(z_div_adv_v(:,jk,:)),&
    &maxval(G_n_c_v(:,jk,:)),minval(G_n_c_v(:,jk,:))
  END DO

  !Final step: calculate new tracer values
  DO jk = 1, n_zlev 
    CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                   & i_startidx_c, i_endidx_c,&
                   & rl_start_c, rl_end_c)
    DO jb = i_startblk_c, i_endblk_c
      DO jc = i_startidx_c, i_endidx_c
        IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

          IF(jk==1)THEN
            delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+p_patch%patch_oce%del_zlev_m(top)
            !delta_z2  = p_os%p_prog(nnew(1))%h(jc,jb)+p_patch%patch_oce%del_zlev_m(top)
          ELSE
            delta_z  = p_patch%patch_oce%del_zlev_m(jk)
           !delta_z2  = p_patch%patch_oce%del_zlev_m(jk) 
          ENDIF
            trac_new(jc,jk,jb) = (trac_new(jc,jk,jb)*delta_z     &
                               & -dtime2*(z_div_adv_v(jc,jk,jb)  &
                               & -G_nimd_c_v(jc,jk,jb)))         &
                               &/(z_h2(jc,jk,jb)+z_tol)
        ENDIF
      END DO
    END DO
  END DO

  G_nm1_c_v(:,:,:)  = G_n_c_v(:,:,:)
  G_nimd_c_v(:,:,:) = 0.0_wp
  G_n_c_v(:,:,:)    = 0.0_wp


ENDIF!(lvertical_diff_implicit)THEN

rl_start_c   = 1
rl_end_c     = min_rlcell
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

  DO jk = 1, n_zlev
  max_val=0.0_wp
  min_val=100.0_wp
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
        IF(trac_new(jc,jk,jb)>max_val)&
          &max_val=trac_new(jc,jk,jb)
        IF(trac_new(jc,jk,jb)<min_val)&
          &min_val=trac_new(jc,jk,jb)

        IF(trac_new(jc,jk,jb)==0.0_wp)&
        & write(*,*)'tracer =0:',jc,jk,jb 
      ENDIF
     END DO
   END DO
   WRITE(*,*)'MAX/MIN tracer after vertical transport:', jk, max_val, min_val
END DO
!    DO jk = 1, n_zlev
!     ctr=0
!     ctr_total=0
!     DO jb = i_startblk_c, i_endblk_c
!       CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
!                        & i_startidx_c, i_endidx_c,&
!                        & rl_start_c, rl_end_c)
!    DO jc = i_startidx_c, i_endidx_c
!      IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary )THEN
!      ctr_total=ctr_total+1
!      ENDIF
!  
! !    IF(trac_new(jc,jk,jb)>t_prof(jk)) THEN
! !    ctr=ctr+1
! !  write(*,*)'indices',jc,jb,jk, p_patch%patch_oce%lsm_oce_c(jc,top,jb),&
! ! ! !&trac_old(jc,jk,jb),
! ! ! &(delta_z/delta_z2)*trac_old(jc,top,jb),&
! !  & dtime*G_nimd_c(jc,top,jb),G_n_c_h(jc,top,jb),z_div_adv_h(jc,top,jb),&
! !  trac_new(jc,top,jb), ctr
! !  ENDIF
!        END DO
!    END DO
!    write(*,*)'counter > initial:',jk, ctr, ctr_total
!  END DO


END SUBROUTINE advect_individual_tracer_ab_old
! !   !-----------------------------------------------------------------------


END MODULE mo_oce_tracer_transport
