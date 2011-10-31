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
USE mo_math_utilities,            ONLY: t_cartesian_coordinates, gc2cc  
USE mo_math_constants,            ONLY: dbl_eps
USE mo_impl_constants,            ONLY: sea_boundary, &
  &                                     min_rlcell, min_rledge, min_rlcell,MIN_DOLIC !, &
!  &                                     max_char_length
USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer, idisc_scheme,    &
                                    &   ab_const, ab_gam, expl_vertical_tracer_diff,&
                                    &   iswm_oce, i_sfc_forcing_form
USE mo_physical_constants,        ONLY: tf
USE mo_parallel_config,           ONLY: nproma
USE mo_dynamics_config,           ONLY: nold, nnew 
USE mo_run_config,                ONLY: dtime
USE mo_oce_state,                 ONLY: t_hydro_ocean_state, v_base, is_initial_timestep
USE mo_model_domain,              ONLY: t_patch
USE mo_exception,                 ONLY: finish !, message_text, message
USE mo_oce_index,                 ONLY: print_mxmn, jkc, jkdim, ipl_src
USE mo_loopindices,               ONLY: get_indices_c, get_indices_e !, get_indices_v
USE mo_oce_boundcond,             ONLY: top_bound_cond_tracer!,&
USE mo_oce_physics
!USE mo_oce_forcing,               ONLY: t_sfc_flx
USE mo_sea_ice,                   ONLY: t_sfc_flx
USE mo_scalar_product,            ONLY:  map_cell2edges,map_edges2cell,map_edges2cell!,&
!                                       & map_cell2edges_upwind,map_cell2edges_upwind2 

USE mo_oce_math_operators,        ONLY: div_oce, grad_fd_norm_oce, grad_fd_norm_oce_2d
!USE mo_interpolation,             ONLY: t_int_state
!USE mo_oce_index,                 ONLY: c_i, c_b, c_k, ne_b, ne_i, nc_b, nc_i, form4ar
USE mo_advection_utils,           ONLY: laxfr_upflux, laxfr_upflux_v
USE mo_oce_diffusion,             ONLY: tracer_diffusion_horz, tracer_diffusion_vert_expl,&
                                      & tracer_diffusion_vert_impl, tracer_diffusion_vert_impl_hom
USE mo_oce_ab_timestepping_mimetic, ONLY: l_STAGGERED_TIMESTEP
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
INTEGER, PARAMETER :: MIMETIC_MIURA= 4
!INTEGER, PARAMETER :: MIN_DOLIC = 2
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
TYPE(t_sfc_flx), INTENT(INOUT)    :: p_sfc_flx
INTEGER                           :: timestep! Actual timestep (to distinghuish initial step from others)
!
!Local variables
INTEGER :: i_no_t
!REAL(wp) :: max_val, min_val
!INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
!CHARACTER(len=max_char_length), PARAMETER :: &
!       & routine = ('mo_tracer_advection:advect_tracer_ab')
!-------------------------------------------------------------------------------

IF(idisc_scheme==1)THEN!Mimetic discretization
  FLUX_CALCULATION = MIMETIC
ELSEIF(idisc_scheme==2)THEN!RBF discretization
  FLUX_CALCULATION = UPWIND
ENDIF

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
REAL(wp) :: trac_tmp(nproma,n_zlev, p_patch%nblks_c)
INTEGER  :: jk
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_tracer_advection:advect_individual_tracer')
!-------------------------------------------------------------------------------
delta_t= dtime
h_tmp    = 0.0_wp
trac_tmp = 0.0_wp

CALL advect_horizontal(p_patch, trac_old,           &
                     & p_os, G_n_c_h,G_nm1_c_h,G_nimd_c_h, &
                     & K_h,                    &
                     & trac_tmp, timestep, delta_t, h_tmp)
  !write(123,*)'-------------timestep---------------',timestep
DO jk = 1, n_zlev
  ipl_src=3  ! output print level (1-5, fix)
!   CALL print_mxmn('adv-horz tracer-old',jk,trac_old(:,:,:),n_zlev, &
!     &              p_patch%nblks_c,'trc',ipl_src)
!   CALL print_mxmn('adv-horz tracer-tmp',jk,trac_tmp(:,:,:),n_zlev, &
!     &              p_patch%nblks_c,'trc',ipl_src)
  write(*,*)'After horizontal max/min old-new tracer:',jk, maxval(trac_old(:,jk,:)),&
                                        & minval(trac_old(:,jk,:)),&
                                        & maxval(trac_tmp(:,jk,:)),&
                                        & minval(trac_tmp(:,jk,:))
 !write(123,*)'After horizontal max/min old-new tracer:',jk, maxval(trac_old(:,jk,:)),&
 !                                      & minval(trac_old(:,jk,:)),&
 !                                      & maxval(trac_tmp(:,jk,:)),&
 !                                      & minval(trac_tmp(:,jk,:))
 END DO

IF( iswm_oce /= 1) THEN


     CALL advect_vertical(p_patch, trac_tmp,              &
                         & p_os,                           &
                         & G_n_c_v, G_nm1_c_v, G_nimd_c_v, &
                         & bc_top_tracer, bc_bot_tracer,   &
                         & A_v,                            &
                         & trac_new, timestep, delta_t)!, h_tmp)
   DO jk = 1, n_zlev 
     write(*,*)'After vertical max/min old-new tracer:',jk, maxval(trac_old(:,jk,:)),&
                                        & minval(trac_old(:,jk,:)),&
                                       & maxval(trac_new(:,jk,:)),&
                                        & minval(trac_new(:,jk,:))
!      ipl_src=3  ! output print level (1-5, fix)
!      CALL print_mxmn('adv-vert tracer-old',jk,trac_old(:,:,:),n_zlev, &
!        &              p_patch%nblks_c,'trc',ipl_src)
!      CALL print_mxmn('adv-vert tracer-new',jk,trac_new(:,:,:),n_zlev, &
!        &              p_patch%nblks_c,'trc',ipl_src)
    END DO
ELSEIF( iswm_oce == 1) THEN

  trac_new=trac_tmp

  DO jk = 1, n_zlev
    ipl_src=3  ! output print level (1-5, fix)
    CALL print_mxmn('adv-vert tracer-old',jk,trac_old(:,:,:),n_zlev, &
      &              p_patch%nblks_c,'trc',ipl_src)
    CALL print_mxmn('adv-vert tracer-new',jk,trac_new(:,:,:),n_zlev, &
      &              p_patch%nblks_c,'trc',ipl_src)
  END DO

ENDIF

DO jk = 1, n_zlev
    ! #slo# - Temperature: tf=-1.9 deg, freezing of sea water, is possible:
! IF (minval(trac_new(:,jk,:))<tf) THEN
  IF (minval(trac_new(:,jk,:))<-4.0_wp) THEN
    write(*,*)'negative tracer', jk, minval(trac_new(:,jk,:)) 
    CALL finish(TRIM('mo_tracer_advection:advect_individual_tracer-h'), &
      &              'Negative tracer values') 
  ENDIF
END DO

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
REAL(wp) :: z_grad_T(nproma,n_zlev,p_patch%nblks_e)
!REAL(wp) :: max_val, min_val, dtime2
!REAL(wp) :: z_tol
LOGICAL  :: ldbg = .TRUE.
!LOGICAL  ::  L_MPDATA_AFTERBURNER
TYPE(t_cartesian_coordinates):: z_vn_c(nproma,n_zlev,p_patch%nblks_c)
TYPE(t_cartesian_coordinates):: z_vn_c2(nproma,n_zlev,p_patch%nblks_c)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_tracer_advection:advect_individual_tracer')
!-------------------------------------------------------------------------------
!z_tol= 1.0E-13

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
z_grad_T=0.0_wp

!Generate 3D array of cell heights and of cell height times tracer content
IF ( iswm_oce /= 1) THEN
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jk = 1, n_zlev
      delta_z = v_base%del_zlev_m(jk)
      DO jc = i_startidx_c, i_endidx_c

        IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          IF (jk == 1) delta_z = v_base%del_zlev_m(jk)&
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
        IF ( v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
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

CASE(UPWIND, MIMETIC)
  !produce weighted transport velocity
  !  z_transport_vn = ab_gam*p_os%p_prog(nnew(1))%vn + (1.0_wp-ab_gam)*p_os%p_prog(nold(1))%vn
  IF(l_STAGGERED_TIMESTEP)THEN
    z_transport_vn = p_os%p_prog(nnew(1))%vn 
  ELSEIF(.NOT.l_STAGGERED_TIMESTEP)THEN 
    z_transport_vn = p_os%p_prog(nold(1))%vn
   ENDIF
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
  !  z_transport_vn = ab_gam*p_os%p_prog(nnew(1))%vn + (1.0_wp-ab_gam)*p_os%p_prog(nold(1))%vn
  IF(l_STAGGERED_TIMESTEP)THEN
    z_transport_vn = p_os%p_prog(nnew(1))%vn 
  ELSEIF(.NOT.l_STAGGERED_TIMESTEP)THEN 
    z_transport_vn = p_os%p_prog(nold(1))%vn
   ENDIF

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
!    CASE(MIMETIC)
  !produce weighted transport velocity
  !  z_transport_vn = ab_gam*p_os%p_prog(nnew(1))%vn + (1.0_wp-ab_gam)*p_os%p_prog(nold(1))%vn
!   IF(l_STAGGERED_TIMESTEP)THEN
!     z_transport_vn = p_os%p_prog(nnew(1))%vn 
!   ELSEIF(.NOT.l_STAGGERED_TIMESTEP)THEN 
!     z_transport_vn = p_os%p_prog(nold(1))%vn
!    ENDIF
!    CALL map_edges2cell( p_patch, z_transport_vn, z_vn_c, p_os%p_diag%h_e)
!    DO jb = i_startblk_c, i_endblk_c
!      CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
!      &                   rl_start_c, rl_end_c)
!      DO jk = 1, n_zlev
!        DO jc = i_startidx_c, i_endidx_c
!           z_vn_c2(jc,jk,jb)%x= z_h_tmp(jc,jk,jb) *z_vn_c(jc,jk,jb)%x
!           z_vn_c(jc,jk,jb)%x = z_trac_c(jc,jk,jb)*z_vn_c(jc,jk,jb)%x
! 
!       END DO
!     END DO
!   END DO
!     CALL map_cell2edges( p_patch, z_vn_c, z_adv_flux_h )
!     CALL map_cell2edges( p_patch, z_vn_c2, z_mass_flux_h )
  CASE(MIMETIC_MIURA)
  !produce weighted transport velocity
  !  z_transport_vn = ab_gam*p_os%p_prog(nnew(1))%vn + (1.0_wp-ab_gam)*p_os%p_prog(nold(1))%vn
  IF(l_STAGGERED_TIMESTEP)THEN
    z_transport_vn = p_os%p_prog(nnew(1))%vn
  ELSEIF(.NOT.l_STAGGERED_TIMESTEP)THEN 
    z_transport_vn = p_os%p_prog(nold(1))%vn
  ENDIF

  !upwind estimate of mass flux
  CALL mimetic_miura_hflux_oce( p_patch,        &
                       & z_h_tmp,        &
                       & z_transport_vn, &    
                       & p_os%p_diag%p_vn,&
                       & z_mass_flux_h )
  !upwind estimate of tracer flux
  CALL mimetic_miura_hflux_oce( p_patch,        &
                       & z_trac_c,       &
                       & z_transport_vn, &    
                       & p_os%p_diag%p_vn,&
                       & z_adv_flux_h )


END SELECT
 
!Step 4) calculate horizontal divergence of mass and tracer flux
CALL div_oce( z_adv_flux_h, p_patch, z_div_adv_h)
CALL div_oce( z_mass_flux_h, p_patch, z_div_mass_flux_h)
 ipl_src=6  ! output print level (1-5, fix)
 DO jk = 1, n_zlev
   CALL print_mxmn('adv-flux',jk,z_adv_flux_h(:,:,:),n_zlev, &
     &              p_patch%nblks_e,'trc',ipl_src)
 END DO
 DO jk = 1, n_zlev
   CALL print_mxmn('div adv-flux',jk,z_div_adv_h(:,:,:),n_zlev, &
     &              p_patch%nblks_c,'trc',ipl_src)
 END DO
 DO jk = 1, n_zlev
   CALL print_mxmn('mass-flux',jk,z_mass_flux_h(:,:,:),n_zlev, &
     &              p_patch%nblks_e,'trc',ipl_src)
 END DO
 DO jk = 1, n_zlev
   CALL print_mxmn('div mass-flux',jk,z_div_mass_flux_h(:,:,:),n_zlev, &
     &              p_patch%nblks_c,'trc',ipl_src)
 END DO
! IF(ldbg)THEN
!  DO jk = 1, n_zlev
!  write(*,*)'max/min adv-flux:',jk,&
!  &maxval(z_adv_flux_h(:,jk,:)),&
!  &minval(z_adv_flux_h(:,jk,:))
!  END DO
!  DO jk = 1, n_zlev
!  write(*,*)'max/min mass-flux:',jk,&
!  &maxval(z_mass_flux_h(:,jk,:)),&
!  &minval(z_mass_flux_h(:,jk,:))
!  END DO
!  ENDIF
!calculate (dummy) height consistent with divergence of mass flux
DO jk=1,n_zlev
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
        z_h(jc,jk,jb) = z_h_tmp(jc,jk,jb) - delta_t*z_div_mass_flux_h(jc,jk,jb)
      ELSE
        z_h(jc,jk,jb) = 0.0_wp
      ENDIF
    END DO
  END DO
!write(*,*)'max/min h_dummy:',jk,maxval(z_h(:,jk,:)),minval(z_h(:,jk,:))
END DO 
!IF(ldbg)THEN
!  DO jk = 1, n_zlev
!   write(*,*)'max/min dummy height:',jk,&
!   &maxval(z_h(:,jk,:)),&
!   &minval(z_h(:,jk,:))
!   END DO
!  ENDIF
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
      IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
        G_n_c_h(jc,jk,jb) = z_div_diff_h(jc,jk,jb)
        !write(*,*)'div',jc,jk,jb,G_n_c_h(jc,jk,jb),z_div_diff_h(jc,jk,jb)
      ENDIF
    END DO
  END DO
END DO

IF(is_initial_timestep(timestep))THEN
    G_nimd_c_h(:,:,:) = G_n_c_h(:,:,:)
  ELSE
    G_nimd_c_h(:,:,:) = (1.5_wp+AB_const) * G_n_c_h(:,:,:)   &
        &             - (0.5_wp+AB_const) * G_nm1_c_h(:,:,:)
ENDIF

ipl_src=5  ! output print level (1-5, fix)
DO jk = 1, n_zlev
  CALL print_mxmn('G-T',jk,G_n_c_h(:,:,:),n_zlev, &
    &              p_patch%nblks_c,'trc',ipl_src)
END DO
DO jk = 1, n_zlev
  CALL print_mxmn('adv flux 2',jk,z_adv_flux_h(:,:,:),n_zlev, &
    &              p_patch%nblks_c,'trc',ipl_src)
END DO
DO jk = 1, n_zlev
  CALL print_mxmn('div adv-flux',jk,z_div_adv_h(:,:,:),n_zlev, &
    &              p_patch%nblks_e,'trc',ipl_src)
END DO
DO jk = 1, n_zlev
  CALL print_mxmn('mass-flux',jk,z_mass_flux_h(:,:,:),n_zlev, &
    &              p_patch%nblks_c,'trc',ipl_src)
END DO
DO jk = 1, n_zlev
  CALL print_mxmn('div mass-flux',jk,z_div_mass_flux_h(:,:,:),n_zlev, &
    &              p_patch%nblks_e,'trc',ipl_src)
END DO
!IF(ldbg)THEN
!  DO jk = 1, n_zlev
!    write(*,*)'max/min G_T:',jk,&
!    &maxval(G_n_c_h(:,jk,:)), minval(G_n_c_h(:,jk,:))
!  END DO
!    DO jk = 1, n_zlev
!      write(*,*)'max/min G_T+1/2:',jk,maxval(G_nimd_c_h(:,jk,:)),&
!                                     &minval(G_nimd_c_h(:,jk,:))
!    END DO
!  DO jk = 1, n_zlev
!    write(*,*)'max/min adv flux & div h:',jk,maxval(z_adv_flux_h(:,jk,:)),&
!                                           & minval(z_adv_flux_h(:,jk,:)),&
!                                           & maxval(z_div_adv_h(:,jk,:)),&
!                                           & minval(z_div_adv_h(:,jk,:))
!  END DO
!  DO jk = 1, n_zlev
!    write(*,*)'max/min mass flux & div h:',jk, maxval(z_mass_flux_h(:,jk,:)),&
!                                             & minval(z_mass_flux_h(:,jk,:)),&
!                                             & maxval(z_div_mass_flux_h(:,jk,:)),&
!                                             & minval(z_div_mass_flux_h(:,jk,:))
!  END DO
!ENDIF

!Final step: calculate new tracer values
!write(123,*)'timestep----------------------',timestep
DO jk = 1, n_zlev
  CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                   & i_startidx_c, i_endidx_c,&
                   & rl_start_c, rl_end_c)
    DO jb = i_startblk_c, i_endblk_c
      DO jc = i_startidx_c, i_endidx_c
        IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

!         IF(jk==1)THEN
!           delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+v_base%del_zlev_m(top)
!           !delta_z2  = p_os%p_prog(nnew(1))%h(jc,jb)+v_base%del_zlev_m(top)
!         ELSE
!           delta_z  = v_base%del_zlev_m(jk)
!          !delta_z2  = v_base%del_zlev_m(jk) 
!         ENDIF

          IF(z_h(jc,jk,jb)/=0.0_wp)THEN
             trac_new(jc,jk,jb) = (z_trac_c(jc,jk,jb)              &
                                & -delta_t*(z_div_adv_h(jc,jk,jb)  &
                                & -G_nimd_c_h(jc,jk,jb)))          &
                                &/z_h(jc,jk,jb)
! write(123,*)'H-adv: tracer_old_ trac_new',jc,jk,jb,trac_old(jc,jk,jb),trac_new(jc,jk,jb)!,&
! !&z_h(jc,jk,jb)
          ELSE
            trac_new(jc,jk,jb) = 0.0_wp
          ENDIF

        ENDIF
      END DO
    END DO
END DO

! IF(L_MPDATA_AFTERBURNER)THEN
! CALL mpdata_afterburner(p_patch, trac_new,           &
!                      & p_os,
!                      & trac_new, timestep, delta_t, h_tmp)
! ENDIF
! DO jk = 1, n_zlev
!   CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
!                    & i_startidx_c, i_endidx_c,&
!                    & rl_start_c, rl_end_c)
!     DO jb = i_startblk_c, i_endblk_c
!       DO jc = i_startidx_c, i_endidx_c
!       IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! 
! !           IF(jk==1)THEN
!             delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+v_base%del_zlev_m(top)
!             !delta_z2  = p_os%p_prog(nnew(1))%h(jc,jb)+v_base%del_zlev_m(top)
! !           ELSE
! !             delta_z  = v_base%del_zlev_m(jk)
! !            !delta_z2  = v_base%del_zlev_m(jk) 
! !           ENDIF
! 
!           IF(z_h(jc,jk,jb)/=0.0_wp)THEN
!             trac_new(jc,jk,jb) = trac_new(jc,jk,jb)             &
!                                & +delta_t*G_nimd_c_h(jc,jk,jb)&
!                                &/z_h(jc,jk,jb)
! ! write(123,*)'trac calc:',jc,jb,jk,&
! ! &trac_new(jc,jk,jb),&
! ! &trac_old(jc,jk,jb),&
! ! &z_trac_c(jc,jk,jb),&
! ! &z_div_adv_h(jc,jk,jb),&
! ! &G_nimd_c_h(jc,jk,jb),&
! ! &z_h(jc,jk,jb),&
! ! &z_trac_c(jc,jk,jb)/z_h(jc,jk,jb),&
! ! &delta_t*(z_div_adv_h(jc,jk,jb)-G_nimd_c_h(jc,jk,jb))/z_h(jc,jk,jb)
!           ELSE
!             trac_new(jc,jk,jb) = 0.0_wp
!           ENDIF
! 
!         ENDIF
!       END DO
!     END DO
! END DO

!   z_trac_c = trac_new
!  
!   CALL elad(p_patch, trac_old, z_trac_c, p_os)
! 
!   DO jk=1,n_zlev
!    write(*,*)'After elad:',maxval(z_trac_c(:,jk,:)),&
!                          & minval(z_trac_c(:,jk,:))
!   END DO
!   trac_new = z_trac_c

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
REAL(wp) :: delta_z
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
!REAL(wp) :: z_inv_h(nproma,n_zlev, p_patch%nblks_c)
REAL(wp) :: z_G_n_c_v   (nproma, n_zlev, p_patch%nblks_c)
REAL(wp) :: z_G_nm1_c_v (nproma, n_zlev, p_patch%nblks_c)
REAL(wp) :: z_G_nimd_c_v(nproma, n_zlev, p_patch%nblks_c)
!REAL(wp) :: max_val, min_val!, dtime2
!REAL(wp) :: z_tol
LOGICAL  :: ldbg = .TRUE.
!TYPE(t_cartesian_coordinates):: z_vn_c(nproma,n_zlev,p_patch%nblks_c)
!TYPE(t_cartesian_coordinates):: z_vn_c2(nproma,n_zlev,p_patch%nblks_c)
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_tracer_advection:advect_individual_tracer')
!-------------------------------------------------------------------------------
!z_tol= 0.0_wp!1.0E-13

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
!z_inv_h       = 0.0_wp
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
    delta_z = v_base%del_zlev_m(jk)
    DO jc = i_startidx_c, i_endidx_c
      IF (jk == 1) THEN
        delta_z = v_base%del_zlev_m(jk)&
              & + p_os%p_prog(nold(1))%h(jc,jb)
      ENDIF
      z_trac_c(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z
    END DO
  END DO
END DO


SELECT CASE(FLUX_CALCULATION)

CASE(UPWIND,MIMETIC)
  !Produce time-weighted vertical transport velocity
  !z_transport_w  = ab_gam*p_os%p_diag%w + (1.0_wp-ab_gam)*p_os%p_diag%w_old
  IF(l_STAGGERED_TIMESTEP)THEN
    z_transport_w = p_os%p_diag%w
  ELSEIF(.NOT.l_STAGGERED_TIMESTEP)THEN 
    z_transport_w = p_os%p_diag%w_old
  ENDIF

  CALL upwind_vflux_oce( p_patch,       &
                       & trac_in,       &
                       & z_transport_w, & 
                       & bc_top_tracer, &
                       & z_adv_flux_v )
CASE(CENTRAL)
  !Produce time-weighted vertical transport velocity
  !z_transport_w  = ab_gam*p_os%p_diag%w + (1.0_wp-ab_gam)*p_os%p_diag%w_old
  IF(l_STAGGERED_TIMESTEP)THEN
    z_transport_w = p_os%p_diag%w
  ELSEIF(.NOT.l_STAGGERED_TIMESTEP)THEN 
    z_transport_w = p_os%p_diag%w_old
  ENDIF

  CALL central_vflux_oce( p_patch,     &
                       & trac_in,      &
                       & z_transport_w,&
                       & z_adv_flux_v )
!   CASE(MIMETIC)
!   z_transport_w  = ab_gam*p_os%p_diag%w + (1.0_wp-ab_gam)*p_os%p_diag%w_old
!   !z_transport_w  = p_os%p_diag%w 
!   IF(l_STAGGERED_TIMESTEP)THEN
!     z_transport_w = p_os%p_diag%w
!   ELSEIF(.NOT.l_STAGGERED_TIMESTEP)THEN 
!     z_transport_w = p_os%p_diag%w_old
!   ENDIF
!   CALL mimetic_vflux_oce( p_patch,      &
!                         & trac_in,       &
!                         & z_transport_w, &
!                         & z_adv_flux_v )
END SELECT
 
!divergence is calculated for advective fluxes
DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
  &                   rl_start_c, rl_end_c)
  DO jc = i_startidx_c, i_endidx_c
    !interior: from one below surface to the ground
    z_dolic = v_base%dolic_c(jc,jb)
    IF(z_dolic>=MIN_DOLIC)THEN
      DO jk = 1, z_dolic
        ! positive vertical divergence in direction of w (upward positive)
        z_div_adv_v(jc,jk,jb) = z_adv_flux_v(jc,jk,jb) - z_adv_flux_v(jc,jk+1,jb) 
      END DO
    ENDIF
  END DO
END DO

!calculate mass flux
DO jb = i_startblk_c, i_endblk_c
  CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
  &                   rl_start_c, rl_end_c)
  DO jc = i_startidx_c, i_endidx_c
    z_dolic = v_base%dolic_c(jc,jb)
    IF(z_dolic>=MIN_DOLIC)THEN
      z_h(jc,1,jb) = (v_base%del_zlev_m(1) + p_os%p_prog(nold(1))%h(jc,jb))&
                   &- dtime*(z_transport_w(jc,1,jb)-z_transport_w(jc,2,jb))

      !z_inv_h(jc,1,jb) = 1.0_wp/z_h(jc,1,jb)
      DO jk = 2, z_dolic
        z_h(jc,jk,jb) = v_base%del_zlev_m(jk)&
                      & - dtime*(z_transport_w(jc,jk,jb)-z_transport_w(jc,jk+1,jb))
        !z_inv_h(jc,jk,jb) = 1.0_wp/z_h(jc,jk,jb)
      END DO
    ELSE
      z_h(jc,1,jb)     =0.0_wp
      !z_inv_h(jc,1,jb) = 0.0_wp
    ENDIF
  END DO
END DO

ipl_src=5  ! output print level (1-5, fix)
DO jk = 1, n_zlev
  CALL print_mxmn('adv flux_v',jk,z_adv_flux_v(:,:,:),n_zlev+1, &
    &              p_patch%nblks_c,'trc',ipl_src)
END DO
DO jk = 1, n_zlev
  CALL print_mxmn('div adv-flux_v',jk,z_div_adv_v(:,:,:),n_zlev, &
    &              p_patch%nblks_c,'trc',ipl_src)
END DO
DO jk = 1, n_zlev
  CALL print_mxmn('vertical mass flux',jk,z_h(:,:,:),n_zlev, &
    &              p_patch%nblks_c,'trc',ipl_src)
END DO

!IF(ldbg)THEN
!  DO jk = 1, n_zlev
!    write(*,*)'max/min adv flux & div:',jk,&
!    &maxval(z_adv_flux_v(:,jk,:)), minval(z_adv_flux_v(:,jk,:)),&
!    &maxval(z_div_adv_v(:,jk,:)),minval(z_div_adv_v(:,jk,:))
!  END DO
!  DO jk = 1, n_zlev
!    write(*,*)'max/min vertical mass flux:',jk,&
!    &maxval(z_h(:,jk,:)),minval(z_h(:,jk,:))
!  END DO
!ENDIF 



!Case: Implicit Vertical diffusion
IF(expl_vertical_tracer_diff==1)THEN

  SELECT CASE(i_sfc_forcing_form) 

  CASE(0)!surface forcing applied as top boundary condition to vertical diffusion


    !Add advective part to old tracer
    DO jb = i_startblk_c, i_endblk_c
      CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                       & i_startidx_c, i_endidx_c,&
                       & rl_start_c, rl_end_c)
        DO jc = i_startidx_c, i_endidx_c
          z_dolic = v_base%dolic_c(jc,jb)
          IF(z_dolic>=MIN_DOLIC)THEN

            DO jk = 1, z_dolic
            !IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
              IF(jk==1)THEN
                delta_z  = v_base%del_zlev_m(top)+p_os%p_prog(nold(1))%h(jc,jb)
              ELSE
                delta_z  = v_base%del_zlev_m(jk)
              ENDIF
              G_n_c_v(jc,jk,jb) = (trac_in(jc,jk,jb)*delta_z      &
                                 & -delta_t*z_div_adv_v(jc,jk,jb))&
                                 &/z_h(jc,jk,jb)
            END DO
          ELSE
            G_n_c_v(jc,:,jb) = 0.0_wp
          ENDIF
      END DO
    END DO
    IF( is_initial_timestep(timestep))THEN
      G_nimd_c_v(:,:,:) = G_n_c_v(:,:,:)
    ELSE
      G_nimd_c_v(:,:,:) = (1.5_wp+AB_const)* G_n_c_v(:,:,:)   &
        &               - (0.5_wp+AB_const)*G_nm1_c_v(:,:,:)
    ENDIF

    ipl_src=5  ! output print level (1-5, fix)
    DO jk = 1, n_zlev
      CALL print_mxmn('bef impl v-trc:G_n',jk,G_n_c_v(:,:,:),n_zlev, &
      &              p_patch%nblks_c,'trc',ipl_src)
    END DO
    DO jk = 1, n_zlev
      CALL print_mxmn('bef.impl v-trc:Gnimd',jk,G_nimd_c_v(:,:,:),n_zlev, &
      &              p_patch%nblks_c,'trc',ipl_src)
    END DO
    DO jk = 1, n_zlev
      CALL print_mxmn('bef.impl adv-flux-v',jk,z_adv_flux_v(:,:,:),n_zlev+1, &
      &              p_patch%nblks_c,'trc',ipl_src)
    END DO

    !IF(ldbg)THEN
    !  DO jk = 1, n_zlev
    ! write(*,*)'before impl: max/min vtracer G_n:G_n+1/2:',jk,&
    ! &maxval(G_n_c_v(:,jk,:)), minval(G_n_c_v(:,jk,:)),& 
    ! &maxval(G_nimd_c_v(:,jk,:)), minval(G_nimd_c_v(:,jk,:))

    ! write(123,*)'before impl: max/min vtracer G_n:G_n+1/2:',jk,&
    ! &maxval(G_n_c_v(:,jk,:)), minval(G_n_c_v(:,jk,:)),& 
    ! &maxval(G_nimd_c_v(:,jk,:)), minval(G_nimd_c_v(:,jk,:))

      !END DO
    !ENDIF

    !calculate vert diffusion impicit: result is stored in trac_out
      CALL tracer_diffusion_vert_impl( p_patch,             &
                                   & G_nimd_c_v,            &!& G_n_c_v,&
                                   & bc_top_tracer,         & 
                                   & bc_bot_tracer,         &
                                   & p_os%p_prog(nold(1))%h,&
                                   & A_v,                   &
                                   & trac_out(:,:,:))
  CASE(1)!=1: surface forcing applied as volume forcing at rhs, i.e.part of explicit term in momentum and tracer eqs.
         !    in this case, top boundary ondition of vertical Laplacians are homogeneous

    !Add advective part and boundary condition to old tracer
    DO jb = i_startblk_c, i_endblk_c
      CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                       & i_startidx_c, i_endidx_c,             &
                       & rl_start_c, rl_end_c)
        DO jc = i_startidx_c, i_endidx_c
          z_dolic = v_base%dolic_c(jc,jb)
          IF(z_dolic>=MIN_DOLIC)THEN

            DO jk = 1, z_dolic
            !IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
              IF(jk==1)THEN
                delta_z  = v_base%del_zlev_m(top)+p_os%p_prog(nold(1))%h(jc,jb)
              ELSE
                delta_z  = v_base%del_zlev_m(jk)
              ENDIF
              G_n_c_v(jc,jk,jb) = (trac_in(jc,jk,jb)*delta_z       &
                                 & -delta_t*z_div_adv_v(jc,jk,jb)  &
                                 & + delta_t*bc_top_tracer(jc,jb)) &
                                 &/z_h(jc,jk,jb)
            END DO
          ELSE
            G_n_c_v(jc,:,jb) = 0.0_wp
          ENDIF
      END DO
    END DO
    IF( is_initial_timestep(timestep))THEN
      G_nimd_c_v(:,:,:) = G_n_c_v(:,:,:)
    ELSE
      G_nimd_c_v(:,:,:) = (1.5_wp+AB_const)* G_n_c_v(:,:,:)   &
        &               - (0.5_wp+AB_const)*G_nm1_c_v(:,:,:)
    ENDIF

    CALL tracer_diffusion_vert_impl_hom( p_patch,         &
                                 & G_nimd_c_v,            &!& G_n_c_v,               &
                                 & p_os%p_prog(nold(1))%h,&
                                 & A_v,                   &
                                 & trac_out(:,:,:))

  END SELECT

!vertival diffusion is calculated explicitely
ELSEIF(expl_vertical_tracer_diff==0)THEN

  CALL tracer_diffusion_vert_expl( p_patch,       &
                                & trac_in,        &! z_trac_c,      &
                                &  z_h,           &
                                &  bc_top_tracer, &
                                &  bc_bot_tracer, & 
                                &  A_v,           &
                                &  z_div_diff_v)
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      !interior: from one below surface to the ground
      z_dolic = v_base%dolic_c(jc,jb)
      IF(z_dolic>=MIN_DOLIC)THEN
        DO jk = 1, z_dolic
           delta_z  = v_base%del_zlev_m(jk)
           ! positive vertical divergence in direction of w (upward positive)
           G_n_c_v(jc,jk,jb) = z_div_adv_v(jc,jk,jb)/z_h(jc,jk,jb) - z_div_diff_v(jc,jk,jb)

        END DO
      ENDIF
    END DO
  END DO

  ipl_src=5  ! output print level (1-5, fix)
! DO jk = 1, n_zlev
!   CALL print_mxmn('diff flux_v',jk,z_diff_flux_v(:,:,:),n_zlev+1, &
!     &              p_patch%nblks_c,'trc',ipl_src)
! END DO
  DO jk = 1, n_zlev
    CALL print_mxmn('div diff-flux_v',jk,z_div_diff_v(:,:,:),n_zlev, &
    &              p_patch%nblks_c,'trc',ipl_src)
  END DO

  IF( is_initial_timestep(timestep))THEN
    G_nimd_c_v(:,:,:) = G_n_c_v(:,:,:)
  ELSE
    G_nimd_c_v(:,:,:) = (1.5_wp+AB_const)* G_n_c_v(:,:,:)   &
      &               - (0.5_wp+AB_const)*G_nm1_c_v(:,:,:)
  ENDIF


  !Add advective and diffusive part to old tracer
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
                     & i_startidx_c, i_endidx_c,&
                     & rl_start_c, rl_end_c)
      DO jc = i_startidx_c, i_endidx_c
        z_dolic = v_base%dolic_c(jc,jb)
        IF(z_dolic>=MIN_DOLIC)THEN
          !top  level
          jk=1
          delta_z = (p_os%p_prog(nold(1))%h(jc,jb)+v_base%del_zlev_m(jk))&
                  &/z_h(jc,jk,jb) 
!           z_trac_c(jc,jk,jb) = (trac_in(jc,jk,jb)*delta_z      &
!                              & -delta_t*z_div_adv_v(jc,jk,jb)) &
!                              &/z_h(jc,jk,jb)
!            trac_out(jc,jk,jb) = z_trac_c(jc,jk,jb)               &
!                               & +delta_t*z_div_diff_v(jc,jk,jb)  
!           trac_out(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z      &
!                              & -delta_t*(z_div_adv_v(jc,jk,jb) &
!                              &/z_h(jc,jk,jb)                  &
!                              &-z_div_diff_v(jc,jk,jb))
          trac_out(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z      &
                             & -delta_t*G_nimd_c_v(jc,jk,jb) 

          !interior down to bottom
          DO jk = 2, z_dolic
            delta_z = v_base%del_zlev_m(jk)/z_h(jc,jk,jb)
!             z_trac_c(jc,jk,jb) = (trac_in(jc,jk,jb)*delta_z     &
!                                & -delta_t*z_div_adv_v(jc,jk,jb))  &
!                                &/z_h(jc,jk,jb) 
!            trac_out(jc,jk,jb) = z_trac_c(jc,jk,jb)                &
!                                 & +delta_t*z_diff_flux_v(jc,jk,jb)
!             trac_out(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z      &
!                                & -delta_t*(z_div_adv_v(jc,jk,jb) &
!                                &/z_h(jc,jk,jb)                  &
!                                &-z_div_diff_v(jc,jk,jb))
            trac_out(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z      &
                               & -delta_t*G_nimd_c_v(jc,jk,jb) 
          END DO
        ELSE
          trac_out(jc,:,jb) = 0.0_wp
        ENDIF
    END DO
  END DO
ENDIF!(lvertical_diff_implicit)THEN

END SUBROUTINE advect_vertical
!-----------------------------------------------------------------------
!
!
!>
!! !  SUBROUTINE advects horizontally the tracers present in the ocean model.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!! 
SUBROUTINE elad(p_patch, trac_old, trac_new, p_os)
!
!
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp) :: trac_old(:,:,:)
REAL(wp) :: trac_new(:,:,:)
TYPE(t_hydro_ocean_state), TARGET :: p_os
!3 arrays for explicit part for tracer in Adams-Bashford  stepping,
!stores information across different timelevels 
!REAL(wp) :: trac_out(:,:,:)                              !new tracer 
!
!Local variables
INTEGER, PARAMETER :: no_cell_edges = 3
!REAL(wp) :: delta_z, delta_z2
!INTEGER  :: ctr, ctr_total
INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, rl_start_c, rl_end_c
INTEGER  :: jc, jk, jb!, jkp1        !< index of edge, vert level, block
INTEGER  :: z_dolic
REAL(wp) :: z_in(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_out(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_pred(nproma,n_zlev,p_patch%nblks_c)

REAL(wp) :: z_up(nproma,n_zlev,p_patch%nblks_c,no_cell_edges)
REAL(wp) :: z_down(nproma,n_zlev,p_patch%nblks_c,no_cell_edges)
REAL(wp) :: z_max(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_min(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_excess(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_grad_excess(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_diff_excess(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_K(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: max_val, min_val!, dtime2, z_tmp
INTEGER  :: il_e(no_cell_edges), ib_e(no_cell_edges)
INTEGER  :: il_c1, ib_c1,il_c2, ib_c2, ie
INTEGER  :: stop_ctr
INTEGER            :: iter
INTEGER, PARAMETER :: i_max_iter = 20
!-----------------------------------------------------------------------
rl_start_c   = 1
rl_end_c     = min_rlcell
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

z_in       = trac_old
z_out      = trac_new
z_pred     = trac_new
z_excess   = 0.0_wp
z_K(:,:,:) = 1.0E12_wp 

! DO jk=1,n_zlev
!   write(*,*)'max/min tracer in:',jk,&
!   &maxval(trac_new(:,jk,:)),minval(trac_new(:,jk,:))
! END DO



  !Step 1: determine minimal and maximal permissible values
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)

    DO jc = i_startidx_c, i_endidx_c
      z_dolic = v_base%dolic_c(jc,jb)
      IF(z_dolic>=MIN_DOLIC)THEN
        DO jk = 1, z_dolic
          DO ie=1,no_cell_edges

            !actual edges of cell c1
            il_e(ie) = p_patch%cells%edge_idx(jc,jb,ie)
            ib_e(ie) = p_patch%cells%edge_blk(jc,jb,ie)

            !get neighbor cells of edge
            il_c1 = p_patch%edges%cell_idx(il_e(ie),ib_e(ie),1)
            ib_c1 = p_patch%edges%cell_blk(il_e(ie),ib_e(ie),1)
            il_c2 = p_patch%edges%cell_idx(il_e(ie),ib_e(ie),2)
            ib_c2 = p_patch%edges%cell_blk(il_e(ie),ib_e(ie),2)
!             IF(z_in(il_c1,jk,ib_c1)>z_in(il_c2,jk,ib_c2))THEN
!               z_up(jc,jk,jb,ie)   = z_in(il_c1,jk,ib_c1) 
!               z_down(jc,jk,jb,ie) = z_in(il_c2,jk,ib_c2)
!             ELSE
!               z_up(jc,jk,jb,ie)   = z_in(il_c2,jk,ib_c2) 
!               z_down(jc,jk,jb,ie) = z_in(il_c1,jk,ib_c1)
!             ENDIF
              IF(p_os%p_diag%ptp_vn(il_e(ie),jk,ib_e(ie))>=0.0_wp)THEN
                z_up(jc,jk,jb,ie) = trac_old(il_c1,jk,ib_c1)
              ELSEIF(p_os%p_diag%ptp_vn(il_e(ie),jk,ib_e(ie))<0.0_wp)THEN
                z_up(jc,jk,jb,ie) = trac_old(il_c2,jk,ib_c2)
              ENDIF
         END DO
         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
           z_max(jc,jk,jb) = maxval(z_up(jc,jk,jb,1:no_cell_edges))
           z_min(jc,jk,jb) = minval(z_up(jc,jk,jb,1:no_cell_edges))
         ELSE
           z_max(jc,jk,jb) = 0.0_wp
           z_min(jc,jk,jb) = 0.0_wp
         ENDIF        
!         z_max(jc,jk,jb) = maxval(z_up(jc,jk,jb,1:no_cell_edges))
!         z_min(jc,jk,jb) = minval(z_down(jc,jk,jb,1:no_cell_edges))

        END DO!level-loop
      ENDIF!(z_dolic>0)
    END DO!idx-loop
  END DO !blk-loop

! DO jk=1,n_zlev
! write(*,*)'admissible bound',jk,&
! &maxval(z_max(:,jk,:)), minval(z_max(:,jk,:)),&
! &maxval(z_min(:,jk,:)), minval(z_min(:,jk,:))
! END DO


ITERATION_LOOP: Do iter = 1, i_max_iter
!write(*,*)'iteration----------',iter
!write(1230,*)'iteration----------',iter
  !Step 1: determine excess field
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)

    DO jc = i_startidx_c, i_endidx_c
      z_dolic = v_base%dolic_c(jc,jb)
      IF(z_dolic>MIN_DOLIC)THEN
        DO jk = 1, z_dolic
         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
           z_excess(jc,jk,jb) = max(z_pred(jc,jk,jb)-z_max(jc,jk,jb),0.0_wp)&
                              &+min(z_pred(jc,jk,jb)-z_min(jc,jk,jb),0.0_wp)
          ELSE
            z_excess(jc,jk,jb) = 0.0_wp
          ENDIF 
 
!  IF(z_excess(jc,jk,jb)/=0.0)THEN
!  write(1230,*)'excess field',jc,jk,jb,z_excess(jc,jk,jb),&
! &z_pred(jc,jk,jb), z_max(jc,jk,jb),z_up(jc,jk,jb,1:no_cell_edges)
!  ENDIF
        END DO!level-loop
      ENDIF!(z_dolic>0)
    END DO!idx-loop
  END DO !blk-loop

! DO jk=1,n_zlev
! write(*,*)'max-min excess',jk,&
! &maxval(z_excess(:,jk,:)), minval(z_excess(:,jk,:))
! END DO


   !Step 3: Calculate diffusion of excess field
   CALL grad_fd_norm_oce( z_excess, p_patch, z_grad_excess)
   z_grad_excess = z_K*z_grad_excess
   CALL div_oce( z_grad_excess, p_patch, z_diff_excess)

! DO jk=1,n_zlev
! write(*,*)'max-min diffusion',jk,&
! &maxval(z_diff_excess(:,jk,:)), minval(z_diff_excess(:,jk,:))
! END DO

  !Step 4
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)

    DO jc = i_startidx_c, i_endidx_c
      z_dolic = v_base%dolic_c(jc,jb)
      IF(z_dolic>=MIN_DOLIC)THEN
        DO jk = 1, z_dolic
          IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
            z_out(jc,jk,jb) = z_pred(jc,jk,jb) - z_excess(jc,jk,jb)+ z_diff_excess(jc,jk,jb)

!  IF(z_excess(jc,jk,jb)/=0.0)THEN
!  write(1230,*)'alg', z_out(jc,jk,jb),z_pred(jc,jk,jb),&
! &z_excess(jc,jk,jb), z_diff_excess(jc,jk,jb), &
! &z_max(jc,jk,jb),z_up(jc,jk,jb,1:no_cell_edges)
!  ENDIF

          ELSE
            z_out(jc,jk,jb) = 0.0_wp
          ENDIF
!           IF(z_diff_excess(jc,1,jb)/=0.0_wp)THEN
!             write(123,*)'correction',jc,jb,z_in(jc,1,jb),&
!             & z_out(jc,1,jb), z_diff_excess(jc,1,jb)
!           ENDIF
        END DO!level-loop
      ENDIF!(z_dolic>0)
    END DO!idx-loop
  END DO !blk-loop

  !Step 4: Stop criterion
  stop_ctr = 0
  DO jk=1,n_zlev

!      write(*,*)'actual state',jk,&
!       & maxval(z_out(:,jk,:)),minval(z_out(:,jk,:))
!      write(*,*)' pred',jk,&
!       & maxval(z_pred(:,jk,:)),minval(z_pred(:,jk,:))

    IF(   maxval(z_excess(:,jk,:))<1.0E-12_wp&
    &.AND.minval(z_excess(:,jk,:))<1.0E-12_wp)THEN

      stop_ctr = stop_ctr +1

    ELSE
    ENDIF
  END DO

  IF(stop_ctr==n_zlev)THEN
    DO jb = i_startblk_c, i_endblk_c
      CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
      &                   rl_start_c, rl_end_c)

      DO jc = i_startidx_c, i_endidx_c
        z_dolic = v_base%dolic_c(jc,jb)
        IF(z_dolic>=MIN_DOLIC)THEN
          DO jk = 1, z_dolic 
           IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
           trac_new(jc,jk,jb) = z_out(jc,jk,jb)
           ENDIF
!            write(1230,*)'before-after',jc,jk,jb,&
!           &trac_new(jc,jk,jb),trac_old(jc,jk,jb),z_excess(jc,jk,jb)
          END DO
        ENDIF
       ENDDO
     END DO

      exit ITERATION_LOOP
   ELSE
       z_pred(:,:,:) = z_out(:,:,:)
  ENDIF


END DO ITERATION_LOOP

!   DO jb = i_startblk_c, i_endblk_c
!     CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
!     &                   rl_start_c, rl_end_c)
!     DO jc = i_startidx_c, i_endidx_c
! write(123,*)'trac old new',trac_old(jc,1,jb),z_old_new(jc,1,jb),&
! &trac_new(jc,1,jb), z_out(jc,1,jb)
! 
!    END DO
!  END DO
END SUBROUTINE elad
!-------------------------------------------------------------------------------
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
          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
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
          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
            pupflux_e(je,jk,jb) =  0.5_wp*pvn_e(je,jk,jb)        &
              &        *( pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)) &
              &          +pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
  END SUBROUTINE central_hflux_oce
!-------------------------------------------------------------------------------
 !>
  !! First order scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using Miura'sapproach with first
  !! order reconstruction. Miuras approach is described in H. Miura, Monthly
  !! Weather Review 135, 2007, 4038-4044. The original method has been formulated for
  !! hexagons and is transfered to triangles and uses the Mimetic capabilities.
  !!
  !! @par Revision History
  !! Developed by P. Korn, (2011).
  !!
  SUBROUTINE mimetic_miura_hflux_oce( ppatch, pvar_c, pvn_e, pvn_cc,pflux_e, opt_slev, opt_elev )
    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch      !< patch on which computation is performed
    REAL(wp), INTENT(INOUT)  :: pvar_c(:,:,:)      !< advected cell centered variable
    REAL(wp), INTENT(INOUT)  :: pvn_e(:,:,:)       !< normal velocity on edges
    TYPE(t_cartesian_coordinates) :: pvn_cc(:,:,:)
    !EAL(wp), INTENT(INOUT), OPTIONAL :: ph_e (:,:)         !< surface elevation on edges
    REAL(wp), INTENT(INOUT)  :: pflux_e(:,:,:)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_slev    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_elev    ! optional vertical end level
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER  :: slev, elev
    INTEGER  :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, rl_start, rl_end
    INTEGER  :: ie,je, jk, jb         !< index of edge, vert level, block 
    INTEGER  :: il_v1, il_v2, ib_v1, ib_v2, il_e, ib_e 
    INTEGER  :: il_c, ib_c
    TYPE(t_cartesian_coordinates) :: u_v1_cc!(nproma,n_zlev,p_patch%nblks_v)
    TYPE(t_cartesian_coordinates) :: u_v2_cc!(nproma,n_zlev,p_patch%nblks_v)
    TYPE(t_cartesian_coordinates) :: u_mean_cc(nproma,n_zlev,ppatch%nblks_e)
    TYPE(t_cartesian_coordinates) :: edge_cc  !(nproma,n_zlev,p_patch%nblks_e)
    TYPE(t_cartesian_coordinates) :: C_e(nproma,n_zlev,ppatch%nblks_e)
    TYPE(t_cartesian_coordinates) :: C_c
    REAL(wp) :: z_thick, z_weight
    REAL(wp) :: z_gradC(nproma,n_zlev,ppatch%nblks_e)
    TYPE(t_cartesian_coordinates) :: z_gradC_cc(nproma,n_zlev,ppatch%nblks_c)
    REAL(wp) :: z_tmpC(nproma,n_zlev,ppatch%nblks_c)  
    REAL(wp)  :: pupflux_e(nproma,n_zlev,ppatch%nblks_e)
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

    i_startblk_e = ppatch%edges%start_blk(rl_start,1)
    i_endblk_e   = ppatch%edges%end_blk(rl_end,1)

!------------------------------------------------------------------

!-------------upwind comparison
    iilc => ppatch%edges%cell_idx
    iibc => ppatch%edges%cell_blk
    ! loop through all patch edges (and blocks)
    DO jb = i_startblk_e, i_endblk_e
      CALL get_indices_e(ppatch, jb, i_startblk_e, i_endblk_e,   &
        &                i_startidx_e, i_endidx_e, rl_start,rl_end)

      DO jk = slev, elev
        DO je = i_startidx_e, i_endidx_e          !
          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
             pupflux_e(je,jk,jb) =  &
             &  laxfr_upflux( pvn_e(je,jk,jb), pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), &
             &                                 pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )
!write(*,*)'upwind miura',je,jk,jb,pupflux_e(je,jk,jb),pflux_e(je,jk,jb), pvn_e(je,jk,jb)!,&
!&pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2))
          ELSE
            pupflux_e(je,jk,jb) = 0.0_wp
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
!---------------------------------------------------------------

    !Step 1: velocity vector (not just normal component) at cell edges.
    !obtained by arithmetic average of two vertex velocity vectors
    LEVEL_LOOP: DO jk = slev, elev
      EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e

        CALL get_indices_e(ppatch, jb, i_startblk_e, i_endblk_e,&
                         & i_startidx_e, i_endidx_e, rl_start, rl_end)

        EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e

          IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN

            !Get indices of two adjacent vertices
            il_v1 = ppatch%edges%vertex_idx(je,jb,1)
            ib_v1 = ppatch%edges%vertex_blk(je,jb,1)
            il_v2 = ppatch%edges%vertex_idx(je,jb,2)
            ib_v2 = ppatch%edges%vertex_blk(je,jb,2)

            u_v1_cc%x  = 0.0_wp
            z_weight   = 0.0_wp
            DO ie=1, ppatch%verts%num_edges(il_v1,ib_v1) !no_vert_edges

              il_e = ppatch%verts%edge_idx(il_v1,ib_v1,ie)
              ib_e = ppatch%verts%edge_blk(il_v1,ib_v1,ie)

              z_thick = v_base%del_zlev_m(jk)
              IF ( iswm_oce == 1 ) THEN
                z_thick    = 1.0_wp!v_base%del_zlev_m(jk) +h_e(il_e,ib_e)
              ELSEIF( iswm_oce /= 1 ) THEN 
                !IF (jk == 1 )THEN
                z_thick = v_base%del_zlev_m(jk) !+ h_e(il_e,ib_e) 
                !ENDIF
              ENDIF 

              IF ( v_base%lsm_oce_e(il_e,jk,ib_e) < sea_boundary ) THEN
                z_weight = z_weight + v_base%variable_dual_vol_norm(il_v1,ib_v1,ie)*z_thick

                u_v1_cc%x = u_v1_cc%x +                      &
                &           v_base%edge2vert_coeff_cc(il_v1,ib_v1,ie)%x * &
                &           pvn_e(il_e,jk,ib_e)*z_thick
                 
              ENDIF
            END DO
            IF(z_weight/=0.0_wp)THEN
              u_v1_cc%x = u_v1_cc%x/(z_weight)!*z_h_approx)
            ENDIF

            u_v2_cc%x   = 0.0_wp
            z_weight    = 0.0_wp
            DO ie=1, ppatch%verts%num_edges(il_v2,ib_v2) !no_vert_edges

              il_e = ppatch%verts%edge_idx(il_v2,ib_v2,ie)
              ib_e = ppatch%verts%edge_blk(il_v2,ib_v2,ie)

              z_thick = v_base%del_zlev_m(jk)
              IF ( iswm_oce == 1 ) THEN
                z_thick    = 1.0_wp!v_base%del_zlev_m(jk) + h_e(il_e,ib_e)
              ELSEIF( iswm_oce /= 1 ) THEN 
               !IF (jk == 1 )THEN
                z_thick = v_base%del_zlev_m(jk) !+ h_e(il_e,ib_e) 
               !ENDIF
              ENDIF 
              IF ( v_base%lsm_oce_e(il_e,jk,ib_e)  < sea_boundary ) THEN
                z_weight = z_weight + v_base%variable_dual_vol_norm(il_v2,ib_v2,ie)*z_thick

                u_v2_cc%x = u_v2_cc%x +                      &
                &           v_base%edge2vert_coeff_cc(il_v2,ib_v2,ie)%x * &
                &           pvn_e(il_e,jk,ib_e)*z_thick
              ENDIF
            END DO
            IF(z_weight/=0.0_wp)THEN
              u_v2_cc%x = u_v2_cc%x/(z_weight)!*z_h_approx)
            ENDIF

           u_mean_cc(je,jk,jb)%x=0.5_wp*( u_v1_cc%x+ u_v2_cc%x)
         ELSE
           u_mean_cc(je,jk,jb)%x= 0.0_wp
       ENDIF
      END DO EDGE_IDX_LOOP
    END DO EDGE_BLK_LOOP
  END DO LEVEL_LOOP

  !Step 2: position of upwind point (this corresponds to the calculation of the
  !point C_i in Miura (cf eq. (9))
  DO jk = slev, elev
    DO jb = i_startblk_e, i_endblk_e

      CALL get_indices_e(ppatch, jb, i_startblk_e, i_endblk_e,&
                       & i_startidx_e, i_endidx_e, rl_start, rl_end)

      DO je =  i_startidx_e, i_endidx_e

        IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
          edge_cc = gc2cc(ppatch%edges%center(je,jb))
       
          C_e(je,jk,jb)%x   = edge_cc%x  -0.5_wp*dtime*u_mean_cc(je,jk,jb)%x

        ENDIF
      END DO 
    END DO
  END DO

    !Step3: Local linear subgridcell distribution of tracer is calculated
    !3a: calculate tracer gradient
    !3b: map tracer gradient from edge to cell. Result is gradient vector at cell centers
    !3c: project gradient vector at cell center in direction of vector that points from
    !    cell center to upwind point C_i from step 2
    !3a:
    CALL grad_fd_norm_oce( pvar_c, ppatch, z_gradC)
    !3b:
    CALL map_edges2cell( ppatch, z_gradC, z_gradC_cc)

    DO jk = slev, elev
      DO jb = i_startblk_e, i_endblk_e

        CALL get_indices_e(ppatch, jb, i_startblk_e, i_endblk_e,&
                         & i_startidx_e, i_endidx_e, rl_start, rl_end)

        DO je =  i_startidx_e, i_endidx_e

         !3c: if vn > 0 (vn < 0), the upwind cell is cell 1 (cell 2)
         IF( pvn_e(je,jk,jb) >0.0_wp)THEN 
           !transform coordinate of cell 1 to cartesian
           il_c = ppatch%edges%cell_idx(je,jb,1)
           ib_c = ppatch%edges%cell_blk(je,jb,1)

           C_c   = gc2cc(ppatch%cells%center(il_c,ib_c)) 
!             pflux_e(je,jk,jb)=&
!             &dot_product(ppatch%edges%primal_cart_normal(je,jb)%x,u_mean_cc(je,jk,jb)%x)&
!             &*(pvar_c(il_c,jk,ib_c)-dot_product(C_e(je,jk,jb)%x-C_c%x,z_gradC_cc(il_c,jk,ib_c)%x))
           pflux_e(je,jk,jb)=pvn_e(je,jk,jb)&
   &*(pvar_c(il_c,jk,ib_c)-dot_product(C_e(je,jk,jb)%x-C_c%x,z_gradC_cc(il_c,jk,ib_c)%x))

         ELSEIF(pvn_e(je,jk,jb) <=0.0_wp)THEN
           !transform coordinate of cell 1 to cartesian
           il_c = ppatch%edges%cell_idx(je,jb,2)
           ib_c = ppatch%edges%cell_blk(je,jb,2)

           C_c   = gc2cc(ppatch%cells%center(il_c,ib_c)) 
!              pflux_e(je,jk,jb)=&
!              &dot_product(ppatch%edges%primal_cart_normal(je,jb)%x,u_mean_cc(je,jk,jb)%x)&
!             &*(pvar_c(il_c,jk,ib_c)-dot_product(C_e(je,jk,jb)%x-C_c%x,z_gradC_cc(il_c,jk,ib_c)%x))
           pflux_e(je,jk,jb)=pvn_e(je,jk,jb)&
   &*(pvar_c(il_c,jk,ib_c)-dot_product(C_e(je,jk,jb)%x-C_c%x,z_gradC_cc(il_c,jk,ib_c)%x))
         ENDIF 

!          C_c   = gc2cc(ppatch%cells%center(il_c,ib_c)) 
!            pflux_e(je,jk,jb)=pvn_e(je,jk,jb)&
!            &*(pvar_c(il_c,jk,ib_c)-dot_product(C_e%x-C_c%x,z_gradC_cc(il_c,jk,ib_c)%x))

!           pflux_e(je,jk,jb)=pvn_e(je,jk,jb)&
!           &*dot_product(C_e%x-C_c%x,z_gradC_cc(il_c,jk,ib_c)%x)
! IF(pvar_c(il_c,jk,ib_c)/=0.0_wp)THEN
! write(*,*)'upwind miura',je,jk,jb,pupflux_e(je,jk,jb),pflux_e(je,jk,jb), pvn_e(je,jk,jb),&
! &dot_product(C_e%x-C_c%x,z_gradC_cc(il_c,jk,ib_c)%x)
! !&pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2))
! ENDIF
        END DO 
      END DO
    END DO

  END SUBROUTINE mimetic_miura_hflux_oce
  !-------------------------------------------------------------------------
  !>
  !! Flux limiter for horizontal advection
  !!
  !! Zalesak Flux-Limiter (Flux corrected transport)
  !! The corrected flux is a weighted average of the low order flux and the
  !! given high order flux. The high order flux is used to the greatest extent
  !! possible without introducing overshoots and undershoots.
  !! Note: This limiter is positive definite and almost monotone (but not strictly).
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2010-03-10)
  !! Modification by Daniel Reinert, DWD (2010-03-25)
  !! - adapted for MPI parallelization
  !!
  SUBROUTINE hflx_limiter_oce_mo( ptr_patch, p_cc, p_mass_flx_e, &
    &                         p_mflx_tracer_h, p_thick_old, p_thick_new, &
    &                         opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &   !< patch on which computation is performed
      &  ptr_patch

    REAL(wp), INTENT(IN) ::     &    !< advected cell centered variable
      &  p_cc(:,:,:)                 !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(in) ::     &    !< contravariant horizontal mass flux
      &  p_mass_flx_e(:,:,:)         !< (provided by dynamical core)
                                     !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(INOUT) ::  &    !< calculated horizontal tracer mass flux
      &  p_mflx_tracer_h(:,:,:)      !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(INOUT) ::  &    !< old surface elevation
      &  p_thick_old(:,:,:)              !< dim: (nproma,n_zlev,nblks_c)

    REAL(wp), INTENT(INOUT) ::  &    !< new surface elevation
      &  p_thick_new(:,:,:)              !< dim: (nproma,n_zlev,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev


!Local variables
    REAL(wp) ::                 &    !< first order tracer mass flux
      &  z_mflx_low(nproma,ptr_patch%nlev,ptr_patch%nblks_e)

    REAL(wp) ::                 &    !< antidiffusive tracer mass flux (F_H - F_L)
      &  z_anti(nproma,ptr_patch%nlev,ptr_patch%nblks_e)

    REAL(wp) ::                 &    !< antidiffusive tracer mass flux (F_H - F_L)
      &  z_mflx_anti(nproma,ptr_patch%nlev,ptr_patch%nblks_c,3) !< (units kg/kg)

    REAL(wp) ::                 &    !< flux divergence at cell center
      &  z_fluxdiv_c(nproma,ptr_patch%nlev,ptr_patch%nblks_c)

    REAL(wp) ::                 &    !< new tracer field after hor. transport,
      &  z_tracer_new_low(nproma,ptr_patch%nlev,ptr_patch%nblks_c) 
                                     !< if the low order fluxes are used

    REAL(wp) ::                 &    !< local maximum of current tracer value and low
      &  z_tracer_max(nproma,ptr_patch%nlev,ptr_patch%nblks_c) !< order update

    REAL(wp) ::                 &    !< local minimum of current tracer value and low
      &  z_tracer_min(nproma,ptr_patch%nlev,ptr_patch%nblks_c) !< order update

    REAL(wp) ::                 &    !< fraction which must multiply all in/out fluxes 
      &  r_p(nproma,ptr_patch%nlev,ptr_patch%nblks_c),&   !< of cell jc to guarantee
      &  r_m(nproma,ptr_patch%nlev,ptr_patch%nblks_c)     !< no overshoot/undershoot

    REAL(wp) :: r_frac !< computed minimum fraction which must multiply
                       !< the flux at the edge

    REAL(wp) :: z_min(nproma,ptr_patch%nlev), & !< minimum/maximum value in cell and neighboring cells
      &         z_max(nproma,ptr_patch%nlev) 
    REAL(wp) :: z_signum             !< sign of antidiffusive velocity
    REAL(wp) :: p_p, p_m             !< sum of antidiffusive fluxes into and out of cell jc

    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices of two
      &  iilc, iibc                          !< neighbor cells (array)
    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices of three
      &  iilnc, iibnc                        !< neighbor cells (array)
    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices (array)
      &  iidx, iblk                          !< of edges

    INTEGER  :: slev, elev             !< vertical start and end level
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend
    INTEGER  :: i_rlstart_e, i_rlend_e, i_rlstart_c, i_rlend_c
    INTEGER  :: je, jk, jb, jc         !< index of edge, vert level, block, cell

!     REAL(wp), POINTER ::  &
!     &  ptr_delp_mc_now(:,:,:) => NULL() !< pointer to old layer thickness
!                                         !< at cell center
!     REAL(wp), POINTER ::  &
!     &  ptr_delp_mc_new(:,:,:) => NULL() !< pointer to new layer thickness
!                                         !< at cell center
  !-------------------------------------------------------------------------
    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = ptr_patch%nlev
    END IF

    i_rlstart = 1
    i_rlend = min_rledge

    ! Set pointers to index-arrays
    ! line and block indices of two neighboring cells
    iilc => ptr_patch%edges%cell_idx
    iibc => ptr_patch%edges%cell_blk
    ! line and block indices of edges as seen from cells
    iidx => ptr_patch%cells%edge_idx
    iblk => ptr_patch%cells%edge_blk
    ! pointers to line and block indices of three neighbor cells
    iilnc => ptr_patch%cells%neighbor_idx
    iibnc => ptr_patch%cells%neighbor_blk
    !
    ! 1. Calculate low (first) order fluxes using the standard upwind scheme and the
    !    antidiffusive fluxes
    !    (not allowed to call upwind_hflux_up directly, due to circular dependency)
    i_rlstart_e  = 1
    i_rlend_e    = min_rledge
    i_startblk   = ptr_patch%edges%start_blk(i_rlstart_e,1)
    i_endblk     = ptr_patch%edges%end_blk(i_rlend_e,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart_e, i_rlend_e)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=5
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif
          z_mflx_low(je,jk,jb) =  &
            &  laxfr_upflux( p_mass_flx_e(je,jk,jb), p_cc(iilc(je,jb,1),jk,iibc(je,jb,1)), &
            &                                        p_cc(iilc(je,jb,2),jk,iibc(je,jb,2)) )

          ! calculate antidiffusive flux for each edge
          z_anti(je,jk,jb)     = p_mflx_tracer_h(je,jk,jb) - z_mflx_low(je,jk,jb)

        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
!$OMP END DO


    i_rlstart_c  = 1
    i_rlend_c    = min_rlcell
    i_startblk   = ptr_patch%cells%start_blk(i_rlstart_c,1)
    i_endblk     = ptr_patch%cells%end_blk(i_rlend_c,1)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk,        &
                         i_startidx, i_endidx, i_rlstart_c, i_rlend_c)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=4
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif
          !
          ! 2. Define "antidiffusive" fluxes A(jc,jk,jb,je) for each cell. It is the difference
          !    between the high order fluxes (given by the FFSL-scheme) and the low order
          !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
          !    - positive for outgoing fluxes
          !    - negative for incoming fluxes
          !    this sign convention is related to the definition of the divergence operator.

          z_mflx_anti(jc,jk,jb,1) =                                                  &
            &     dtime * v_base%geofac_div(jc,1,jb) / p_thick_new(jc,jk,jb)  &
            &   * z_anti(iidx(jc,jb,1),jk,iblk(jc,jb,1))

          z_mflx_anti(jc,jk,jb,2) =                                                  &
            &     dtime *  v_base%geofac_div(jc,2,jb) / p_thick_new(jc,jk,jb)  &
            &   * z_anti(iidx(jc,jb,2),jk,iblk(jc,jb,2))

          z_mflx_anti(jc,jk,jb,3) =                                                  &
            &     dtime * v_base%geofac_div(jc,3,jb) / p_thick_new(jc,jk,jb)  &
            &   * z_anti(iidx(jc,jb,3),jk,iblk(jc,jb,3))

          !  compute also divergence of low order fluxes
          z_fluxdiv_c(jc,jk,jb) =  &
            & z_mflx_low(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * v_base%geofac_div(jc,1,jb) + &
            & z_mflx_low(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * v_base%geofac_div(jc,2,jb) + &
            & z_mflx_low(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * v_base%geofac_div(jc,3,jb)

        ENDDO
      ENDDO
      !
      ! 3. Compute the updated low order solution z_tracer_new_low
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          z_tracer_new_low(jc,jk,jb) =                           &
            &      ( p_cc(jc,jk,jb) * p_thick_old(jc,jk,jb)  &
            &      - dtime * z_fluxdiv_c(jc,jk,jb) )           &
            &      / p_thick_new(jc,jk,jb)

          ! precalculate local maximum of current tracer value and low order
          ! updated value
          z_tracer_max(jc,jk,jb) = MAX(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))

          ! precalculate local minimum of current tracer value and low order
          ! updated value
          z_tracer_min(jc,jk,jb) = MIN(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO

    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.
    i_rlstart_c  = 1
    i_rlend_c    = min_rlcell
    i_startblk   = ptr_patch%cells%start_blk(i_rlstart_c,1)
    i_endblk     = ptr_patch%cells%end_blk(i_rlend_c,1)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_max,p_p,z_min,p_m)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart_c, i_rlend_c)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=2
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif
          ! max value of cell and its neighbors
          ! also look back to previous time step
          z_max(jc,jk) = MAX( z_tracer_max(jc,jk,jb),                          &
            &                 z_tracer_max(iilnc(jc,jb,1),jk,iibnc(jc,jb,1)),  &
            &                 z_tracer_max(iilnc(jc,jb,2),jk,iibnc(jc,jb,2)),  &
            &                 z_tracer_max(iilnc(jc,jb,3),jk,iibnc(jc,jb,3)) )
          ! min value of cell and its neighbors
          ! also look back to previous time step
          z_min(jc,jk) = MIN( z_tracer_min(jc,jk,jb),                          &
            &                 z_tracer_min(iilnc(jc,jb,1),jk,iibnc(jc,jb,1)),  &
            &                 z_tracer_min(iilnc(jc,jb,2),jk,iibnc(jc,jb,2)),  &
            &                 z_tracer_min(iilnc(jc,jb,3),jk,iibnc(jc,jb,3)) )
        ENDDO
      ENDDO

      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          ! Sum of all incoming antidiffusive fluxes into cell jc
          p_p =  -1._wp * (MIN(0._wp,z_mflx_anti(jc,jk,jb,1))   &
                         + MIN(0._wp,z_mflx_anti(jc,jk,jb,2))   &
                         + MIN(0._wp,z_mflx_anti(jc,jk,jb,3)) )
          ! Sum of all outgoing antidiffusive fluxes out of cell jc
          p_m =  MAX(0._wp,z_mflx_anti(jc,jk,jb,1))  &
            &  + MAX(0._wp,z_mflx_anti(jc,jk,jb,2))  &
            &  + MAX(0._wp,z_mflx_anti(jc,jk,jb,3))

          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of q
          r_m(jc,jk,jb) = (z_tracer_new_low(jc,jk,jb) - z_min(jc,jk))/(p_m + dbl_eps)
          ! fraction which must multiply all fluxes into cell jc to guarantee no
          ! overshoot
          ! Nominator: maximum allowable increase of q
          r_p(jc,jk,jb) = (z_max(jc,jk) - z_tracer_new_low(jc,jk,jb))/(p_p + dbl_eps)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ! Synchronize r_m and r_p
    !CALL sync_patch_array_mult(SYNC_C1, ptr_patch, 2, r_m, r_p)

    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !    At the end, compute new, limited fluxes which are then passed to the main
    !    program. Note that p_mflx_tracer_h now denotes the LIMITED flux.
    i_startblk = ptr_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = ptr_patch%edges%end_blk(i_rlend,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,r_frac,z_signum)
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk,    &
                         i_startidx, i_endidx, i_rlstart, i_rlend)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=3
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif
          z_signum = SIGN(1._wp,z_anti(je,jk,jb))

          ! This does the same as an IF (z_signum > 0) THEN ... ELSE ... ENDIF,
          ! but is computationally more efficient
          r_frac = 0.5_wp*( (1._wp+z_signum)*                &
             &     MIN(r_m(iilc(je,jb,1),jk,iibc(je,jb,1)),  &
             &         r_p(iilc(je,jb,2),jk,iibc(je,jb,2)))  &
             &     +  (1._wp-z_signum)*                      &
             &     MIN(r_m(iilc(je,jb,2),jk,iibc(je,jb,2)),  &
             &         r_p(iilc(je,jb,1),jk,iibc(je,jb,1)))  )

          ! Limited flux
          p_mflx_tracer_h(je,jk,jb) = z_mflx_low(je,jk,jb)               &
            &                       + MIN(1._wp,r_frac) * z_anti(je,jk,jb)

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE hflx_limiter_oce_mo
  !-------------------------------------------------------------------------
  !>
  !! Monotonous slope limiter for MIURA advection scheme
  !!
  !! Gradient limitation using monotonous Barth-Jesperson limiter
  !! to avoid over/undershoots in the advected field. The limiter is modified
  !! in the sense, that the points where the reconstruction is evaluated
  !! (the Gauss points) have been made time dependent. This wass necessary in
  !! order to achieve positive definiteness.
  !!
  !! @par Literature
  !! - T. J. Barth and D. C. Jespersen. The design and application of upwind
  !!   schemes on unstructured meshes. In 27th Aerospace Sciences Meeting, Jan.
  !!   1989. AIAA paper 89-0366.
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2010-02-26)
  !!
  SUBROUTINE h_miura_slimiter_oce_mo( p_patch, p_cc, ptr_cc, p_distv_gausspoint, p_grad, &
    &                              opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &      !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::    &       !< advected cell centered variable
      &  p_cc(:,:,:)                   !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &       !< different from p_cc if conservative lsq
      &  ptr_cc(:,:,:)                 !< reconstruction is chosen
                                       !< In this case ptr_cc contains coefficient
                                       !< c0 from the lsq reconstruction
                                       !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &       !< distance vectors cell center -->
      &  p_distv_gausspoint(:,:,:,:,:) !< shifted Gauss-points
                                       !< dim: (nproma,nlev,p_patch%nblks_c,3,2)

    REAL(wp), INTENT(INOUT) :: &    !< reconstructed gradients, to be limited
      &  p_grad(:,:,:,:)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    !--- Local variables ---!
    REAL(wp) ::            &
      &  z_min(nproma),    &        !< minimum of cell centered value
      &  z_max(nproma)              !< maximum of cell centered value

    REAL(wp) :: z_lext_val_1,  &    !< linear extrapolation value of cell
      &  z_lext_val_2, z_lext_val_3 !< for edges 1, 2, 3

    REAL(wp) :: z_limit             !< value of gradient limiter

    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices (array)
      &  iilc, iibc

    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER  :: jk, jb, jc         !< index of vert level, block, cell
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend

  !-------------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = p_patch%nlev
    END IF

      i_rlstart = 1
      i_rlend = min_rlcell


    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,1)

    ! pointers to line and block indices of three neighbor cells
    iilc => p_patch%cells%neighbor_idx
    iibc => p_patch%cells%neighbor_blk


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_min,z_max,z_lext_val_1,   &
!$OMP            z_lext_val_2,z_lext_val_3,z_limit)

    ! Calculate reconstructed tracer values at Gauss-points provided and limit
    ! the gradient appropriately. Note that one needs to distinguish beween
    ! the prognostic field (p_cc) and the constant c0 (ptr_cc) from the lsq
    ! reconstruction. These fields differ in the case of a conservative least
    ! squares reconstruction.
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev

        DO jc = i_startidx, i_endidx
#endif

          ! min value of cell and its neighbors
          z_min(jc) = MIN( p_cc(jc,jk,jb), p_cc(iilc(jc,jb,1),jk,iibc(jc,jb,1)), &
            &              p_cc(iilc(jc,jb,2),jk,iibc(jc,jb,2)),                 &
            &              p_cc(iilc(jc,jb,3),jk,iibc(jc,jb,3)) )
          ! max value of cell and its neighbors
          z_max(jc) = MAX( p_cc(jc,jk,jb), p_cc(iilc(jc,jb,1),jk,iibc(jc,jb,1)), &
            &              p_cc(iilc(jc,jb,2),jk,iibc(jc,jb,2)),                 &
            &              p_cc(iilc(jc,jb,3),jk,iibc(jc,jb,3)) )


          ! value for linear extrapolation of cell centered value to Gauss-point 1
          z_lext_val_1 = ( ptr_cc(jc,jk,jb)                                      &
            &          + p_grad(jc,jk,jb,1) * p_distv_gausspoint(jc,jk,jb,1,1)   &
            &          + p_grad(jc,jk,jb,2) * p_distv_gausspoint(jc,jk,jb,1,2) ) &
            &          - p_cc(jc,jk,jb)

          ! value for linear extrapolation of cell centered value to Gauss-point 2
          z_lext_val_2 = ( ptr_cc(jc,jk,jb)                                      &
            &          + p_grad(jc,jk,jb,1) * p_distv_gausspoint(jc,jk,jb,2,1)   &
            &          + p_grad(jc,jk,jb,2) * p_distv_gausspoint(jc,jk,jb,2,2) ) &
            &          - p_cc(jc,jk,jb)

          ! value for linear extrapolation of cell centered value to Gauss-point 3
          z_lext_val_3 = ( ptr_cc(jc,jk,jb)                                      &
            &          + p_grad(jc,jk,jb,1) * p_distv_gausspoint(jc,jk,jb,3,1)   &
            &          + p_grad(jc,jk,jb,2) * p_distv_gausspoint(jc,jk,jb,3,2) ) &
            &          - p_cc(jc,jk,jb)



          z_limit = HUGE(1.0_wp)


          ! calculation of limiters
          IF ( z_lext_val_1 > 0._wp ) THEN
            z_limit = MIN( z_limit,                         &
              &  MIN( 1.0_wp, (z_max(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_1) + dbl_eps),z_lext_val_1 ) ) )
          ELSE IF ( z_lext_val_1 < 0._wp ) THEN
            z_limit = MIN( z_limit,                         &
              &  MIN( 1.0_wp, (z_min(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_1) + dbl_eps),z_lext_val_1 ) ) )
          ELSE
            z_limit = MIN( z_limit, 1.0_wp )
          END IF

          IF ( z_lext_val_2 > 0._wp ) THEN
            z_limit = MIN( z_limit,                         &
              &  MIN( 1.0_wp, (z_max(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_2) + dbl_eps),z_lext_val_2 ) ) )
          ELSE IF ( z_lext_val_2 < 0._wp ) THEN
            z_limit = MIN( z_limit,                         &
              &  MIN( 1.0_wp, (z_min(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_2) + dbl_eps),z_lext_val_2 ) ) )
          ELSE
            z_limit = MIN( z_limit, 1.0_wp )
          END IF

          IF ( z_lext_val_3 > 0._wp ) THEN
           z_limit = MIN( z_limit,                          &
              &  MIN( 1.0_wp, (z_max(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_3) + dbl_eps),z_lext_val_3 ) ) )
          ELSE IF ( z_lext_val_3 < 0._wp ) THEN
            z_limit = MIN( z_limit,                         &
              &  MIN( 1.0_wp, (z_min(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_3) + dbl_eps),z_lext_val_3 ) ) )
          ELSE
            z_limit = MIN( z_limit, 1.0_wp )
          END IF


          !
          ! limited value of gradient at cell center
          !
          p_grad(jc,jk,jb,1:2) = z_limit * p_grad(jc,jk,jb,1:2)

          END DO ! end loop over cells

        END DO

      END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE h_miura_slimiter_oce_mo




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
  SUBROUTINE upwind_vflux_oce( ppatch, pvar_c, pw_c,top_bc_t, pupflux_i )

    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch      !< patch on which computation is performed
    REAL(wp), INTENT(INOUT)  :: pvar_c(:,:,:)      !< advected cell centered variable
    REAL(wp), INTENT(INOUT)  :: pw_c(:,:,:)        !< vertical velocity on cells 
    REAL(wp), INTENT(INOUT)  :: top_bc_t(:,:)        !< top boundary condition traver!vertical velocity on cells
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

    !fluxes at first layer
    pupflux_i(:,1,:) = pvar_c(:,1,:)*pw_c(:,1,:) !0.0!
    !pupflux_i(:,1,:) = pvar_c(:,1,:)*(pw_c(:,1,:)+pw_c(:,2,:))*0.5_wp !0.0_wp
    !pupflux_i(:,1,:) = pvar_c(:,1,:)*top_bc_t  ! max(abs(pw_c(:,1,:)),abs(pw_c(:,2,:)))
    !
    DO jb =  i_startblk, i_endblk
      CALL get_indices_c(ppatch, jb, i_startblk, i_endblk,&
                   & i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
!          IF(abs(pw_c(jc,1,jb))>=abs(pw_c(jc,2,jb)))THEN
!            pupflux_i(jc,1,jb) = pvar_c(jc,1,jb)*pw_c(jc,1,jb)
!          ELSE
!            pupflux_i(jc,1,jb) = pvar_c(jc,1,jb)*pw_c(jc,2,jb)
!          ENDIF
        z_dolic = v_base%dolic_c(jc,jb)
        IF(z_dolic>=MIN_DOLIC)THEN
          DO jk = 2, z_dolic
            jkm1 = jk - 1
            ! calculate vertical tracer flux using upwind method
             pupflux_i(jc,jk,jb) =                 &
               &  laxfr_upflux_v( pw_c(jc,jk,jb),  &
               &                  pvar_c(jc,jkm1,jb), pvar_c(jc,jk,jb), zcoeff_grid )
          END DO
          ! no fluxes at bottom boundary
          pupflux_i(jc,z_dolic+1,jb) = 0.0_wp
        ENDIF
      END DO
    END DO 
    !pupflux_i(:,1,:) = 0.5_wp*(pupflux_i(:,2,:)+pvar_c(:,1,:)*pw_c(:,1,:))
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
    !fluxes at first layer   
    pupflux_i(:,1,:) = pvar_c(:,1,:)*pw_c(:,1,:)
    !pupflux_i(:,1,:) =  pvar_c(:,1,:)*(pw_c(:,1,:)+pw_c(:,2,:))*0.5_wp!pvar_c(:,1,:)*pw_c(:,1,:) !0.0_wp
   

    DO jb =  i_startblk, i_endblk
      CALL get_indices_c(ppatch, jb, i_startblk, i_endblk,&
                   & i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          z_dolic = v_base%dolic_c(jc,jb)
          IF(z_dolic>=MIN_DOLIC)THEN
            DO jk = 2, z_dolic -1
              jkm1 = jk - 1
              jkp1 = jk + 1

              w_ave(jk)=0.5_wp*(pw_c(jc,jk,jb)+pw_c(jc,jkm1,jb))&
                    &*pvar_c(jc,jk,jb)

             w_avep1(jk)=0.5_wp*(pw_c(jc,jk,jb)+pw_c(jc,jkp1,jb))&
                    &*pvar_c(jc,jkp1,jb)

            pupflux_i(jc,jk,jb) = 0.5_wp*(w_ave(jk) &! *v_base%del_zlev_m(jk) &
                                &        +w_avep1(jk))!*v_base%del_zlev_m(jkp1))&
                                !&        /(v_base%del_zlev_m(jk)+v_base%del_zlev_m(jkp1))
            END DO 
            w_ave(z_dolic)=0.5_wp*(pw_c(jc,z_dolic,jb)+pw_c(jc,z_dolic-1,jb))&
                    &*pvar_c(jc,z_dolic,jb)
            pupflux_i(jc,z_dolic,jb) = w_ave(z_dolic)*v_base%del_zlev_m(z_dolic)&
                                      &/v_base%del_zlev_i(z_dolic)
            ! no fluxes at bottom boundary
            pupflux_i(jc,z_dolic+1,jb) = 0.0_wp
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
!            pupflux_i(jc,jk,jb) = 0.5_wp*(w_ave(jk)  *v_base%del_zlev_m(jk) &
!                                &        +w_avep1(jk)*v_base%del_zlev_m(jkp1))/v_base%del_zlev_i(jk)
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

    !fluxes at first layer
    !pupflux_i(:,1,:) =  pvar_c(:,1,:)*(pw_c(:,1,:)+pw_c(:,2,:))*0.5_wp!
    pupflux_i(:,1,:) =pvar_c(:,1,:)* pw_c(:,1,:)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(ppatch, jb, i_startblk, i_endblk,&
                      & i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          z_dolic = v_base%dolic_c(jc,jb)
          IF(z_dolic>=MIN_DOLIC)THEN
            DO jk = 2, z_dolic
              ! index of top half level
              jkm1 = jk - 1

              ! calculate vertical tracer flux using upwind method
              pupflux_i(jc,jk,jb) = 0.5_wp*pw_c(jc,jk,jb) &
                &  * (pvar_c(jc,jkm1,jb)+ pvar_c(jc,jk,jb) )
! IF(jb==959.and.jc==10)THEN 
! IF(v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary)THEN
! write(*,*)'vert fluxes:',jc,jk,jb,pupflux_i(jc,jk,jb),pw_c(jc,jk,jb),&
! & pvar_c(jc,jkm1,jb), pvar_c(jc,jk,jb) 
! ENDIF
! ENDIF 
          END DO 
          ! no fluxes at bottom boundary
          pupflux_i(jc,z_dolic+1,jb) = 0.0_wp


        ENDIF
      END DO
    END DO
  END SUBROUTINE central_vflux_oce
!------------------------------------------------------------------------
! ! 
! ! 
! ! !!>
! ! ! !  SUBROUTINE advects the tracers present in the ocean model.
! ! !
! ! ! @par Revision History
! ! ! Developed  by  Peter Korn, MPI-M (2010).
! ! ! 
! ! SUBROUTINE advect_individual_tracer_ab_old(p_patch, trac_old,           &
! !                                      & p_os, G_n_c_h,G_nm1_c_h,G_nimd_c_h, &
! !                                      & G_n_c_v, G_nm1_c_v, G_nimd_c_v, &
! !                                      & bc_top_tracer, bc_bot_tracer,&
! !                                      & K_h, A_v,                    &
! !                                      & trac_new, timestep)
! ! 
! ! 
! ! TYPE(t_patch), TARGET, INTENT(in) :: p_patch
! ! REAL(wp)                          :: trac_old(:,:,:)
! ! TYPE(t_hydro_ocean_state), TARGET :: p_os
! ! REAL(wp) :: G_n_c_h   (nproma, n_zlev,   p_patch%nblks_c)  !G^n
! ! REAL(wp) :: G_nm1_c_h (nproma, n_zlev,   p_patch%nblks_c)  !G^(n-1)
! ! REAL(wp) :: G_nimd_c_h(nproma, n_zlev,   p_patch%nblks_c)  !G^(n+1/2)
! ! REAL(wp) :: G_n_c_v   (nproma, n_zlev,   p_patch%nblks_c)  !G^n
! ! REAL(wp) :: G_nm1_c_v (nproma, n_zlev,   p_patch%nblks_c)  !G^(n-1)
! ! REAL(wp) :: G_nimd_c_v(nproma, n_zlev,   p_patch%nblks_c)  !G^(n+1/2)
! ! REAL(wp) :: bc_top_tracer(nproma, p_patch%nblks_c)
! ! REAL(wp) :: bc_bot_tracer(nproma, p_patch%nblks_c)
! ! REAL(wp) :: K_h(:,:,:)                                  !horizontal mixing coeff
! ! REAL(wp) :: A_v(:,:,:)                                   !vertical mixing coeff
! ! REAL(wp) :: trac_new(:,:,:)                              !new tracer 
! ! INTEGER  :: timestep                                     ! Actual timestep (to distinghuish initial step from others)
! ! 
! ! !Local variables
! ! REAL(wp) :: delta_z!, delta_z2
! ! !INTEGER  :: ctr, ctr_total
! ! INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, rl_start_c, rl_end_c
! ! INTEGER  :: jc, jk, jb, jkp1        !< index of edge, vert level, block
! ! INTEGER  :: z_dolic
! ! REAL(wp) :: z_adv_flux_h (nproma, n_zlev,   p_patch%nblks_e)  ! horizontal advective tracer flux
! ! REAL(wp) :: z_adv_flux_v (nproma, n_zlev+1, p_patch%nblks_c)  ! vertical advective tracer flux
! ! REAL(wp) :: z_div_adv_h (nproma, n_zlev,  p_patch%nblks_c)        ! horizontal tracer divergence
! ! REAL(wp) :: z_div_adv_v (nproma, n_zlev,p_patch%nblks_c)        ! vertical tracer divergence
! ! REAL(wp) :: z_div_diff_h (nproma, n_zlev,  p_patch%nblks_c)        ! horizontal tracer divergence
! ! REAL(wp) :: z_div_diff_v (nproma, n_zlev,p_patch%nblks_c)        ! vertical tracer divergence
! ! 
! ! REAL(wp) :: z_diff_flux_h(nproma, n_zlev,  p_patch%nblks_e)   ! horizontal diffusive tracer flux
! ! REAL(wp) :: z_diff_flux_v(nproma, n_zlev+1,p_patch%nblks_c)   ! vertical diffusive tracer flux
! ! 
! ! 
! ! REAL(wp) :: z_transport_vn(nproma,n_zlev, p_patch%nblks_e)  ! horizontal transport velocity
! ! REAL(wp) :: z_transport_w(nproma,n_zlev+1,p_patch%nblks_c)  ! vertical transport velocity
! ! 
! ! REAL(wp) :: z_mass_flux_h(nproma,n_zlev, p_patch%nblks_e)
! ! REAL(wp) :: z_trac_c(nproma,n_zlev, p_patch%nblks_c)
! ! REAL(wp) :: z_div_mass_flux_h(nproma,n_zlev, p_patch%nblks_c)
! ! REAL(wp) :: z_h(nproma,n_zlev, p_patch%nblks_c)
! ! REAL(wp) :: z_h2(nproma,n_zlev, p_patch%nblks_c)
! ! REAL(wp) :: z_h_tmp(nproma,n_zlev, p_patch%nblks_c)
! ! 
! ! REAL(wp) :: z_G_n_c_v   (nproma, n_zlev, p_patch%nblks_c)
! ! REAL(wp) :: z_G_nm1_c_v (nproma, n_zlev, p_patch%nblks_c)
! ! REAL(wp) :: z_G_nimd_c_v(nproma, n_zlev, p_patch%nblks_c)
! ! REAL(wp) :: max_val, min_val, dtime2
! ! 
! ! INTEGER            :: FLUX_CALCULATION = 3 
! ! INTEGER, PARAMETER :: UPWIND = 1
! ! INTEGER, PARAMETER :: CENTRAL= 2
! ! INTEGER, PARAMETER :: MIMETIC= 3
! ! REAL(wp)           :: z_tol
! ! LOGICAL :: ldbg = .false.
! ! TYPE(t_cartesian_coordinates):: z_vn_c(nproma,n_zlev,p_patch%nblks_c)
! ! TYPE(t_cartesian_coordinates):: z_vn_c2(nproma,n_zlev,p_patch%nblks_c)
! ! ! CHARACTER(len=max_char_length), PARAMETER :: &
! ! !        & routine = ('mo_tracer_advection:advect_individual_tracer')
! ! !-------------------------------------------------------------------------------
! ! z_tol= 0.0_wp!1.0E-13
! ! 
! ! dtime2= 1.0_wp*dtime
! ! 
! ! rl_start_c   = 1
! ! rl_end_c     = min_rlcell
! ! i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
! ! i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
! ! 
! ! ! ! rl_start_e   = 1
! ! ! ! rl_end_e     = min_rledge
! ! ! ! i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
! ! ! ! i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
! ! 
! ! z_adv_flux_h  = 0.0_wp
! ! z_div_adv_h   = 0.0_wp
! ! z_div_diff_h  = 0.0_wp
! ! z_diff_flux_h = 0.0_wp
! ! z_adv_flux_v  = 0.0_wp
! ! z_div_adv_v   = 0.0_wp
! ! z_div_diff_v  = 0.0_wp
! ! z_diff_flux_v = 0.0_wp
! ! z_mass_flux_h = 0.0_wp
! ! z_h           = 0.0_wp
! ! z_h2          = 0.0_wp
! ! z_h_tmp       =0.0_wp
! ! z_transport_w =0.0_wp
! ! z_trac_c      =0.0_wp
! ! z_div_mass_flux_h=0.0_wp
! ! z_G_n_c_v        = 0.0_wp
! ! z_G_nm1_c_v      = 0.0_wp
! ! z_G_nimd_c_v     = 0.0_wp
! ! 
! ! DO jb = i_startblk_c, i_endblk_c
! !   CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
! !   &                   rl_start_c, rl_end_c)
! !   DO jk = 1, n_zlev
! !     delta_z = v_base%del_zlev_m(jk)
! !     DO jc = i_startidx_c, i_endidx_c
! !       IF (jk == 1) delta_z = v_base%del_zlev_m(jk)&
! !                          & + p_os%p_prog(nold(1))%h(jc,jb)
! !       z_h_tmp(jc,jk,jb) = delta_z 
! !       z_trac_c(jc,jk,jb)=trac_old(jc,jk,jb)*delta_z
! !     END DO
! !   END DO
! ! END DO
! ! 
! ! !Step 1) Horizontal advection and diffusion
! ! SELECT CASE(FLUX_CALCULATION)
! ! 
! ! CASE(UPWIND)
! !   z_transport_vn = ab_gam*p_os%p_prog(nnew(1))%vn + (1.0_wp-ab_gam)*p_os%p_prog(nold(1))%vn
! ! 
! !   CALL upwind_hflux_oce( p_patch,        &
! !                        & z_h_tmp,        &
! !                        & z_transport_vn, &
! !                        & z_mass_flux_h )
! !   CALL upwind_hflux_oce( p_patch,        &
! !                        & z_trac_c,       &
! !                        & z_transport_vn, &
! !                        & z_adv_flux_h )
! ! CASE(CENTRAL)
! !   z_transport_vn = ab_gam*p_os%p_prog(nnew(1))%vn + (1.0_wp-ab_gam)*p_os%p_prog(nold(1))%vn
! ! 
! !   CALL central_hflux_oce( p_patch,       &
! !                        & z_h_tmp,        &
! !                        & z_transport_vn, &
! !                        & z_mass_flux_h )
! !   CALL central_hflux_oce( p_patch,       &
! !                        & z_trac_c,       &
! !                        & z_transport_vn, &
! !                        & z_adv_flux_h )
! ! CASE(MIMETIC)
! ! 
! !   DO jb = i_startblk_c, i_endblk_c
! !     CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
! !     &                   rl_start_c, rl_end_c)
! !     DO jk = 1, n_zlev
! !       DO jc = i_startidx_c, i_endidx_c
! ! 
! !         z_vn_c(jc,jk,jb)%x(:) =0.0_wp
! !         z_vn_c2(jc,jk,jb)%x(:)=0.0_wp
! !         IF (jk == 1)THEN
! !           delta_z = p_os%p_prog(nold(1))%h(jc,jb) +v_base%del_zlev_m(jk)!+p_os%p_prog(nold(1))%h(jc,jb)
! ! 
! !           z_vn_c(jc,jk,jb)%x  = trac_old(jc,jk,jb)*delta_z*p_os%p_diag%p_vn(jc,jk,jb)%x&
! !           &*p_patch%cells%area(jc,jb)
! !           z_vn_c2(jc,jk,jb)%x = delta_z*p_os%p_diag%p_vn(jc,jk,jb)%x*p_patch%cells%area(jc,jb)
! !          ELSE
! !           z_vn_c(jc,jk,jb)%x = trac_old(jc,jk,jb)*p_os%p_diag%p_vn(jc,jk,jb)%x
! !           z_vn_c2(jc,jk,jb)%x= p_os%p_diag%p_vn(jc,jk,jb)%x
! !          ENDIF
! !       END DO
! !     END DO
! !   END DO
! ! 
! !    CALL map_cell2edges( p_patch, z_vn_c, z_adv_flux_h )
! !    CALL map_cell2edges( p_patch, z_vn_c2, z_mass_flux_h )
! ! END SELECT
! !  
! ! !Step 4) calculate horizontal divergence
! ! CALL div_oce( z_adv_flux_h, p_patch, z_div_adv_h)
! ! CALL div_oce( z_mass_flux_h, p_patch, z_div_mass_flux_h)
! ! 
! ! DO jk = 1, n_zlev
! ! write(*,*)'max/min adv-flux:',jk,&
! ! &maxval(z_adv_flux_h(:,jk,:)),&
! ! &minval(z_adv_flux_h(:,jk,:))
! ! END DO
! ! 
! ! DO jk=1,n_zlev
! !   DO jb = i_startblk_c, i_endblk_c
! !     CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
! !     &                   rl_start_c, rl_end_c)
! !     DO jc = i_startidx_c, i_endidx_c
! !       z_h(jc,jk,jb) = z_h_tmp(jc,jk,jb) - dtime2*z_div_mass_flux_h(jc,jk,jb)
! !     END DO
! !   END DO!write(*,*)'max/min h_dummy:',maxval(z_h(:,jk,:)),minval(z_h(:,jk,:))
! ! END DO
! ! 
! ! !The diffusion part
! ! !Calculate horizontal diffusive flux
! ! CALL tracer_diffusion_horz( p_patch, &
! !  &                          trac_old,&
! !  &                          p_os,    &
! !  &                          K_h,     & 
! !  &                          z_diff_flux_h)
! ! 
! ! 
! ! CALL div_oce( z_diff_flux_h, p_patch, z_div_diff_h)
! ! 
! ! 
! ! DO jb = i_startblk_c, i_endblk_c
! !   CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
! !                    & i_startidx_c, i_endidx_c,&
! !                    & rl_start_c, rl_end_c)
! !   DO jk = 1, n_zlev
! !     DO jc = i_startidx_c, i_endidx_c
! !       G_n_c_h(jc,jk,jb) = z_div_diff_h(jc,jk,jb)
! !     END DO
! !   END DO
! ! END DO
! ! 
! ! IF(is_initial_timestep(timestep))THEN
! !     G_nimd_c_h(:,:,:) = G_n_c_h(:,:,:)
! !   ELSE
! !     G_nimd_c_h(:,:,:) = (1.5_wp+AB_const)* G_n_c_h(:,:,:)   &
! !         &             - (0.5_wp+AB_const)*G_nm1_c_h(:,:,:)
! ! ENDIF
! ! 
! ! IF(ldbg)THEN
! !   DO jk = 1, n_zlev
! !     write(*,*)'max/min G_T:',jk,&
! !     &maxval(G_n_c_h(:,jk,:)), minval(G_n_c_h(:,jk,:))
! !   END DO
! ! 
! !   DO jk = 1, n_zlev
! !     write(*,*)'max/min G_T+1/2:',jk,maxval(G_nimd_c_h(:,jk,:)),&
! !                                    &minval(G_nimd_c_h(:,jk,:))
! !   END DO
! !   DO jk = 1, n_zlev
! !     write(*,*)'max/min adv flux & div h:',jk,maxval(z_adv_flux_h(:,jk,:)),&
! !                                            & minval(z_adv_flux_h(:,jk,:)),&
! !                                            & maxval(z_div_adv_h(:,jk,:)),&
! !                                            & minval(z_div_adv_h(:,jk,:))
! !   END DO
! !   DO jk = 1, n_zlev
! !     write(*,*)'max/min mass flux & div h:',jk, maxval(z_mass_flux_h(:,jk,:)),&
! !                                               & minval(z_mass_flux_h(:,jk,:)),&
! !                                              & maxval(z_div_mass_flux_h(:,jk,:)),&
! !                                              & minval(z_div_mass_flux_h(:,jk,:))
! !   END DO
! ! ENDIF
! ! 
! ! !Final step: calculate new tracer values
! ! DO jk = 1, n_zlev
! !   CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
! !                    & i_startidx_c, i_endidx_c,&
! !                    & rl_start_c, rl_end_c)
! !     DO jb = i_startblk_c, i_endblk_c
! !       DO jc = i_startidx_c, i_endidx_c
! !         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! ! 
! !           IF(jk==1)THEN
! !             delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+v_base%del_zlev_m(top)
! !             !delta_z2  = p_os%p_prog(nnew(1))%h(jc,jb)+v_base%del_zlev_m(top)
! !           ELSE
! !             delta_z  = v_base%del_zlev_m(jk)
! !            !delta_z2  = v_base%del_zlev_m(jk) 
! !           ENDIF
! !           trac_new(jc,jk,jb) = (z_trac_c(jc,jk,jb)             &
! !                              & -dtime2*(z_div_adv_h(jc,jk,jb)  &
! !                              &            -G_nimd_c_h(jc,jk,jb)))&
! !                              &/(z_h(jc,jk,jb)+z_tol)
! !         ENDIF
! !       END DO
! !     END DO
! ! END DO
! ! 
! !   G_nm1_c_h(:,:,:)  = G_n_c_h(:,:,:)
! !   G_nimd_c_h(:,:,:) = 0.0_wp
! !   G_n_c_h(:,:,:)    = 0.0_wp
! ! 
! !   DO jk = 1, n_zlev
! !   max_val=0.0_wp
! !   min_val=100.0_wp
! !   DO jb = i_startblk_c, i_endblk_c
! !     CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
! !     &                   rl_start_c, rl_end_c)
! !     DO jc = i_startidx_c, i_endidx_c
! !       IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! !         IF(trac_new(jc,jk,jb)>max_val)&
! !           &max_val=trac_new(jc,jk,jb)
! !         IF(trac_new(jc,jk,jb)<min_val)&
! !           &min_val=trac_new(jc,jk,jb)
! ! 
! !         IF(trac_new(jc,jk,jb)==0.0_wp)&
! !         & write(*,*)'tracer =0:',jc,jk,jb 
! !       ENDIF
! !      END DO
! !    END DO
! !    WRITE(*,*)'MAX/MIN tracer after horizontal transport:', jk, max_val, min_val
! ! END DO
! ! 
! ! !-----------------------------------------------------------------------------------------------------------------------------
! ! !Step 2) Vertical advection and diffusion
! ! 
! ! 
! ! z_transport_w  = ab_gam*p_os%p_diag%w + (1.0_wp-ab_gam)*p_os%p_diag%w_old
! ! 
! ! SELECT CASE(FLUX_CALCULATION)
! ! 
! ! CASE(UPWIND)
! !   CALL upwind_vflux_oce( p_patch,      &
! !                        & trac_new,     &
! !                        & z_transport_w,&
! !                        & z_adv_flux_v )
! ! CASE(CENTRAL)
! !   CALL central_vflux_oce( p_patch,     &
! !                        & trac_new,     &
! !                        & z_transport_w,&
! !                        & z_adv_flux_v )
! ! CASE(MIMETIC)
! !   CALL mimetic_vflux_oce( p_patch,      &
! !                        & trac_new,     &
! !                        & z_transport_w,&
! !                        & z_adv_flux_v )
! ! END SELECT
! ! 
! ! !calculate mass flux
! ! DO jb = i_startblk_c, i_endblk_c
! !   CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
! !   &                   rl_start_c, rl_end_c)
! !   DO jc = i_startidx_c, i_endidx_c
! !     z_dolic = v_base%dolic_c(jc,jb)
! !     IF(z_dolic/=0)THEN
! !       DO jk = 1, z_dolic
! !         z_h2(jc,jk,jb) = z_h_tmp(jc,jk,jb)&
! !         & - (z_transport_w(jc,jk,jb) -z_transport_w(jc,jk+1,jb))
! !       END DO
! !       z_h2(jc,z_dolic,jb) = z_h_tmp(jc,z_dolic,jb)
! !    ENDIF
! !   END DO
! ! END DO
! ! 
! ! !The diffusion part
! ! IF(expl_vertical_tracer_diff==1)THEN
! ! 
! !   !vertical divergence is calculated for vertical advective fluxes
! !   DO jb = i_startblk_c, i_endblk_c
! !     CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
! !     &                   rl_start_c, rl_end_c)
! !     DO jc = i_startidx_c, i_endidx_c
! !       !interior: from one below surface to the ground
! !       z_dolic = v_base%dolic_c(jc,jb)
! !       IF(z_dolic>0)THEN
! !         DO jk = 1, z_dolic
! !           jkp1 = jk + 1
! !           !positive vertical divergence in direction of w (upward positive)
! !           z_div_adv_v(jc,jk,jb) = &
! !           & (z_adv_flux_v(jc,jk,jb)-z_adv_flux_v(jc,jkp1,jb))
! ! 
! !         END DO!level-loop
! !       ENDIF!(z_dolic>0)
! !     END DO!idx-loop
! !   END DO !blk-loop
! ! 
! ! !Final step: calculate new tracer values
! ! DO jk = 1, n_zlev 
! !   CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
! !                    & i_startidx_c, i_endidx_c,&
! !                    & rl_start_c, rl_end_c)
! !     DO jb = i_startblk_c, i_endblk_c
! !       DO jc = i_startidx_c, i_endidx_c
! !         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! !           IF(jk==1)THEN
! !             delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+v_base%del_zlev_m(top)
! !             !delta_z2  = p_os%p_prog(nnew(1))%h(jc,jb)+v_base%del_zlev_m(top)
! !           ELSE
! !             delta_z  = v_base%del_zlev_m(jk)  !delta_z2  = v_base%del_zlev_m(jk) 
! !           ENDIF
! !             z_trac_c(jc,jk,jb) = (trac_new(jc,jk,jb)*delta_z     &
! !                                & -dtime2*z_div_adv_v(jc,jk,jb))  &
! !                                &/(z_h2(jc,jk,jb)+z_tol)
! !         ENDIF
! !       END DO
! !     END DO
! ! END DO
! ! 
! !   CALL tracer_diffusion_vert_impl( p_patch,        &
! !                               & z_trac_c(:,:,:),&
! !                               & bc_top_tracer,  & 
! !                               & bc_bot_tracer,  &
! !                               & A_v,            &
! !                               & trac_new(:,:,:))
! ! 
! ! ELSEIF(expl_vertical_tracer_diff==0)THEN
! ! 
! !   !calculation of vertical diffusive fluxes
! !   CALL tracer_diffusion_vert_expl( p_patch,      &
! !                                 &  trac_old, z_h, &
! !                                 &  bc_top_tracer, &
! !                                 &  bc_bot_tracer,  & 
! !                                 &  A_v,           &
! !                                 &  z_diff_flux_v)
! ! 
! !   DO jk = 1, n_zlev+1
! !     write(*,*)'max/min adv & diff tracer flux v:',jk,&
! !     &maxval(z_adv_flux_v(:,jk,:)), minval(z_adv_flux_v(:,jk,:)),&
! !     &maxval(z_diff_flux_v(:,jk,:)),minval(z_diff_flux_v(:,jk,:))
! !   END DO
! !   DO jk = 1, n_zlev
! !     write(*,*)'max/min vertical mass flux:',jk,&
! !     &maxval(z_h2(:,jk,:)),minval(z_h2(:,jk,:))
! !   END DO
! ! 
! !   !vertical divergence is calculated for vertical advective & diffusive fluxes
! !   DO jb = i_startblk_c, i_endblk_c
! !     CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
! !     &                   rl_start_c, rl_end_c)
! !     DO jc = i_startidx_c, i_endidx_c
! !       !interior: from one below surface to the ground
! !       z_dolic = v_base%dolic_c(jc,jb)
! !       IF(z_dolic>0)THEN
! !         DO jk = 1, z_dolic
! !           jkp1 = jk + 1 
! ! 
! !           IF(jk==1)THEN
! !             delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+v_base%del_zlev_m(top)
! !           ELSE
! !             delta_z  = v_base%del_zlev_m(jk)
! !           ENDIF
! !           !positive vertical divergence in direction of w (upward positive)
! !           z_div_adv_v(jc,jk,jb) = &
! !           & (z_adv_flux_v(jc,jk,jb)-z_adv_flux_v(jc,jkp1,jb))
! ! 
! !           G_n_c_v(jc,jk,jb) = &   !z_div_diff_v(jc,jk,jb) = &
! !           & (z_diff_flux_v(jc,jk,jb)-z_diff_flux_v(jc,jkp1,jb))
! !         END DO!level-loop
! !       ENDIF!(z_dolic>0)
! !     END DO!idx-loop
! !   END DO !blk-loop
! ! 
! !   IF(is_initial_timestep(timestep))THEN
! !     G_nimd_c_v(:,:,:) = G_n_c_v(:,:,:)
! !   ELSE
! !     G_nimd_c_v(:,:,:) = (1.5_wp+AB_const)* G_n_c_v(:,:,:)   &
! !       &               - (0.5_wp+AB_const)*G_nm1_c_v(:,:,:)
! !   ENDIF
! ! 
! ! 
! ! 
! !   DO jk = 1, n_zlev
! !     write(*,*)'max/min div adv & diffusive  flux v:',&
! !     &jk,maxval(z_div_adv_v(:,jk,:)), minval(z_div_adv_v(:,jk,:)),&
! !     &maxval(G_n_c_v(:,jk,:)),minval(G_n_c_v(:,jk,:))
! !   END DO
! ! 
! !   !Final step: calculate new tracer values
! !   DO jk = 1, n_zlev 
! !     CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
! !                    & i_startidx_c, i_endidx_c,&
! !                    & rl_start_c, rl_end_c)
! !     DO jb = i_startblk_c, i_endblk_c
! !       DO jc = i_startidx_c, i_endidx_c
! !         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! ! 
! !           IF(jk==1)THEN
! !             delta_z  = p_os%p_prog(nold(1))%h(jc,jb)+v_base%del_zlev_m(top)
! !             !delta_z2  = p_os%p_prog(nnew(1))%h(jc,jb)+v_base%del_zlev_m(top)
! !           ELSE
! !             delta_z  = v_base%del_zlev_m(jk)
! !            !delta_z2  = v_base%del_zlev_m(jk) 
! !           ENDIF
! !             trac_new(jc,jk,jb) = (trac_new(jc,jk,jb)*delta_z     &
! !                                & -dtime2*(z_div_adv_v(jc,jk,jb)  &
! !                                & -G_nimd_c_v(jc,jk,jb)))         &
! !                                &/(z_h2(jc,jk,jb)+z_tol)
! !         ENDIF
! !       END DO
! !     END DO
! !   END DO
! ! 
! !   G_nm1_c_v(:,:,:)  = G_n_c_v(:,:,:)
! !   G_nimd_c_v(:,:,:) = 0.0_wp
! !   G_n_c_v(:,:,:)    = 0.0_wp
! ! 
! ! 
! ! ENDIF!(lvertical_diff_implicit)THEN
! ! 
! ! rl_start_c   = 1
! ! rl_end_c     = min_rlcell
! ! i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
! ! i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
! ! 
! !   DO jk = 1, n_zlev
! !   max_val=0.0_wp
! !   min_val=100.0_wp
! !   DO jb = i_startblk_c, i_endblk_c
! !     CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
! !     &                   rl_start_c, rl_end_c)
! !     DO jc = i_startidx_c, i_endidx_c
! !       IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! !         IF(trac_new(jc,jk,jb)>max_val)&
! !           &max_val=trac_new(jc,jk,jb)
! !         IF(trac_new(jc,jk,jb)<min_val)&
! !           &min_val=trac_new(jc,jk,jb)
! ! 
! !         IF(trac_new(jc,jk,jb)==0.0_wp)&
! !         & write(*,*)'tracer =0:',jc,jk,jb 
! !       ENDIF
! !      END DO
! !    END DO
! !    WRITE(*,*)'MAX/MIN tracer after vertical transport:', jk, max_val, min_val
! ! END DO
! ! !    DO jk = 1, n_zlev
! ! !     ctr=0
! ! !     ctr_total=0
! ! !     DO jb = i_startblk_c, i_endblk_c
! ! !       CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c,&
! ! !                        & i_startidx_c, i_endidx_c,&
! ! !                        & rl_start_c, rl_end_c)
! ! !    DO jc = i_startidx_c, i_endidx_c
! ! !      IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary )THEN
! ! !      ctr_total=ctr_total+1
! ! !      ENDIF
! ! !  
! ! ! !    IF(trac_new(jc,jk,jb)>t_prof(jk)) THEN
! ! ! !    ctr=ctr+1
! ! ! !  write(*,*)'indices',jc,jb,jk, v_base%lsm_oce_c(jc,top,jb),&
! ! ! ! ! !&trac_old(jc,jk,jb),
! ! ! ! ! &(delta_z/delta_z2)*trac_old(jc,top,jb),&
! ! ! !  & dtime*G_nimd_c(jc,top,jb),G_n_c_h(jc,top,jb),z_div_adv_h(jc,top,jb),&
! ! ! !  trac_new(jc,top,jb), ctr
! ! ! !  ENDIF
! ! !        END DO
! ! !    END DO
! ! !    write(*,*)'counter > initial:',jk, ctr, ctr_total
! ! !  END DO
! ! 
! ! 
! ! END SUBROUTINE advect_individual_tracer_ab_old
! ! ! !   !-----------------------------------------------------------------------


END MODULE mo_oce_tracer_transport
