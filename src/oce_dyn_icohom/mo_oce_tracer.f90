!>
!! Contains the implementation of the tracer transport routines for the ICON ocean model.
!! This comprises advection and diffusion in horizontal and vertical direction.
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
MODULE mo_oce_tracer
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
  &                                     min_rlcell, min_rledge, min_rlcell,MIN_DOLIC
USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer, idisc_scheme,    &
  &                                     irelax_3d_T, relax_3d_mon_T, irelax_3d_S, relax_3d_mon_S, &
  &                                     ab_const, ab_gam, expl_vertical_tracer_diff,&
  &                                     iswm_oce
USE mo_physical_constants,        ONLY: tf
USE mo_parallel_config,           ONLY: nproma
USE mo_dynamics_config,           ONLY: nold, nnew 
USE mo_run_config,                ONLY: dtime
USE mo_oce_state,                 ONLY: t_hydro_ocean_state, v_base, is_initial_timestep
USE mo_model_domain,              ONLY: t_patch
USE mo_exception,                 ONLY: finish !, message_text, message
USE mo_oce_index,                 ONLY: print_mxmn, jkc, jkdim, ipl_src
USE mo_loopindices,               ONLY: get_indices_c, get_indices_e !, get_indices_v
USE mo_oce_boundcond,             ONLY: top_bound_cond_tracer
USE mo_oce_physics
USE mo_sea_ice,                   ONLY: t_sfc_flx
!USE mo_scalar_product,            ONLY:  map_cell2edges,map_edges2cell,map_edges2cell
!USE mo_oce_math_operators,        ONLY: div_oce, grad_fd_norm_oce, grad_fd_norm_oce_2d
!USE mo_advection_utils,           ONLY: laxfr_upflux, laxfr_upflux_v
USE mo_oce_diffusion,             ONLY: tracer_diffusion_horz, tracer_diffusion_vert_expl,&
                                      & tracer_diffusion_vert_impl_hom
!USE mo_oce_ab_timestepping_mimetic,ONLY: l_STAGGERED_TIMESTEP
USE mo_intp_data_strc,             ONLY: p_int_state
USE mo_oce_tracer_transport_horz,    ONLY: advect_horizontal
USE mo_oce_tracer_transport_vert,    ONLY: advect_vertical

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
PRIVATE :: advect_individual_tracer_ab


INTEGER, PARAMETER  :: top=1
INTEGER, PARAMETER :: UPWIND = 1
INTEGER, PARAMETER :: CENTRAL= 2
INTEGER, PARAMETER :: MIMETIC= 3
INTEGER, PARAMETER :: MIMETIC_MIURA= 4


INTEGER  :: FLUX_CALCULATION_HORZ=MIMETIC
INTEGER  :: FLUX_CALCULATION_VERT=MIMETIC

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
INTEGER  :: i_no_t, jk
REAL(wp) :: z_relax
REAL(wp) :: z_c(nproma,n_zlev,p_patch%nblks_c)
!-------------------------------------------------------------------------------

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
!                                & p_os%p_aux%g_n_c_h(:,:,:,i_no_t),          &
!                                & p_os%p_aux%g_nm1_c_h(:,:,:,i_no_t),        &
!                                & p_os%p_aux%g_nimd_c_h(:,:,:,i_no_t),       &
!                                & p_os%p_aux%g_n_c_v(:,:,:,i_no_t),          &
!                                & p_os%p_aux%g_nm1_c_v(:,:,:,i_no_t),        &
!                                & p_os%p_aux%g_nimd_c_v(:,:,:,i_no_t),       &
                               & p_os%p_aux%bc_top_tracer(:,:,i_no_t),      &
                               & p_os%p_aux%bc_bot_tracer(:,:,i_no_t),      &
                               & p_param%K_tracer_h(:,:,:,i_no_t ),         &
                               & p_param%A_tracer_v(:,:,:, i_no_t),         &
                               & p_os%p_prog(nnew(1))%tracer(:,:,:,i_no_t), &
                               & timestep )
END DO

! Final step: 3-dim temperature relaxation
!  - strict time constant, i.e. independent of layer thickness 
!  - additional forcing Term F_T = -1/tau(T-T*) [ K/s ]
!    when using the sign convention
!      dT/dt = Operators + F_T
!    i.e. F_T <0 for  T-T* >0 (i.e. decreasing temperature if it is warmer than relaxation data) 
!  - discretized: 
!    tracer = tracer - 1/(relax_3d_mon_T[months]) * (tracer(1)-relax_3d_data_T)
IF (no_tracer>=1 .AND. irelax_3d_T >0) THEN

  ! calculate relaxation term
  z_relax = 1.0_wp/(relax_3d_mon_T*2.592e6_wp)
  p_os%p_aux%relax_3d_forc_T(:,:,:) = -z_relax* &
    &  ( p_os%p_prog(nnew(1))%tracer(:,:,:,1) - p_os%p_aux%relax_3d_data_T(:,:,:))

  ! add relaxation term to new temperature
  p_os%p_prog(nnew(1))%tracer(:,:,:,1) = p_os%p_prog(nnew(1))%tracer(:,:,:,1) - &
    &                                    p_os%p_aux%relax_3d_forc_T(:,:,:) * dtime

  DO jk = 1, n_zlev 
    ipl_src=3  ! output print level (1-5, fix)
    CALL print_mxmn('3d_relax_T: forc',jk, p_os%p_aux%relax_3d_forc_T(:,:,:),n_zlev, &
      &              p_patch%nblks_c,'trc',ipl_src)
    CALL print_mxmn('3d_relax_T: data',jk, p_os%p_aux%relax_3d_data_T(:,:,:),n_zlev, &
      &              p_patch%nblks_c,'trc',ipl_src)
    ipl_src=2  ! output print level (1-5, fix)
    z_c(:,:,:) =  p_os%p_prog(nnew(1))%tracer(:,:,:,1)
    CALL print_mxmn('3d_relax_T: trac',jk,z_c(:,:,:),n_zlev, &
      &              p_patch%nblks_c,'trc',ipl_src)
  END DO

END IF

! Final step: 3-dim salinity relaxation
!  - additional forcing Term F_S = -1/tau(S-S*) [ psu/s ]
!    when using the sign convention
!      dS/dt = Operators + F_S
!    i.e. F_S <0 for  S-S* >0 (i.e. decreasing salinity if it is larger than relaxation data) 
!    note that freshwater flux is positive to decrease salinity, i.e. freshening water
!  - discretized: 
!    tracer = tracer - 1/(relax_3d_mon_T[months]) * (tracer(1)-relax_3d_data_T)
IF (no_tracer==2 .AND. irelax_3d_S >0) THEN

  ! calculate relaxation term
  z_relax = 1.0_wp/(relax_3d_mon_S*2.592e6_wp)
  p_os%p_aux%relax_3d_forc_S(:,:,:) = -z_relax* &
    &  ( p_os%p_prog(nnew(1))%tracer(:,:,:,2) - p_os%p_aux%relax_3d_data_S(:,:,:))

  ! add relaxation term to new salinity
  p_os%p_prog(nnew(1))%tracer(:,:,:,2) = p_os%p_prog(nnew(1))%tracer(:,:,:,2) + &
    &                                    p_os%p_aux%relax_3d_forc_S(:,:,:) * dtime

  DO jk = 1, n_zlev 
    ipl_src=3  ! output print level (1-5, fix)
    CALL print_mxmn('3d_relax_S: forc',jk, p_os%p_aux%relax_3d_forc_S(:,:,:),n_zlev, &
      &              p_patch%nblks_c,'trc',ipl_src)
    CALL print_mxmn('3d_relax_S: data',jk, p_os%p_aux%relax_3d_data_S(:,:,:),n_zlev, &
      &              p_patch%nblks_c,'trc',ipl_src)
    ipl_src=2  ! output print level (1-5, fix)
    z_c(:,:,:) =  p_os%p_prog(nnew(1))%tracer(:,:,:,2)
    CALL print_mxmn('3d_relax_S: trac',jk,z_c(:,:,:),n_zlev, &
      &              p_patch%nblks_c,'trc',ipl_src)
  END DO
END IF

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
                                     & p_os, &!G_n_c_h,G_nm1_c_h,G_nimd_c_h, &
                                     !& G_n_c_v, G_nm1_c_v, G_nimd_c_v, &
                                     & bc_top_tracer, bc_bot_tracer,&
                                     & K_h, A_v,                    &
                                     & trac_new, timestep)
!
!
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp)                          :: trac_old(:,:,:)
TYPE(t_hydro_ocean_state), TARGET :: p_os
! REAL(wp) :: G_n_c_h   (nproma, n_zlev,   p_patch%nblks_c)  !G^n
! REAL(wp) :: G_nm1_c_h (nproma, n_zlev,   p_patch%nblks_c)  !G^(n-1)
! REAL(wp) :: G_nimd_c_h(nproma, n_zlev,   p_patch%nblks_c)  !G^(n+1/2)
! REAL(wp) :: G_n_c_v   (nproma, n_zlev,   p_patch%nblks_c)  !G^n
! REAL(wp) :: G_nm1_c_v (nproma, n_zlev,   p_patch%nblks_c)  !G^(n-1)
! REAL(wp) :: G_nimd_c_v(nproma, n_zlev,   p_patch%nblks_c)  !G^(n+1/2)
REAL(wp) :: bc_top_tracer(nproma, p_patch%nblks_c)
REAL(wp) :: bc_bot_tracer(nproma, p_patch%nblks_c)
REAL(wp) :: K_h(:,:,:)                                  !horizontal mixing coeff
REAL(wp) :: A_v(:,:,:)                                   !vertical mixing coeff
REAL(wp) :: trac_new(:,:,:)                              !new tracer 
INTEGER  :: timestep                                     ! Actual timestep (to distinghuish initial step from others)
!
!Local variables
REAL(wp) :: delta_t
REAL(wp) :: dummy_h_c(nproma,n_zlev, p_patch%nblks_c)
REAL(wp) :: trac_tmp(nproma,n_zlev, p_patch%nblks_c)
INTEGER  :: jk
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_tracer_advection:advect_individual_tracer')
!-------------------------------------------------------------------------------
delta_t= dtime
dummy_h_c= 0.0_wp
trac_tmp = 0.0_wp

ipl_src=1  ! output print level (1-5, fix)
CALL print_mxmn('on entry - Tracer',1,trac_old(:,:,:),n_zlev, p_patch%nblks_c,'trc',ipl_src)

CALL advect_horizontal(p_patch, trac_old,           &
                     & p_os,&! G_n_c_h,G_nm1_c_h,G_nimd_c_h, &
                     & K_h,                    &
                     & trac_tmp, timestep, delta_t,p_os%p_diag%h_e, dummy_h_c,&
                     & FLUX_CALCULATION_HORZ)
  !write(123,*)'-------------timestep---------------',timestep
! DO jk = 1, n_zlev
! !   ipl_src=3  ! output print level (1-5, fix)
! !     CALL print_mxmn('adv-horz trac-old',jk,trac_old(:,:,:),n_zlev, &
! !       &              p_patch%nblks_c,'trc',ipl_src)
! !     CALL print_mxmn('adv-horz trac-tmp',jk,trac_tmp(:,:,:),n_zlev, &
! !       &              p_patch%nblks_c,'trc',ipl_src)
!   write(*,*)'After horizontal max/min old-new tracer:',jk, maxval(trac_old(:,jk,:)),&
!                                         & minval(trac_old(:,jk,:)),&
!                                         & maxval(trac_tmp(:,jk,:)),&
!                                         & minval(trac_tmp(:,jk,:))
!  END DO

IF( iswm_oce /= 1) THEN

     CALL advect_vertical(p_patch, trac_tmp,              &
                         & p_os,                           &
!                         & G_n_c_v, G_nm1_c_v, G_nimd_c_v, &
                         & bc_top_tracer, bc_bot_tracer,   &
                         & A_v,                            &
                         & trac_new, timestep, delta_t, dummy_h_c,&
                         & FLUX_CALCULATION_VERT)
!     DO jk = 1, n_zlev 
!       write(*,*)'After vertical max/min old-new tracer:',jk,&
!       & maxval(trac_old(:,jk,:)), minval(trac_old(:,jk,:)),&
!       & maxval(trac_new(:,jk,:)),minval(trac_new(:,jk,:))
! ! !        ipl_src=3  ! output print level (1-5, fix)
! ! !        CALL print_mxmn('adv-vert trac-old',jk,trac_old(:,:,:),n_zlev, &
! ! !          &              p_patch%nblks_c,'trc',ipl_src)
! ! !        CALL print_mxmn('adv-vert trac-new',jk,trac_new(:,:,:),n_zlev, &
! ! !          &              p_patch%nblks_c,'trc',ipl_src)
!      END DO
ELSEIF( iswm_oce == 1) THEN

  trac_new=trac_tmp

  DO jk = 1, n_zlev
    ipl_src=3  ! output print level (1-5, fix)
    CALL print_mxmn('adv-vert trac-old',jk,trac_old(:,:,:),n_zlev, &
      &              p_patch%nblks_c,'trc',ipl_src)
    CALL print_mxmn('adv-vert trac-new',jk,trac_new(:,:,:),n_zlev, &
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
 
END MODULE mo_oce_tracer
