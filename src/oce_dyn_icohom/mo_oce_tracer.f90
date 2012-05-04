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
  &                                     ab_const, ab_gam, expl_vertical_tracer_diff, &
  &                                     iswm_oce, &
  &                                     FLUX_CALCULATION_HORZ, FLUX_CALCULATION_VERT, &
  &                                     UPWIND, CENTRAL, MIMETIC, MIMETIC_MIURA
USE mo_physical_constants,        ONLY: tf
USE mo_math_constants,            ONLY: pi
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
USE mo_scalar_product,            ONLY:  map_cell2edges,map_edges2cell,map_edges2cell
USE mo_oce_math_operators,        ONLY: div_oce_3D
!USE mo_advection_utils,           ONLY: laxfr_upflux, laxfr_upflux_v
USE mo_oce_diffusion,             ONLY: tracer_diffusion_horz, tracer_diffusion_vert_expl,&
                                      & tracer_diffusion_vert_impl_hom
!USE mo_oce_ab_timestepping_mimetic,ONLY: l_STAGGERED_TIMESTEP
USE mo_intp_data_strc,            ONLY: p_int_state
USE mo_oce_tracer_transport_horz, ONLY: advect_diffuse_horizontal
USE mo_oce_tracer_transport_vert, ONLY: advect_diffuse_vertical
USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
USE mo_sync,                      ONLY: SYNC_C, SYNC_C1, SYNC_E, SYNC_V, sync_patch_array
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
PRIVATE :: prepare_tracer_transport

CONTAINS
!-------------------------------------------------------------------------
!
!>
!! !  SUBROUTINE advects the tracers present in the ocean model.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!!
SUBROUTINE advect_tracer_ab(p_patch, p_os, p_param, p_sfc_flx,p_op_coeff, timestep)
  TYPE(t_patch),             TARGET, INTENT(in)    :: p_patch
  TYPE(t_hydro_ocean_state), TARGET                :: p_os
  TYPE(t_ho_params),                 INTENT(INOUT) :: p_param
  TYPE(t_sfc_flx),                   INTENT(INOUT) :: p_sfc_flx
  TYPE(t_operator_coeff),            INTENT(INOUT) :: p_op_coeff
  INTEGER                                          :: timestep

  !Local variables
  INTEGER  :: i_no_t, jk
  REAL(wp) :: z_relax
  INTEGER  :: iloc(2)
  REAL(wp) :: zlat, zlon
  REAL(wp) :: z_cellthick_intmed(nproma,n_zlev, p_patch%nblks_c)
  REAL(wp) :: z_c(nproma,n_zlev,p_patch%nblks_c)
  !-------------------------------------------------------------------------------
!p_os%p_prog(nnew(1))%tracer=p_os%p_prog(nold(1))%tracer
!return
  z_cellthick_intmed=0.0_wp
  CALL prepare_tracer_transport( p_patch,           &
                               & p_os,              &
                               & p_param,           &
                               & p_sfc_flx,         &
                               & p_op_coeff,        &
                               & z_cellthick_intmed,&
                               & timestep)

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
                                 & p_os, p_op_coeff,  z_cellthick_intmed,     &
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



  DO jk = 1, n_zlev
    ! Abort if tracer is below treshold
    ! Temperature: tf<-1.9 deg, may be possible, limit set to lower value
    IF (minval(p_os%p_prog(nnew(1))%tracer(:,jk,:,1))<-4.0_wp) THEN
      write(0,*) ' TEMPERATURE BELOW TRESHOLD:'
      iloc(:) = minloc(p_os%p_prog(nnew(1))%tracer(:,jk,:,1))
      zlat    = p_patch%cells%center(iloc(1),iloc(2))%lat * 180.0_wp / pi
      zlon    = p_patch%cells%center(iloc(1),iloc(2))%lon * 180.0_wp / pi
      write(0,*) ' negative temperature at jk =', jk, &
      &minval(p_os%p_prog(nnew(1))%tracer(:,jk,:,1))
      write(0,*) ' location is at    idx =',iloc(1),' blk=',iloc(2)
      write(0,*) ' lat/lon  is at    lat =',zlat   ,' lon=',zlon
      CALL finish(TRIM('mo_tracer_advection:advect_individual_tracer-h'), &
        &              'Temperatur below treshold')
    ENDIF
  END DO

  IF (no_tracer>=2)THEN
    DO jk = 1, n_zlev
      ! Abort if salinity is negative:
      IF (minval(p_os%p_prog(nnew(1))%tracer(:,jk,:,2))<0.0_wp) THEN
        write(0,*) ' SALINITY NEGATIVE:'
        iloc(:) = minloc(p_os%p_prog(nnew(1))%tracer(:,jk,:,2))
        zlat    = p_patch%cells%center(iloc(1),iloc(2))%lat * 180.0_wp / pi
        zlon    = p_patch%cells%center(iloc(1),iloc(2))%lon * 180.0_wp / pi
        write(0,*) ' negative temperature at jk =', jk, &
        &minval(p_os%p_prog(nnew(1))%tracer(:,jk,:,2))
        write(0,*) ' location is at    idx =',iloc(1),' blk=',iloc(2)
        write(0,*) ' lat/lon  is at    lat =',zlat   ,' lon=',zlon
        CALL finish(TRIM('mo_tracer_advection:advect_individual_tracer-h'), &
        &              'SALINITY NEGATIVE')
      ENDIF
    END DO
  ENDIF
END SUBROUTINE advect_tracer_ab
!-------------------------------------------------------------------------
!
!
!>
!! !  SUBROUTINE prepares next tracer transport step. Currently needed in horizontal
!!    flux-scheme "MIMETIC-Miura". Geometric quantities are updated according to
!!    actual velocity. This information is required by MIURA-scheme and is identical
!!    for all tracers.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2012).
!!
!! mpi parallelized, sync required
SUBROUTINE prepare_tracer_transport(p_patch, p_os, p_param, p_sfc_flx, p_op_coeff,&
                                   &z_cellthick_intmed, timestep)
  TYPE(t_patch), TARGET, INTENT(in)    :: p_patch
  TYPE(t_hydro_ocean_state), TARGET    :: p_os
  TYPE(t_ho_params), INTENT(INOUT)     :: p_param
  TYPE(t_sfc_flx), INTENT(INOUT)       :: p_sfc_flx
  TYPE(t_operator_coeff),INTENT(INOUT) :: p_op_coeff
  REAL(wp),INTENT(INOUT)               :: z_cellthick_intmed(nproma,n_zlev, p_patch%nblks_c)
  INTEGER                              :: timestep
  !
  !Local variables
  REAL(wp) :: z_relax
  REAL(wp) :: z_c(nproma,n_zlev,p_patch%nblks_c)
  INTEGER  :: slev, elev
  INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, rl_start_c, rl_end_c
  INTEGER  :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, rl_start_e, rl_end_e
  INTEGER  :: je, jk, jb,jc         !< index of edge, vert level, block
  INTEGER  :: il_v1, il_v2, ib_v1, ib_v2!, il_e, ib_e
  INTEGER  :: il_c, ib_c
  REAL(wp) :: delta_z
  TYPE(t_cartesian_coordinates):: z_vn_c(nproma,n_zlev,p_patch%nblks_c)
  !TYPE(t_cartesian_coordinates) :: u_mean_cc(nproma,n_zlev,p_patch%nblks_e)
  !-------------------------------------------------------------------------------
  TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
  !-------------------------------------------------------------------------------
  cells_in_domain => p_patch%cells%in_domain
  edges_in_domain => p_patch%edges%in_domain

  slev = 1
  elev = n_zlev

   p_os%p_diag%w_time_weighted=ab_gam*p_os%p_diag%w + (1.0_wp-ab_gam)*p_os%p_diag%w_old
   !p_os%p_diag%w_time_weighted=p_os%p_diag%w

  IF(FLUX_CALCULATION_HORZ==MIMETIC_MIURA.OR.FLUX_CALCULATION_HORZ==MIMETIC)THEN
    DO jk = slev, elev
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        DO je = i_startidx_e, i_endidx_e
          IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN

            !Get indices of two adjacent vertices
            il_v1 = p_patch%edges%vertex_idx(je,jb,1)
            ib_v1 = p_patch%edges%vertex_blk(je,jb,1)
            il_v2 = p_patch%edges%vertex_idx(je,jb,2)
            ib_v2 = p_patch%edges%vertex_blk(je,jb,2)

            p_os%p_diag%p_vn_mean(je,jk,jb)%x = 0.5_wp*&
              &(p_os%p_diag%p_vn_dual(il_v1,jk,ib_v1)%x+p_os%p_diag%p_vn_dual(il_v2,jk,ib_v2)%x)

            p_op_coeff%moved_edge_position_cc(je,jk,jb)%x = &
              & p_op_coeff%edge_position_cc(je,jk,jb)%x     &
              &  - 0.5_wp*dtime*p_os%p_diag%p_vn_mean(je,jk,jb)%x

            IF ( p_os%p_diag%vn_time_weighted(je,jk,jb) > 0.0_wp ) THEN
              il_c = p_patch%edges%cell_idx(je,jb,1)
              ib_c = p_patch%edges%cell_blk(je,jb,1)
            ELSE  ! p_os%p_diag%vn_time_weighted <= 0.0
              il_c = p_patch%edges%cell_idx(je,jb,2)
              ib_c = p_patch%edges%cell_blk(je,jb,2)
            ENDIF

            p_op_coeff%upwind_cell_idx(je,jk,jb) = il_c
            p_op_coeff%upwind_cell_blk(je,jk,jb) = ib_c

            p_op_coeff%upwind_cell_position_cc(je,jk,jb)%x = &
              & p_op_coeff%cell_position_cc(il_c,jk,ib_c)%x

          ENDIF
        END DO
      END DO
    END DO
  ENDIF

  !CALL sync_patch_array(SYNC_E, p_patch,p_os%p_diag%h_e )
  !CALL sync_patch_array(SYNC_E, p_patch,p_os%p_diag%vn_time_weighted )
  ! CALL sync_patch_array(SYNC_E, p_patch,v_base%del_zlev_m )

!   !Calculation of mass flux and related quantities that are identical for all tracers
!   CALL map_edges2cell( p_patch, p_os%p_diag%vn_time_weighted, z_vn_c)
!   DO jb = cells_in_domain%start_block, cells_in_domain%end_block
!     CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
!     DO jc = i_startidx_c, i_endidx_c
!       DO jk = 1, n_zlev
!         delta_z = v_base%del_zlev_m(jk)
!           IF(jk==1)THEN
!             delta_z=v_base%del_zlev_m(jk)+p_os%p_prog(nold(1))%h(jc,jb)
!           ENDIF
!           z_vn_c(jc,jk,jb)%x = delta_z*z_vn_c(jc,jk,jb)%x
!       END DO
!     END DO
!   END DO
!   CALL map_cell2edges( p_patch, z_vn_c, p_os%p_diag%mass_flx_e)
  DO jb = edges_in_domain%start_block, edges_in_domain%end_block
    CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
    DO jk = 1, n_zlev
      delta_z = v_base%del_zlev_m(jk)
      DO je = i_startidx_e, i_endidx_e

      IF(jk==1) delta_z=v_base%del_zlev_m(jk) + p_os%p_diag%h_e(je,jb)

        !IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
          p_os%p_diag%mass_flx_e(je,jk,jb)  = delta_z*p_os%p_diag%vn_time_weighted(je,jk,jb)          
         !ENDIF
      END DO
    END DO
  END DO

  ! sync the result, is it necessary? this is not used elsewhere
  ! keep it only for checking
  CALL sync_patch_array(SYNC_E, p_patch,p_os%p_diag%mass_flx_e )
  
  CALL div_oce_3d(  p_os%p_diag%mass_flx_e,&
                 & p_patch,               &
                 & p_op_coeff%div_coeff,&
                 & p_os%p_diag%div_mass_flx_c,subset_range=cells_in_domain)

  !calculate (dummy) height consistent with divergence of mass fluxx (UPPERMOST LEVEL)
  jk = 1
  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
    DO jc = i_startidx_c, i_endidx_c
      !IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
      delta_z = v_base%del_zlev_m(jk)+p_os%p_prog(nold(1))%h(jc,jb)&
                &*v_base%wet_c(jc,jk,jb)

      p_os%p_diag%depth_c(jc,jk,jb) = delta_z

      z_cellthick_intmed(jc,jk,jb)&
      & = (delta_z-dtime*p_os%p_diag%div_mass_flx_c(jc,jk,jb))&
      &*v_base%wet_c(jc,jk,jb)
      !ENDIF
    END DO
  END DO

  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
    DO jk = 2, n_zlev
      delta_z = v_base%del_zlev_m(jk)
      DO jc = i_startidx_c, i_endidx_c
        !IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
        p_os%p_diag%depth_c(jc,jk,jb) = delta_z

        z_cellthick_intmed(jc,jk,jb)&
        & =(delta_z-dtime*p_os%p_diag%div_mass_flx_c(jc,jk,jb))&
        &*v_base%wet_c(jc,jk,jb)
        !ENDIF
      END DO
    END DO
  END DO

END SUBROUTINE prepare_tracer_transport
!-------------------------------------------------------------------------
!
!
!>
!! !  SUBROUTINE advects the tracers present in the ocean model.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!!
SUBROUTINE advect_individual_tracer_ab(p_patch, trac_old,                  &
                                     & p_os, p_op_coeff,z_cellthick_intmed,&
                                     & bc_top_tracer, bc_bot_tracer,       &
                                     & K_h, A_v,                           &
                                     & trac_new, timestep)

  TYPE(t_patch), TARGET, INTENT(in)    :: p_patch
  REAL(wp), INTENT(in)                 :: trac_old(nproma,n_zlev, p_patch%nblks_c)
  TYPE(t_hydro_ocean_state), TARGET    :: p_os
  TYPE(t_operator_coeff),INTENT(INOUT) :: p_op_coeff
  REAL(wp), INTENT(in)                 :: z_cellthick_intmed(nproma,n_zlev, p_patch%nblks_c)
  REAL(wp) :: bc_top_tracer(nproma, p_patch%nblks_c)
  REAL(wp) :: bc_bot_tracer(nproma, p_patch%nblks_c)
  REAL(wp) :: K_h(:,:,:)                                  !horizontal mixing coeff
  REAL(wp) :: A_v(:,:,:)                                   !vertical mixing coeff
  REAL(wp) :: trac_new(:,:,:)                              !new tracer
  INTEGER  :: timestep                                     ! Actual timestep (to distinghuish initial step from others)

  !Local variables
  REAL(wp) :: delta_t
  REAL(wp) :: trac_tmp(nproma,n_zlev, p_patch%nblks_c)
  INTEGER  :: jk
  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_tracer_advection:advect_individual_tracer')
  !-------------------------------------------------------------------------------
  delta_t  = dtime
  trac_tmp = 0.0_wp

  ipl_src=1  ! output print level (1-5, fix)
  CALL print_mxmn('on entry - Tracer',1,trac_old(:,:,:),n_zlev, p_patch%nblks_c,'trc',ipl_src)

  CALL advect_diffuse_horizontal(p_patch, trac_old,       &
                               & p_os,p_op_coeff,         &
                               & K_h,                     &
                               & trac_tmp,                &
                               & timestep, delta_t,       &
                               & p_os%p_diag%h_e,         &
                               & z_cellthick_intmed,      &
                               & FLUX_CALCULATION_HORZ)

!     DO jk = 1, n_zlev
!        write(*,*)'After horizontal max/min old-new tracer:',jk, maxval(trac_old(:,jk,:)),&
!                                              & minval(trac_old(:,jk,:)),&
!                                              & maxval(trac_tmp(:,jk,:)),&
!                                              & minval(trac_tmp(:,jk,:))
!     END DO

  IF( iswm_oce /= 1) THEN

    CALL advect_diffuse_vertical(p_patch, trac_tmp,       &
                        & p_os,                           &
                        & bc_top_tracer, bc_bot_tracer,   &
                        & A_v,                            &
                        & trac_new, timestep, delta_t,    &
                        & z_cellthick_intmed,             &
                        & FLUX_CALCULATION_VERT)
!    DO jk = 1, n_zlev
!      write(*,*)'After vertical max/min old-new tracer:',jk,&
!      & maxval(trac_old(:,jk,:)), minval(trac_old(:,jk,:)),&
!      & maxval(trac_new(:,jk,:)),minval(trac_new(:,jk,:))
!           ipl_src=3  ! output print level (1-5, fix)
!           CALL print_mxmn('adv-vert trac-old',jk,trac_old(:,:,:),n_zlev, &
!             &              p_patch%nblks_c,'trc',ipl_src)
!           CALL print_mxmn('adv-vert trac-new',jk,trac_new(:,:,:),n_zlev, &
!             &              p_patch%nblks_c,'trc',ipl_src)
!      END DO
  ELSEIF( iswm_oce == 1) THEN

    trac_new = trac_tmp

    DO jk = 1, n_zlev
      ipl_src=3  ! output print level (1-5, fix)
      CALL print_mxmn('adv-vert trac-old',jk,trac_old(:,:,:),n_zlev, &
        &              p_patch%nblks_c,'trc',ipl_src)
      CALL print_mxmn('adv-vert trac-new',jk,trac_new(:,:,:),n_zlev, &
        &              p_patch%nblks_c,'trc',ipl_src)
    END DO

  ENDIF


END SUBROUTINE advect_individual_tracer_ab

END MODULE mo_oce_tracer
