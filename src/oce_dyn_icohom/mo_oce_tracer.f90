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
!USE mo_math_utilities,            ONLY: t_cartesian_coordinates
USE mo_impl_constants,            ONLY: sea_boundary, sea
USE mo_math_constants,            ONLY: pi
USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer, &
  &                                     irelax_3d_T, relax_3d_mon_T, irelax_3d_S, relax_3d_mon_S, &
  &                                     expl_vertical_tracer_diff, iswm_oce,&! l_edge_based,    &
  &                                     FLUX_CALCULATION_HORZ, &!FLUX_CALCULATION_VERT, &
  &                                     MIMETIC_MIURA
USE mo_util_dbg_prnt,             ONLY: dbg_print
USE mo_parallel_config,           ONLY: nproma
USE mo_dynamics_config,           ONLY: nold, nnew
USE mo_run_config,                ONLY: dtime, ltimer
USE mo_oce_state,                 ONLY: t_hydro_ocean_state!, v_base
USE mo_model_domain,              ONLY: t_patch, t_patch_3D_oce
USE mo_exception,                 ONLY: finish !, message_text, message
!USE mo_oce_index,                 ONLY: print_mxmn, jkc, jkdim, ipl_src
!USE mo_loopindices,               ONLY: get_indices_c, get_indices_e !, get_indices_v
USE mo_oce_boundcond,             ONLY: top_bound_cond_tracer
USE mo_oce_physics
USE mo_sea_ice_types,             ONLY: t_sfc_flx
!USE mo_scalar_product,            ONLY:  map_cell2edges_3D,map_edges2cell_3D
!USE mo_oce_math_operators,        ONLY: div_oce_3D
USE mo_oce_diffusion,             ONLY: tracer_diffusion_vert_impl_hom!, tracer_diffusion_vert_expl,&
USE mo_oce_tracer_transport_horz, ONLY: advect_diffuse_flux_horz
USE mo_oce_tracer_transport_vert, ONLY: advect_flux_vertical 
USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
USE mo_sync,                      ONLY: SYNC_C, SYNC_E, sync_patch_array
USE mo_timer,                     ONLY: timer_start, timer_stop,&
  &                                     timer_dif_vert
IMPLICIT NONE

PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'
CHARACTER(len=12)           :: str_module = 'oceTracer   '  ! Output of module for 1 line debug
INTEGER                     :: idt_src    = 1               ! Level of detail for 1 line debug

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
SUBROUTINE advect_tracer_ab(p_patch_3D, p_os, p_param, p_sfc_flx,p_op_coeff, timestep)
  TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)      :: p_patch_3D
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
  REAL(wp) :: z_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)
  TYPE(t_patch), POINTER         :: p_patch 
  !-------------------------------------------------------------------------------
  p_patch => p_patch_3D%p_patch_2D(1)

  CALL prepare_tracer_transport( p_patch_3D,&
                               & p_os,              &
                               & p_op_coeff)

  DO i_no_t = 1,no_tracer
!CALL sync_patch_array(SYNC_C, p_patch,p_os%p_prog(nold(1))%tracer(:,:,:,i_no_t))
    !First tracer is temperature
    !Second tracer is salinity
    IF( iswm_oce /= 1) THEN
      CALL top_bound_cond_tracer( p_patch,  &
                                & p_os,     &
                                & i_no_t,   &
                                & p_sfc_flx,&
                                & p_os%p_aux%bc_top_tracer)
    ENDIF

    CALL advect_individual_tracer_ab( p_patch_3D,                             &
                                 & p_os%p_prog(nold(1))%tracer(:,:,:,i_no_t), &
                                 & p_os, p_op_coeff,                          &
                                 & p_os%p_aux%bc_top_tracer(:,:,i_no_t),      &
                                 & p_os%p_aux%bc_bot_tracer(:,:,i_no_t),      &
                                 & p_param%K_tracer_h(:,:,:,i_no_t ),         &
                                 & p_param%A_tracer_v(:,:,:,i_no_t),          &
                                 & p_os%p_prog(nnew(1))%tracer(:,:,:,i_no_t), &
                                 & i_no_t )
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

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('3d_rel: AdvTracT forc'  ,p_os%p_aux%relax_3d_forc_T   ,str_module,idt_src)
    CALL dbg_print('3d_rel: AdvTracT data'  ,p_os%p_aux%relax_3d_data_T   ,str_module,idt_src)
    idt_src=2  ! output print level (1-5, fix)
    z_c(:,:,:) =  p_os%p_prog(nnew(1))%tracer(:,:,:,1)
    CALL dbg_print('3d_relax: AdvTracT trac'  ,z_c                          ,str_module,idt_src)
    !---------------------------------------------------------------------

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
    p_os%p_aux%relax_3d_forc_S(1:nproma,jk,1:p_patch%nblks_c) = -z_relax* &
    &  ( p_os%p_prog(nnew(1))%tracer(1:nproma,jk,1:p_patch%nblks_c,2) - p_os%p_aux%relax_3d_data_S(1:nproma,jk,1:p_patch%nblks_c))

    ! add relaxation term to new salinity
    p_os%p_prog(nnew(1))%tracer(1:nproma,jk,1:p_patch%nblks_c,2) = p_os%p_prog(nnew(1))%tracer(1:nproma,jk,1:p_patch%nblks_c,2) + &
      &                                    p_os%p_aux%relax_3d_forc_S(1:nproma,jk,1:p_patch%nblks_c) * dtime

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('3d_rel: AdvTracS forc'  ,p_os%p_aux%relax_3d_forc_S   ,str_module,idt_src)
    CALL dbg_print('3d_rel: AdvTracS data'  ,p_os%p_aux%relax_3d_data_S   ,str_module,idt_src)
    idt_src=2  ! output print level (1-5, fix)
    z_c(:,:,:) =  p_os%p_prog(nnew(1))%tracer(:,:,:,2)
    CALL dbg_print('3d_relax: AdvTracS trac'  ,z_c                          ,str_module,idt_src)
    !---------------------------------------------------------------------

  END IF
!!Commented out because of NAG-compiler, PK 
  IF (no_tracer>=1) THEN
  DO jk = 1, n_zlev
    ! Abort if tracer is below threshold
    ! Temperature: <-1.9 deg, may be possible, limit set to lower value
    IF (minval(p_os%p_prog(nnew(1))%tracer(1:nproma,jk,1:p_patch%nblks_c,1))<-4.0_wp) THEN
      write(0,*) ' TEMPERATURE BELOW THRESHOLD:'
      iloc(:) = minloc(p_os%p_prog(nnew(1))%tracer(:,jk,:,1))
      zlat    = p_patch%cells%center(iloc(1),iloc(2))%lat * 180.0_wp / pi
      zlon    = p_patch%cells%center(iloc(1),iloc(2))%lon * 180.0_wp / pi
      write(0,*) ' negative temperature at jk =', jk, &
      &minval(p_os%p_prog(nnew(1))%tracer(:,jk,:,1))
      write(0,*) ' location is at    idx =',iloc(1),' blk=',iloc(2)
      write(0,*) ' lat/lon  is at    lat =',zlat   ,' lon=',zlon
      CALL finish(TRIM('mo_tracer_advection:advect_tracer'), &
        &              'Temperature below threshold')
    ENDIF
  END DO
  ENDIF
  IF (no_tracer>=2)THEN
    DO jk = 1, n_zlev
      ! Abort if salinity is negative:
      IF (minval(p_os%p_prog(nnew(1))%tracer(1:nproma,jk,1:p_patch%nblks_c,2))<0.0_wp) THEN
        write(0,*) ' SALINITY NEGATIVE:'
        iloc(:) = minloc(p_os%p_prog(nnew(1))%tracer(:,jk,:,2))
        zlat    = p_patch%cells%center(iloc(1),iloc(2))%lat * 180.0_wp / pi
        zlon    = p_patch%cells%center(iloc(1),iloc(2))%lon * 180.0_wp / pi
        write(0,*) ' negative salinity at jk =', jk, &
        &minval(p_os%p_prog(nnew(1))%tracer(:,jk,:,2))
        write(0,*) ' location is at    idx =',iloc(1),' blk=',iloc(2)
        write(0,*) ' lat/lon  is at    lat =',zlat   ,' lon=',zlon
        CALL finish(TRIM('mo_tracer_advection:advect_tracer'), &
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
SUBROUTINE prepare_tracer_transport(p_patch_3D, p_os, p_op_coeff)

  TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
  TYPE(t_hydro_ocean_state), TARGET    :: p_os
  TYPE(t_operator_coeff),INTENT(INOUT) :: p_op_coeff
  !
  !Local variables
  INTEGER  :: slev, elev
  INTEGER  :: i_startidx_e, i_endidx_e
  INTEGER  :: il_c1, il_c2, ib_c1, ib_c2
  INTEGER  :: il_c, ib_c
  INTEGER  :: je, jk, jb!,jc         !< index of edge, vert level, block
  TYPE(t_patch), POINTER :: p_patch 
  TYPE(t_subset_range), POINTER :: edges_in_domain
  !-------------------------------------------------------------------------------
  p_patch         => p_patch_3D%p_patch_2D(1)
  edges_in_domain => p_patch%edges%in_domain

  slev = 1
  elev = n_zlev

  p_os%p_diag%w_time_weighted(1:nproma,1:n_zlev+1,1:p_patch%nblks_c)= 0.0_wp

  IF( FLUX_CALCULATION_HORZ==MIMETIC_MIURA)THEN
    DO jk = slev, elev
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        DO je = i_startidx_e, i_endidx_e
          !IF(v_base%lsm_oce_e(je,jk,jb) == sea)THEN
          IF (p_patch_3D%lsm_oce_e(je,jk,jb) == sea) THEN

            !Get indices of two adjacent vertices
!             il_v1 = p_patch%edges%vertex_idx(je,jb,1)
!             ib_v1 = p_patch%edges%vertex_blk(je,jb,1)
!             il_v2 = p_patch%edges%vertex_idx(je,jb,2)
!             ib_v2 = p_patch%edges%vertex_blk(je,jb,2)
            !Get indices of two adjacent cells
            il_c1 = p_patch%edges%cell_idx(je,jb,1)
            ib_c1 = p_patch%edges%cell_blk(je,jb,1)
            il_c2 = p_patch%edges%cell_idx(je,jb,2)
            ib_c2 = p_patch%edges%cell_blk(je,jb,2)

            !  p_os%p_diag%p_vn_mean(je,jk,jb)%x = 0.5_wp*&
            !    &(p_os%p_diag%p_vn_dual(il_v1,jk,ib_v1)%x+p_os%p_diag%p_vn_dual(il_v2,jk,ib_v2)%x)
            p_os%p_diag%p_vn_mean(je,jk,jb)%x = 0.5_wp*&
              &(p_os%p_diag%p_vn(il_c1,jk,ib_c1)%x+p_os%p_diag%p_vn(il_c2,jk,ib_c2)%x)

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
      CALL sync_patch_array(SYNC_E,p_patch,p_op_coeff%upwind_cell_position_cc(1:nproma,jk,1:p_patch%nblks_e)%x(1))
      CALL sync_patch_array(SYNC_E,p_patch,p_op_coeff%upwind_cell_position_cc(1:nproma,jk,1:p_patch%nblks_e)%x(2))
      CALL sync_patch_array(SYNC_E,p_patch,p_op_coeff%upwind_cell_position_cc(1:nproma,jk,1:p_patch%nblks_e)%x(3))
    END DO

  ENDIF

  CALL sync_patch_array(SYNC_E, p_patch,p_os%p_diag%vn_time_weighted ) 
  CALL sync_patch_array(SYNC_C, p_patch,p_os%p_diag%w_time_weighted )


  !---------DEBUG DIAGNOSTICS-------------------------------------------
!   idt_src=4  ! output print level (1-5, fix)
!   CALL dbg_print('PrepTrans: prism thickness'  ,p_os%p_diag%prism_thick_c      ,str_module,idt_src)
  !---------------------------------------------------------------------

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
SUBROUTINE advect_individual_tracer_ab(p_patch_3D, trac_old,               &
                                     & p_os, p_op_coeff,                   &
                                     & bc_top_tracer, bc_bot_tracer,       &
                                     & K_h, A_v,                           &
                                     & trac_new, tracer_id)

  TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
  REAL(wp), INTENT(INOUT)              :: trac_old(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%nblks_c)
  TYPE(t_hydro_ocean_state), TARGET    :: p_os
  TYPE(t_operator_coeff),INTENT(INOUT) :: p_op_coeff
  REAL(wp), INTENT(in)                 :: bc_top_tracer(nproma, p_patch_3D%p_patch_2D(1)%nblks_c)
  REAL(wp), INTENT(in)                 :: bc_bot_tracer(nproma, p_patch_3D%p_patch_2D(1)%nblks_c)
  REAL(wp), INTENT(in)                 :: K_h(:,:,:)       !horizontal mixing coeff
  REAL(wp), INTENT(in)                 :: A_v(:,:,:)       !vertical mixing coeff
  REAL(wp), INTENT(INOUT)              :: trac_new(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_c)  !new tracer
  INTEGER,  INTENT(IN)                 :: tracer_id
  !Local variables
  REAL(wp) :: delta_t, delta_z
  REAL(wp) :: flux_horz(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%nblks_c)
  REAL(wp) :: flux_vert(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%nblks_c)
  REAL(wp) :: z_temp(nproma,n_zlev,    p_patch_3D%p_patch_2D(1)%nblks_c)
  INTEGER  :: jc,jk,jb
  INTEGER  :: z_dolic
  INTEGER  :: i_startidx_c, i_endidx_c
  TYPE(t_subset_range), POINTER :: cells_in_domain
  TYPE(t_patch), POINTER         :: p_patch 
  ! CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_tracer_advection:advect_individual_tracer')
  !-------------------------------------------------------------------------------_
  p_patch => p_patch_3D%p_patch_2D(1)
  delta_t = dtime

  z_temp(1:nproma,1:n_zlev,1:p_patch%nblks_c)     =0.0_wp
  flux_horz(1:nproma,1:n_zlev,1:p_patch%nblks_c)  =0.0_wp
  flux_vert(1:nproma,1:n_zlev,1:p_patch%nblks_c)  =0.0_wp

  trac_new(1:nproma,1:n_zlev,1:p_patch%nblks_c)=0.0_wp

  cells_in_domain => p_patch%cells%in_domain

  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=1  ! output print level (1-5, fix)
  CALL dbg_print('on entry: IndTrac: trac_old',trac_old(1:nproma,1:n_zlev,1:p_patch%nblks_c) ,str_module,idt_src)
  !---------------------------------------------------------------------

  CALL advect_diffuse_flux_horz( p_patch_3D,       &
                               & trac_old,         &
                               & p_os,             &
                               & p_op_coeff,       &
                               & K_h,              &
                               & flux_horz)

  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=3  ! output print level (1-5, fix)
  CALL dbg_print('after AdvDiffHorz: flux horz',flux_horz,str_module,idt_src)
  !---------------------------------------------------------------------

  !Shallow water is done with horizontal advection
  IF(iswm_oce == 1) THEN

    jk=1
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
          !IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          IF ( p_patch_3D%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
            delta_z =  p_patch_3D%p_patch_1D(1)%del_zlev_m(1)!v_base%del_zlev_m(1)
            trac_new(jc,jk,jb)= trac_old(jc,jk,jb)&
            & +(delta_t/delta_z)*flux_horz(jc,jk,jb)
          ENDIF
        END DO
      END DO

  !The 3D-case: first vertical fluxes than preliminary tracer value and
  !finally implicit vertical diffusion
  ELSEIF( iswm_oce /= 1) THEN
!write(*,*)'dolic_c21',maxval(p_patch_3D%p_patch_1D(1)%dolic_c),minval(p_patch_3D%p_patch_1D(1)%dolic_c)
!write(*,*)'dolic_e21',maxval(p_patch_3D%p_patch_1D(1)%dolic_e),minval(p_patch_3D%p_patch_1D(1)%dolic_e)
    CALL advect_flux_vertical( p_patch_3D,        &
                             & trac_old,          &
                             & p_os,              &
                             & bc_top_tracer,     &
                             & bc_bot_tracer,     &
                             & flux_vert,         &
                             & tracer_id)
!       DO jk=1,n_zlev
!        write(*,*)'FLX:HORZ:VERT',jk,&
!        &minval(flux_horz(:,jk,:)),maxval(flux_horz(:,jk,:)),&
!        &minval(flux_vert(:,jk,:)),maxval(flux_vert(:,jk,:))
!        END DO

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('after AdvDiffVert: flux vert',flux_vert,str_module,idt_src)
    !---------------------------------------------------------------------

    !Case: Implicit Vertical diffusion
    IF(expl_vertical_tracer_diff==1)THEN

      !Calculate preliminary tracer value out of horizontal advective and
      !diffusive fluxes and vertical advective fluxes, plus surface forcing.
      !Surface forcing applied as volume forcing at rhs, i.e.part of explicit term 
      !in tracer (and also momentum) eqs. In this case, top boundary condition of 
      !vertical Laplacians are homogeneous
      jk=1
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
          DO jc = i_startidx_c, i_endidx_c
            z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)
            !IF(z_dolic>=MIN_DOLIC)THEN
            !IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
            IF ( p_patch_3D%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
              !delta_z = p_os%p_diag%inv_prism_thick_c(jc,jk,jb)!v_base%del_zlev_m(jk)!+p_os%p_prog(nold(1))%h(jc,jb)
              delta_z = p_patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,jk,jb)

                z_temp(jc,jk,jb)= trac_old(jc,jk,jb) &
                & -(delta_t*delta_z)*(flux_vert(jc,jk,jb)-flux_horz(jc,jk,jb))

               z_temp(jc,jk,jb)=z_temp(jc,jk,jb)+(delta_t*delta_z)*bc_top_tracer(jc,jb)
            ENDIF
        END DO
      END DO

      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
          DO jc = i_startidx_c, i_endidx_c
            z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)! v_base%dolic_c(jc,jb)
            !IF(z_dolic>=MIN_DOLIC)THEN
            DO jk = 2, z_dolic
              !IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
              IF ( p_patch_3D%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
                !delta_z = p_os%p_diag%inv_prism_thick_c(jc,jk,jb)
                delta_z = p_patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,jk,jb)
                z_temp(jc,jk,jb)= trac_old(jc,jk,jb)&
                & -(delta_t*delta_z)*(flux_vert(jc,jk,jb)-flux_horz(jc,jk,jb))

              ENDIF
            ENDDO
        END DO
      END DO
      CALL sync_patch_array(SYNC_C, p_patch, z_temp)

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('BefImplDiff: trac_old',trac_old, str_module,idt_src)
      CALL dbg_print('BefImplDiff: z_temp'  ,z_temp  , str_module,idt_src)
      !---------------------------------------------------------------------
       !DO jk=1,n_zlev
       !write(*,*)'BEFORE DIFF: TRACER old:new',jk,&
       !&minval(trac_old(1:nproma,jk,1:p_patch%nblks_c)),maxval(trac_old(1:nproma,jk,1:p_patch%nblks_c)),&
       !&minval(z_temp(1:nproma,jk,1:p_patch%nblks_c)),maxval(z_temp(1:nproma,jk,1:p_patch%nblks_c)),&
      !&minval(A_v(:,jk,:)),maxval(A_v(:,jk,:)) 
       !END DO

      IF (ltimer) CALL timer_start(timer_dif_vert)

      !calculate vert diffusion impicit: result is stored in trac_out
      CALL tracer_diffusion_vert_impl_hom( p_patch_3D,                                        &
                                         & z_temp(1:nproma,1:n_zlev,1:p_patch%nblks_c),       &
                                         & p_os%p_prog(nnew(1))%h(1:nproma,1:p_patch%nblks_c),&
                                         & A_v,                                               &
                                         & p_op_coeff,                                        &
                                         & trac_new(1:nproma,1:n_zlev,1:p_patch%nblks_c))

       !DO jk=1,n_zlev
       !write(*,*)'AFTER DIFF: TRACER old:new',jk,&
       !&minval(trac_old(1:nproma,jk,1:p_patch%nblks_c)),maxval(trac_old(1:nproma,jk,1:p_patch%nblks_c)),&
       !&minval(trac_new(1:nproma,jk,1:p_patch%nblks_c)),maxval(trac_new(1:nproma,jk,1:p_patch%nblks_c))
       !END DO


      CALL sync_patch_array(SYNC_C, p_patch, trac_new)
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('AftImplDiff: z_temp'       ,z_temp                   ,str_module,idt_src)
      CALL dbg_print('AftImplDiff: trac_new'     ,trac_new                 ,str_module,idt_src)
      !---------------------------------------------------------------------

    !vertival diffusion is calculated explicitely
    ELSEIF(expl_vertical_tracer_diff==0)THEN
      CALL finish("advect_diffuse_vertical", &
      &"xplicit vertical tracer mixing currently not supported -  change namlist option")
    ENDIF ! lvertical_diff_implicit
    IF (ltimer) CALL timer_stop(timer_dif_vert)
  ENDIF
    !DO jk=1,n_zlev
    ! write(*,*)'TRACER old:new',jk,&
    ! &minval(trac_old(1:nproma,jk,1:p_patch%nblks_c)),maxval(trac_old(1:nproma,jk,1:p_patch%nblks_c)),&
    ! &minval(trac_new(1:nproma,jk,1:p_patch%nblks_c)),maxval(trac_new(1:nproma,jk,1:p_patch%nblks_c))
    ! END DO
  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=2  ! output print level (1-5, fix)
  CALL dbg_print('aft. AdvIndivTrac: trac_old',trac_old                 ,str_module,idt_src)
  CALL dbg_print('aft. AdvIndivTrac: trac_new',trac_new                 ,str_module,idt_src)
  !---------------------------------------------------------------------

END SUBROUTINE advect_individual_tracer_ab

END MODULE mo_oce_tracer
