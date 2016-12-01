!>
!! Contains the implementation of the tracer transport routines for the ICON ocean model.
!! This comprises advection and diffusion in horizontal and vertical direction.
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/01)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!----------------------------
#include "omp_definitions.inc"
#include "icon_definitions.inc"
!----------------------------
MODULE mo_ocean_tracer
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  USE mo_impl_constants,            ONLY: sea_boundary, sea, min_dolic
  USE mo_math_constants,            ONLY: pi
  USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer,              &
    & threshold_min_t, threshold_max_t, threshold_min_s, threshold_max_s, &
    & type_3dimrelax_temp, para_3dimrelax_temp,                           &
    & type_3dimrelax_salt, para_3dimrelax_salt,                           &
    & iswm_oce, l_edge_based,             &
    & flux_calculation_horz, flux_calculation_vert, miura_order1,         &
    & l_with_vert_tracer_diffusion, l_with_vert_tracer_advection,         &
    & l_skip_tracer,                                       &! , use_ThermoExpansion_Correction
    & GMRedi_configuration,GMRedi_combined,  GM_only,Redi_only ,          &
    & Cartesian_Mixing, tracer_threshold_min, tracer_threshold_max,       &
    & namelist_tracer_name, nbgcadv,                                      &
    & GMREDI_COMBINED_DIAGNOSTIC,GM_INDIVIDUAL_DIAGNOSTIC,REDI_INDIVIDUAL_DIAGNOSTIC
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer, debug_check_level
  USE mo_ocean_types,               ONLY: t_hydro_ocean_state, t_ocean_tracer !, v_base
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish !, message_text, message
  USE mo_ocean_boundcond,           ONLY: top_bound_cond_tracer
  USE mo_ocean_physics_types,       ONLY: t_ho_params
  USE mo_sea_ice_types,             ONLY: t_sfc_flx
  USE mo_ocean_diffusion,             ONLY: tracer_diffusion_vertical_implicit, tracer_diffusion_vert_explicit,tracer_diffusion_horz
  USE mo_ocean_tracer_transport_horz, ONLY: advect_horz, diffuse_horz
  USE mo_ocean_tracer_transport_vert, ONLY: advect_flux_vertical
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_c, sync_e, sync_patch_array
  USE mo_timer,                     ONLY: timer_start, timer_stop, timers_level, timer_dif_vert, timer_extra30
  USE mo_statistics,                ONLY: global_minmaxmean, print_value_location
  USE mo_mpi,                       ONLY: my_process_is_stdio !global_mpi_barrier
  USE mo_ocean_GM_Redi,             ONLY: calc_ocean_physics, prepare_ocean_physics
  USE mo_ocean_math_operators,      ONLY: div_oce_3d, verticalDiv_scalar_onFullLevels! !verticalDiv_scalar_midlevel
  USE mo_scalar_product,            ONLY: map_edges2edges_viacell_3d_const_z
  USE mo_physical_constants,        ONLY: clw, rho_ref
  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=12)           :: str_module = 'oceTracer   '  ! Output of module for 1 line debug
  INTEGER :: idt_src    = 1               ! Level of detail for 1 line debug

  PUBLIC  :: advect_ocean_tracers
  PUBLIC  :: advect_diffuse_tracer
  PRIVATE :: advect_diffuse_individual_tracer

CONTAINS
  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE advects the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE advect_ocean_tracers(patch_3d, p_os, p_param, p_sfc_flx,p_op_coeff, timestep)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)      :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: p_os
    TYPE(t_ho_params),                 INTENT(inout) :: p_param
    TYPE(t_sfc_flx),                   INTENT(inout) :: p_sfc_flx
    TYPE(t_operator_coeff),            INTENT(inout) :: p_op_coeff
    INTEGER :: timestep

    !Local variables
    INTEGER :: tracer_index, TracerDiffusion_coeff_index
    INTEGER :: start_cell_index, end_cell_index, jc, jb, level
    REAL(wp) :: z_relax!, delta_z
    INTEGER :: iloc(2)
    REAL(wp) :: zlat, zlon
    REAL(wp) :: z_c(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: minmaxmean(3)
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !-------------------------------------------------------------------------------
    patch_2D => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    start_detail_timer(timer_extra30,6)

    !alculate some information that is used for all tracers
    CALL prepare_tracer_transport( patch_3d, &
      & p_os,       &
      & p_op_coeff)
    !calculation of isopycnical slopes and tapering
    IF(GMRedi_configuration/=Cartesian_Mixing)THEN

      CALL prepare_ocean_physics(patch_3d, &
        & p_os,    &
        & p_param, &
        & p_op_coeff)

    ENDIF

    IF( iswm_oce /= 1) THEN
      DO tracer_index = 1, no_tracer
        CALL top_bound_cond_tracer( patch_2D,            &
                                  & p_os,               &
                                  & tracer_index,       &
                                  & p_sfc_flx,          &
                                  & p_os%p_aux%bc_top_tracer)
      ENDDO
    ENDIF

    stop_detail_timer(timer_extra30,6)

    DO tracer_index = 1, no_tracer+nbgcadv
      TracerDiffusion_coeff_index = MIN(tracer_index, 2) ! use the salinity coeeficient for tracers > 2      
      CALL advect_diffuse_individual_tracer( patch_3d,         &
        & p_os%p_prog(nold(1))%ocean_tracers(tracer_index), &
        & p_os, p_op_coeff,                                 &
        & p_os%p_aux%bc_top_tracer(:,:,tracer_index),       &
        & p_os%p_aux%bc_bot_tracer(:,:,tracer_index),       &
        & p_param,                                          &
        & p_param%k_tracer_h(:,:,:,TracerDiffusion_coeff_index),          &
        & p_param%a_tracer_v(:,:,:,TracerDiffusion_coeff_index),          &
        & p_os%p_prog(nnew(1))%ocean_tracers(tracer_index), &
        & tracer_index )

    END DO

    ! Final step: 3-dim temperature relaxation
    !  - strict time constant, i.e. independent of layer thickness
    !  - additional forcing Term F_T = -1/tau(T-T*) [ K/s ]
    !    when using the sign convention
    !      dT/dt = Operators + F_T
    !    i.e. F_T <0 for  T-T* >0 (i.e. decreasing temperature if it is warmer than relaxation data)
    !  - discretized:
    !    tracer = tracer - 1/(para_3dimRelax_Temp[months]) * (tracer(1)-data_3dimRelax_Temp)
    IF (no_tracer>=1 .AND. type_3dimrelax_temp >0) THEN

      ! calculate relaxation term
      z_relax = 1.0_wp/(para_3dimrelax_temp*2.592e6_wp)
      p_os%p_aux%forc_3dimrelax_temp(:,:,:) = -z_relax * p_os%p_aux%relax_3dim_coefficient(:,:,:) &
        & * ( p_os%p_prog(nnew(1))%tracer(:,:,:,1) - p_os%p_aux%data_3dimrelax_temp(:,:,:))

      ! add relaxation term to new temperature
      p_os%p_prog(nnew(1))%tracer(:,:,:,1) = p_os%p_prog(nnew(1))%tracer(:,:,:,1) + &
        &                                    p_os%p_aux%forc_3dimRelax_Temp(:,:,:) * dtime

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('3d_relax: AdvTracT forc', p_os%p_aux%forc_3dimRelax_Temp, str_module,idt_src, in_subset=cells_in_domain)
      CALL dbg_print('3d_relax: AdvTracT data', p_os%p_aux%data_3dimRelax_Temp, str_module,idt_src, in_subset=cells_in_domain)
      idt_src=2  ! output print level (1-5, fix)
      z_c(:,:,:) =  p_os%p_prog(nnew(1))%tracer(:,:,:,1)
      CALL dbg_print('3d_relax: AdvTracT trac', z_c, str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

    END IF

    ! Final step: 3-dim salinity relaxation
    !  - additional forcing Term F_S = -1/tau(S-S*) [ psu/s ]
    !    when using the sign convention
    !      dS/dt = Operators + F_S
    !    i.e. F_S <0 for  S-S* >0 (i.e. decreasing salinity if it is larger than relaxation data)
    !  - discretized:
    !    tracer = tracer - 1/(para_3dimRelax_Temp[months]) * (tracer(1)-data_3dimRelax_Temp)
    IF (no_tracer==2 .AND. type_3dimrelax_salt >0) THEN

      ! calculate relaxation term
      z_relax = 1.0_wp/(para_3dimrelax_salt*2.592e6_wp)
      p_os%p_aux%forc_3dimrelax_salt(1:nproma,:,1:patch_2D%nblks_c) = -z_relax* &
        & ( p_os%p_prog(nnew(1))%tracer(1:nproma,:,1:patch_2D%nblks_c,2) -       &
        & p_os%p_aux%forc_3dimrelax_salt(1:nproma,:,1:patch_2D%nblks_c))

      ! add relaxation term to new salinity
      p_os%p_prog(nnew(1))%tracer(1:nproma,:,1:patch_2D%nblks_c,2) = &
        & p_os%p_prog(nnew(1))%tracer(1:nproma,:,1:patch_2D%nblks_c,2) + &
        & p_os%p_aux%forc_3dimrelax_salt(1:nproma,:,1:patch_2D%nblks_c) * dtime

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('3d_relax: AdvTracS forc'  ,p_os%p_aux%forc_3dimrelax_salt,str_module,idt_src, in_subset=cells_in_domain)
      CALL dbg_print('3d_relax: AdvTracS data'  ,p_os%p_aux%data_3dimrelax_salt,str_module,idt_src, in_subset=cells_in_domain)
      idt_src=2  ! output print level (1-5, fix)
      z_c(:,:,:) =  p_os%p_prog(nnew(1))%tracer(:,:,:,2)
      CALL dbg_print('3d_relax: AdvTracS trac'  ,z_c                           ,str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

    END IF
    !!Commented out because of NAG-compiler, PK
    !TODO review IF statements concerning tracer transport


    !! apply additional volume flux to surface endLevelation - add to h_new after tracer advection
    !IF (l_forc_freshw) THEN
    !  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    !    CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
    !    DO jc = start_cell_index, end_cell_index
    !      p_os%p_prog(nnew(1))%h(jc,jb) = p_os%p_prog(nnew(1))%h(jc,jb) + p_sfc_flx%forc_fwfx(jc,jb)*dtime
    !    END DO
    !  END DO
    !END IF
    !CALL dbg_print('aft. AdvTracer: h-new (fwf)',p_os%p_prog(nnew(1))%h   ,str_module,idt_src)

!DO tracer_index = 1, no_tracer
!CALL dbg_print('OLD:',  p_os%p_prog(nold(1))%ocean_tracers(tracer_index),str_module,idt_src)
!END DO
  END SUBROUTINE advect_ocean_tracers
  !-------------------------------------------------------------------------


    !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE advect_diffuse_individual_tracer(patch_3d, old_ocean_tracer,       &
    & p_os, p_op_coeff,                      &
    & bc_top_tracer, bc_bot_tracer, p_param, &
    & k_h, a_v,                              &
    & new_ocean_tracer, tracer_index)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: old_ocean_tracer
    TYPE(t_ocean_tracer), TARGET :: new_ocean_tracer

    TYPE(t_hydro_ocean_state), TARGET :: p_os
    TYPE(t_operator_coeff),INTENT(inout) :: p_op_coeff
    REAL(wp), INTENT(in)                 :: bc_top_tracer(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                 :: bc_bot_tracer(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_ho_params),INTENT(inout)      :: p_param
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)       !horizontal mixing coeff
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)       !vertical mixing coeff, in
    INTEGER,  INTENT(in)                 :: tracer_index
!     REAL(wp), INTENT(inout), OPTIONAL :: horizontally_diffused_tracer(:,:,:)

    !Local variables
    REAL(wp) :: delta_t, delta_z,delta_z_new, delta_z1,delta_z_new1
    REAL(wp) :: div_adv_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flx(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    REAL(wp), POINTER :: trac_old(:,:,:), trac_new(:,:,:) ! temporary pointers to the concentration arrays

    INTEGER :: jc,level,jb,ic,ib
    INTEGER :: z_dolic
    INTEGER :: start_cell_index, end_cell_index
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_advection:advect_diffuse_tracer')
    !-------------------------------------------------------------------------------_
    trac_old => old_ocean_tracer%concentration
    trac_new => new_ocean_tracer%concentration

    patch_2D => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    delta_t = dtime

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('on entry: IndTrac: trac_old',trac_old(:,:,:) ,str_module,idt_src, in_subset=patch_2D%cells%owned)
    !---------------------------------------------------------------------
    IF (l_skip_tracer) THEN

      !   trac_new(1:nproma,1:n_zlev,1:patch_2D%nblks_c) = trac_old(1:nproma,1:n_zlev,1:patch_2D%nblks_c)
      new_ocean_tracer%concentration(1:nproma,1:n_zlev,1:patch_2D%nblks_c) = &
        & old_ocean_tracer%concentration(1:nproma,1:n_zlev,1:patch_2D%nblks_c)

      RETURN
    ENDIF
    !---------------------------------------------------------------------
    
    !Shallow water is done with horizontal advection
    IF(iswm_oce == 1) THEN

      div_adv_flux_horz   (1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks) = 0.0_wp
      div_adv_flux_vert   (1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks) = 0.0_wp      
      div_diff_flux_horz  (1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks) = 0.0_wp

      !---------------------------------------------------------------------
      CALL advect_horz( patch_3d,         &
        & old_ocean_tracer%concentration, &
        & p_os,                           &
        & p_op_coeff,                     &
        & k_h,                            &
        & p_os%p_prog(nold(1))%h,         &
        & p_os%p_prog(nnew(1))%h,         &
        & div_adv_flux_horz,              &
        & div_adv_flux_vert,              &
        & tracer_index)

      CALL diffuse_horz( patch_3d,         &
        & old_ocean_tracer%concentration, &
        & p_os,                           &
        & p_op_coeff,                     &
        & k_h,                            &
        & p_os%p_prog(nold(1))%h,         &
        & p_os%p_prog(nnew(1))%h,         &
        & div_diff_flux_horz)


      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('DivAdvFlx',div_adv_flux_horz,str_module,idt_src, in_subset=cells_in_domain)
      CALL dbg_print('DivDiffFlx',div_diff_flux_horz,str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

      !level=1
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          DO level = 1, 1!MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1

            !delta_z =  patch_3d%p_patch_1d(1)%del_zlev_m(1)
            delta_z = p_os%p_prog(nnew(1))%h(jc,jb) !- p_ext_data%oce%bathymetry_c(jc,jb)
            new_ocean_tracer%concentration(jc,level,jb)= old_ocean_tracer%concentration(jc,level,jb) - &
              & (delta_t/delta_z) * (div_adv_flux_horz(jc,level,jb)-div_diff_flux_horz(jc,level,jb))

          END DO
        END DO
      END DO

      CALL sync_patch_array(sync_c, patch_2D, new_ocean_tracer%concentration)
        
    !The 3D-case
    ELSE ! IF( iswm_oce /= 1) THEN

         CALL advect_diffuse_tracer( patch_3d, &
           & old_ocean_tracer,&
           & p_os,            &
           & p_op_coeff,      &
           & bc_top_tracer,   &
           & bc_bot_tracer,   &
           & p_param,         &
           & k_h,             &
           & a_v,             &
           & new_ocean_tracer,&
           & tracer_index)

    ENDIF

!    !---------DEBUG DIAGNOSTICS-------------------------------------------
!    CALL dbg_print('aft. AdvIndivTrac: trac_old', trac_old, str_module, 3, in_subset=cells_in_domain)
!    CALL dbg_print('aft. AdvIndivTrac: trac_new', trac_new, str_module, 1, in_subset=cells_in_domain)
!    !---------------------------------------------------------------------
! DO level = 1, 4
! write(0,*)'Trac new:old',level, maxval(trac_new(:,level,:)) ,minval(trac_new(:,level,:)),&
! &maxval(trac_old(:,level,:)) ,minval(trac_old(:,level,:))
! END DO

  END SUBROUTINE advect_diffuse_individual_tracer
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE advect_diffuse_tracer(patch_3d, old_ocean_tracer,       &
    & p_os, p_op_coeff,                    &
    & bc_top_tracer, bc_bot_tracer,p_param,&
    & k_h, a_v,                            &
    & new_ocean_tracer, tracer_index)!,        &
    ! & horizontally_diffused_tracer        )

    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: old_ocean_tracer
    TYPE(t_ocean_tracer), TARGET :: new_ocean_tracer

    TYPE(t_hydro_ocean_state), TARGET :: p_os
    TYPE(t_operator_coeff),INTENT(inout) :: p_op_coeff
    REAL(wp), INTENT(in)                 :: bc_top_tracer(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                 :: bc_bot_tracer(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_ho_params), INTENT(inout)     :: p_param
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)       !horizontal mixing coeff
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)       !vertical mixing coeff, in
    INTEGER,  INTENT(in)                 :: tracer_index
!     REAL(wp), INTENT(inout), OPTIONAL :: horizontally_diffused_tracer(:,:,:)

    !Local variables
    REAL(wp) :: delta_t, delta_z,delta_z_new, delta_z1,delta_z_new1
    REAL(wp) :: div_adv_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: flux_horz(nproma,n_zlev, patch_3d%p_patch_2D(1)%nblks_e)
    REAL(wp) :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flx_vert(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)        
    REAL(wp), POINTER :: trac_old(:,:,:), trac_new(:,:,:) ! temporary pointers to the concentration arrays
    TYPE(t_ocean_tracer) :: temp_tracer_before!(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)       
    TYPE(t_ocean_tracer) :: temp_tracer_after!(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)           
    INTEGER :: jc,level,jb, je
    INTEGER :: z_dolic
    INTEGER :: start_cell_index, end_cell_index
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D

    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_advection:advect_diffuse_tracer')
    !-------------------------------------------------------------------------------
    trac_old => old_ocean_tracer%concentration
    trac_new => new_ocean_tracer%concentration

    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain
    delta_t = dtime

!     div_adv_flux_horz (1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
!     div_diff_flux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
!     flux_horz(1:nproma,1:n_zlev,1: 1:patch_3d%p_patch_2d(1)%nblks_e)=0.0_wp
    !---------------------------------------------------------------------
 

    !---------------------------------------------------------------------
    IF ( l_with_vert_tracer_advection ) THEN

      CALL advect_flux_vertical( patch_3d,&
        & old_ocean_tracer%concentration, &
        & p_os,                           &
        & p_op_coeff,                     &
        & div_adv_flux_vert,              &
        & tracer_index)

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('aft. AdvFluxVert:divfluxvert',div_adv_flux_vert          ,str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

    ELSE
      div_adv_flux_vert(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks) = 0.0_wp
    ENDIF  ! l_with_vert_tracer_advection

    !---------------------------------------------------------------------
    CALL advect_horz( patch_3d,       &
      & old_ocean_tracer%concentration, &
      & p_os,                           &
      & p_op_coeff,                     &
      & k_h,                            &
      & p_os%p_prog(nold(1))%h,         &
      & p_os%p_prog(nnew(1))%h,         &
      & div_adv_flux_horz,              &
      & div_adv_flux_vert,              &
      & tracer_index )
    !---------------------------------------------------------------------

    IF(GMRedi_configuration==Cartesian_Mixing)THEN
      !horizontal diffusion, vertical is handled implicitely below
      CALL diffuse_horz( patch_3d,      &
      & old_ocean_tracer%concentration, &
      & p_os,                           &
      & p_op_coeff,                     &
      & k_h,                            &
      & p_os%p_prog(nold(1))%h,         &
      & p_os%p_prog(nnew(1))%h,         &
      & div_diff_flux_horz)

      div_diff_flx_vert(:,:,:) = 0.0_wp

    ELSEIF(GMRedi_configuration/=Cartesian_Mixing)THEN
    
      !calculate horizontal and vertical Redi and GM fluxes
      CALL calc_ocean_physics(patch_3d, p_os, p_param,p_op_coeff, tracer_index)
    
      !calculate horizontal divergence of diffusive flux
      div_diff_flux_horz(:,:,:)=0.0_wp
      CALL div_oce_3d( p_os%p_diag%GMRedi_flux_horz(:,:,:,tracer_index),&
                   &   patch_3d, &
                   &   p_op_coeff%div_coeff, &
                   &   div_diff_flux_horz )
                                     
      !vertical div of explicit part of vertical GMRedi-flux
      div_diff_flx_vert(:,:,:) = 0.0_wp
      CALL verticalDiv_scalar_onFullLevels( patch_3d, &
        & p_os%p_diag%GMRedi_flux_vert(:,:,:,tracer_index), &
        & div_diff_flx_vert)
                   
    END IF   
    !
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=1  ! output print level (1-5, fix)
    !CALL dbg_print('AftGMRedi: GMRediflux_h',p_os%p_diag%GMRedi_flux_horz(:,:,:,tracer_index),&
    !&str_module,idt_src, in_subset=edges_in_domain)
    Do level=1,n_zlev
      !Keep in mind that the div below lacks the division by prism thickness (cf. the equation
      !for the tracer update below),
      !
      CALL dbg_print('AftGMRedi: divGMRediflux_h',div_diff_flux_horz(:,level,:),&
      &str_module,idt_src, in_subset=cells_in_domain)
    END DO
    !This is only non-zero for GMRedi switched on
    IF(GMRedi_configuration/=Cartesian_Mixing)THEN  
      Do level=1,n_zlev
        CALL dbg_print('AftGMRedi: divGMRediflux_v',div_diff_flx_vert(:,level,:),&
        & str_module, idt_src, in_subset=cells_in_domain)      
      END DO
    ENDIF   
   !---------------------------------------------------------------------
    
                   
    IF(GMREDI_COMBINED_DIAGNOSTIC)THEN
    
      IF(tracer_index == 1) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, &
!ICON_OMP level) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index
      
            !top level
            level=1
            delta_z     = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+p_os%p_prog(nold(1))%h(jc,jb)
            delta_z_new = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+p_os%p_prog(nnew(1))%h(jc,jb)

            !This contains the divergence of the diffusive temperature fluxes,
            !either due to the GMRedi-scheme or due to the cartesian mixing  
            !The implicit contribution is not included.               
            p_os%p_diag%opottempGMRedi(jc,level,jb)&
            &=(div_diff_flux_horz(jc,level,jb)+div_diff_flx_vert(jc,level,jb))/delta_z_new
            !& * clw *rho_ref

            !This contains sum of advective and diffusive fluxes, i.e the whole tendency,
            !except for the implicit contribution.
            !This is overwritten each time a new tracer is calculated
            p_os%p_diag%div_of_GMRedi_flux(jc,level,jb)&
            &=-(div_adv_flux_horz(jc,level,jb) +div_adv_flux_vert(jc,level,jb)&
            & - div_diff_flux_horz(jc,level,jb)-div_diff_flx_vert(jc,level,jb))/delta_z_new
            
            DO level = 2, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

              !This contains the divergence of the diffusive temperature fluxes,
              !either due to the GMRedi-scheme or due to the cartesian mixing  
              !The implicit contribution is not included.   
              p_os%p_diag%opottempGMRedi(jc,level,jb)&
              &=(div_diff_flux_horz(jc,level,jb)+div_diff_flx_vert(jc,level,jb))/delta_z_new
              !& * clw *rho_ref
        
             !This contains sum of advective and diffusive fluxes, i.e the whole tendency,
             !except for the implicit contribution.
             !This is overwritten each time a new tracer is calculated
              p_os%p_diag%div_of_GMRedi_flux(jc,level,jb)&
              &=-(div_adv_flux_horz(jc,level,jb) +div_adv_flux_vert(jc,level,jb)&
              & - div_diff_flux_horz(jc,level,jb)-div_diff_flx_vert(jc,level,jb))/delta_z

            ENDDO
          END DO
        END DO
        
        DO level=1,n_zlev
          !CALL dbg_print('AftGMRedi: divGMRediflux',p_os%p_diag%div_of_GMRedi_flux(:,level,:),&
          !&str_module,idt_src, in_subset=cells_in_domain)
          CALL dbg_print('AftGMRedi: opottempGMRedi',p_os%p_diag%opottempGMRedi(:,level,:),&
          & str_module, idt_src, in_subset=cells_in_domain)      
        END DO
        
!ICON_OMP_END_PARALLEL_DO
     
      ELSEIF(tracer_index == 2) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, &
!ICON_OMP level) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index
      
            !top level
            level=1
            delta_z     = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+p_os%p_prog(nold(1))%h(jc,jb)
            delta_z_new = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+p_os%p_prog(nnew(1))%h(jc,jb)
              
            !This contains the divergence of the diffusive temperature fluxes,
            !either due to the GMRedi-scheme or due to the cartesian mixing   
            !The implicit contribution is not included.                         
            p_os%p_diag%osaltGMRedi(jc,level,jb)&
            &=-(div_adv_flux_horz(jc,level,jb) +div_adv_flux_vert(jc,level,jb)&
            &  - div_diff_flux_horz(jc,level,jb)-div_diff_flx_vert(jc,level,jb))/delta_z_new
                     
            !This contains sum of advective and diffusive fluxes, i.e the whole tendency,
            !except for the implicit contribution.
            !This is overwritten each time a new tracer is calculated
            p_os%p_diag%div_of_GMRedi_flux(jc,level,jb)&
            &=(div_diff_flux_horz(jc,level,jb)+div_diff_flx_vert(jc,level,jb))/delta_z_new

        
            DO level = 2, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
        
              !This contains the divergence of the diffusive temperature fluxes,
              !either due to the GMRedi-scheme or due to the cartesian mixing.
              !The implicit contribution is not included.   
              p_os%p_diag%osaltGMRedi(jc,level,jb)&
              &=(div_diff_flux_horz(jc,level,jb)+div_diff_flx_vert(jc,level,jb))/delta_z

             !This contains sum of advective and diffusive fluxes, i.e the whole tendency,
             !except for the implicit contribution.
             !This is overwritten each time a new tracer is calculated
              p_os%p_diag%div_of_GMRedi_flux(jc,level,jb)&
              &=-(div_adv_flux_horz(jc,level,jb) +div_adv_flux_vert(jc,level,jb)&
              & - div_diff_flux_horz(jc,level,jb)-div_diff_flx_vert(jc,level,jb))/delta_z

            ENDDO
          END DO
        END DO
        
        DO level=1,n_zlev
          !CALL dbg_print('AftGMRedi: divGMRediflux_h',div_diff_flux_horz(:,level,:),&
          !&str_module,idt_src, in_subset=cells_in_domain)
          CALL dbg_print('AftGMRedi: opotsaltGMRedi',p_os%p_diag%osaltGMRedi(:,level,:),&
          & str_module, idt_src, in_subset=cells_in_domain)      
        END DO
        
!ICON_OMP_END_PARALLEL_DO
      ENDIF!(tracer_index == 1)
    ENDIF!(GMREDI_COMBINED_DIAGNOSTIC)THEN  
      
      
    !Case: Implicit Vertical diffusion
    start_timer(timer_dif_vert,4)

    !Calculate preliminary tracer value out of horizontal advective and
    !diffusive fluxes and vertical advective fluxes, plus surface forcing.
    !Surface forcing applied as volume forcing at rhs, i.e.part of explicit term
    !in tracer (and also momentum) eqs. In this case, top boundary condition of
    !vertical Laplacians are homogeneous

!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, level, &
!ICON_OMP delta_z, delta_z_new) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        !TODO check algorithm: inv_prism_thick_c vs. del_zlev_m | * vs. /
        DO level = 1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1

          delta_z     = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+p_os%p_prog(nold(1))%h(jc,jb)
          delta_z_new = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+p_os%p_prog(nnew(1))%h(jc,jb)


          new_ocean_tracer%concentration(jc,level,jb)= &
            & (old_ocean_tracer%concentration(jc,level,jb) * delta_z &
            & - delta_t * (&
            &  div_adv_flux_horz(jc,level,jb) +div_adv_flux_vert(jc,level,jb)&
            & -div_diff_flux_horz(jc,level,jb)-div_diff_flx_vert(jc,level,jb))) / delta_z_new

          new_ocean_tracer%concentration(jc,level,jb) =         &
            & ( new_ocean_tracer%concentration(jc,level,jb) +   &
            & (delta_t  / delta_z_new) * bc_top_tracer(jc,jb))
            
            
          !p_os%p_diag%osaltGMRedi(jc,level,jb)&
          !&= (new_ocean_tracer%concentration(jc,level,jb)- old_ocean_tracer%concentration(jc,level,jb))*delta_z/delta_t
          
          !p_os%p_diag%opottempGMRedi(jc,level,jb)=-(div_adv_flux_horz(jc,level,jb) +div_adv_flux_vert(jc,level,jb)&
          !  &  - div_diff_flux_horz(jc,level,jb)-div_diff_flx_vert(jc,level,jb))
          


        ENDDO

        DO level = 2, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

          new_ocean_tracer%concentration(jc,level,jb) =                          &
            &  old_ocean_tracer%concentration(jc,level,jb) -                     &
            &  (delta_t /  patch_3d%p_patch_1D(1)%prism_thick_c(jc,level,jb))    &
            & * (div_adv_flux_horz(jc,level,jb) +div_adv_flux_vert(jc,level,jb)&
            &  - div_diff_flux_horz(jc,level,jb)-div_diff_flx_vert(jc,level,jb))

          !p_os%p_diag%osaltGMRedi(jc,level,jb)&
          !&= (new_ocean_tracer%concentration(jc,level,jb)- old_ocean_tracer%concentration(jc,level,jb))*&
          !&patch_3d%p_patch_1D(1)%prism_thick_c(jc,level,jb)/delta_t
                   
          !p_os%p_diag%opottempGMRedi(jc,level,jb)=-(div_adv_flux_horz(jc,level,jb) +div_adv_flux_vert(jc,level,jb)&
          !  &  - div_diff_flux_horz(jc,level,jb)-div_diff_flx_vert(jc,level,jb))

        ENDDO

      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO


    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('BefImplDiff: div_adv_flux_vert',div_adv_flux_vert, str_module,idt_src, in_subset=cells_in_domain)
    CALL dbg_print('BefImplDiff: trac_inter', new_ocean_tracer%concentration, str_module,idt_src, in_subset=cells_in_domain)
    !---------------------------------------------------------------------

    !calculate vert diffusion impicit: result is stored in trac_out
    ! no sync because of columnwise computation
    IF ( l_with_vert_tracer_diffusion ) THEN
    
      IF(GMREDI_COMBINED_DIAGNOSTIC)THEN
        ALLOCATE(temp_tracer_before%concentration(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks))
        ALLOCATE(temp_tracer_after%concentration(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks))      
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, &
!ICON_OMP level) ICON_OMP_DEFAULT_SCHEDULE
 
        !Store new tracer concentration in two arrays
        DO jb = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index
            DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)     

              temp_tracer_before%concentration(jc,level,jb)=new_ocean_tracer%concentration(jc,level,jb)            
              temp_tracer_after%concentration(jc,level,jb) =new_ocean_tracer%concentration(jc,level,jb)

            END DO
          END DO
        ENDDO 
!ICON_OMP_END_PARALLEL_DO         
      ENDIF!IF(GMREDI_COMBINED_DIAGNOSTIC)
      
      !Vertical mixing: implicit and with coefficient a_v
      !that is the sum of PP-coeff and implicit part of Redi-scheme
!DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)           
!write(0,*)'vert coeff',level,maxval(a_v(:,level,:)),minval(a_v(:,level,:))
!END DO
      CALL tracer_diffusion_vertical_implicit( &
          & patch_3d,                        &
          & new_ocean_tracer,                &
          & a_v,                             &
          & p_op_coeff)
          
      IF(GMREDI_COMBINED_DIAGNOSTIC)THEN
      
        !This is vertical mixing with of tracer with implicit part of GMRedi only.
        !It repeats the application of implicit vertical diffusion for diagnsotic reasons. 
        CALL tracer_diffusion_vertical_implicit(     &
          & patch_3d,                                &
          & temp_tracer_before,                      &
          & p_os%p_diag%vertical_mixing_coeff_GMRedi_implicit,&
          & p_op_coeff)
               
      
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, &
!ICON_OMP level ) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index

            DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)       

              p_os%p_diag%div_of_GMRedi_flux(jc,level,jb)&
              &=p_os%p_diag%div_of_GMRedi_flux(jc,level,jb)&
              &+(temp_tracer_after%concentration(jc,level,jb)-temp_tracer_before%concentration(jc,level,jb))/dtime            
            END DO
          END DO
        ENDDO 
!ICON_OMP_END_PARALLEL_DO    

        IF(tracer_index == 1) THEN     
          CALL dbg_print('AftGMRedi: temp. complete divofGMRediflux',p_os%p_diag%div_of_GMRedi_flux(:,:,:),&
          & str_module, idt_src, in_subset=cells_in_domain) 
        ELSEIF(tracer_index == 2) THEN
          CALL dbg_print('AftGMRedi: sal complete divofGMRediflux',p_os%p_diag%div_of_GMRedi_flux(:,:,:),&
          & str_module, idt_src, in_subset=cells_in_domain)        
        ENDIF     
        DEALLOCATE(temp_tracer_before%concentration)
        DEALLOCATE(temp_tracer_after%concentration) 
        
      ENDIF!IF(GMREDI_COMBINED_DIAGNOSTIC)THEN          
      

      IF(tracer_index == 1) THEN           
      !ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, &
!ICON_OMP level ) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index

            DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)       

              p_os%p_diag%opottemptend(jc,level,jb)&
              &=(new_ocean_tracer%concentration(jc,level,jb)&
              &- old_ocean_tracer%concentration(jc,level,jb))/dtime            
            END DO
          END DO
        ENDDO 
!ICON_OMP_END_PARALLEL_DO 

      ELSEIF(tracer_index == 2) THEN 
     !ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, &
!ICON_OMP level ) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index

            DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)       

              p_os%p_diag%osalttend(jc,level,jb)&
              &=(new_ocean_tracer%concentration(jc,level,jb)&
              &- old_ocean_tracer%concentration(jc,level,jb))/dtime            
            END DO
          END DO
        ENDDO 
!ICON_OMP_END_PARALLEL_DO 

      ENDIF!IF(tracer_index == 1)   
           
          
    ENDIF!IF ( l_with_vert_tracer_diffusion )

    CALL sync_patch_array(sync_c, patch_2D, new_ocean_tracer%concentration)

    stop_timer(timer_dif_vert,4)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('AftImplDiff: trac_new', trac_new, str_module, idt_src, in_subset=cells_in_domain)
    !---------------------------------------------------------------------

    IF (debug_check_level > 1 .and. tracer_index < 3 ) &
      CALL check_min_max_tracer(info_text="After advect_diffuse_tracer", tracer=new_ocean_tracer%concentration,     &
        & min_tracer=tracer_threshold_min(tracer_index), max_tracer=tracer_threshold_max(tracer_index), &
        & tracer_name=namelist_tracer_name(tracer_index), in_subset=cells_in_domain)
  
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('aft. AdvIndivTrac: trac_old', trac_old, str_module, 3, in_subset=cells_in_domain)
    CALL dbg_print('aft. AdvIndivTrac: trac_new', trac_new, str_module, 3, in_subset=cells_in_domain)
    IF(tracer_index == 1) THEN
      DO level=1,n_zlev
        CALL dbg_print('after trac: temp chg', p_os%p_diag%opottemptend(:,level,:), str_module, 3, in_subset=cells_in_domain)
      END DO
     ELSEIF(tracer_index == 2) THEN
      DO level=1,n_zlev
        CALL dbg_print('after trac: salt chg', p_os%p_diag%osalttend(:,level,:), str_module, 3, in_subset=cells_in_domain)
      END DO
    ENDIF  
      
    !---------------------------------------------------------------------
  

  END SUBROUTINE advect_diffuse_tracer
  !-------------------------------------------------------------------------

 
 
  !-------------------------------------------------------------------------
  !>
  !!    SUBROUTINE prepares next tracer transport step. Currently needed in horizontal
  !!    flux-scheme "MIMETIC-Miura". Geometric quantities are updated according to
  !!    actual velocity. This information is required by MIURA-scheme and is identical
  !!    for all tracers.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2012).
!<Optimize:inUse>
  SUBROUTINE prepare_tracer_transport(patch_3d, p_os, p_op_coeff)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: p_os
    TYPE(t_operator_coeff),INTENT(inout) :: p_op_coeff
    !
    !Local variables
    INTEGER :: startLevel, fin_level
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: je, level, jb,jc         !< index of edge, vert level, block
    INTEGER :: edge_cell_index(2), edge_cell_block(2)
!     INTEGER :: edge_vert_index(2), edge_vert_block(2)
    INTEGER :: upwind_index
    REAL(wp) :: delta_z, half_time
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc
    TYPE(t_cartesian_coordinates):: flux_sum
    !-------------------------------------------------------------------------------
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells
    !-------------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    all_cells       => patch_2d%cells%all
    edges_in_domain => patch_2d%edges%in_domain

    startLevel = 1
    half_time = 0.5_wp * dtime

!ICON_OMP_PARALLEL
    ! This should be changed
    ! just moving data around should not take place
!ICON_OMP_DO SCHEDULE(static)
    DO jb = all_cells%start_block, all_cells%end_block
      p_os%p_diag%w_time_weighted(1:nproma, 1:n_zlev+1, jb) = &
        & p_os%p_diag%w(1:nproma, 1:n_zlev+1, jb)
    ENDDO
!ICON_OMP_END_DO
    !In case of shallow water we have to to this here, for 3D fluid its done within vertical velocity calculation
    IF(iswm_oce==1)THEN
     CALL map_edges2edges_viacell_3d_const_z( patch_3d, p_os%p_diag%vn_time_weighted, p_op_coeff, &
        & p_os%p_diag%mass_flx_e)
    ENDIF    
 
    ! p_diag%w is compouted in_domain cells
    ! CALL sync_patch_array(SYNC_C, patch_2d,p_os%p_diag%w_time_weighted )

    ! This is already synced on edges_in_domain !
    ! CALL sync_patch_array(SYNC_E, patch_2d,p_os%p_diag%vn_time_weighted )

    !IF( .NOT.l_edge_based .OR. flux_calculation_horz==MIMETIC_MIURA) THEN
    ! default is flux_calculation_horz = fct_horz
 !   IF( flux_calculation_horz == miura_order1 ) THEN
!ICON_OMP_DO PRIVATE(start_edge_index, end_edge_index, je, edge_cell_index, edge_cell_block, &
!ICON_OMP fin_level, level, upwind_index) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_edge_index, end_edge_index)
      DO je = start_edge_index, end_edge_index
        !Get indices of two adjacent cells
        edge_cell_index(1) = patch_2d%edges%cell_idx(je,jb,1)
        edge_cell_block(1) = patch_2d%edges%cell_blk(je,jb,1)
        edge_cell_index(2) = patch_2d%edges%cell_idx(je,jb,2)
        edge_cell_block(2) = patch_2d%edges%cell_blk(je,jb,2)
!         edge_vert_index(1) = patch_2d%edges%vertex_idx(je,jb,1)
!         edge_vert_block(1) = patch_2d%edges%vertex_blk(je,jb,1)
!         edge_vert_index(2) = patch_2d%edges%vertex_idx(je,jb,2)
!         edge_vert_block(2) = patch_2d%edges%vertex_blk(je,jb,2)

        fin_level  = patch_3d%p_patch_1d(1)%dolic_e(je,jb)

        DO level = startLevel, fin_level
                        p_os%p_diag%p_vn_mean(je,level,jb)%x = 0.5_wp *                          &
                          & (p_os%p_diag%p_vn(edge_cell_index(1), level, edge_cell_block(1))%x + &
                          &  p_os%p_diag%p_vn(edge_cell_index(2), level, edge_cell_block(2))%x)
          !p_os%p_diag%p_vn_mean(je,level,jb)%x = 0.5_wp *                               &
          !  & (p_os%p_diag%p_vn_dual(edge_vert_index(1), level, edge_vert_block(1))%x + &
          !  & p_os%p_diag%p_vn_dual(edge_vert_index(2), level, edge_vert_block(2))%x)

          ! this is specific to the miura_order1_hflux_oce
          p_op_coeff%moved_edge_position_cc(je,level,jb)%x =      &
            & p_op_coeff%edge_position_cc(je,level,jb)%x          &
            & - half_time * p_os%p_diag%p_vn_mean(je,level,jb)%x

        END DO

        DO level = startLevel, fin_level
          upwind_index = MERGE(1, 2, p_os%p_diag%vn_time_weighted(je,level,jb) > 0.0_wp)

          p_op_coeff%upwind_cell_idx(je,level,jb) = edge_cell_index(upwind_index)
          p_op_coeff%upwind_cell_blk(je,level,jb) = edge_cell_block(upwind_index)

          p_op_coeff%upwind_cell_position_cc(je,level,jb)%x = &
            & patch_2d%cells%cartesian_center(edge_cell_index(upwind_index), edge_cell_block(upwind_index))%x
          ! & p_op_coeff%cell_position_cc(edge_cell_index(upwind_index), level, edge_cell_block(upwind_index))%x
        END DO

      END DO
    END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL

  !  ENDIF

  END SUBROUTINE prepare_tracer_transport
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE check_min_max_tracer(info_text, tracer, min_tracer, max_tracer, tracer_name, in_subset)
    CHARACTER(*) :: info_text
    REAL(wp), POINTER :: tracer(:,:,:)
    REAL(wp), INTENT(in) :: min_tracer, max_tracer
    CHARACTER(*) :: tracer_name
    TYPE(t_subset_range), POINTER :: in_subset

!     INTEGER  :: level
    REAL(wp) :: minmaxmean(3)
!     REAL(wp) :: lon, lat

      minmaxmean(:) = global_minmaxmean(values = tracer(:,:,:), in_subset=in_subset)
        IF (minmaxmean(1) < min_tracer) THEN
          WRITE(0,*) TRIM(tracer_name), ' too low:', minmaxmean(1)
          CALL print_value_location(tracer(:,:,:), minmaxmean(1), in_subset)
          CALL finish(TRIM(info_text), 'tracer below threshold')
        ENDIF

        IF (minmaxmean(2) > max_tracer) THEN
          WRITE(0,*) TRIM(tracer_name), ' too high:', minmaxmean(2)
          CALL print_value_location(tracer(:,:,:), minmaxmean(2), in_subset)
          CALL finish(TRIM(info_text), 'tracer above threshold')
        ENDIF

  END SUBROUTINE check_min_max_tracer
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  FUNCTION tracer_content(patch_3d, tracer, height) result(content)
    TYPE(t_patch_3d), TARGET, INTENT(in)    :: patch_3d
    REAL(wp), INTENT(in) :: tracer(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in) :: height(nproma,        patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) content
    REAL(wp) :: delta_z
    INTEGER :: jc,level,jb
    INTEGER :: z_dolic
    INTEGER :: start_cell_index, end_cell_index
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch
    !-------------------------------------------------------------------------------
    patch => patch_3d%p_patch_2d(1)
    cells_in_domain => patch%cells%in_domain
    !-------------------------------------------------------------------------------
    content = 0.0_wp
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          content = content + tracer(jc,level,jb)*patch%cells%area(jc,jb)*patch_3d%p_patch_1d(1)%prism_thick_c(jc,level,jb)
        END DO
      END DO
    END DO
  END FUNCTION tracer_content
  !-------------------------------------------------------------------------

END MODULE mo_ocean_tracer


