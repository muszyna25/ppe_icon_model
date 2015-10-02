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
    & vertical_tracer_diffusion_type, iswm_oce, l_edge_based,             &
    & flux_calculation_horz, flux_calculation_vert, miura_order1,         &
    & l_with_vert_tracer_diffusion, l_with_vert_tracer_advection,         &
    & tracer_update_mode, i_post_step, i_during_step,      &
    & l_skip_tracer, implicit_diffusion,explicit_diffusion,               &! , use_ThermoExpansion_Correction
    & GMRedi_configuration,GMRedi_combined,  GM_only,Redi_only ,Cartesian_Mixing
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer, debug_check_level
  USE mo_ocean_types,                 ONLY: t_hydro_ocean_state, t_ocean_tracer !, v_base
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish !, message_text, message
  USE mo_ocean_boundcond,             ONLY: top_bound_cond_tracer
  USE mo_ocean_physics
  USE mo_sea_ice_types,             ONLY: t_sfc_flx
  USE mo_ocean_diffusion,             ONLY: tracer_diffusion_vertical_implicit, tracer_diffusion_vert_explicit,tracer_diffusion_horz
  USE mo_ocean_tracer_transport_horz, ONLY: advect_horz, diffuse_horz
  USE mo_ocean_tracer_transport_vert, ONLY: advect_flux_vertical
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_c, sync_e, sync_patch_array
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_dif_vert, timer_extra30
  USE mo_statistics,                ONLY: global_minmaxmean
  USE mo_mpi,                       ONLY: my_process_is_stdio !global_mpi_barrier
  USE mo_ocean_GM_Redi,               ONLY: calc_ocean_physics, prepare_ocean_physics
  USE mo_ocean_math_operators,        ONLY: div_oce_3d, verticalDiv_scalar_midlevel
  USE mo_scalar_product,            ONLY: map_edges2edges_viacell_3d_const_z
  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=12)           :: str_module = 'oceTracer   '  ! Output of module for 1 line debug
  INTEGER :: idt_src    = 1               ! Level of detail for 1 line debug

  PUBLIC :: advect_tracer_ab

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
  SUBROUTINE advect_tracer_ab(p_patch_3d, p_os, p_param, p_sfc_flx,p_op_coeff, timestep)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)      :: p_patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: p_os
    TYPE(t_ho_params),                 INTENT(inout) :: p_param
    TYPE(t_sfc_flx),                   INTENT(inout) :: p_sfc_flx
    TYPE(t_operator_coeff),            INTENT(inout) :: p_op_coeff
    INTEGER :: timestep

    !Local variables
    INTEGER :: tracer_index, level
    INTEGER :: start_cell_index, end_cell_index, jc, jb
    REAL(wp) :: z_relax!, delta_z
    INTEGER :: iloc(2)
    REAL(wp) :: zlat, zlon
    REAL(wp) :: z_c(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: minmaxmean(3)
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: p_patch
    !-------------------------------------------------------------------------------
    p_patch => p_patch_3d%p_patch_2d(1)
    cells_in_domain => p_patch%cells%in_domain

    CALL timer_start(timer_extra30)

    !alculate some information that is used for all tracers
    CALL prepare_tracer_transport( p_patch_3d, &
      & p_os,       &
      & p_op_coeff)

    !calculation of isopycnical slopes and tapering
    IF(GMRedi_configuration/=Cartesian_Mixing)THEN

      CALL prepare_ocean_physics(p_patch_3d, &
        & p_os,    &
        & p_param, &
        & p_op_coeff)

    ENDIF

    IF( iswm_oce /= 1) THEN
      DO tracer_index = 1, no_tracer
        CALL top_bound_cond_tracer( p_patch,            &
                                  & p_os,               &
                                  & tracer_index,       &
                                  & p_sfc_flx,          &
                                  & p_os%p_aux%bc_top_tracer)
      ENDDO
    ENDIF

    CALL timer_stop(timer_extra30)

    DO tracer_index = 1, no_tracer

      CALL advect_individual_tracer_ab( p_patch_3d,         &
        & p_os%p_prog(nold(1))%ocean_tracers(tracer_index), &
        & p_os, p_op_coeff,                                 &
        & p_os%p_aux%bc_top_tracer(:,:,tracer_index),       &
        & p_os%p_aux%bc_bot_tracer(:,:,tracer_index),       &
        & p_param,                                          &
        & p_param%k_tracer_h(:,:,:,tracer_index ),          &
        & p_param%a_tracer_v(:,:,:, tracer_index),          &
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
      p_os%p_aux%forc_3dimrelax_salt(1:nproma,:,1:p_patch%nblks_c) = -z_relax* &
        & ( p_os%p_prog(nnew(1))%tracer(1:nproma,:,1:p_patch%nblks_c,2) -       &
        & p_os%p_aux%forc_3dimrelax_salt(1:nproma,:,1:p_patch%nblks_c))

      ! add relaxation term to new salinity
      p_os%p_prog(nnew(1))%tracer(1:nproma,:,1:p_patch%nblks_c,2) = &
        & p_os%p_prog(nnew(1))%tracer(1:nproma,:,1:p_patch%nblks_c,2) + &
        & p_os%p_aux%forc_3dimrelax_salt(1:nproma,:,1:p_patch%nblks_c) * dtime

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

    IF (debug_check_level > 5) THEN
      IF (no_tracer>=1) THEN
        DO level = 1, n_zlev


          minmaxmean(:) = global_minmaxmean(values = p_os%p_prog(nnew(1))%tracer(:,level,:,1), &
            & in_subset=cells_in_domain)
          IF (my_process_is_stdio()) THEN
            ! Abort if tracer is below or above threshold, read from namelist
            ! Temperature: <-1.9 deg C, may be possible, limit set to lower value
            IF (minmaxmean(1) < threshold_min_t) THEN
              WRITE(0,*) ' TEMPERATURE BELOW THRESHOLD = ', threshold_min_t
              WRITE(0,*) ' too low temperature at level =', level, minmaxmean(1)
              CALL finish(TRIM('mo_tracer_advection:advect_tracer'), &
                & 'Temperature below threshold')
            ENDIF

            ! Temperature: >100 deg C aborts
            IF (minmaxmean(2) > threshold_max_t) THEN
              WRITE(0,*) ' TEMPERATURE ABOVE THRESHOLD = ', threshold_max_t
              WRITE(0,*) ' too high temperature at level =', level, minmaxmean(2)
              CALL finish(TRIM('mo_tracer_advection:advect_tracer'), &
                & 'Temperature above threshold')
            ENDIF
          ENDIF

        END DO
      ENDIF

      IF (no_tracer>=2)THEN
        DO level = 1, n_zlev
          minmaxmean(:) = global_minmaxmean(values = p_os%p_prog(nnew(1))%tracer(:,level,:,2), &
            & in_subset=cells_in_domain)

          IF (my_process_is_stdio()) THEN
            ! Abort if salinity is negative:
            IF (minmaxmean(1) < threshold_min_s) THEN
              WRITE(0,*) ' SALINITY BELOW THRESHOLD = ', threshold_min_s
              WRITE(0,*) ' too low salinity at level =', level, minmaxmean(1)
              CALL finish(TRIM('mo_tracer_advection:advect_tracer'), &
                & 'SALINITY NEGATIVE')
            ENDIF
            ! Abort if salinity is >60 psu:
            IF (minmaxmean(2) > threshold_max_s) THEN
              WRITE(0,*) ' SALINITY ABOVE THRESHOLD = ', threshold_max_s
              WRITE(0,*) ' too high salinity at level =', level, minmaxmean(2)
              CALL finish(TRIM('mo_tracer_advection:advect_tracer'), &
                & 'TOO LARGE SALINITY')
            ENDIF
          ENDIF

        END DO
      ENDIF
    ENDIF !     IF (debug_check_level > 5)

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

  END SUBROUTINE advect_tracer_ab
  !-------------------------------------------------------------------------


    !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE advect_individual_tracer_ab(p_patch_3d, old_ocean_tracer,       &
    & p_os, p_op_coeff,                      &
    & bc_top_tracer, bc_bot_tracer, p_param, &
    & k_h, a_v,                              &
    & new_ocean_tracer, tracer_index)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: p_patch_3d
    TYPE(t_ocean_tracer), TARGET :: old_ocean_tracer
    TYPE(t_ocean_tracer), TARGET :: new_ocean_tracer

    TYPE(t_hydro_ocean_state), TARGET :: p_os
    TYPE(t_operator_coeff),INTENT(inout) :: p_op_coeff
    REAL(wp), INTENT(in)                 :: bc_top_tracer(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                 :: bc_bot_tracer(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_ho_params),INTENT(inout)      :: p_param
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)       !horizontal mixing coeff
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)       !vertical mixing coeff, in
    INTEGER,  INTENT(in)                 :: tracer_index
!     REAL(wp), INTENT(inout), OPTIONAL :: horizontally_diffused_tracer(:,:,:)

    !Local variables
    REAL(wp) :: delta_t, delta_z,delta_z_new, delta_z1,delta_z_new1
    REAL(wp) :: div_adv_flux_horz(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: flux_vert(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flx(nproma, n_zlev,p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    REAL(wp), POINTER :: trac_old(:,:,:), trac_new(:,:,:) ! temporary pointers to the concentration arrays

    INTEGER :: jc,level,jb,ic,ib
    INTEGER :: z_dolic
    INTEGER :: start_cell_index, end_cell_index
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: p_patch
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_advection:advect_individual_tracer')
    !-------------------------------------------------------------------------------_
    trac_old => old_ocean_tracer%concentration
    trac_new => new_ocean_tracer%concentration

    p_patch => p_patch_3d%p_patch_2d(1)
    cells_in_domain => p_patch%cells%in_domain
    delta_t = dtime



    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('on entry: IndTrac: trac_old',trac_old(:,:,:) ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------
    IF (l_skip_tracer) THEN

      !   trac_new(1:nproma,1:n_zlev,1:p_patch%nblks_c) = trac_old(1:nproma,1:n_zlev,1:p_patch%nblks_c)
      new_ocean_tracer%concentration(1:nproma,1:n_zlev,1:p_patch%nblks_c) = &
        & old_ocean_tracer%concentration(1:nproma,1:n_zlev,1:p_patch%nblks_c)

      RETURN
    ENDIF
    !---------------------------------------------------------------------
    
    !Shallow water is done with horizontal advection
    IF(iswm_oce == 1) THEN
      div_adv_flux_horz   (1:nproma,1:n_zlev,1:p_patch%alloc_cell_blocks) = 0.0_wp
      div_diff_flux_horz  (1:nproma,1:n_zlev,1:p_patch%alloc_cell_blocks) = 0.0_wp

      !---------------------------------------------------------------------
      CALL advect_horz( p_patch_3d,       &
        & old_ocean_tracer%concentration, &
        & p_os,                           &
        & p_op_coeff,                     &
        & k_h,                            &
        & p_os%p_prog(nold(1))%h,         &
        & p_os%p_prog(nnew(1))%h,         &
        & div_adv_flux_horz)

      CALL diffuse_horz( p_patch_3d,      &
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
          DO level = 1, 1!MIN(p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1

            !delta_z =  p_patch_3d%p_patch_1d(1)%del_zlev_m(1)
            delta_z = p_os%p_prog(nnew(1))%h(jc,jb) !- p_ext_data%oce%bathymetry_c(jc,jb)
            new_ocean_tracer%concentration(jc,level,jb)= old_ocean_tracer%concentration(jc,level,jb) - &
              & (delta_t/delta_z) * (div_adv_flux_horz(jc,level,jb)-div_diff_flux_horz(jc,level,jb))
!write(123,*)'data',new_ocean_tracer%concentration(jc,level,jb), old_ocean_tracer%concentration(jc,level,jb),&
!&delta_z,div_adv_flux_horz(jc,level,jb),div_diff_flux_horz(jc,level,jb),p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          END DO
        END DO
      END DO

      CALL sync_patch_array(sync_c, p_patch, new_ocean_tracer%concentration)
        
   !The 3D-case: first vertical fluxes than preliminary tracer value and
    !finally implicit vertical diffusion
    ELSE ! IF( iswm_oce /= 1) THEN
       IF(tracer_update_mode == i_post_step)THEN
         ! default setup
         CALL advect_individual_tracer_ab_post_step( &
           & p_patch_3d, old_ocean_tracer,           &
           & p_os, p_op_coeff,                       &
           & bc_top_tracer, bc_bot_tracer,p_param,   &
           & k_h, a_v,                               &
           & new_ocean_tracer, tracer_index)

       ELSEIF(tracer_update_mode == i_during_step)THEN

         CALL advect_individual_tracer_ab_during_step( &
           & p_patch_3d, old_ocean_tracer,         &
           & p_os, p_op_coeff,                     &
           & bc_top_tracer, bc_bot_tracer, p_param,&
           & k_h, a_v,                             &
           & new_ocean_tracer, tracer_index)

       ENDIF
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('aft. AdvIndivTrac: trac_old', trac_old, str_module, 3, in_subset=cells_in_domain)
    CALL dbg_print('aft. AdvIndivTrac: trac_new', trac_new, str_module, 1, in_subset=cells_in_domain)
    !---------------------------------------------------------------------
! DO level = 1, 4
! write(0,*)'Trac new:old',level, maxval(trac_new(:,level,:)) ,minval(trac_new(:,level,:)),&
! &maxval(trac_old(:,level,:)) ,minval(trac_old(:,level,:))
! END DO

  END SUBROUTINE advect_individual_tracer_ab
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE advect_individual_tracer_ab_post_step(p_patch_3d, old_ocean_tracer,       &
    & p_os, p_op_coeff,                    &
    & bc_top_tracer, bc_bot_tracer,p_param,&
    & k_h, a_v,                            &
    & new_ocean_tracer, tracer_index)!,        &
    ! & horizontally_diffused_tracer        )

    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: p_patch_3d
    TYPE(t_ocean_tracer), TARGET :: old_ocean_tracer
    TYPE(t_ocean_tracer), TARGET :: new_ocean_tracer

    TYPE(t_hydro_ocean_state), TARGET :: p_os
    TYPE(t_operator_coeff),INTENT(inout) :: p_op_coeff
    REAL(wp), INTENT(in)                 :: bc_top_tracer(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                 :: bc_bot_tracer(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_ho_params), INTENT(inout)     :: p_param
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)       !horizontal mixing coeff
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)       !vertical mixing coeff, in
    INTEGER,  INTENT(in)                 :: tracer_index
!     REAL(wp), INTENT(inout), OPTIONAL :: horizontally_diffused_tracer(:,:,:)

    !Local variables
    REAL(wp) :: delta_t, delta_z,delta_z_new, delta_z1,delta_z_new1
    REAL(wp) :: div_adv_flux_horz(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flux_horz2(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: flux_horz(nproma,n_zlev, p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: div_adv_flux_vert(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flx_vert(nproma, n_zlev,p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    REAL(wp), POINTER :: trac_old(:,:,:), trac_new(:,:,:) ! temporary pointers to the concentration arrays

    INTEGER :: jc,level,jb, je
    INTEGER :: z_dolic
    INTEGER :: start_cell_index, end_cell_index
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: p_patch

    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_advection:advect_individual_tracer')
    !-------------------------------------------------------------------------------
    trac_old => old_ocean_tracer%concentration
    trac_new => new_ocean_tracer%concentration

    p_patch         => p_patch_3d%p_patch_2d(1)
    cells_in_domain => p_patch%cells%in_domain
    edges_in_domain => p_patch%edges%in_domain
    delta_t = dtime

    div_adv_flux_horz (1:nproma,1:n_zlev,1:p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    div_diff_flux_horz(1:nproma,1:n_zlev,1:p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    div_diff_flx_vert (1:nproma,1:n_zlev,1:p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp

    div_diff_flux_horz2(1:nproma,1:n_zlev, 1:p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    flux_horz(1:nproma,1:n_zlev,1: 1:p_patch_3d%p_patch_2d(1)%nblks_e)=0.0_wp
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    CALL advect_horz( p_patch_3d,       &
      & old_ocean_tracer%concentration, &
      & p_os,                           &
      & p_op_coeff,                     &
      & k_h,                            &
      & p_os%p_prog(nold(1))%h,         &
      & p_os%p_prog(nnew(1))%h,         &
      & div_adv_flux_horz)
    !---------------------------------------------------------------------

    IF ( l_with_vert_tracer_advection ) THEN

      CALL advect_flux_vertical( p_patch_3d,                     &
        & old_ocean_tracer%concentration, &
        & p_os,                           &
        & p_op_coeff,                     &
!         & bc_top_tracer,                  &
!         & bc_bot_tracer,                  &
        & div_adv_flux_vert,                      &
        & tracer_index)

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('aft. AdvFluxVert:divfluxvert',div_adv_flux_vert          ,str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

    ELSE
      div_adv_flux_vert(1:nproma,1:n_zlev,1:p_patch%alloc_cell_blocks) = 0.0_wp
    ENDIF  ! l_with_vert_tracer_advection

    IF(GMRedi_configuration/=Cartesian_Mixing)THEN
      !calculate horizontal and vertical Redi and GM fluxes
      CALL calc_ocean_physics(p_patch_3d, p_os, p_param,p_op_coeff, tracer_index)
    ENDIF


    IF(GMRedi_configuration==Cartesian_Mixing)THEN
      !horizontal diffusion, vertical is handled implicitely below
      CALL diffuse_horz( p_patch_3d,      &
      & old_ocean_tracer%concentration, &
      & p_os,                           &
      & p_op_coeff,                     &
      & k_h,                            &
      & p_os%p_prog(nold(1))%h,         &
      & p_os%p_prog(nnew(1))%h,         &
      & div_diff_flux_horz)

    ELSE
      !horizontal
      CALL div_oce_3d( p_os%p_diag%GMRedi_flux_horz(:,:,:,tracer_index),&
                   &   p_patch_3D, &
                   &   p_op_coeff%div_coeff, &
                   &   div_diff_flux_horz )
      !vertical div of GMRedi-flux
      CALL verticalDiv_scalar_midlevel( p_patch_3d, &
                                      & p_os%p_diag%GMRedi_flux_vert(:,:,:,tracer_index), &
                                      & div_diff_flx_vert)
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('AftGMRedi: GMRediflux_h',p_os%p_diag%GMRedi_flux_horz(:,:,:,tracer_index),&
      &str_module,idt_src, in_subset=edges_in_domain)
      CALL dbg_print('AftGMRedi: GMRediflux_v',p_os%p_diag%GMRedi_flux_vert(:,:,:,tracer_index),&
      & str_module, idt_src, in_subset=cells_in_domain)
      CALL dbg_print('AftGMRedi: divGMRediflux_h',div_diff_flux_horz(:,:,:),&
      &str_module,idt_src, in_subset=cells_in_domain)
      CALL dbg_print('AftGMRedi: divGMRediflux_v',div_diff_flx_vert(:,:,:),&
      & str_module, idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

    END IF

    !Case: Implicit Vertical diffusion
    IF ( vertical_tracer_diffusion_type == implicit_diffusion ) THEN

      IF (ltimer) CALL timer_start(timer_dif_vert)

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
          DO level = 1, MIN(p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1

            delta_z     = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nold(1))%h(jc,jb)
            delta_z_new = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nnew(1))%h(jc,jb)
            !  not yet changed
            !  delta_z     = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+p_os%p_prog(nold(1))%h(jc,jb)
            !  delta_z_new = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+p_os%p_prog(nnew(1))%h(jc,jb)

            new_ocean_tracer%concentration(jc,level,jb)= &
              & (old_ocean_tracer%concentration(jc,level,jb) * delta_z &
              & - delta_t * (div_adv_flux_vert(jc,level,jb)-&
              &  (div_diff_flux_horz(jc,level,jb)-div_adv_flux_horz(jc,level,jb)))) / delta_z_new

            new_ocean_tracer%concentration(jc,level,jb) =         &
              & ( new_ocean_tracer%concentration(jc,level,jb) +   &
              & (delta_t  / delta_z_new) * bc_top_tracer(jc,jb))

          ENDDO

          DO level = 2, p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

            ! delta_z = p_patch_3d%p_patch_1d(1)%del_zlev_m(level)
            ! delta_z = p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb)

            new_ocean_tracer%concentration(jc,level,jb) =                          &
              &  old_ocean_tracer%concentration(jc,level,jb) -                     &
              &  (delta_t /  p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb))  &
              &    * (div_adv_flux_vert(jc,level,jb) - (div_diff_flux_horz(jc,level,jb)-div_adv_flux_horz(jc,level,jb)))
            !   test
            !   IF( delta_z/= delta_z1)THEN
            !     write(0,*)'no agreement',level,jc,jb,&
            !     &p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb), p_patch_3D%p_patch_1D(1)%del_zlev_m(level)
            !     &p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb), p_patch_3D%p_patch_1D(1)%del_zlev_m(level)
            !   ENDIF
          ENDDO
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('BefImplDiff: div_adv_flux_vert',div_adv_flux_vert, str_module,idt_src, in_subset=cells_in_domain)
      CALL dbg_print('BefImplDiff: trac_inter', new_ocean_tracer%concentration,  str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

      !calculate vert diffusion impicit: result is stored in trac_out
      ! no sync because of columnwise computation
      IF ( l_with_vert_tracer_diffusion ) THEN
        CALL tracer_diffusion_vertical_implicit( &
            & p_patch_3d,                      &
            & new_ocean_tracer,                &
            & a_v,                             &
            & p_op_coeff)
      ENDIF

      CALL sync_patch_array(sync_c, p_patch, new_ocean_tracer%concentration)

      IF (ltimer) CALL timer_stop(timer_dif_vert)

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('AftImplDiff: trac_new', trac_new, str_module, idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

            !vertival diffusion is calculated explicitely
    ELSEIF ( vertical_tracer_diffusion_type == explicit_diffusion ) THEN
      ! notInUse
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('BefExplDiff: trac_old', trac_old,  str_module,idt_src, in_subset=cells_in_domain)
      CALL dbg_print('BefExplDiff: div_adv_flux_vert',div_adv_flux_vert, str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

      SELECT CASE(GMRedi_configuration)

      CASE (Cartesian_Mixing)
        !this is the vertical Laplacian
        CALL tracer_diffusion_vert_explicit( p_patch_3d,            &
          & old_ocean_tracer%concentration, &
          & bc_top_tracer,                  &
          & a_v,                            &
          & div_diff_flx_vert)
!  DO level = 1, n_zlev
!    write(0,*)'After Cartesian div h/v', tracer_index,level,&
!    &maxval(div_diff_flux_horz(:,level,:)), minval(div_diff_flux_horz(:,level,:)),&
!    &maxval(div_diff_flx_vert(:,level,:)), minval(div_diff_flx_vert(:,level,:))
!  END DO

      END SELECT

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('AftExplDiff: div_diffflx_horz',div_diff_flux_horz,str_module,idt_src, in_subset=cells_in_domain)
      CALL dbg_print('AftExplDiff: div_diffflx_vert',div_diff_flx_vert, str_module, idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------


      ! top layer
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          DO level = 1, MIN(p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1
            ! NOTE: this will not wotk with partial cells on the top level
            delta_z     = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nold(1))%h(jc,jb)
            delta_z_new = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nnew(1))%h(jc,jb)

            new_ocean_tracer%concentration(jc,level,jb) = (old_ocean_tracer%concentration(jc,level,jb) * delta_z &
            & -delta_t*(div_adv_flux_vert(jc,level,jb)-(div_diff_flux_horz(jc,level,jb)-div_adv_flux_horz(jc,level,jb))&
            &-div_diff_flx_vert(jc,level,jb)))/delta_z_new

            new_ocean_tracer%concentration(jc,level,jb) = new_ocean_tracer%concentration(jc,level,jb) + &
            & (delta_t/delta_z) * bc_top_tracer(jc,jb)

          ENDDO

          DO level = 2, p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
            delta_z = p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb)

            new_ocean_tracer%concentration(jc,level,jb)= old_ocean_tracer%concentration(jc,level,jb)   &
            & - (delta_t/delta_z) * (div_adv_flux_vert(jc,level,jb)-&
            &(div_diff_flux_horz(jc,level,jb)-div_adv_flux_horz(jc,level,jb))&
            &-div_diff_flx_vert(jc,level,jb))

          ENDDO
        END DO
      END DO
      
      CALL sync_patch_array(sync_c, p_patch, new_ocean_tracer%concentration)
       
    ELSE!wrong tracer diffusion configuration
      CALL finish('TRIM(advect tracer)',"This wrong tracer diffusion configuration: neither explicit nor implicit")

    ENDIF ! implicit or explicit
!  Do level=1,n_zlev
!  write(0,*)'tracer old-new',level,&
!  &maxval(old_ocean_tracer%concentration(:,level,:)),minval(old_ocean_tracer%concentration(:,level,:)),&
!  &maxval(new_ocean_tracer%concentration(:,level,:)),minval(new_ocean_tracer%concentration(:,level,:))
!  End do
  END SUBROUTINE advect_individual_tracer_ab_post_step
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
  SUBROUTINE advect_individual_tracer_ab_during_step(p_patch_3d, old_ocean_tracer,       &
    & p_os, p_op_coeff,                    &
    & bc_top_tracer, bc_bot_tracer,p_param,&
    & k_h, a_v,                            &
    & new_ocean_tracer, tracer_index)!,        &
    ! & horizontally_diffused_tracer        )

    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: p_patch_3d
    TYPE(t_ocean_tracer), TARGET :: old_ocean_tracer
    TYPE(t_ocean_tracer), TARGET :: new_ocean_tracer

    TYPE(t_hydro_ocean_state), TARGET :: p_os
    TYPE(t_operator_coeff),INTENT(inout) :: p_op_coeff
    REAL(wp), INTENT(in)                 :: bc_top_tracer(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                 :: bc_bot_tracer(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_ho_params), INTENT(inout)     :: p_param
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)       !horizontal mixing coeff
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)       !vertical mixing coeff, in
    INTEGER,  INTENT(in)                 :: tracer_index
!     REAL(wp), INTENT(inout), OPTIONAL :: horizontally_diffused_tracer(:,:,:)

    !Local variables
    REAL(wp) :: delta_t, delta_z,delta_z_new, delta_z1,delta_z_new1
    REAL(wp) :: div_adv_flux_horz(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    REAL(wp) :: trac_old_plus_forc(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_adv_flux_vert(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flx(nproma, n_zlev,p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    REAL(wp), POINTER :: trac_old(:,:,:), trac_new(:,:,:) ! temporary pointers to the concentration arrays

    INTEGER :: jc,level,jb,ic,ib
    INTEGER :: z_dolic
    INTEGER :: start_cell_index, end_cell_index
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: p_patch

    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_advection:advect_individual_tracer')
    !-------------------------------------------------------------------------------_!

   !-------------------------------------------------------------------------------_!
    CALL finish(TRIM('mo_tracer_advection:advect_individual_tracer_ab_during_step'), 'This not tested and does probably not work !')
    !-------------------------------------------------------------------------------_!


    trac_old => old_ocean_tracer%concentration
    trac_new => new_ocean_tracer%concentration

    p_patch => p_patch_3d%p_patch_2d(1)
    cells_in_domain => p_patch%cells%in_domain
    delta_t = dtime

    div_adv_flux_horz(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) =0.0_wp
    div_diff_flux_horz(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    trac_old_plus_forc(nproma,n_zlev, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    !---------------------------------------------------------------------


    !update old tracer with surface forcing
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, level, &
!ICON_OMP delta_z, delta_z_new) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        !TODO check algorithm: inv_prism_thick_c vs. del_zlev_m | * vs. /
        DO level = 1, MIN(p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1

          delta_z     = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nold(1))%h(jc,jb)
          delta_z_new = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nnew(1))%h(jc,jb)
          !  not yet changed
          !  delta_z     = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+p_os%p_prog(nold(1))%h(jc,jb)
          !  delta_z_new = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+p_os%p_prog(nnew(1))%h(jc,jb)

          trac_old(jc,level,jb) = trac_old(jc,level,jb) +   &
            & (delta_t/delta_z_new ) * bc_top_tracer(jc,jb)

        ENDDO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

    !---------------------------------------------------------------------
    CALL diffuse_horz( p_patch_3d,      &
      & trac_old,                       &
      & p_os,                           &
      & p_op_coeff,                     &
      & k_h,                            &
      & p_os%p_prog(nold(1))%h,         &
      & p_os%p_prog(nnew(1))%h,         &
      & div_diff_flux_horz)!,                      &
      ! & horizontally_diffused_tracer )

!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, level, &
!ICON_OMP delta_z, delta_z_new) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index
            !TODO check algorithm: inv_prism_thick_c vs. del_zlev_m | * vs. /
            DO level = 1, MIN(p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1

              delta_z     = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nold(1))%h(jc,jb)
              delta_z_new = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nnew(1))%h(jc,jb)

              trac_old(jc,level,jb)= &
                & (trac_old(jc,level,jb) * delta_z         &
                & - delta_t * (- &
                &  (div_diff_flux_horz(jc,level,jb)))) / delta_z_new

            ENDDO

            DO level = 2, p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

              trac_old(jc,level,jb) =                          &
                &  trac_old(jc,level,jb) -                                           &
                &  (delta_t /  p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb))  &
                &    * (- div_diff_flux_horz(jc,level,jb))
            ENDDO
          END DO
        END DO
!ICON_OMP_END_PARALLEL_DO


    CALL advect_horz( p_patch_3d,       &
      & trac_old,                       &
      & p_os,                           &
      & p_op_coeff,                     &
      & k_h,                            &
      & p_os%p_prog(nold(1))%h,         &
      & p_os%p_prog(nnew(1))%h,         &
      & div_adv_flux_horz)!,                      &
      ! & horizontally_diffused_tracer )

!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, level, &
!ICON_OMP delta_z, delta_z_new) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index
            !TODO check algorithm: inv_prism_thick_c vs. del_zlev_m | * vs. /
            DO level = 1, MIN(p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1

              delta_z     = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nold(1))%h(jc,jb)
              delta_z_new = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nnew(1))%h(jc,jb)

              trac_old(jc,level,jb)= &
                & (trac_old(jc,level,jb) * delta_z         &
                & - delta_t * div_adv_flux_horz(jc,level,jb))/ delta_z_new

            ENDDO

            DO level = 2, p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

              trac_old(jc,level,jb) =                          &
                &  trac_old(jc,level,jb) -                                           &
                &  (delta_t /  p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb))  &
                &    * div_adv_flux_horz(jc,level,jb)
            ENDDO
          END DO
        END DO
!ICON_OMP_END_PARALLEL_DO


    !---------------------------------------------------------------------

      IF ( l_with_vert_tracer_advection ) THEN

        CALL advect_flux_vertical( p_patch_3d,&
          & trac_old,                         &
          & p_os,                             &
          & p_op_coeff,                       &
!           & bc_top_tracer,                    &
!           & bc_bot_tracer,                    &
          & div_adv_flux_vert,                    &
          & tracer_index)

        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=3  ! output print level (1-5, fix)
         CALL dbg_print('aft. AdvFluxVert:divfluxvert',div_adv_flux_vert          ,str_module,idt_src, in_subset=cells_in_domain)
        !---------------------------------------------------------------------

      ELSE
        div_adv_flux_vert(1:nproma,1:n_zlev,1:p_patch%alloc_cell_blocks) = 0.0_wp
      ENDIF  ! l_with_vert_tracer_advection


      !Case: Implicit Vertical diffusion
      IF(vertical_tracer_diffusion_type == implicit_diffusion)THEN
        IF (ltimer) CALL timer_start(timer_dif_vert)

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
            DO level = 1, MIN(p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1

              delta_z     = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nold(1))%h(jc,jb)
              delta_z_new = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nnew(1))%h(jc,jb)
              !  not yet changed
              !  delta_z     = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+p_os%p_prog(nold(1))%h(jc,jb)
              !  delta_z_new = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+p_os%p_prog(nnew(1))%h(jc,jb)

              new_ocean_tracer%concentration(jc,level,jb)= &
                & (trac_old(jc,level,jb) * delta_z         &
                & - delta_t * (div_adv_flux_vert(jc,level,jb))) / delta_z_new

            ENDDO

            DO level = 2, p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

              ! delta_z = p_patch_3d%p_patch_1d(1)%del_zlev_m(level)
              ! delta_z = p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb)

              new_ocean_tracer%concentration(jc,level,jb) =                          &
                &  trac_old(jc,level,jb) -                                           &
                &  (delta_t /  p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb))  &
                &    * (div_adv_flux_vert(jc,level,jb))
              !   test
              !   IF( delta_z/= delta_z1)THEN
              !     write(0,*)'no agreement',level,jc,jb,&
              !     &p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb), p_patch_3D%p_patch_1D(1)%del_zlev_m(level)
              !     &p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb), p_patch_3D%p_patch_1D(1)%del_zlev_m(level)
              !   ENDIF
            ENDDO
          END DO
        END DO
!ICON_OMP_END_PARALLEL_DO
        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=3  ! output print level (1-5, fix)
        CALL dbg_print('BefImplDiff: div_adv_flux_vert',div_adv_flux_vert, str_module,idt_src, in_subset=cells_in_domain)
         CALL dbg_print('BefImplDiff: trac_inter', new_ocean_tracer%concentration,  str_module,idt_src, in_subset=cells_in_domain)
        !---------------------------------------------------------------------
        !CALL finish(TRIM('mo_tracer_advection:before vertical diffusion'), 'test of thickness')

        !calculate vert diffusion impicit: result is stored in trac_out
        ! no sync because of columnwise computation
        IF ( l_with_vert_tracer_diffusion ) THEN
          CALL tracer_diffusion_vertical_implicit( &
            & p_patch_3d,                      &
            & new_ocean_tracer,                &
            & a_v,                             &
            & p_op_coeff)!,                      &
          !&                                  trac_new(:,:,:))
        ENDIF

        CALL sync_patch_array(sync_c, p_patch, new_ocean_tracer%concentration)

        IF (ltimer) CALL timer_stop(timer_dif_vert)

        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=3  ! output print level (1-5, fix)
        CALL dbg_print('AftImplDiff: trac_new', trac_new, str_module, idt_src, in_subset=cells_in_domain)
        !---------------------------------------------------------------------

        !vertival diffusion is calculated explicitely
      ELSE ! IF(vertical_tracer_diffusion_type == explicit)THEN
        ! notInUse
        IF ( l_with_vert_tracer_diffusion ) THEN

          !---------DEBUG DIAGNOSTICS-------------------------------------------
          idt_src=3  ! output print level (1-5, fix)
          CALL dbg_print('BefExplDiff: trac_old', trac_old,  str_module,idt_src, in_subset=cells_in_domain)
          CALL dbg_print('BefExplDiff: div_adv_flux_vert',div_adv_flux_vert, str_module,idt_src, in_subset=cells_in_domain)
          !CALL dbg_print('BefExplDiff: flux_horz',flux_horz, str_module,idt_src, in_subset=cells_in_domain)
          !---------------------------------------------------------------------

          CALL tracer_diffusion_vert_explicit( p_patch_3d,            &
            & old_ocean_tracer%concentration, &
            & bc_top_tracer,                  &
            & a_v,                            &
            & div_diff_flx)
        ELSE
          div_diff_flx(1:nproma,1:n_zlev,1:p_patch%alloc_cell_blocks) = 0.0_wp
        ENDIF

        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=3  ! output print level (1-5, fix)
        CALL dbg_print('AftExplDiff: div_diff_flx', div_diff_flx, str_module, idt_src, in_subset=cells_in_domain)
        !---------------------------------------------------------------------

        ! top layer
        DO jb = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index
            DO level = 1, MIN(p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1
              ! NOTE: this will not wotk with partial cells on the top level
              delta_z     = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nold(1))%h(jc,jb)
              delta_z_new = p_patch_3d%p_patch_1d(1)%del_zlev_m(level) + p_os%p_prog(nnew(1))%h(jc,jb)

                new_ocean_tracer%concentration(jc,level,jb) = (old_ocean_tracer%concentration(jc,level,jb) * delta_z &
                  & -delta_t*(div_adv_flux_vert(jc,level,jb)-(div_diff_flux_horz(jc,level,jb)-div_adv_flux_horz(jc,level,jb))&
                  &-div_diff_flx(jc,level,jb)))/delta_z_new

                new_ocean_tracer%concentration(jc,level,jb) = new_ocean_tracer%concentration(jc,level,jb) + &
                  & (delta_t/delta_z) * bc_top_tracer(jc,jb)

            ENDDO

            DO level = 2, p_patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
              delta_z = p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,level,jb) ! p_patch_3d%p_patch_1d(1)%del_zlev_m(level)

              new_ocean_tracer%concentration(jc,level,jb)= old_ocean_tracer%concentration(jc,level,jb)   &
                & - (delta_t/delta_z) * (div_adv_flux_vert(jc,level,jb)-&
                &(div_diff_flux_horz(jc,level,jb)-div_adv_flux_horz(jc,level,jb))&
                &-div_diff_flx(jc,level,jb))

            ENDDO
          END DO
        END DO

        CALL sync_patch_array(sync_c, p_patch, new_ocean_tracer%concentration)

      ENDIF ! lvertical_diff_implicit

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('aft. AdvIndivTrac: trac_old', trac_old, str_module, 3, in_subset=cells_in_domain)
    CALL dbg_print('aft. AdvIndivTrac: trac_new', trac_new, str_module, 1, in_subset=cells_in_domain)
    !---------------------------------------------------------------------

  END SUBROUTINE advect_individual_tracer_ab_during_step
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


