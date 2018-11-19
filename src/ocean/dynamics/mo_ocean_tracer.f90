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
#include "iconfor_dsl_definitions.inc"
!----------------------------
MODULE mo_ocean_tracer
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  USE mo_math_types,                ONLY: t_cartesian_coordinates
  USE mo_impl_constants,            ONLY: sea_boundary, sea, min_dolic
  USE mo_math_constants,            ONLY: pi
  USE mo_ocean_nml,                 ONLY: n_zlev,                         &
    & threshold_min_t, threshold_max_t, threshold_min_s, threshold_max_s, &
    & iswm_oce,                 use_none,                                 &
    & flux_calculation_horz, flux_calculation_vert, miura_order1,         &
    & l_with_vert_tracer_diffusion, l_with_vert_tracer_advection,         &
    & GMRedi_configuration,GMRedi_combined,  GM_only,Redi_only ,          &
    & Cartesian_Mixing, tracer_threshold_min, tracer_threshold_max,       &
    & namelist_tracer_name, tracer_update_mode, use_none,                 &
    & GMREDI_COMBINED_DIAGNOSTIC,GM_INDIVIDUAL_DIAGNOSTIC,REDI_INDIVIDUAL_DIAGNOSTIC
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer, debug_check_level
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish !, message_text, message
  USE mo_ocean_boundcond,           ONLY: top_bound_cond_tracer
  USE mo_ocean_surface_types,       ONLY: t_ocean_surface
  USE mo_ocean_tracer_diffusion,    ONLY: tracer_diffusion_vertical_implicit, tracer_diffusion_vert_explicit,tracer_diffusion_horz
  USE mo_ocean_tracer_transport_horz, ONLY: advect_horz, diffuse_horz
  USE mo_ocean_tracer_transport_vert, ONLY: advect_flux_vertical
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_c, sync_e, sync_patch_array
  USE mo_timer,                     ONLY: timer_start, timer_stop, timers_level, timer_dif_vert, timer_extra30
  USE mo_statistics,                ONLY: global_minmaxmean, print_value_location
  USE mo_mpi,                       ONLY: my_process_is_stdio !global_mpi_barrier
!   USE mo_ocean_GM_Redi,             ONLY: calc_ocean_physics, prepare_ocean_physics
  USE mo_ocean_math_operators,      ONLY: div_oce_3d, verticalDiv_scalar_onFullLevels! !verticalDiv_scalar_midlevel
  USE mo_scalar_product,            ONLY: map_edges2edges_viacell_3d_const_z
  USE mo_physical_constants,        ONLY: clw, rho_ref,sitodbar
  USE mo_ocean_thermodyn,           ONLY: calculate_density, calc_potential_density
  USE mo_ocean_pp_scheme,           ONLY: calculate_rho4GMRedi
  USE mo_ocean_types,               ONLY: t_hydro_ocean_state
  USE mo_ocean_tracer_transport_types,  ONLY: t_ocean_tracer, t_tracer_collection, t_ocean_transport_state
  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=12)           :: str_module = 'oceTracer   '  ! Output of module for 1 line debug
  INTEGER :: idt_src    = 1               ! Level of detail for 1 line debug

  PUBLIC :: advect_ocean_tracers
  PUBLIC :: prepare_tracer_transport

CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE advect_ocean_tracers(old_tracers, new_tracers, transport_state,operators_coeff)
    TYPE(t_tracer_collection), INTENT(inout)      :: old_tracers
    TYPE(t_tracer_collection), INTENT(inout)      :: new_tracers
    TYPE(t_ocean_transport_state), TARGET         :: transport_state
    TYPE(t_operator_coeff),            INTENT(inout) :: operators_coeff

    !Local variables
    TYPE(t_patch_3d ), POINTER     :: patch_3d
    INTEGER :: tracer_index
    !-------------------------------------------------------------------------------
    patch_3d => old_tracers%patch_3d

    DO tracer_index = 1, old_tracers%no_of_tracers
      IF ( old_tracers%tracer(tracer_index)%is_advected) THEN
        CALL advect_diffuse_individual_tracer( patch_3d,    &
          & old_tracers%tracer(tracer_index),               &
          & transport_state, operators_coeff,                   &
          & new_tracers%tracer(tracer_index))
      ENDIF
    END DO

  END SUBROUTINE advect_ocean_tracers
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE copy_individual_tracer_ab(patch_3d, old_tracer, new_tracer)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: old_tracer
    TYPE(t_ocean_tracer), TARGET :: new_tracer


    INTEGER :: jc,level,jb
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp), POINTER :: trac_old(:,:,:), trac_new(:,:,:) ! temporary pointers to the concentration arrays
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2D
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_advection:advect_individual_tracer')
    !-------------------------------------------------------------------------------_
    trac_old => old_tracer%concentration
    trac_new => new_tracer%concentration
    patch_2D => patch_3d%p_patch_2d(1)
    all_cells => patch_2D%cells%all

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        DO level = 1, n_zlev
          new_tracer%concentration(jc,level,jb)= old_tracer%concentration(jc,level,jb)
        END DO
      END DO
    END DO
  END SUBROUTINE copy_individual_tracer_ab
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE advect_diffuse_individual_tracer(patch_3d, old_tracer,       &
    & transport_state, operators_coeff, new_tracer)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: old_tracer
    TYPE(t_ocean_tracer), TARGET :: new_tracer

    TYPE(t_ocean_transport_state), TARGET :: transport_state
    TYPE(t_operator_coeff),INTENT(inout) :: operators_coeff
!     REAL(wp), INTENT(inout), OPTIONAL :: horizontally_diffused_tracer(:,:,:)

    !Local variables

    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_advection:advect_diffuse_tracer')
    !-------------------------------------------------------------------------------_

    patch_2D => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('on entry: IndTrac: trac_old',old_tracer%concentration(:,:,:), &
      & str_module,idt_src, in_subset=patch_2D%cells%owned)
    !---------------------------------------------------------------------
    IF(tracer_update_mode == use_none ) THEN
      CALL copy_individual_tracer_ab( patch_3d,            &
        & old_tracer,new_tracer)
      RETURN
    ENDIF
   
    !Shallow water is done with horizontal advection
    IF(iswm_oce == 1) THEN
      CALL advect_diffuse_SW_tracer(patch_3d, old_tracer,       &
        & transport_state, operators_coeff,                      &
        & old_tracer%hor_diffusion_coeff,        &
        & new_tracer)
        
    !The 3D-case
    ELSE ! IF( iswm_oce /= 1) THEN

         CALL advect_diffuse_tracer( patch_3d, &
           & old_tracer,&
           & transport_state,            &
           & operators_coeff,      &
           & old_tracer%hor_diffusion_coeff,  &
           & old_tracer%ver_diffusion_coeff,  &
           & new_tracer)

    ENDIF

  END SUBROUTINE advect_diffuse_individual_tracer
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers for shallow water case
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE advect_diffuse_SW_tracer(patch_3d, old_tracer,       &
    & transport_state, operators_coeff,                      &
    & k_h,                                   &
    & new_tracer)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: old_tracer
    TYPE(t_ocean_transport_state), TARGET :: transport_state
    TYPE(t_operator_coeff),INTENT(inout) :: operators_coeff
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)       !horizontal mixing coeff
    TYPE(t_ocean_tracer), TARGET :: new_tracer
 
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
    trac_old => old_tracer%concentration
    trac_new => new_tracer%concentration

    patch_2D => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    delta_t = dtime

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('on entry: IndTrac: trac_old',trac_old(:,:,:) ,str_module,idt_src, in_subset=patch_2D%cells%owned)
    !---------------------------------------------------------------------
    
    !Shallow water is done with horizontal advection
      div_adv_flux_horz   (1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks) = 0.0_wp
      div_adv_flux_vert   (1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks) = 0.0_wp      
      div_diff_flux_horz  (1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks) = 0.0_wp

      !---------------------------------------------------------------------
      CALL advect_horz( patch_3d,         &
        & old_tracer%concentration, &
        & transport_state,                           &
        & operators_coeff,                     &
        & k_h,                            &
        & transport_state%h_old,         &
        & transport_state%h_new,         &
        & div_adv_flux_horz,              &
        & div_adv_flux_vert)

      CALL diffuse_horz( patch_3d,         &
        & old_tracer%concentration, &
        & transport_state,                           &
        & operators_coeff,                     &
        & k_h,                            &
        & transport_state%h_old,         &
        & transport_state%h_new,         &
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
            delta_z = transport_state%h_new(jc,jb) !- p_ext_data%oce%bathymetry_c(jc,jb)
            new_tracer%concentration(jc,level,jb)= old_tracer%concentration(jc,level,jb) - &
              & (delta_t/delta_z) * (div_adv_flux_horz(jc,level,jb)-div_diff_flux_horz(jc,level,jb))

          END DO
        END DO
      END DO

      CALL sync_patch_array(sync_c, patch_2D, new_tracer%concentration)
        
  END SUBROUTINE advect_diffuse_SW_tracer
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE advect_diffuse_tracer(patch_3d, old_tracer,       &
    & transport_state, operators_coeff,           &
    & k_h, a_v,                            &
    & new_tracer)!,        &
    ! & horizontally_diffused_tracer        )

    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: old_tracer
    TYPE(t_ocean_transport_state), TARGET :: transport_state
    TYPE(t_operator_coeff),INTENT(inout) :: operators_coeff
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)       !horizontal mixing coeff
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)       !vertical mixing coeff, in
    TYPE(t_ocean_tracer), TARGET :: new_tracer
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
    REAL(wp) :: top_bc(nproma)


    CHARACTER(len=*), PARAMETER :: method_name = 'mo_ocean_tracer:advect_diffuse_tracer'
    !-------------------------------------------------------------------------------
    trac_old => old_tracer%concentration
    trac_new => new_tracer%concentration

    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain
    delta_t = dtime

!     div_adv_flux_horz (1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
!     div_diff_flux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
!     flux_horz(1:nproma,1:n_zlev,1: 1:patch_3d%p_patch_2d(1)%nblks_e)=0.0_wp
    !---------------------------------------------------------------------
 
    ! these are probably not necessary
    div_diff_flx_vert = 0.0_wp
    div_adv_flux_vert = 0.0_wp
    div_adv_flux_horz = 0.0_wp
    div_diff_flux_horz = 0.0_wp
    !---------------------------------------------------------------------
    IF ( l_with_vert_tracer_advection ) THEN

      CALL advect_flux_vertical( patch_3d,&
        & old_tracer%concentration, &
        & transport_state,                           &
        & operators_coeff,                     &
        & div_adv_flux_vert)

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('aft. AdvFluxVert:divfluxvert',div_adv_flux_vert          ,str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

    ENDIF  ! l_with_vert_tracer_advection

    !---------------------------------------------------------------------
    CALL advect_horz( patch_3d,       &
      & old_tracer%concentration, &
      & transport_state,                           &
      & operators_coeff,                     &
      & k_h,                            &
      & transport_state%h_old,         &
      & transport_state%h_new,         &
      & div_adv_flux_horz,              &
      & div_adv_flux_vert)
    !---------------------------------------------------------------------

    IF(GMRedi_configuration==Cartesian_Mixing)THEN
      !horizontal diffusion, vertical is handled implicitely below
      CALL diffuse_horz( patch_3d,      &
      & old_tracer%concentration, &
      & transport_state,                           &
      & operators_coeff,                     &
      & k_h,                            &
      & transport_state%h_old,         &
      & transport_state%h_new,         &
      & div_diff_flux_horz)

    ELSE
      CALL finish(method_name, "wrong GMredi call")
    ENDIF  
      
    !Case: Implicit Vertical diffusion
    start_timer(timer_dif_vert,4)

    !Calculate preliminary tracer value out of horizontal advective and
    !diffusive fluxes and vertical advective fluxes, plus surface forcing.
    !Surface forcing applied as volume forcing at rhs, i.e.part of explicit term
    !in tracer (and also momentum) eqs. In this case, top boundary condition of
    !vertical Laplacians are homogeneous

!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, level, &
!ICON_OMP delta_z, delta_z_new, top_bc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
      IF (ASSOCIATED(old_tracer%top_bc)) THEN
        top_bc(:) = old_tracer%top_bc(:,jb)
      ELSE
        top_bc(:) = 0.0_wp
      ENDIF
        
      DO jc = start_cell_index, end_cell_index
        !TODO check algorithm: inv_prism_thick_c vs. del_zlev_m | * vs. /
        DO level = 1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1

          delta_z     = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+transport_state%h_old(jc,jb)
          delta_z_new = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)+transport_state%h_new(jc,jb)

          new_tracer%concentration(jc,level,jb)= &
            & (old_tracer%concentration(jc,level,jb) * delta_z &
            & - delta_t * (&
            &  div_adv_flux_horz(jc,level,jb) +div_adv_flux_vert(jc,level,jb)&
            & -div_diff_flux_horz(jc,level,jb)-div_diff_flx_vert(jc,level,jb))) / delta_z_new

          new_tracer%concentration(jc,level,jb) =         &
            & ( new_tracer%concentration(jc,level,jb) +   &
            & (delta_t  / delta_z_new) * top_bc(jc))

        ENDDO

        DO level = 2, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

          new_tracer%concentration(jc,level,jb) =                          &
            &  old_tracer%concentration(jc,level,jb) -                     &
            &  (delta_t /  patch_3d%p_patch_1D(1)%prism_thick_c(jc,level,jb))    &
            & * (div_adv_flux_horz(jc,level,jb)  +div_adv_flux_vert(jc,level,jb)&
            &  - div_diff_flux_horz(jc,level,jb)- div_diff_flx_vert(jc,level,jb))
          !   test
          !   IF( delta_z/= delta_z1)THEN
          !     write(0,*)'no agreement',level,jc,jb,&
          !     &patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb), patch_3d%p_patch_1D(1)%del_zlev_m(level)
          !     &patch_3d%p_patch_1D(1)%prism_thick_c(jc,level,jb), patch_3d%p_patch_1D(1)%del_zlev_m(level)
          !   ENDIF

        ENDDO

      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO


    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('BefImplDiff: div_adv_flux_vert',div_adv_flux_vert, str_module,idt_src, in_subset=cells_in_domain)
    CALL dbg_print('BefImplDiff: trac_inter', new_tracer%concentration,  str_module,idt_src, in_subset=cells_in_domain)
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
            DO level = 1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1      
             temp_tracer_before%concentration(jc,level,jb)=new_tracer%concentration(jc,level,jb)            
             temp_tracer_after%concentration(jc,level,jb)=new_tracer%concentration(jc,level,jb)                         
            END DO
          END DO
        ENDDO 
!ICON_OMP_END_PARALLEL_DO         
      ENDIF!IF(GMREDI_COMBINED_DIAGNOSTIC)
      
      !Vertical mixing: implicit and with coefficient a_v
      !that is the sum of PP-coeff and implicit part of Redi-scheme
      
      CALL tracer_diffusion_vertical_implicit( &
          & patch_3d,                   &
          & new_tracer,                 &
          & a_v,                        &
          & operators_coeff)           
          
    ENDIF!IF ( l_with_vert_tracer_diffusion )


    CALL sync_patch_array(sync_c, patch_2D, new_tracer%concentration)

    stop_timer(timer_dif_vert,4)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('AftImplDiff: trac_new', trac_new, str_module, idt_src, in_subset=cells_in_domain)
    !---------------------------------------------------------------------

!     IF (debug_check_level > 5 ) &
!       CALL check_min_max_tracer(info_text="After advect_diffuse_tracer", tracer=new_tracer%concentration,     &
!         & min_tracer=tracer_threshold_min(tracer_index), max_tracer=tracer_threshold_max(tracer_index), &
!         & tracer_name=namelist_tracer_name(tracer_index), in_subset=cells_in_domain)
  
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('aft. AdvIndivTrac: trac_old', trac_old, str_module, 3, in_subset=cells_in_domain)
    CALL dbg_print('aft. AdvIndivTrac: trac_new', trac_new, str_module, 3, in_subset=cells_in_domain)
!     CALL dbg_print('tracer tend', (new_tracer%concentration(:,:,:)&
!             &- old_tracer%concentration(:,:,:))/dtime, str_module, 4, in_subset=cells_in_domain)
      
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
  SUBROUTINE prepare_tracer_transport(patch_3d, ocean_state, operators_coeff)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_operator_coeff),INTENT(inout) :: operators_coeff
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
      ocean_state%p_diag%w_time_weighted(1:nproma, 1:n_zlev+1, jb) = &
        & ocean_state%p_diag%w(1:nproma, 1:n_zlev+1, jb)
    ENDDO
!ICON_OMP_END_DO
    !In case of shallow water we have to to this here, for 3D fluid its done within vertical velocity calculation
    IF(iswm_oce==1)THEN
     CALL map_edges2edges_viacell_3d_const_z( patch_3d, ocean_state%p_diag%vn_time_weighted, operators_coeff, &
        & ocean_state%p_diag%mass_flx_e)
    ENDIF    
 
    ! p_diag%w is compouted in_domain cells
    ! CALL sync_patch_array(SYNC_C, patch_2d,ocean_state%p_diag%w_time_weighted )

    ! This is already synced on edges_in_domain !
    ! CALL sync_patch_array(SYNC_E, patch_2d,ocean_state%p_diag%vn_time_weighted )

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
                        ocean_state%p_diag%p_vn_mean(je,level,jb)%x = 0.5_wp *                          &
                          & (ocean_state%p_diag%p_vn(edge_cell_index(1), level, edge_cell_block(1))%x + &
                          &  ocean_state%p_diag%p_vn(edge_cell_index(2), level, edge_cell_block(2))%x)
          !ocean_state%p_diag%p_vn_mean(je,level,jb)%x = 0.5_wp *                               &
          !  & (ocean_state%p_diag%p_vn_dual(edge_vert_index(1), level, edge_vert_block(1))%x + &
          !  & ocean_state%p_diag%p_vn_dual(edge_vert_index(2), level, edge_vert_block(2))%x)

          ! this is specific to the miura_order1_hflux_oce
          operators_coeff%moved_edge_position_cc(je,level,jb)%x =      &
            & operators_coeff%edge_position_cc(je,level,jb)%x          &
            & - half_time * ocean_state%p_diag%p_vn_mean(je,level,jb)%x

        END DO

        DO level = startLevel, fin_level
          upwind_index = MERGE(1, 2, ocean_state%p_diag%vn_time_weighted(je,level,jb) > 0.0_wp)

          operators_coeff%upwind_cell_idx(je,level,jb) = edge_cell_index(upwind_index)
          operators_coeff%upwind_cell_blk(je,level,jb) = edge_cell_block(upwind_index)

          operators_coeff%upwind_cell_position_cc(je,level,jb)%x = &
            & patch_2d%cells%cartesian_center(edge_cell_index(upwind_index), edge_cell_block(upwind_index))%x
          ! & operators_coeff%cell_position_cc(edge_cell_index(upwind_index), level, edge_cell_block(upwind_index))%x
        END DO

      END DO
    END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL

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


