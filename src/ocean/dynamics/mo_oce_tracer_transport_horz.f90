!>
!! Contains the implementation of the horizontal tracer transport routines for the ICON ocean model.
!! This comprises horizontal advection and diffusion
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/01)
!!
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
MODULE mo_oce_tracer_transport_horz
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  USE mo_math_constants,            ONLY: dbl_eps
  USE mo_impl_constants,            ONLY: sea_boundary
  USE mo_ocean_nml,                 ONLY: n_zlev, l_edge_based, ab_gam,                   &
    & upwind, central,lax_friedrichs, fct_horz, miura_order1, flux_calculation_horz,      &
    & fct_high_order_flux,  fct_low_order_flux,FCT_Limiter_horz, fct_limiter_horz_zalesak,&
    & fct_limiter_horz_minmod, fct_limiter_horz_posdef, l_with_horz_tracer_diffusion, l_with_horz_tracer_advection
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma, p_test_run
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_adv_horz, timer_hflx_lim, &
    & timer_dif_horz
  USE mo_oce_types,                 ONLY: t_hydro_ocean_state, t_ocean_tracer  
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish !, message_text, message
  USE mo_oce_physics
  USE mo_scalar_product,            ONLY: map_cell2edges_3d,map_edges2cell_3d, &
    & map_edges2edges_viacell_3d_const_z
  USE mo_oce_math_operators,        ONLY: div_oce_3d, grad_fd_norm_oce_3d
  USE mo_oce_diffusion,             ONLY: tracer_diffusion_horz
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff, no_primal_edges
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_c, sync_c1, sync_e, sync_patch_array, sync_patch_array_mult
  USE mo_mpi,                       ONLY: my_process_is_mpi_parallel
  
  
  
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=12)           :: str_module    = 'oceTracHorz '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  !
  ! PUBLIC INTERFACE
  !
  PUBLIC :: advect_diffuse_horz
  ! Private implemenation
  !
  PRIVATE :: advect_diffuse_edge_based
  PRIVATE :: advect_diffuse_cell_based
  PRIVATE :: fct_high_res
  PRIVATE :: flux_corr_transport_cell  
  PRIVATE :: flux_corr_transport_edge 
  PRIVATE :: hflx_limiter_oce_zalesak
  
  PRIVATE :: ratio_consecutive_gradients
  PRIVATE :: calculate_limiter_function

  PRIVATE :: upwind_hflux_oce  
  PRIVATE :: upwind_hflux_oce_mimetic
  PRIVATE :: central_hflux_oce  
  PRIVATE :: central_hflux_oce_mimetic
  PRIVATE :: miura_order1_hflux_oce  
  PRIVATE :: lax_friedrichs_hflux_oce 
  PRIVATE :: laxfr_upflux
  !PRIVATE :: elad
  !PRIVATE :: hflx_limiter_oce_mo
  
  
  INTEGER, PARAMETER :: top=1
  
CONTAINS
  !-----------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE advect_diffuse_horz( patch_3d,          &
    & trac_old,            &
    & p_os,                &
    & p_op_coeff,          &
    & k_h,                 &
    & h_old,               &
    & h_new,               &
    & flux_horz,           &
    & horizontally_diffused_tracer)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp)                               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET      :: p_os
    TYPE(t_operator_coeff), INTENT(in)     :: p_op_coeff
    REAL(wp), INTENT(in)                   :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                   :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                   :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)                :: flux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout), OPTIONAL      :: horizontally_diffused_tracer(:,:,:)
    !
    !
    !Local variables
   TYPE(t_subset_range), POINTER :: cells_in_domain!,edges_in_domain
    !-------------------------------------------------------------------------------
    cells_in_domain => patch_3d%p_patch_2d(1)%cells%in_domain
    !-------------------------------------------------------------------------------
    
    
    IF(l_edge_based)THEN
      CALL advect_diffuse_edge_based( patch_3d, &
      & trac_old,            &
      & p_os,                &
      & p_op_coeff,          &
      & k_h,                 &
      & h_old,               &
      & h_new,               &
      & flux_horz)
    
    ELSE
      CALL advect_diffuse_cell_based( patch_3d,          &
      & trac_old,            &
      & p_os,                &
      & p_op_coeff,          &
      & k_h,                 &
      & h_old,               &
      & h_new,               &
      & flux_horz,           &
      & horizontally_diffused_tracer)
    
    ENDIF
     
  END SUBROUTINE advect_diffuse_horz
  !-------------------------------------------------------------------------------


  !-----------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE advect_diffuse_cell_based( patch_3d,          &
    & trac_old,            &
    & p_os,                &
    & p_op_coeff,          &
    & k_h,                 &
    & h_old,               &
    & h_new,               &
    & flux_horz,           &
    & horizontally_diffused_tracer)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp)                             :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET    :: p_os
    TYPE(t_operator_coeff), INTENT(in)   :: p_op_coeff
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                 :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                 :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)              :: flux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout), OPTIONAL    :: horizontally_diffused_tracer(:,:,:)
    !
    !
    !Local variables
    INTEGER :: start_index, end_index
    INTEGER :: jc, jk, jb
    REAL(wp) :: z_adv_flux_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_div_adv_h  (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp) :: z_div_diff_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_diff_flux_h(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: p_vn_c(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-------------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------------
    z_adv_flux_h (1:nproma,1:n_zlev,1:patch_2d%nblks_e)=0.0_wp
    z_div_adv_h  (1:nproma,1:n_zlev,1:patch_2d%nblks_c)=0.0_wp
    z_div_diff_h (1:nproma,1:n_zlev,1:patch_2d%nblks_c)=0.0_wp
    z_diff_flux_h(1:nproma,1:n_zlev,1:patch_2d%nblks_e)=0.0_wp
    flux_horz    (1:nproma,1:n_zlev,1:patch_2d%nblks_c)=0.0_wp

    !Calculate tracer fluxes at edges
    !This step takes already the edge length into account
    !but not the edge height
    IF ( l_with_horz_tracer_advection ) THEN
      
      ! Initialize timer for horizontal advection
      IF (ltimer) CALL timer_start(timer_adv_horz)
      
      SELECT CASE(flux_calculation_horz)
      
      CASE(upwind)
      
        ! Note: this should be done only once, not for each tracer
        CALL map_edges2cell_3d(patch_3D, &
        & p_os%p_diag%vn_time_weighted,  &
        & p_op_coeff,                    &
        & p_vn_c)
        
        !!ICON_OMP_PARALLEL
        !!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, jb, start_index, end_index)
          DO jc = start_index, end_index
            DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
              p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x * trac_old(jc,jk,jb)
            END DO
          END DO
        END DO
        !!ICON_OMP_END_DO NOWAIT
        !!ICON_OMP_END_PARALLEL
        CALL upwind_hflux_oce_mimetic( patch_3d,&
        & p_vn_c, p_op_coeff,&
        & z_adv_flux_h )        
        
      CASE(central)  
      
        CALL map_edges2cell_3d(patch_3D, &
        &p_os%p_diag%vn_time_weighted,   &
        & p_op_coeff,                    &
        & p_vn_c)
        !!ICON_OMP_PARALLEL
        !!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, jb, start_index, end_index)
          DO jc = start_index, end_index
            DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
              p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x * trac_old(jc,jk,jb)
            END DO
          END DO
        END DO
        !!ICON_OMP_END_DO NOWAIT
        !!ICON_OMP_END_PARALLEL        
        CALL central_hflux_oce_mimetic( patch_3d,&
        & p_vn_c,                                &
        & z_adv_flux_h )        
        
      CASE(miura_order1)
        
        !MIURA-type upwind-biased estimate of tracer flux
        CALL miura_order1_hflux_oce( patch_3d,   &
          & trac_old,                            &
          & p_os%p_diag%mass_flx_e,              &
          & p_op_coeff,                          &
          & z_adv_flux_h )
        
      CASE(fct_horz)
      
       CALL flux_corr_transport_cell( patch_3d,  &
        & trac_old,                              &
        & p_os,                                  &
        & p_op_coeff,                            &
        & k_h,                                   &
        & h_old,                                 &
        & h_new,                                 &
        & z_adv_flux_h)              
              
      CASE default
        CALL finish('TRIM(advect_diffuse_flux_horz)',"This flux option is not supported")
        
      END SELECT
      CALL sync_patch_array(sync_e, patch_2d, z_adv_flux_h)
      ! Stop timer for horizontal advection
      IF (ltimer) CALL timer_stop(timer_adv_horz)      
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=5  ! output print level (1-5, fix)
      CALL dbg_print('aft. AdvHorz: adv_flux_h',z_adv_flux_h,str_module,idt_src,patch_2d%edges%owned)
      !---------------------------------------------------------------------
      
      !Calculate divergence of advective fluxes
      CALL div_oce_3d( z_adv_flux_h, patch_3D, p_op_coeff%div_coeff, z_div_adv_h )
      
    ELSE
      z_div_adv_h   (:,:,:) = 0.0_wp
    ENDIF ! l_with_horz_tracer_advection
    
    !The diffusion part: calculate horizontal diffusive flux
    IF ( l_with_horz_tracer_diffusion) THEN
      IF (ltimer) CALL timer_start(timer_dif_horz)
      
      CALL tracer_diffusion_horz( patch_3d,     &
        & trac_old,     &
        & p_os,         &
        & k_h,          &
        & z_diff_flux_h,&
        & subset_range = edges_in_domain)
      
      !Calculate divergence of diffusive fluxes
      CALL div_oce_3d( z_diff_flux_h, patch_3D, p_op_coeff%div_coeff, z_div_diff_h )
      IF (ltimer) CALL timer_stop(timer_dif_horz)
      
    ELSE
      z_diff_flux_h (:,:,:) = 0.0_wp
      z_div_diff_h  (:,:,:) = 0.0_wp
    ENDIF
    
    !Final step: calculate sum of advective and diffusive horizontal fluxes
    !!ICON_OMP_PARALLEL
    !!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
           flux_horz(jc,jk,jb) = z_div_diff_h(jc,jk,jb)-z_div_adv_h(jc,jk,jb)
        END DO
      END DO
    END DO
    !!ICON_OMP_END_DO NOWAIT
    !!ICON_OMP_END_PARALLEL
    
    CALL sync_patch_array(sync_c, patch_2d, flux_horz)
    
    IF (PRESENT(horizontally_diffused_tracer)) THEN
      horizontally_diffused_tracer(:,:,:) = z_div_diff_h(:,:,:)
      CALL sync_patch_array(sync_c, patch_2d, horizontally_diffused_tracer)

      CALL dbg_print('diffusion: t_diff', horizontally_diffused_tracer, "" , 3, &
        & patch_2D%cells%owned )

    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('AdvDifHorz: adv_flux_h'     ,z_adv_flux_h ,str_module,idt_src,patch_2d%edges%owned)
    CALL dbg_print('AdvDifHorz: div adv_flux_h' ,z_div_adv_h  ,str_module,idt_src,cells_in_domain)
    CALL dbg_print('AdvDifHorz: div diff_flux_h',z_div_diff_h ,str_module,idt_src,cells_in_domain)
    CALL dbg_print('AdvDifHorz: flux_horz'      ,flux_horz    ,str_module,idt_src,cells_in_domain)
    !---------------------------------------------------------------------
    
  END SUBROUTINE advect_diffuse_cell_based
  !-------------------------------------------------------------------------------
  
! !   
  !-----------------------------------------------------------------------
  SUBROUTINE advect_diffuse_edge_based( patch_3d,          &
    & trac_old,            &
    & p_os,                &
    & p_op_coeff,          &
    & k_h,                 &
    & h_old,               &
    & h_new,               &
    & flux_horz)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp)                             :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET    :: p_os
    TYPE(t_operator_coeff), INTENT(in)   :: p_op_coeff
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                 :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                 :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)              :: flux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !
    !
    !Local variables
    INTEGER :: start_index, end_index
    INTEGER :: jc, jk, jb    
    REAL(wp) :: z_adv_flux_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_div_adv_h  (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_div_diff_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_diff_flux_h(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_transport_horz:advect_diffuse_flux_horz')
    TYPE(t_patch), POINTER :: patch_2d
    !-------------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------------
    z_adv_flux_h (1:nproma,1:n_zlev,1:patch_2d%nblks_e)=0.0_wp
    z_div_adv_h  (1:nproma,1:n_zlev,1:patch_2d%nblks_c)=0.0_wp
    z_div_diff_h (1:nproma,1:n_zlev,1:patch_2d%nblks_c)=0.0_wp
    z_diff_flux_h(1:nproma,1:n_zlev,1:patch_2d%nblks_e)=0.0_wp
    flux_horz    (1:nproma,1:n_zlev,1:patch_2d%nblks_c)=0.0_wp

    !Calculate tracer fluxes at edges
    !This step takes already the edge length into account
    !but not the edge height
    IF ( l_with_horz_tracer_advection ) THEN

      ! Initialize timer for horizontal advection
      IF (ltimer) CALL timer_start(timer_adv_horz)
      
      SELECT CASE(flux_calculation_horz)
      
      CASE(upwind)
        
        !upwind estimate of tracer flux
        CALL upwind_hflux_oce( patch_3d,   &
           & trac_old,                     &
           & p_os%p_diag%mass_flx_e,       &
           & z_adv_flux_h )           
  
      CASE(central)
        
        !central estimate of tracer flux
        CALL central_hflux_oce( patch_3d, &
          & trac_old,                     &
          & p_os%p_diag%mass_flx_e,       &
          & z_adv_flux_h )
          
      CASE(lax_friedrichs)
      
        CALL lax_friedrichs_hflux_oce( patch_3d, &
          & trac_old,                            &
          & p_os%p_diag%mass_flx_e,              &
          & z_adv_flux_h )

      CASE(miura_order1)
        
        !MIURA-type upwind-biased estimate of tracer flux
        CALL miura_order1_hflux_oce( patch_3d, &
          & trac_old,                          &
          & p_os%p_diag%mass_flx_e,            &
          & p_op_coeff,                        &
          & z_adv_flux_h  )
                   
      CASE(fct_horz)
      
        CALL flux_corr_transport_edge( patch_3d,&
          & trac_old,                           &
          & p_os,                               &
          & p_op_coeff,                         &
          & k_h,                                &
          & h_old,                              &
          & h_new,                              &
          & z_adv_flux_h)      
              
      CASE default
        CALL finish('TRIM(advect_diffuse_flux_horz)',"This flux option is not supported")
        
      END SELECT     
      CALL sync_patch_array(sync_e, patch_2d, z_adv_flux_h)      
      
      !---------------------------------------------------------------------       
      ! Stop timer for horizontal advection
      IF (ltimer) CALL timer_stop(timer_adv_horz)      
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=5  ! output print level (1-5, fix)
      CALL dbg_print('aft. AdvHorz: adv_flux_h',z_adv_flux_h,str_module,idt_src,patch_2d%edges%owned)
      !---------------------------------------------------------------------               
      
      !Calculate divergence of advective fluxes
      CALL div_oce_3d( z_adv_flux_h, patch_3D, p_op_coeff%div_coeff, z_div_adv_h )
      
    ELSE
      z_div_adv_h   (:,:,:) = 0.0_wp
    ENDIF ! l_with_horz_tracer_diffusion
    
    !The diffusion part: calculate horizontal diffusive flux
    IF ( l_with_horz_tracer_diffusion)THEN 

      IF (ltimer) CALL timer_start(timer_dif_horz)
      CALL tracer_diffusion_horz( patch_3d,     &
        & trac_old,     &
        & p_os,         &
        & k_h,          &
        & z_diff_flux_h,&
        & subset_range = edges_in_domain)
      
      !Calculate divergence of diffusive fluxes
      CALL div_oce_3d( z_diff_flux_h, patch_3D, p_op_coeff%div_coeff, z_div_diff_h)
      IF (ltimer) CALL timer_stop(timer_dif_horz)
      
    ELSE
      z_diff_flux_h (:,:,:) = 0.0_wp
      z_div_diff_h  (:,:,:) = 0.0_wp
    ENDIF
    
    !Final step: calculate sum of advective and diffusive horizontal fluxes
    !!ICON_OMP_PARALLEL
    !!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          flux_horz(jc,jk,jb) = z_div_diff_h(jc,jk,jb)-z_div_adv_h(jc,jk,jb)
        END DO
      END DO
    END DO
    !!ICON_OMP_END_DO NOWAIT
    !!ICON_OMP_END_PARALLEL
    
    CALL sync_patch_array(sync_c, patch_2d, flux_horz)
    
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('AdvDifHorz: adv_flux_h'     ,z_adv_flux_h ,str_module,idt_src,patch_2d%edges%owned)
    CALL dbg_print('AdvDifHorz: div adv_flux_h' ,z_div_adv_h  ,str_module,idt_src,cells_in_domain)
    CALL dbg_print('AdvDifHorz: div diff_flux_h',z_div_diff_h ,str_module,idt_src,cells_in_domain)
    CALL dbg_print('AdvDifHorz: flux_horz'      ,flux_horz    ,str_module,idt_src,cells_in_domain)
    !---------------------------------------------------------------------   
  END SUBROUTINE advect_diffuse_edge_based
  !-------------------------------------------------------------------------------


  SUBROUTINE flux_corr_transport_edge( patch_3d,&
    & trac_old,                                 &
    & p_os,                                     &
    & p_op_coeff,                               &
    & k_h,                                      &
    & h_old,                                    &
    & h_new,                                    &
    & adv_flux_h)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp)                               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET      :: p_os
    TYPE(t_operator_coeff), INTENT(in)     :: p_op_coeff
    REAL(wp), INTENT(in)                   :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                   :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                   :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), TARGET, INTENT(inout)        :: adv_flux_h(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)!< variable in which the upwind flux is stored
    !Local Variables
    REAL(wp)         :: z_adv_flux_high(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp)         :: z_adv_flux_low (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    TYPE(t_patch), POINTER        :: patch_2d
    !-------------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------------

    !1) provide high- & low order tracer flux
    
    !low order flux: currently upwind
    SELECT CASE(fct_low_order_flux)
    
      CASE(upwind)
      
      !upwind serves as low order flux 
      CALL upwind_hflux_oce( patch_3d,  &
      & trac_old,                       &
      & p_os%p_diag%mass_flx_e,         & 
      & z_adv_flux_low ) 
      
    CASE DEFAULT
      CALL finish('TRIM(flux_corr_transport_edge)',"This low-order  option is not supported")   
    END SELECT
    
    !high order flux options
    SELECT CASE(fct_high_order_flux)
    
    CASE(central)
    
      !central as high order flux
      CALL central_hflux_oce( patch_3d,   &
        & trac_old,                       &
        & p_os%p_diag%mass_flx_e,         &
        & z_adv_flux_high) 
      
    CASE(lax_friedrichs)
    
      CALL lax_friedrichs_hflux_oce( patch_3d, &
        & trac_old,                            &
        & p_os%p_diag%mass_flx_e,              &
        & z_adv_flux_high)
    
    CASE(miura_order1)
    
      CALL miura_order1_hflux_oce( patch_3d, &
        & trac_old,                          &
        & p_os%p_diag%mass_flx_e,            &
        & p_op_coeff,                        &
        & z_adv_flux_high )
    
    CASE DEFAULT
      CALL finish('TRIM(flux_corr_transport_edge)',"This high-order  option is not supported")
    END SELECT
    


    !2)call limiter 
    SELECT CASE(fct_limiter_horz)
    
    CASE(fct_limiter_horz_minmod)
    
      CALL fct_high_res( patch_3d, &
      & trac_old,                              &
      & p_os%p_diag%vn_time_weighted,          &
      & z_adv_flux_low,                        &    
      & z_adv_flux_high,                       &      
      & p_op_coeff,                            &
      & k_h,                                   &
      & adv_flux_h )

      !CALL hflx_limiter_oce_posdef( patch_3d,&
      !& trac_old,                    &
      !& z_adv_flux_high,              &    
      !& p_op_coeff )      
      
    CASE(fct_limiter_horz_zalesak)

       CALL hflx_limiter_oce_zalesak( patch_3d, &
       & trac_old,                              &
       & p_os%p_diag%mass_flx_e,                &
       & z_adv_flux_low,                        &
       & z_adv_flux_high,                       &    
       & adv_flux_h,                            &
       & p_op_coeff,                            &
       & h_old,                                 &
       & h_new              )       
       
    CASE(fct_limiter_horz_posdef)  
      CALL hflx_limiter_oce_posdef( patch_3d,  &
      & trac_old,                              &
      & z_adv_flux_high,                       &    
      & p_op_coeff )      

    CASE DEFAULT
      CALL finish('TRIM(flux_corr_transport_edge)',"This limiter_type option is not supported")
    END SELECT

    
  END SUBROUTINE flux_corr_transport_edge
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE flux_corr_transport_cell( patch_3d, &
    & trac_old,                                   & 
    & p_os,                                       &
    & p_op_coeff,                                 &
    & k_h,                                        &
    & h_old,                                      &
    & h_new,                                      &
    & adv_flux_h)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp)                               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET      :: p_os
    TYPE(t_operator_coeff), INTENT(in)     :: p_op_coeff
    REAL(wp), INTENT(in)                   :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                   :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                   :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET, INTENT(inout)         :: adv_flux_h(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)!< variable in which the upwind flux is stored
    !Local Variables
    REAL(wp) :: z_adv_flux_high(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_adv_flux_low (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-------------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------------
    !1) provide high- & low order tracer flux
    !upwind serves as low order flux 
    CALL upwind_hflux_oce( patch_3d,  &
    & trac_old,                       &
    & p_os%p_diag%mass_flx_e,         &
    & z_adv_flux_low )     
  
    !mimetic fluc calculation high order flux
     CALL map_edges2edges_viacell_3d_const_z( patch_3d, &
     & p_os%p_diag%vn_time_weighted,                    &
     & p_op_coeff,                                      &
     & z_adv_flux_high,                                 &
     & trac_old)  
    
    !2)call limiter 
    SELECT CASE(fct_limiter_horz)
    
    CASE(fct_Limiter_horz_minmod)
    
      CALL fct_high_res( patch_3d, &
      & trac_old,                              &
      & p_os%p_diag%vn_time_weighted,          &
      & z_adv_flux_low,                        &    
      & z_adv_flux_high,                       &      
      & p_op_coeff,                            &
      & k_h,                                   &
      & adv_flux_h )
      
    CASE(fct_limiter_horz_zalesak)

      CALL hflx_limiter_oce_zalesak( patch_3d, &
      & trac_old,                              &
      & p_os%p_diag%mass_flx_e,                &
      & z_adv_flux_low,                        & 
      & z_adv_flux_high,                       &            
      & adv_flux_h,                            &
      & p_op_coeff,                            &
      & h_old,                                 &
      & h_new              )  
      
    CASE(fct_limiter_horz_posdef)  
      CALL hflx_limiter_oce_posdef( patch_3d,  &
      & trac_old,                              &
      & z_adv_flux_high,                       &    
      & p_op_coeff )      
      
    CASE DEFAULT
     CALL finish('TRIM(flux_corr_transport_h)',"This limiter_type option is not supported")
    END SELECT
    
  END SUBROUTINE flux_corr_transport_cell
  !-------------------------------------------------------------------------------
 
  
  
  !-------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2013).
  SUBROUTINE fct_high_res( patch_3d, trac_old, vn_e,adv_flux_low,adv_flux_high, &
    & p_op_coeff, k_h, adv_flux_h)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)           :: trac_old(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< advected cell centered variable
    REAL(wp), INTENT(inout)           :: vn_e (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) !< normal velocity on edges
    REAL(wp), INTENT(IN)              :: adv_flux_low (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(IN)              :: adv_flux_high(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff)            :: p_op_coeff
    REAL(wp), INTENT(in)              :: k_h(:,:,:)         !horizontal mixing coeff
    
    REAL(wp), INTENT(OUT) :: adv_flux_h(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)!< variable in which the upwind flux is stored

    REAL(wp) :: z_diff_flux_h(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)    
    REAL(wp) :: z_consec_grad(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_mflux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_limit_phi(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_limit_psi(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_limit_sigma(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    INTEGER :: start_index_e, end_index_e
    INTEGER :: jk, jb, je
    INTEGER, DIMENSION(:,:,:), POINTER :: cell_of_edge_idx, cell_of_edge_blk
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    TYPE(t_patch), POINTER        :: patch_2d
    !-------------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !------------------------------------------------------------------------------   
      
    ! Set pointers to index-arrays
    ! line and block indices of two neighboring cells
    cell_of_edge_idx  => patch_2d%edges%cell_idx
    cell_of_edge_blk  => patch_2d%edges%cell_blk
    
    z_consec_grad(1:nproma,1:n_zlev,1:patch_2d%nblks_e) = 0.0_wp
    z_mflux(1:nproma,1:n_zlev,1:patch_2d%nblks_e)       = 0.0_wp
    z_limit_phi(1:nproma,1:n_zlev,1:patch_2d%nblks_e)   = 0.0_wp
    z_limit_psi(1:nproma,1:n_zlev,1:patch_2d%nblks_e)   = 0.0_wp
    z_limit_sigma(1:nproma,1:n_zlev,1:patch_2d%nblks_e) = 0.0_wp   
    
    !1) calculate phi-part of limiter
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index_e, end_index_e)
      DO jk = 1, n_zlev
        DO je = start_index_e, end_index_e
        
           !calculate the difference between adjacent tracer values, without
           !dividing by dual edge length.
           z_diff_flux_h(je,jk,jb) = trac_old(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2))&
                                   &-trac_old(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1))
           
           !outflow case
           !if vn_e(je,jk,jb) is positive then the flow is from cell 1 towards cell 2
           z_mflux(je,jk,jb) = abs(vn_e(je,jk,jb))

          !This corresponds to eq. (14) in Casulli-Zanolli.
          !Their "little" phi is here called sigma. We have
          !canceled down their primal edge-length and thickness 
          !in numerator and denominator of sigma.
           IF (z_mflux(je,jk,jb) == 0.0_wp) THEN
            z_limit_sigma(je,jk,jb) = 1.0_wp   !  z_mflux can be zero!
          ELSE
            z_limit_sigma(je,jk,jb) = MIN(1.0_wp, &
              & 2.0_wp * k_h(je,jk,jb)/(patch_2d%edges%dual_edge_length(je,jb)*z_mflux(je,jk,jb)))
          ENDIF

        END DO
      END DO
    END DO
    CALL sync_patch_array(sync_e, patch_2d, z_diff_flux_h)    
    CALL sync_patch_array(sync_e, patch_2d, z_mflux)
    CALL sync_patch_array(sync_e, patch_2d, z_limit_sigma)        

    !This corresponds to (16) in Casulli-Zanolli.
     CALL ratio_consecutive_gradients(patch_3d,p_op_coeff,vn_e, trac_old,z_consec_grad)

    !3) calculate psi-part of limiter (see (17)-(19) in Casulli-Zanolli).
    CALL calculate_limiter_function(patch_3d, z_limit_sigma, z_consec_grad,z_limit_phi)

    !4) Calculate limited advective flux
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index_e, end_index_e)
      DO jk = 1, n_zlev
        DO je = start_index_e, end_index_e

          IF(z_consec_grad(je,jk,jb)>0)THEN
            z_limit_psi(je,jk,jb) = z_limit_phi(je,jk,jb)-z_limit_sigma(je,jk,jb)
          ELSE
            z_limit_psi(je,jk,jb) = 0.0_wp
          ENDIF
          !only low gives upwind, low +mflux*diff and without limiter gives central
          adv_flux_h(je,jk,jb) =           &
          & adv_flux_low (je,jk,jb)&
          & +0.5_wp*z_mflux(je,jk,jb)*z_diff_flux_h(je,jk,jb)*z_limit_psi(je,jk,jb)
          
        END DO
      END DO
    END DO
    CALL sync_patch_array(sync_e, patch_2d, z_limit_psi)        
    CALL sync_patch_array(sync_e, patch_2d, adv_flux_h)        
    idt_src=5  ! output print level (1-5, fix)
    CALL dbg_print('AdvDifHorzHR: limit_phi ', z_limit_phi,  str_module,idt_src, patch_2d%edges%owned)
    CALL dbg_print('AdvDifHorzHR: limit_psi ', z_limit_psi,  str_module,idt_src, patch_2d%edges%owned)

  END SUBROUTINE fct_high_res
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates ratio of consecutive gradients following Casulli-Zanolli.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2013).
  !!
   SUBROUTINE ratio_consecutive_gradients( patch_3d, p_op_coeff,      &
    & vn_time_weighted, &
    & trac_old,         &
    & consec_grad)
    
    TYPE(t_patch_3d),TARGET, INTENT(in):: patch_3d
    TYPE(t_operator_coeff)             :: p_op_coeff    
    REAL(wp), INTENT(in)               :: vn_time_weighted(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%nblks_c)
    REAL(wp), INTENT(inout)            :: consec_grad(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%nblks_e)
    !
    !Local variables
    !REAL(wp) :: delta_z
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: jc, jk, jb!, je!, ic,ib
    INTEGER :: i_edge, ii_e, ib_e
    REAL(wp) :: z_diff_trac    
    REAL(wp) :: z_cellsum_mass_flx_in(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_c)
    REAL(wp) :: z_cellsum_tracdiff_flx_in(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_c)    
    REAL(wp) :: z_mflux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_c,1:no_primal_edges)
    INTEGER, DIMENSION(:,:,:), POINTER :: cell_of_edge_idx, cell_of_edge_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    cells_in_domain => patch_2D%cells%in_domain
    !-------------------------------------------------------------------------------
    consec_grad              (1:nproma,1:n_zlev,1:patch_2D%nblks_e) = 0.0_wp
    z_cellsum_mass_flx_in    (1:nproma,1:n_zlev,1:patch_2d%nblks_c) = 0.0_wp
    z_cellsum_tracdiff_flx_in(1:nproma,1:n_zlev,1:patch_2d%nblks_c) = 0.0_wp
    z_mflux(1:nproma,1:n_zlev,1:patch_2d%nblks_c,1:no_primal_edges) = 0.0_wp
      
    ! Set pointers to index-arrays
    ! line and block indices of two neighboring cells
    cell_of_edge_idx  => patch_2D%edges%cell_idx
    cell_of_edge_blk  => patch_2D%edges%cell_blk
    edge_of_cell_idx  => patch_2D%cells%edge_idx
    edge_of_cell_blk  => patch_2D%cells%edge_blk
    
#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jk = 1, n_zlev
        DO jc = i_startidx_c, i_endidx_c
          
          DO i_edge=1,no_primal_edges
            ii_e = patch_2D%cells%edge_idx(jc,jb,i_edge)
            ib_e = patch_2D%cells%edge_blk(jc,jb,i_edge)
            
            !calculate mass flux of actual edge.
            z_mflux(jc,jk,jb, i_edge)=p_op_coeff%div_coeff(jc,jk,jb,i_edge)*patch_2D%cells%area(jc,jb)  &
            & * vn_time_weighted(edge_of_cell_idx(jc,jb,i_edge),jk,edge_of_cell_blk(jc,jb,i_edge))&
            & * patch_3d%p_patch_1d(1)%prism_thick_e(ii_e,jk,ib_e)
            
            !consider here inflow edges only
            IF( z_mflux(jc,jk,jb, i_edge)<0.0)THEN
            
            !IF( z_mflux(jc,jk,jb, i_edge)>0.0)THEN                                    
              !sum up mass fluxes over edge
              z_cellsum_mass_flx_in(jc,jk,jb) = z_cellsum_mass_flx_in(jc,jk,jb)&
              &+z_mflux(jc,jk,jb, i_edge)
                  
              !one of the two terms below is zero.
              z_diff_trac = &
              & (trac_old(cell_of_edge_idx(ii_e,ib_e,2),jk, cell_of_edge_blk(ii_e,ib_e,2))-trac_old(jc,jk,jb))&
              &+(trac_old(cell_of_edge_idx(ii_e,ib_e,1),jk, cell_of_edge_blk(ii_e,ib_e,1))-trac_old(jc,jk,jb) )
 
              !This is in Casulli Zanolli, eq. (16) the numerator of
              !the consecutive gradient
              z_cellsum_tracdiff_flx_in(jc,jk,jb) &
              &= z_cellsum_tracdiff_flx_in(jc,jk,jb)+z_mflux(jc,jk,jb, i_edge)*z_diff_trac!&
              !&*patch%cells%edge_orientation(jc,jb,i_edge)  !diff_flux_h(iii_e,jk,iib_e)!
            END IF
            !ENDIF
          END DO
        END DO
      END DO
    END DO
    DO i_edge=1,no_primal_edges
      CALL sync_patch_array(sync_c, patch_2D, z_mflux(:,:,:,i_edge))    
    END DO  
    CALL sync_patch_array(sync_c, patch_2D, z_cellsum_mass_flx_in)
    CALL sync_patch_array(sync_c, patch_2D, z_cellsum_tracdiff_flx_in)    
    
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jk = 1, n_zlev
        DO jc = i_startidx_c, i_endidx_c
          
          DO i_edge=1,no_primal_edges
            
            !For positive velocities the flow is directed from cell1 to cell2
             IF(z_mflux(jc,jk,jb, i_edge)>=0.0)THEN                   
             !IF(z_mflux(jc,jk,jb, i_edge)<0.0)THEN                               

                ii_e = patch_2D%cells%edge_idx(jc,jb,i_edge)
                ib_e = patch_2D%cells%edge_blk(jc,jb,i_edge)        
                
                z_diff_trac = &
                & ( trac_old(cell_of_edge_idx(ii_e,ib_e,2),jk, cell_of_edge_blk(ii_e,ib_e,2))-trac_old(jc,jk,jb))&
                &+( trac_old(cell_of_edge_idx(ii_e,ib_e,1),jk, cell_of_edge_blk(ii_e,ib_e,1))-trac_old(jc,jk,jb))
            
              
              IF(z_cellsum_mass_flx_in(jc,jk,jb)/=0.0_wp.AND.z_diff_trac/=0.0_wp)THEN
                
                consec_grad(edge_of_cell_idx(jc,jb,i_edge),jk,edge_of_cell_blk(jc,jb,i_edge))&
                & = z_cellsum_tracdiff_flx_in(jc,jk,jb)/z_cellsum_mass_flx_in(jc,jk,jb)
                                
                consec_grad(edge_of_cell_idx(jc,jb,i_edge),jk,edge_of_cell_blk(jc,jb,i_edge))&
                & = (consec_grad(ii_e,jk,ib_e)/z_diff_trac)!*patch%cells%edge_orientation(jc,jb,i_edge) 
              ENDIF
            END IF
          END DO
        END DO
      END DO
    END DO    
    CALL sync_patch_array(sync_e, patch_2D, consec_grad)
    
  END SUBROUTINE ratio_consecutive_gradients

  !-------------------------------------------------------------------------------
  !>
  !! Calculates ratio of consecutive gradients following Casulli-Zanolli.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2013).
  !!
  SUBROUTINE calculate_limiter_function( p_patch_3d,   &
    & limit_sigma,                            &
    & consec_grad,                            &
    & limit_phi)
    
    TYPE(t_patch_3d),TARGET, INTENT(in):: p_patch_3d
    REAL(wp), INTENT(in)               :: limit_sigma(1:nproma,1:n_zlev,1:p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)               :: consec_grad(1:nproma,1:n_zlev,1:p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)            :: limit_phi(1:nproma,1:n_zlev,1:p_patch_3d%p_patch_2d(1)%nblks_e)
    !
    !Local variables
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: jk, jb, je
    
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    TYPE(t_patch), POINTER        :: patch_2D
    !-------------------------------------------------------------------------------
    patch_2D        => p_patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    cells_in_domain => patch_2D%cells%in_domain
    !-------------------------------------------------------------------------------
    limit_phi(1:nproma,1:n_zlev,1:patch_2D%nblks_e) = 0.0_wp
    
    !SELECT CASE(limiter_type)    
    !CASE(minmod)
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        DO jk = 1, n_zlev
          DO je = i_startidx_e, i_endidx_e
            
            limit_phi(je,jk,jb)&
              & =MAX(limit_sigma(je,jk,jb), MIN(1.0_wp,consec_grad(je,jk,jb)))
          END DO
        END DO
      END DO
      
!     CASE(superbee)
!       DO jb = edges_in_domain%start_block, edges_in_domain%end_block
!         CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
!         DO jk = 1, n_zlev
!           DO je = i_startidx_e, i_endidx_e
!             limit_phi(je,jk,jb) = MAX(limit_sigma(je,jk,jb),&
!               & MIN(1.0_wp,2.0_wp*consec_grad(je,jk,jb)),&
!               & MIN(2.0_wp,consec_grad(je,jk,jb)))
!           END DO
!         END DO
!       END DO
!     CASE(van_leer)
!       DO jb = edges_in_domain%start_block, edges_in_domain%end_block
!         CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
!         DO jk = 1, n_zlev
!           DO je = i_startidx_e, i_endidx_e
!             limit_phi(je,jk,jb) = MAX(limit_sigma(je,jk,jb),&
!               & (consec_grad(je,jk,jb)+ABS(consec_grad(je,jk,jb))&
!               & /(1.0_wp+ABS(consec_grad(je,jk,jb)))))
!               IF(limit_phi(je,jk,jb)>2.0)THEN
!               write(0,*)'Something wrong: limiter exceeds 2.0'
!               ENDIF
!           END DO
!         END DO
!       END DO
!     END SELECT
     CALL sync_patch_array(sync_e, patch_2D, limit_phi)
!write(0,*)'limiter-function: max-min',maxval(limit_phi(:,1,:)),minval(limit_phi(:,1,:))   
  END SUBROUTINE calculate_limiter_function
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !>
  !! First order upwind scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !!
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  SUBROUTINE upwind_hflux_oce_mimetic( patch_3d, flux_cc,p_op_coeff,&
    & edge_upwind_flux, opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    TYPE(t_cartesian_coordinates)      :: flux_cc(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_operator_coeff), INTENT(in) :: p_op_coeff
    REAL(wp), INTENT(inout)            :: edge_upwind_flux(nproma,n_zlev, patch_3d%p_patch_2d(1)%nblks_e)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level    ! optional vertical end level
    ! local variables
    !INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc
    INTEGER :: start_level, end_level
    INTEGER :: start_index, end_index
    INTEGER :: je, jk, jb        !< index of edge, vert level, block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d         => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_start_level) ) THEN
      start_level = opt_start_level
    ELSE
      start_level = 1
    END IF
    IF ( PRESENT(opt_end_level) ) THEN
      end_level = opt_end_level
    ELSE
      end_level = n_zlev
    END IF
    
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    ! for no-slip boundary conditions, boundary treatment for tracer (zero at leteral walls)
    !is implicit done via velocity boundary conditions
    !
    ! line and block indices of two neighboring cells
    ! loop through all patch edges (and blocks)
    
    !!ICON_OMP_PARALLEL
    !!ICON_OMP_DO PRIVATE(start_index, end_index, je, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      edge_upwind_flux(:,:,jb) = 0.0_wp
      DO je = start_index, end_index
        DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,jb), end_level)
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          edge_upwind_flux(je,jk,jb) = &
            & DOT_PRODUCT(flux_cc(&
            & p_op_coeff%upwind_cell_idx(je,jk,jb), jk, p_op_coeff%upwind_cell_blk(je,jk,jb))%x, &
            & patch_2d%edges%primal_cart_normal(je,jb)%x)
          
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
    !!ICON_OMP_END_DO NOWAIT
    !!ICON_OMP_END_PARALLEL
    
    
  END SUBROUTINE upwind_hflux_oce_mimetic
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  !! First order upwind scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !!
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  SUBROUTINE central_hflux_oce_mimetic( patch_3d, flux_cc,&
    & edge_upwind_flux, opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d                                  !< patch on which computation is performed
    TYPE(t_cartesian_coordinates)     :: flux_cc  (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)           :: edge_upwind_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level    ! optional vertical end level
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER :: start_level, end_level
    INTEGER :: start_index, end_index
    INTEGER :: je, jk, jb      !< index of edge, vert level, block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_cartesian_coordinates):: flux_mean_e
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d         => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_start_level) ) THEN
      start_level = opt_start_level
    ELSE
      start_level = 1
    END IF
    IF ( PRESENT(opt_end_level) ) THEN
      end_level = opt_end_level
    ELSE
      end_level = n_zlev
    END IF
    flux_mean_e%x= 0.0_wp
    !
    ! for no-slip boundary conditions, boundary treatment for tracer (zero at leteral walls)
    !is implicit done via velocity boundary conditions
    !
    ! line and block indices of two neighboring cells
    iilc => patch_2d%edges%cell_idx
    iibc => patch_2d%edges%cell_blk
    
    ! loop through all patch edges (and blocks)
    !!ICON_OMP_PARALLEL
    !!ICON_OMP_DO PRIVATE(start_index, end_index, je, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      edge_upwind_flux(:,:,jb) = 0.0_wp
      DO je = start_index, end_index
        DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,jb), end_level)
          !
          ! compute the central flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          flux_mean_e%x=0.5_wp*(flux_cc(iilc(je,jb,1),jk,iibc(je,jb,1))%x&
            & +flux_cc(iilc(je,jb,2),jk,iibc(je,jb,2))%x)
          
          edge_upwind_flux(je,jk,jb)=DOT_PRODUCT(flux_mean_e%x,&
            & patch_2d%edges%primal_cart_normal(je,jb)%x)
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
    !!ICON_OMP_END_DO NOWAIT
    !!ICON_OMP_END_PARALLEL
    
  END SUBROUTINE central_hflux_oce_mimetic
  !-------------------------------------------------------------------------------
  
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
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  !<Optimize:inUse>
  SUBROUTINE upwind_hflux_oce( patch_3d, pvar_c, pvn_e, edge_upwind_flux, opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)              :: pvar_c   (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)      !< advected cell centered variable
    REAL(wp), INTENT(in)              :: pvn_e    (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)       !< normal velocity on edges
    REAL(wp), INTENT(inout)             :: edge_upwind_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level    ! optional vertical end level
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER :: start_level, end_level
    INTEGER :: start_index, end_index
    INTEGER :: je, jk, jb         !< index of edge, vert level, block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d         => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_start_level) ) THEN
      start_level = opt_start_level
    ELSE
      start_level = 1
    END IF
    IF ( PRESENT(opt_end_level) ) THEN
      end_level = opt_end_level
    ELSE
      end_level = n_zlev
    END IF
    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    !
    ! for no-slip boundary conditions, boundary treatment for tracer (zero at lateral walls)
    !is implicit done via velocity boundary conditions
    !
    ! line and block indices of two neighboring cells
    iilc => patch_2d%edges%cell_idx
    iibc => patch_2d%edges%cell_blk
    
    ! loop through all patch edges (and blocks)
    !!ICON_OMP_PARALLEL
    !!ICON_OMP_DO PRIVATE(start_index, end_index, je, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      edge_upwind_flux(:,:,jb) = 0.0_wp
      DO je = start_index, end_index
        DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,jb), end_level)
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          edge_upwind_flux(je,jk,jb) =  &
            & laxfr_upflux( pvn_e(je,jk,jb), pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), &
            & pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )
          
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
    !!ICON_OMP_END_DO NOWAIT
    !!ICON_OMP_END_PARALLEL
    
  END SUBROUTINE upwind_hflux_oce
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  !! Central scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using central fluxes
  !!
  !! @par Revision History
  !! Peter korn, MPI-M, 2011
  !!
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  SUBROUTINE central_hflux_oce( patch_3d, pvar_c, pvn_e, edge_flux )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)           :: pvar_c   (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    !< advected cell centered variable
    REAL(wp), INTENT(inout)           :: pvn_e    (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     !< normal velocity on edges
    REAL(wp), INTENT(inout)           :: edge_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) !< variable in which the upwind flux is stored
    
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER :: start_index, end_index
    INTEGER :: je, jk, jb
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d         => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    ! line and block indices of two neighboring cells
    iilc => patch_2d%edges%cell_idx
    iibc => patch_2d%edges%cell_blk
    
    ! loop through all patch edges (and blocks)
    !!ICON_OMP_PARALLEL
    !!ICON_OMP_DO PRIVATE(start_index, end_index, je, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      edge_flux(:,:,jb) = 0.0_wp
      DO je = start_index, end_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
          !
          ! compute the central flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          edge_flux(je,jk,jb) =  0.5_wp * pvn_e(je,jk,jb)    &
            & * (   pvar_c(iilc(je,jb,1), jk, iibc(je,jb,1)) &
            & + pvar_c(iilc(je,jb,2), jk, iibc(je,jb,2)))
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
    !!ICON_OMP_END_DO NOWAIT
    !!ICON_OMP_END_PARALLEL
    
  END SUBROUTINE central_hflux_oce
  !-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Central scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using Lax-Friedrichs fluxes
  !!
  !! @par Revision History
  !! Peter korn, MPI-M, 2011
  !!
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  SUBROUTINE lax_friedrichs_hflux_oce( patch_3d, pvar_c, pvn_e, edge_flux )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)              :: pvar_c   (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< advected cell centered variable
    REAL(wp), INTENT(inout)              :: pvn_e    (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)  
    REAL(wp), INTENT(inout)              :: edge_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER :: start_index, end_index
    INTEGER :: je, jk, jb
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d         => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    ! line and block indices of two neighboring cells
    iilc => patch_2d%edges%cell_idx
    iibc => patch_2d%edges%cell_blk
    
    ! loop through all patch edges (and blocks)
    !!ICON_OMP_PARALLEL
    !!ICON_OMP_DO PRIVATE(start_index, end_index, je, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      edge_flux(:,:,jb) = 0.0_wp
      DO je = start_index, end_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
          !
          edge_flux(je,jk,jb) =  0.5_wp * pvn_e(je,jk,jb)    &
            & * (   pvar_c(iilc(je,jb,1), jk, iibc(je,jb,1)) &
            & + pvar_c(iilc(je,jb,2), jk, iibc(je,jb,2)))    &
            &          +0.5_wp*pvn_e(je,jk,jb)*pvn_e(je,jk,jb)*dtime&
            &          * patch_2d%edges%inv_dual_edge_length(je,jb)   &
            &        *( pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2))      &
            &          -pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)))
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
    !!ICON_OMP_END_DO NOWAIT
    !!ICON_OMP_END_PARALLEL
    
  END SUBROUTINE lax_friedrichs_hflux_oce
  !-------------------------------------------------------------------------------
  
  
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
  !!  mpi parallelized, the result is NOT synced. Should be done in the calling method if required
  !!
  SUBROUTINE miura_order1_hflux_oce( patch_3d,pvar_c, pvn_e, &
    & p_op_coeff, edge_upwind_flux, &
    & opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)           :: pvar_c(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< advected cell centered variable
    REAL(wp), INTENT(inout)           :: pvn_e (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) !< normal velocity on edges
    TYPE(t_operator_coeff)            :: p_op_coeff
    REAL(wp), INTENT(inout)           :: edge_upwind_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)!< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level      ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level      ! optional vertical end level
    
    ! local variables
    INTEGER :: start_level, end_level
    INTEGER :: start_index, end_index
    INTEGER :: je, jk, jb
    !INTEGER  :: il_v1, il_v2, ib_v1, ib_v2
    INTEGER :: il_c, ib_c
    REAL(wp)                      :: z_gradc   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: z_gradc_cc(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d         => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_start_level) ) THEN
      start_level = opt_start_level
    ELSE
      start_level = 1
    END IF
    IF ( PRESENT(opt_end_level) ) THEN
      end_level = opt_end_level
    ELSE
      end_level = n_zlev
    END IF
    !------------------------------------------------------------------
    ! ! !-------------upwind comparison
    ! !     iilc => patch_2d%edges%cell_idx
    ! !     iibc => patch_2d%edges%cell_blk
    ! !     ! loop through all patch edges (and blocks)
    ! !     DO jb = i_startblk_e, i_endblk_e
    ! !       CALL get_indices_e(patch_2d, jb, i_startblk_e, i_endblk_e,   &
    ! !         &                start_index, end_index, rl_start,rl_end)
    ! !
    ! !       DO jk = start_level, end_level
    ! !         DO je = start_index, end_index          !
    ! !           IF ( v_base%lsm_e(je,jk,jb) <= sea_boundary ) THEN
    ! !              edge_upwind_flux(je,jk,jb) =  &
    ! !              &  laxfr_upflux( pvn_e(je,jk,jb), pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), &
    ! !              &                                 pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )
    ! ! !write(*,*)'upwind miura',je,jk,jb,edge_upwind_flux(je,jk,jb),edge_upwind_flux(je,jk,jb), pvn_e(je,jk,jb)!,&
    ! ! !&pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2))
    ! !           ELSE
    ! !             edge_upwind_flux(je,jk,jb) = 0.0_wp
    ! !           ENDIF
    ! !         END DO  ! end loop over edges
    ! !       END DO  ! end loop over levels
    ! !     END DO  ! end loop over blocks
    ! ! !---------------------------------------------------------------
    
    !Step3: Local linear subgridcell distribution of tracer is calculated
    !3a: calculate tracer gradient
    !3b: map tracer gradient from edge to cell. Result is gradient vector at cell centers
    !3c: project gradient vector at cell center in direction of vector that points from
    !    cell center to upwind point C_i from step 2
    !3a:
    
    z_gradc(1:nproma,1:n_zlev,1:patch_2d%nblks_e) = 0.0_wp
    
    z_gradc_cc(1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks)%x(1) = 0.0_wp
    z_gradc_cc(1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks)%x(2) = 0.0_wp
    z_gradc_cc(1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks)%x(3) = 0.0_wp
    
    CALL grad_fd_norm_oce_3d( pvar_c,                 &
      & patch_3d,    &
      & p_op_coeff%grad_coeff,  &
      & z_gradc)
    
    CALL sync_patch_array(sync_e, patch_2d, z_gradc)
    !3b:
    CALL map_edges2cell_3d( patch_3d, z_gradc,p_op_coeff, z_gradc_cc)
    
    !!ICON_OMP_PARALLEL
    !!ICON_OMP_DO PRIVATE(start_index, end_index, je, jk, il_c, ib_c)  ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      edge_upwind_flux(:,:,jb) = 0.0_wp
      DO je = start_index, end_index
        DO jk = start_level, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
        
          il_c = p_op_coeff%upwind_cell_idx(je,jk,jb)
          ib_c = p_op_coeff%upwind_cell_blk(je,jk,jb)
          
          edge_upwind_flux(je,jk,jb) = pvn_e(je,jk,jb)&
            & *( pvar_c(il_c, jk, ib_c) +          &
            & DOT_PRODUCT(                    &
            & p_op_coeff%moved_edge_position_cc(je,jk,jb)%x -    &
            & p_op_coeff%upwind_cell_position_cc(je,jk,jb)%x,  &
            & z_gradc_cc(il_c, jk, ib_c)%x))
        END DO
      END DO
    END DO
    !!ICON_OMP_END_DO NOWAIT
    !!ICON_OMP_END_PARALLEL
    
  END SUBROUTINE miura_order1_hflux_oce
  !-------------------------------------------------------------------------

  
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
  !!  mpi note: computed on domain edges. Results is not synced.
  !!
  !<Optimize:inUse>
  SUBROUTINE hflx_limiter_oce_zalesak( patch_3d,        &
    & tracer,              &
    & p_mass_flx_e,      &
    & flx_tracer_low,    &    
    & flx_tracer_high,   &
    & flx_tracer_final,  &
    & p_op_coeff,        &
    & h_old,             &
    & h_new              )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
    REAL(wp), INTENT(inout)             :: tracer           (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: p_mass_flx_e     (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_low   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: flx_tracer_high  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp), INTENT(inout)             :: flx_tracer_final (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    TYPE(t_operator_coeff),INTENT(in)   :: p_op_coeff
    REAL(wp), INTENT(in)                :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    
    !Local variables
    !REAL(wp)              :: flx_tracer_high2  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp) :: z_mflx_anti(patch_3d%p_patch_2d(1)%cells%max_connectivity)
    REAL(wp) :: z_fluxdiv_c     !< flux divergence at cell center
    REAL(wp) :: z_anti          (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)          !< antidiffusive tracer mass flux (F_H - F_L)    
    REAL(wp) :: z_tracer_new_low(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after transport, if low order fluxes are used
    REAL(wp) :: z_tracer_max    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local maximum of current tracer value and low order update
    REAL(wp) :: z_tracer_min    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local minimum of current tracer value and low order update
    REAL(wp) :: r_p             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< fraction which must multiply all in/out fluxes of cell jc to guarantee
    REAL(wp) :: r_m             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< no overshoot/undershoot
    REAL(wp) :: r_frac          !< computed minimum fraction which must multiply< the flux at the edge
    REAL(wp) :: z_min, z_max    !< minimum/maximum value in cell and neighboring cells
    REAL(wp) :: z_signum        !< sign of antidiffusive velocity
    REAL(wp) :: p_p, p_m        !< sum of antidiffusive fluxes into and out of cell jc
    REAL(wp) :: prism_thick_old(n_zlev), inv_prism_thick_new(n_zlev)
    INTEGER, DIMENSION(:,:,:), POINTER ::  cell_of_edge_idx, cell_of_edge_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    INTEGER :: start_level, end_level            
    INTEGER :: start_index, end_index
    INTEGER :: je, jk, jb, jc,  cell_connect
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d    
    !-------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------
    start_level = 1
    end_level   = n_zlev
    cell_of_edge_idx  => patch_2d%edges%cell_idx
    cell_of_edge_blk  => patch_2d%edges%cell_blk
    edge_of_cell_idx  => patch_2d%cells%edge_idx
    edge_of_cell_blk  => patch_2d%cells%edge_blk
    neighbor_cell_idx => patch_2d%cells%neighbor_idx
    neighbor_cell_blk => patch_2d%cells%neighbor_blk
    
    ! 1. Calculate low (first) order fluxes using the standard upwind scheme and the
    !    antidiffusive fluxes
    !    (not allowed to call upwind_hflux_up directly, due to circular dependency)
!     IF ( p_test_run ) THEN
! !       z_tracer_max    (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
! !       z_tracer_min    (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
!       r_m             (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
!       r_p             (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
!     ENDIF

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, je, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      
      z_anti(:,:,jb)     = 0.0_wp
      DO je = start_index, end_index
        DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,jb), end_level)
          
          ! calculate antidiffusive flux for each edge
          z_anti(je,jk,jb) = flx_tracer_high(je,jk,jb) - flx_tracer_low(je,jk,jb)

        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
!ICON_OMP_END_DO
    

#ifdef NAGFOR
 z_tracer_max(:,:,:) = 0.0_wp
 z_tracer_min(:,:,:) = 0.0_wp
 r_m(:,:,:)          = 0.0_wp
 r_p(:,:,:)          = 0.0_wp
#endif
    
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk, inv_prism_thick_new, prism_thick_old, &
!ICON_OMP z_fluxdiv_c ) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      
      z_tracer_new_low(:,:,jb)  = 0.0_wp
      
      DO jc = start_index, end_index
        
        ! get prism thickness
        inv_prism_thick_new(start_level) = 1.0_wp / (patch_3d%p_patch_1d(1)%del_zlev_m(start_level) + h_new(jc,jb))
        prism_thick_old(start_level)     = patch_3d%p_patch_1d(1)%del_zlev_m(start_level)           + h_old(jc,jb)

        DO jk = start_level+1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb), end_level)
          prism_thick_old (jk)    = patch_3d%p_patch_1d(1)%del_zlev_m(jk)
          inv_prism_thick_new(jk) = patch_3d%p_patch_1d(1)%inv_del_zlev_m(jk)
        ENDDO
        
        DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb), end_level)
          
          !  compute also divergence of low order fluxes
!           z_fluxdiv_c =  &
!             & flx_tracer_low(edge_of_cell_idx(jc,jb,1),jk,edge_of_cell_blk(jc,jb,1))   * &
!             & p_op_coeff%div_coeff(jc,jk,jb,1) &
!             & + flx_tracer_low(edge_of_cell_idx(jc,jb,2),jk,edge_of_cell_blk(jc,jb,2)) * &
!             & p_op_coeff%div_coeff(jc,jk,jb,2)  &
!             & + flx_tracer_low(edge_of_cell_idx(jc,jb,3),jk,edge_of_cell_blk(jc,jb,3)) * &
!             & p_op_coeff%div_coeff(jc,jk,jb,3)
          z_fluxdiv_c = 0
          DO cell_connect = 1, patch_2d%cells%num_edges(jc,jb)
            z_fluxdiv_c =  z_fluxdiv_c + &
              & flx_tracer_low(edge_of_cell_idx(jc,jb,cell_connect),jk,edge_of_cell_blk(jc,jb,cell_connect))   * &
              & p_op_coeff%div_coeff(jc,jk,jb,cell_connect)
          ENDDO
          
          !
          ! 3. Compute the updated low order solution z_tracer_new_low
          z_tracer_new_low(jc,jk,jb) = (tracer(jc,jk,jb) * prism_thick_old(jk)     &
            & - dtime * z_fluxdiv_c) * inv_prism_thick_new(jk)
          
          ! precalculate local maximum/minimum of current tracer value and low order
          ! updated value
!            z_tracer_max(jc,jk,jb) = MAX(tracer(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
!            z_tracer_min(jc,jk,jb) = MIN(tracer(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
          
        ENDDO
      ENDDO
      
!       z_tracer_max(start_index:end_index,start_level:end_level,jb) =            &
!         & MAX(          tracer(start_index:end_index,start_level:end_level,jb), &
!         &     z_tracer_new_low(start_index:end_index,start_level:end_level,jb))
!       z_tracer_min(start_index:end_index,start_level:end_level,jb) =            &
!         & MIN(          tracer(start_index:end_index,start_level:end_level,jb), &
!         &     z_tracer_new_low(start_index:end_index,start_level:end_level,jb))
      z_tracer_max(:,:,jb) =            &
        & MAX(          tracer(:,:,jb), &
        &     z_tracer_new_low(:,:,jb))
      z_tracer_min(:,:,jb) =            &
        & MIN(          tracer(:,:,jb), &
        &     z_tracer_new_low(:,:,jb))

!      write(0,*) jb, ":", z_tracer_max(start_index:end_index,start_level:end_level,jb)
!      write(0,*) jb, ":", z_tracer_min(start_index:end_index,start_level:end_level,jb)
    ENDDO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, z_tracer_max, z_tracer_min)
    
    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.    
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk, inv_prism_thick_new, &
!ICON_OMP z_mflx_anti, z_max, z_min, cell_connect, p_p, p_m) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block

      ! this is only needed for the parallel test setups
      ! it will try  tocheck the uninitialized (land) parts
      r_m(:,:,jb) = 0.0_wp
      r_p(:,:,jb) = 0.0_wp
        
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        
        ! get prism thickness
        inv_prism_thick_new(start_level) = 1.0_wp / (patch_3d%p_patch_1d(1)%del_zlev_m(start_level) + h_new(jc,jb))
        DO jk = start_level+1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb), end_level)
          inv_prism_thick_new(jk) = patch_3d%p_patch_1d(1)%inv_del_zlev_m(jk)
        ENDDO
        
        DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb), end_level)
          
          ! 2. Define "antidiffusive" fluxes A(jc,jk,jb,je) for each cell. It is the difference
          !    between the high order fluxes (given by the FFSL-scheme) and the low order
          !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
          !    - positive for outgoing fluxes
          !    - negative for incoming fluxes
          !    this sign convention is related to the definition of the divergence operator.

          
!           z_mflx_anti(1) =                                                        &
!             & dtime * p_op_coeff%div_coeff(jc,jk,jb,1) * inv_prism_thick_new(jk)  &
!             & * z_anti(edge_of_cell_idx(jc,jb,1),jk,edge_of_cell_blk(jc,jb,1))
!           
!           z_mflx_anti(2) =                                                         &
!             & dtime *  p_op_coeff%div_coeff(jc,jk,jb,2) * inv_prism_thick_new(jk)  &
!             & * z_anti(edge_of_cell_idx(jc,jb,2),jk,edge_of_cell_blk(jc,jb,2))
!           
!           z_mflx_anti(3) =                                                         &
!             & dtime * p_op_coeff%div_coeff(jc,jk,jb,3) * inv_prism_thick_new(jk)   &
!             & * z_anti(edge_of_cell_idx(jc,jb,3),jk,edge_of_cell_blk(jc,jb,3))
          
          z_mflx_anti(:) = 0.0_wp
          z_max = z_tracer_max(jc,jk,jb)
          z_min = z_tracer_min(jc,jk,jb)
          p_p = 0.0_wp
          p_m = 0_wp
          DO cell_connect = 1, patch_2d%cells%num_edges(jc,jb)
            IF (patch_3d%p_patch_1d(1)% &
              & dolic_c(neighbor_cell_idx(jc,jb,cell_connect), neighbor_cell_blk(jc,jb,cell_connect)) >= jk) THEN
              
              z_max = MAX(z_max, &
                & z_tracer_max(neighbor_cell_idx(jc,jb,cell_connect),jk,neighbor_cell_blk(jc,jb,cell_connect)))
              z_min = MIN(z_min, &
                & z_tracer_min(neighbor_cell_idx(jc,jb,cell_connect),jk,neighbor_cell_blk(jc,jb,cell_connect)))

              z_mflx_anti(cell_connect) =                                                        &
                & dtime * p_op_coeff%div_coeff(jc,jk,jb,cell_connect) * inv_prism_thick_new(jk)  &
                & * z_anti(edge_of_cell_idx(jc,jb,cell_connect),jk,edge_of_cell_blk(jc,jb,cell_connect))

              ! Sum of all incoming antidiffusive fluxes into cell jc
              p_p = p_p - MIN(0._wp, z_mflx_anti(cell_connect))
              ! Sum of all outgoing antidiffusive fluxes out of cell jc
              p_m = p_m + MAX(0._wp, z_mflx_anti(cell_connect))
              
            ENDIF
          ENDDO
          
!           ! Sum of all incoming antidiffusive fluxes into cell jc
!           p_p = - (  MIN(0._wp,z_mflx_anti(1))   &
!             &      + MIN(0._wp,z_mflx_anti(2))   &
!             &      + MIN(0._wp,z_mflx_anti(3)) )
!           ! Sum of all outgoing antidiffusive fluxes out of cell jc
!           p_m =  MAX(0._wp,z_mflx_anti(1))  &
!             &  + MAX(0._wp,z_mflx_anti(2))  &
!             &  + MAX(0._wp,z_mflx_anti(3))
          
          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of q
          r_m(jc,jk,jb) = (z_tracer_new_low(jc,jk,jb) - z_min ) / (p_m + dbl_eps)
          ! fraction which must multiply all fluxes into cell jc to guarantee no
          ! overshoot
          ! Nominator: maximum allowable increase of q
          r_p(jc,jk,jb) = (z_max - z_tracer_new_low(jc,jk,jb)) / (p_p + dbl_eps)
          
        ENDDO
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    
    ! Synchronize r_m and r_p
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, r_m, r_p)
    
    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !    At the end, compute new, limited fluxes which are then passed to the main
    !    program. Note that flx_tracer_high now denotes the LIMITED flux.
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, je, jk, z_signum, r_frac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      DO je = start_index, end_index
        DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,jb), end_level)
          ! IF( patch_3D%lsm_e(je,jk,jb) <= sea_boundary ) THEN
          z_signum = SIGN(1._wp, z_anti(je,jk,jb))
                    
          ! This does the same as an IF (z_signum > 0) THEN ... ELSE ... ENDIF,
          ! but is computationally more efficient
          r_frac = 0.5_wp * (                                                        &
            & (1._wp + z_signum) *                                           &
            & MIN(r_m(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)),  &
            & r_p(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)))  &
            & +  (1._wp - z_signum) *                                           &
            & MIN(r_m(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)),  &
            & r_p(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)))  )
          
          ! Limited flux
          !flx_tracer_high(je,jk,jb) = flx_tracer_low(je,jk,jb)               &
          !  & + MIN(1._wp,r_frac) * z_anti(je,jk,jb)
            
          flx_tracer_final(je,jk,jb) = flx_tracer_low(je,jk,jb)               &
            & + MIN(1._wp,r_frac) * z_anti(je,jk,jb)            
          ! ELSE
          !   flx_tracer_high(je,jk,jb)= 0.0_wp
          ! ENDIF
        END DO
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    
  END SUBROUTINE hflx_limiter_oce_zalesak
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Positive definite flux limiter for horizontal advection
  !!
  !! Positive definite Zalesak Flux-Limiter (Flux corrected transport).
  !! Only outward fluxes are re-scaled, in order to maintain positive
  !! definiteness.
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !! - Harris, L. M. and P. H. Lauritzen (2010): A flux-form version of the
  !!   Conservative Semi-Lagrangian Multi-tracer transport scheme (CSLAM) on
  !!   the cubed sphere grid. JCP, in press
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2010-10-06)
  !!
   SUBROUTINE hflx_limiter_oce_posdef( patch_3d,&
    & tracer,                             &
    & flx_tracer,                       &    
    & p_op_coeff )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
    REAL(wp), INTENT(inout)             :: tracer      (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: flx_tracer(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
    TYPE(t_operator_coeff),INTENT(in)   :: p_op_coeff

    !Local variables
    REAL(wp) :: z_mflx(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)   
    !REAL(wp) :: r_p             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< fraction which must multiply all in/out fluxes of cell jc to guarantee
    REAL(wp) :: r_m             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< no overshoot/undershoot
    REAL(wp) :: z_signum        !< sign of antidiffusive velocity
    REAL(wp) :: p_m        !< sum of antidiffusive fluxes into and out of cell jc
    INTEGER, DIMENSION(:,:,:), POINTER :: cell_of_edge_idx, cell_of_edge_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    INTEGER :: start_level, end_level            
    INTEGER :: start_index, end_index
    INTEGER :: je, jk, jb, jc!,  cell_connect
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d    
    !-------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------
    start_level = 1
    end_level   = n_zlev
    
    ! Set pointers to index-arrays
    ! line and block indices of two neighboring cells
    cell_of_edge_idx  => patch_2d%edges%cell_idx
    cell_of_edge_blk  => patch_2d%edges%cell_blk
    edge_of_cell_idx  => patch_2d%cells%edge_idx
    edge_of_cell_blk  => patch_2d%cells%edge_blk
    neighbor_cell_idx => patch_2d%cells%neighbor_idx
    neighbor_cell_blk => patch_2d%cells%neighbor_blk
    
    !
    ! 1. Reformulate all fluxes in terms of the total mass [kg m^-3]
    !    that crosses each of the CV-edges and store them in a cell-based structure.
    !
    !    z_mflx > 0: outward
    !    z_mflx < 0: inward
    !    
    !!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk, inv_prism_thick_new, prism_thick_old, &
    !!ICON_OMP z_fluxdiv_c ) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      
      DO jc = start_index, end_index
                
        DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb), end_level)          
        
          z_mflx(jc,jk,1) = p_op_coeff%div_coeff(jc,jk,jb,1) * dtime&
          &*patch_2d%cells%area(jc,jb)*patch_2D%edges%primal_edge_length(edge_of_cell_idx(jc,jb,1),edge_of_cell_blk(jc,jb,1)) &
            &             * flx_tracer(edge_of_cell_idx(jc,jb,1),jk,edge_of_cell_blk(jc,jb,1))

          z_mflx(jc,jk,2) = p_op_coeff%div_coeff(jc,jk,jb,2) * dtime &
          &*patch_2d%cells%area(jc,jb)*patch_2D%edges%primal_edge_length(edge_of_cell_idx(jc,jb,2),edge_of_cell_blk(jc,jb,2)) &          
            &             * flx_tracer(edge_of_cell_idx(jc,jb,2),jk,edge_of_cell_blk(jc,jb,2))
  
          z_mflx(jc,jk,3) = p_op_coeff%div_coeff(jc,jk,jb,3) * dtime &
          &*patch_2d%cells%area(jc,jb)*patch_2D%edges%primal_edge_length(edge_of_cell_idx(jc,jb,3),edge_of_cell_blk(jc,jb,3)) &          
            &             * flx_tracer(edge_of_cell_idx(jc,jb,3),jk,edge_of_cell_blk(jc,jb,3))
          
        ENDDO
      ENDDO
    ENDDO
    !!ICON_OMP_END_DO NOWAIT
    !!ICON_OMP_END_PARALLEL
    
          
    ! 2. Compute total outward mass
    !
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb), end_level)

          ! Sum of all outgoing fluxes out of cell jc
          p_m =  MAX(0._wp,z_mflx(jc,jk,1))  &
            &  + MAX(0._wp,z_mflx(jc,jk,2))  &
            &  + MAX(0._wp,z_mflx(jc,jk,3))

          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of \rho q
          r_m(jc,jk,jb) = MIN(1._wp, (tracer(jc,jk,jb)*patch_3d%p_patch_1d(1)%prism_thick_c(jc,jk,jb)) &
            &                        /(p_m + dbl_eps) )

        ENDDO
      ENDDO
    ENDDO
   ! synchronize r_m
    CALL sync_patch_array(SYNC_C1, patch_2d, r_m)

    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      DO je = start_index, end_index
        DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,jb), end_level)

          ! p_mflx_tracer_h > 0: flux directed from cell 1 -> 2
          ! p_mflx_tracer_h < 0: flux directed from cell 2 -> 1
          z_signum = SIGN(1._wp,flx_tracer(je,jk,jb))

          flx_tracer(je,jk,jb) = flx_tracer(je,jk,jb) * 0.5_wp  &
            & *( (1._wp + z_signum) * r_m(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)) &
            &   +(1._wp - z_signum) * r_m(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)) )

        END DO
      ENDDO
    ENDDO
        
  END SUBROUTINE hflx_limiter_oce_posdef
  
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Lax Friedrichs first order upwind flux,.
  !!
  !! Lax Friedrichs first order upwind flux,
  !! used in conservative advection routines.
  !! For passive advection, equivalent to
  !! any other first order upwind flux.
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura  (2004).
  !!
!<Optimize:inUse>
  FUNCTION laxfr_upflux( p_vn, p_psi1, p_psi2 )  result(p_upflux)
    !
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: p_vn
    REAL(wp), INTENT(in) :: p_psi1, p_psi2
    REAL(wp)             :: p_upflux
    !-----------------------------------------------------------------------
    p_upflux = 0.5_wp * (        p_vn  *( p_psi1 + p_psi2 )    &
      & - ABS( p_vn )*( p_psi2 - p_psi1 ) )
    
  END FUNCTION laxfr_upflux
  
  !-----------------------------------------------------------------------------------------
  ! ! !>
  ! ! !! !  SUBROUTINE advects horizontally the tracers present in the ocean model.
  ! ! !!
  ! ! !! @par Revision History
  ! ! !! Developed  by  Peter Korn, MPI-M (2011).
  ! ! !!
  ! ! SUBROUTINE elad(patch_2d, trac_old, trac_new,p_op_coeff, p_os)
  ! ! !
  ! ! !
  ! ! TYPE(t_patch), TARGET, INTENT(in) :: patch_2d
  ! ! REAL(wp), INTENT(IN)  :: trac_old(:,:,:)
  ! ! REAL(wp), INTENT(inout) :: trac_new(:,:,:)
  ! ! TYPE(t_operator_coeff), INTENT(IN)        :: p_op_coeff
  ! ! TYPE(t_hydro_ocean_state), TARGET :: p_os
  ! ! !3 arrays for explicit part for tracer in Adams-Bashford  stepping,
  ! ! !stores information across different timelevels
  ! ! !REAL(wp) :: trac_out(:,:,:)                              !new tracer
  ! ! !
  ! ! !Local variables
  ! ! INTEGER, PARAMETER :: no_cell_edges = 3
  ! ! !REAL(wp) :: delta_z, delta_z2
  ! ! !INTEGER  :: ctr, ctr_total
  ! ! INTEGER  :: i_startblk_c, i_endblk_c, start_index, end_index, rl_start_c, rl_end_c
  ! ! INTEGER  :: jc, jk, jb!, jkp1        !< index of edge, vert level, block
  ! ! INTEGER  :: z_dolic
  ! ! REAL(wp) :: z_in(nproma,n_zlev,patch_2d%alloc_cell_blocks)
  ! ! REAL(wp) :: z_out(nproma,n_zlev,patch_2d%alloc_cell_blocks)
  ! ! REAL(wp) :: z_pred(nproma,n_zlev,patch_2d%alloc_cell_blocks)
  ! !
  ! ! REAL(wp) :: z_up(nproma,n_zlev,patch_2d%alloc_cell_blocks,no_cell_edges)
  ! ! !REAL(wp) :: z_down(nproma,n_zlev,patch_2d%alloc_cell_blocks,no_cell_edges)
  ! ! REAL(wp) :: z_max(nproma,n_zlev,patch_2d%alloc_cell_blocks)
  ! ! REAL(wp) :: z_min(nproma,n_zlev,patch_2d%alloc_cell_blocks)
  ! ! REAL(wp) :: z_excess(nproma,n_zlev,patch_2d%alloc_cell_blocks)
  ! ! REAL(wp) :: z_grad_excess(nproma,n_zlev,patch_2d%nblks_e)
  ! ! REAL(wp) :: z_diff_excess(nproma,n_zlev,patch_2d%alloc_cell_blocks)
  ! ! REAL(wp) :: z_K(nproma,n_zlev,patch_2d%nblks_e)
  ! ! !REAL(wp) :: max_val, min_val!, dtime2, z_tmp
  ! ! INTEGER  :: il_e(no_cell_edges), ib_e(no_cell_edges)
  ! ! INTEGER  :: il_c1, ib_c1,il_c2, ib_c2, ie
  ! ! INTEGER  :: stop_ctr
  ! ! INTEGER            :: iter
  ! ! INTEGER, PARAMETER :: i_max_iter = 20
  ! ! !-----------------------------------------------------------------------
  ! !   IF (my_process_is_mpi_parallel()) &
  ! !     & CALL finish("elad", "is not mpi parallelized")
  ! !
  ! ! rl_start_c   = 1
  ! ! rl_end_c     = min_rlcell
  ! ! i_startblk_c = patch_2d%cells%start_blk(rl_start_c,1)
  ! ! i_endblk_c   = patch_2d%cells%end_blk(rl_end_c,1)
  ! !
  ! ! z_in       = trac_old
  ! ! z_out      = trac_new
  ! ! z_pred     = trac_new
  ! ! z_excess   = 0.0_wp
  ! ! z_K(:,:,:) = 1.0E12_wp
  ! !
  ! ! ! DO jk=1,n_zlev
  ! ! !   write(*,*)'max/min tracer in:',jk,&
  ! ! !   &maxval(trac_new(:,jk,:)),minval(trac_new(:,jk,:))
  ! ! ! END DO
  ! !
  ! !   !Step 1: determine minimal and maximal permissible values
  ! !   DO jb = i_startblk_c, i_endblk_c
  ! !     CALL get_indices_c( patch_2d, jb, i_startblk_c, i_endblk_c, start_index, end_index, &
  ! !     &                   rl_start_c, rl_end_c)
  ! !
  ! !     DO jc = start_index, end_index
  ! !       z_dolic = v_base%dolic_c(jc,jb)
  ! !       IF(z_dolic>=MIN_DOLIC)THEN
  ! !         DO jk = 1, z_dolic
  ! !           DO ie=1,no_cell_edges
  ! !
  ! !             !actual edges of cell c1
  ! !             il_e(ie) = patch_2d%cells%edge_idx(jc,jb,ie)
  ! !             ib_e(ie) = patch_2d%cells%edge_blk(jc,jb,ie)
  ! !
  ! !             !get neighbor cells of edge
  ! !             il_c1 = patch_2d%edges%cell_idx(il_e(ie),ib_e(ie),1)
  ! !             ib_c1 = patch_2d%edges%cell_blk(il_e(ie),ib_e(ie),1)
  ! !             il_c2 = patch_2d%edges%cell_idx(il_e(ie),ib_e(ie),2)
  ! !             ib_c2 = patch_2d%edges%cell_blk(il_e(ie),ib_e(ie),2)
  ! ! !             IF(z_in(il_c1,jk,ib_c1)>z_in(il_c2,jk,ib_c2))THEN
  ! ! !               z_up(jc,jk,jb,ie)   = z_in(il_c1,jk,ib_c1)
  ! ! !               z_down(jc,jk,jb,ie) = z_in(il_c2,jk,ib_c2)
  ! ! !             ELSE
  ! ! !               z_up(jc,jk,jb,ie)   = z_in(il_c2,jk,ib_c2)
  ! ! !               z_down(jc,jk,jb,ie) = z_in(il_c1,jk,ib_c1)
  ! ! !             ENDIF
  ! !               IF(p_os%p_diag%ptp_vn(il_e(ie),jk,ib_e(ie))>=0.0_wp)THEN
  ! !                 z_up(jc,jk,jb,ie) = trac_old(il_c1,jk,ib_c1)
  ! !               ELSEIF(p_os%p_diag%ptp_vn(il_e(ie),jk,ib_e(ie))<0.0_wp)THEN
  ! !                 z_up(jc,jk,jb,ie) = trac_old(il_c2,jk,ib_c2)
  ! !               ENDIF
  ! !          END DO
  ! !          IF ( v_base%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
  ! !            z_max(jc,jk,jb) = maxval(z_up(jc,jk,jb,1:no_cell_edges))
  ! !            z_min(jc,jk,jb) = minval(z_up(jc,jk,jb,1:no_cell_edges))
  ! !          ELSE
  ! !            z_max(jc,jk,jb) = 0.0_wp
  ! !            z_min(jc,jk,jb) = 0.0_wp
  ! !          ENDIF
  ! ! !         z_max(jc,jk,jb) = maxval(z_up(jc,jk,jb,1:no_cell_edges))
  ! ! !         z_min(jc,jk,jb) = minval(z_down(jc,jk,jb,1:no_cell_edges))
  ! !
  ! !         END DO!level-loop
  ! !       ENDIF!(z_dolic>0)
  ! !     END DO!idx-loop
  ! !   END DO !blk-loop
  ! !
  ! ! ! DO jk=1,n_zlev
  ! ! ! write(*,*)'admissible bound',jk,&
  ! ! ! &maxval(z_max(:,jk,:)), minval(z_max(:,jk,:)),&
  ! ! ! &maxval(z_min(:,jk,:)), minval(z_min(:,jk,:))
  ! ! ! END DO
  ! !
  ! !
  ! ! ITERATION_LOOP: Do iter = 1, i_max_iter
  ! ! !write(*,*)'iteration----------',iter
  ! ! !write(1230,*)'iteration----------',iter
  ! !   !Step 1: determine excess field
  ! !   DO jb = i_startblk_c, i_endblk_c
  ! !     CALL get_indices_c( patch_2d, jb, i_startblk_c, i_endblk_c, start_index, end_index, &
  ! !     &                   rl_start_c, rl_end_c)
  ! !
  ! !     DO jc = start_index, end_index
  ! !       z_dolic = v_base%dolic_c(jc,jb)
  ! !       IF(z_dolic>MIN_DOLIC)THEN
  ! !         DO jk = 1, z_dolic
  ! !          IF ( v_base%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
  ! !            z_excess(jc,jk,jb) = max(z_pred(jc,jk,jb)-z_max(jc,jk,jb),0.0_wp)&
  ! !                               &+min(z_pred(jc,jk,jb)-z_min(jc,jk,jb),0.0_wp)
  ! !           ELSE
  ! !             z_excess(jc,jk,jb) = 0.0_wp
  ! !           ENDIF
  ! !
  ! ! !  IF(z_excess(jc,jk,jb)/=0.0)THEN
  ! ! !  write(1230,*)'excess field',jc,jk,jb,z_excess(jc,jk,jb),&
  ! ! ! &z_pred(jc,jk,jb), z_max(jc,jk,jb),z_up(jc,jk,jb,1:no_cell_edges)
  ! ! !  ENDIF
  ! !         END DO!level-loop
  ! !       ENDIF!(z_dolic>0)
  ! !     END DO!idx-loop
  ! !   END DO !blk-loop
  ! !
  ! ! ! DO jk=1,n_zlev
  ! ! ! write(*,*)'max-min excess',jk,&
  ! ! ! &maxval(z_excess(:,jk,:)), minval(z_excess(:,jk,:))
  ! ! ! END DO
  ! !
  ! !
  ! !    !Step 3: Calculate diffusion of excess field
  ! !    !CALL grad_fd_norm_oce( z_excess, patch_2d, z_grad_excess)
  ! !     CALL grad_fd_norm_oce_3D( z_excess,               &
  ! !            &                  patch_2d,                &
  ! !            &                  p_op_coeff%grad_coeff,  &
  ! !            &                  z_grad_excess)
  ! !    z_grad_excess = z_K*z_grad_excess
  ! !    CALL div_oce_3D( z_grad_excess, patch_2d,p_op_coeff%div_coeff, z_diff_excess)
  ! !
  ! ! ! DO jk=1,n_zlev
  ! ! ! write(*,*)'max-min diffusion',jk,&
  ! ! ! &maxval(z_diff_excess(:,jk,:)), minval(z_diff_excess(:,jk,:))
  ! ! ! END DO
  ! !
  ! !   !Step 4
  ! !   DO jb = i_startblk_c, i_endblk_c
  ! !     CALL get_indices_c( patch_2d, jb, i_startblk_c, i_endblk_c, start_index, end_index, &
  ! !     &                   rl_start_c, rl_end_c)
  ! !
  ! !     DO jc = start_index, end_index
  ! !       z_dolic = v_base%dolic_c(jc,jb)
  ! !       IF(z_dolic>=MIN_DOLIC)THEN
  ! !         DO jk = 1, z_dolic
  ! !           IF ( v_base%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
  ! !             z_out(jc,jk,jb) = z_pred(jc,jk,jb) - z_excess(jc,jk,jb)+ z_diff_excess(jc,jk,jb)
  ! !
  ! ! !  IF(z_excess(jc,jk,jb)/=0.0)THEN
  ! ! !  write(1230,*)'alg', z_out(jc,jk,jb),z_pred(jc,jk,jb),&
  ! ! ! &z_excess(jc,jk,jb), z_diff_excess(jc,jk,jb), &
  ! ! ! &z_max(jc,jk,jb),z_up(jc,jk,jb,1:no_cell_edges)
  ! ! !  ENDIF
  ! !
  ! !           ELSE
  ! !             z_out(jc,jk,jb) = 0.0_wp
  ! !           ENDIF
  ! ! !           IF(z_diff_excess(jc,1,jb)/=0.0_wp)THEN
  ! ! !             write(123,*)'correction',jc,jb,z_in(jc,1,jb),&
  ! ! !             & z_out(jc,1,jb), z_diff_excess(jc,1,jb)
  ! ! !           ENDIF
  ! !         END DO!level-loop
  ! !       ENDIF!(z_dolic>0)
  ! !     END DO!idx-loop
  ! !   END DO !blk-loop
  ! !
  ! !   !Step 4: Stop criterion
  ! !   stop_ctr = 0
  ! !   DO jk=1,n_zlev
  ! !
  ! ! !      write(*,*)'actual state',jk,&
  ! ! !       & maxval(z_out(:,jk,:)),minval(z_out(:,jk,:))
  ! ! !      write(*,*)' pred',jk,&
  ! ! !       & maxval(z_pred(:,jk,:)),minval(z_pred(:,jk,:))
  ! !
  ! !     IF(   maxval(z_excess(:,jk,:))<1.0E-12_wp&
  ! !     &.AND.minval(z_excess(:,jk,:))<1.0E-12_wp)THEN
  ! !
  ! !       stop_ctr = stop_ctr +1
  ! !
  ! !     ELSE
  ! !     ENDIF
  ! !   END DO
  ! !
  ! !   IF(stop_ctr==n_zlev)THEN
  ! !     DO jb = i_startblk_c, i_endblk_c
  ! !       CALL get_indices_c( patch_2d, jb, i_startblk_c, i_endblk_c, start_index, end_index, &
  ! !       &                   rl_start_c, rl_end_c)
  ! !
  ! !       DO jc = start_index, end_index
  ! !         z_dolic = v_base%dolic_c(jc,jb)
  ! !         IF(z_dolic>=MIN_DOLIC)THEN
  ! !           DO jk = 1, z_dolic
  ! !            IF ( v_base%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
  ! !            trac_new(jc,jk,jb) = z_out(jc,jk,jb)
  ! !            ENDIF
  ! ! !            write(1230,*)'before-after',jc,jk,jb,&
  ! ! !           &trac_new(jc,jk,jb),trac_old(jc,jk,jb),z_excess(jc,jk,jb)
  ! !           END DO
  ! !         ENDIF
  ! !        ENDDO
  ! !      END DO
  ! !       exit ITERATION_LOOP
  ! !    ELSE
  ! !        z_pred(:,:,:) = z_out(:,:,:)
  ! !   ENDIF
  ! ! END DO ITERATION_LOOP
  ! !
  ! ! !   DO jb = i_startblk_c, i_endblk_c
  ! ! !     CALL get_indices_c( patch_2d, jb, i_startblk_c, i_endblk_c, start_index, end_index, &
  ! ! !     &                   rl_start_c, rl_end_c)
  ! ! !     DO jc = start_index, end_index
  ! ! ! write(123,*)'trac old new',trac_old(jc,1,jb),z_old_new(jc,1,jb),&
  ! ! ! &trac_new(jc,1,jb), z_out(jc,1,jb)
  ! ! !
  ! ! !    END DO
  ! ! !  END DO
  ! ! END SUBROUTINE elad
  ! !   !-------------------------------------------------------------------------
   !-----------------------------------------------------------------------
! !   !>
! !   !! !  SUBROUTINE advects horizontally the tracers present in the ocean model.
! !   !!
! !   !! @par Revision History
! !   !! Developed  by  Peter Korn, MPI-M (2011).
! !   !!
! !   ! Not used currently, no openmp
! !   SUBROUTINE advect_diffuse_horz_high_res( p_patch_3d,          &
! !     & trac_old,            &
! !     & p_os,                &
! !     & p_op_coeff,          &
! !     & k_h,                 &
! !     & h_old,               &
! !     & h_new,               &
! !     & flux_horz)
! !     
! !     TYPE(t_patch_3d ),TARGET, INTENT(in)   :: p_patch_3d
! !     REAL(wp)                               :: trac_old(1:nproma,1:n_zlev,1:p_patch_3d%p_patch_2d(1)%nblks_c)
! !     TYPE(t_hydro_ocean_state), TARGET :: p_os
! !     TYPE(t_operator_coeff), INTENT(in)     :: p_op_coeff
! !     REAL(wp), INTENT(in)                   :: k_h(:,:,:)         !horizontal mixing coeff
! !     REAL(wp), INTENT(in)                   :: h_old(1:nproma,1:p_patch_3d%p_patch_2d(1)%nblks_c)
! !     REAL(wp), INTENT(in)                   :: h_new(1:nproma,1:p_patch_3d%p_patch_2d(1)%nblks_c)
! !     REAL(wp), INTENT(out)                  :: flux_horz(1:nproma,1:n_zlev,1:p_patch_3d%p_patch_2d(1)%nblks_c)
! !     !
! !     !Local variables
! !     !REAL(wp) :: delta_z
! !     INTEGER :: i_startidx_c, i_endidx_c
! !     INTEGER :: i_startidx_e, i_endidx_e
! !     INTEGER :: jc, jk, jb, je,ictr,ictr1
! !     INTEGER :: i_edge, ii_edge, ii_e, ib_e, iii_e, iib_e
! !     INTEGER :: ii_c1, ib_c1, ii_c2, ib_c2,ii_c3, ib_c3
! !     INTEGER :: ii_e1, ib_e1, ii_e2, ib_e2,ii_e3, ib_e3
! !     INTEGER :: in_idx(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e,2)
! !     INTEGER :: in_blk(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e,2)
! !     INTEGER :: out_idx(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e,2)
! !     INTEGER :: out_blk(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e,2)
! !     
! !     REAL(wp) :: z_adv_flux_h (nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e)  ! horizontal advective tracer flux
! !     REAL(wp) :: z_adv_flux_high (nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e)  ! horizontal advective tracer flux
! !     REAL(wp) :: z_adv_flux_low (nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e)  ! horizontal advective tracer flux
! !     REAL(wp) :: z_adv_flux_h2 (nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e)  ! horizontal advective tracer flux
! !     REAL(wp) :: z_diff_trac, z_sum_flux_diff,z_sum_tmp_mflux,z_tmp_flux, z_tmp
! !     !REAL(wp) :: z_difference_h(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)  ! horizontal advective tracer flux
! !     REAL(wp) :: z_div_adv_h  (nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_c)   ! horizontal tracer divergence
! !     REAL(wp) :: z_div_adv_h2  (nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_c)   ! horizontal tracer divergence
! !     
! !     REAL(wp) :: z_div_diff_h (nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_c)  ! horizontal tracer divergence
! !     REAL(wp) :: z_diff_flux_h(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e) ! horizontal diffusive tracer flux
! !     
! !     REAL(wp) :: z_tflux_out(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e)
! !     REAL(wp) :: z_tflux_in(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e)
! !     REAL(wp) :: z_mflux(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e)
! !     
! !     REAL(wp) :: z_adv_plus_diff_flux(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e)
! !     REAL(wp) :: z_consec_grad(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e)
! !     REAL(wp) :: z_consec_grad2(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e)
! !     REAL(wp) :: z_limit_phi(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e)
! !     REAL(wp) :: z_limit_psi(nproma,n_zlev,p_patch_3d%p_patch_2d(1)%nblks_e)
! !     TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
! !     TYPE(t_patch), POINTER :: p_patch
! !     !-------------------------------------------------------------------------------
! !     p_patch         => p_patch_3d%p_patch_2d(1)
! !     edges_in_domain => p_patch%edges%in_domain
! !     cells_in_domain => p_patch%cells%in_domain
! !     !-------------------------------------------------------------------------------
! !     !z_vn         (1:nproma,1:n_zlev,1:p_patch%nblks_e)=0.0_wp
! !     z_adv_flux_h (1:nproma,1:n_zlev,1:p_patch%nblks_e)=0.0_wp
! !     z_adv_flux_high (1:nproma,1:n_zlev,1:p_patch%nblks_e)=0.0_wp
! !     z_adv_flux_low (1:nproma,1:n_zlev,1:p_patch%nblks_e)=0.0_wp
! !     z_adv_flux_h2 (1:nproma,1:n_zlev,1:p_patch%nblks_e)=0.0_wp
! !     !z_difference_h(1:nproma,1:n_zlev,1:p_patch%nblks_e)=0.0_wp
! !     z_div_adv_h  (1:nproma,1:n_zlev,1:p_patch%nblks_c)=0.0_wp
! !     z_div_adv_h2  (1:nproma,1:n_zlev,1:p_patch%nblks_c)=0.0_wp
! !     z_div_diff_h (1:nproma,1:n_zlev,1:p_patch%nblks_c)=0.0_wp
! !     z_diff_flux_h(1:nproma,1:n_zlev,1:p_patch%nblks_e)=0.0_wp
! !     flux_horz    (1:nproma,1:n_zlev,1:p_patch%nblks_c)=0.0_wp
! !     
! !     in_idx(1:nproma,1:n_zlev,1:p_patch%nblks_e,1:2) =0
! !     in_blk(1:nproma,1:n_zlev,1:p_patch%nblks_e,1:2) =0
! !     out_idx(1:nproma,1:n_zlev,1:p_patch%nblks_e,1:2)=0
! !     out_blk(1:nproma,1:n_zlev,1:p_patch%nblks_e,1:2)=0
! !     
! !     z_tflux_out(1:nproma,1:n_zlev,1:p_patch%nblks_e)  = 0.0_wp
! !     z_tflux_in(1:nproma,1:n_zlev,1:p_patch%nblks_e)   = 0.0_wp
! !     z_mflux(1:nproma,1:n_zlev,1:p_patch%nblks_e)     = 0.0_wp
! !     
! !     z_consec_grad(1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp
! !     z_consec_grad2(1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp
! !     z_limit_phi(1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp
! !     z_limit_psi(1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp
! !     z_adv_plus_diff_flux(1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp
! !     !Does not make sense to implement this as an option to the
! !     !available (horz/vert) transport routines!
! !     !Its still possible to come out of the sbr with a
! !     !horz flux and leave the caller sbr on oce-tracer unchanged.
! !     !
! !     !Keep in mind mass-tracer consistency (what happens to identical constant tracer)
! !     !in this case the "r" below equals 0 and psi should become 0  as well (this constitutes a test)
! !     !
! !     !this sbr changes the effective diffusion (small and captial D in Cassulli-Zanolli)
! !     !
! !     !horizontal:
! !     !input data for limiter sbr: grid, tracer grad, fluxes Q
! !     !output of limiter sbr: a) the term psi times abs(Q)times grad tracer
! !     !                       b) effective diffusion coeff
! !     !write(*,*)'INSIDE input tracer',maxval(trac_old(:,1,:)), minval(trac_old(:,1,:))
! !     !      CALL upwind_hflux_oce( p_patch_3D,                   &
! !     !                             & trac_old,                     &
! !     !                             & p_os%p_diag%vn_time_weighted, &
! !     !                             & z_adv_flux_h2 )
! !     ! !      CALL central_hflux_oce( p_patch_3D,                   &
! !     ! !                             & trac_old,                     &
! !     ! !                             & p_os%p_diag%vn_time_weighted, &
! !     ! !                             & z_adv_flux_h2 )
! !     ! CALL miura_order1_hflux_oce( p_patch_3D,                    &
! !     !                              & trac_old,                     &
! !     !                              & p_os%p_diag%vn_time_weighted, &
! !     !                              & p_op_coeff,                   &
! !     !                              & z_adv_flux_h2 )
! !     !
! !     !
! !     !  DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !     !       CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
! !     !       DO jk = 1, n_zlev
! !     !         DO je = i_startidx_e, i_endidx_e
! !     !           z_adv_flux_h2(je,jk,jb) = z_adv_flux_h2(je,jk,jb)*p_patch_3D%p_patch_1D(1)%prism_thick_e(je,jk,jb)
! !     !         END DO
! !     !       END DO
! !     !     END DO
! !     !     CALL div_oce_3D( z_adv_flux_h2, p_patch,p_op_coeff%div_coeff, z_div_adv_h2,&
! !     !     & subset_range=cells_in_domain)
! !     
! !     !1) calculate phi-part of limiter
! !     DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !       CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
! ! #ifdef __SX__
! ! !CDIR UNROLL=6
! ! #endif
! !       DO jk = 1, n_zlev
! !         DO je = i_startidx_e, i_endidx_e
! !           
! !           !Determine edges with inflow (can later be done in prep-sbr).
! !           !This depends on sign of velocity: positive = outflow, negtive = inflow
! !           
! !           !get adjacent cells and calculate tracer differrence
! !           !for later use.
! !           ii_c1 = p_patch%edges%cell_idx(je,jb,1)
! !           ib_c1 = p_patch%edges%cell_blk(je,jb,1)
! !           
! !           ii_c2 = p_patch%edges%cell_idx(je,jb,2)
! !           ib_c2 = p_patch%edges%cell_blk(je,jb,2)
! !           
! !           !z_difference_h(je,jk,jb)=(trac_old(ii_c2,jk,ib_c2)-trac_old(ii_c1,jk,ib_c1))
! !           z_diff_flux_h(je,jk,jb)=trac_old(ii_c2,jk,ib_c2)-trac_old(ii_c1,jk,ib_c1)
! !           
! !           !outflow case
! !           IF(p_os%p_diag%vn_time_weighted(je,jk,jb)>=0.0)THEN
! !             
! !             z_tflux_out(je,jk,jb) = p_os%p_diag%vn_time_weighted(je,jk,jb)*trac_old(ii_c1,jk,ib_c1)
! !             z_mflux(je,jk,jb)     = p_os%p_diag%vn_time_weighted(je,jk,jb)
! !             
! !             !inflow case
! !           ELSEIF(p_os%p_diag%vn_time_weighted(je,jk,jb)<0.0)THEN
! !             
! !             z_tflux_in(je,jk,jb) =-p_os%p_diag%vn_time_weighted(je,jk,jb)*trac_old(ii_c2,jk,ib_c2)
! !             z_mflux(je,jk,jb)    =-p_os%p_diag%vn_time_weighted(je,jk,jb)
! !             
! !           ENDIF
! !           z_adv_flux_low(je,jk,jb) =(z_tflux_out(je,jk,jb)-z_tflux_in(je,jk,jb))
! !           
! !           z_adv_flux_high(je,jk,jb) = z_mflux(je,jk,jb)*z_diff_flux_h(je,jk,jb)
! !           
! !           z_tmp = 2.0_wp*k_h(je,jk,jb)/p_patch%edges%dual_edge_length(je,jb)
! !           
! !           !This corresponds to eq. (14) in Casulli-Zanolli
! !           z_limit_phi(je,jk,jb) = MIN(1.0_wp,z_tmp/z_mflux(je,jk,jb))
! !           
! !           !write(345,*)'phi-limit',je,jk,jb, z_limit_phi(je,jk,jb), z_tmp, z_mflux(je,jk,jb)
! !           !_adv_flux_h(je,jk,jb)=&
! !           !&z_flux_out(je,jk,jb)-z_flux_in(je,jk,jb)&
! !           !&+z_mflux(je,jk,jb)*z_diff_h(je,jk,jb)
! !           !z_adv_flux_h(je,jk,jb)=z_adv_flux_h(je,jk,jb) * p_patch_3D%p_patch_1D(1)%prism_thick_e(je,jk,jb)
! !         END DO
! !       END DO
! !     END DO
! !     
! !     CALL ratio_consecutive_gradients(p_patch_3d,p_os%p_diag%vn_time_weighted,trac_old,z_consec_grad)
! !     
! !     !3) calculate psi-part of limiter
! !     CALL calculate_limiter(p_patch_3d,z_limit_phi,z_consec_grad,z_limit_psi)
! !     
! !     !ictr=0
! !     !ictr1=0
! !     !4) Calculate limited advective flux
! !     DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !       CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
! !       DO jk = 1, n_zlev
! !         DO je = i_startidx_e, i_endidx_e
! !           
! !           !only low gives upwind, high+low without limiter gives central
! !           z_adv_flux_h(je,jk,jb) =           &
! !             & z_adv_flux_low (je,jk,jb)&
! !             & +0.5_wp*z_adv_flux_high(je,jk,jb)&
! !             & *(z_limit_phi(je,jk,jb)-z_limit_psi(je,jk,jb))
! !           
! !           z_diff_flux_h(je,jk,jb) = &
! !             & z_diff_flux_h(je,jk,jb) &
! !             & *(k_h(je,jk,jb)/p_patch%edges%dual_edge_length(je,jb))
! !           
! !           z_adv_plus_diff_flux(je,jk,jb)=&
! !             & (z_diff_flux_h(je,jk,jb)-z_adv_flux_h(je,jk,jb))&
! !             & *p_patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,jb)
! !           
! !         END DO
! !       END DO
! !     END DO
! !     ! write(*,*)'max-min LIMITER',&
! !     ! maxval(0.5_wp*(z_limit_phi+z_limit_psi)),minval(0.5_wp*(z_limit_phi+z_limit_psi)),ictr,ictr1,&
! !     ! &  (edges_in_domain%end_block-edges_in_domain%start_block+1)*(i_endidx_e-i_startidx_e+1)
! !     !write(*,*)'max-min ADV-DIFF-SUM',&
! !     !maxval(z_adv_flux_h),minval(z_adv_flux_h),  &
! !     !maxval(z_diff_flux_h),minval(z_diff_flux_h),&
! !     !maxval(z_adv_plus_diff_flux),minval(z_adv_plus_diff_flux)
! !     
! !     !---------DEBUG DIAGNOSTICS-------------------------------------------
! !     idt_src=4  ! output print level (1-5, fix)
! !     CALL dbg_print('AdvDifHorzHR: adv_flux_h', z_adv_flux_h, str_module,idt_src, p_patch_3d%p_patch_2d(1)%edges%owned)
! !     CALL dbg_print('AdvDifHorzHR: adv+diff_h', z_adv_flux_h, str_module,idt_src, p_patch_3d%p_patch_2d(1)%edges%owned)
! !     idt_src=5  ! output print level (1-5, fix)
! !     CALL dbg_print('AdvDifHorzHR: limit_phi ', z_limit_phi,  str_module,idt_src, p_patch_3d%p_patch_2d(1)%edges%owned)
! !     CALL dbg_print('AdvDifHorzHR: limit_psi ', z_limit_psi,  str_module,idt_src, p_patch_3d%p_patch_2d(1)%edges%owned)
! !     !---------------------------------------------------------------------
! !     
! !     
! !     
! !     ! write(*,*)'max-min LIMITER PHI-PSI-CONSEC-GRAD',&
! !     ! maxval(z_limit_phi),minval(z_limit_phi),    &
! !     ! maxval(z_limit_psi),minval(z_limit_psi),&
! !     ! maxval(z_consec_grad),minval(z_consec_grad)
! !     !4) calculate divergence of advective and diffusive fluxes
! !     !CALL div_oce_3D( z_adv_flux_h, p_patch,p_op_coeff%div_coeff, z_div_adv_h,&
! !     ! & subset_range=cells_in_domain)
! !     !CALL div_oce_3D( z_diff_flux_h, p_patch,p_op_coeff%div_coeff, z_div_diff_h,&
! !     ! & subset_range=cells_in_domain)
! !     
! !     
! !     CALL div_oce_3d( z_adv_plus_diff_flux, p_patch,p_op_coeff%div_coeff, flux_horz,&
! !       & subset_range=cells_in_domain)
! !     
! !   END SUBROUTINE advect_diffuse_horz_high_res
! !   !-----------------------------------------------------------------------
! !  
! !   !-------------------------------------------------------------------------
! !   !>
! !   !! Flux limiter for horizontal advection
! !   !!
! !   !! Zalesak Flux-Limiter (Flux corrected transport)
! !   !! The corrected flux is a weighted average of the low order flux and the
! !   !! given high order flux. The high order flux is used to the greatest extent
! !   !! possible without introducing overshoots and undershoots.
! !   !! Note: This limiter is positive definite and almost monotone (but not strictly).
! !   !!
! !   !! @par Literature:
! !   !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
! !   !!   Algorithms for Fluids. JCP, 31, 335-362
! !   !!
! !   !! @par Revision History
! !   !! - Inital revision by Daniel Reinert, DWD (2010-03-10)
! !   !! Modification by Daniel Reinert, DWD (2010-03-25)
! !   !! - adapted for MPI parallelization
! !   !!
! !   !!  mpi note: computed on domain edges. Results is not synced.
! !   !!
! !   SUBROUTINE hflx_limiter_oce_mo( patch_3d,        &
! !     & tracer,              &
! !     & p_mass_flx_e,      &
! !     & p_mflx_tracer_h,   &
! !     & p_op_coeff,        &
! !     & h_old,             &
! !     & h_new,             &
! !     & opt_start_level, opt_end_level )
! !     
! !     TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
! !     REAL(wp), INTENT(inout)           :: tracer             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks) !< advected cell centered variable
! !     REAL(wp), INTENT(inout)           :: p_mass_flx_e     (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) !< horizontal mass flux
! !     REAL(wp), INTENT(inout)           :: p_mflx_tracer_h  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) !< calculated horizontal tracer mass flux
! !     TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
! !     REAL(wp), INTENT(in)              :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
! !     REAL(wp), INTENT(in)              :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
! !     INTEGER, INTENT(in), OPTIONAL :: opt_start_level !< optional vertical start level
! !     INTEGER, INTENT(in), OPTIONAL :: opt_end_level !< optional vertical end level
! !     
! !     
! !     !Local variables
! !     REAL(wp) :: z_mflx_low      (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)      !< first order tracer mass flux
! !     REAL(wp) :: z_anti          (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)          !< antidiffusive tracer mass flux (F_H - F_L)
! !     REAL(wp) :: z_mflx_anti     (3)
! !     REAL(wp) :: z_fluxdiv_c        !< flux divergence at cell center
! !     REAL(wp) :: z_tracer_new_low(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after hor. transport,
! !     !< if the low order fluxes are used
! !     REAL(wp) :: z_tracer_max(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    !< local maximum of current tracer value and low
! !     !< order update
! !     REAL(wp) :: z_tracer_min(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    !< local minimum of current tracer value and low
! !     !< order update
! !     REAL(wp) :: r_p(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)             !< fraction which must multiply all in/out fluxes
! !     !< of cell jc to guarantee
! !     REAL(wp) :: r_m(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)             !< no overshoot/undershoot
! !     REAL(wp) :: r_frac                                         !< computed minimum fraction which must multiply< the flux at the edge
! !     REAL(wp) :: z_min, z_max                         !< minimum/maximum value in cell and neighboring cells
! !     REAL(wp) :: z_signum                                       !< sign of antidiffusive velocity
! !     REAL(wp) :: p_p, p_m                                       !< sum of antidiffusive fluxes into and out of cell jc
! !     REAL(wp) :: prism_thick_old(n_zlev), inv_prism_thick_new(n_zlev)
! !     INTEGER, DIMENSION(:,:,:), POINTER ::  cell_of_edge_idx, cell_of_edge_blk   !< Pointer to line and block indices of two
! !     !< neighbor cells (array)
! !     INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk  !< Pointer to line and block indices of three
! !     !< neighbor cells (array)
! !     INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk    !< Pointer to line and block indices (array)
! !     !< of edges
! !     INTEGER :: start_level, end_level             !< vertical start and end level
! !     INTEGER :: start_index, end_index
! !     INTEGER :: je, jk, jb, jc,  cell_connect        !< index of edge, vert level, block, cell
! !     TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
! !     TYPE(t_patch), POINTER :: patch_2d
! !     
! !     !-------------------------------------------------------------------------
! !     patch_2d        => patch_3d%p_patch_2d(1)
! !     edges_in_domain => patch_2d%edges%in_domain
! !     cells_in_domain => patch_2d%cells%in_domain
! !     !-------------------------------------------------------------------------
! !     ! Check for optional arguments
! !     IF ( PRESENT(opt_start_level) ) THEN
! !       start_level = opt_start_level
! !     ELSE
! !       start_level = 1
! !     END IF
! !     IF ( PRESENT(opt_end_level) ) THEN
! !       end_level = opt_end_level
! !     ELSE
! !       end_level = n_zlev
! !     END IF
! !     
! !     ! Set pointers to index-arrays
! !     ! line and block indices of two neighboring cells
! !     cell_of_edge_idx => patch_2d%edges%cell_idx
! !     cell_of_edge_blk => patch_2d%edges%cell_blk
! !     ! line and block indices of edges as seen from cells
! !     edge_of_cell_idx => patch_2d%cells%edge_idx
! !     edge_of_cell_blk => patch_2d%cells%edge_blk
! !     ! pointers to line and block indices of three neighbor cells
! !     neighbor_cell_idx => patch_2d%cells%neighbor_idx
! !     neighbor_cell_blk => patch_2d%cells%neighbor_blk
! !     !
! !     ! 1. Calculate low (first) order fluxes using the standard upwind scheme and the
! !     !    antidiffusive fluxes
! !     !    (not allowed to call upwind_hflux_up directly, due to circular dependency)
! !     IF ( p_test_run ) THEN
! !       ! z_tracer_new_low(1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
! !       z_tracer_max    (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
! !       z_tracer_min    (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
! !       r_m             (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
! !       r_p             (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
! !       !    z_anti          (1:nproma,1:n_zlev,1:patch_2d%nblks_e) = 0.0_wp
! !       !    z_mflx_low      (1:nproma,1:n_zlev,1:patch_2d%nblks_e) = 0.0_wp
! !     ENDIF
! !     
! !     !ICON_OMP_PARALLEL
! !     !ICON_OMP_DO PRIVATE(start_index, end_index, je, jk) ICON_OMP_DEFAULT_SCHEDULE
! !     DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !       CALL get_index_range(edges_in_domain, jb, start_index, end_index)
! !       z_mflx_low(:,:,jb) = 0.0_wp
! !       z_anti(:,:,jb)     = 0.0_wp
! !       DO je = start_index, end_index
! !         DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,jb), end_level)
! !           !IF( patch_3D%lsm_e(je,jk,jb) <= sea_boundary ) THEN
! !           z_mflx_low(je,jk,jb) = &
! !             & laxfr_upflux( p_mass_flx_e(je,jk,jb), &
! !             & tracer(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)), &
! !             & tracer(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)) )
! !           
! !           ! calculate antidiffusive flux for each edge
! !           z_anti(je,jk,jb) = p_mflx_tracer_h(je,jk,jb) - z_mflx_low(je,jk,jb)
! !           !ENDIF
! !         END DO  ! end loop over edges
! !       END DO  ! end loop over levels
! !     END DO  ! end loop over blocks
! !     !ICON_OMP_END_DO
! !     
! !     ! no need to sync since we compute on cells_in_domain
! !     ! CALL sync_patch_array(SYNC_E, patch_2d,z_mflx_low)
! !     ! CALL sync_patch_array(SYNC_E, patch_2d, z_anti)
! !     
! !     !ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk, inv_prism_thick_new, prism_thick_old, &
! !     !ICON_OMP z_fluxdiv_c ) ICON_OMP_DEFAULT_SCHEDULE
! !     DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !       CALL get_index_range(cells_in_domain, jb, start_index, end_index)
! !       
! !       ! the number of levels are checked when using them
! !       ! so no zeroeing is required
! !       ! z_tracer_max    (:,:,jb) = 0.0_wp
! !       ! z_tracer_min    (:,:,jb) = 0.0_wp
! !       
! !       DO jc = start_index, end_index
! !         
! !         ! get prism thickness
! !         inv_prism_thick_new(start_level) = 1.0_wp / (patch_3d%p_patch_1d(1)%del_zlev_m(start_level) + h_new(jc,jb))
! !         prism_thick_old(start_level)     = patch_3d%p_patch_1d(1)%del_zlev_m(start_level)           + h_old(jc,jb)
! !         DO jk = start_level+1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb), end_level)
! !           prism_thick_old (jk)    = patch_3d%p_patch_1d(1)%del_zlev_m(jk)
! !           !            inv_prism_thick_new = 1.0_wp / patch_3D%p_patch_1D(1)%del_zlev_m(jk) ! should be calclulated only once
! !           inv_prism_thick_new(jk) = patch_3d%p_patch_1d(1)%inv_del_zlev_m(jk)
! !         ENDDO
! !         
! !         
! !         DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb), end_level)
! !           
! !           ! IF( patch_3D%lsm_c(jc,jk,jb) > sea_boundary ) &
! !           !& CALL finish("","patch_3D%lsm_c(jc,jk,jb) > sea_boundar")
! !           
! !           !  compute also divergence of low order fluxes
! !           z_fluxdiv_c =  &
! !             & z_mflx_low(edge_of_cell_idx(jc,jb,1),jk,edge_of_cell_blk(jc,jb,1)) * &
! !             & p_op_coeff%div_coeff(jc,jk,jb,1) &
! !             & + z_mflx_low(edge_of_cell_idx(jc,jb,2),jk,edge_of_cell_blk(jc,jb,2)) * &
! !             & p_op_coeff%div_coeff(jc,jk,jb,2)  &
! !             & + z_mflx_low(edge_of_cell_idx(jc,jb,3),jk,edge_of_cell_blk(jc,jb,3)) * &
! !             & p_op_coeff%div_coeff(jc,jk,jb,3)
! !           
! !           !
! !           ! 3. Compute the updated low order solution z_tracer_new_low
! !           ! IF( patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
! !           z_tracer_new_low(jc,jk,jb) = (tracer(jc,jk,jb) * prism_thick_old(jk)     &
! !             & - dtime * z_fluxdiv_c) * inv_prism_thick_new(jk)
! !           
! !           ! precalculate local maximum/minimum of current tracer value and low order
! !           ! updated value
! !           z_tracer_max(jc,jk,jb) = MAX(tracer(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
! !           z_tracer_min(jc,jk,jb) = MIN(tracer(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
! !           !          ELSE
! !           !            z_tracer_new_low(jc,jk,jb) = 0.0_wp
! !           !            z_tracer_max(jc,jk,jb)     = 0.0_wp
! !           !            z_tracer_min(jc,jk,jb)     = 0.0_wp
! !           !          ENDIF
! !           
! !         ENDDO
! !       ENDDO
! !     ENDDO
! !     !ICON_OMP_END_DO NOWAIT
! !     !ICON_OMP_END_PARALLEL
! !     
! !     ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
! !     !    field is free of any new extrema.
! !     CALL sync_patch_array_mult(sync_c1, patch_2d, 2, z_tracer_max, z_tracer_min)
! !     
! !     !ICON_OMP_PARALLEL
! !     !ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk, inv_prism_thick_new, &
! !     !ICON_OMP z_mflx_anti, z_max, z_min, cell_connect, p_p, p_m) ICON_OMP_DEFAULT_SCHEDULE
! !     DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !       
! !       ! The dolic_e is the min of the neigboring dolic_c, so
! !       ! zeroeing of r_m, r_p is not required
! !       ! r_m(:,:,jb) = 0.0_wp
! !       ! r_p(:,:,jb) = 0.0_wp
! !       
! !       CALL get_index_range(cells_in_domain, jb, start_index, end_index)
! !       DO jc = start_index, end_index
! !         
! !         ! get prism thickness
! !         inv_prism_thick_new(start_level) = 1.0_wp / (patch_3d%p_patch_1d(1)%del_zlev_m(start_level) + h_new(jc,jb))
! !         ! prism_thick_old(start_level)     = patch_3d%p_patch_1d(1)%del_zlev_m(start_level)           + h_old(jc,jb)
! !         DO jk = start_level+1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb), end_level)
! !           ! prism_thick_old (jk)    = patch_3d%p_patch_1d(1)%del_zlev_m(jk)
! !           !            inv_prism_thick_new = 1.0_wp / patch_3D%p_patch_1D(1)%del_zlev_m(jk) ! should be calclulated only once
! !           inv_prism_thick_new(jk) = patch_3d%p_patch_1d(1)%inv_del_zlev_m(jk)
! !         ENDDO
! !         
! !         DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb), end_level)
! !           
! !           ! 2. Define "antidiffusive" fluxes A(jc,jk,jb,je) for each cell. It is the difference
! !           !    between the high order fluxes (given by the FFSL-scheme) and the low order
! !           !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
! !           !    - positive for outgoing fluxes
! !           !    - negative for incoming fluxes
! !           !    this sign convention is related to the definition of the divergence operator.
! !           
! !           ! IF( patch_3D%lsm_c(jc,jk,jb) > sea_boundary ) &
! !           !& CALL finish("","patch_3D%lsm_c(jc,jk,jb) > sea_boundar")
! !           
! !           z_mflx_anti(1) =                                                        &
! !             & dtime * p_op_coeff%div_coeff(jc,jk,jb,1) * inv_prism_thick_new(jk)  &
! !             & * z_anti(edge_of_cell_idx(jc,jb,1),jk,edge_of_cell_blk(jc,jb,1))
! !           
! !           z_mflx_anti(2) =                                                         &
! !             & dtime *  p_op_coeff%div_coeff(jc,jk,jb,2) * inv_prism_thick_new(jk)  &
! !             & * z_anti(edge_of_cell_idx(jc,jb,2),jk,edge_of_cell_blk(jc,jb,2))
! !           
! !           z_mflx_anti(3) =                                                         &
! !             & dtime * p_op_coeff%div_coeff(jc,jk,jb,3) * inv_prism_thick_new(jk)   &
! !             & * z_anti(edge_of_cell_idx(jc,jb,3),jk,edge_of_cell_blk(jc,jb,3))
! !           
! !           !  IF( patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
! !           !          ! max value of cell and its neighbors
! !           !          ! also look back to previous time step
! !           !          z_max = MAX( z_tracer_max(jc,jk,jb),                          &
! !           !            & z_tracer_max(neighbor_cell_idx(jc,jb,1),jk,neighbor_cell_blk(jc,jb,1)),  &
! !           !            & z_tracer_max(neighbor_cell_idx(jc,jb,2),jk,neighbor_cell_blk(jc,jb,2)),  &
! !           !            & z_tracer_max(neighbor_cell_idx(jc,jb,3),jk,neighbor_cell_blk(jc,jb,3)) )
! !           !          ! min value of cell and its neighbors
! !           !          ! also look back to previous time step
! !           !          z_min = MIN( z_tracer_min(jc,jk,jb),                          &
! !           !            & z_tracer_min(neighbor_cell_idx(jc,jb,1),jk,neighbor_cell_blk(jc,jb,1)),  &
! !           !            & z_tracer_min(neighbor_cell_idx(jc,jb,2),jk,neighbor_cell_blk(jc,jb,2)),  &
! !           !            & z_tracer_min(neighbor_cell_idx(jc,jb,3),jk,neighbor_cell_blk(jc,jb,3)) )
! !           !  ENDIF
! !           z_max = z_tracer_max(jc,jk,jb)
! !           z_min = z_tracer_min(jc,jk,jb)
! !           DO cell_connect = 1, 3
! !             !            IF (neighbor_cell_idx(jc,jb,cell_connect) > 0) THEN
! !             IF (patch_3d%p_patch_1d(1)% &
! !               & dolic_c(neighbor_cell_idx(jc,jb,cell_connect), neighbor_cell_blk(jc,jb,cell_connect)) >= jk) THEN
! !               z_max = MAX(z_max, &
! !                 & z_tracer_max(neighbor_cell_idx(jc,jb,cell_connect),jk,neighbor_cell_blk(jc,jb,cell_connect)))
! !               z_min = MIN(z_min, &
! !                 & z_tracer_min(neighbor_cell_idx(jc,jb,cell_connect),jk,neighbor_cell_blk(jc,jb,cell_connect)))
! !             ENDIF
! !             !            ENDIF
! !           ENDDO
! !           
! !           ! Sum of all incoming antidiffusive fluxes into cell jc
! !           p_p = - (MIN(0._wp,z_mflx_anti(1))   &
! !             & + MIN(0._wp,z_mflx_anti(2))   &
! !             & + MIN(0._wp,z_mflx_anti(3)) )
! !           ! Sum of all outgoing antidiffusive fluxes out of cell jc
! !           p_m =  MAX(0._wp,z_mflx_anti(1))  &
! !             & + MAX(0._wp,z_mflx_anti(2))  &
! !             & + MAX(0._wp,z_mflx_anti(3))
! !           
! !           ! fraction which must multiply all fluxes out of cell jc to guarantee no
! !           ! undershoot
! !           ! Nominator: maximum allowable decrease of q
! !           r_m(jc,jk,jb) = (z_tracer_new_low(jc,jk,jb) - z_min ) / (p_m + dbl_eps)
! !           ! fraction which must multiply all fluxes into cell jc to guarantee no
! !           ! overshoot
! !           ! Nominator: maximum allowable increase of q
! !           r_p(jc,jk,jb) = (z_max - z_tracer_new_low(jc,jk,jb)) / (p_p + dbl_eps)
! !           
! !         ENDDO
! !       ENDDO
! !     ENDDO
! !     !ICON_OMP_END_DO NOWAIT
! !     !ICON_OMP_END_PARALLEL
! !     
! !     ! Synchronize r_m and r_p
! !     CALL sync_patch_array_mult(sync_c1, patch_2d, 2, r_m, r_p)
! !     
! !     ! 5. Now loop over all edges and determine the minimum fraction which must
! !     !    multiply the antidiffusive flux at the edge.
! !     !    At the end, compute new, limited fluxes which are then passed to the main
! !     !    program. Note that p_mflx_tracer_h now denotes the LIMITED flux.
! !     !ICON_OMP_PARALLEL
! !     !ICON_OMP_DO PRIVATE(start_index, end_index, je, jk, z_signum, r_frac) ICON_OMP_DEFAULT_SCHEDULE
! !     DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !       CALL get_index_range(edges_in_domain, jb, start_index, end_index)
! !       DO je = start_index, end_index
! !         DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,jb), end_level)
! !           ! IF( patch_3D%lsm_e(je,jk,jb) <= sea_boundary ) THEN
! !           z_signum = SIGN(1._wp, z_anti(je,jk,jb))
! !           
! !           ! The dolic_e is the min of the neigboring dolic_c, so
! !           ! zeroeing of r_m, r_p is not required
! !           
! !           !          IF (patch_3D%lsm_e(je,jk,jb) > sea_boundary) &
! !           !            CALL finish("","patch_3D%lsm_e(je,jk,jb) > sea_boundary")
! !           !          IF (patch_3D%lsm_c(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)) > sea_boundary) &
! !           !            CALL finish("","patch_3D%lsm_c(1) > sea_boundary")
! !           !          IF (patch_3D%lsm_c(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)) > sea_boundary) &
! !           !            CALL finish("","patch_3D%lsm_c(2) > sea_boundary")
! !           
! !           ! This does the same as an IF (z_signum > 0) THEN ... ELSE ... ENDIF,
! !           ! but is computationally more efficient
! !           r_frac = 0.5_wp * (                                                        &
! !             & (1._wp + z_signum) *                                           &
! !             & MIN(r_m(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)),  &
! !             & r_p(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)))  &
! !             & +  (1._wp - z_signum) *                                           &
! !             & MIN(r_m(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)),  &
! !             & r_p(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)))  )
! !           
! !           ! Limited flux
! !           p_mflx_tracer_h(je,jk,jb) = z_mflx_low(je,jk,jb)               &
! !             & + MIN(1._wp,r_frac) * z_anti(je,jk,jb)
! !           ! ELSE
! !           !   p_mflx_tracer_h(je,jk,jb)= 0.0_wp
! !           ! ENDIF
! !         END DO
! !       ENDDO
! !     ENDDO
! !     !ICON_OMP_END_DO NOWAIT
! !     !ICON_OMP_END_PARALLEL
! !     
! !   END SUBROUTINE hflx_limiter_oce_mo
! !   !-------------------------------------------------------------------------
END MODULE mo_oce_tracer_transport_horz
