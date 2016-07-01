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
#include "icon_definitions.inc"
!----------------------------
MODULE mo_ocean_tracer_transport_horz
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  USE mo_math_constants,            ONLY: dbl_eps
  USE mo_impl_constants,            ONLY: sea_boundary, SEA
  USE mo_ocean_nml,                 ONLY: n_zlev, l_edge_based, ab_gam,                   &
    & upwind, central,lax_friedrichs, horz_flux_twisted_vec_recon, miura_order1, flux_calculation_horz,      &
    & fct_high_order_flux,  fct_low_order_flux,FCT_Limiter_horz, fct_limiter_horz_zalesak,&
    & fct_limiter_horz_minmod, fct_limiter_horz_posdef, l_with_horz_tracer_diffusion, l_with_horz_tracer_advection,&
    &l_LAX_FRIEDRICHS, l_GRADIENT_RECONSTRUCTION, k_pot_temp_h, tracer_HorizontalAdvection_type, &
    & edge_based, cell_based
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma, p_test_run
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer
  USE mo_timer,                     ONLY: timer_start, timer_stop, timers_level, timer_adv_horz, timer_hflx_lim, &
    & timer_dif_horz, timer_extra10, timer_extra11, timer_extra12, timer_extra13, timer_extra15
  USE mo_ocean_types,               ONLY: t_hydro_ocean_state, t_ocean_tracer
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish !, message_text, message
  USE mo_ocean_physics
  USE mo_scalar_product,            ONLY: map_cell2edges_3d,map_edges2cell_3d, &
    & map_edges2edges_viacell_3d_const_z
  USE mo_ocean_math_operators,      ONLY: div_oce_3d, grad_fd_norm_oce_3d, grad_fd_norm_oce_3d_onBlock
  USE mo_ocean_diffusion,           ONLY: tracer_diffusion_horz
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff, no_primal_edges
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_c, sync_c1, sync_e, sync_patch_array, sync_patch_array_mult
  USE mo_mpi,                       ONLY: global_mpi_barrier
  USE mo_ocean_tracer_transport_vert, ONLY: upwind_vflux_oce
  USE mo_ocean_limiter,             ONLY: limiter_ocean_zalesak_horizontal,limiter_ocean_posdef_horizontal
  
  
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=12)           :: str_module    = 'oceTracHorz '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  !
  ! PUBLIC INTERFACE
  !
  PUBLIC :: advect_horz
  PUBLIC :: diffuse_horz  
  ! Private implemenation
  !
  PRIVATE :: advect_edge_based
  PRIVATE :: diffuse_edge_based  
  PRIVATE :: advect_cell_based
  PRIVATE :: diffuse_cell_based  
  PRIVATE :: fct_high_res
  PRIVATE :: flux_corr_transport_cell  
  PRIVATE :: flux_corr_transport_edge 
 
  PRIVATE :: upwind_hflux_oce  
  PRIVATE :: upwind_hflux_oce_mimetic
  PRIVATE :: central_hflux_oce  
  PRIVATE :: central_hflux_oce_mimetic
  PRIVATE :: miura_order1_hflux_oce  
  PRIVATE :: lax_friedrichs_hflux_oce 
  PRIVATE :: laxfr_upflux
 
  PRIVATE :: ratio_consecutive_gradients
  PRIVATE :: calculate_limiter_function
 
  !PRIVATE :: elad
  !PRIVATE :: hflx_limiter_oce_mo
  !PRIVATE :: hflx_limiter_oce_zalesak
  
  
  INTEGER, PARAMETER :: top=1
  
CONTAINS
  !-----------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE advect_horz( patch_3d,&
    & trac_old,              &
    & p_os,                  &
    & operators_coefficients,&
    & k_h,                   &
    & h_old,                 &
    & h_new,                 &
    & div_flux_horz,         &
    & div_flux_vert,         &
    & tracer_index,          &
    & horizontally_diffused_tracer)

    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp)                               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET      :: p_os
    TYPE(t_operator_coeff), INTENT(in)     :: operators_coefficients
    REAL(wp), INTENT(in)                   :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                   :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                   :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)                :: div_flux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)                :: div_flux_vert(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)       
    INTEGER,  INTENT(in)                   :: tracer_index
    REAL(wp), INTENT(inout), OPTIONAL      :: horizontally_diffused_tracer(:,:,:)    
    !
    !-------------------------------------------------------------------------------    
    start_timer(timer_adv_horz,2)

    SELECT CASE(tracer_HorizontalAdvection_type)
    CASE(edge_based)
      CALL advect_edge_based( patch_3d, &
      & trac_old,              &
      & p_os,                  &
      & operators_coefficients,&
      & k_h,                   &
      & h_old,                 &
      & h_new,                 &
      & div_flux_horz,         &
      & div_flux_vert,         &
      & tracer_index)
    
    CASE(cell_based)

      CALL advect_cell_based( patch_3d, &
      & trac_old,              &
      & p_os,                  &
      & operators_coefficients,&
      & k_h,                   &
      & h_old,                 &
      & h_new,                 &
      & div_flux_horz,         &
      & div_flux_vert,         &
      & tracer_index) 
    CASE default
      CALL finish("advect_horz","uknown tracer_HorizontalAdvection_type")
    END SELECT

    stop_timer(timer_adv_horz,2)
     
  END SUBROUTINE advect_horz
  !-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE diffuse_horz( patch_3d,          &
    & trac_old,            &
    & p_os,                &
    & operators_coefficients,          &
    & k_h,                 &
    & h_old,               &
    & h_new,               &
    & flux_horz,           &
    & horizontally_diffused_tracer)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp)                               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET      :: p_os
    TYPE(t_operator_coeff), INTENT(in)     :: operators_coefficients
    REAL(wp), INTENT(in)                   :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                   :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                   :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)                :: flux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout), OPTIONAL      :: horizontally_diffused_tracer(:,:,:)
    !
    !-------------------------------------------------------------------------------    
    start_timer(timer_dif_horz,3)
    
    IF(l_edge_based)THEN
      CALL diffuse_edge_based( patch_3d, &
      & trac_old,            &
      & p_os,                &
      & operators_coefficients,          &
      & k_h,                 &
      & h_old,               &
      & h_new,               &
      & flux_horz)
    
    ELSE
      CALL diffuse_cell_based( patch_3d,          &
      & trac_old,            &
      & p_os,                &
      & operators_coefficients,          &
      & k_h,                 &
      & h_old,               &
      & h_new,               &
      & flux_horz) !,           &
      ! & horizontally_diffused_tracer)
    
    ENDIF

    stop_timer(timer_dif_horz,3)
     
  END SUBROUTINE diffuse_horz
  !-------------------------------------------------------------------------------


  !-----------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE advect_cell_based( patch_3d,          &
    & trac_old,            &
    & p_os,                &
    & operators_coefficients,          &
    & k_h,                 &
    & h_old,               &
    & h_new,               &
    & div_advflux_horz,    &
    & div_advflux_vert,    &
    &tracer_index)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp)                             :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET    :: p_os
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                 :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                 :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)              :: div_advflux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)              :: div_advflux_vert(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)   
    INTEGER, INTENT(in)                  :: tracer_index
    !
    !Local variables
    INTEGER :: start_index, end_index
    INTEGER :: jc, level, blockNo
    REAL(wp) :: z_adv_flux_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: p_vn_c(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-------------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------------
    !Calculate tracer fluxes at edges
    !This step takes already the edge length into account
    !but not the edge height

    IF ( l_with_horz_tracer_advection ) THEN
      
      ! Initialize timer for horizontal advection
      
      SELECT CASE(flux_calculation_horz)
      
      CASE(upwind)

        CALL upwind_hflux_oce( patch_3d,  &
        & trac_old,                       &
        & p_os%p_diag%mass_flx_e,         & 
        & z_adv_flux_h)

        
      CASE(central)  
      
        CALL central_hflux_oce( patch_3d, &
        & trac_old,                       &
        & p_os%p_diag%mass_flx_e,         &
        & z_adv_flux_h)    
         
        
      CASE(miura_order1)
        
        !MIURA-type upwind-biased estimate of tracer flux
        z_adv_flux_h = 0.0_wp
        CALL miura_order1_hflux_oce( patch_3d,   &
          & trac_old,                            &
          & p_os%p_diag%mass_flx_e,              &
          & operators_coefficients,              &
          & z_adv_flux_h )
        
      CASE(horz_flux_twisted_vec_recon)
        ! inUse

        CALL flux_corr_transport_cell( patch_3d,  &
          & trac_old,                             &
          & p_os,                                 &
          & operators_coefficients,               &
          & k_h,                                  &
          & h_old,                                &
          & h_new,                                &
          & z_adv_flux_h,                         &
          & div_advflux_vert,                     &
          & tracer_index)
              
      CASE default
        CALL finish('TRIM(advect_diffuse_flux_horz)',"This flux option is not supported")
        
      END SELECT
      !---------------------------------------------------------------------
      
      !Calculate divergence of advective fluxes
      CALL div_oce_3d( z_adv_flux_h, patch_3D, operators_coefficients%div_coeff, div_advflux_horz, subset_range=cells_in_domain )

      
    ELSE ! no l_with_horz_tracer_advection
      div_advflux_horz (:,:,:) = 0.0_wp
    ENDIF ! l_with_horz_tracer_advection


    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('TracAdvHorz: adv_flux_h'     ,z_adv_flux_h ,str_module,idt_src,patch_2d%edges%owned)
    CALL dbg_print('TracAdvHorz: div adv_flux_h' ,div_advflux_horz,str_module,idt_src,cells_in_domain)
    !---------------------------------------------------------------------
    
  END SUBROUTINE advect_cell_based
  !-------------------------------------------------------------------------------


  !-----------------------------------------------------------------------
  SUBROUTINE diffuse_cell_based( patch_3d,          &
    & trac_old,            &
    & p_os,                &
    & operators_coefficients,          &
    & k_h,                 &
    & h_old,               &
    & h_new,               &
    & div_diff_flux_horz)!,           &
    ! & horizontally_diffused_tracer)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp)                             :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET    :: p_os
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                 :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                 :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)              :: div_diff_flux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    ! REAL(wp), INTENT(inout), OPTIONAL    :: horizontally_diffused_tracer(:,:,:)
    !
    !
    !Local variables
    INTEGER :: start_index, end_index
    INTEGER :: jc, level, blockNo
    REAL(wp) :: z_diff_flux_h(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: p_vn_c(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: diff_option=1
    INTEGER :: diff_option_standard=1
    !-------------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------------
    
    !The diffusion part: calculate horizontal diffusive flux

    IF ( l_with_horz_tracer_diffusion) THEN
      
      IF(diff_option==diff_option_standard)THEN
        CALL tracer_diffusion_horz( patch_3d,     &
        & trac_old,     &
        & p_os,         &
        & z_diff_flux_h,&
        & k_h,          &        
        & subset_range = edges_in_domain)
        
      ELSE!this is only for testing
        CALL tracer_diffusion_horz_ptp( patch_3d,     &
        & trac_old,     &
        & p_os, operators_coefficients,        &
        & k_h,z_diff_flux_h)
      
      ENDIF
      
    ELSE
      div_diff_flux_horz(:,:,:) = 0.0_wp
    ENDIF
    
    !Calculate divergence of diffusive fluxes
    CALL div_oce_3d( z_diff_flux_h, patch_3D, operators_coefficients%div_coeff, div_diff_flux_horz )

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('TracDiffHorz: diff_flux_h',       z_diff_flux_h ,str_module,idt_src,cells_in_domain)
    CALL dbg_print('TracDiffHorz: div_diff_flux_horz' ,div_diff_flux_horz    ,str_module,idt_src,cells_in_domain)
    !---------------------------------------------------------------------
    
  END SUBROUTINE diffuse_cell_based
  !-------------------------------------------------------------------------------



  !-----------------------------------------------------------------------
  SUBROUTINE tracer_diffusion_horz_ptp( patch_3d,          &
    & trac_old,            &
    & p_os,                &
    & operators_coefficients,          &
    & k_h,                 &
    & diff_flux_h)!,           &
    ! & horizontally_diffused_tracer)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp)                             :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET    :: p_os
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(inout)              :: diff_flux_h(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%nblks_e)
    ! REAL(wp), INTENT(inout), OPTIONAL    :: horizontally_diffused_tracer(:,:,:)
    !
    !
    !Local variables
    INTEGER :: start_index, end_index, start_edge_index, end_edge_index
    INTEGER :: il_1, ib_c1,il_2, ib_c2
    INTEGER :: jc, jk, blockNo
    REAL(wp) :: z_div_diff_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: z_diff_flux_h(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: p_vn_c(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    REAL(wp) :: grad_T_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: flux_vec_horz_center(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: grad_T_vec_horz(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    !-------------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain 
    grad_T_horz(1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e)=0.0_wp
    
     grad_T_vec_horz(1:nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)%x(1)=0.0_wp
     grad_T_vec_horz(1:nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)%x(2)=0.0_wp
     grad_T_vec_horz(1:nproma,1:n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)%x(3)=0.0_wp

    !-------------------------------------------------------------------------------
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      
      !1a) calculate horizontal gradient of tracer
      CALL grad_fd_norm_oce_3d_onBlock ( &
        & trac_old, &
        & patch_3D,                           &
        & operators_coefficients%grad_coeff(:,:,blockNo), &
        & grad_T_horz(:,:,blockNo),           &
        & start_edge_index, end_edge_index, blockNo)

    END DO ! blocks
   CALL sync_patch_array(sync_e, patch_2D, grad_T_horz)

   CALL map_edges2cell_3d(patch_3D, &
        & grad_T_horz,               &
        & operators_coefficients,                &
        & grad_T_vec_horz)

   DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
   
     flux_vec_horz_center(1:nproma,1:n_zlev,blockNo)%x(1)=0.0_wp
     flux_vec_horz_center(1:nproma,1:n_zlev,blockNo)%x(2)=0.0_wp
     flux_vec_horz_center(1:nproma,1:n_zlev,blockNo)%x(3)=0.0_wp  
     
     CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      
     DO jc = start_index, end_index

        
       DO jk = 1, patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)
         flux_vec_horz_center(jc,jk,blockNo)%x=k_pot_temp_h*grad_T_vec_horz(jc,jk,blockNo)%x
       END DO
     END DO
   END DO
 
   CALL map_cell2edges_3D( patch_3D,flux_vec_horz_center, diff_flux_h,operators_coefficients)
  END SUBROUTINE tracer_diffusion_horz_ptp
  !-------------------------------------------------------------------------------



  !-----------------------------------------------------------------------
  SUBROUTINE advect_edge_based( patch_3d,&
    & trac_old,               &
    & p_os,                   &
    & operators_coefficients, &
    & k_h,                    &
    & h_old,                  &
    & h_new,                  &
    & div_advflux_horz,       &
    & div_advflux_vert,       &
    & tracer_index)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp)                             :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET    :: p_os
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                 :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                 :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)              :: div_advflux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)              :: div_advflux_vert(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    INTEGER, INTENT(in) :: tracer_index
    !
    !
    !Local variables
    INTEGER :: start_index_e, end_index_e
    INTEGER :: je, level, blockNo
    REAL(wp) :: z_adv_flux_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_div_adv_h  (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    
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
    z_div_adv_h  (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks)=0.0_wp

    !Calculate tracer fluxes at edges
    !This step takes already the edge length into account
    !but not the edge height
    IF ( l_with_horz_tracer_advection ) THEN
      
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
          & operators_coefficients,            &
          & z_adv_flux_h  )
         
                   
      CASE(horz_flux_twisted_vec_recon)
      
        CALL flux_corr_transport_edge( patch_3d,&
          & trac_old,                           &
          & p_os,                               &
          & operators_coefficients,             &
          & k_h,                                &
          & h_old,                              &
          & h_new,                              &
          & z_adv_flux_h,                       &
          & div_advflux_vert,                   &
          & tracer_index)      
              
      CASE default
        CALL finish('TRIM(advect_diffuse_flux_horz)',"This flux option is not supported")
        
      END SELECT     
      CALL sync_patch_array(sync_e, patch_2d, z_adv_flux_h)      
      
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=5  ! output print level (1-5, fix)
      CALL dbg_print('aft. AdvHorz: adv_flux_h',z_adv_flux_h,str_module,idt_src,patch_2d%edges%owned)
      !---------------------------------------------------------------------               
      
      !Calculate divergence of advective fluxes
      CALL div_oce_3d( z_adv_flux_h, patch_3D, operators_coefficients%div_coeff, div_advflux_horz )

      
    ELSE
      div_advflux_horz(:,:,:) = 0.0_wp
    ENDIF ! l_with_horz_tracer_diffusion
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('AdvDifHorz: adv_flux_h'     ,z_adv_flux_h ,str_module,idt_src,patch_2d%edges%owned)
    CALL dbg_print('AdvDifHorz: div adv_flux_h' ,div_advflux_horz  ,str_module,idt_src,cells_in_domain)
    !---------------------------------------------------------------------   
  END SUBROUTINE advect_edge_based
  !-------------------------------------------------------------------------------


  !-----------------------------------------------------------------------
  SUBROUTINE diffuse_edge_based( patch_3d,          &
    & trac_old,            &
    & p_os,                &
    & operators_coefficients,          &
    & k_h,                 &
    & h_old,               &
    & h_new,               &
    & div_diffflux_horz)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp)                             :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET    :: p_os
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                 :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                 :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)              :: div_diffflux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !
    !
    !Local variables
    INTEGER :: start_index, end_index
    INTEGER :: jc, level, blockNo
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
    z_div_diff_h (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks)=0.0_wp
    z_diff_flux_h(1:nproma,1:n_zlev,1:patch_2d%nblks_e)=0.0_wp

    
    !The diffusion part: calculate horizontal diffusive flux
    IF ( l_with_horz_tracer_diffusion)THEN 

      CALL tracer_diffusion_horz( patch_3d,     &
        & trac_old,     &
        & p_os,         &
        & z_diff_flux_h,&
        & k_h,          &
        & subset_range = edges_in_domain)
      
      !Calculate divergence of diffusive fluxes
      CALL div_oce_3d( z_diff_flux_h, patch_3D, operators_coefficients%div_coeff, div_diffflux_horz)
      
    ELSE
      div_diffflux_horz(:,:,:) = 0.0_wp
    ENDIF
    
    
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)    
    CALL dbg_print('AdvDifHorz: div diff_flux_h',div_diffflux_horz ,str_module,idt_src,cells_in_domain)
    !---------------------------------------------------------------------   
  END SUBROUTINE diffuse_edge_based
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE flux_corr_transport_edge( patch_3d,&
    & trac_old,                                 &
    & p_os,                                     &
    & operators_coefficients,                   &
    & k_h,                                      &
    & h_old,                                    &
    & h_new,                                    &
    & adv_flux_h,                               &
    & div_advflux_vert,                         &
    & tracer_index)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp)                               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET      :: p_os
    TYPE(t_operator_coeff), INTENT(in)     :: operators_coefficients
    REAL(wp), INTENT(in)                   :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                   :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                   :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), TARGET, INTENT(inout)        :: adv_flux_h(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)!< variable in which the upwind flux is stored
    REAL(wp), INTENT(inout)                :: div_advflux_vert(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)     
    INTEGER, INTENT(in)                    :: tracer_index
    !Local Variables
    INTEGER          :: je, blockNo, level, start_index_e, end_index_e
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
        & operators_coefficients,                        &
        & z_adv_flux_high )
    
    CASE DEFAULT
      CALL finish('TRIM(flux_corr_transport_edge)',"This high-order  option is not supported")
    END SELECT
    
 
    !2)call limiter 
    SELECT CASE(fct_limiter_horz)
    
    CASE(fct_limiter_horz_minmod)
    
      CALL fct_high_res( patch_3d, &
      & trac_old,                              &
      & p_os%p_diag%mass_flx_e,                &
      & z_adv_flux_low,                        &    
      & z_adv_flux_high,                       &      
      & operators_coefficients,                            &
      & k_h,                                   &
      & adv_flux_h )

      !CALL hflx_limiter_oce_posdef( patch_3d,&
      !& trac_old,                    &
      !& z_adv_flux_high,              &    
      !& operators_coefficients )
      
    CASE(fct_limiter_horz_zalesak)

       CALL limiter_ocean_zalesak_horizontal( patch_3d, &
       & p_os%p_diag%w_time_weighted,           &
       & trac_old,                              &
       & p_os%p_diag%mass_flx_e,                &
       & z_adv_flux_low,                        &
       & z_adv_flux_high,                       &    
       & adv_flux_h,                            &
       & div_advflux_vert,                      &     
       & operators_coefficients,                &
       & h_old,                                 &
       & h_new,                                 &
       & p_os%p_diag%zlim,tracer_index             )       
       
    CASE(fct_limiter_horz_posdef)  
      CALL limiter_ocean_posdef_horizontal( patch_3d,  &
      & trac_old,                              &
      & z_adv_flux_high,                       &    
      & adv_flux_h,                            &       
      & operators_coefficients )

    CASE DEFAULT
      CALL finish('TRIM(flux_corr_transport_edge)',"This limiter_type option is not supported")
    END SELECT
    
  END SUBROUTINE flux_corr_transport_edge
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE flux_corr_transport_cell( patch_3d, &
    & trac_old,                                  & 
    & p_os,                                      &
    & operators_coefficients,                    &
    & k_h,                                       &
    & h_old,                                     &
    & h_new,                                     &
    & adv_flux_h,                                &
    & div_advflux_vert,                          &
    & tracer_index)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp)                               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET      :: p_os
    TYPE(t_operator_coeff), INTENT(in)     :: operators_coefficients
    REAL(wp), INTENT(in)                   :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                   :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                   :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET, INTENT(inout)         :: adv_flux_h(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)!< variable in which the upwind flux is stored
    REAL(wp), INTENT(inout)                :: div_advflux_vert(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)        
    INTEGER, INTENT(in) :: tracer_index
    !Local Variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices    
    INTEGER  :: je,blockNo,level,start_index_e, end_index_e, edge_index, jc
    INTEGER  :: start_index, end_index
    REAL(wp) :: z_adv_flux_high(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_adv_flux_gradrecon(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_adv_flux_low (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_cartesian_coordinates) :: p_vn_c(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)  
    REAL(wp) :: grad_C_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)       
    !-------------------------------------------------------------------------------
    
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    
!     z_adv_flux_high(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%nblks_e)=0.0_wp
!     z_adv_flux_low(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%nblks_e)=0.0_wp
!     z_adv_flux_gradrecon(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%nblks_e)=0.0_wp    
!     p_vn_c(1:nproma, 1:n_zlev, 1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)%x(1)=0.0_wp        
!     p_vn_c(1:nproma, 1:n_zlev, 1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)%x(2)=0.0_wp        
!     p_vn_c(1:nproma, 1:n_zlev, 1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)%x(3)=0.0_wp   
    !-------------------------------------------------------------------------------
    !1) provide high- & low order tracer flux
    z_adv_flux_high = 0.0_wp
    z_adv_flux_low  = 0.0_wp 
    SELECT CASE(fct_low_order_flux)
	
    CASE(upwind)
	
      !upwind serves as low order flux
      start_detail_timer(timer_extra11,5)
      CALL upwind_hflux_oce( patch_3d,  &
        & trac_old,                       &
        & p_os%p_diag%mass_flx_e,         &
        & z_adv_flux_low )
      stop_detail_timer(timer_extra11,5)
        
    CASE(miura_order1)
      CALL miura_order1_hflux_oce( patch_3d,   &
        & trac_old,                            &
        & p_os%p_diag%mass_flx_e,              &
        & operators_coefficients,                          &
        & z_adv_flux_low )      
    
    CASE DEFAULT
      CALL finish('TRIM(flux_corr_transport_edge)',"This low-order  option is not supported")   
    END SELECT    



    SELECT CASE(fct_high_order_flux)
    CASE(central)
    
      !central as high order flux
      start_detail_timer(timer_extra12,5)
      CALL central_hflux_oce( patch_3d,   &
        & trac_old,                       &
        & p_os%p_diag%mass_flx_e,         &
        & z_adv_flux_high)    
      stop_detail_timer(timer_extra12,5)
 
    CASE(horz_flux_twisted_vec_recon)
      !mimetic fluc calculation high order flux
      CALL map_edges2edges_viacell_3d_const_z( patch_3d, &
        & p_os%p_diag%vn_time_weighted,          &
        & operators_coefficients,                &
        & z_adv_flux_high,                       &
        & trac_old)
        
      IF(l_GRADIENT_RECONSTRUCTION)THEN
    
        grad_C_horz(1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e)= 0.0_wp
        
        CALL grad_fd_norm_oce_3d ( &
          & trac_old, &
          & patch_3D,              &
          & operators_coefficients%grad_coeff, &
          & grad_C_horz)


        CALL map_edges2edges_viacell_3d_const_z( patch_3d, &
          & grad_C_horz,                           &
          & operators_coefficients,                            &
          & z_adv_flux_gradrecon)  
    
        ! loop through all patch edges (and blocks)
        !!ICON_OMP_PARALLEL
        !!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDUL
        DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, blockNo, start_index_e, end_index_e)
          DO edge_index = start_index_e, end_index_e
            DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo)
          
              z_adv_flux_high(edge_index,level,blockNo) &
              &= z_adv_flux_high(edge_index,level,blockNo) + z_adv_flux_gradrecon(edge_index,level,blockNo)
            END DO  ! end loop over edges
          END DO  ! end loop over levels
        END DO  ! end loop over blocks
        !!ICON_OMP_END_DO NOWAIT
        !!ICON_OMP_END_PARALLEL
      ENDIF

      IF(l_LAX_FRIEDRICHS)THEN
        iilc => patch_2d%edges%cell_idx
        iibc => patch_2d%edges%cell_blk

        ! loop through all patch edges (and blocks)
        !!ICON_OMP_PARALLEL
        !!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDUL
        DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, blockNo, start_index_e, end_index_e)
          DO edge_index = start_index_e, end_index_e
            DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo)

            z_adv_flux_high(edge_index,level,blockNo) =  z_adv_flux_high(edge_index,level,blockNo)    &
              & + 0.5_wp * p_os%p_diag%vn_time_weighted(edge_index,level,blockNo)                     &
              &          * p_os%p_diag%vn_time_weighted(edge_index,level,blockNo) * dtime             &
              &          * patch_2d%edges%inv_dual_edge_length(edge_index,blockNo)                    &
              &          * ( trac_old(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2))    &
              &              -trac_old(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1)))
              
            END DO  ! end loop over edges
          END DO  ! end loop over levels
        END DO  ! end loop over blocks
        !!ICON_OMP_END_DO NOWAIT
        !!ICON_OMP_END_PARALLEL
      ENDIF!l_LAX_FRIEDRICHS
    END SELECT
    !-----------------------------------------------------------------------
    
    !2)call limiter
    SELECT CASE(fct_limiter_horz)
    CASE(fct_Limiter_horz_minmod)
    
      CALL fct_high_res( patch_3d, &
        & trac_old,                              &
        & p_os%p_diag%mass_flx_e,                &
        & z_adv_flux_low,                        &
        & z_adv_flux_high,                       &
        & operators_coefficients,                            &
        & k_h,                                   &
        & adv_flux_h )
      
    CASE(fct_limiter_horz_zalesak)
	
      ! inUse
     ! adv_flux_h=z_adv_flux_high
      start_detail_timer(timer_extra13,4)
      CALL limiter_ocean_zalesak_horizontal( patch_3d,   &
        & p_os%p_diag%w_time_weighted,           &
        & trac_old,                              &
        & p_os%p_diag%mass_flx_e,                &
        & z_adv_flux_low,                        &
        & z_adv_flux_high,                       &
        & adv_flux_h,                            &
        & div_advflux_vert,                      &            
        & operators_coefficients,                &
        & h_old,                                 &
        & h_new,                                 &
        & p_os%p_diag%zlim,                      &
        & tracer_index   )       
      stop_detail_timer(timer_extra13,4)
      
    CASE(fct_limiter_horz_posdef)  
      CALL limiter_ocean_posdef_horizontal( patch_3d,    &
        & trac_old,                              &
        & z_adv_flux_high,                       &
        & adv_flux_h,                            &        
        & operators_coefficients )
  
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
    & operators_coefficients, k_h, adv_flux_h)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)           :: trac_old(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< advected cell centered variable
    REAL(wp), INTENT(inout)           :: vn_e (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) !< normal velocity on edges
    REAL(wp), INTENT(IN)              :: adv_flux_low (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(IN)              :: adv_flux_high(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff)            :: operators_coefficients
    REAL(wp), INTENT(in)              :: k_h(:,:,:)         !horizontal mixing coeff
    
    REAL(wp), INTENT(OUT) :: adv_flux_h(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)!< variable in which the upwind flux is stored

    REAL(wp) :: z_diff_flux_h(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)    
    REAL(wp) :: z_consec_grad(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_mflux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_limit_phi(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_limit_psi(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_limit_sigma(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    INTEGER :: start_index_e, end_index_e
    INTEGER :: level, blockNo, edge_index
    INTEGER, DIMENSION(:,:,:), POINTER :: cellOfEdge_idx, cellOfEdge_blk
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    TYPE(t_patch), POINTER        :: patch_2d
    !-------------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !------------------------------------------------------------------------------   
      
    ! Set pointers to index-arrays
    ! line and block indices of two neighboring cells
    cellOfEdge_idx  => patch_2d%edges%cell_idx
    cellOfEdge_blk  => patch_2d%edges%cell_blk
    
    z_consec_grad(1:nproma,1:n_zlev,1:patch_2d%nblks_e) = 0.0_wp
    z_mflux(1:nproma,1:n_zlev,1:patch_2d%nblks_e)       = 0.0_wp
    z_limit_phi(1:nproma,1:n_zlev,1:patch_2d%nblks_e)   = 0.0_wp
    z_limit_psi(1:nproma,1:n_zlev,1:patch_2d%nblks_e)   = 0.0_wp
    z_limit_sigma(1:nproma,1:n_zlev,1:patch_2d%nblks_e) = 0.0_wp   
    
    !1) calculate phi-part of limiter
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index_e, end_index_e)
      DO level = 1, n_zlev
        DO edge_index = start_index_e, end_index_e
        
           !calculate the difference between adjacent tracer values, without
           !dividing by dual edge length.
           z_diff_flux_h(edge_index,level,blockNo) = &
             & trac_old(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2))   &
             & - trac_old(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1))
           
           !outflow case
           !if vn_e(edge_index,level,blockNo) is positive then the flow is from cell 1 towards cell 2
           z_mflux(edge_index,level,blockNo) = vn_e(edge_index,level,blockNo)

          !This corresponds to eq. (14) in Casulli-Zanolli.
          !Their "little" phi is here called sigma. We have
          !canceled down their primal edge-length and thickness 
          !in numerator and denominator of sigma.
           IF (z_mflux(edge_index,level,blockNo) == 0.0_wp) THEN
            z_limit_sigma(edge_index,level,blockNo) = 1.0_wp   !  z_mflux can be zero!
          ELSE
            z_limit_sigma(edge_index,level,blockNo) = MIN(1.0_wp, &
              & 2.0_wp * k_h(edge_index,level,blockNo) / &
              & (patch_2d%edges%dual_edge_length(edge_index,blockNo) * &
              & z_mflux(edge_index,level,blockNo)))
          ENDIF

        END DO
      END DO
    END DO
    CALL sync_patch_array(sync_e, patch_2d, z_diff_flux_h)    
    CALL sync_patch_array(sync_e, patch_2d, z_mflux)
    CALL sync_patch_array(sync_e, patch_2d, z_limit_sigma)        

    !This corresponds to (16) in Casulli-Zanolli.
     CALL ratio_consecutive_gradients(patch_3d,operators_coefficients,vn_e, trac_old,z_consec_grad)

    !3) calculate psi-part of limiter (see (17)-(19) in Casulli-Zanolli).
    CALL calculate_limiter_function(patch_3d, z_limit_sigma, z_consec_grad,z_limit_phi)

    !4) Calculate limited advective flux
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index_e, end_index_e)
      DO level = 1, n_zlev
        DO edge_index = start_index_e, end_index_e

          IF(z_consec_grad(edge_index,level,blockNo)>0)THEN
            z_limit_psi(edge_index,level,blockNo) = &
              & z_limit_phi(edge_index,level,blockNo)-z_limit_sigma(edge_index,level,blockNo)
          ELSE
            z_limit_psi(edge_index,level,blockNo) = 0.0_wp
          ENDIF
          !only low gives upwind, low +mflux*diff and without limiter gives central
          adv_flux_h(edge_index,level,blockNo) =             &
            & adv_flux_low (edge_index,level,blockNo)        &
            & + 0.5_wp * z_mflux(edge_index,level,blockNo)   &
            &   * z_diff_flux_h(edge_index,level,blockNo)    &
            &   * z_limit_psi(edge_index,level,blockNo)
          
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
   SUBROUTINE ratio_consecutive_gradients( patch_3d, operators_coefficients,      &
    & vn_time_weighted, &
    & trac_old,         &
    & consec_grad)
    
    TYPE(t_patch_3d),TARGET, INTENT(in):: patch_3d
    TYPE(t_operator_coeff)             :: operators_coefficients
    REAL(wp), INTENT(in)               :: vn_time_weighted(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)            :: consec_grad(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%nblks_e)
    !
    !Local variables
    !REAL(wp) :: delta_z
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: jc, level, blockNo!, edge_index!, ic,ib
    INTEGER :: i_edge, ii_e, ib_e
    REAL(wp) :: z_diff_trac    
    REAL(wp) :: z_cellsum_mass_flx_in(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_cellsum_tracdiff_flx_in(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_mflux(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks,1:no_primal_edges)
    INTEGER, DIMENSION(:,:,:), POINTER :: cellOfEdge_idx, cellOfEdge_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    cells_in_domain => patch_2D%cells%in_domain
    !-------------------------------------------------------------------------------
    consec_grad              (1:nproma,1:n_zlev,1:patch_2D%nblks_e) = 0.0_wp
    z_cellsum_mass_flx_in    (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
    z_cellsum_tracdiff_flx_in(1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
    z_mflux(1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks,1:no_primal_edges) = 0.0_wp
      
    ! Set pointers to index-arrays
    ! line and block indices of two neighboring cells
    cellOfEdge_idx  => patch_2D%edges%cell_idx
    cellOfEdge_blk  => patch_2D%edges%cell_blk
    edge_of_cell_idx  => patch_2D%cells%edge_idx
    edge_of_cell_blk  => patch_2D%cells%edge_blk
    
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, i_startidx_c, i_endidx_c)
      DO level = 1, n_zlev
        DO jc = i_startidx_c, i_endidx_c
          
          DO i_edge=1,no_primal_edges
            ii_e = patch_2D%cells%edge_idx(jc,blockNo,i_edge)
            ib_e = patch_2D%cells%edge_blk(jc,blockNo,i_edge)
            
            !calculate mass flux of actual edge.
            z_mflux(jc,level,blockNo, i_edge) = &
              & operators_coefficients%div_coeff(jc,level,blockNo,i_edge)*patch_2D%cells%area(jc,blockNo)                     &
              & * vn_time_weighted(edge_of_cell_idx(jc,blockNo,i_edge),level,edge_of_cell_blk(jc,blockNo,i_edge)) &
              & * patch_3d%p_patch_1d(1)%prism_thick_e(ii_e,level,ib_e)
            
            !consider here inflow edges only
            IF( z_mflux(jc,level,blockNo, i_edge)<0.0)THEN
            
            !IF( z_mflux(jc,level,blockNo, i_edge)>0.0)THEN
              !sum up mass fluxes over edge
              z_cellsum_mass_flx_in(jc,level,blockNo) = z_cellsum_mass_flx_in(jc,level,blockNo)&
              &+z_mflux(jc,level,blockNo, i_edge)
                  
              !one of the two terms below is zero.
              z_diff_trac = &
              & (trac_old(cellOfEdge_idx(ii_e,ib_e,2),level, cellOfEdge_blk(ii_e,ib_e,2))-trac_old(jc,level,blockNo))&
              &+(trac_old(cellOfEdge_idx(ii_e,ib_e,1),level, cellOfEdge_blk(ii_e,ib_e,1))-trac_old(jc,level,blockNo) )
 
              !This is in Casulli Zanolli, eq. (16) the numerator of
              !the consecutive gradient
              z_cellsum_tracdiff_flx_in(jc,level,blockNo) &
              &= z_cellsum_tracdiff_flx_in(jc,level,blockNo)+z_mflux(jc,level,blockNo, i_edge)*z_diff_trac!&
              !&*patch%cells%edge_orientation(jc,blockNo,i_edge)  !diff_flux_h(iii_e,level,iib_e)!
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
    
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, i_startidx_c, i_endidx_c)
      DO level = 1, n_zlev
        DO jc = i_startidx_c, i_endidx_c
          
          DO i_edge=1,no_primal_edges
            
            !For positive velocities the flow is directed from cell1 to cell2
             IF(z_mflux(jc,level,blockNo, i_edge)>=0.0)THEN
             !IF(z_mflux(jc,level,blockNo, i_edge)<0.0)THEN

                ii_e = patch_2D%cells%edge_idx(jc,blockNo,i_edge)
                ib_e = patch_2D%cells%edge_blk(jc,blockNo,i_edge)
                
                z_diff_trac = &
                  &  ( trac_old(cellOfEdge_idx(ii_e,ib_e,2),level,cellOfEdge_blk(ii_e,ib_e,2))  &
                  &     -trac_old(jc,level,blockNo))                                                &
                  &+ ( trac_old(cellOfEdge_idx(ii_e,ib_e,1),level, cellOfEdge_blk(ii_e,ib_e,1)) &
                  &    - trac_old(jc,level,blockNo))
            
              
              IF(z_cellsum_mass_flx_in(jc,level,blockNo)/=0.0_wp.AND.z_diff_trac/=0.0_wp)THEN
                
                consec_grad(edge_of_cell_idx(jc,blockNo,i_edge),level,edge_of_cell_blk(jc,blockNo,i_edge))&
                & = z_cellsum_tracdiff_flx_in(jc,level,blockNo)/z_cellsum_mass_flx_in(jc,level,blockNo)
                                
                consec_grad(edge_of_cell_idx(jc,blockNo,i_edge),level,edge_of_cell_blk(jc,blockNo,i_edge))&
                & = (consec_grad(ii_e,level,ib_e)/z_diff_trac)!*patch%cells%edge_orientation(jc,blockNo,i_edge)
              ENDIF
            END IF
          END DO
        END DO
      END DO
    END DO    
    CALL sync_patch_array(sync_e, patch_2D, consec_grad)
    
  END SUBROUTINE ratio_consecutive_gradients
  !-------------------------------------------------------------------------------

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
    INTEGER :: level, blockNo, edge_index
    
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
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, i_startidx_e, i_endidx_e)
        DO level = 1, n_zlev
          DO edge_index = i_startidx_e, i_endidx_e
            
            limit_phi(edge_index,level,blockNo)&
              & =MAX(limit_sigma(edge_index,level,blockNo), MIN(1.0_wp,consec_grad(edge_index,level,blockNo)))
          END DO
        END DO
      END DO
      
!     CASE(superbee)
!       DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!         CALL get_index_range(edges_in_domain, blockNo, i_startidx_e, i_endidx_e)
!         DO level = 1, n_zlev
!           DO edge_index = i_startidx_e, i_endidx_e
!             limit_phi(edge_index,level,blockNo) = MAX(limit_sigma(edge_index,level,blockNo),&
!               & MIN(1.0_wp,2.0_wp*consec_grad(edge_index,level,blockNo)),&
!               & MIN(2.0_wp,consec_grad(edge_index,level,blockNo)))
!           END DO
!         END DO
!       END DO
!     CASE(van_leer)
!       DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!         CALL get_index_range(edges_in_domain, blockNo, i_startidx_e, i_endidx_e)
!         DO level = 1, n_zlev
!           DO edge_index = i_startidx_e, i_endidx_e
!             limit_phi(edge_index,level,blockNo) = MAX(limit_sigma(edge_index,level,blockNo),&
!               & (consec_grad(edge_index,level,blockNo)+ABS(consec_grad(edge_index,level,blockNo))&
!               & /(1.0_wp+ABS(consec_grad(edge_index,level,blockNo)))))
!               IF(limit_phi(edge_index,level,blockNo)>2.0)THEN
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
  SUBROUTINE upwind_hflux_oce_mimetic( patch_3d, flux_cc,operators_coefficients,&
    & edge_upwind_flux, opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    TYPE(t_cartesian_coordinates)      :: flux_cc(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_operator_coeff), INTENT(in) :: operators_coefficients
    REAL(wp), INTENT(inout)            :: edge_upwind_flux(nproma,n_zlev, patch_3d%p_patch_2d(1)%nblks_e)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level    ! optional vertical end level
    ! local variables
    !INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc
    INTEGER :: start_level, end_level
    INTEGER :: start_index, end_index
    INTEGER :: edge_index, level, blockNo        !< index of edge, vert level, block
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
    !!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      edge_upwind_flux(:,:,blockNo) = 0.0_wp
      DO edge_index = start_index, end_index
        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          edge_upwind_flux(edge_index,level,blockNo) = &
            & DOT_PRODUCT(flux_cc(    &
            & operators_coefficients%upwind_cell_idx(edge_index,level,blockNo), level, &
            & operators_coefficients%upwind_cell_blk(edge_index,level,blockNo))%x,     &
            & patch_2d%edges%primal_cart_normal(edge_index,blockNo)%x)
          
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
    INTEGER :: edge_index, level, blockNo      !< index of edge, vert level, block
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
    !!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      edge_upwind_flux(:,:,blockNo) = 0.0_wp
      DO edge_index = start_index, end_index
        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)
          !
          ! compute the central flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          flux_mean_e%x=0.5_wp*(flux_cc(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1))%x&
            & +flux_cc(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2))%x)
          
          edge_upwind_flux(edge_index,level,blockNo)=DOT_PRODUCT(flux_mean_e%x,&
            & patch_2d%edges%primal_cart_normal(edge_index,blockNo)%x)
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
  SUBROUTINE upwind_hflux_oce( patch_3d, cell_value, edge_vn, edge_upwind_flux, opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)              :: cell_value   (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)      !< advected cell centered variable
    REAL(wp), INTENT(in)              :: edge_vn    (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)       !< normal velocity on edges
    REAL(wp), INTENT(inout)             :: edge_upwind_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level    ! optional vertical end level
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER :: start_level, end_level
    INTEGER :: start_index, end_index
    INTEGER :: edge_index, level, blockNo         !< index of edge, vert level, block
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
!ICON_OMP_PARALLEL PRIVATE(iilc, iibc)
    ! line and block indices of two neighboring cells
    iilc => patch_2d%edges%cell_idx
    iibc => patch_2d%edges%cell_blk
    
    ! loop through all patch edges (and blocks)
!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      edge_upwind_flux(:,:,blockNo) = 0.0_wp
      DO edge_index = start_index, end_index
        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          edge_upwind_flux(edge_index,level,blockNo) =  &
             0.5_wp * (        edge_vn(edge_index,level,blockNo)  *           &
               & ( cell_value(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1)) + &
               &   cell_value(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2)) ) &
               &   - ABS( edge_vn(edge_index,level,blockNo) ) *               &
               & ( cell_value(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2)) - &
               &   cell_value(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1)) ) )
           ! inlined above
!          FUNCTION laxfr_upflux( p_vn, p_psi1, p_psi2 )  result(p_upflux)
!             & laxfr_upflux( edge_vn(edge_index,level,blockNo), cell_value(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1)), &
!             & cell_value(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2)) )
          
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    
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
!<Optimize:inUse>
  SUBROUTINE central_hflux_oce( patch_3d, cell_value, edge_vn, edge_flux )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)           :: cell_value   (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    !< advected cell centered variable
    REAL(wp), INTENT(inout)           :: edge_vn    (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     !< normal velocity on edges
    REAL(wp), INTENT(inout)           :: edge_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) !< variable in which the upwind flux is stored
    
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER :: start_index, end_index
    INTEGER :: edge_index, level, blockNo
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
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      edge_flux(:,:,blockNo) = 0.0_wp
      DO edge_index = start_index, end_index
        DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo)
          !
          ! compute the central flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          edge_flux(edge_index,level,blockNo) =  0.5_wp * edge_vn(edge_index,level,blockNo)    &
            & * (   cell_value(iilc(edge_index,blockNo,1), level, iibc(edge_index,blockNo,1)) &
            &     + cell_value(iilc(edge_index,blockNo,2), level, iibc(edge_index,blockNo,2)))
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
!ICON_OMP_END_PARALLE_DO
    
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
  SUBROUTINE lax_friedrichs_hflux_oce( patch_3d, cell_value, edge_vn, edge_flux )    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)              :: cell_value   (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< advected cell centered variable
    REAL(wp), INTENT(inout)              :: edge_vn    (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)              :: edge_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER :: start_index, end_index
    INTEGER :: edge_index, level, blockNo
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
    !!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      edge_flux(:,:,blockNo) = 0.0_wp
      DO edge_index = start_index, end_index
        DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo)
          !
          edge_flux(edge_index,level,blockNo) =  0.5_wp * edge_vn(edge_index,level,blockNo)    &
            & * (   cell_value(iilc(edge_index,blockNo,1), level, iibc(edge_index,blockNo,1)) &
            & + cell_value(iilc(edge_index,blockNo,2), level, iibc(edge_index,blockNo,2)))    &
            &          +0.5_wp*edge_vn(edge_index,level,blockNo)*edge_vn(edge_index,level,blockNo)*dtime&
            &          * patch_2d%edges%inv_dual_edge_length(edge_index,blockNo)   &
            &        *( cell_value(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2))      &
            &          -cell_value(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1)))
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
  SUBROUTINE miura_order1_hflux_oce( patch_3d,cell_value, edge_vn, &
    & operators_coefficients, edge_upwind_flux, &
    & opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)           :: cell_value(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< advected cell centered variable
    REAL(wp), INTENT(inout)           :: edge_vn (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) !< normal velocity on edges
    TYPE(t_operator_coeff)            :: operators_coefficients
    REAL(wp), INTENT(inout)           :: edge_upwind_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)!< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level      ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level      ! optional vertical end level
    
    ! local variables
    INTEGER :: start_level, end_level
    INTEGER :: start_index, end_index
    INTEGER :: edge_index, level, blockNo
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
    ! !     DO blockNo = i_startblk_e, i_endblk_e
    ! !       CALL get_indices_e(patch_2d, blockNo, i_startblk_e, i_endblk_e,   &
    ! !         &                start_index, end_index, rl_start,rl_end)
    ! !
    ! !       DO level = start_level, end_level
    ! !         DO edge_index = start_index, end_index          !
    ! !           IF ( v_base%lsm_e(edge_index,level,blockNo) <= sea_boundary ) THEN
    ! !              edge_upwind_flux(edge_index,level,blockNo) =  &
    ! !              &  laxfr_upflux( edge_vn(edge_index,level,blockNo), cell_value(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1)), &
    ! !              &                                 cell_value(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2)) )
    ! ! !write(*,*)'upwind miura',edge_index,level,blockNo,edge_upwind_flux(edge_index,level,blockNo),edge_upwind_flux(edge_index,level,blockNo), edge_vn(edge_index,level,blockNo)!,&
    ! ! !&cell_value(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1)), cell_value(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2))
    ! !           ELSE
    ! !             edge_upwind_flux(edge_index,level,blockNo) = 0.0_wp
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
    
    CALL grad_fd_norm_oce_3d( cell_value,                 &
      & patch_3d,    &
      & operators_coefficients%grad_coeff,  &
      & z_gradc)
    
    CALL sync_patch_array(sync_e, patch_2d, z_gradc)
    !3b:
    CALL map_edges2cell_3d( patch_3d, z_gradc,operators_coefficients, z_gradc_cc)
    
    !!ICON_OMP_PARALLEL
    !!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level, il_c, ib_c)  ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      edge_upwind_flux(:,:,blockNo) = 0.0_wp
      DO edge_index = start_index, end_index
        DO level = start_level, patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo)
        
          il_c = operators_coefficients%upwind_cell_idx(edge_index,level,blockNo)
          ib_c = operators_coefficients%upwind_cell_blk(edge_index,level,blockNo)
          
          edge_upwind_flux(edge_index,level,blockNo) = edge_vn(edge_index,level,blockNo)&
            & *( cell_value(il_c, level, ib_c) +          &
            & DOT_PRODUCT(                    &
            & operators_coefficients%moved_edge_position_cc(edge_index,level,blockNo)%x -    &
            & operators_coefficients%upwind_cell_position_cc(edge_index,level,blockNo)%x,  &
            & z_gradc_cc(il_c, level, ib_c)%x))
        END DO
      END DO
    END DO
    !!ICON_OMP_END_DO NOWAIT
    !!ICON_OMP_END_PARALLEL
    
  END SUBROUTINE miura_order1_hflux_oce
  !-------------------------------------------------------------------------


  
!  !-------------------------------------------------------------------------
!  !>
!  !! Flux limiter for horizontal advection
!  !!
!  !! Zalesak Flux-Limiter (Flux corrected transport)
!  !! The corrected flux is a weighted average of the low order flux and the
!  !! given high order flux. The high order flux is used to the greatest extent
!  !! possible without introducing overshoots and undershoots.
!  !! In vicinity of a lateral boundary only the low order flux is used: The criterion 
!  !! is that at least one of the edges of the two neighboring cells of
!  !! a central edges is a boundary edge.
!  !! Note: This limiter is positive definite and almost monotone (but not strictly).
!  !!
!  !! @par Literature:
!  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
!  !!   Algorithms for Fluids. JCP, 31, 335-362
!  !!
!  !! @par Revision History
!  !! - Inital revision by Daniel Reinert, DWD (2010-03-10)
!  !! Modification by Daniel Reinert, DWD (2010-03-25)
!  !! - adapted for MPI parallelization
!  !! - adapted for ocean use by P. Korn (2012)
!  !! - optimized by L. Linardakis (2015)
!  !! - criterion for switch to low order scheme near boundaries added P. Korn (2015)  
!  !!
!  !!  mpi note: computed on domain edges. Results is not synced.
!  !!
!!<Optimize:inUse>
!  SUBROUTINE hflx_limiter_oce_zalesak( patch_3d, vert_velocity, &
!    & tracer,              &
!    & p_mass_flx_e,      &
!    & flx_tracer_low,    &    
!    & flx_tracer_high,   &
!    & flx_tracer_final,  &
!    & operators_coefficients,        &
!    & h_old,             &
!    & h_new,             &
!    & zlim,tracer_index  )
!    
!    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
!    REAL(wp),INTENT(inout)              :: vert_velocity(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
!    REAL(wp), INTENT(inout)             :: tracer           (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!    REAL(wp), INTENT(inout)             :: p_mass_flx_e     (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
!    REAL(wp), INTENT(inout)             :: flx_tracer_low   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
!    REAL(wp), INTENT(inout)             :: flx_tracer_high  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
!    REAL(wp), INTENT(inout)             :: flx_tracer_final (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
!    TYPE(t_operator_coeff),INTENT(in)   :: operators_coefficients
!    REAL(wp), INTENT(in)                :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!    REAL(wp), INTENT(in)                :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!    REAL(wp), INTENT(inout)             :: zlim(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
!    INTEGER, INTENT(in) :: tracer_index
!    !Local variables
!    !REAL(wp)              :: flx_tracer_high2  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
!    REAL(wp) :: z_mflx_anti(patch_3d%p_patch_2d(1)%cells%max_connectivity)
!    REAL(wp) :: z_fluxdiv_c     !< flux divergence at cell center
!    REAL(wp) :: z_anti          (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)          !< antidiffusive tracer mass flux (F_H - F_L)    
!    REAL(wp) :: z_tracer_new_low(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after transport, if low order fluxes are used
!    REAL(wp) :: z_tracer_max    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local maximum of current tracer value and low order update
!    REAL(wp) :: z_tracer_min    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local minimum of current tracer value and low order update
!    REAL(wp) :: r_p             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< fraction which must multiply all in/out fluxes of cell jc to guarantee
!    REAL(wp) :: r_m             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< no overshoot/undershoot
!    REAL(wp) :: r_frac          !< computed minimum fraction which must multiply< the flux at the edge
!    REAL(wp) :: z_min, z_max    !< minimum/maximum value in cell and neighboring cells
!    REAL(wp) :: z_signum        !< sign of antidiffusive velocity
!    REAL(wp) :: p_p, p_m        !< sum of antidiffusive fluxes into and out of cell jc
!    REAL(wp) :: prism_thick_old(n_zlev), inv_prism_thick_new(n_zlev)
!    REAL(wp) :: delta_z_new, delta_z
!    INTEGER, DIMENSION(:,:,:), POINTER ::  cellOfEdge_idx, cellOfEdge_blk
!    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk
!    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
!    INTEGER :: start_level, end_level            
!    INTEGER :: start_index, end_index
!!     INTEGER :: cell_1_idx,cell_2_idx,cell_1_blk,cell_2_blk
!!     INTEGER :: cell_1_edge_1_idx,cell_1_edge_2_idx,cell_1_edge_3_idx
!!     INTEGER :: cell_2_edge_1_idx,cell_2_edge_2_idx,cell_2_edge_3_idx        
!!     INTEGER :: cell_1_edge_1_blk,cell_1_edge_2_blk,cell_1_edge_3_blk
!!     INTEGER :: cell_2_edge_1_blk,cell_2_edge_2_blk,cell_2_edge_3_blk      
!    INTEGER :: edge_index, level, blockNo, jc,  cell_connect, sum_lsm_quad_edge
!!    INTEGER :: all_water_edges 
!    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
!    TYPE(t_patch), POINTER :: patch_2d    
!    REAL(wp) :: z_adv_flux_v(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)  ! upwind flux
!    REAL(wp) :: flux_div_vert !new tracer
!    !-------------------------------------------------------------------------
!    patch_2d        => patch_3d%p_patch_2d(1)
!    edges_in_domain => patch_2d%edges%in_domain
!    cells_in_domain => patch_2d%cells%in_domain
!    !-------------------------------------------------------------------------
!    start_level = 1
!    end_level   = n_zlev
!    cellOfEdge_idx  => patch_2d%edges%cell_idx
!    cellOfEdge_blk  => patch_2d%edges%cell_blk
!    edge_of_cell_idx  => patch_2d%cells%edge_idx
!    edge_of_cell_blk  => patch_2d%cells%edge_blk
!    neighbor_cell_idx => patch_2d%cells%neighbor_idx
!    neighbor_cell_blk => patch_2d%cells%neighbor_blk
!    
!            
!    ! 1. Calculate low (first) order fluxes using the standard upwind scheme and the
!    !    antidiffusive fluxes
!    !    (not allowed to call upwind_hflux_up directly, due to circular dependency)
!!     IF ( p_test_run ) THEN
!! !       z_tracer_max    (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
!! !       z_tracer_min    (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
!!       r_m             (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
!!       r_p             (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
!!     ENDIF
!#ifdef NAGFOR
!    z_tracer_max(:,:,:) = 0.0_wp
!    z_tracer_min(:,:,:) = 0.0_wp
!    r_m(:,:,:)          = 0.0_wp
!    r_p(:,:,:)          = 0.0_wp
!#endif
!
!      CALL upwind_vflux_oce( patch_3d, &
!        & tracer,                      &
!        & vert_velocity,               &
!        & z_adv_flux_v )
!        
!      !CALL sync_patch_array(sync_c, patch_2d, z_adv_flux_v)
!   
!   
!!ICON_OMP_PARALLEL
!! !ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE        
!!       DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
!!       
!!         flux_div_vert(:,:,blockNo) = 0.0_wp     
!!         
!!         CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
!!             
!!         DO jc = start_index, end_index
!!    
!!           DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)        
!!             ! positive vertical divergence in direction of w (upward positive)
!!             flux_div_vert(jc,level,blockNo) = z_adv_flux_v(jc, level, blockNo) &
!!             & - z_adv_flux_v(jc, level+1, blockNo)
!!           ENDDO
!!         END DO
!!       END DO
!! !ICON_OMP_END_DO      
!! !       CALL sync_patch_array(sync_c, patch_2D, flux_div_vert)
! 
!
!!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
!    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
!      
!      z_anti(:,:,blockNo)     = 0.0_wp
!      DO edge_index = start_index, end_index
!        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)
!          
!          ! calculate antidiffusive flux for each edge
!          z_anti(edge_index,level,blockNo) = flx_tracer_high(edge_index,level,blockNo)&
!                                          &- flx_tracer_low(edge_index,level,blockNo)
!        END DO  ! end loop over edges
!      END DO  ! end loop over levels
!    END DO  ! end loop over blocks
!!ICON_OMP_END_DO
!
!    
!!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, delta_z, delta_z_new, &
!!ICON_OMP z_fluxdiv_c, flux_div_vert) ICON_OMP_DEFAULT_SCHEDULE
!    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
!      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
!      
!      z_tracer_new_low(:,:,blockNo)  = 0.0_wp
!      
!      DO jc = start_index, end_index
!        IF (patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo) < start_level) CYCLE 
!        ! get prism thickness
!        !inv_prism_thick_new(start_level) = 1.0_wp / (patch_3d%p_patch_1d(1)%del_zlev_m(start_level) + h_new(jc,blockNo))
!        !prism_thick_old(start_level)     = patch_3d%p_patch_1d(1)%del_zlev_m(start_level)           + h_old(jc,blockNo)
!        !DO level = start_level+1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)
!        !  prism_thick_old (level)    = patch_3d%p_patch_1d(1)%del_zlev_m(level)
!        !  inv_prism_thick_new(level) = patch_3d%p_patch_1d(1)%inv_del_zlev_m(level)
!        !ENDDO
!        
!        ! 3. Compute the complete (with horizontal and vertical divergence) updated low order solution z_tracer_new_low
!        ! First at top level than in fluid interior
!        !       
!        level = start_level
!        !  compute divergence of low order fluxes
!        z_fluxdiv_c = 0
!        DO cell_connect = 1, patch_2d%cells%num_edges(jc,blockNo)
!          z_fluxdiv_c =  z_fluxdiv_c + &
!            & flx_tracer_low(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect)) * &
!            & operators_coefficients%div_coeff(jc,level,blockNo,cell_connect)
!        ENDDO
!
!        ! note: the is = - z_adv_flux_v(jc, 2, blockNo), as the the interface is 0.
!        flux_div_vert = z_adv_flux_v(jc, level, blockNo) &
!          & - z_adv_flux_v(jc, level+1, blockNo)
!
!
!        delta_z = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)&
!             &  + h_old(jc,blockNo)
!        delta_z_new = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)&
!             &  + h_new(jc,blockNo)
!             
!        ! Low order flux at top level
!        z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
!          & - dtime * (z_fluxdiv_c+flux_div_vert))/delta_z_new
!        !
!!         CALL global_mpi_barrier()
!!         write(0,*) "Over the barrier"
!       !Fluid interior       
!        DO level = start_level+1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)       
!          !  compute divergence of low order fluxes
!          z_fluxdiv_c = 0
!          DO cell_connect = 1, patch_2d%cells%num_edges(jc,blockNo)
!            z_fluxdiv_c =  z_fluxdiv_c + &
!              & flx_tracer_low(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect)) * &
!              & operators_coefficients%div_coeff(jc,level,blockNo,cell_connect)
!          ENDDO
!
!          flux_div_vert = z_adv_flux_v(jc, level, blockNo) &
!            & - z_adv_flux_v(jc, level+1, blockNo)
!
!          delta_z     = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)
!          delta_z_new = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)
!          !
!          ! low order flux in flow interior
!          !z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * prism_thick_old(level)     &
!          !  & - dtime * (z_fluxdiv_c+flux_div_vert(jc,level,blockNo))) * inv_prism_thick_new(level)
!          z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
!            & - dtime * (z_fluxdiv_c+flux_div_vert))/delta_z_new
!        ENDDO
!      ENDDO
!      
!      ! precalculate local maximum/minimum of current tracer value and low order
!      ! updated value
!      z_tracer_max(:,:,blockNo) =            &
!        & MAX(          tracer(:,:,blockNo), &
!        &     z_tracer_new_low(:,:,blockNo))
!      z_tracer_min(:,:,blockNo) =            &
!        & MIN(          tracer(:,:,blockNo), &
!        &     z_tracer_new_low(:,:,blockNo))
!
!!      write(0,*) blockNo, ":", z_tracer_max(start_index:end_index,start_level:end_level,blockNo)
!!      write(0,*) blockNo, ":", z_tracer_min(start_index:end_index,start_level:end_level,blockNo)
!    ENDDO
!!ICON_OMP_END_DO
!
!
!!ICON_OMP_MASTER
!    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, z_tracer_max, z_tracer_min)
!!ICON_OMP_END_MASTER    
!!ICON_OMP_BARRIER
!    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
!    !    field is free of any new extrema.    
!!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, inv_prism_thick_new, &
!!ICON_OMP z_mflx_anti, z_max, z_min, cell_connect, p_p, p_m) ICON_OMP_DEFAULT_SCHEDULE
!    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
!
!      ! this is only needed for the parallel test setups
!      ! it will try  tocheck the uninitialized (land) parts
!      r_m(:,:,blockNo) = 0.0_wp
!      r_p(:,:,blockNo) = 0.0_wp
!        
!      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
!      DO jc = start_index, end_index
!        
!        ! get prism thickness
!        inv_prism_thick_new(start_level) = 1.0_wp / (patch_3d%p_patch_1d(1)%del_zlev_m(start_level)+ h_new(jc,blockNo))
!        DO level = start_level+1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)
!          inv_prism_thick_new(level) = patch_3D%p_patch_1d(1)%inv_prism_thick_c(jc,level,blockNo)
!        ENDDO
!        
!        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)         
!          ! 2. Define "antidiffusive" fluxes A(jc,level,blockNo,edge_index) for each cell. It is the difference
!          !    between the high order fluxes (given by the FFSL-scheme) and the low order
!          !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
!          !    - positive for outgoing fluxes
!          !    - negative for incoming fluxes
!          !    this sign convention is related to the definition of the divergence operator.          
!          z_mflx_anti(:) = 0.0_wp
!          z_max = z_tracer_max(jc,level,blockNo)
!          z_min = z_tracer_min(jc,level,blockNo)
!          p_p = 0.0_wp
!          p_m = 0_wp
!          DO cell_connect = 1, patch_2d%cells%num_edges(jc,blockNo)
!            IF (patch_3d%p_patch_1d(1)% &
!              & dolic_c(neighbor_cell_idx(jc,blockNo,cell_connect), neighbor_cell_blk(jc,blockNo,cell_connect)) >= level) THEN
!              
!              z_max = MAX(z_max, &
!                & z_tracer_max(neighbor_cell_idx(jc,blockNo,cell_connect),level,neighbor_cell_blk(jc,blockNo,cell_connect)))
!              z_min = MIN(z_min, &
!                & z_tracer_min(neighbor_cell_idx(jc,blockNo,cell_connect),level,neighbor_cell_blk(jc,blockNo,cell_connect)))
!
!              z_mflx_anti(cell_connect) =                                                        &
!                & dtime * operators_coefficients%div_coeff(jc,level,blockNo,cell_connect) * inv_prism_thick_new(level)  &
!                & * z_anti(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect))
!
!              ! Sum of all incoming antidiffusive fluxes into cell jc
!              ! outgoing fluxes carry a positive sign, incoming a negative
!              p_p = p_p - MIN(0._wp, z_mflx_anti(cell_connect))
!              ! Sum of all outgoing antidiffusive fluxes out of cell jc
!              p_m = p_m + MAX(0._wp, z_mflx_anti(cell_connect))
!            ENDIF
!          ENDDO                
!          ! fraction which must multiply all fluxes out of cell jc to guarantee no
!          ! undershoot
!          ! Nominator: maximum allowable decrease of tracer
!          r_m(jc,level,blockNo) = (z_tracer_new_low(jc,level,blockNo) - z_min ) / (p_m + dbl_eps)!&
!          !
!          ! fraction which must multiply all fluxes into cell jc to guarantee no
!          ! overshoot
!          ! Nominator: maximum allowable increase of tracer
!          r_p(jc,level,blockNo) = (z_max - z_tracer_new_low(jc,level,blockNo)) / (p_p + dbl_eps)!&
!          !
!          !update old tracer with low-order flux
!          tracer(jc,level,blockNo)=z_tracer_new_low(jc,level,blockNo)         
!        ENDDO
!      ENDDO
!    ENDDO
!!ICON_OMP_END_DO
!    
!!ICON_OMP_MASTER
!    ! Synchronize r_m and r_p
!    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, r_m, r_p)
!!ICON_OMP_END_MASTER
!!ICON_OMP_BARRIER   
!
!    ! 5. Now loop over all edges and determine the minimum fraction which must
!    !    multiply the antidiffusive flux at the edge.
!    !    At the end, compute new, limited fluxes which are then passed to the main
!!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level, z_signum, r_frac) ICON_OMP_DEFAULT_SCHEDULE
!    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
!      flx_tracer_final(:,:,blockNo) = 0.0_wp
!      DO edge_index = start_index, end_index
!      
!        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)
!        
!          IF( operators_coefficients%edges_SeaBoundaryLevel(edge_index,level,blockNo) > -2)THEN! edge < 2nd order boundary
!          
!            flx_tracer_final(edge_index,level,blockNo) = flx_tracer_low(edge_index,level,blockNo)
!            
!          ELSE!IF(sum_lsm_quad_edge==all_water_edges)THEN
!          
!            !z_anti>0 returns  1: here z_anti is outgoing, i.e. flux_high>flux_low
!            !z_anti<0 returns -1: here z_anti is ingoing, i.e. flux_high<flux_low
!            z_signum = SIGN(1._wp, z_anti(edge_index,level,blockNo))
!                    
!          ! This does the same as an IF (z_signum > 0) THEN ... ELSE ... ENDIF,
!          ! but is computationally more efficient
!          r_frac = 0.5_wp * (       &
!            & (1._wp + z_signum) * & !<- active for z_signum=1
!            & MIN(r_m(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)),  &
!            &     r_p(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)))  &
!            &+(1._wp - z_signum) * & !<- active for z_signum=-1
!            & MIN(r_m(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)),  &
!            &     r_p(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)))  )
!
!          if (tracer_index == 1) zlim(edge_index,level,blockNo)=r_frac
!          
!          ! Limited flux
!          flx_tracer_final(edge_index,level,blockNo) = flx_tracer_low(edge_index,level,blockNo)&
!           & + MIN(1.0_wp,r_frac) *z_anti(edge_index,level,blockNo)      
!                  
!            ENDIF
!        END DO
!       ENDDO
!    ENDDO
!!ICON_OMP_END_DO NOWAIT
!!ICON_OMP_END_PARALLEL
!
!
!  END SUBROUTINE hflx_limiter_oce_zalesak
!  !-------------------------------------------------------------------------
!
!
!  !-------------------------------------------------------------------------
!  !>
!  !! Positive definite flux limiter for horizontal advection
!  !!
!  !! Positive definite Zalesak Flux-Limiter (Flux corrected transport).
!  !! Only outward fluxes are re-scaled, in order to maintain positive
!  !! definiteness.
!  !!
!  !! @par Literature:
!  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
!  !!   Algorithms for Fluids. JCP, 31, 335-362
!  !! - Harris, L. M. and P. H. Lauritzen (2010): A flux-form version of the
!  !!   Conservative Semi-Lagrangian Multi-tracer transport scheme (CSLAM) on
!  !!   the cubed sphere grid. JCP, in press
!  !!
!  !! @par Revision History
!  !! - Inital revision by Daniel Reinert, DWD (2010-10-06)
!  !!
!  SUBROUTINE hflx_limiter_oce_posdef( patch_3d,&
!    & tracer,                             &
!    & flx_tracer,                         &
!    & flx_tracer_limit,                   &
!    & operators_coefficients )
!    
!    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
!    REAL(wp), INTENT(inout)             :: tracer      (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!    REAL(wp), INTENT(inout)             :: flx_tracer(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
!    REAL(wp), INTENT(inout)             :: flx_tracer_limit(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
!    TYPE(t_operator_coeff),INTENT(in)   :: operators_coefficients
!
!    !Local variables
!    REAL(wp) :: z_mflx(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)   
!    REAL(wp) :: r_m             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< no overshoot/undershoot
!    REAL(wp) :: z_signum        !< sign of antidiffusive velocity
!    REAL(wp) :: p_m        !< sum of antidiffusive fluxes into and out of cell jc
!    INTEGER, DIMENSION(:,:,:), POINTER :: cellOfEdge_idx, cellOfEdge_blk
!    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk
!    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
!    INTEGER :: start_level, end_level            
!    INTEGER :: start_index, end_index
!    INTEGER :: edge_index, level, blockNo, jc!,  cell_connect
!    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
!    TYPE(t_patch), POINTER :: patch_2d    
!    REAL(wp) :: z_adv_flux_v(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)  ! upwind flux
!    REAL(wp) :: flux_div_vert(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !new tracer        
!    !-------------------------------------------------------------------------
!    patch_2d        => patch_3d%p_patch_2d(1)
!    edges_in_domain => patch_2d%edges%in_domain
!    cells_in_domain => patch_2d%cells%in_domain
!    !-------------------------------------------------------------------------
!    start_level = 1
!    end_level   = n_zlev
!    
!    flx_tracer_limit(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%nblks_e)          = 0.0_wp
!    z_mflx          (1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)= 0.0_wp 
!    r_m             (1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)= 0.0_wp  
!    
!    ! Set pointers to index-arrays
!    ! line and block indices of two neighboring cells
!    cellOfEdge_idx  => patch_2d%edges%cell_idx
!    cellOfEdge_blk  => patch_2d%edges%cell_blk
!    edge_of_cell_idx  => patch_2d%cells%edge_idx
!    edge_of_cell_blk  => patch_2d%cells%edge_blk
!    neighbor_cell_idx => patch_2d%cells%neighbor_idx
!    neighbor_cell_blk => patch_2d%cells%neighbor_blk
!
!    
!    !
!    ! 1. Reformulate all fluxes in terms of the total mass [kg m^-3]
!    !    that crosses each of the CV-edges and store them in a cell-based structure.
!    !
!    !    z_mflx > 0: outward
!    !    z_mflx < 0: inward
!    !    
!    !!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, inv_prism_thick_new, prism_thick_old, &
!    !!ICON_OMP z_fluxdiv_c ) ICON_OMP_DEFAULT_SCHEDULE
!    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
!      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
!      
!      DO jc = start_index, end_index
!                
!        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)
!        
!          z_mflx(jc,level,1) = dtime*operators_coefficients%div_coeff(jc,level,blockNo,1) &
!            & * flx_tracer(edge_of_cell_idx(jc,blockNo,1),level,edge_of_cell_blk(jc,blockNo,1))
!            
!
!          z_mflx(jc,level,2) = dtime*operators_coefficients%div_coeff(jc,level,blockNo,2)&
!            & * flx_tracer(edge_of_cell_idx(jc,blockNo,2),level,edge_of_cell_blk(jc,blockNo,2))
!
!  
!          z_mflx(jc,level,3) = dtime*operators_coefficients%div_coeff(jc,level,blockNo,3)&
!            & * flx_tracer(edge_of_cell_idx(jc,blockNo,3),level,edge_of_cell_blk(jc,blockNo,3))
!
!          
!        ENDDO
!      ENDDO
!    ENDDO
!    !!ICON_OMP_END_DO NOWAIT
!    !!ICON_OMP_END_PARALLEL
!    CALL sync_patch_array(SYNC_C1, patch_2d, z_mflx)
!    ! 2. Compute total outward mass
!    !
!    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
!      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
!      DO jc = start_index, end_index
!        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)
!
!          ! Sum of all outgoing fluxes out of cell jc
!          p_m =  MAX(0._wp,z_mflx(jc,level,1))  &
!            &  + MAX(0._wp,z_mflx(jc,level,2))  &
!            &  + MAX(0._wp,z_mflx(jc,level,3))
!
!          ! fraction which must multiply all fluxes out of cell jc to guarantee no
!          ! undershoot
!          ! Nominator: maximum allowable decrease of \rho q
!          r_m(jc,level,blockNo) = MIN(1._wp, &
!            & (tracer(jc,level,blockNo) &
!            &  * patch_3d%p_patch_1d(1)%prism_thick_c(jc,level,blockNo)) &
!            &  /(p_m + dbl_eps) )
!
!        ENDDO
!      ENDDO
!    ENDDO
!   ! synchronize r_m
!    CALL sync_patch_array(SYNC_C1, patch_2d, r_m)
!
!    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
!      DO edge_index = start_index, end_index
!        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)
!
!          ! p_mflx_tracer_h > 0: flux directed from cell 1 -> 2
!          ! p_mflx_tracer_h < 0: flux directed from cell 2 -> 1
!          z_signum = SIGN(1._wp,flx_tracer(edge_index,level,blockNo))
!          flx_tracer_limit(edge_index,level,blockNo) = flx_tracer(edge_index,level,blockNo) * 0.5_wp  &
!            & * ( (1._wp + z_signum) &
!            &      * r_m(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)) &
!            &    + (1._wp - z_signum) &
!            &      * r_m(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)) )
!
!        END DO
!      ENDDO
!    ENDDO
!    CALL sync_patch_array(sync_e, patch_2d, flx_tracer_limit)
!  !---------DEBUG DIAGNOSTICS-------------------------------------------
!    idt_src=3  ! output print level (1-5, fix)
!    !CALL dbg_print('LimitPosDef: scalefac'   ,r_m,str_module,idt_src,patch_2d%cells%owned)    
!    CALL dbg_print('LimitPosDef: flux_h'     ,flx_tracer_limit,str_module,idt_src,patch_2d%edges%owned)
!    !---------------------------------------------------------------------
!            
!  END SUBROUTINE hflx_limiter_oce_posdef
!  !-------------------------------------------------------------------------
!  
  !-------------------------------------------------------------------------
  !>
  !! Lax Friedrichs first order upwind flux,
  !! used in conservative advection routines.
  !! For passive advection, equivalent to
  !! any other first order upwind flux.
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura  (2004).
  !!
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
  ! ! SUBROUTINE elad(patch_2d, trac_old, trac_new,operators_coefficients, p_os)
  ! ! !
  ! ! !
  ! ! TYPE(t_patch), TARGET, INTENT(in) :: patch_2d
  ! ! REAL(wp), INTENT(IN)  :: trac_old(:,:,:)
  ! ! REAL(wp), INTENT(inout) :: trac_new(:,:,:)
  ! ! TYPE(t_operator_coeff), INTENT(IN)        :: operators_coefficients
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
  ! ! INTEGER  :: jc, level, blockNo!, jkp1        !< index of edge, vert level, block
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
  ! ! ! DO level=1,n_zlev
  ! ! !   write(*,*)'max/min tracer in:',level,&
  ! ! !   &maxval(trac_new(:,level,:)),minval(trac_new(:,level,:))
  ! ! ! END DO
  ! !
  ! !   !Step 1: determine minimal and maximal permissible values
  ! !   DO blockNo = i_startblk_c, i_endblk_c
  ! !     CALL get_indices_c( patch_2d, blockNo, i_startblk_c, i_endblk_c, start_index, end_index, &
  ! !     &                   rl_start_c, rl_end_c)
  ! !
  ! !     DO jc = start_index, end_index
  ! !       z_dolic = v_base%dolic_c(jc,blockNo)
  ! !       IF(z_dolic>=MIN_DOLIC)THEN
  ! !         DO level = 1, z_dolic
  ! !           DO ie=1,no_cell_edges
  ! !
  ! !             !actual edges of cell c1
  ! !             il_e(ie) = patch_2d%cells%edge_idx(jc,blockNo,ie)
  ! !             ib_e(ie) = patch_2d%cells%edge_blk(jc,blockNo,ie)
  ! !
  ! !             !get neighbor cells of edge
  ! !             il_c1 = patch_2d%edges%cell_idx(il_e(ie),ib_e(ie),1)
  ! !             ib_c1 = patch_2d%edges%cell_blk(il_e(ie),ib_e(ie),1)
  ! !             il_c2 = patch_2d%edges%cell_idx(il_e(ie),ib_e(ie),2)
  ! !             ib_c2 = patch_2d%edges%cell_blk(il_e(ie),ib_e(ie),2)
  ! ! !             IF(z_in(il_c1,level,ib_c1)>z_in(il_c2,level,ib_c2))THEN
  ! ! !               z_up(jc,level,blockNo,ie)   = z_in(il_c1,level,ib_c1)
  ! ! !               z_down(jc,level,blockNo,ie) = z_in(il_c2,level,ib_c2)
  ! ! !             ELSE
  ! ! !               z_up(jc,level,blockNo,ie)   = z_in(il_c2,level,ib_c2)
  ! ! !               z_down(jc,level,blockNo,ie) = z_in(il_c1,level,ib_c1)
  ! ! !             ENDIF
  ! !               IF(p_os%p_diag%ptp_vn(il_e(ie),level,ib_e(ie))>=0.0_wp)THEN
  ! !                 z_up(jc,level,blockNo,ie) = trac_old(il_c1,level,ib_c1)
  ! !               ELSEIF(p_os%p_diag%ptp_vn(il_e(ie),level,ib_e(ie))<0.0_wp)THEN
  ! !                 z_up(jc,level,blockNo,ie) = trac_old(il_c2,level,ib_c2)
  ! !               ENDIF
  ! !          END DO
  ! !          IF ( v_base%lsm_c(jc,level,blockNo) <= sea_boundary ) THEN
  ! !            z_max(jc,level,blockNo) = maxval(z_up(jc,level,blockNo,1:no_cell_edges))
  ! !            z_min(jc,level,blockNo) = minval(z_up(jc,level,blockNo,1:no_cell_edges))
  ! !          ELSE
  ! !            z_max(jc,level,blockNo) = 0.0_wp
  ! !            z_min(jc,level,blockNo) = 0.0_wp
  ! !          ENDIF
  ! ! !         z_max(jc,level,blockNo) = maxval(z_up(jc,level,blockNo,1:no_cell_edges))
  ! ! !         z_min(jc,level,blockNo) = minval(z_down(jc,level,blockNo,1:no_cell_edges))
  ! !
  ! !         END DO!level-loop
  ! !       ENDIF!(z_dolic>0)
  ! !     END DO!idx-loop
  ! !   END DO !blk-loop
  ! !
  ! ! ! DO level=1,n_zlev
  ! ! ! write(*,*)'admissible bound',level,&
  ! ! ! &maxval(z_max(:,level,:)), minval(z_max(:,level,:)),&
  ! ! ! &maxval(z_min(:,level,:)), minval(z_min(:,level,:))
  ! ! ! END DO
  ! !
  ! !
  ! ! ITERATION_LOOP: Do iter = 1, i_max_iter
  ! ! !write(*,*)'iteration----------',iter
  ! ! !write(1230,*)'iteration----------',iter
  ! !   !Step 1: determine excess field
  ! !   DO blockNo = i_startblk_c, i_endblk_c
  ! !     CALL get_indices_c( patch_2d, blockNo, i_startblk_c, i_endblk_c, start_index, end_index, &
  ! !     &                   rl_start_c, rl_end_c)
  ! !
  ! !     DO jc = start_index, end_index
  ! !       z_dolic = v_base%dolic_c(jc,blockNo)
  ! !       IF(z_dolic>MIN_DOLIC)THEN
  ! !         DO level = 1, z_dolic
  ! !          IF ( v_base%lsm_c(jc,level,blockNo) <= sea_boundary ) THEN
  ! !            z_excess(jc,level,blockNo) = max(z_pred(jc,level,blockNo)-z_max(jc,level,blockNo),0.0_wp)&
  ! !                               &+min(z_pred(jc,level,blockNo)-z_min(jc,level,blockNo),0.0_wp)
  ! !           ELSE
  ! !             z_excess(jc,level,blockNo) = 0.0_wp
  ! !           ENDIF
  ! !
  ! ! !  IF(z_excess(jc,level,blockNo)/=0.0)THEN
  ! ! !  write(1230,*)'excess field',jc,level,blockNo,z_excess(jc,level,blockNo),&
  ! ! ! &z_pred(jc,level,blockNo), z_max(jc,level,blockNo),z_up(jc,level,blockNo,1:no_cell_edges)
  ! ! !  ENDIF
  ! !         END DO!level-loop
  ! !       ENDIF!(z_dolic>0)
  ! !     END DO!idx-loop
  ! !   END DO !blk-loop
  ! !
  ! ! ! DO level=1,n_zlev
  ! ! ! write(*,*)'max-min excess',level,&
  ! ! ! &maxval(z_excess(:,level,:)), minval(z_excess(:,level,:))
  ! ! ! END DO
  ! !
  ! !
  ! !    !Step 3: Calculate diffusion of excess field
  ! !    !CALL grad_fd_norm_oce( z_excess, patch_2d, z_grad_excess)
  ! !     CALL grad_fd_norm_oce_3D( z_excess,               &
  ! !            &                  patch_2d,                &
  ! !            &                  operators_coefficients%grad_coeff,  &
  ! !            &                  z_grad_excess)
  ! !    z_grad_excess = z_K*z_grad_excess
  ! !    CALL div_oce_3D( z_grad_excess, patch_2d,operators_coefficients%div_coeff, z_diff_excess)
  ! !
  ! ! ! DO level=1,n_zlev
  ! ! ! write(*,*)'max-min diffusion',level,&
  ! ! ! &maxval(z_diff_excess(:,level,:)), minval(z_diff_excess(:,level,:))
  ! ! ! END DO
  ! !
  ! !   !Step 4
  ! !   DO blockNo = i_startblk_c, i_endblk_c
  ! !     CALL get_indices_c( patch_2d, blockNo, i_startblk_c, i_endblk_c, start_index, end_index, &
  ! !     &                   rl_start_c, rl_end_c)
  ! !
  ! !     DO jc = start_index, end_index
  ! !       z_dolic = v_base%dolic_c(jc,blockNo)
  ! !       IF(z_dolic>=MIN_DOLIC)THEN
  ! !         DO level = 1, z_dolic
  ! !           IF ( v_base%lsm_c(jc,level,blockNo) <= sea_boundary ) THEN
  ! !             z_out(jc,level,blockNo) = z_pred(jc,level,blockNo) - z_excess(jc,level,blockNo)+ z_diff_excess(jc,level,blockNo)
  ! !
  ! ! !  IF(z_excess(jc,level,blockNo)/=0.0)THEN
  ! ! !  write(1230,*)'alg', z_out(jc,level,blockNo),z_pred(jc,level,blockNo),&
  ! ! ! &z_excess(jc,level,blockNo), z_diff_excess(jc,level,blockNo), &
  ! ! ! &z_max(jc,level,blockNo),z_up(jc,level,blockNo,1:no_cell_edges)
  ! ! !  ENDIF
  ! !
  ! !           ELSE
  ! !             z_out(jc,level,blockNo) = 0.0_wp
  ! !           ENDIF
  ! ! !           IF(z_diff_excess(jc,1,blockNo)/=0.0_wp)THEN
  ! ! !             write(123,*)'correction',jc,blockNo,z_in(jc,1,blockNo),&
  ! ! !             & z_out(jc,1,blockNo), z_diff_excess(jc,1,blockNo)
  ! ! !           ENDIF
  ! !         END DO!level-loop
  ! !       ENDIF!(z_dolic>0)
  ! !     END DO!idx-loop
  ! !   END DO !blk-loop
  ! !
  ! !   !Step 4: Stop criterion
  ! !   stop_ctr = 0
  ! !   DO level=1,n_zlev
  ! !
  ! ! !      write(*,*)'actual state',level,&
  ! ! !       & maxval(z_out(:,level,:)),minval(z_out(:,level,:))
  ! ! !      write(*,*)' pred',level,&
  ! ! !       & maxval(z_pred(:,level,:)),minval(z_pred(:,level,:))
  ! !
  ! !     IF(   maxval(z_excess(:,level,:))<1.0E-12_wp&
  ! !     &.AND.minval(z_excess(:,level,:))<1.0E-12_wp)THEN
  ! !
  ! !       stop_ctr = stop_ctr +1
  ! !
  ! !     ELSE
  ! !     ENDIF
  ! !   END DO
  ! !
  ! !   IF(stop_ctr==n_zlev)THEN
  ! !     DO blockNo = i_startblk_c, i_endblk_c
  ! !       CALL get_indices_c( patch_2d, blockNo, i_startblk_c, i_endblk_c, start_index, end_index, &
  ! !       &                   rl_start_c, rl_end_c)
  ! !
  ! !       DO jc = start_index, end_index
  ! !         z_dolic = v_base%dolic_c(jc,blockNo)
  ! !         IF(z_dolic>=MIN_DOLIC)THEN
  ! !           DO level = 1, z_dolic
  ! !            IF ( v_base%lsm_c(jc,level,blockNo) <= sea_boundary ) THEN
  ! !            trac_new(jc,level,blockNo) = z_out(jc,level,blockNo)
  ! !            ENDIF
  ! ! !            write(1230,*)'before-after',jc,level,blockNo,&
  ! ! !           &trac_new(jc,level,blockNo),trac_old(jc,level,blockNo),z_excess(jc,level,blockNo)
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
  ! ! !   DO blockNo = i_startblk_c, i_endblk_c
  ! ! !     CALL get_indices_c( patch_2d, blockNo, i_startblk_c, i_endblk_c, start_index, end_index, &
  ! ! !     &                   rl_start_c, rl_end_c)
  ! ! !     DO jc = start_index, end_index
  ! ! ! write(123,*)'trac old new',trac_old(jc,1,blockNo),z_old_new(jc,1,blockNo),&
  ! ! ! &trac_new(jc,1,blockNo), z_out(jc,1,blockNo)
  ! ! !
  ! ! !    END DO
  ! ! !  END DO
  ! ! END SUBROUTINE elad
  ! !   !-------------------------------------------------------------------------
   !-----------------------------------------------------------------------
! !   !-------------------------------------------------------------------------
END MODULE mo_ocean_tracer_transport_horz
