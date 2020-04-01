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
  USE mo_math_types,                ONLY: t_cartesian_coordinates
  USE mo_ocean_nml,                 ONLY: n_zlev, l_edge_based, ab_gam,                   &
    & upwind, central,lax_friedrichs, horz_flux_twisted_vec_recon, miura_order1, flux_calculation_horz,      &
    & fct_high_order_flux,  fct_low_order_flux,FCT_Limiter_horz, fct_limiter_horz_zalesak,&
    & fct_limiter_horz_minmod, fct_limiter_horz_posdef, l_with_horz_tracer_diffusion, l_with_horz_tracer_advection,&
    &l_LAX_FRIEDRICHS, l_GRADIENT_RECONSTRUCTION, Tracer_HorizontalDiffusion_PTP_coeff, tracer_HorizontalAdvection_type, &
    & edge_based, cell_based
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma
  USE mo_run_config,                ONLY: dtime, ltimer
  USE mo_timer,                     ONLY: timer_start, timer_stop, timers_level, timer_adv_horz, timer_hflx_lim, &
    & timer_dif_horz, timer_extra10, timer_extra11, timer_extra12, timer_extra13, timer_extra15
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish !, message_text, message
!   USE mo_ocean_physics
  USE mo_scalar_product,            ONLY: map_cell2edges_3d,map_edges2cell_3d, &
    & map_edges2edges_viacell_3d_const_z
  USE mo_ocean_math_operators,      ONLY: div_oce_3d
  USE mo_ocean_tracer_diffusion,    ONLY: tracer_diffusion_horz
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff, no_primal_edges
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_c, sync_c1, sync_e, sync_patch_array, sync_patch_array_mult
  USE mo_mpi,                       ONLY: global_mpi_barrier
  USE mo_ocean_limiter,             ONLY: limiter_ocean_zalesak_horizontal
  USE mo_ocean_tracer_transport_types,  ONLY: t_ocean_transport_state
  
  
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=12)           :: str_module    = 'oceTracHorz '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  !
  ! PUBLIC INTERFACE
  !
  PUBLIC :: advect_horz
  PUBLIC :: diffuse_horz  
  
  INTEGER, PARAMETER :: top=1
  
CONTAINS
  !-----------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE advect_horz( patch_3d,&
    & trac_old,              &
    & transport_state,                  &
    & operators_coefficients,&
    & k_h,                   &
    & h_old,                 &
    & h_new,                 &
    & div_flux_horz,         &
    & div_flux_vert,         &
    & horizontally_diffused_tracer)

    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp)                               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_ocean_transport_state), TARGET      :: transport_state
    TYPE(t_operator_coeff), INTENT(in)     :: operators_coefficients
    REAL(wp), INTENT(in)                   :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                   :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                   :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)                :: div_flux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)                :: div_flux_vert(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)       
    REAL(wp), INTENT(inout), OPTIONAL      :: horizontally_diffused_tracer(:,:,:)    
    !
    !-------------------------------------------------------------------------------    
    start_timer(timer_adv_horz,2)

    SELECT CASE(tracer_HorizontalAdvection_type)
    
    CASE(cell_based)

      CALL advect_cell_based( patch_3d, &
      & trac_old,              &
      & transport_state,                  &
      & operators_coefficients,&
      & k_h,                   &
      & h_old,                 &
      & h_new,                 &
      & div_flux_horz,         &
      & div_flux_vert) 
    CASE default
      CALL finish("advect_horz","uknown tracer_HorizontalAdvection_type")
    END SELECT

    stop_timer(timer_adv_horz,2)
     
  END SUBROUTINE advect_horz
  !-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE diffuse_horz( patch_3d,          &
    & trac_old,            &
    & transport_state,                &
    & operators_coefficients,          &
    & k_h,                 &
    & h_old,               &
    & h_new,               &
    & flux_horz,           &
    & horizontally_diffused_tracer)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp)                               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_ocean_transport_state), TARGET      :: transport_state
    TYPE(t_operator_coeff), INTENT(in)     :: operators_coefficients
    REAL(wp), INTENT(in)                   :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                   :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                   :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)                :: flux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout), OPTIONAL      :: horizontally_diffused_tracer(:,:,:)
    !
    !-------------------------------------------------------------------------------    
    start_timer(timer_dif_horz,3)
    
    CALL diffuse_cell_based( patch_3d,          &
      & trac_old,            &
      & transport_state,                &
      & operators_coefficients,          &
      & k_h,                 &
      & h_old,               &
      & h_new,               &
      & flux_horz) !,           &
      ! & horizontally_diffused_tracer)
    
    stop_timer(timer_dif_horz,3)
     
  END SUBROUTINE diffuse_horz
  !-------------------------------------------------------------------------------


  !-----------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE advect_cell_based( patch_3d,          &
    & trac_old,            &
    & transport_state,                &
    & operators_coefficients,          &
    & k_h,                 &
    & h_old,               &
    & h_new,               &
    & div_advflux_horz,    &
    & div_advflux_vert)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp)                             :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_ocean_transport_state), TARGET    :: transport_state
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                 :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                 :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)              :: div_advflux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)              :: div_advflux_vert(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)   
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
      
!       CASE(upwind)
! 
!         CALL upwind_hflux_oce( patch_3d,  &
!         & trac_old,                       &
!         & transport_state%mass_flux_e,         & 
!         & z_adv_flux_h)
! 
!         
!       CASE(central)  
!       
!         CALL central_hflux_oce( patch_3d, &
!         & trac_old,                       &
!         & transport_state%mass_flux_e,         &
!         & z_adv_flux_h)    
                 
      CASE(horz_flux_twisted_vec_recon)
        ! inUse

        CALL flux_corr_transport_cell( patch_3d,  &
          & trac_old,                             &
          & transport_state,                                 &
          & operators_coefficients,               &
          & k_h,                                  &
          & h_old,                                &
          & h_new,                                &
          & z_adv_flux_h,                         &
          & div_advflux_vert)
              
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
    & transport_state,                &
    & operators_coefficients,         &
    & k_h,                 &
    & h_old,               &
    & h_new,               &
    & div_diff_flux_horz)!,           &
    ! & horizontally_diffused_tracer)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp)                             :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_ocean_transport_state), TARGET    :: transport_state
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
        & z_diff_flux_h,&
        & k_h,          &        
        & subset_range = edges_in_domain)
              
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

  !-------------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE flux_corr_transport_cell( patch_3d, &
    & trac_old,                                  & 
    & transport_state,                                      &
    & operators_coefficients,                    &
    & k_h,                                       &
    & h_old,                                     &
    & h_new,                                     &
    & adv_flux_h,                                &
    & div_advflux_vert)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp)                               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_ocean_transport_state), TARGET      :: transport_state
    TYPE(t_operator_coeff), INTENT(in)     :: operators_coefficients
    REAL(wp), INTENT(in)                   :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                   :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                   :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET, INTENT(inout)         :: adv_flux_h(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)!< variable in which the upwind flux is stored
    REAL(wp), INTENT(inout)                :: div_advflux_vert(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)        
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
    
    CHARACTER(len=*), PARAMETER :: method_name = 'flux_corr_transport_cell'
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
        & transport_state%mass_flux_e,         &
        & z_adv_flux_low )
      stop_detail_timer(timer_extra11,5)
            
    CASE DEFAULT
      CALL finish(method_name,"This low-order  option is not supported")   
    END SELECT    


    SELECT CASE(fct_high_order_flux) 
    CASE(horz_flux_twisted_vec_recon)
      !mimetic fluc calculation high order flux
      ! in_use
      CALL map_edges2edges_viacell_3d_const_z( patch_3d, &
        & transport_state%vn,          &
        & operators_coefficients,                &
        & z_adv_flux_high,                       &
        & trac_old)
        
    END SELECT
    !-----------------------------------------------------------------------
    
    !2)call limiter
    SELECT CASE(fct_limiter_horz)
      
    CASE(fct_limiter_horz_zalesak)

      ! inUse
     ! adv_flux_h=z_adv_flux_high
      start_detail_timer(timer_extra13,4)
      CALL limiter_ocean_zalesak_horizontal( patch_3d,   &
        & transport_state%w,           &
        & trac_old,                              &
        & transport_state%mass_flux_e,                &
        & z_adv_flux_low,                        &
        & z_adv_flux_high,                       &
        & adv_flux_h,                            &
        & div_advflux_vert,                      &            
        & operators_coefficients,                &
        & h_old,                                 &
        & h_new)       
      stop_detail_timer(timer_extra13,4)
        
    CASE DEFAULT
     CALL finish('TRIM(flux_corr_transport_h)',"This limiter_type option is not supported")
    END SELECT
    
  END SUBROUTINE flux_corr_transport_cell
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
    
  
END MODULE mo_ocean_tracer_transport_horz
