!>
!! Contains the implementation of the vertical tracer transport routines for the ICON ocean model.
!! This comprises vertical advection.
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/01)
!!
!! @par Copyright and License
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
MODULE mo_ocean_tracer_transport_vert
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  !USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  USE mo_impl_constants,            ONLY: sea_boundary, max_char_length
  USE mo_math_constants,            ONLY: dbl_eps
  USE mo_ocean_nml,                 ONLY: n_zlev,                    &
    & upwind, central,fct_vert_adpo,fct_vert_ppm,fct_vert_zalesak,  &
    & fct_vert_minmod, l_adpo_flowstrength,                         &
    & flux_calculation_vert, fast_performance_level
  USE mo_parallel_config,           ONLY: nproma
  USE mo_run_config,                ONLY: dtime, ltimer
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_adv_vert, timer_ppm_slim, &
    & timer_adpo_vert, timers_level  !, timer_dif_vert,
  USE mo_ocean_types,                 ONLY: t_hydro_ocean_state, t_verticalAdvection_ppm_coefficients, &
    & t_operator_coeff
  USE mo_model_domain,              ONLY: t_patch,t_patch_3d, t_patch_vert
  USE mo_exception,                 ONLY: finish !, message_text, message
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_ocean_physics
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_c, sync_patch_array
  USE mo_ocean_limiter,             ONLY: v_ppm_slimiter_mo, v_ppm_slimiter_mo_onblock 

IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=12)           :: str_module    = 'oceTracVert '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  PUBLIC :: advect_flux_vertical
  PUBLIC :: advect_flux_vertical_high_res
  PUBLIC :: adpo_vtrac_oce
  PUBLIC :: upwind_vflux_oce
 
CONTAINS
  
  !-------------------------------------------------------------------------
  !! SUBROUTINE advects vertically the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
  !! mpi parallelized, sync required: trac_out
!<Optimize:inUse>
  SUBROUTINE advect_flux_vertical( patch_3d,   &
    & trac_old,             &
    & ocean_state,          &
    & operators_coeff,      &
!     & bc_top_tracer,        &
!     & bc_bot_tracer,        &
    & flux_div_vert,        &
    & tracer_id)
    
    TYPE(t_patch_3d ),TARGET :: patch_3d
    REAL(wp), INTENT(inout)           :: trac_old(:,:,:) ! (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_operator_coeff), TARGET    :: operators_coeff
!     REAL(wp)                          :: bc_top_tracer(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!     REAL(wp)                          :: bc_bot_tracer(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)           :: flux_div_vert(:,:,:) ! (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !new tracer
    INTEGER, INTENT(in)               :: tracer_id
    
    !Local variables
    REAL(wp) :: deriv_fst, deriv_sec, adpo_weight_upw_cntr, adpo_weight_cntr_upw
    REAL(wp) :: prism_volume, transport_in, transport_out, adpo_r1, adpo_r2
    INTEGER :: startIndex, endIndex
    INTEGER :: jc, jk, jb
    ! vertical advective tracer fluxes:
    REAL(wp), ALLOCATABLE :: z_adv_flux_v2(:,:,:) ! resulting flux
    REAL(wp), ALLOCATABLE :: z_adv_flux_v (:,:,:)  ! resulting flux
    REAL(wp), ALLOCATABLE :: z_adv_flux_vu(:,:,:)  ! upwind flux
    
    REAL(wp), ALLOCATABLE :: adpo_weight(:,:,:)
    REAL(wp) :: z_flux_div_upw, z_flux_div_cnt
    REAL(wp), ALLOCATABLE :: a_v(:,:,:)
    
    TYPE(t_patch), POINTER :: patch_2D
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_tracer_advection:advect_flux_vertical')
    !-------------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------------
    start_timer(timer_adv_vert,2)

    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    
    ! CALL sync_patch_array(sync_c, patch_2D, trac_old)
    
    ! This is already synced in  edges_in_domain !
    ! CALL sync_patch_array(SYNC_C, patch_2D, ocean_state%p_diag%w_time_weighted)
    

    IF (flux_calculation_vert == fct_vert_ppm) THEN

      ! Vertical advection scheme: piecewise parabolic method (ppm) inUse
      CALL upwind_vflux_ppm( patch_3d,              &
        & trac_old,                                 &
        & ocean_state%p_diag%w_time_weighted,              &
        & dtime, 1,                                &
        & patch_3d%p_patch_1d(1)%prism_thick_c,     &
        & patch_3d%p_patch_1d(1)%inv_prism_thick_c, &
        & operators_coeff%verticalAdvectionPPMcoeffs, &
        & flux_div_vert)

        stop_timer(timer_adv_vert,2)
        RETURN

    ENDIF

    !  The rest is for not upwind_vflux_ppm cases, notInUse
    ALLOCATE(z_adv_flux_v2 (nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks),  & ! resulting flux
              z_adv_flux_v (nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks),  &  ! resulting flux
              z_adv_flux_vu(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks),  &  ! upwind flux
              adpo_weight(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks),      &
              a_v(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks))

    z_adv_flux_v2 (1:nproma, 1:n_zlev+1, 1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp

    z_adv_flux_v (1:nproma, 1:n_zlev+1, 1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    z_adv_flux_vu(1:nproma, 1:n_zlev+1, 1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    adpo_weight  (1:nproma, 1:n_zlev  , 1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp

    SELECT CASE(flux_calculation_vert)

    CASE(upwind)
      !IF (flux_calculation_vert == upwind .OR. flux_calculation_vert == fct_vert_adpo) THEN

      CALL upwind_vflux_oce( patch_3d,                  &
        & trac_old,                    &
        & ocean_state%p_diag%w_time_weighted, &
        & z_adv_flux_v )

      !z_adv_flux_v (1:nproma, 1:n_zlev+1, 1:patch_3D%p_patch_2D(1)%alloc_cell_blocks) = &
      !&  z_adv_flux_vu(1:nproma, 1:n_zlev+1, 1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    CASE(fct_vert_zalesak)
      CALL finish(TRIM(routine), ' This option for vertical tracer transport is not supported')
      !       CALL central_vflux_oce(patch_3D,                 &
      !         &                    trac_old,                   &
      !         &                    ocean_state%p_diag%w_time_weighted,&
      !         &                    z_adv_flux_v)
      !
      !       CALL vflx_limiter_pd_oce( patch_3D,              &
      !         &                       dtime,                   &
      !         &                       trac_old,                &
      !         &                       patch_3D%p_patch_1D(1)%prism_thick_c, &
      !         &                       z_adv_flux_v)


      ! Vertical advection scheme: fct_vert_adpo, adapted from MPIOM (Ernst Meier-Reimer)
    CASE(fct_vert_adpo)

      CALL upwind_vflux_oce( patch_3d,          &
        & trac_old,                             &
        & ocean_state%p_diag%w_time_weighted,          &
        & z_adv_flux_vu )

      CALL adpo_vtrac_oce( patch_3d,            &
        & trac_old,                             &
        & ocean_state%p_diag%w_time_weighted,          &
        & dtime,                                &
        & patch_3d%p_patch_1d(1)%prism_thick_c, &
        & adpo_weight)

      CALL sync_patch_array(sync_c, patch_2D, adpo_weight)

      ! fct_vert_adpo weights as well as flux divergence are defined at cell centers
      ! calculate resulting advective flux divergence using weighting factors for upwind and central
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
        DO jc = startIndex, endIndex
          DO jk = 2,patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

            !z_adv_flux_v(jc,jk,jb)=z_adv_flux_vu(jc,jk,jb)*(1.0_wp-adpo_weight(jc,jk,jb))&
            !&+ z_adv_flux_v2(jc,jk,jb)*adpo_weight(jc,jk,jb)
            !
            !without the weighting factors (i.e. 1-adpo_weight=0=adpo_weigth) the added flux is identical to central flux
            z_adv_flux_v(jc,jk,jb) = z_adv_flux_vu(jc,jk,jb)&
              & +0.5_wp*ABS(ocean_state%p_diag%w_time_weighted(jc,jk,jb))&
              & *(trac_old(jc,jk-1,jb) - trac_old(jc,jk,jb))*adpo_weight(jc,jk,jb)

          ENDDO
        END DO
      END DO

    CASE(fct_vert_minmod)
      CALL advect_flux_vertical_high_res( patch_3d,  &
        & trac_old,           &
        & ocean_state,               &
        & a_v,                &
        & z_adv_flux_v)

    CASE default
      CALL finish(TRIM(routine), ' This option for vertical tracer transport is not supported')
    END SELECT

    ! CALL sync_patch_array(SYNC_C, patch_2D, z_adv_flux_v)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      DO jc = startIndex, endIndex

        DO jk = 1,patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          ! positive vertical divergence in direction of w (upward positive)
          flux_div_vert(jc,jk,jb) = z_adv_flux_v(jc, jk, jb) &
            & - z_adv_flux_v(jc, jk+1, jb)
        ENDDO
      END DO
    END DO

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print thisLevel (1-5, fix)
    CALL dbg_print('AdvVert: adv_flux_v'     ,z_adv_flux_v                ,str_module,idt_src, in_subset=cells_in_domain)
    !---------------------------------------------------------------------
 
    DEALLOCATE(z_adv_flux_v2, & ! resulting flux
              z_adv_flux_v,    &  ! resulting flux
              z_adv_flux_vu,   &  ! upwind flux
              adpo_weight,     &
              a_v)
    
    CALL sync_patch_array(sync_c, patch_2D, flux_div_vert)
    stop_timer(timer_adv_vert,2)
    
  END SUBROUTINE advect_flux_vertical
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !! SUBROUTINE advects vertically the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  SUBROUTINE advect_flux_vertical_high_res( &
    & patch_3d,         &  
    & trac_old,         &
    & ocean_state,      &
    & a_v,              &
    & adv_flux_v)
    
    !TYPE(t_patch), TARGET, INTENT(IN) :: patch_2D
    TYPE(t_patch_3d ),TARGET :: patch_3d
    REAL(wp), INTENT(inout)           :: trac_old(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    REAL(wp)                          :: a_v(:,:,:)
    REAL(wp), INTENT(inout)           :: adv_flux_v(nproma,n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !new tracer
    
    !Local variables
    REAL(wp) :: z_tmp
    INTEGER :: startidx_c, endidx_c
    INTEGER :: jc, jk, jb
    INTEGER :: z_dolic
    !REAL(wp) :: z_adv_flux_v (nproma, n_zlev+1, patch_3D%p_patch_2D(1)%alloc_cell_blocks)  ! vertical advective tracer flux
    REAL(wp) :: z_adv_flux_v2 (nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)  ! vertical advective tracer flux
    REAL(wp) :: z_diff_flux_v(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_mflux(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_consec_grad(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    
    REAL(wp) :: z_limit_phi(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_limit_psi(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !REAL(wp) :: z_adv_flux_high (nproma,n_zlev+1,patch_3D%p_patch_2D(1)%alloc_cell_blocks)  ! vert advective tracer flux
    REAL(wp) :: z_adv_flux_low (nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)  ! vert advective tracer flux
    REAL(wp) :: z_limit_sigma(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_patch_vert), POINTER :: patch_1d
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------------
    patch_2d         => patch_3d%p_patch_2d(1)
    patch_1d         => patch_3d%p_patch_1d(1)
    cells_in_domain  => patch_2d%cells%in_domain
    
    z_limit_psi   (1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    z_limit_sigma (1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    z_mflux       (1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    z_adv_flux_low(1:nproma,1:n_zlev+1,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    adv_flux_v(1:nproma,1:n_zlev+1,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    z_adv_flux_v2(1:nproma,1:n_zlev+1,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    
    CALL sync_patch_array(sync_c, patch_2d, trac_old)
    
    !0) For testing
    CALL upwind_vflux_oce( patch_3d,                  &
      & trac_old,                   &
      & ocean_state%p_diag%w_time_weighted,&
      & z_adv_flux_low )
    
    CALL sync_patch_array(sync_c, patch_2d, z_adv_flux_low)
    CALL sync_patch_array(sync_c, patch_2d, z_adv_flux_v2)
    
    
    !1) calculate phi-part of limiter
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startidx_c, endidx_c)
      
      DO jk = 2, n_zlev
        DO jc = startidx_c, endidx_c
          
          z_diff_flux_v(jc,jk,jb)=trac_old(jc,jk-1,jb)-trac_old(jc,jk,jb)
          
          z_mflux(jc,jk,jb)     = ABS(ocean_state%p_diag%w_time_weighted(jc,jk,jb))
          
          IF (z_mflux(jc,jk,jb) == 0.0_wp) THEN
            z_limit_sigma(jc,jk,jb) = 1.0_wp   !  z_mflux can be zero!
          ELSE
            z_limit_sigma(jc,jk,jb) = MIN(1.0_wp, &
              & 2.0_wp * a_v(jc,jk,jb)/(patch_1d%del_zlev_i(jk)*z_mflux(jc,jk,jb)))
          ENDIF
          
        END DO
      END DO
    END DO
    CALL sync_patch_array(sync_c, patch_2d, z_limit_sigma)
    CALL sync_patch_array(sync_c, patch_2d, z_diff_flux_v)
    CALL sync_patch_array(sync_c, patch_2d, z_mflux)
    
    !This corresponds to (16) in Casulli-Zanolli.
    CALL ratio_consecutive_gradients(patch_3d,ocean_state%p_diag%w_time_weighted, trac_old,z_consec_grad)
    
    !3) calculate psi-part of limiter (see (17)-(19) in Casulli-Zanolli).
    CALL calculate_limiter_function(patch_3d, z_limit_sigma, z_consec_grad,z_limit_phi)
    
    !4) Calculate limited advective flux
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startidx_c, endidx_c)
      DO jk = 2, n_zlev
        DO jc = startidx_c, endidx_c
          
          IF(z_consec_grad(jc,jk,jb)>0)THEN
            z_limit_psi(jc,jk,jb) = z_limit_phi(jc,jk,jb)-z_limit_sigma(jc,jk,jb)
          ELSE
            z_limit_psi(jc,jk,jb) = 0.0_wp
          ENDIF
          !only low gives upwind, low +mflux*diff and without limiter gives central
          adv_flux_v(jc,jk,jb) =           &
            & z_adv_flux_low (jc,jk,jb)        &
            & +0.5_wp*z_mflux(jc,jk,jb)*z_diff_flux_v(jc,jk,jb)*z_limit_psi(jc,jk,jb)
        END DO
      END DO
    END DO
    CALL sync_patch_array(sync_c, patch_2d, z_limit_psi)
    CALL sync_patch_array(sync_c, patch_2d, adv_flux_v)
    
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print thisLevel (1-5, fix)
    CALL dbg_print('AdvDifVert: w_time_weighted',ocean_state%p_diag%w_time_weighted ,str_module,idt_src)
    !CALL dbg_print('AdvDifVert: div_mass_flx_c' ,ocean_state%p_diag%div_mass_flx_c  ,str_module,idt_src)
    CALL dbg_print('AdvVert: adv_flux_v'     ,adv_flux_v                ,str_module,idt_src)
    !CALL dbg_print('AdvVert: flux_div_vert'  ,flux_div_vert                ,str_module,idt_src)
    !---------------------------------------------------------------------
  END SUBROUTINE advect_flux_vertical_high_res
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates ratio of consecutive gradients following Casulli-Zanolli.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2013).
  !!
  SUBROUTINE ratio_consecutive_gradients( patch_3d, w_time_weighted, &
    & trac_old,         &
    & consec_grad)
    
    TYPE(t_patch_3d),TARGET, INTENT(in):: patch_3d
    REAL(wp), INTENT(in)               :: w_time_weighted(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)            :: consec_grad(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !
    !Local variables
    INTEGER :: startIndex, endIndex
    INTEGER :: jc, jk, jb!, je!, ic,ib
    INTEGER :: i_edge, ii_e, ib_e
    REAL(wp) :: z_diff_trac
    REAL(wp) :: z_cellsum_mass_flx_in(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_cellsum_tracdiff_flx_in(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_mflux(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks,1:2)
    INTEGER, DIMENSION(:,:,:), POINTER :: cell_of_edge_idx, cell_of_edge_blk
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_patch_vert), POINTER :: patch_1d
    !-------------------------------------------------------------------------------
    patch_1d        => patch_3d%p_patch_1d(1)
    patch_2d        => patch_3d%p_patch_2d(1)
    ! edges_in_domain => patch_2D%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------------
    consec_grad              (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
    z_cellsum_mass_flx_in    (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
    z_cellsum_tracdiff_flx_in(1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
    z_mflux(1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks,1:2) = 0.0_wp
    
    
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      DO jk = 2, n_zlev-1
        DO jc = startIndex, endIndex
          
          !we are on tracer cell (jc,jk,jb)
          !vertical velocity at top of the cell is w(jc,jk,jb)
          !vertical velocity at bottom of the cell is w(jc,jk+1,jb)
          
          !consider outflow at prism top
          IF( w_time_weighted(jc,jk,jb)>0.0)THEN
            
            IF( w_time_weighted(jc,jk+1,jb)<=0.0)THEN
              z_cellsum_mass_flx_in(jc,jk,jb) =ABS(w_time_weighted(jc,jk+1,jb))
              
              !This is in Casulli Zanolli, eq. (16) the numerator of
              !the consecutive gradient
              z_cellsum_tracdiff_flx_in(jc,jk,jb) &
                & = z_cellsum_mass_flx_in(jc,jk,jb)*(trac_old(jc,jk,jb)-trac_old(jc,jk+1,jb))
              
            ENDIF
            
            IF(z_cellsum_mass_flx_in(jc,jk,jb)/=0.0_wp.AND.(trac_old(jc,jk-1,jb)-trac_old(jc,jk,jb))/=0.0_wp)THEN
              consec_grad(jc,jk,jb)=z_cellsum_tracdiff_flx_in(jc,jk,jb)/&
                & (z_cellsum_mass_flx_in(jc,jk,jb)*(trac_old(jc,jk-1,jb)-trac_old(jc,jk,jb)))
            ELSE
              consec_grad(jc,jk,jb)=0.0_wp
            ENDIF
          END IF
          
          !consider outflow at prism bottom
          IF( w_time_weighted(jc,jk+1,jb)>0.0)THEN
            IF( w_time_weighted(jc,jk,jb)<=0.0)THEN
              z_cellsum_mass_flx_in(jc,jk,jb) =ABS(w_time_weighted(jc,jk,jb))
              
              !This is in Casulli Zanolli, eq. (16) the numerator of
              !the consecutive gradient
              z_cellsum_tracdiff_flx_in(jc,jk,jb) &
                & = z_cellsum_mass_flx_in(jc,jk,jb)*(trac_old(jc,jk,jb)-trac_old(jc,jk+1,jb))
            ENDIF
            IF(z_cellsum_mass_flx_in(jc,jk,jb)/=0.0_wp.AND.(trac_old(jc,jk-1,jb)-trac_old(jc,jk,jb))/=0.0_wp)THEN
              consec_grad(jc,jk,jb)=z_cellsum_tracdiff_flx_in(jc,jk,jb)/&
                & (z_cellsum_mass_flx_in(jc,jk,jb)*(trac_old(jc,jk-1,jb)-trac_old(jc,jk,jb)))
            ELSE
              consec_grad(jc,jk,jb)=0.0_wp
            ENDIF
            
          ENDIF
        END DO
      END DO
    END DO
    CALL sync_patch_array(sync_c, patch_2d, consec_grad)
  END SUBROUTINE ratio_consecutive_gradients
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  !>
  !!  SUBROUTINE calculates ratio of consecutive gradients following Casulli-Zanolli.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2013).
  !!
  SUBROUTINE calculate_limiter_function( patch_3d,   &
    & limit_sigma,                            &
    & consec_grad,                            &
    & limit_phi)
    
    TYPE(t_patch_3d),TARGET, INTENT(in):: patch_3d
    REAL(wp), INTENT(in)               :: limit_sigma(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)               :: consec_grad(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)            :: limit_phi(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !
    !Local variables
    INTEGER :: startIndex, endIndex
    INTEGER :: jk, jb, jc
    
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-------------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    !edges_in_domain => patch_2D%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------------
    limit_phi(1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
    
    !SELECT CASE(limiter_type)
    !CASE(minmod)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      DO jk = 1, n_zlev
        DO jc = startIndex, endIndex
          
          limit_phi(jc,jk,jb)&
            & =MAX(limit_sigma(jc,jk,jb), MIN(1.0_wp,consec_grad(jc,jk,jb)))
        END DO
      END DO
    END DO
    
    !     CASE(superbee)
    !       DO jb = edges_in_domain%start_block, edges_in_domain%end_block
    !         CALL get_index_range(edges_in_domain, jb, startIndex_e, endIndex_e)
    !         DO jk = 1, n_zlev
    !           DO je = startIndex_e, endIndex_e
    !             limit_phi(je,jk,jb) = MAX(limit_sigma(je,jk,jb),&
    !               & MIN(1.0_wp,2.0_wp*consec_grad(je,jk,jb)),&
    !               & MIN(2.0_wp,consec_grad(je,jk,jb)))
    !           END DO
    !         END DO
    !       END DO
    !     CASE(van_leer)
    !       DO jb = edges_in_domain%start_block, edges_in_domain%end_block
    !         CALL get_index_range(edges_in_domain, jb, startIndex_e, endIndex_e)
    !         DO jk = 1, n_zlev
    !           DO je = startIndex_e, endIndex_e
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
    CALL sync_patch_array(sync_c, patch_2d, limit_phi)
    !write(0,*)'limiter-function: max-min',maxval(limit_phi(5,:,5)),minval(limit_phi(5,:,5))
  END SUBROUTINE calculate_limiter_function
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
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
  !! mpi parallelized, no sync
  SUBROUTINE upwind_vflux_oce( patch_3d, pvar_c, pw_c, pupflux_i )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)           :: pvar_c(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)     !< advected cell centered variable
    REAL(wp), INTENT(inout)           :: pw_c(nproma,n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)     !< in: vertical velocity on cells
    REAL(wp), INTENT(inout)           :: pupflux_i(nproma,n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !< variable in which the upwind flux is stored
    ! local variables
    ! height based but reversed (downward increasing depth) coordinate system,
    ! grid coefficient is negative (same as pressure based atmospheric coordinate system
    !REAL(wp), PARAMETER :: zcoeff_grid = -1.0_wp
    !INTEGER             :: z_dolic
    INTEGER :: startIndex, endIndex
    INTEGER :: jc, jk, jb               !< index of cell, vertical thisLevel and block
    
    TYPE(t_patch), POINTER :: patch_2D
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    
!ICON_OMP_PARALLEL_DO PRIVATE(startIndex, endIndex, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      pupflux_i(:,:,jb) = 0.0_wp
      DO jc = startIndex, endIndex
        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          ! calculate vertical tracer flux using upwind method
          pupflux_i(jc,jk,jb) =                 &
            & laxfr_upflux_v( pw_c(jc,jk,jb),  &
            & pvar_c(jc,jk - 1,jb), pvar_c(jc,jk,jb) )
        ENDDO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO
    
  END SUBROUTINE upwind_vflux_oce
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!
  !! Calculation of central vertical tracer fluxes
  !!
  !! @par Revision History
  !! Peter Korn, MPI-M
  !!
  SUBROUTINE central_vflux_oce( patch_3d, pvar_c, pw_c, c_flux_i )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)  :: pvar_c(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)     !< advected cell centered variable
    REAL(wp), INTENT(inout)  :: pw_c(nproma,n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)     !< in: vertical velocity on cells
    REAL(wp), INTENT(inout)  :: c_flux_i(nproma,n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !< variable in which the central flux is stored
    ! local variables
    INTEGER :: startIndex, endIndex
    INTEGER :: jc, jk, jb               !< index of cell, vertical thisLevel and block
    INTEGER :: jkm1                     !< jk - 1
    INTEGER :: z_dolic
    TYPE(t_patch), POINTER :: patch_2D
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      DO jc = startIndex, endIndex
        z_dolic = patch_3d%p_patch_1d(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)
        !IF(z_dolic>=MIN_DOLIC)THEN
        DO jk = 2, z_dolic
          IF ( patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            ! index of top half thisLevel
            jkm1 = jk - 1
            ! calculate vertical tracer flux using central method
            c_flux_i(jc,jk,jb) = 0.5_wp*pw_c(jc,jk,jb)              &
              & * (pvar_c(jc,jkm1,jb)+ pvar_c(jc,jk,jb) )
          ENDIF
        ENDDO
      END DO
    END DO
    
  END SUBROUTINE central_vflux_oce
  
  !-------------------------------------------------------------------------
  !>
  !!
  !! Calculation of vertical tracer advection using ADPO-scheme of MPIOM
  !!
  !! This scheme uses a weighted average of upwind and central tracer advection schemes.
  !! The weight 'r' is dependent on the first and second derivative
  !! of the tracer, the weight 'R' reflects the strength of the flow, according
  !! to the MPIOM documentation
  !!
  !! The code is adopted from the MPIOM-model, written by Ernst Meier-Reimer
  !! The implementation follows the MPIOM Draft documentation, section 5.2.13 describing ocadpo.f90
  !! The adpo weights are located in the center of the cells
  !!
  !! LITERATURE
  !!  - MPIOM Draft Documentation (Wetzel, Haak, Jungclaus, Meier-Reimer)
  !!  - Sweby (1984), SIAM-JNA, 21, 995-1011
  !!
  !! @par Revision History
  !! Initial revision by Stephan Lorenz, MPI-M (2013-10)
  !!
  !! vertical column only, no sync
  !!
  SUBROUTINE adpo_vtrac_oce ( patch_3d,  &
    & p_cc,        & ! advected cell centered variable
    & p_w,         & ! vertical velocity
    & p_dtime,     & ! time step
    & p_thick_c,   & ! layer thickness
    & p_adpo_w)      ! output: weights for upwind and central part
    
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)  :: p_cc        (nproma,n_zlev  , patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)  :: p_w         (nproma,n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)     :: p_dtime
    REAL(wp), INTENT(inout)  :: p_thick_c   (nproma,n_zlev,   patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)  :: p_adpo_w    (nproma,n_zlev,   patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    
    ! local variables
    INTEGER :: jc, jk, jb, startIndex, endIndex
    REAL(wp) :: deriv_fst, deriv_sec
    REAL(wp) :: prism_volume, new_volume, adpo_r1, adpo_r2
    REAL(wp) :: advvol_up, advvol_down, wupw_up, wupw_down
    REAL(wp) :: z_adpo_r1     (nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_adpo_vol_up  (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: z_adpo_vol_down(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_patch),        POINTER :: patch_2D
    REAL(wp),             POINTER :: cell_area(:,:)
    TYPE(t_subset_range), POINTER :: cells_in_domain
    
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    cell_area       => patch_2D%cells%area
    
    z_adpo_r1      (1:nproma, 1:n_zlev+1, 1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    z_adpo_vol_up  (1:nproma, 1:n_zlev,   1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    z_adpo_vol_down(1:nproma, 1:n_zlev,   1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    
    ! special treatment of top
    !  - boundary condition for top and bottom: pure upwind; r=R=0.0
    !  - top: no calculation necessary
    z_adpo_r1(1:nproma,1,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    p_adpo_w (1:nproma,1,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    
    ! special treatment of bottom: r=R=0.0
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      DO jc = startIndex, endIndex
        jk = patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
        IF ( jk > 1 ) THEN  !  sea-point
          z_adpo_r1(jc,jk,jb) = 0.0_wp
          p_adpo_w (jc,jk,jb) = 0.0_wp
        ENDIF
      ENDDO
    ENDDO
    
    ! calculate weighting factors for part of upwind and central schemes: main vertical loop
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      DO jc = startIndex, endIndex
        DO jk = 2,patch_3d%p_patch_1d(1)%dolic_c(jc,jb)-1
          
          ! calculation of first and second spatial derivative
          prism_volume  = p_thick_c(jc,jk,jb) * cell_area(jc,jb)
          
          !upwind biased first derivative
          !IF(p_w(jc,jk,jb)<0)THEN
          !  deriv_fst = ABS(p_cc(jc,jk-1,jb)-p_cc(jc,jk,jb))
          !ELSE
          !  deriv_fst = ABS(p_cc(jc,jk,jb)-p_cc(jc,jk+1,jb))
          !ENDIF
          deriv_fst = ABS(p_cc(jc,jk-1,jb)-p_cc(jc,jk+1,jb))
          deriv_sec = ABS(p_cc(jc,jk-1,jb)+p_cc(jc,jk+1,jb) - 2.0_wp*p_cc(jc,jk,jb))
          ! first weighting factor corresponding to "r" in MPIOM documentation
          !  - second derivative is small (r->1) -> central
          !  - with small first but large second derivative there is an extremum (r->0) -> upwind
          adpo_r1 = MAX(0.0_wp, (deriv_fst-deriv_sec) / (deriv_fst + dbl_eps))
          z_adpo_r1(jc,jk,jb) = adpo_r1  ! for checks
          
          
          IF(l_adpo_flowstrength) THEN
            ! second weighting factor corresponding to "R" in MPIOM documentation
            !  - weak flow or large r -> central; strong flow or small r -> upwind
            !  - this calculation is independent of tracers, see ocadpo_base in MPIOM
            !  - z_adpo_vol_in  is water transport wtp in mpiom, is U_in  in MPIOM documentation
            !  - z_adpo_vol_out is water transport wtm in mpiom, is U_out in MPIOM documentation
            !  - these transports can be replaced by upwind flux z_adv_flux_vu
            
            !transport upward and downward of gridcell (jc,jk,jb)
            !
            !if vertical velocity at cell top jk is positive -> flow goes from cell (jc,jk,jb) to (jc,jk+1,jb)
            !flow moves upward -> wupw_up >0
            !if vertical velocity at cell top jk is negative -> flow goes from cell (jc,jk-1,jb) to (jc,jk,jb)
            !flow moves downward wupw_up =0
            !if vertical velocity at cell bottom jk+1 is positive -> flow goes from cell (jc,jk+1,jb) to (jc,jk,jb)
            !flow moves upward -> wupw_down=0
            !if vertical velocity at cell bottom jk+1 is negative -> flow goes from cell (jc,jk,jb) to (jc,jk+1,jb)
            !flow moves downward -> wupw_down>0
            !
            !this part of the volume moves up
            wupw_up    =     p_w(jc,jk  ,jb)  + ABS(p_w(jc,jk ,jb))
            !this part of the volume moves down
            wupw_down   = ABS(p_w(jc,jk+1,jb)) -  p_w(jc,jk+1,jb)
            
            ! multiplying with factors - ventilation time of the cell
            !            z_adpo_vol_up (jc,jk,jb) = wupw_up   * 0.5_wp * p_dtime * cell_area(jc,jb)
            !            z_adpo_vol_down(jc,jk,jb) = wupw_down * 0.5_wp * p_dtime * cell_area(jc,jb)
            !            p_adpo_w(jc,jk,jb) = MIN(1.0_wp, &
            !              &       adpo_r1*prism_volume / (z_adpo_vol_up(jc,jk,jb) + z_adpo_vol_down(jc,jk,jb) + 1.E-20_wp))! * &
            !              !&       patch_3D%wet_c(jc,jk+1,jb)  ! set zero on land
            
            z_adpo_vol_up (jc,jk,jb)  = wupw_up   * 0.5_wp * p_dtime
            z_adpo_vol_down(jc,jk,jb) = wupw_down * 0.5_wp * p_dtime
            
            p_adpo_w(jc,jk,jb) = MIN(1.0_wp, &
              & adpo_r1*p_thick_c(jc,jk,jb) / (z_adpo_vol_up(jc,jk,jb) + z_adpo_vol_down(jc,jk,jb) + dbl_eps))! * &
            !&       patch_3D%wet_c(jc,jk+1,jb)  ! set zero on land
            
          ELSEIF(.NOT.l_adpo_flowstrength) THEN
            ! only the first condition is active, this is when cpp-key SMOADV in MPIOM exists
            ! r1 is reduced to a maximum of 1; why is wet(jk+1) instead of wet(jk) involved?
            p_adpo_w(jc,jk,jb) = adpo_r1!MIN(adpo_r1,patch_3D%wet_c(jc,jk+1,jb))
          ENDIF
          
        ENDDO
      ENDDO
    ENDDO
    
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print thisLevel (1-5, fix)
    !IF (tracer_id == 1) THEN
    CALL dbg_print('AdvVertAdpo: weight: r1 '    ,z_adpo_r1,     str_module,idt_src, in_subset=cells_in_domain)
    CALL dbg_print('AdvVertAdpo: adpo_weight'    ,p_adpo_w,      str_module,idt_src, in_subset=cells_in_domain)
    CALL dbg_print('AdvVertAdpo: adpo_vol_up'    ,z_adpo_vol_up ,str_module,idt_src, in_subset=cells_in_domain)
    CALL dbg_print('AdvVertAdpo: adpo_vol_down'  ,z_adpo_vol_down,str_module,idt_src, in_subset=cells_in_domain)
    !ENDIF
    !---------------------------------------------------------------------
    
  END SUBROUTINE adpo_vtrac_oce
  !------------------------------------------------------------------------
  
  !------------------------------------------------------------------------
  !! The third order PPM scheme
  !!
  !! Calculation of time averaged vertical tracer fluxes using the third
  !! order PPM scheme.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2009-08-12)
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized to height based vertical coordinate systems. Included
  !!   parameter coeff_grid, in order to apply the same code to either a
  !!   pressure based or height based vertical coordinate system.
  !! Modification by Daniel Reinert, DWD (2011-01-17)
  !! - added optional parameter opt_lout_edge which will provide the
  !!   reconstructed 'edge' value.
  !!
  !! mpi parallelized, sync required: p_upflux
  !
  ! !LITERATURE
  ! - Colella and Woodward (1984), JCP, 54, 174-201
  ! - Carpenter et al. (1989), MWR, 118, 586-612
  ! - Lin and Rood (1996), MWR, 124, 2046-2070
  !
!<Optimize:inUse>
  SUBROUTINE upwind_vflux_ppm( patch_3d, p_cc,        &
    & p_w, p_dtime, p_itype_vlimit,&
    & p_cellhgt_mc_now, cell_invheight, &
    & verticalAdvection_ppm_coefficients, &
    & flux_div_vert)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)           :: p_cc(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< advected cell centered variable
    REAL(wp), INTENT(inout)           :: p_w(nproma,n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !<  in, vertical velocity
    REAL(wp), INTENT(in)              :: p_dtime  !< time step
    REAL(wp), INTENT(inout)           :: p_cellhgt_mc_now(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< layer thickness at cell center at time n
    REAL(wp), INTENT(inout)           :: cell_invheight(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< layer thickness at cell center at time n
    TYPE(t_verticalAdvection_ppm_coefficients), POINTER :: verticalAdvection_ppm_coefficients(:)
    REAL(wp), INTENT(inout)           :: flux_div_vert(:,:,:) ! (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !new tracer
    INTEGER, INTENT(in)               :: p_itype_vlimit                                  !< parameter to select limiter
    !
    !local variables
    REAL(wp) :: p_upflux(nproma,n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)      !<  tracer flux
    REAL(wp) :: z_face(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)   !< face values of transported field
    REAL(wp) :: z_face_up(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< face value (upper face)
    REAL(wp) :: z_face_low(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks) !< face value (lower face)
    REAL(wp) :: z_lext_1(nproma,n_zlev+1)                 !< linear extrapolation value 1
    REAL(wp) :: z_lext_2(nproma,n_zlev+1)                 !< linear extrapolation value 2
    REAL(wp) :: z_cfl_m(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< CFL number (weta>0, w<0)
    REAL(wp) :: z_cfl_p(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< CFL number (weta<0, w>0)
    REAL(wp) :: z_slope(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< monotonized slope
    REAL(wp) :: z_slope_u, z_slope_l                            !< one-sided slopes
    REAL(wp) :: z_delta_m, z_delta_p                            !< difference between lower and upper face value
    !< for weta >0 and weta <0
    REAL(wp) :: z_a11, z_a12                                    !< 1/6 * a6,i (see Colella and Woodward (1984))
    REAL(wp) :: z_weta_dt                                       !< weta times p_dtime
    INTEGER :: firstLevel, secondLevel                                    !< vertical start thisLevel and start thisLevel +1
    INTEGER :: nlevp1                                          !< number of full and half levels
    INTEGER :: ikm1, ikp1, ikp1_ic,level2Below                        !< vertical thisLevel minus and plus one, plus two
    INTEGER :: startIndex, endIndex
    INTEGER :: jc, jk, jb                                      !< index of cell, vertical thisLevel and block
    INTEGER, PARAMETER :: islopel_vsm = 1
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-----------------------------------------------------------------------
    IF (fast_performance_level > 20) THEN
      CALL upwind_vflux_ppm_fast( patch_3d, p_cc,  &
        & p_w, p_dtime, p_itype_vlimit,            &
        & p_cellhgt_mc_now, cell_invheight,        &
        & verticalAdvection_ppm_coefficients,      &
        & flux_div_vert)
      RETURN
    ENDIF
    !-----------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    
    firstLevel   = 1
    secondLevel = 2
    nlevp1 = n_zlev+1
    
    z_cfl_m   (1:nproma,1:n_zlev+1,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    z_cfl_p   (1:nproma,1:n_zlev+1,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    z_face    (1:nproma,1:n_zlev+1,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    z_face_low(1:nproma,1:n_zlev,  1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    z_face_up (1:nproma,1:n_zlev,  1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    z_slope   (1:nproma,1:n_zlev,  1:patch_3d%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    z_lext_1  (1:nproma,1:n_zlev+1)=0.0_wp
    z_lext_2  (1:nproma,1:n_zlev+1)=0.0_wp
    
    ! this is not needed
    ! CALL sync_patch_array(SYNC_C, patch_2D, p_cc)
    ! CALL sync_patch_array(SYNC_C, patch_2D, p_cellhgt_mc_now)
    ! CALL sync_patch_array(SYNC_C, patch_2D, p_w)
    
    ! advection is done with an upwind scheme and a piecwise parabolic
    ! approx. of the subgrid distribution is used.
    ! 3 options:  standard without limiter
    !             standard with semi-monotone or monotone limiter
    !             special version with limiter which handles CFL >1
    !
    ! 1. Calculate Courant number for weta>0 (w<0) and weta<0 (w>0)
    !
    ! !$OMP PARALLEL
    ! !$OMP DO PRIVATE(jb,jk,jc,startIndex,endIndex,ikm1,z_weta_dt,ikp1_ic,ikp1, &
    ! !$OMP            z_slope_u,z_slope_l,level2Below)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      
      ! Courant number at top and bottom
      z_cfl_p(startIndex:endIndex,firstLevel,jb)   = 0._wp
      z_cfl_m(startIndex:endIndex,firstLevel,jb)   = 0._wp
      z_cfl_p(startIndex:endIndex,nlevp1,jb) = 0._wp
      z_cfl_m(startIndex:endIndex,nlevp1,jb) = 0._wp
      
      DO jk = secondLevel, n_zlev
        ! index of top half thisLevel
        ikm1 = jk - 1
        DO jc = startIndex, endIndex
          ! Calculate local Courant number at half levels
          ! z_cfl_m for weta >0 (w <0)
          ! z_cfl_p for weta <0 (w >0)
          !z_weta_dt = 0.0_wp
          
          IF ( patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            z_weta_dt = ABS(p_w(jc,jk,jb)) * p_dtime
            
            z_cfl_m(jc,jk,jb) = z_weta_dt &
              & / p_cellhgt_mc_now(jc,ikm1,jb) !patch_3D%p_patch_1D(1)%del_zlev_m(ikm1)!
            z_cfl_p(jc,jk,jb) = z_weta_dt &
              & / p_cellhgt_mc_now(jc,jk,jb)!patch_3D%p_patch_1D(1)%del_zlev_m(jk)!
          ENDIF
        END DO ! end loop over cells
      ENDDO ! end loop over vertical levels
      !
      ! 2. Calculate monotonized slope
      !
      z_slope(startIndex:endIndex,firstLevel,jb) = 0._wp
      
      DO jk = secondLevel, n_zlev
        ! index of top half thisLevel
        ikm1    = jk - 1
        ! index of bottom half thisLevel
        ikp1_ic = jk + 1
        ikp1    = MIN( ikp1_ic, n_zlev )
        
        DO jc = startIndex, endIndex
          IF ( patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            z_slope_u = 2._wp * (p_cc(jc,jk,jb) - p_cc(jc,ikm1,jb))
          ELSE
            z_slope_u = 0.0_wp
          ENDIF
          IF ( patch_3d%lsm_c(jc,ikp1,jb) <= sea_boundary ) THEN
            z_slope_l = 2._wp * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))
          ELSE
            z_slope_l = 0.0_wp
          ENDIF
          
          IF ((z_slope_u * z_slope_l) > 0._wp) THEN
            
            z_slope(jc,jk,jb) = ( p_cellhgt_mc_now(jc,jk,jb)                            &
              & / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)            &
              & + p_cellhgt_mc_now(jc,ikp1,jb)) )                                       &
              & * ( (2._wp * p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)) &
              & / (p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb))           &
              & * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))                                   &
              & + (p_cellhgt_mc_now(jc,jk,jb) + 2._wp * p_cellhgt_mc_now(jc,ikp1,jb))   &
              & / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb))           &
              & * (p_cc(jc,jk,jb) - p_cc(jc,ikm1,jb)) )
            
            z_slope(jc,jk,jb) = SIGN(                                            &
              & MIN( ABS(z_slope(jc,jk,jb)), ABS(z_slope_u), ABS(z_slope_l) ),  &
              & z_slope(jc,jk,jb))
          ELSE
            z_slope(jc,jk,jb) = 0._wp
          ENDIF
        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      
      !
      ! 3. reconstruct face values at vertical half-levels
      !
      ! Boundary values for two highest and lowest half-levels
      !
      ! for faces k=secondLevel and k=nlevp1-1 reconstructed face values are calculated by
      ! interpolating a quadratic (instead of quartic) polynomial through 3
      ! values of the indefinite integral A=\int_{\eta_{0}}^{\eta}q\,\mathrm{d}\eta
      !
      ! for faces k=firstLevel and k=nlevp1 a zero gradient condition is assumed and the
      ! face values are set to the tracer values of the corresponding cell centers
      !
      DO jc = startIndex, endIndex
        IF ( patch_3d%lsm_c(jc,secondLevel,jb) <= sea_boundary ) THEN
          z_face(jc,secondLevel,jb) = p_cc(jc,firstLevel,jb)*(1._wp - (p_cellhgt_mc_now(jc,firstLevel,jb)&
            & / p_cellhgt_mc_now(jc,secondLevel,jb))) + (p_cellhgt_mc_now(jc,firstLevel,jb)  &
            & /(p_cellhgt_mc_now(jc,firstLevel,jb) + p_cellhgt_mc_now(jc,secondLevel,jb)))   &
            & * ((p_cellhgt_mc_now(jc,firstLevel,jb) / p_cellhgt_mc_now(jc,secondLevel,jb))  &
            & * p_cc(jc,firstLevel,jb) + p_cc(jc,secondLevel,jb))
        ELSE
          z_face(jc,secondLevel,jb) = 0.0_wp
        ENDIF
        IF ( patch_3d%lsm_c(jc,n_zlev,jb) <= sea_boundary ) THEN
          z_face(jc,n_zlev,jb) = p_cc(jc,n_zlev-1,jb)*( 1._wp                               &
            & - (p_cellhgt_mc_now(jc,n_zlev-1,jb) / p_cellhgt_mc_now(jc,n_zlev,jb)))  &
            & + (p_cellhgt_mc_now(jc,n_zlev-1,jb)/(p_cellhgt_mc_now(jc,n_zlev-1,jb)   &
            & + p_cellhgt_mc_now(jc,n_zlev,jb))) * ((p_cellhgt_mc_now(jc,n_zlev-1,jb) &
            & / p_cellhgt_mc_now(jc,n_zlev,jb)) * p_cc(jc,n_zlev-1,jb)                &
            & + p_cc(jc,n_zlev,jb))
        ELSE
          z_face(jc,n_zlev,jb) = 0.0_wp
        ENDIF
        IF ( patch_3d%lsm_c(jc,firstLevel,jb) <= sea_boundary ) THEN
          z_face(jc,firstLevel,jb) = p_cc(jc,firstLevel,jb)
        ELSE
          z_face(jc,firstLevel,jb) = 0.0_wp
        ENDIF
        IF ( patch_3d%lsm_c(jc,n_zlev,jb) <= sea_boundary ) THEN
          z_face(jc,nlevp1,jb) = p_cc(jc,n_zlev,jb)
        ELSE
          z_face(jc,nlevp1,jb) = 0.0_wp
        ENDIF
      ENDDO
      
      
      DO jk = secondLevel, n_zlev-2
        ! index of top half thisLevel
        ikm1 = jk - 1
        ! index of bottom half thisLevel
        ikp1 = jk + 1
        level2Below = jk + 2
        DO jc = startIndex, endIndex
          IF ( patch_3d%lsm_c(jc,level2Below,jb) <= sea_boundary ) THEN
            z_face(jc,ikp1,jb) = p_cc(jc,jk,jb) &
              & + (p_cellhgt_mc_now(jc,jk,jb)                                          &
              & / (p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb)))         &
              & * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))                                  &
              & + (1._wp/(p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)    &
              & + p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,level2Below,jb)))        &
              & * ( (2._wp * p_cellhgt_mc_now(jc,ikp1,jb) * p_cellhgt_mc_now(jc,jk,jb) &
              & / (p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb)))         &
              & * ( (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb))        &
              & / (2._wp*p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb))    &
              & - (p_cellhgt_mc_now(jc,level2Below,jb) + p_cellhgt_mc_now(jc,ikp1,jb))        &
              & / (2._wp*p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb)) )  &
              & * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb)) - p_cellhgt_mc_now(jc,jk,jb)     &
              & * z_slope(jc,ikp1,jb) * (p_cellhgt_mc_now(jc,ikm1,jb)                  &
              & + p_cellhgt_mc_now(jc,jk,jb)) / (2._wp*p_cellhgt_mc_now(jc,jk,jb)      &
              & + p_cellhgt_mc_now(jc,ikp1,jb)) + p_cellhgt_mc_now(jc,ikp1,jb)         &
              & * z_slope(jc,jk,jb) * (p_cellhgt_mc_now(jc,ikp1,jb)                    &
              & + p_cellhgt_mc_now(jc,level2Below,jb)) / (p_cellhgt_mc_now(jc,jk,jb)          &
              & + 2._wp*p_cellhgt_mc_now(jc,ikp1,jb)) )
          ELSE
            z_face(jc,ikp1,jb) = 0.0_wp
          ENDIF
        END DO ! end loop over cells
      END DO ! end loop over vertical levels
      
    END DO
    ! !$OMP END DO
    ! !$OMP END PARALLEL
    
    CALL sync_patch_array(sync_c, patch_2D, z_face)
    CALL sync_patch_array(sync_c, patch_2D, z_face_up)
    CALL sync_patch_array(sync_c, patch_2D, z_face_low)
    !CALL sync_patch_array(SYNC_C, patch_2D, z_slope)
    
    !
    ! 4. Limitation of first guess parabola (which is based on z_face)
    ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
    ! the limitation procedure.
    ! Therefore 2 additional fields z_face_up and z_face_low are
    ! introduced.
    !
    IF (p_itype_vlimit == islopel_vsm) THEN
      !     ! monotonic (mo) limiter
      
      CALL v_ppm_slimiter_mo( patch_3d, &
        & p_cc,       &
        & z_face,     &
        & z_slope,    &
        & z_face_up,  &
        & z_face_low)
      
    ELSE
      !      ! simply copy face values to 'face_up' and 'face_low' arrays
      ! !$OMP PARALLEL
      ! !$OMP DO PRIVATE(jk,ikp1,jb,startIndex,endIndex)
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
        DO jc = startIndex, endIndex
          DO jk = secondLevel, n_zlev-1
            ! index of bottom half thisLevel
            ikp1 = jk + 1
            IF ( patch_3d%lsm_c(jc,ikp1,jb) <= sea_boundary ) THEN
              z_face_up(startIndex:endIndex,jk,jb)  = z_face(startIndex:endIndex,jk,jb)
              z_face_low(startIndex:endIndex,jk,jb) = z_face(startIndex:endIndex,ikp1,jb)
            ELSE
              z_face_up(startIndex:endIndex,jk,jb)  = 0.0_wp
              z_face_low(startIndex:endIndex,jk,jb) = 0.0_wp
            ENDIF
          ENDDO
        END DO
      ENDDO
      ! !$OMP ENDDO
      ! !$OMP END PARALLEL
    ENDIF  !  p_ityp_vlimit
    
    
    ! !$OMP PARALLEL
    ! 5b. extrapolation using piecewise parabolic approx. of the transported
    ! quantity to the edge and finally, calculation of the upwind fluxes
    !
    
    ! !$OMP DO PRIVATE(jb,jk,jc,startIndex,endIndex,z_lext_1,z_lext_2,ikm1,z_delta_m, &
    ! !$OMP            z_delta_p,z_a11,z_a12)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      
      z_lext_1(startIndex:endIndex,firstLevel)   = p_cc(startIndex:endIndex,firstLevel,jb)
      z_lext_2(startIndex:endIndex,firstLevel)   = p_cc(startIndex:endIndex,firstLevel,jb)
      z_lext_1(startIndex:endIndex,nlevp1) = p_cc(startIndex:endIndex,n_zlev,jb)
      z_lext_2(startIndex:endIndex,nlevp1) = p_cc(startIndex:endIndex,n_zlev,jb)
      p_upflux(:,:,jb) = 0.0_wp
      
      DO jk = secondLevel, n_zlev
        ! index of top half thisLevel
        ikm1 = jk -1
        DO jc = startIndex, endIndex
          IF ( patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            ! linear extrapolated values
            ! for the height based coordinate system multiplication by coeff_grid
            ! is not necessary due to compensating (-) signs.
            ! first (of cell above) (case of w < 0; weta > 0)
            z_delta_m = z_face_low(jc,ikm1,jb) - z_face_up(jc,ikm1,jb)
            z_a11     = p_cc(jc,ikm1,jb)                                  &
              & - 0.5_wp * (z_face_low(jc,ikm1,jb) + z_face_up(jc,ikm1,jb))
            
            z_lext_1(jc,jk) = p_cc(jc,ikm1,jb)                            &
              & + (0.5_wp * z_delta_m * (1._wp - z_cfl_m(jc,jk,jb)))     &
              & - z_a11*(1._wp - 3._wp*z_cfl_m(jc,jk,jb)                 &
              & + 2._wp*z_cfl_m(jc,jk,jb)*z_cfl_m(jc,jk,jb))
            
            ! second (of cell below) (case of w > 0; weta < 0)
            z_delta_p = z_face_low(jc,jk,jb) - z_face_up(jc,jk,jb)
            z_a12     = p_cc(jc,jk,jb)                                    &
              & - 0.5_wp * (z_face_low(jc,jk,jb) + z_face_up(jc,jk,jb))
            
            z_lext_2(jc,jk) = p_cc(jc,jk,jb)                              &
              & - (0.5_wp * z_delta_p * (1._wp - z_cfl_p(jc,jk,jb)))     &
              & - z_a12*(1._wp - 3._wp*z_cfl_p(jc,jk,jb)                 &
              & + 2._wp*z_cfl_p(jc,jk,jb)*z_cfl_p(jc,jk,jb))
            !
            ! calculate vertical tracer flux
            !
            p_upflux(jc,jk,jb) =                                  &
              & laxfr_upflux_v( p_w(jc,jk,jb),       &
              & z_lext_1(jc,jk), z_lext_2(jc,jk))
            
            
          ENDIF
        END DO ! end loop over cells
      ENDDO ! end loop over vertical levels
      !
      ! set lower boundary condition
      !
      p_upflux(startIndex:endIndex,nlevp1,jb) = 0.0_wp
      !
    ENDDO ! end loop over blocks
    ! !$OMP END DO
    ! !$OMP END PARALLEL
    
    ! jc=5
    ! jb=5
    ! Do jk=1,patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
    ! write(0,*)'INSIDE PPM:',jc,jk,jb,p_upflux(jc,jk,jb)
    ! END DO
    
    ! CALL sync_patch_array(sync_c, patch_2D, p_upflux)
    
    
    ! 6. If desired, apply a flux limiter to limit computed fluxes.
    !    These flux limiters are based on work by Zalesak (1979)
    !    positive-definite (pd) limiter      !
    !    IF (p_itype_vlimit == ifluxl_vpd) THEN
    !     IF (l_vert_limiter_advection) &
    !CALL vflx_limiter_pd_oce( patch_3D, p_dtime, p_cc, p_cellhgt_mc_now, p_upflux)
    !   ENDIF
    
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      DO jc = startIndex, endIndex

        DO jk = 1,patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          ! positive vertical divergence in direction of w (upward positive)
          flux_div_vert(jc,jk,jb) = p_upflux(jc, jk, jb) &
            & - p_upflux(jc, jk+1, jb)
        ENDDO
      END DO
    END DO

  END SUBROUTINE upwind_vflux_ppm
  !-------------------------------------------------------------------------
  
  !------------------------------------------------------------------------
  !! Otpimized version of the third order PPM scheme
  !!
!<Optimize:inUse>
  SUBROUTINE upwind_vflux_ppm_fast(       &
    & patch_3d, tracer,                   &
    & w, dtime, vertical_limiter_type,    &
    & cell_thickeness,  cell_invheight,   &
    & verticalAdvection_ppm_coefficients, &
    & flux_div_vert)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d   
    REAL(wp), INTENT(inout)           :: tracer(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< in: advected cell centered variable
    REAL(wp), INTENT(inout)           :: w(nproma,n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !<  in : vertical velocity
    REAL(wp), INTENT(in)              :: dtime  !< time step
    REAL(wp), INTENT(inout)           :: cell_thickeness(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< in: layer thickness at cell center at time n
    REAL(wp), INTENT(inout)           :: cell_invheight(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< layer thickness at cell center at time n
    TYPE(t_verticalAdvection_ppm_coefficients) :: verticalAdvection_ppm_coefficients(patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)           :: flux_div_vert(:,:,:) ! (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !new tracer
    INTEGER, INTENT(in)               :: vertical_limiter_type                                  !< parameter to select limiter
    !
    !-----------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    INTEGER                       :: startIndex, endIndex, jb
    !-----------------------------------------------------------------------
    cells_in_domain => patch_3d%p_patch_2d(1)%cells%in_domain
#ifdef NAGFOR
    flux_div_vert(:,:,:) = 0.0_wp
#endif
   
!ICON_OMP_PARALLEL_DO PRIVATE(startIndex, endIndex) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      CALL  upwind_vflux_ppm_onBlock(         &
            & tracer(:,:,jb),                 &
            & w(:,:,jb),                      &
            & dtime, vertical_limiter_type,   &
            & cell_thickeness(:,:,jb),        &
            & cell_invheight(:,:,jb),         &
            & verticalAdvection_ppm_coefficients(jb), &
            & flux_div_vert(:,:,jb),          & ! out
            & startIndex, endIndex,           &
            & patch_3d%p_patch_1d(1)%dolic_c(:,jb))
    ENDDO ! end loop over blocks
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE upwind_vflux_ppm_fast
  !-------------------------------------------------------------------------

  !------------------------------------------------------------------------
  !! Otpimized version of the third order PPM scheme
  !!
  !! Calculation of time averaged vertical tracer fluxes using the third
  !! order PPM scheme.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2009-08-12)
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized to height based vertical coordinate systems. Included
  !!   parameter coeff_grid, in order to apply the same code to either a
  !!   pressure based or height based vertical coordinate system.
  !! Modification by Daniel Reinert, DWD (2011-01-17)
  !! - added optional parameter opt_lout_edge which will provide the
  !!   reconstructed 'edge' value.
  !!
  !! mpi parallelized, sync required: p_upflux
  !
  ! !LITERATURE
  ! - Colella and Woodward (1984), JCP, 54, 174-201
  ! - Carpenter et al. (1989), MWR, 118, 586-612
  ! - Lin and Rood (1996), MWR, 124, 2046-2070
  !
  ! Optimized version of upwind_vflux_ppm
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE upwind_vflux_ppm_onBlock(  &
    & tracer,                            &
    & w, dtime, vertical_limiter_type,   &
    & cell_thickeness, cell_invheight,   &
    & ppmCoeffs,                         &
    & flux_div_vert,                     &
    & startIndex, endIndex, cells_noOfLevels)

    REAL(wp), INTENT(in)           :: tracer(nproma,n_zlev)      !< advected cell centered variable
    REAL(wp), INTENT(in)           :: w(nproma,n_zlev+1)         !<  in : vertical velocity
    REAL(wp), INTENT(in)           :: dtime                      !< time step
    REAL(wp), INTENT(in)           :: cell_thickeness(nproma,n_zlev) !< layer thickness at cell center at time n
    REAL(wp), INTENT(in)           :: cell_invheight(nproma,n_zlev)!< layer thickness at cell center at time n
    TYPE(t_verticalAdvection_ppm_coefficients) :: ppmCoeffs
    INTEGER, INTENT(in)            :: vertical_limiter_type                    !< parameter to select limiter
    REAL(wp), INTENT(inout)        :: flux_div_vert(nproma, n_zlev) !new tracer
    INTEGER, INTENT(in)            :: startIndex, endIndex
    INTEGER, INTENT(in)            :: cells_noOfLevels(nproma)
    !
    REAL(wp) :: upward_tracer_flux(nproma,n_zlev+1)      !< tracer flux
    REAL(wp) :: z_face(nproma,n_zlev+1)   !< face values of transported field
    REAL(wp) :: z_face_up(nproma,n_zlev)  !< face value (upper face)
    REAL(wp) :: z_face_low(nproma,n_zlev) !< face value (lower face)
    REAL(wp) :: z_lext_1(n_zlev+1)                 !< linear extrapolation value 1
    REAL(wp) :: z_lext_2(n_zlev+1)                 !< linear extrapolation value 2
    REAL(wp) :: z_cfl_m, z_cfl_p !< CFL number (weta>0, w<0), CFL number (weta<0, w>0)
    REAL(wp) :: z_slope(nproma,n_zlev+1)  !< monotonized slope
    REAL(wp) :: z_slope_u, z_slope_l                            !< one-sided slopes
    !< for weta >0 and weta <0
    REAL(wp) :: z_a11, z_a12                                    !< 1/6 * a6,i (see Colella and Woodward (1984))
    REAL(wp) :: z_weta_dt                                       !< weta times dtime
    INTEGER :: firstLevel, secondLevel                                    !< vertical start thisLevel and start thisLevel +1
    INTEGER :: levelAbove, levelBelow,level2Below                        !< vertical thisLevel minus and plus one, plus two
    INTEGER :: jc, thisLevel, cell_levels
    !LOGICAL  :: opt_lout_edge !< optional: output edge value (.TRUE.),
    !                          !< or the flux across the edge   !< (.FALSE./not specified)
    !REAL(wp) :: opt_topflx_tra(nproma,patch_3D%p_patch_2D(1)%alloc_cell_blocks)  !< vertical tracer flux at upper boundary
    INTEGER, PARAMETER :: islopel_vsm = 1
    
    REAL(wp), POINTER ::  cellHeightRatio_This_toBelow(:,:)
    REAL(wp), POINTER ::  cellHeightRatio_This_toThisBelow(:,:)
    REAL(wp), POINTER ::  cellHeight_2xBelow_x_RatioThis_toThisBelow(:,:)
    REAL(wp), POINTER ::  cellHeightRatio_This_toThisAboveBelow(:,:)
    REAL(wp), POINTER ::  cellHeightRatio_2xAboveplusThis_toThisBelow(:,:)
    REAL(wp), POINTER ::  cellHeightRatio_2xBelowplusThis_toThisAbove(:,:)
    REAL(wp), POINTER ::  cellHeightRatio_ThisAbove_to2xThisplusBelow(:,:)
    REAL(wp), POINTER ::  cellHeightRatio_ThisBelow_to2xThisplusAbove(:,:)
    REAL(wp), POINTER ::  cellHeight_inv_ThisAboveBelow2Below(:,:)
    !-----------------------------------------------------------------------
    cellHeightRatio_This_toBelow                   =>  ppmCoeffs%cellHeightRatio_This_toBelow
    cellHeightRatio_This_toThisBelow               =>  ppmCoeffs%cellHeightRatio_This_toThisBelow
    cellHeight_2xBelow_x_RatioThis_toThisBelow     =>  ppmCoeffs%cellHeight_2xBelow_x_RatioThis_toThisBelow
    cellHeightRatio_This_toThisAboveBelow          =>  ppmCoeffs%cellHeightRatio_This_toThisAboveBelow
    cellHeightRatio_2xAboveplusThis_toThisBelow    =>  ppmCoeffs%cellHeightRatio_2xAboveplusThis_toThisBelow
    cellHeightRatio_2xBelowplusThis_toThisAbove    =>  ppmCoeffs%cellHeightRatio_2xBelowplusThis_toThisAbove
    cellHeightRatio_ThisAbove_to2xThisplusBelow    =>  ppmCoeffs%cellHeightRatio_ThisAbove_to2xThisplusBelow
    cellHeightRatio_ThisBelow_to2xThisplusAbove    =>  ppmCoeffs%cellHeightRatio_ThisBelow_to2xThisplusAbove
    cellHeight_inv_ThisAboveBelow2Below            =>  ppmCoeffs%cellHeight_inv_ThisAboveBelow2Below

    firstLevel  = 1
    secondLevel = 2

    ! advection is done with an upwind scheme and a piecwise parabolic
    ! approx. of the subgrid distribution is used.
    ! 3 options:  standard without limiter
    !             standard with semi-monotone or monotone limiter
    !             special version with limiter which handles CFL >1
    !
    !------------------------------------------------
    ! 1. Calculate Courant number for weta>0 (w<0) and weta<0 (w>0)
    ! ..............

    ! 2. Calculate monotonized slope
    !
    z_slope(:, :) = 0._wp

! !CDIR NODEP
    DO jc = startIndex, endIndex

! !CDIR NODEP
      DO thisLevel = secondLevel, cells_noOfLevels(jc) - 1

        levelAbove    = thisLevel - 1                 ! index of top half thisLevel
        levelBelow    = thisLevel + 1  ! index of bottom half thisLevel

        z_slope_u = 2._wp * (tracer(jc,thisLevel)  - tracer(jc,levelAbove))
        z_slope_l = 2._wp * (tracer(jc,levelBelow) - tracer(jc,thisLevel))

        IF ((z_slope_u * z_slope_l) > 0._wp) THEN

          z_slope(jc,thisLevel) = &
            & ( cellHeightRatio_This_toThisAboveBelow(jc,thisLevel)  )           &
            & * ( &
            & (cellHeightRatio_2xAboveplusThis_toThisBelow(jc,thisLevel))        &
            & &
            & * (tracer(jc,levelBelow) - tracer(jc,thisLevel))                   &
            & &
            & + (cellHeightRatio_2xBelowplusThis_toThisAbove(jc,thisLevel))      &
            & &
            & * (tracer(jc,thisLevel) - tracer(jc,levelAbove)) )

          z_slope(jc,thisLevel) = SIGN(                                           &
            & MIN( ABS(z_slope(jc,thisLevel)), ABS(z_slope_u), ABS(z_slope_l) ),  &
            & z_slope(jc,thisLevel))

        ENDIF

      END DO ! jc = startIndex, endIndex
    END DO ! end loop over vertical levels

    !
    ! 3. reconstruct face values at vertical half-levels
    !
    ! Boundary values for two highest and lowest half-levels
    !
    ! for faces k=secondLevel and k=nlevp1-1 reconstructed face values are calculated by
    ! interpolating a quadratic (instead of quartic) polynomial through 3
    ! values of the indefinite integral A=\int_{\eta_{0}}^{\eta}q\,\mathrm{d}\eta
    !
    ! for faces k=firstLevel and k=nlevp1 a zero gradient condition is assumed and the
    ! face values are set to the tracer values of the corresponding cell centers
    !
    z_face(:, :)  = 0.0_wp

! !CDIR NODEP
    DO jc = startIndex, endIndex
      cell_levels = cells_noOfLevels(jc)

! !CDIR NODEP
      DO thisLevel = secondLevel, cell_levels - 2
        levelAbove  = thisLevel - 1
        levelBelow  = thisLevel + 1
        level2Below = thisLevel + 2

        z_face(jc,levelBelow) = tracer(jc,thisLevel) &
          & + (cellHeightRatio_This_toThisBelow(jc,thisLevel))                   &
          & &
          & * (tracer(jc,levelBelow) - tracer(jc,thisLevel))                     &
          & &
          & + cellHeight_inv_ThisAboveBelow2Below(jc,thisLevel)                  &
          & * &
          & ( &
          &   (cellHeight_2xBelow_x_RatioThis_toThisBelow(jc,thisLevel))         &
          & &
          &  * (cellHeightRatio_ThisAbove_to2xThisplusBelow(jc,thisLevel)        &
          &     -  cellHeightRatio_ThisBelow_to2xThisplusAbove(jc,thisLevel))    &
          & &
          &  * (tracer(jc,levelBelow) - tracer(jc,thisLevel))                    &
          &  - z_slope(jc,levelBelow)                                            &
          &  * cell_thickeness(jc,thisLevel)                                     &
          &  * cellHeightRatio_ThisAbove_to2xThisplusBelow(jc,thisLevel)         &
          & &
          &  +  z_slope(jc,thisLevel) *                                          &
          & &
          &  cell_thickeness(jc,levelBelow) * &
          &  cellHeightRatio_ThisBelow_to2xThisplusAbove(jc,levelBelow)          &
          & )

      END DO ! end loop over vertical levels

      ! compute top 2 levels
      IF ( cells_noOfLevels(jc) >= firstLevel  ) THEN
      
        z_face(jc,firstLevel) = tracer(jc,firstLevel)
        
        IF ( cell_levels >= secondLevel ) THEN

          z_face(jc,secondLevel) = &
            & tracer(jc,firstLevel) *                                        &
            & (1._wp - cellHeightRatio_This_toBelow(jc,firstLevel))          &
            & &
            & + (cellHeightRatio_This_toThisBelow(jc,firstLevel))            &
            & * &
            & ( cellHeightRatio_This_toBelow(jc,firstLevel)                  &
            & * tracer(jc,firstLevel)                                        & 
            & + tracer(jc,secondLevel))
            
        ENDIF
      ENDIF

      ! compute bottom thisLevel
      IF ( cells_noOfLevels(jc) > secondLevel ) THEN

        z_face(jc, cell_levels) =                                                       &
          & tracer(jc,cell_levels-1) *                                                  &
          & ( 1._wp - cellHeightRatio_This_toBelow(jc,cell_levels-1))                   &
          & +  &
          & (cell_thickeness(jc,cell_levels-1) / (cell_thickeness(jc,cell_levels-1)     &
          &   + cell_thickeness(jc,cell_levels)))                                       &
          &  * &
          &  (cellHeightRatio_This_toBelow(jc,cell_levels-1)                            &
          &  * tracer(jc,cell_levels-1)                                                 &
          &  + tracer(jc,cell_levels))

      ENDIF

    ENDDO

    ! 4. Limitation of first guess parabola (which is based on z_face)
    ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
    ! the limitation procedure.
    ! Therefore 2 additional fields z_face_up and z_face_low are
    ! introduced.
    z_face_low(1:nproma,1:n_zlev) = 0.0_wp
    z_face_up (1:nproma,1:n_zlev) = 0.0_wp

    IF (vertical_limiter_type == islopel_vsm) THEN
      !     ! monotonic (mo) limiter
      CALL v_ppm_slimiter_mo_onBlock( &
        & tracer(:,:),       &
        & z_face(:,:),       &
        & z_slope(:,:),      &
        & z_face_up,         &
        & z_face_low,        &
        & startIndex,        &
        & endIndex,          &
        & cells_noOfLevels)

    ELSE
        ! simply copy face values to 'face_up' and 'face_low' arrays
        DO jc = startIndex, endIndex
          DO thisLevel = secondLevel, cells_noOfLevels(jc)-1
            z_face_up(jc, thisLevel)  = z_face(jc, thisLevel    )
            z_face_low(jc,thisLevel)  = z_face(jc, thisLevel + 1)
          ENDDO
        END DO

    ENDIF  !  p_ityp_vlimit

    upward_tracer_flux(:,:) = 0.0_wp

! !CDIR NODEP
    DO jc = startIndex, endIndex
! !CDIR NODEP
      DO thisLevel = secondLevel, cells_noOfLevels(jc)
        ! index of top half thisLevel
        levelAbove = thisLevel - 1
        ! linear extrapolated values
        ! for the height based coordinate system multiplication by coeff_grid
        ! is not necessary due to compensating (-) signs.
        ! first (of cell above) (case of w < 0; weta > 0)
        ! z_delta_m = (z_face_low(jc,levelAbove) - z_face_up(jc,levelAbove))
        z_a11     = tracer(jc,levelAbove)                                  &
          & - 0.5_wp * (z_face_low(jc,levelAbove) + z_face_up(jc,levelAbove))

        ! Calculate local Courant number at half levels
        ! z_cfl_m for weta >0 (w <0)
        ! z_cfl_p for weta <0 (w >0)
        ! z_weta_dt = 0.0_wp
        z_weta_dt = ABS(w(jc,thisLevel)) * dtime
        z_cfl_p = z_weta_dt * cell_invheight(jc, thisLevel)
        z_cfl_m = z_weta_dt * cell_invheight(jc, levelAbove)

        z_lext_1(thisLevel) = tracer(jc,levelAbove)                             &
          & + 0.5_wp * (z_face_low(jc,levelAbove) - z_face_up(jc,levelAbove))   &
          &          * ( 1.0_wp - z_cfl_m )                                     &
          & - z_a11 - z_a11 * z_cfl_m * ( -3._wp  + 2._wp * z_cfl_m)

        ! second (of cell below) (case of w > 0; weta < 0)
        ! z_delta_p = (z_face_low(jc,thisLevel) - z_face_up(jc,thisLevel))
        z_a12     = tracer(jc,thisLevel)                                      &
          & - 0.5_wp * (z_face_low(jc,thisLevel) + z_face_up(jc,thisLevel))

        z_lext_2(thisLevel) = tracer(jc,thisLevel)                            &
          & - 0.5_wp * (z_face_low(jc,thisLevel) - z_face_up(jc,thisLevel))   &
          &          * ( 1.0_wp - z_cfl_p )                                   &
          & - z_a12 + z_a12 * z_cfl_p * (- 3._wp + 2._wp * z_cfl_p)
        !
        ! calculate vertical tracer flux
        !
        !upward_tracer_flux(jc,thisLevel) =                                    &
        !  & laxfr_upflux_v( w(jc,thisLevel),                                  &
        !  & z_lext_1(jc,thisLevel), z_lext_2(jc,thisLevel))
        ! copy of   FUNCTION laxfr_upflux_v( p_vn, p_psi1, p_psi2 )  RESULT(upward_tracer_flux)
        upward_tracer_flux(jc,thisLevel) = 0.5_wp *  &
          & ( w(jc,thisLevel)  * ( z_lext_1(thisLevel) + z_lext_2(thisLevel) )           &
          & +  ABS( w(jc,thisLevel) ) * ( z_lext_2(thisLevel) - z_lext_1(thisLevel) )    &
          & )

      END DO ! end loop over cells
    ENDDO ! end loop over vertical levels
    !
    ! set lower boundary condition
    !
    ! upward_tracer_flux(startIndex:endIndex,nlevp1) = 0.0_wp

    DO jc = startIndex, endIndex
! !CDIR NODEP
      DO thisLevel = firstLevel, cells_noOfLevels(jc)
        ! positive vertical divergence in direction of w (upward positive)
        flux_div_vert(jc,thisLevel) = upward_tracer_flux(jc, thisLevel) &
          & - upward_tracer_flux(jc, thisLevel+1)
      ENDDO
      DO thisLevel = cells_noOfLevels(jc)+1, n_zlev
        ! positive vertical divergence in direction of w (upward positive)
        flux_div_vert(jc,thisLevel) = 0.0_wp
      ENDDO
    END DO

  END SUBROUTINE upwind_vflux_ppm_onBlock
  !-------------------------------------------------------------------------
 

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Lax Friedrichs first order upwind flux for vertical advection,.
  !!
  !! Generalized Lax Friedrichs first order upwind flux,
  !! used in conservative vertical advection routines.
  !! For passive advection, equivalent to any other first
  !! order upwind flux.
  !! Applicable to both pressure based and height based vertical
  !! coordinate systems. Depending on the coordinate system chosen,
  !! the sign of the second term in the flux equation changes.
  !! - (-) for pressure based vertical coordinate systems
  !! - (+) for height based coordinate systems
  !! In order to get the correct sign, the variable p_coeff_grid
  !! has been introduced which is =1 for pressure based and =-1
  !! for height based coordinate systems.
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura  (2004).
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized for p- and z-based vertical coordinate systems
  !!
  ELEMENTAL FUNCTION laxfr_upflux_v( p_vn, p_psi1, p_psi2 )  result(p_upflux)
    
    REAL(wp), INTENT(in) :: p_vn
    REAL(wp), INTENT(in) :: p_psi1, p_psi2
    
    REAL(wp) :: p_upflux
    
    !-----------------------------------------------------------------------
    !p_upflux = 0.5_wp * (                       p_vn  *( p_psi1 + p_psi2 )    &
    !  &                   - p_coeff_grid * ABS( p_vn )*( p_psi2 - p_psi1 ) )
    p_upflux = 0.5_wp * (                       p_vn  *( p_psi1 + p_psi2 )    &
      & +  ABS( p_vn )*( p_psi2 - p_psi1 ) )
  END FUNCTION laxfr_upflux_v
  
END MODULE mo_ocean_tracer_transport_vert
!  !------------------------------------------------------------------------


  !------------------------------------------------------------------------
!   !>
!   !! First order upwind scheme for vertical tracer advection
!   !!
!   !! Calculation of vertical tracer fluxes
!   !!
!   !! @par Revision History
!   !!!! mpi parallelized, no sync
!   SUBROUTINE mimetic_vflux_oce( patch_2D, pvar_c, pw_c, pupflux_i, tracer_id )
!
!     TYPE(t_patch), TARGET, INTENT(IN) :: patch_2D           !< patch on which computation is performed
!     REAL(wp), INTENT(INOUT)           :: pvar_c(:,:,:)    !< advected cell centered variable
!     REAL(wp), INTENT(INOUT)           :: pw_c(:,:,:)      !< vertical velocity on cells
!     REAL(wp), INTENT(INOUT)           :: pupflux_i(:,:,:) !< variable in which the upwind flux is stored
!                                                           !< dim: (nproma,n_zlev+1,alloc_cell_blocks)
!     INTEGER, INTENT(IN)               :: tracer_id
!     ! local variables
!     INTEGER  :: startIndex, endIndex
!     INTEGER  :: jc, jk, jb               !< index of cell, vertical thisLevel and block
!     INTEGER  :: jkm1, jkp1, z_dolic                    !< jk - 1
!     REAL(wp) :: w_ave(n_zlev)
!     REAL(wp) :: w_avep1(n_zlev)
!     !-------------------------------------------------------------------------
!     TYPE(t_subset_range), POINTER :: cells_in_domain
!     !-------------------------------------------------------------------------
!     cells_in_domain => patch_2D%cells%in_domain
!
!     w_ave(:)  = 0.0_wp
!     w_avep1(:)= 0.0_wp
!
!     DO jb = cells_in_domain%start_block, cells_in_domain%end_block
!       CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
!         DO jc = startIndex, endIndex
!           z_dolic = v_base%dolic_c(jc,jb)
!           IF(z_dolic>=MIN_DOLIC)THEN
!             DO jk = 2, z_dolic
!               jkm1 = jk - 1
!               jkp1 = jk + 1
!               IF(pw_c(jc,jk,jb)>pw_c(jc,jkm1,jb))THEN
!
!                 w_ave(jk)=pw_c(jc,jk,jb)*pvar_c(jc,jk,jb)
!               ELSE
!                 w_ave(jk)=pw_c(jc,jk,jb)*pvar_c(jc,jkm1,jb)
!               ENDIF
!              !w_avep1(jk)=0.5_wp*(pw_c(jc,jk,jb)+pw_c(jc,jkp1,jb))&
!              !       &*pvar_c(jc,jk,jb)
!
!             !pupflux_i(jc,jk,jb) = 0.5_wp*(w_ave(jk) +w_avep1(jk))
!             !pupflux_i(jc,jk,jb) = (w_ave(jk)*v_base%del_zlev_m(jk) &
!             !                    & +w_avep1(jk)*v_base%del_zlev_m(jkp1))&
!             !                   & /(v_base%del_zlev_m(jk)+v_base%del_zlev_m(jkp1))
!             END DO
!             DO jk=1,z_dolic-1
!               jkm1 = jk - 1
!               jkp1 = jk + 1
!               pupflux_i(jc,jk,jb) = (w_ave(jk)*v_base%del_zlev_m(jk) &
!                                  & +w_avep1(jk)*v_base%del_zlev_m(jkp1))&
!                                  & /(v_base%del_zlev_m(jk)+v_base%del_zlev_m(jkp1))
!             END DO
!             !w_ave(z_dolic)=0.5_wp*(pw_c(jc,z_dolic,jb)+pw_c(jc,z_dolic-1,jb))&
!             !        &*pvar_c(jc,z_dolic,jb)
!             !pupflux_i(jc,z_dolic,jb) = w_ave(z_dolic)*v_base%del_zlev_m(z_dolic)&
!             !                          &/v_base%del_zlev_i(z_dolic)
!             ! no fluxes at bottom boundary
!             pupflux_i(jc,z_dolic+1,jb) = 0.0_wp
!         ENDIF
!       END DO
!     END DO
! !     DO jb = all_cells%start_block, all_cells%end_block
! !       CALL get_index_range(all_cells, jb, startIndex, endIndex)
! !         DO jc = startIndex, endIndex
! !           DO jk = 2, n_zlev
! !            ! index of top & bottom half thisLevel
! !            jkm1 = jk - 1
! !            jkp1 = jk + 1
! !              w_ave(jk)=0.5_wp*(pw_c(jc,jk,jb)+pw_c(jc,jkm1,jb))&
! !                    &*pvar_c(jc,jk,jb)
! !            END DO ! end cell loop
! !        END DO ! end thisLevel loop
! !      END DO ! end block loop
! !     DO jb = all_cells%start_block, all_cells%end_block
! !       CALL get_index_range(all_cells, jb, startIndex, endIndex)
! !         DO jc = startIndex, endIndex
! !          DO jk = 1, n_zlev-1
! !            ! index of top & bottom half thisLevel
! !            jkm1 = jk - 1
! !            jkp1 = jk + 1
! !             w_avep1(jk)=0.5_wp*(pw_c(jc,jkp1,jb)+pw_c(jc,jk,jb))&
! !                    &*pvar_c(jc,jkp1,jb)
! !            END DO ! end cell loop
! !        END DO ! end thisLevel loop
! !      END DO ! end block loop
! !     DO jb = all_cells%start_block, all_cells%end_block
! !       CALL get_index_range(all_cells, jb, startIndex, endIndex)
! !         DO jc = startIndex, endIndex
! !          DO jk = 1, n_zlev-1
! !            ! index of top & bottom half thisLevel
! !            jkm1 = jk - 1
! !            jkp1 = jk + 1
! !            pupflux_i(jc,jk,jb) = 0.5_wp*(w_ave(jk)  *v_base%del_zlev_m(jk) &
! !                                &        +w_avep1(jk)*v_base%del_zlev_m(jkp1))/v_base%del_zlev_i(jk)
! !            END DO ! end cell loop
! !        END DO ! end thisLevel loop
! !      END DO ! end block loop
!   END SUBROUTINE mimetic_vflux_oce



!  SUBROUTINE upwind_vflux_ppm_onBlock(  &
!     & tracer,                           &
!     & w, dtime, vertical_limiter_type,  &
!     & cell_thickeness, cell_invheight,      &
!     & flux_div_vert,               &
!     & startIndex, endIndex, cells_noOfLevels)
!
!     REAL(wp), INTENT(in)           :: tracer(nproma,n_zlev)      !< advected cell centered variable
!     REAL(wp), INTENT(in)           :: w(nproma,n_zlev+1)         !<  in : vertical velocity
!     REAL(wp), INTENT(in)           :: dtime                      !< time step
!     REAL(wp), INTENT(in)           :: cell_thickeness(nproma,n_zlev) !< layer thickness at cell center at time n
!     REAL(wp), INTENT(in)           :: cell_invheight(nproma,n_zlev)!< layer thickness at cell center at time n
!     INTEGER, INTENT(in)            :: vertical_limiter_type                    !< parameter to select limiter
!     REAL(wp), INTENT(inout)        :: flux_div_vert(nproma, n_zlev) !new tracer
!     INTEGER, INTENT(in)            :: startIndex, endIndex
!     INTEGER, INTENT(in)            :: cells_noOfLevels(nproma)
!     !
!     REAL(wp) :: upward_tracer_flux(nproma,n_zlev+1)      !< tracer flux
!     REAL(wp) :: z_face(nproma,n_zlev+1)   !< face values of transported field
!     REAL(wp) :: z_face_up(nproma,n_zlev)  !< face value (upper face)
!     REAL(wp) :: z_face_low(nproma,n_zlev) !< face value (lower face)
!     REAL(wp) :: z_lext_1(n_zlev+1)                 !< linear extrapolation value 1
!     REAL(wp) :: z_lext_2(n_zlev+1)                 !< linear extrapolation value 2
!     REAL(wp) :: z_cfl_m, z_cfl_p !< CFL number (weta>0, w<0), CFL number (weta<0, w>0)
!     REAL(wp) :: z_slope(nproma,n_zlev+1)  !< monotonized slope
!     REAL(wp) :: z_slope_u, z_slope_l                            !< one-sided slopes
!     !< for weta >0 and weta <0
!     REAL(wp) :: z_a11, z_a12                                    !< 1/6 * a6,i (see Colella and Woodward (1984))
!     REAL(wp) :: z_weta_dt                                       !< weta times dtime
!     INTEGER :: firstLevel, secondLevel                                    !< vertical start thisLevel and start thisLevel +1
!     INTEGER :: levelAbove, levelBelow,level2Below                        !< vertical thisLevel minus and plus one, plus two
!     INTEGER :: jc, thisLevel, cell_levels
!     !LOGICAL  :: opt_lout_edge !< optional: output edge value (.TRUE.),
!     !                          !< or the flux across the edge   !< (.FALSE./not specified)
!     !REAL(wp) :: opt_topflx_tra(nproma,patch_3D%p_patch_2D(1)%alloc_cell_blocks)  !< vertical tracer flux at upper boundary
!     INTEGER, PARAMETER :: islopel_vsm = 1
!     !-----------------------------------------------------------------------
!
!     firstLevel  = 1
!     secondLevel = 2
!
!     ! advection is done with an upwind scheme and a piecwise parabolic
!     ! approx. of the subgrid distribution is used.
!     ! 3 options:  standard without limiter
!     !             standard with semi-monotone or monotone limiter
!     !             special version with limiter which handles CFL >1
!     !
!     !------------------------------------------------
!     ! 1. Calculate Courant number for weta>0 (w<0) and weta<0 (w>0)
!     ! ..............
!
!     ! 2. Calculate monotonized slope
!     !
!     z_slope(:, :) = 0._wp
!
! ! !CDIR NODEP
!     DO jc = startIndex, endIndex
!
! ! !CDIR NODEP
!       DO thisLevel = secondLevel, cells_noOfLevels(jc) - 1
!
!         levelAbove    = thisLevel - 1                 ! index of top half thisLevel
!         levelBelow    = thisLevel + 1  ! index of bottom half thisLevel
!
!         z_slope_u = 2._wp * (tracer(jc,thisLevel)      - tracer(jc,levelAbove))
!         z_slope_l = 2._wp * (tracer(jc,levelBelow) - tracer(jc,thisLevel))
!
!         IF ((z_slope_u * z_slope_l) > 0._wp) THEN
!
!           z_slope(jc,thisLevel) = &
!             & ( cell_thickeness(jc,thisLevel)                             &
!             & / (cell_thickeness(jc,levelAbove) + cell_thickeness(jc,thisLevel)            &
!             & + cell_thickeness(jc,levelBelow)) )                                       &
!             & * ( &
!             & (2._wp * cell_thickeness(jc,levelAbove) + cell_thickeness(jc,thisLevel)) &
!             & / (cell_thickeness(jc,levelBelow) + cell_thickeness(jc,thisLevel))           &
!             & &
!             & * (tracer(jc,levelBelow) - tracer(jc,thisLevel))                         &
!             & &
!             & + (cell_thickeness(jc,thisLevel) + 2._wp * cell_thickeness(jc,levelBelow))   &
!             & / (cell_thickeness(jc,levelAbove) + cell_thickeness(jc,thisLevel))           &
!             & &
!             & * (tracer(jc,thisLevel) - tracer(jc,levelAbove)) )
!
!           z_slope(jc,thisLevel) = SIGN(                                            &
!             & MIN( ABS(z_slope(jc,thisLevel)), ABS(z_slope_u), ABS(z_slope_l) ),  &
!             & z_slope(jc,thisLevel))
!
!         ENDIF
!
!       END DO ! jc = startIndex, endIndex
!     END DO ! end loop over vertical levels
!
!     !
!     ! 3. reconstruct face values at vertical half-levels
!     !
!     ! Boundary values for two highest and lowest half-levels
!     !
!     ! for faces k=secondLevel and k=nlevp1-1 reconstructed face values are calculated by
!     ! interpolating a quadratic (instead of quartic) polynomial through 3
!     ! values of the indefinite integral A=\int_{\eta_{0}}^{\eta}q\,\mathrm{d}\eta
!     !
!     ! for faces k=firstLevel and k=nlevp1 a zero gradient condition is assumed and the
!     ! face values are set to the tracer values of the corresponding cell centers
!     !
!     z_face(:, :)  = 0.0_wp
!
! ! !CDIR NODEP
!     DO jc = startIndex, endIndex
!       cell_levels = cells_noOfLevels(jc)
!
! ! !CDIR NODEP
!       DO thisLevel = secondLevel, cell_levels - 2
!         levelAbove  = thisLevel - 1
!         levelBelow  = thisLevel + 1
!         level2Below = thisLevel + 2
!
!         z_face(jc,levelBelow) = tracer(jc,thisLevel) &
!           & + (cell_thickeness(jc,thisLevel)                                           &
!           & / (cell_thickeness(jc,thisLevel) + cell_thickeness(jc,levelBelow)))         &
!           & &
!           & * (tracer(jc,levelBelow) - tracer(jc,thisLevel))                    &
!           & &
!           & + (1._wp/(cell_thickeness(jc,levelAbove) + cell_thickeness(jc,thisLevel)    &
!           & + cell_thickeness(jc,levelBelow) + cell_thickeness(jc,level2Below)))    &
!           & &
!           & *( &
!           & (2._wp * cell_thickeness(jc,levelBelow) * cell_thickeness(jc,thisLevel)     &
!           & / (cell_thickeness(jc,thisLevel) + cell_thickeness(jc,levelBelow)))         &
!           & &
!           & * ( (cell_thickeness(jc,levelAbove) + cell_thickeness(jc,thisLevel))        &
!           & / (2._wp*cell_thickeness(jc,thisLevel) + cell_thickeness(jc,levelBelow))    &
!           & - (cell_thickeness(jc,level2Below) + cell_thickeness(jc,levelBelow))    &
!           & / (2._wp*cell_thickeness(jc,levelBelow) + cell_thickeness(jc,thisLevel)) )  &
!           & &
!           & * (tracer(jc,levelBelow) - tracer(jc,thisLevel)) -                  &
!           & z_slope(jc,levelBelow) *                                        &
!           & &
!           & cell_thickeness(jc,thisLevel) *                                         &
!           & (cell_thickeness(jc,levelAbove) + cell_thickeness(jc,thisLevel)) /          &
!           & (2._wp*cell_thickeness(jc,thisLevel) + cell_thickeness(jc,levelBelow))      &
!           & &
!           & +  z_slope(jc,thisLevel) *                                                &
!           & &
!           & cell_thickeness(jc,levelBelow) * &
!           & (cell_thickeness(jc,levelBelow)  + cell_thickeness(jc,level2Below)) &
!           & / (cell_thickeness(jc,thisLevel) + 2._wp*cell_thickeness(jc,levelBelow)) )
!
!       END DO ! end loop over vertical levels
!
!       ! compute top 2 levels
!       IF ( cells_noOfLevels(jc) >= firstLevel  ) THEN
!         z_face(jc,firstLevel) = tracer(jc,firstLevel)
!         IF ( cell_levels >= secondLevel ) THEN
!
!           z_face(jc,secondLevel) = &
!             & tracer(jc,firstLevel) *                           &
!             & (1._wp - (cell_thickeness(jc,firstLevel) / cell_thickeness(jc,secondLevel)))          &
!             & &
!             & + (cell_thickeness(jc,firstLevel)  &
!             & /(cell_thickeness(jc,firstLevel) + cell_thickeness(jc,secondLevel)))   &
!             & * &
!             & ((cell_thickeness(jc,firstLevel) / cell_thickeness(jc,secondLevel))  &
!             & * tracer(jc,firstLevel) &
!             & + tracer(jc,secondLevel))
!
!         ENDIF
!       ENDIF
!
!       ! compute bottom thisLevel
!       IF ( cells_noOfLevels(jc) > secondLevel ) THEN
!
!         z_face(jc, cell_levels) = &
!           & tracer(jc,cell_levels-1) * &
!           & ( 1._wp - (cell_thickeness(jc,cell_levels-1) / cell_thickeness(jc,cell_levels)))  &
!           & &
!           & + (cell_thickeness(jc,cell_levels-1)/(cell_thickeness(jc,cell_levels-1)   &
!           & + cell_thickeness(jc,cell_levels))) &
!           & * &
!           & ((cell_thickeness(jc,cell_levels-1) / cell_thickeness(jc,cell_levels)) &
!           & * tracer(jc,cell_levels-1)                &
!           & + tracer(jc,cell_levels))
!
!       ENDIF
!
!     ENDDO
!
!     ! 4. Limitation of first guess parabola (which is based on z_face)
!     ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
!     ! the limitation procedure.
!     ! Therefore 2 additional fields z_face_up and z_face_low are
!     ! introduced.
!     z_face_low(1:nproma,1:n_zlev) = 0.0_wp
!     z_face_up (1:nproma,1:n_zlev) = 0.0_wp
!
!     IF (vertical_limiter_type == islopel_vsm) THEN
!       !     ! monotonic (mo) limiter
!       CALL v_ppm_slimiter_mo_onBlock( &
!         & tracer(:,:),       &
!         & z_face(:,:),       &
!         & z_slope(:,:),      &
!         & z_face_up,         &
!         & z_face_low,        &
!         & startIndex,        &
!         & endIndex,          &
!         & cells_noOfLevels)
!
!     ELSE
!         ! simply copy face values to 'face_up' and 'face_low' arrays
!         DO jc = startIndex, endIndex
!           DO thisLevel = secondLevel, cells_noOfLevels(jc)-1
!             z_face_up(jc, thisLevel)  = z_face(jc, thisLevel    )
!             z_face_low(jc,thisLevel)  = z_face(jc, thisLevel + 1)
!           ENDDO
!         END DO
!
!     ENDIF  !  p_ityp_vlimit
!
!
!     upward_tracer_flux(:,:) = 0.0_wp
!
! ! !CDIR NODEP
!     DO jc = startIndex, endIndex
! ! !CDIR NODEP
!       DO thisLevel = secondLevel, cells_noOfLevels(jc)
!         ! index of top half thisLevel
!         levelAbove = thisLevel - 1
!         ! linear extrapolated values
!         ! for the height based coordinate system multiplication by coeff_grid
!         ! is not necessary due to compensating (-) signs.
!         ! first (of cell above) (case of w < 0; weta > 0)
!         ! z_delta_m = (z_face_low(jc,levelAbove) - z_face_up(jc,levelAbove))
!         z_a11     = tracer(jc,levelAbove)                                  &
!           & - 0.5_wp * (z_face_low(jc,levelAbove) + z_face_up(jc,levelAbove))
!
!         ! Calculate local Courant number at half levels
!         ! z_cfl_m for weta >0 (w <0)
!         ! z_cfl_p for weta <0 (w >0)
!         ! z_weta_dt = 0.0_wp
!         z_weta_dt = ABS(w(jc,thisLevel)) * dtime
!         z_cfl_p = z_weta_dt * cell_invheight(jc, thisLevel)
!         z_cfl_m = z_weta_dt * cell_invheight(jc, levelAbove)
!
!         z_lext_1(thisLevel) = tracer(jc,levelAbove)                           &
!           & + 0.5_wp * (z_face_low(jc,levelAbove) - z_face_up(jc,levelAbove))   &
!           &          * ( 1.0_wp - z_cfl_m )                                     &
!           & - z_a11 - z_a11 * z_cfl_m * ( -3._wp  + 2._wp * z_cfl_m)
!
!         ! second (of cell below) (case of w > 0; weta < 0)
!         ! z_delta_p = (z_face_low(jc,thisLevel) - z_face_up(jc,thisLevel))
!         z_a12     = tracer(jc,thisLevel)                                      &
!           & - 0.5_wp * (z_face_low(jc,thisLevel) + z_face_up(jc,thisLevel))
!
!         z_lext_2(thisLevel) = tracer(jc,thisLevel)                            &
!           & - 0.5_wp * (z_face_low(jc,thisLevel) - z_face_up(jc,thisLevel))         &
!           &          * ( 1.0_wp - z_cfl_p )                                 &
!           & - z_a12 + z_a12 * z_cfl_p * (- 3._wp + 2._wp * z_cfl_p)
!         !
!         ! calculate vertical tracer flux
!         !
!         !upward_tracer_flux(jc,thisLevel) =                                    &
!         !  & laxfr_upflux_v( w(jc,thisLevel),                                  &
!         !  & z_lext_1(jc,thisLevel), z_lext_2(jc,thisLevel))
!         ! copy of   FUNCTION laxfr_upflux_v( p_vn, p_psi1, p_psi2 )  RESULT(upward_tracer_flux)
!         upward_tracer_flux(jc,thisLevel) = 0.5_wp *  &
!           & ( w(jc,thisLevel)  * ( z_lext_1(thisLevel) + z_lext_2(thisLevel) )           &
!           & +  ABS( w(jc,thisLevel) ) * ( z_lext_2(thisLevel) - z_lext_1(thisLevel) )    &
!           & )
!
!       END DO ! end loop over cells
!     ENDDO ! end loop over vertical levels
!     !
!     ! set lower boundary condition
!     !
!     ! upward_tracer_flux(startIndex:endIndex,nlevp1) = 0.0_wp
!
!     DO jc = startIndex, endIndex
! ! !CDIR NODEP
!       DO thisLevel = firstLevel, cells_noOfLevels(jc)
!         ! positive vertical divergence in direction of w (upward positive)
!         flux_div_vert(jc,thisLevel) = upward_tracer_flux(jc, thisLevel) &
!           & - upward_tracer_flux(jc, thisLevel+1)
!       ENDDO
!     END DO
!
!   END SUBROUTINE upwind_vflux_ppm_onBlock
!   !-------------------------------------------------------------------------
