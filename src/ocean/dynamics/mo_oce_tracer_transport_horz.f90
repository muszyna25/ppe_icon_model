!>
!! Contains the implementation of the horizontal tracer transport routines for the ICON ocean model.
!! This comprises horizontal advection and diffusion
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/01)
!!
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
!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_oce_tracer_transport_horz
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  USE mo_math_constants,            ONLY: dbl_eps
  USE mo_impl_constants,            ONLY: sea_boundary
  USE mo_ocean_nml,                 ONLY: n_zlev, l_edge_based, ab_gam,            &
    & upwind, central,mimetic, mimetic_miura,  &
    & flux_calculation_horz,                   &
    & l_with_horz_tracer_diffusion,            &
    & l_with_horz_tracer_advection,            &
    & l_horz_limiter_advection
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_adv_horz, timer_hflx_lim, &
    & timer_dif_horz
  USE mo_oce_state,                 ONLY: t_hydro_ocean_state
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish !, message_text, message
  USE mo_oce_physics
  USE mo_scalar_product,            ONLY: map_cell2edges_3d,map_edges2cell_3d, &
    & map_edges2edges_viacell_3d_const_z
  USE mo_oce_math_operators,        ONLY: div_oce_3d, grad_fd_norm_oce_3d
  USE mo_oce_diffusion,             ONLY: tracer_diffusion_horz
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_c, sync_c1, sync_e, sync_patch_array, sync_patch_array_mult
  USE mo_mpi,                       ONLY: my_process_is_mpi_parallel
  
  IMPLICIT NONE
  
  PRIVATE
  
  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  CHARACTER(LEN=12)           :: str_module    = 'oceTracHorz '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  !
  ! PUBLIC INTERFACE
  !
  PUBLIC :: advect_diffuse_flux_horz
  ! Private implemenation
  !
  
  PRIVATE :: upwind_hflux_oce
  PRIVATE :: central_hflux_oce
  PRIVATE :: upwind_hflux_oce_mimetic
  PRIVATE :: central_hflux_oce_mimetic
  PRIVATE :: mimetic_miura_hflux_oce
  PRIVATE :: hflx_limiter_oce_mo
  PRIVATE :: laxfr_upflux
  !PRIVATE :: elad
  
  
  INTEGER, PARAMETER :: top=1
CONTAINS
  !-----------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects horizontally the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !!
  !!  mpi parallelized LL
  SUBROUTINE advect_diffuse_flux_horz( patch_3d,          &
    & trac_old,            &
    & p_os,                &
    & p_op_coeff,          &
    & k_h,                 &
    & h_old,               &
    & h_new,               &
    & flux_horz)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp)                               :: trac_old(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET :: p_os
    TYPE(t_operator_coeff), INTENT(in)     :: p_op_coeff
    REAL(wp), INTENT(in)                   :: k_h(:,:,:)         !horizontal mixing coeff
    REAL(wp), INTENT(in)                   :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                   :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)                :: flux_horz(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !
    !Local variables
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: jc, jk, jb, je
    
    REAL(wp) :: z_adv_flux_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)  ! horizontal advective tracer flux
    REAL(wp) :: z_div_adv_h  (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)   ! horizontal tracer divergence
    REAL(wp) :: z_div_diff_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)  ! horizontal tracer divergence
    REAL(wp) :: z_diff_flux_h(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e) ! horizontal diffusive tracer flux
    
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_transport_horz:advect_diffuse_flux_horz')
    TYPE(t_patch), POINTER :: patch_2d
    !-------------------------------------------------------------------------------
    patch_2d         => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------------
    !z_vn         (1:nproma,1:n_zlev,1:patch_2d%nblks_e)=0.0_wp
    z_adv_flux_h  (:,:,:) = 0.0_wp
    z_div_adv_h   (:,:,:) = 0.0_wp
    z_div_diff_h  (:,:,:) = 0.0_wp
    z_diff_flux_h (:,:,:) = 0.0_wp
    
    CALL sync_patch_array(sync_e, patch_2d, p_os%p_diag%vn_time_weighted)
    
    !Calculate tracer fluxes at edges
    !This step takes already the edge length into account
    !but not the edge height
    IF ( l_with_horz_tracer_advection ) THEN
      
      ! Initialize timer for horizontal advection
      IF (ltimer) CALL timer_start(timer_adv_horz)
      
      SELECT CASE(flux_calculation_horz)
      
      CASE(upwind)
        
        !upwind estimate of tracer flux
        CALL upwind_hflux_oce( patch_3d,                   &
          & trac_old,                     &
          & p_os%p_diag%vn_time_weighted, &
          & z_adv_flux_h )
      CASE(central)
        
        !central estimate of tracer flux
        CALL central_hflux_oce( patch_3d,                   &
          & trac_old,                     &
          & p_os%p_diag%vn_time_weighted, &
          & z_adv_flux_h )
        
      CASE(mimetic_miura)
        
        !MIURA-type upwind-biased estimate of tracer flux
        CALL mimetic_miura_hflux_oce( patch_3d,                 &
          & trac_old,                     &
          & p_os%p_diag%vn_time_weighted, &
          & p_op_coeff,                   &
          & z_adv_flux_h )
      CASE(mimetic)!For non-edged-based option this is handled completely below
        IF(l_edge_based)THEN
          CALL finish('TRIM(advect_diffuse_flux_horz)',"Incompatible options: Mimetic flux together with edge-based")
        ENDIF
        
      CASE default
        CALL finish('TRIM(advect_diffuse_flux_horz)',"This flux option is not supported")
        
      END SELECT
      
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=5  ! output print level (1-5, fix)
      CALL dbg_print('aft. AdvHorz: adv_flux_h',z_adv_flux_h,str_module,idt_src,patch_2d%edges%owned)
      !---------------------------------------------------------------------
      
      !Multiply fluxes with edge height
      ! !-------------------------------------------------------------------------------
      IF( l_edge_based)THEN
        ! !-------------------------------------------------------------------------------
        !compute new edge
        DO jb = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
          DO je = i_startidx_e, i_endidx_e
            DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
              z_adv_flux_h(je,jk,jb) = z_adv_flux_h(je,jk,jb) * patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,jb)
            END DO
          END DO
        END DO
        
        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=5  ! output print level (1-5, fix)
        CALL dbg_print('AdvH*thick_e: adv_flux_h', z_adv_flux_h,str_module, idt_src, patch_2d%edges%owned)
        !---------------------------------------------------------------------
        
        ! !-------------------------------------------------------------------------------
      ELSEIF(.NOT. l_edge_based)THEN
        ! !-------------------------------------------------------------------------------
        !Identical procedure of depth multiplication for all options
        !except the mimetic one
        !    IF(    FLUX_CALCULATION_HORZ==UPWIND &
        !      &.OR.FLUX_CALCULATION_HORZ==CENTRAL&
        !      &.OR.FLUX_CALCULATION_HORZ==MIMETIC_MIURA)THEN
        !       !Multiply with height at edges
        !       !Non-time-dependent contribution first (already at edges)
        !       DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        !         CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        !         DO jk = 1, n_zlev
        !           DO je = i_startidx_e, i_endidx_e
        !            z_adv_flux_h(je,jk,jb) = z_adv_flux_h(je,jk,jb)*patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_e(je,jk,jb)
        !             !ENDIF
        !           END DO
        !         END DO
        !       END DO
        !       !Now time-dependent contribution (at cell centers)
        !       DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        !         CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        !         DO jc = i_startidx_c, i_endidx_c
        !           !Mass flux is already multiplied by thickness
        !           z_vn_c(jc,jb)%x = trac_old(jc,1,jb)*p_os%p_diag%p_mass_flux_sfc_cc(jc,jb)%x
        !           !ENDIF
        !         END DO
        !       END DO
        !       CALL sync_patch_array(SYNC_C,patch_2d,z_vn_c(1:nproma,1:patch_2d%nblks_c)%x(1))
        !       CALL sync_patch_array(SYNC_C,patch_2d,z_vn_c(1:nproma,1:patch_2d%nblks_c)%x(2))
        !       CALL sync_patch_array(SYNC_C,patch_2d,z_vn_c(1:nproma,1:patch_2d%nblks_c)%x(3))
        !       !Map time-dependent contribution back to edges.
        !       !This is a 2D-mapping
        !       CALL map_cell2edges_3D( patch_3D,     &
        !                             & z_vn_c(:,:),    &
        !                             & z_flux_2D(:,:), &
        !                             & p_op_coeff,     &
        !                             & level=1)
        !       CALL sync_patch_array(SYNC_E,patch_2d,z_flux_2D)
        !       DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        !         CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        !         DO je = i_startidx_e, i_endidx_e
        !           IF ( patch_3D%lsm_e(je,1,jb) <= sea_boundary ) THEN
        !              z_adv_flux_h(je,1,jb) =  z_adv_flux_h(je,1,jb) + z_flux_2D(je,jb)
        !           ENDIF
        !         END DO
        !       END DO
        !    ELSEIF(FLUX_CALCULATION_HORZ==MIMETIC)THEN
        
        CALL map_edges2edges_viacell_3d_const_z( patch_3d,                  &
          & p_os%p_diag%vn_time_weighted,&
          & p_op_coeff,                  &
          & z_adv_flux_h,                &
          & trac_old)
        CALL sync_patch_array(sync_e, patch_2d, z_adv_flux_h)
        
        !ENDIF
        ! !------------------------------------------------------------------------------
      ENDIF!l_edge_based
      ! !-------------------------------------------------------------------------------
      ! Stop timer for horizontal advection
      IF (ltimer) CALL timer_stop(timer_adv_horz)
      
      ! Flux limiting process, dependent on tracer configuration
      !IF(FLUX_CALCULATION_HORZ/=UPWIND)THEN
      IF (l_horz_limiter_advection) THEN
        
        ! Start timer for horizontal flux limitation
        IF (ltimer) CALL timer_start(timer_hflx_lim)
        
        CALL hflx_limiter_oce_mo( patch_3d,               &
          & trac_old,               &
          & p_os%p_diag%mass_flx_e, &
          & z_adv_flux_h,           &
          & p_op_coeff,             &
          & h_old,h_new)
        IF (ltimer) CALL timer_stop(timer_hflx_lim)
        
        !---------DEBUG DIAGNOSTICS-------------------------------------------
        idt_src=3  ! output print level (1-5, fix)
        CALL dbg_print('aft. HorzLim: adv_flux_h',z_adv_flux_h,str_module,idt_src,patch_2d%edges%owned)
        !---------------------------------------------------------------------
        
      ENDIF
      
      !Calculate divergence of advective fluxes
      CALL div_oce_3d( z_adv_flux_h, patch_2d,p_op_coeff%div_coeff, z_div_adv_h,&
        & subset_range=cells_in_domain)
      
    ELSE
      z_div_adv_h   (:,:,:) = 0.0_wp
    ENDIF ! l_with_horz_tracer_diffusion
    
    !The diffusion part: calculate horizontal diffusive flux
    IF ( l_with_horz_tracer_diffusion ) THEN
      IF (ltimer) CALL timer_start(timer_dif_horz)
      CALL tracer_diffusion_horz( patch_3d,     &
        & trac_old,     &
        & p_os,         &
        & k_h,          &
        & z_diff_flux_h,&
        & subset_range = edges_in_domain)
      
      !Calculate divergence of diffusive fluxes
      CALL div_oce_3d( z_diff_flux_h, patch_2d,p_op_coeff%div_coeff, z_div_diff_h, &
        & subset_range=cells_in_domain)
      IF (ltimer) CALL timer_stop(timer_dif_horz)
      
    ELSE
      z_diff_flux_h (:,:,:) = 0.0_wp
      z_div_diff_h  (:,:,:) = 0.0_wp
    ENDIF
    
    !Final step: calculate sum of advective and diffusive horizontal fluxes
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          flux_horz(jc,jk,jb) = z_div_diff_h(jc,jk,jb)-z_div_adv_h(jc,jk,jb)
        END DO
      END DO
    END DO
    
    CALL sync_patch_array(sync_c, patch_2d, flux_horz)
    
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('AdvDifHorz: adv_flux_h'     ,z_adv_flux_h ,str_module,idt_src,patch_2d%edges%owned)
    CALL dbg_print('AdvDifHorz: div adv_flux_h' ,z_div_adv_h  ,str_module,idt_src,cells_in_domain)
    CALL dbg_print('AdvDifHorz: div diff_flux_h',z_div_diff_h ,str_module,idt_src,cells_in_domain)
    CALL dbg_print('AdvDifHorz: flux_horz'      ,flux_horz    ,str_module,idt_src,cells_in_domain)
    !---------------------------------------------------------------------
    
  END SUBROUTINE advect_diffuse_flux_horz
  !-------------------------------------------------------------------------------
  !>
  !! First order upwind scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !!
  !!  mpi parallelized LL
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  SUBROUTINE upwind_hflux_oce_mimetic( patch_3d, flux_cc,p_op_coeff,&
    & pupflux_e, opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    TYPE(t_cartesian_coordinates)      :: flux_cc(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_operator_coeff), INTENT(in) :: p_op_coeff
    REAL(wp), INTENT(inout)            :: pupflux_e(nproma,n_zlev, patch_3d%p_patch_2d(1)%nblks_e)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level    ! optional vertical end level
    ! local variables
    !INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc
    INTEGER :: start_level, end_level
    INTEGER :: start_index, end_index
    INTEGER :: je, jk, jb, il_c, ib_c         !< index of edge, vert level, block
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
    !iilc => patch_2d%edges%cell_idx
    !iibc => patch_2d%edges%cell_blk
    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    ! for no-slip boundary conditions, boundary treatment for tracer (zero at leteral walls)
    !is implicit done via velocity boundary conditions
    !
    ! line and block indices of two neighboring cells
    ! loop through all patch edges (and blocks)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      
#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO jk = start_level, end_level
        DO je = start_index, end_index
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          IF ( patch_3d%lsm_e(je,jk,jb) <= sea_boundary ) THEN
            
            il_c=p_op_coeff%upwind_cell_idx(je,jk,jb)
            ib_c=p_op_coeff%upwind_cell_blk(je,jk,jb)
            
            pupflux_e(je,jk,jb)=&
              & DOT_PRODUCT(flux_cc(il_c,jk,ib_c)%x,patch_2d%edges%primal_cart_normal(je,jb)%x)
            
            !    write(999,*)'upwind flux -h',je,jk,jb,z_pub_flux_e_up(je,jk,jb),pupflux_e(je,jk,jb)&
            !    &, pvn_e(je,jk,jb),&
            !    &pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2))
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
    
    
  END SUBROUTINE upwind_hflux_oce_mimetic
  !-----------------------------------------------------------------------
  !>
  !! First order upwind scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !!
  !!  mpi parallelized LL
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  SUBROUTINE central_hflux_oce_mimetic( patch_3d, flux_cc,&
    & pupflux_e, opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d                                  !< patch on which computation is performed
    TYPE(t_cartesian_coordinates)     :: flux_cc  (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)           :: pupflux_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)   !< variable in which the upwind flux is stored
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
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO jk = start_level, end_level
        DO je = start_index, end_index
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          IF ( patch_3d%lsm_e(je,jk,jb) <= sea_boundary ) THEN
            flux_mean_e%x=0.5_wp*(flux_cc(iilc(je,jb,1),jk,iibc(je,jb,1))%x&
              & +flux_cc(iilc(je,jb,2),jk,iibc(je,jb,2))%x)
            
            pupflux_e(je,jk,jb)=DOT_PRODUCT(flux_mean_e%x,&
              & patch_2d%edges%primal_cart_normal(je,jb)%x)
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
    
  END SUBROUTINE central_hflux_oce_mimetic
  
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
  !!  mpi parallelized LL
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  SUBROUTINE upwind_hflux_oce( patch_3d, pvar_c, pvn_e, pupflux_e, opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)              :: pvar_c   (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)      !< advected cell centered variable
    REAL(wp), INTENT(in)              :: pvn_e    (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)       !< normal velocity on edges
    REAL(wp), INTENT(inout)             :: pupflux_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)   !< variable in which the upwind flux is stored
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
    ! for no-slip boundary conditions, boundary treatment for tracer (zero at leteral walls)
    !is implicit done via velocity boundary conditions
    !
    ! line and block indices of two neighboring cells
    iilc => patch_2d%edges%cell_idx
    iibc => patch_2d%edges%cell_blk
    
    ! loop through all patch edges (and blocks)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      
#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO je = start_index, end_index
        DO jk = start_level, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          pupflux_e(je,jk,jb) =  &
            & laxfr_upflux( pvn_e(je,jk,jb), pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), &
            & pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
    
  END SUBROUTINE upwind_hflux_oce
  !-----------------------------------------------------------------------
  !>
  !! Central scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using central fluxes
  !!
  !! @par Revision History
  !! Peter korn, MPI-M, 2011
  !!
  !!  mpi parallelized LL
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  SUBROUTINE central_hflux_oce( patch_3d, pvar_c, pvn_e, pupflux_e )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)           :: pvar_c   (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    !< advected cell centered variable
    REAL(wp), INTENT(inout)           :: pvn_e    (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     !< normal velocity on edges
    REAL(wp), INTENT(inout)           :: pupflux_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) !< variable in which the upwind flux is stored
    
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
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO jk = 1, n_zlev
        DO je = start_index, end_index
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          IF ( patch_3d%lsm_e(je,jk,jb) <= sea_boundary ) THEN
            pupflux_e(je,jk,jb) =  0.5_wp*pvn_e(je,jk,jb)             &
              & *( pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1))      &
              & +pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)))     !&
            !               &          +0.5_wp*pvn_e(je,jk,jb)*pvn_e(je,jk,jb)*dtime&
            !               &          * patch_2d%edges%inv_dual_edge_length(je,jb)   &
            !               &        *( pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2))      &
            !               &          -pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)))
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
  END SUBROUTINE central_hflux_oce
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
  SUBROUTINE mimetic_miura_hflux_oce( patch_3d,pvar_c, pvn_e, &
    & p_op_coeff, pflux_e,&
    & opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)           :: pvar_c(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< advected cell centered variable
    REAL(wp), INTENT(inout)           :: pvn_e (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) !< normal velocity on edges
    TYPE(t_operator_coeff)            :: p_op_coeff
    REAL(wp), INTENT(inout)           :: pflux_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)!< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level      ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level      ! optional vertical end level
    
    ! local variables
    INTEGER :: start_level, end_level
    INTEGER :: i_startidx_e, i_endidx_e
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
    ! !         &                i_startidx_e, i_endidx_e, rl_start,rl_end)
    ! !
    ! !       DO jk = start_level, end_level
    ! !         DO je = i_startidx_e, i_endidx_e          !
    ! !           IF ( v_base%lsm_e(je,jk,jb) <= sea_boundary ) THEN
    ! !              pupflux_e(je,jk,jb) =  &
    ! !              &  laxfr_upflux( pvn_e(je,jk,jb), pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), &
    ! !              &                                 pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )
    ! ! !write(*,*)'upwind miura',je,jk,jb,pupflux_e(je,jk,jb),pflux_e(je,jk,jb), pvn_e(je,jk,jb)!,&
    ! ! !&pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2))
    ! !           ELSE
    ! !             pupflux_e(je,jk,jb) = 0.0_wp
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
    
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
      DO jk = start_level, end_level
        DO je =  i_startidx_e, i_endidx_e
          IF ( patch_3d%lsm_e(je,jk,jb) <= sea_boundary ) THEN
            il_c=p_op_coeff%upwind_cell_idx(je,jk,jb)
            ib_c=p_op_coeff%upwind_cell_blk(je,jk,jb)
            
            pflux_e(je,jk,jb)=pvn_e(je,jk,jb)&
              & *(pvar_c(il_c,jk,ib_c)&
              & +DOT_PRODUCT(p_op_coeff%moved_edge_position_cc(je,jk,jb)%x  &
              & -p_op_coeff%upwind_cell_position_cc(je,jk,jb)%x,&
              & z_gradc_cc(il_c,jk,ib_c)%x))
            
          ENDIF
        END DO
      END DO
    END DO
    
  END SUBROUTINE mimetic_miura_hflux_oce
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
  SUBROUTINE hflx_limiter_oce_mo( patch_3d,        &
    & p_cc,              &
    & p_mass_flx_e,      &
    & p_mflx_tracer_h,   &
    & p_op_coeff,        &
    & h_old,             &
    & h_new,             &
    & opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)           :: p_cc             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks) !< advected cell centered variable
    REAL(wp), INTENT(inout)           :: p_mass_flx_e     (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) !< horizontal mass flux
    REAL(wp), INTENT(inout)           :: p_mflx_tracer_h  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) !< calculated horizontal tracer mass flux
    TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
    REAL(wp), INTENT(in)              :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)              :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level !< optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level !< optional vertical end level
    
    
    !Local variables
    REAL(wp) :: z_mflx_low      (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)      !< first order tracer mass flux
    REAL(wp) :: z_anti          (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)          !< antidiffusive tracer mass flux (F_H - F_L)
    REAL(wp) :: z_mflx_anti     (3)
    REAL(wp) :: z_fluxdiv_c        !< flux divergence at cell center
    REAL(wp) :: z_tracer_new_low(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after hor. transport,
    !< if the low order fluxes are used
    REAL(wp) :: z_tracer_max(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    !< local maximum of current tracer value and low
    !< order update
    REAL(wp) :: z_tracer_min(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    !< local minimum of current tracer value and low
    !< order update
    REAL(wp) :: r_p(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)             !< fraction which must multiply all in/out fluxes
    !< of cell jc to guarantee
    REAL(wp) :: r_m(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)             !< no overshoot/undershoot
    REAL(wp) :: r_frac                                         !< computed minimum fraction which must multiply< the flux at the edge
    REAL(wp) :: z_min, z_max                         !< minimum/maximum value in cell and neighboring cells
    REAL(wp) :: z_signum                                       !< sign of antidiffusive velocity
    REAL(wp) :: p_p, p_m                                       !< sum of antidiffusive fluxes into and out of cell jc
    REAL(wp) :: prism_thick_old, inv_prism_thick_new
    INTEGER, DIMENSION(:,:,:), POINTER ::  cell_of_edge_idx, cell_of_edge_blk   !< Pointer to line and block indices of two
    !< neighbor cells (array)
    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk  !< Pointer to line and block indices of three
    !< neighbor cells (array)
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk    !< Pointer to line and block indices (array)
    !< of edges
    INTEGER :: start_level, end_level             !< vertical start and end level
    INTEGER :: start_index, end_index
    INTEGER :: je, jk, jb, jc,  cell_connect        !< index of edge, vert level, block, cell
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    
    ! CALL timer_start(timer_hflx_limiter)
    !-------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------
    ! Check for optional arguments
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
    
    ! Set pointers to index-arrays
    ! line and block indices of two neighboring cells
    cell_of_edge_idx => patch_2d%edges%cell_idx
    cell_of_edge_blk => patch_2d%edges%cell_blk
    ! line and block indices of edges as seen from cells
    edge_of_cell_idx => patch_2d%cells%edge_idx
    edge_of_cell_blk => patch_2d%cells%edge_blk
    ! pointers to line and block indices of three neighbor cells
    neighbor_cell_idx => patch_2d%cells%neighbor_idx
    neighbor_cell_blk => patch_2d%cells%neighbor_blk
    !
    ! 1. Calculate low (first) order fluxes using the standard upwind scheme and the
    !    antidiffusive fluxes
    !    (not allowed to call upwind_hflux_up directly, due to circular dependency)
    z_tracer_new_low(1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
    z_tracer_max    (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
    z_tracer_min    (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
    r_m             (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
    r_p             (1:nproma,1:n_zlev,1:patch_2d%alloc_cell_blocks) = 0.0_wp
    !    z_anti          (1:nproma,1:n_zlev,1:patch_2d%nblks_e) = 0.0_wp
    !    z_mflx_low      (1:nproma,1:n_zlev,1:patch_2d%nblks_e) = 0.0_wp
    
    
    !ICON_OMP_PARALLEL
    !ICON_OMP_DO PRIVATE(jb,start_index, end_index, je, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      DO je = start_index, end_index
        DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,jb), end_level)
          !IF( patch_3D%lsm_e(je,jk,jb) <= sea_boundary ) THEN
          z_mflx_low(je,jk,jb) = &
            & laxfr_upflux( p_mass_flx_e(je,jk,jb), &
            & p_cc(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)), &
            & p_cc(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)) )
          
          ! calculate antidiffusive flux for each edge
          z_anti(je,jk,jb) = p_mflx_tracer_h(je,jk,jb) - z_mflx_low(je,jk,jb)
          !ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
    !ICON_OMP_END_DO
    !ICON_OMP_END_PARALLEL
    
    ! no need to sync since we compute on cells_in_domain
    ! CALL sync_patch_array(SYNC_E, patch_2d,z_mflx_low)
    ! CALL sync_patch_array(SYNC_E, patch_2d, z_anti)
    
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        
        DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb), end_level)
          
          ! replace with merge
          IF (jk == start_level) THEN
            inv_prism_thick_new = 1.0_wp / (patch_3d%p_patch_1d(1)%del_zlev_m(start_level)+h_new(jc,jb))
            prism_thick_old     = patch_3d%p_patch_1d(1)%del_zlev_m(start_level) + h_old(jc,jb)
          ELSE
            prism_thick_old     = patch_3d%p_patch_1d(1)%del_zlev_m(jk)
            !            inv_prism_thick_new = 1.0_wp / patch_3D%p_patch_1D(1)%del_zlev_m(jk) ! should be calclulated only once
            inv_prism_thick_new = patch_3d%p_patch_1d(1)%inv_del_zlev_m(jk)
          ENDIF
          !
          
          ! IF( patch_3D%lsm_c(jc,jk,jb) > sea_boundary ) &
          !& CALL finish("","patch_3D%lsm_c(jc,jk,jb) > sea_boundar")
          
          !  compute also divergence of low order fluxes
          z_fluxdiv_c =  &
            & z_mflx_low(edge_of_cell_idx(jc,jb,1),jk,edge_of_cell_blk(jc,jb,1)) * &
            & p_op_coeff%div_coeff(jc,jk,jb,1) &
            & + z_mflx_low(edge_of_cell_idx(jc,jb,2),jk,edge_of_cell_blk(jc,jb,2)) * &
            & p_op_coeff%div_coeff(jc,jk,jb,2)  &
            & + z_mflx_low(edge_of_cell_idx(jc,jb,3),jk,edge_of_cell_blk(jc,jb,3)) * &
            & p_op_coeff%div_coeff(jc,jk,jb,3)
          
          !
          ! 3. Compute the updated low order solution z_tracer_new_low
          ! IF( patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
          z_tracer_new_low(jc,jk,jb) = (p_cc(jc,jk,jb) * prism_thick_old     &
            & - dtime * z_fluxdiv_c) * inv_prism_thick_new
          
          ! precalculate local maximum/minimum of current tracer value and low order
          ! updated value
          z_tracer_max(jc,jk,jb) = MAX(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
          z_tracer_min(jc,jk,jb) = MIN(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
          !          ELSE
          !            z_tracer_new_low(jc,jk,jb) = 0.0_wp
          !            z_tracer_max(jc,jk,jb)     = 0.0_wp
          !            z_tracer_min(jc,jk,jb)     = 0.0_wp
          !          ENDIF
          
        ENDDO
      ENDDO
    ENDDO
    ! !$OMP END DO
    ! !$OMP END PARALLEL
    
    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.
    CALL sync_patch_array(sync_c1, patch_2d,  z_tracer_max)
    CALL sync_patch_array(sync_c1, patch_2d,  z_tracer_min)
    
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        DO jk = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb), end_level)
          
          ! replace with merge
          IF (jk == start_level) THEN
            inv_prism_thick_new = 1.0_wp / (patch_3d%p_patch_1d(1)%del_zlev_m(start_level)+h_new(jc,jb))
            prism_thick_old     = patch_3d%p_patch_1d(1)%del_zlev_m(start_level) + h_old(jc,jb)
          ELSE
            prism_thick_old     = patch_3d%p_patch_1d(1)%del_zlev_m(jk)
            !            inv_prism_thick_new = 1.0_wp / patch_3D%p_patch_1D(1)%del_zlev_m(jk) ! should be calclulated only once
            inv_prism_thick_new = patch_3d%p_patch_1d(1)%inv_del_zlev_m(jk)
          ENDIF
          !
          ! 2. Define "antidiffusive" fluxes A(jc,jk,jb,je) for each cell. It is the difference
          !    between the high order fluxes (given by the FFSL-scheme) and the low order
          !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
          !    - positive for outgoing fluxes
          !    - negative for incoming fluxes
          !    this sign convention is related to the definition of the divergence operator.
          
          ! IF( patch_3D%lsm_c(jc,jk,jb) > sea_boundary ) &
          !& CALL finish("","patch_3D%lsm_c(jc,jk,jb) > sea_boundar")
          
          z_mflx_anti(1) =                                                        &
            & dtime * p_op_coeff%div_coeff(jc,jk,jb,1) * inv_prism_thick_new  &
            & * z_anti(edge_of_cell_idx(jc,jb,1),jk,edge_of_cell_blk(jc,jb,1))
          
          z_mflx_anti(2) =                                                         &
            & dtime *  p_op_coeff%div_coeff(jc,jk,jb,2) * inv_prism_thick_new  &
            & * z_anti(edge_of_cell_idx(jc,jb,2),jk,edge_of_cell_blk(jc,jb,2))
          
          z_mflx_anti(3) =                                                         &
            & dtime * p_op_coeff%div_coeff(jc,jk,jb,3) * inv_prism_thick_new   &
            & * z_anti(edge_of_cell_idx(jc,jb,3),jk,edge_of_cell_blk(jc,jb,3))
          
          IF( patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            !          ! max value of cell and its neighbors
            !          ! also look back to previous time step
            !          z_max = MAX( z_tracer_max(jc,jk,jb),                          &
            !            & z_tracer_max(neighbor_cell_idx(jc,jb,1),jk,neighbor_cell_blk(jc,jb,1)),  &
            !            & z_tracer_max(neighbor_cell_idx(jc,jb,2),jk,neighbor_cell_blk(jc,jb,2)),  &
            !            & z_tracer_max(neighbor_cell_idx(jc,jb,3),jk,neighbor_cell_blk(jc,jb,3)) )
            !          ! min value of cell and its neighbors
            !          ! also look back to previous time step
            !          z_min = MIN( z_tracer_min(jc,jk,jb),                          &
            !            & z_tracer_min(neighbor_cell_idx(jc,jb,1),jk,neighbor_cell_blk(jc,jb,1)),  &
            !            & z_tracer_min(neighbor_cell_idx(jc,jb,2),jk,neighbor_cell_blk(jc,jb,2)),  &
            !            & z_tracer_min(neighbor_cell_idx(jc,jb,3),jk,neighbor_cell_blk(jc,jb,3)) )
          ENDIF
          z_max = z_tracer_max(jc,jk,jb)
          z_min = z_tracer_min(jc,jk,jb)
          DO cell_connect = 1, 3
            IF (neighbor_cell_idx(jc,jb,cell_connect) > 0) THEN
              IF (patch_3d%p_patch_1d(1)% &
                & dolic_c(neighbor_cell_idx(jc,jb,cell_connect), neighbor_cell_blk(jc,jb,cell_connect)) >= jk) THEN
                z_max = MAX(z_max, &
                  & z_tracer_max(neighbor_cell_idx(jc,jb,cell_connect),jk,neighbor_cell_blk(jc,jb,cell_connect)))
                z_min = MIN(z_min, &
                  & z_tracer_min(neighbor_cell_idx(jc,jb,cell_connect),jk,neighbor_cell_blk(jc,jb,cell_connect)))
              ENDIF
            ENDIF
          ENDDO
          
          ! Sum of all incoming antidiffusive fluxes into cell jc
          p_p = - (MIN(0._wp,z_mflx_anti(1))   &
            & + MIN(0._wp,z_mflx_anti(2))   &
            & + MIN(0._wp,z_mflx_anti(3)) )
          ! Sum of all outgoing antidiffusive fluxes out of cell jc
          p_m =  MAX(0._wp,z_mflx_anti(1))  &
            & + MAX(0._wp,z_mflx_anti(2))  &
            & + MAX(0._wp,z_mflx_anti(3))
          
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
    
    ! Synchronize r_m and r_p
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, r_m, r_p)
    
    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !    At the end, compute new, limited fluxes which are then passed to the main
    !    program. Note that p_mflx_tracer_h now denotes the LIMITED flux.
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
          p_mflx_tracer_h(je,jk,jb) = z_mflx_low(je,jk,jb)               &
            & + MIN(1._wp,r_frac) * z_anti(je,jk,jb)
          ! ELSE
          !   p_mflx_tracer_h(je,jk,jb)= 0.0_wp
          ! ENDIF
        END DO
      ENDDO
    ENDDO
    
    ! CALL timer_stop(timer_hflx_limiter)
    
  END SUBROUTINE hflx_limiter_oce_mo
  !-------------------------------------------------------------------------
  
  
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
  ! ! INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, rl_start_c, rl_end_c
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
  ! !     CALL get_indices_c( patch_2d, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
  ! !     &                   rl_start_c, rl_end_c)
  ! !
  ! !     DO jc = i_startidx_c, i_endidx_c
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
  ! !     CALL get_indices_c( patch_2d, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
  ! !     &                   rl_start_c, rl_end_c)
  ! !
  ! !     DO jc = i_startidx_c, i_endidx_c
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
  ! !     CALL get_indices_c( patch_2d, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
  ! !     &                   rl_start_c, rl_end_c)
  ! !
  ! !     DO jc = i_startidx_c, i_endidx_c
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
  ! !       CALL get_indices_c( patch_2d, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
  ! !       &                   rl_start_c, rl_end_c)
  ! !
  ! !       DO jc = i_startidx_c, i_endidx_c
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
  ! ! !     CALL get_indices_c( patch_2d, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
  ! ! !     &                   rl_start_c, rl_end_c)
  ! ! !     DO jc = i_startidx_c, i_endidx_c
  ! ! ! write(123,*)'trac old new',trac_old(jc,1,jb),z_old_new(jc,1,jb),&
  ! ! ! &trac_new(jc,1,jb), z_out(jc,1,jb)
  ! ! !
  ! ! !    END DO
  ! ! !  END DO
  ! ! END SUBROUTINE elad
  ! !   !-------------------------------------------------------------------------
  
END MODULE mo_oce_tracer_transport_horz
