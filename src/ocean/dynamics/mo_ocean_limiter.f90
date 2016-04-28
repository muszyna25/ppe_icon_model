!>
!! Contains the implementation of the limter of the ocean model
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2016)
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
MODULE mo_ocean_limiter
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  USE mo_math_constants,            ONLY: dbl_eps
  USE mo_impl_constants,            ONLY: sea_boundary, SEA
  USE mo_ocean_nml,                 ONLY: n_zlev
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma, p_test_run
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer
  USE mo_timer,                     ONLY: timer_start, timer_stop, timers_level, timer_adv_horz, timer_hflx_lim, &
    & timer_dif_horz, timer_extra10, timer_extra11, timer_extra12, timer_extra13, timer_extra15
  USE mo_ocean_types,               ONLY: t_hydro_ocean_state, t_ocean_tracer
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish !, message_text, message
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff, no_primal_edges
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_c, sync_c1, sync_e, sync_patch_array, sync_patch_array_mult
  USE mo_mpi,                       ONLY: global_mpi_barrier
  
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=12)           :: str_module    = 'oceTracHorz '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  !
  ! PUBLIC INTERFACE
  !
  PUBLIC :: limiter_ocean_zalesak_horizontal
  PUBLIC :: limiter_ocean_posdef_horizontal 
  PUBLIC :: v_ppm_slimiter_mo
  PUBLIC :: v_ppm_slimiter_mo_onBlock
  PUBLIC :: vflx_limiter_pd_oce
  
  
  
  INTEGER, PARAMETER :: top=1
  
CONTAINS

  
  !-------------------------------------------------------------------------
  !>
  !! Flux limiter for horizontal advection
  !!
  !! Zalesak Flux-Limiter (Flux corrected transport)
  !! The corrected flux is a weighted average of the low order flux and the
  !! given high order flux. The high order flux is used to the greatest extent
  !! possible without introducing overshoots and undershoots.
  !! In vicinity of a lateral boundary only the low order flux is used: The criterion 
  !! is that at least one of the edges of the two neighboring cells of
  !! a central edges is a boundary edge.
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
  !! - adapted for ocean use by P. Korn (2012)
  !! - optimized by L. Linardakis (2015)
  !! - criterion for switch to low order scheme near boundaries added P. Korn (2015)  
  !!
  !!  mpi note: computed on domain edges. Results is not synced.
  !!
!<Optimize:inUse>
  SUBROUTINE limiter_ocean_zalesak_horizontal( patch_3d,&
    & vert_velocity,          &
    & tracer,                 &
    & p_mass_flx_e,           &
    & flx_tracer_low,         &    
    & flx_tracer_high,        &
    & flx_tracer_final,       &
    & div_adv_flux_vert,      &   
    & operators_coefficients, &
    & h_old,                  &
    & h_new, zlim, tracer_index )                  
    
    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
    REAL(wp),INTENT(inout)              :: vert_velocity(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp), INTENT(inout)             :: tracer           (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: p_mass_flx_e     (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_low   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: flx_tracer_high  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp), INTENT(inout)             :: flx_tracer_final (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    TYPE(t_operator_coeff),INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(in)                :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: zlim(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    INTEGER, INTENT(in)                 :: tracer_index

    
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
    REAL(wp) :: z_tracer_update_horz(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after transport, if low order fluxes are used    
    REAL(wp) :: r_frac          !< computed minimum fraction which must multiply< the flux at the edge
    REAL(wp) :: z_min, z_max    !< minimum/maximum value in cell and neighboring cells
    REAL(wp) :: z_signum        !< sign of antidiffusive velocity
    REAL(wp) :: p_p, p_m        !< sum of antidiffusive fluxes into and out of cell jc
    REAL(wp) :: prism_thick_old(n_zlev), inv_prism_thick_new(n_zlev)
    REAL(wp) :: delta_z_new, delta_z
    INTEGER, DIMENSION(:,:,:), POINTER ::  cellOfEdge_idx, cellOfEdge_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    INTEGER :: start_level, end_level            
    INTEGER :: start_index, end_index
!     INTEGER :: cell_1_idx,cell_2_idx,cell_1_blk,cell_2_blk
!     INTEGER :: cell_1_edge_1_idx,cell_1_edge_2_idx,cell_1_edge_3_idx
!     INTEGER :: cell_2_edge_1_idx,cell_2_edge_2_idx,cell_2_edge_3_idx        
!     INTEGER :: cell_1_edge_1_blk,cell_1_edge_2_blk,cell_1_edge_3_blk
!     INTEGER :: cell_2_edge_1_blk,cell_2_edge_2_blk,cell_2_edge_3_blk      
    INTEGER :: edge_index, level, blockNo, jc,  cell_connect, sum_lsm_quad_edge, ctr
!    INTEGER :: all_water_edges 
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d    
    !-------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------
    start_level = 1
    end_level   = n_zlev
    cellOfEdge_idx  => patch_2d%edges%cell_idx
    cellOfEdge_blk  => patch_2d%edges%cell_blk
    edge_of_cell_idx  => patch_2d%cells%edge_idx
    edge_of_cell_blk  => patch_2d%cells%edge_blk
    neighbor_cell_idx => patch_2d%cells%neighbor_idx
    neighbor_cell_blk => patch_2d%cells%neighbor_blk
    
            

#ifdef NAGFOR
    z_tracer_max(:,:,:) = 0.0_wp
    z_tracer_min(:,:,:) = 0.0_wp
    r_m(:,:,:)          = 0.0_wp
    r_p(:,:,:)          = 0.0_wp
#endif

   
!ICON_OMP_PARALLEL
! !ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE        
!       DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block      
!         flux_div_vert(:,:,blockNo) = 0.0_wp     
!         CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)            
!         DO jc = start_index, end_index    
!           DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)        
!             ! positive vertical divergence in direction of w (upward positive)
!             flux_div_vert(jc,level,blockNo) = z_adv_flux_v(jc, level, blockNo) &
!             & - z_adv_flux_v(jc, level+1, blockNo)
!           ENDDO
!         END DO
!       END DO
! !ICON_OMP_END_DO      
! !       CALL sync_patch_array(sync_c, patch_2D, flux_div_vert)
 

!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      
      z_anti(:,:,blockNo)     = 0.0_wp
      DO edge_index = start_index, end_index
        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)
          
          ! calculate antidiffusive flux for each edge
          z_anti(edge_index,level,blockNo) = flx_tracer_high(edge_index,level,blockNo)&
                                          &- flx_tracer_low(edge_index,level,blockNo)
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
!ICON_OMP_END_DO

    
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, delta_z, delta_z_new, &
!ICON_OMP z_fluxdiv_c, div_adv_flux_vert) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      
      z_tracer_new_low(:,:,blockNo)    = 0.0_wp
      z_tracer_update_horz(:,:,blockNo)= 0.0_wp           
      DO jc = start_index, end_index
        
        ! get prism thickness
        !inv_prism_thick_new(start_level) = 1.0_wp / (patch_3d%p_patch_1d(1)%del_zlev_m(start_level) + h_new(jc,blockNo))
        !prism_thick_old(start_level)     = patch_3d%p_patch_1d(1)%del_zlev_m(start_level)           + h_old(jc,blockNo)
        !DO level = start_level+1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)
        !  prism_thick_old (level)    = patch_3d%p_patch_1d(1)%del_zlev_m(level)
        !  inv_prism_thick_new(level) = patch_3d%p_patch_1d(1)%inv_del_zlev_m(level)
        !ENDDO
        
        ! 3. Compute the complete (with horizontal and vertical divergence) updated low order solution z_tracer_new_low
        ! First at top level than in fluid interior
        !       
        level = start_level
        !  compute divergence of low order fluxes
        z_fluxdiv_c = 0
        DO cell_connect = 1, patch_2d%cells%num_edges(jc,blockNo)
          z_fluxdiv_c =  z_fluxdiv_c + &
            & flx_tracer_low(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect)) * &
            & operators_coefficients%div_coeff(jc,level,blockNo,cell_connect)
        ENDDO
       
        delta_z = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)&
             &  + h_old(jc,blockNo)
        delta_z_new = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)&
             &  + h_new(jc,blockNo)
             
        ! Low order flux at top level
        !z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
        !  & - dtime * (z_fluxdiv_c+flux_div_vert))/delta_z_new
        z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
          & - dtime * (z_fluxdiv_c+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new

         z_tracer_update_horz(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
         & - dtime * (div_adv_flux_vert(jc,level,blockNo)))/delta_z_new

        !write(0,*) "z_tracer_new_low is done"
        !
        !Fluid interior       
        DO level = start_level+1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)       
          !  compute divergence of low order fluxes
          z_fluxdiv_c = 0
          DO cell_connect = 1, patch_2d%cells%num_edges(jc,blockNo)
            z_fluxdiv_c =  z_fluxdiv_c + &
              & flx_tracer_low(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect)) * &
              & operators_coefficients%div_coeff(jc,level,blockNo,cell_connect)
          ENDDO


          delta_z     = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)
          delta_z_new = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)
          !
          ! low order flux in flow interior
          !z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * prism_thick_old(level)     &
          !  & - dtime * (z_fluxdiv_c+flux_div_vert(jc,level,blockNo))) * inv_prism_thick_new(level)
          ! z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
          !   & - dtime * (z_fluxdiv_c+flux_div_vert))/delta_z_new
          !z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
          !    & - dtime * (z_fluxdiv_c+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new
              
          z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
            & - dtime * (z_fluxdiv_c+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new
            
          !z_tracer_update_horz(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
          !  & - dtime * (div_adv_flux_vert(jc,level,blockNo)))/delta_z_new
        ENDDO
      ENDDO
      
      ! precalculate local maximum/minimum of current tracer value and low order
      ! updated value
      z_tracer_max(:,:,blockNo) =            &
        & MAX(          tracer(:,:,blockNo), &
        &     z_tracer_new_low(:,:,blockNo))
      z_tracer_min(:,:,blockNo) =            &
        & MIN(          tracer(:,:,blockNo), &
        &     z_tracer_new_low(:,:,blockNo))

!      write(0,*) blockNo, ":", z_tracer_max(start_index:end_index,start_level:end_level,blockNo)
!      write(0,*) blockNo, ":", z_tracer_min(start_index:end_index,start_level:end_level,blockNo)
    ENDDO
!ICON_OMP_END_DO

! DO level = 1, 4
! write(0,*)'tracer bounds',level, maxval(z_tracer_max(:,level,:)) ,minval(z_tracer_max(:,level,:)),&
! &maxval(z_tracer_min(:,level,:)) ,minval(z_tracer_min(:,level,:))
! END DO

!ICON_OMP_MASTER
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, z_tracer_max, z_tracer_min)
!ICON_OMP_END_MASTER    
!ICON_OMP_BARRIER
    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.    
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, inv_prism_thick_new, &
!ICON_OMP z_mflx_anti, z_max, z_min, cell_connect, p_p, p_m) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block

      ! this is only needed for the parallel test setups
      ! it will try  tocheck the uninitialized (land) parts
      r_m(:,:,blockNo) = 0.0_wp
      r_p(:,:,blockNo) = 0.0_wp
        
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        
        ! get prism thickness
        inv_prism_thick_new(start_level) = 1.0_wp / (patch_3d%p_patch_1d(1)%del_zlev_m(start_level)+ h_new(jc,blockNo))
        DO level = start_level+1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)
          inv_prism_thick_new(level) = patch_3D%p_patch_1d(1)%inv_prism_thick_c(jc,level,blockNo)
        ENDDO
        
        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)         
          ! 2. Define "antidiffusive" fluxes A(jc,level,blockNo,edge_index) for each cell. It is the difference
          !    between the high order fluxes (given by the FFSL-scheme) and the low order
          !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
          !    - positive for outgoing fluxes
          !    - negative for incoming fluxes
          !    this sign convention is related to the definition of the divergence operator.          
          z_mflx_anti(:) = 0.0_wp
          z_max = z_tracer_max(jc,level,blockNo)
          z_min = z_tracer_min(jc,level,blockNo)
          p_p = 0.0_wp
          p_m = 0_wp
          DO cell_connect = 1, patch_2d%cells%num_edges(jc,blockNo)
            IF (patch_3d%p_patch_1d(1)% &
              & dolic_c(neighbor_cell_idx(jc,blockNo,cell_connect), neighbor_cell_blk(jc,blockNo,cell_connect)) >= level) THEN
              
              z_max = MAX(z_max, &
                & z_tracer_max(neighbor_cell_idx(jc,blockNo,cell_connect),level,neighbor_cell_blk(jc,blockNo,cell_connect)))
              z_min = MIN(z_min, &
                & z_tracer_min(neighbor_cell_idx(jc,blockNo,cell_connect),level,neighbor_cell_blk(jc,blockNo,cell_connect)))

              z_mflx_anti(cell_connect) =                                                        &
                & dtime * operators_coefficients%div_coeff(jc,level,blockNo,cell_connect) * inv_prism_thick_new(level)  &
                & * z_anti(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect))

              ! Sum of all incoming antidiffusive fluxes into cell jc
              ! outgoing fluxes carry a positive sign, incoming a negative
              p_p = p_p - MIN(0._wp, z_mflx_anti(cell_connect))
              ! Sum of all outgoing antidiffusive fluxes out of cell jc
              p_m = p_m + MAX(0._wp, z_mflx_anti(cell_connect))
            ENDIF
          ENDDO                
          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of tracer
          r_m(jc,level,blockNo) = (z_tracer_new_low(jc,level,blockNo) - z_min ) / (p_m + dbl_eps)!&
          !
          ! fraction which must multiply all fluxes into cell jc to guarantee no
          ! overshoot
          ! Nominator: maximum allowable increase of tracer
          r_p(jc,level,blockNo) = (z_max - z_tracer_new_low(jc,level,blockNo)) / (p_p + dbl_eps)!&
          !
          !update old tracer with low-order flux
          !!tracer(jc,level,blockNo)=z_tracer_new_low(jc,level,blockNo) 
          !tracer(jc,level,blockNo)=z_tracer_update_horz(jc,level,blockNo)         
        ENDDO
      ENDDO
    ENDDO
!ICON_OMP_END_DO

! DO level = 1, 10
! write(0,*)'updated:trac',level, maxval(tracer(:,level,:)) ,minval(tracer(:,level,:)),&
! &maxval(z_tracer_new_low(:,level,:)) ,minval(z_tracer_new_low(:,level,:))
! END DO

!ctr=0

    
!ICON_OMP_MASTER
    ! Synchronize r_m and r_p
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, r_m, r_p)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER   

    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !    At the end, compute new, limited fluxes which are then passed to the main
!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level, z_signum, r_frac) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      flx_tracer_final(:,:,blockNo) = 0.0_wp
      DO edge_index = start_index, end_index
      
        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)
        
          IF( operators_coefficients%edges_SeaBoundaryLevel(edge_index,level,blockNo) > -2)THEN! edge < 2nd order boundary
          
            flx_tracer_final(edge_index,level,blockNo) = flx_tracer_low(edge_index,level,blockNo)
            
          ELSE!IF(sum_lsm_quad_edge==all_water_edges)THEN
          
            !z_anti>0 returns  1: here z_anti is outgoing, i.e. flux_high>flux_low
            !z_anti<0 returns -1: here z_anti is ingoing, i.e. flux_high<flux_low
            z_signum = SIGN(1._wp, z_anti(edge_index,level,blockNo))
                    
          ! This does the same as an IF (z_signum > 0) THEN ... ELSE ... ENDIF,
          ! but is computationally more efficient
          r_frac = 0.5_wp * (       &
            & (1._wp + z_signum) * & !<- active for z_signum=1
            & MIN(r_m(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)),  &
            &     r_p(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)))  &
            &+(1._wp - z_signum) * & !<- active for z_signum=-1
            & MIN(r_m(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)),  &
            &     r_p(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)))  )
          
          ! Limited flux
          flx_tracer_final(edge_index,level,blockNo) = flx_tracer_low(edge_index,level,blockNo)&
           & + MIN(1.0_wp,r_frac) *z_anti(edge_index,level,blockNo)      

            ENDIF
        END DO
       ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL


  END SUBROUTINE limiter_ocean_zalesak_horizontal
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
  SUBROUTINE limiter_ocean_posdef_horizontal( patch_3d,&
    & tracer,                             &
    & flx_tracer,                         &
    & flx_tracer_limit,                   &
    & operators_coefficients )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
    REAL(wp), INTENT(inout)             :: tracer      (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: flx_tracer(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp), INTENT(inout)             :: flx_tracer_limit(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    TYPE(t_operator_coeff),INTENT(in)   :: operators_coefficients

    !Local variables
    REAL(wp) :: z_mflx(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)   
    REAL(wp) :: r_m             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< no overshoot/undershoot
    REAL(wp) :: z_signum        !< sign of antidiffusive velocity
    REAL(wp) :: p_m        !< sum of antidiffusive fluxes into and out of cell jc
    INTEGER, DIMENSION(:,:,:), POINTER :: cellOfEdge_idx, cellOfEdge_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    INTEGER :: start_level, end_level            
    INTEGER :: start_index, end_index
    INTEGER :: edge_index, level, blockNo, jc!,  cell_connect
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d    
    REAL(wp) :: z_adv_flux_v(nproma, n_zlev+1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)  ! upwind flux
    !REAL(wp) :: flux_div_vert(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !new tracer        
    !-------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------
    start_level = 1
    end_level   = n_zlev
    
    flx_tracer_limit(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%nblks_e)          = 0.0_wp
    z_mflx          (1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)= 0.0_wp 
    r_m             (1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)= 0.0_wp  
    
    ! Set pointers to index-arrays
    ! line and block indices of two neighboring cells
    cellOfEdge_idx  => patch_2d%edges%cell_idx
    cellOfEdge_blk  => patch_2d%edges%cell_blk
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
    !!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, inv_prism_thick_new, prism_thick_old, &
    !!ICON_OMP z_fluxdiv_c ) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      
      DO jc = start_index, end_index
                
        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)
        
          z_mflx(jc,level,1) = dtime*operators_coefficients%div_coeff(jc,level,blockNo,1) &
            & * flx_tracer(edge_of_cell_idx(jc,blockNo,1),level,edge_of_cell_blk(jc,blockNo,1))
            

          z_mflx(jc,level,2) = dtime*operators_coefficients%div_coeff(jc,level,blockNo,2)&
            & * flx_tracer(edge_of_cell_idx(jc,blockNo,2),level,edge_of_cell_blk(jc,blockNo,2))

  
          z_mflx(jc,level,3) = dtime*operators_coefficients%div_coeff(jc,level,blockNo,3)&
            & * flx_tracer(edge_of_cell_idx(jc,blockNo,3),level,edge_of_cell_blk(jc,blockNo,3))

          
        ENDDO
      ENDDO
    ENDDO
    !!ICON_OMP_END_DO NOWAIT
    !!ICON_OMP_END_PARALLEL
    CALL sync_patch_array(SYNC_C1, patch_2d, z_mflx)
    ! 2. Compute total outward mass
    !
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)

          ! Sum of all outgoing fluxes out of cell jc
          p_m =  MAX(0._wp,z_mflx(jc,level,1))  &
            &  + MAX(0._wp,z_mflx(jc,level,2))  &
            &  + MAX(0._wp,z_mflx(jc,level,3))

          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of \rho q
          r_m(jc,level,blockNo) = MIN(1._wp, &
            & (tracer(jc,level,blockNo) &
            &  * patch_3d%p_patch_1d(1)%prism_thick_c(jc,level,blockNo)) &
            &  /(p_m + dbl_eps) )

        ENDDO
      ENDDO
    ENDDO
   ! synchronize r_m
    CALL sync_patch_array(SYNC_C1, patch_2d, r_m)

    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      DO edge_index = start_index, end_index
        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)

          ! p_mflx_tracer_h > 0: flux directed from cell 1 -> 2
          ! p_mflx_tracer_h < 0: flux directed from cell 2 -> 1
          z_signum = SIGN(1._wp,flx_tracer(edge_index,level,blockNo))
          flx_tracer_limit(edge_index,level,blockNo) = flx_tracer(edge_index,level,blockNo) * 0.5_wp  &
            & * ( (1._wp + z_signum) &
            &      * r_m(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)) &
            &    + (1._wp - z_signum) &
            &      * r_m(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)) )

        END DO
      ENDDO
    ENDDO
    CALL sync_patch_array(sync_e, patch_2d, flx_tracer_limit)
  !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    !CALL dbg_print('LimitPosDef: scalefac'   ,r_m,str_module,idt_src,patch_2d%cells%owned)    
    CALL dbg_print('LimitPosDef: flux_h'     ,flx_tracer_limit,str_module,idt_src,patch_2d%edges%owned)
    !---------------------------------------------------------------------
            
  END SUBROUTINE limiter_ocean_posdef_horizontal
  !-------------------------------------------------------------------------
 
 
   !-------------------------------------------------------------------------
  !>
  !! Limiter for PPM (3rd order) vertical advection (monotone version)
  !!
  !! Removes over- and undershoots in first guess parabola by resetting the
  !! upper or lower interface values.
  !! Avoids non-physical over/undershoots in advected fields.
  !!
  !! Note that this limiter was coded assuming a pressure based vertical
  !! coordinate system. Nevertheless this limiter works for a height based
  !! vertical system, too. This is due to a 'wrong' computation of z_delta
  !! in the case of a height based coordinate system (i.e. z_delta is
  !! implicity multiplied by -1)
  !!
  !! Literature
  !! Lin and Rood (1996), MWR, 124, 2046-2070
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-02-04)
  !!
  SUBROUTINE v_ppm_slimiter_mo( patch_3d, p_cc, p_face, p_slope, p_face_up, p_face_low )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)           :: p_cc(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)      !< advected cell centered variable
    REAL(wp), INTENT(inout)           :: p_face(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)  !< reconstructed face values of the advected field
    REAL(wp), INTENT(inout)           :: p_slope(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks) !< monotonized slope
    REAL(wp), INTENT(inout)           :: p_face_up(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks) !< final face value (upper face, height based)
    REAL(wp), INTENT(inout)           :: p_face_low(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< final face value (lower face, height based)
    
    ! locals
    INTEGER :: nlev                      !< number of full levels
    INTEGER :: firstLevel                      !< vertical start thisLevel
    INTEGER :: jc, jk, jb                !< index of cell, vertical thisLevel and block
    INTEGER :: startIndex, endIndex
    INTEGER :: ikp1                      !< vertical thisLevel plus one
    INTEGER :: i_dolic
    REAL(wp) :: z_delta                   !< lower minus upper face value
    REAL(wp) :: z_a6i                     !< curvature of parabola
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-----------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    
    firstLevel = 1
    nlev = n_zlev
    ! !$OMP PARALLEL
    ! !$OMP DO PRIVATE(jb,jk,jc,startIndex,endIndex,ikp1,z_delta,z_a6i)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      DO jc = startIndex, endIndex
        
        i_dolic=patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
        
        DO jk = firstLevel, i_dolic
          ! index of bottom half thisLevel
          ikp1 = jk + 1
          
          IF ( patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            z_delta   = p_face(jc,ikp1,jb) - p_face(jc,jk,jb)
            z_a6i     = 6._wp * (p_cc(jc,jk,jb)                           &
              & - 0.5_wp * (p_face(jc,jk,jb) + p_face(jc,ikp1,jb)))
            
            IF ( p_slope(jc,jk,jb) == 0._wp) THEN
              p_face_up(jc,jk,jb)  = p_cc(jc,jk,jb)
              p_face_low(jc,jk,jb) = p_cc(jc,jk,jb)
              
            ELSE IF (z_delta * z_a6i > z_delta * z_delta) THEN
              p_face_up(jc,jk,jb)  = 3._wp*p_cc(jc,jk,jb) - 2._wp*p_face(jc,ikp1,jb)
              p_face_low(jc,jk,jb) = p_face(jc,ikp1,jb)
              
            ELSE IF (z_delta * z_a6i < -1._wp * (z_delta * z_delta)) THEN
              p_face_up(jc,jk,jb)  = p_face(jc,jk,jb)
              p_face_low(jc,jk,jb) = 3._wp*p_cc(jc,jk,jb) - 2._wp*p_face(jc,jk,jb)
              
            ELSE
              p_face_up(jc,jk,jb)  = p_face(jc,jk,jb)
              p_face_low(jc,jk,jb) = p_face(jc,ikp1,jb)
            ENDIF
          ENDIF
        END DO
      END DO
    END DO
    ! !$OMP END DO
    ! !$OMP END PARALLEL
    
  END SUBROUTINE v_ppm_slimiter_mo
  !-------------------------------------------------------------------------
  !>
  !! Limiter for PPM (3rd order) vertical advection (monotone version)
  !!
  !! Removes over- and undershoots in first guess parabola by resetting the
  !! upper or lower interface values.
  !! Avoids non-physical over/undershoots in advected fields.
  !!
  !! Note that this limiter was coded assuming a pressure based vertical
  !! coordinate system. Nevertheless this limiter works for a height based
  !! vertical system, too. This is due to a 'wrong' computation of z_delta
  !! in the case of a height based coordinate system (i.e. z_delta is
  !! implicity multiplied by -1)
  !!
  !! Literature
  !! Lin and Rood (1996), MWR, 124, 2046-2070
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-02-04)
  !!
  !! mpi parallelized, only cells_in_domain are computed, no sync
!<Optimize:inUse>
  SUBROUTINE v_ppm_slimiter_mo_onBlock( p_cc, p_face, p_slope, p_face_up, p_face_low, &
    & startIndex, endIndex, cells_noOfLevels )

    REAL(wp), INTENT(in)           :: p_cc(nproma,n_zlev)      !< advected cell centered variable
    REAL(wp), INTENT(in)           :: p_face(nproma,n_zlev+1)  !< reconstructed face values of the advected field
    REAL(wp), INTENT(in)           :: p_slope(nproma,n_zlev+1) !< monotonized slope
    REAL(wp), INTENT(inout)        :: p_face_up(nproma,n_zlev) !< final face value (upper face, height based)
    REAL(wp), INTENT(inout)        :: p_face_low(nproma,n_zlev)!< final face value (lower face, height based)
    INTEGER,  INTENT(in)           :: startIndex, endIndex
    INTEGER,  INTENT(in)           :: cells_noOfLevels(nproma)

    ! locals
    INTEGER :: nlev                      !< number of full levels
    INTEGER :: firstLevel                      !< vertical start thisLevel
    INTEGER :: jc, jk                   !< index of cell, vertical thisLevel
    INTEGER :: ikp1                      !< vertical thisLevel plus one
    REAL(wp) :: z_delta                   !< lower minus upper face value
    REAL(wp) :: z_a6i                     !< curvature of parabola
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------

    firstLevel = 1
    nlev = n_zlev

! !CDIR NODEP
    DO jc = startIndex, endIndex

! !CDIR NODEP
      DO jk = firstLevel, cells_noOfLevels(jc)
        ! index of bottom half thisLevel
        ikp1 = jk + 1

        z_delta   = p_face(jc,ikp1) - p_face(jc,jk)
        z_a6i     = 6._wp * (p_cc(jc,jk)                           &
          & - 0.5_wp * (p_face(jc,jk) + p_face(jc,ikp1)))

        IF ( p_slope(jc,jk) == 0._wp) THEN
          p_face_up(jc,jk)  = p_cc(jc,jk)
          p_face_low(jc,jk) = p_cc(jc,jk)

        ELSE IF (z_delta * z_a6i > z_delta * z_delta) THEN
          p_face_up(jc,jk)  = 3._wp*p_cc(jc,jk) - 2._wp*p_face(jc,ikp1)
          p_face_low(jc,jk) = p_face(jc,ikp1)

        ELSE IF (z_delta * z_a6i < -1._wp * (z_delta * z_delta)) THEN
          p_face_up(jc,jk)  = p_face(jc,jk)
          p_face_low(jc,jk) = 3._wp*p_cc(jc,jk) - 2._wp*p_face(jc,jk)

        ELSE
          p_face_up(jc,jk)  = p_face(jc,jk)
          p_face_low(jc,jk) = p_face(jc,ikp1)
        ENDIF

      END DO
        
    END DO

  END SUBROUTINE v_ppm_slimiter_mo_onBlock
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Positive definite flux limiter for vertical advection
  !!
  !! Positive definite Zalesak Flux-Limiter (Flux corrected transport).
  !! for the hydrostatic core. Only outward fluxes are re-scaled, in
  !! order to maintain positive definiteness.
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !! - Harris, L. M. and P. H. Lauritzen (2010): A flux-form version of the
  !!   Conservative Semi-Lagrangian Multi-tracer transport scheme (CSLAM) on
  !!   the cubed sphere grid.  J. Comput. Phys., 230, 1215-1237
  !! - Smolarkiewicz, P. K., 1989: Comment on "A positive definite advection
  !!   scheme obtained by nonlinear renormalization of the advective fluxes.",
  !!   Mon. Wea. Rev., 117, 2626-2632
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2011-01-07)
  !!
  SUBROUTINE vflx_limiter_pd_oce( patch_3d, p_dtime, p_cc, p_cellhgt_mc_now, p_flx_tracer_v, &
    & opt_slev, opt_elev )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)          :: p_cc(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)           !< advected cell centered variable at time (n)
    REAL(wp), INTENT(inout)          :: p_cellhgt_mc_now(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)             :: p_dtime
    REAL(wp), INTENT(inout)          :: p_flx_tracer_v(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks) !< calculated vertical tracer mass flux
    INTEGER,  INTENT(in), OPTIONAL :: opt_slev
    INTEGER,  INTENT(in), OPTIONAL :: opt_elev
    !
    !Local variables
    REAL(wp) :: r_m(nproma,n_zlev)      !< fraction which must multiply all
    !< outgoing fluxes of cell jc
    !< to guarantee positive definiteness
    REAL(wp) :: p_m(nproma)            !< sum of fluxes out of cell
    REAL(wp) :: z_signum(nproma)       !< sign of mass flux       !< >0: upward; <0: downward
    
    INTEGER :: firstLevel, elev             !< vertical start and end thisLevel
    INTEGER :: startIndex, endIndex
    INTEGER :: jk, jb, jc            !< index of edge, vert thisLevel, block, cell
    INTEGER :: jkp1, jkm1
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-----------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    ! Do jk=1,n_zlev
    ! write(0,*)'profile before limit',jk,p_flx_tracer_v(5,jk,5)
    ! END DO
    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      firstLevel = opt_slev
    ELSE
      firstLevel = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF
    
    r_m(1:nproma,1:n_zlev)= 0.0_wp
    p_m(1:nproma)         = 0.0_wp
    z_signum(1:nproma)    = 0.0_wp
    
    
    ! !$OMP PARALLEL
    ! !$OMP DO PRIVATE(jb,jk,jc,startIndex,endIndex,jkp1,p_m,r_m,jkm1,z_signum)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, startIndex, endIndex)
      !
      ! 1. Compute total outward mass
      !
      DO jk = firstLevel, elev
        jkp1 = jk+1
        
        DO jc = startIndex, endIndex
          
          IF ( patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            ! Sum of all outgoing fluxes out of cell jk
            p_m(jc) = p_dtime                                &
              & * (MAX(0._wp,p_flx_tracer_v(jc,jk,jb))  &  ! upper half thisLevel
              & - MIN(0._wp,p_flx_tracer_v(jc,jkp1,jb)) )     ! lower half thisLevel
            ! fraction which must multiply the fluxes out of cell jk to guarantee no
            ! undershoot
            ! Nominator: maximum allowable decrease \rho^n q^n
            r_m(jc,jk) = MIN(1._wp, (p_cc(jc,jk,jb)*p_cellhgt_mc_now(jc,jk,jb)) &
              & /(p_m(jc) + dbl_eps) )
          ENDIF
        ENDDO
        
      ENDDO
      !
      ! 2. Limit outward fluxes (loop over half levels)
      !    Choose r_m depending on the sign of p_mflx_tracer_v
      !
      DO jk = firstLevel+1, elev
        jkm1 = jk-1
        DO jc = startIndex, endIndex
          
          IF ( patch_3d%lsm_c(jc,jk,jb) <= sea_boundary ) THEN
            ! p_mflx_tracer_v(k-1/2) > 0: flux directed from cell k   -> k-1
            ! p_mflx_tracer_v(k-1/2) < 0: flux directed from cell k-1 -> k
            z_signum(jc) = SIGN(1._wp,p_flx_tracer_v(jc,jk,jb))
            
            p_flx_tracer_v(jc,jk,jb) =  p_flx_tracer_v(jc,jk,jb)  * 0.5_wp    &
              & * ( (1._wp + z_signum(jc)) * r_m(jc,jk) &
              & +   (1._wp - z_signum(jc)) * r_m(jc,jkm1) )
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    
    ! !$OMP END DO NOWAIT
    ! !$OMP END PARALLEL
  END SUBROUTINE vflx_limiter_pd_oce
  
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
  
  !------------------------------------------------------------------------
 
 
 
  
!  !-------------------------------------------------------------------------
!  !>
!  !! Lax Friedrichs first order upwind flux,
!  !! used in conservative advection routines.
!  !! For passive advection, equivalent to
!  !! any other first order upwind flux.
!  !!
!  !! @par Revision History
!  !! Developed  by L.Bonaventura  (2004).
!  !!
!  FUNCTION laxfr_upflux( p_vn, p_psi1, p_psi2 )  result(p_upflux)
!    !
!    IMPLICIT NONE
!    REAL(wp), INTENT(in) :: p_vn
!    REAL(wp), INTENT(in) :: p_psi1, p_psi2
!    REAL(wp)             :: p_upflux
!    !-----------------------------------------------------------------------
!    p_upflux = 0.5_wp * (        p_vn  *( p_psi1 + p_psi2 )    &
!      & - ABS( p_vn )*( p_psi2 - p_psi1 ) )
!    
!  END FUNCTION laxfr_upflux
!  
! !   !-------------------------------------------------------------------------
END MODULE mo_ocean_limiter
