!>
!! Testbed for modifications requiresd to enable zstar 
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_ocean_testbed_zstar
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: max_char_length, sea_boundary, zero_coriolis
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d,t_subset_range
  USE mo_grid_config,            ONLY: n_dom, grid_sphere_radius, grid_angular_velocity
  USE mo_math_constants,         ONLY: pi, pi_2, rad2deg, deg2rad, dbl_eps
  USE mo_math_types,             ONLY: t_cartesian_coordinates
  USE mo_ocean_nml,              ONLY: n_zlev, GMRedi_configuration, GMRedi_combined, Cartesian_Mixing, &
    & atmos_flux_analytical_type, no_tracer, OceanReferenceDensity, l_with_vert_tracer_advection, &
    & tracer_update_mode, use_none, l_edge_based
  USE mo_sea_ice_nml,            ONLY: init_analytic_conc_param, t_heat_base
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_run_config,             ONLY: nsteps, dtime, output_mode, test_mode !, test_param
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ext_data_types,         ONLY: t_external_data
  !USE mo_io_units,               ONLY: filename_max
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_total
  USE mo_ocean_ab_timestepping,  ONLY: update_time_indices , &
     & solve_free_surface_eq_ab, calc_vert_velocity
  USE mo_random_util,            ONLY: add_random_noise_global
  USE mo_ocean_types,            ONLY: t_hydro_ocean_state, t_operator_coeff, t_solvercoeff_singleprecision
  USE mo_hamocc_types,           ONLY: t_hamocc_state
  USE mo_restart,                ONLY: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor
  USE mo_restart_attributes,     ONLY: t_RestartAttributeList, getAttributesForRestarting
  USE mo_io_config,              ONLY: n_checkpoints, write_last_restart
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff, no_primal_edges
  USE mo_ocean_tracer,           ONLY: advect_ocean_tracers
  USE mo_ocean_surface_refactor, ONLY: update_ocean_surface_refactor
  USE mo_ocean_surface_types,    ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_sea_ice,                ONLY: salt_content_in_surface, energy_content_in_surface
  USE mo_sea_ice_types,          ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_ice_diagnostics,        ONLY: energy_in_surface, salt_in_surface
  USE mo_physical_constants,     ONLY: rhoi, rhos, clw, alf, Tf
  USE mo_ocean_physics_types,    ONLY: t_ho_params
  USE mo_master_config,          ONLY: isRestart
  USE mo_ocean_GM_Redi,          ONLY: prepare_ocean_physics,calc_ocean_physics
  USE mo_ocean_diagnostics,      ONLY: calc_fast_oce_diagnostics, calc_psi
  USE mo_ocean_thermodyn,        ONLY: calc_potential_density, calculate_density,&
  &                                    calc_neutralslope_coeff_func_onColumn,calc_neutralslope_coeff_func_onColumn_UNESCO
  USE mo_time_config,            ONLY: time_config
  USE mo_statistics
  USE mo_util_dbg_prnt,          ONLY: dbg_print
  USE mo_ocean_statistics
  USE mo_ocean_output
  USE mo_parallel_config,        ONLY: nproma
  USE mo_statistics
  USE mo_ocean_testbed_vertical_diffusion
  USE mo_ocean_math_operators
  USE mo_grid_subset,            ONLY: t_subset_range, get_index_range 
  USE mo_scalar_product,         ONLY: calc_scalar_product_veloc_3d, &
      & map_edges2cell_3d, map_scalar_center2prismtop, map_vec_prismtop2center_on_block, &
      & map_cell2edges_3D, map_edges2edges_viacell_3d_const_z
  USE mo_ocean_tracer_transport_horz, ONLY: diffuse_horz
  USE mo_hydro_ocean_run
  USE mo_var_list
  USE mo_linked_list
  USE mo_cdi
  use mo_cdi_constants
  use mo_zaxis_type
  use mo_cf_convention
  use mo_grib2

  USE mtime,                     ONLY: datetime, newDatetime, deallocateDatetime, datetimeToString, &
       &                               timedelta, newTimedelta, deallocateTimedelta,                &
       &                               MAX_DATETIME_STR_LEN, newDatetime,                           &
       &                               MAX_MTIME_ERROR_STR_LEN, no_error, mtime_strerror,           &
       &                               OPERATOR(-), OPERATOR(+), OPERATOR(>), OPERATOR(*),          &
       &                               ASSIGNMENT(=), OPERATOR(==), OPERATOR(>=), OPERATOR(/=),     &
       &                               event, eventGroup, newEvent,                                 &
       &                               addEventToEventGroup, isCurrentEventActive
  USE mo_event_manager,          ONLY: initEventManager, addEventGroup, getEventGroup, printEventGroup

  USE mo_hamocc_types,          ONLY: t_hamocc_state
  USE mo_ocean_physics,         ONLY: update_ho_params

  USE mo_ocean_tracer_transport_types

  !! Needed to test advection of velocity 
  USE mo_ocean_velocity_advection, ONLY: veloc_adv_horz_mimetic, &
      & veloc_adv_vert_mimetic

  USE mo_ocean_tracer_transport_horz, ONLY: advect_horz, diffuse_horz
  USE mo_ocean_tracer_transport_vert, ONLY: advect_flux_vertical
  USE mo_sync,                        ONLY: sync_c, sync_c1, sync_patch_array, sync_patch_array_mult
  
  !-------------------------------------------------------------------------
    IMPLICIT NONE
  PRIVATE

  PUBLIC :: ocean_test_zstar_advection
  
!  CHARACTER(len=12)           :: debug_string = 'testbed     '  ! Output of module for 1 line debug
  
  !-------------------------------------------------------------------------
CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! FIXME: Notes and overview
  !! Routines for testing advection with zstar are below
  !! Routines have been copied from other files and modified for zstar
  !! 1. Ideally, coefficients for the modified discretization should be calculated at one place
  !! and passed as arguments to repeat calculations
  !! 2. One routine can be used to calculate both low and high order flux for speedup
  !! 3. Variable to calculate depth needs to be clarified
  !! 4. This can be converted to using generalized vertical co-ordinates by using only 
  !! the coefficient of dz as variables to be modified
  !-------------------------------------------------------------------------

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
  !! Adapted for zstar
  !! FIXME: The limiter assumes no knowledge of eta for next time step which 
  !! would be required if the formulation were to be correct
  !!
  SUBROUTINE limiter_ocean_zalesak_horz_zstar( patch_3d,&
    & vert_velocity,          &
    & tracer,                 &
    & p_mass_flx_e,           &
    & flx_tracer_low,         &    
    & flx_tracer_high,        &
    & div_adv_flux_vert,      &   
    & eta,                    &   
    & operators_coefficients, &
    & flx_tracer_final )       
    
    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
    REAL(wp),INTENT(inout)              :: vert_velocity(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp), INTENT(inout)             :: tracer           (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: p_mass_flx_e     (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_low   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: flx_tracer_high  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp), INTENT(inout)             :: flx_tracer_final (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp), INTENT(in)                 :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !Surface ht 
    TYPE(t_operator_coeff),INTENT(in)   :: operators_coefficients
    
    !Local variables
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
    INTEGER :: edge_index, level, blockNo, jc,  cell_connect, sum_lsm_quad_edge, ctr
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d    
    
    INTEGER  :: bt_lev 
    REAL(wp) :: H_l, eta_l
    REAL(wp) :: coeff_l(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)      
    TYPE(t_subset_range), POINTER :: all_cells 

    !-------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    all_cells       => patch_2d%cells%ALL
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
 
    !-----------------------------------------------------------------------
    
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        bt_lev = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)      
 
        H_l   = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, blockNo)

        eta_l = eta(jc, blockNo) 

        coeff_l(jc, blockNo) = (H_l + eta_l)/H_l
        
      END DO
    END DO
    
    !-----------------------------------------------------------------------

  
!ICON_OMP_PARALLEL

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
!ICON_OMP z_fluxdiv_c) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      
      z_tracer_new_low(:,:,blockNo)    = 0.0_wp
      z_tracer_update_horz(:,:,blockNo)= 0.0_wp
      z_tracer_max(:,:,blockNo)        = 0.0_wp
      z_tracer_min(:,:,blockNo)        = 0.0_wp

      DO jc = start_index, end_index
        IF (patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo) < 1) CYCLE
        
        ! 3. Compute the complete (with horizontal and vertical divergence) updated low order solution z_tracer_new_low
        DO level = start_level  , MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)       
          !  compute divergence of low order fluxes
          z_fluxdiv_c = 0
          DO cell_connect = 1, patch_2d%cells%num_edges(jc,blockNo)
            z_fluxdiv_c =  z_fluxdiv_c + &
              & flx_tracer_low(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect)) * &
              & operators_coefficients%div_coeff(jc,level,blockNo,cell_connect)
          ENDDO

          delta_z     = coeff_l(jc, blockNo)*patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)
          delta_z_new = coeff_l(jc, blockNo)*patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)
          !
              
          z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
            & - dtime * (z_fluxdiv_c+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new
            
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

    ENDDO
!ICON_OMP_END_DO

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
        DO level = start_level  , MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)
          inv_prism_thick_new(level) = patch_3D%p_patch_1d(1)%inv_prism_thick_c(jc,level,blockNo) &
            & /coeff_l(jc, blockNo)
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
        ENDDO
      ENDDO
    ENDDO
!ICON_OMP_END_DO

    
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


  END SUBROUTINE limiter_ocean_zalesak_horz_zstar
  !-------------------------------------------------------------------------

  

  !-----------------------------------------------------------------------------
  ! the map_edges2edges_viacell_3d_mlev_constZs optimized for triangles
  ! modified for zstar 
  SUBROUTINE map_edges2edges_zstar( patch_3d, vn_e, eta, scalar, operators_coefficients, out_vn_e)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(in)                 :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)                 :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !Surface ht 
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(inout)              :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)                 :: scalar(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !Local variables
    INTEGER :: startLevel, endLevel
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_subset_range), POINTER :: all_edges 
    TYPE(t_subset_range), POINTER :: all_cells 
    TYPE(t_patch), POINTER :: patch_2d
    REAL(wp), POINTER :: coeffs(:,:,:,:)

    ! omp private
    ! defined in a struct for omp clearness
    TYPE omp_local_private
    
      INTEGER :: start_edge_index, end_edge_index
      INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
      INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
      INTEGER :: edge_11_block, edge_12_block, edge_13_block
      INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
      INTEGER :: edge_21_block, edge_22_block, edge_23_block
      
    END TYPE omp_local_private

    REAL(wp) :: coeff_e(nproma, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: H1, H2, H_e, eta_e 
    INTEGER  :: bt_lev 
    REAL(wp) :: H_l, eta_l
    REAL(wp) :: coeff_l(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)      
    INTEGER  :: start_index, end_index

    INTEGER, DIMENSION(:,:,:), POINTER :: idx, blk
    INTEGER  :: id1, id2, bl1, bl2 

!     TYPE(omp_local_private) :: omp_this 
    INTEGER :: jc, je, blockNo, level
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
    INTEGER :: edge_11_block, edge_12_block, edge_13_block
    INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
    INTEGER :: edge_21_block, edge_22_block, edge_23_block
   !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')
    
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    all_edges  => patch_2d%edges%ALL
    all_cells  => patch_2d%cells%ALL
    edges_in_domain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = n_zlev
    coeffs => operators_coefficients%edge2edge_viacell_coeff
    idx      => patch_3D%p_patch_2D(1)%edges%cell_idx
    blk      => patch_3D%p_patch_2D(1)%edges%cell_blk
 
    !-----------------------------------------------------------------------
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)
      DO je =  start_edge_index, end_edge_index
        !! FIXME: Assuming levels are the same at the bottom
        bt_lev = patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)      
 
        id1 = idx(je, blockNo, 1)
        id2 = idx(je, blockNo, 2)
        bl1 = blk(je, blockNo, 1)
        bl2 = blk(je, blockNo, 2)
 
        H1  = patch_3d%p_patch_1d(1)%depth_CellInterface(id1, bt_lev + 1, bl1)
        H2  = patch_3d%p_patch_1d(1)%depth_CellInterface(id2, bt_lev + 1, bl2)

        H_e = 0.5_wp*(H1 + H2)

        eta_e = 0.5_wp*( eta(id1, bl1) + eta(id2, bl2) )

        coeff_e(je, blockNo) = (H_e + eta_e)/H_e
        
      END DO
    END DO

    !-----------------------------------------------------------------------
    
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        bt_lev = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)      
 
        H_l   = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, blockNo)

        eta_l = eta(jc, blockNo) 

        coeff_l(jc, blockNo) = (H_l + eta_l)/H_l
        
      END DO
    END DO
    
    !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------
    
!ICON_OMP_PARALLEL_DO PRIVATE( je, level, start_edge_index, end_edge_index, &
!ICON_OMP cell_1_index, cell_2_index, cell_1_block, cell_2_block, &
!ICON_OMP edge_11_index, edge_12_index, edge_13_index, &
!ICON_OMP edge_11_block, edge_12_block, edge_13_block, &
!ICON_OMP edge_21_index, edge_22_index, edge_23_index, &
!ICON_OMP edge_21_block, edge_22_block, edge_23_block) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      
      DO je =  start_edge_index, end_edge_index
        
        cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
        cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)
        
        edge_11_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
        edge_12_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
        edge_13_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
        edge_11_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
        edge_12_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
        edge_13_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)
        
        edge_21_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
        edge_22_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
        edge_23_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
        edge_21_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
        edge_22_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
        edge_23_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)
                
        ! levels
        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          
          out_vn_e(je, level, blockNo) =  &
            & (  vn_e(edge_11_index, level, edge_11_block) * coeffs(je, level, blockNo, 1)      &
            &    * patch_3d%p_patch_1d(1)%prism_thick_e(edge_11_index, level, edge_11_block)    &
            & * coeff_e(edge_11_index, edge_11_block)                                           &
            &  + vn_e(edge_12_index, level, edge_12_block) * coeffs(je, level, blockNo, 2)      &
            &    * patch_3d%p_patch_1d(1)%prism_thick_e(edge_12_index, level, edge_12_block)    &
            & * coeff_e(edge_12_index, edge_12_block)                                           &
            &  + vn_e(edge_13_index, level, edge_13_block) * coeffs(je, level, blockNo, 3)      &
            &    * patch_3d%p_patch_1d(1)%prism_thick_e(edge_13_index, level, edge_13_block)    &
            & * coeff_e(edge_13_index, edge_13_block)                                           &
            & ) * scalar(cell_1_index, level, cell_1_block)                                     &
            & + &
            & (  vn_e(edge_21_index, level, edge_21_block) * coeffs(je, level, blockNo, 4)           &
            &   * patch_3d%p_patch_1d(1)%prism_thick_e(edge_21_index, level, edge_21_block)     &
            & * coeff_e(edge_21_index, edge_21_block)                                           &
            &  + vn_e(edge_22_index, level, edge_22_block) * coeffs(je, level, blockNo, 5)      &
            &  * patch_3d%p_patch_1d(1)%prism_thick_e(edge_22_index, level, edge_22_block)      &
            & * coeff_e(edge_22_index, edge_22_block)                                           &
            &  + vn_e(edge_23_index, level, edge_23_block) * coeffs(je, level, blockNo, 6)      &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_23_index, level, edge_23_block)       &
            & * coeff_e(edge_23_index, edge_23_block)                                           &
            & ) * scalar(cell_2_index, level, cell_2_block)                                      
      
        END DO !levels
    
      END DO
      
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO
     
  END SUBROUTINE map_edges2edges_zstar
  !-----------------------------------------------------------------------------
  




  
  SUBROUTINE upwind_zstar_hflux_oce( patch_3d, cell_value, eta, edge_vn, edge_upwind_flux, opt_start_level, opt_end_level )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)              :: cell_value   (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)      !< advected cell centered variable
    REAL(wp), INTENT(in)              :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !Surface ht 
    REAL(wp), INTENT(in)              :: edge_vn    (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)       !< normal velocity on edges
    REAL(wp), INTENT(inout)           :: edge_upwind_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level    ! optional vertical end level
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER, DIMENSION(:,:,:), POINTER :: idx, blk
    INTEGER  :: start_level, end_level
    INTEGER  :: start_index, end_index
    INTEGER  :: edge_index, level, blockNo         !< index of edge, vert level, block
    INTEGER  :: id1, id2, bl1, bl2 
    INTEGER  :: bt_level 
    REAL(wp) :: H1, H2, coeff1, coeff2, eta1, eta2 
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    idx             => patch_3D%p_patch_2D(1)%edges%cell_idx
    blk             => patch_3D%p_patch_2D(1)%edges%cell_blk
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
        bt_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo)      

        id1 = idx(edge_index, blockNo, 1)
        id2 = idx(edge_index, blockNo, 2)
        bl1 = blk(edge_index, blockNo, 1)
        bl2 = blk(edge_index, blockNo, 2)
 
        !! Get height of the cell center at mid point from the bottom
        H1  = patch_3d%p_patch_1d(1)%depth_CellInterface(id1,bt_level+1,bl1)
        H2  = patch_3d%p_patch_1d(1)%depth_CellInterface(id2,bt_level+1,bl2) 
        
        eta1 = eta(id1, bl1) 
        eta2 = eta(id2, bl2) 

        !! Transform to z from zstar
        coeff1  = (H1 + eta1)/H1 
        coeff2  = (H2 + eta2)/H2 

        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          edge_upwind_flux(edge_index,level,blockNo) =  &
             0.5_wp * (        edge_vn(edge_index,level,blockNo)  *           &
               & ( coeff1*cell_value(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1)) + &
               &   coeff2*cell_value(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2)) ) &
               &   - ABS( edge_vn(edge_index,level,blockNo) ) *               &
               & ( coeff2*cell_value(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2)) - &
               &   coeff1*cell_value(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1)) ) )
          
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    
  END SUBROUTINE upwind_zstar_hflux_oce
  !-----------------------------------------------------------------------
   


  !-------------------------------------------------------------------------
  !>
  !! Computation of new vertical velocity using continuity equation
  !! Calculate diagnostic vertical velocity from horizontal velocity using the
  !! incommpressibility condition in the continuity equation.
  !! vertical velocity is integrated from bottom to topLevel
  !! vertical velocity is negative for positive divergence
  !! of horizontal velocity
  !!
  SUBROUTINE calc_vert_velocity_zstar( patch_3d, ocean_state, op_coeffs, eta)
    TYPE(t_patch_3d), TARGET :: patch_3d       ! patch on which computation is performed
    TYPE(t_hydro_ocean_state) :: ocean_state
    TYPE(t_operator_coeff), INTENT(in) :: op_coeffs
    REAL(wp), INTENT(in)               :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !Surface ht 
    ! Local variables
    INTEGER :: jc, jk, blockNo, je, z_dolic, start_index, end_index
    REAL(wp) :: z_c(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks), z_abort
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells, cells_owned
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp), POINTER :: vertical_velocity(:,:,:)
    REAL(wp) :: H_l, eta_l, coeff_l 
    INTEGER  :: bt_level 
    REAL(wp) :: z_adv_flux_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: div_z_c(nproma,n_zlev)
    REAL(wp) :: div_z_depth_int_c(nproma)
    REAL(wp) :: temp(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), POINTER :: cell_thickness(:,:,:)

    CHARACTER(len=*), PARAMETER :: method_name='mo_ocean_ab_timestepping_mimetic:alc_vert_velocity_mim_bottomup'
    !-----------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain  => patch_2D%cells%in_domain
    cells_owned      => patch_2D%cells%owned
    all_cells        => patch_2D%cells%all
    edges_in_domain  => patch_2D%edges%in_domain
    vertical_velocity=> ocean_state%p_diag%w
    cell_thickness   => patch_3D%p_patch_1d(1)%prism_thick_c
    ! due to nag -nan compiler-option:
    !------------------------------------------------------------------
    ! Step 1) Calculate divergence of horizontal velocity at all levels
    !------------------------------------------------------------------
    !-------------------------------------------------------------------------------
 
    CALL map_edges2edges_viacell_3d_const_z( patch_3d, ocean_state%p_prog(nold(1))%vn, &
      & op_coeffs, ocean_state%p_diag%mass_flx_e)

    !! Use trick of using constant coefficient to get mass flux 
    temp = 1.0_wp
  
    CALL map_edges2edges_zstar( patch_3d, ocean_state%p_prog(nold(1))%vn, eta, &
      & temp, op_coeffs, z_adv_flux_h)

    vertical_velocity = 0.0_wp

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      
      CALL div_oce_3D_onTriangles_onBlock(z_adv_flux_h, patch_3D, op_coeffs%div_coeff, &
        & div_z_c(:,:), blockNo=blockNo, start_index=start_index, &
        & end_index=end_index, start_level=1, end_level=n_zlev)

      DO jc = start_index, end_index
        !! Get summation over depth of divergence for RHS of sfc equation 
        div_z_depth_int_c(jc) = SUM(div_z_c(jc, 1:patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)))
       
        bt_level = patch_3d%p_patch_1d(1)%dolic_c(jc, blockNo)      

        !! Get height of the cell center at mid point from the bottom
        H_l    = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_level+1, blockNo)
        
        eta_l  = eta(jc, blockNo) 

        !! Transform to z from zstar
        coeff_l = (H_l + eta_l)/H_l 

        !! Trick here is that dw/dz = -div(v) or w2-w1=-dz.div(v)=-div(dz.v)
        !! In z*, it becomes d(J.w)/dz = -div(J.v) or w2-w1=-dz*.div(J.v)/J=-div(dz*.J.v)/J
        !! J being the vertical co-ordinate Jacobian
        !! d_t eta = -div_z_depth_int_c(jc)/H_l
        DO jk = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), 1, -1
        
          vertical_velocity(jc,jk,blockNo) = vertical_velocity(jc,jk+1,blockNo) - &
            & ( div_z_c(jc,jk) - &
            & cell_thickness(jc, jk, blockNo)*div_z_depth_int_c(jc)/H_l)/coeff_l
          !! Switch to the below lines to use w^* for d_t eta = 0
          !! This will give incorrect results for topmost layer where a zero
          !! flux boundary condn is used in vertical advection
!          vertical_velocity(jc,jk,blockNo) = vertical_velocity(jc,jk+1,blockNo) - &
!            &  div_z_c(jc,jk)/coeff_l


        END DO
      END DO
    END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO
    
    CALL sync_patch_array(sync_c,patch_2D,vertical_velocity)

  END SUBROUTINE calc_vert_velocity_zstar
  !-------------------------------------------------------------------------
 


   !------------------------------------------------------------------------
  SUBROUTINE tracer_diffusion_vertical_implicit_zstar( &
    & patch_3d,                  &
    & ocean_tracer,              &
    & a_v,     &
    & eta)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: ocean_tracer
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)
    REAL(wp), INTENT(in)                 :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !Surface ht 
    !
    INTEGER :: cell_block, start_index, end_index
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d

    !-----------------------------------------------------------------------
    cells_in_domain       =>  patch_3d%p_patch_2d(1)%cells%in_domain
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO cell_block = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, cell_block, start_index, end_index)

      CALL tracer_diffusion_vertical_implicit_zstar_onBlock( &
        & patch_3d,                  &
        & ocean_tracer,              &
        & a_v, eta,                  &
        & cell_block, start_index, end_index)

    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE tracer_diffusion_vertical_implicit_zstar
  !------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!Subroutine implements implicit vertical diffusion for scalar fields.
  !>
  !! The result ocean_tracer%concetration is calculated on domain_cells
  !-------------------------------------------------------------------------
  SUBROUTINE tracer_diffusion_vertical_implicit_zstar_onBlock( &
    & patch_3d,                &
    & ocean_tracer,            &
    & a_v, eta,                &
    & blockNo, start_index, end_index) !,  &
    ! & diff_column)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: ocean_tracer
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)
    INTEGER, INTENT(in)                  :: blockNo, start_index, end_index
    REAL(wp), INTENT(in)                 :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !Surface ht 
    !
    !
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)! , nb(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    REAL(wp) :: column_tracer(1:n_zlev)
    REAL(wp) :: dt_inv, diagonal_product
    REAL(wp), POINTER :: field_column(:,:,:)
    INTEGER :: bottom_level
    INTEGER :: cell_index, level
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    
    REAL(wp) :: H_l, eta_l, str_l, inv_str_l 

    !-----------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2d%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
    dt_inv = 1.0_wp/dtime
    
    DO cell_index = start_index, end_index
      bottom_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)

      !! Get height of the cell center at mid point from the bottom
      H_l    = patch_3d%p_patch_1d(1)%depth_CellInterface(cell_index, bottom_level+1, blockNo)
      
      eta_l  = eta(cell_index, blockNo) 

      !! Transform to z from zstar
      str_l      = (H_l + eta_l)/H_l 
      inv_str_l  = 1._wp/str_l 

 
      IF (bottom_level < 2 ) CYCLE ! nothing to diffuse

      DO level=1,bottom_level
        inv_prism_thickness(level)        = inv_str_l*patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_index,level,blockNo)
        inv_prisms_center_distance(level) = inv_str_l*patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(cell_index,level,blockNo)
      ENDDO

      !------------------------------------
      ! Fill triangular matrix
      ! b is diagonal, a is the upper diagonal, c is the lower
      !   top level
      a(1) = 0.0_wp
      c(1) = -a_v(cell_index,2,blockNo) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
      b(1) = dt_inv - c(1)
      DO level = 2, bottom_level-1
        a(level) = - a_v(cell_index,level,blockNo)   * inv_prism_thickness(level) * inv_prisms_center_distance(level)
        c(level) = - a_v(cell_index,level+1,blockNo) * inv_prism_thickness(level) * inv_prisms_center_distance(level+1)
        b(level) = dt_inv - a(level) - c(level)
      END DO
      ! bottom
      a(bottom_level) = -a_v(cell_index,bottom_level,blockNo) * &
        & inv_prism_thickness(bottom_level) * inv_prisms_center_distance(bottom_level)
      b(bottom_level) = dt_inv - a(bottom_level)

      ! precondition: set diagonal equal to diagonal_product
      diagonal_product = PRODUCT(b(1:bottom_level))

      DO level = 1, bottom_level
        fact(level) = diagonal_product / b(level)
        a(level)  = a(level)  * fact(level)
        b(level)  = diagonal_product
        c(level)  = dt_inv * fact(level) - a(level) - b(level)
 
        column_tracer(level) = field_column(cell_index,level,blockNo) * dt_inv * fact(level)

      ENDDO
      c(bottom_level) = 0.0_wp

      !------------------------------------
      ! solver from lapack
      !
      ! eliminate lower diagonal
      DO level=bottom_level-1, 1, -1
        fact(level+1)  = c( level ) / b( level+1 )
        b( level ) = b( level ) - fact(level+1) * a( level +1 )
        column_tracer( level ) = column_tracer( level ) - fact(level+1) * column_tracer( level+1 )
      ENDDO

      !     Back solve with the matrix U from the factorization.
      column_tracer( 1 ) = column_tracer( 1 ) / b( 1 )
      DO level =  2, bottom_level
        column_tracer( level ) = ( column_tracer( level ) - a( level ) * column_tracer( level-1 ) ) / b( level )
      ENDDO

      DO level = 1, bottom_level
        ocean_tracer%concentration(cell_index,level,blockNo) = column_tracer(level)
      ENDDO
    
    ENDDO ! cell_index
    
  END SUBROUTINE tracer_diffusion_vertical_implicit_zstar_onBlock
  !------------------------------------------------------------------------
    

  !-------------------------------------------------------------------------
  !> Setup a test case that advects tracers for testing with zstar
  !  Should start by using the low order horizontal advection
  SUBROUTINE ocean_test_zstar_advection( patch_3d, ocean_state, &
    & this_datetime, ocean_surface, physics_parameters,             &
    & ocean_ice,operators_coefficients)
    
    TYPE(t_patch_3d), POINTER, INTENT(in)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_ocean_surface)                            :: ocean_surface
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(in)          :: operators_coefficients
    
    ! local variables
    TYPE (t_hamocc_state)        :: hamocc_State
    INTEGER :: jstep, jg
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: jstep0 ! start counter for time loop
    INTEGER :: i
    INTEGER :: tracer_index 
    TYPE(timedelta), POINTER :: model_time_step => NULL()
    
    !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
    TYPE(t_tracer_collection) , POINTER              :: old_tracer_collection, new_tracer_collection
    TYPE(t_ocean_transport_state)                    :: transport_state
    
    INTEGER  :: jb, jc, je, level 
    REAL(wp) :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: delta_t, delta_z,delta_z_new, delta_z1,delta_z_new1
    REAL(wp) :: div_adv_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: flux_horz(nproma,n_zlev, patch_3d%p_patch_2D(1)%nblks_e)
    REAL(wp) :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: top_bc(nproma)
    INTEGER  :: start_index, end_index
    INTEGER  :: start_cell_index, end_cell_index
    REAL(wp) :: z_adv_flux_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_adv_low (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_adv_high(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z2(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: H_l, eta_l, coeff_l 
    INTEGER  :: bt_level 
    
    TYPE(t_ocean_tracer), POINTER :: new_tracer
    TYPE(t_ocean_tracer), POINTER :: old_tracer

    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    
    REAL(wp) :: temp(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & method_name = 'mo_ocean_testbed_modules:ocean_test_zstar_advection'
    !------------------------------------------------------------------
    
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2d%edges%in_domain
    CALL datetimeToString(this_datetime, datestring)

    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)
    
    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(method_name), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    !! sea surface height type 201 set explicitly since we don't want
    !! the grid to change
    eta = 0.  
    ! Initialize eta for zstar
    ! #slo#: simple elevation between 30W and 30E (pi/3.)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        IF ( patch_3d%lsm_c(jc, 1, jb) <= sea_boundary ) THEN
           eta(jc, jb) = 10.0_wp * &
            & SIN(patch_2d%cells%center(jc, jb)%lon * 6.0_wp) &
            & * COS(patch_2d%cells%center(jc, jb)%lat * 3.0_wp)
        ENDIF
      END DO
    END DO

!    !---------------------------------------------------------------------
!    !-FIXME: test divergence of constant fn
!    !---------------------------------------------------------------------
!  
!    ! calc_vert_vel uses vn_time_weighter instead of vn
!    ocean_state(jg)%p_diag%vn_time_weighted = ocean_state(jg)%p_prog(nold(1))%vn
!    
!    !! Update mass_flux and w 
!    CALL calc_vert_velocity_zstar( patch_3d, ocean_state(jg),operators_coefficients, eta)
!
!    ! fill transport_state
!    transport_state%patch_3d    => patch_3d
!    transport_state%h_old       => ocean_state(jg)%p_prog(nold(1))%h
!    transport_state%h_new       => ocean_state(jg)%p_prog(nnew(1))%h
!    transport_state%vn          => ocean_state(jg)%p_prog(nold(1))%vn
!    transport_state%mass_flux_e => ocean_state(jg)%p_diag%mass_flx_e
!    transport_state%w           => ocean_state(jg)%p_diag%w
!
!    temp = 1.0_wp
!
!    CALL advect_flux_vertical( patch_3d,&
!        & temp, &
!        & transport_state,                           &
!        & operators_coefficients,                     &
!        & div_adv_flux_vert)
!
!    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
!      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
!      DO jc = start_cell_index, end_cell_index
!        bt_level = patch_3d%p_patch_1d(1)%dolic_c(jc, jb)      
!
!        !! Get height of the cell center at mid point from the bottom
!        H_l    = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_level+1, jb)
!        
!        eta_l  = eta(jc, jb) 
!
!        !! Transform to z from zstar
!        coeff_l = (H_l + eta_l)/H_l 
!
!        div_adv_flux_vert(jc, :, jb) = coeff_l*div_adv_flux_vert(jc, :, jb)
!      END DO
!    END DO
!
!    !---------------------------------------------------------------------
!    !-Horizontal  advection
!    !---------------------------------------------------------------------
!    CALL upwind_zstar_hflux_oce( patch_3d,  &
!      & temp, &
!      & eta, &
!      & transport_state%mass_flux_e,         &
!      & z_adv_flux_h)                         
! 
!    !Calculate divergence of advective fluxes
!    CALL div_oce_3d( z_adv_flux_h, patch_3D, operators_coefficients%div_coeff, &
!      & div_adv_flux_horz, subset_range=cells_in_domain )
!
!    CALL map_edges2edges_zstar( patch_3d, transport_state%vn, eta, &
!      & temp, operators_coefficients, z2)
!
!    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
!      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
!      DO jc = start_cell_index, end_cell_index
!        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
!
!          if ( (jb == 4) .AND. (jc == 10) ) THEN
!!            write(*, *) level, div_adv_flux_horz(jc,level,jb) , temp(jc,level,jb)
!!              & transport_state%w(jb, level + 1, jc) - transport_state%w(jb, level, jc), &
!!              & patch_3d%p_patch_1d(1)%prism_thick_c(jb, level, jc)
!          ENDIF
!    
!        ENDDO
!      ENDDO
!    ENDDO
! 
!    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
!      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
!      DO je = start_index, end_index
!        DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
!           if ( (jb == 4) .AND. (je == 10) ) THEN
!            write(*, *) level, z_adv_flux_h(je,level,jb) , z2(je,level,jb)
!          ENDIF
!    
!        ENDDO
!      ENDDO
!    ENDDO
!
!
!    !---------------------------------------------------------------------
!    !-FIXME: test end 
!    !---------------------------------------------------------------------
 
    jstep0 = 0
    DO jstep = (jstep0+1), (jstep0+nsteps)

        CALL datetimeToString(this_datetime, datestring)
        WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
        CALL message (TRIM(method_name), message_text)
 
        old_tracer_collection => ocean_state(jg)%p_prog(nold(1))%tracer_collection
        new_tracer_collection => ocean_state(jg)%p_prog(nnew(1))%tracer_collection

        IF (no_tracer>=1) THEN

          ! calc_vert_vel uses vn_time_weighter instead of vn
          ocean_state(jg)%p_diag%vn_time_weighted = ocean_state(jg)%p_prog(nold(1))%vn

          !! Update w and mass_flx_e for tracer advection
          CALL calc_vert_velocity_zstar( patch_3d, ocean_state(jg),operators_coefficients, eta)

          ! fill transport_state
          transport_state%patch_3d    => patch_3d
          transport_state%h_old       => ocean_state(jg)%p_prog(nold(1))%h
          transport_state%h_new       => ocean_state(jg)%p_prog(nnew(1))%h
          transport_state%vn          => ocean_state(jg)%p_prog(nold(1))%vn
          transport_state%w           => ocean_state(jg)%p_diag%w
          transport_state%mass_flux_e => ocean_state(jg)%p_diag%mass_flx_e

          ! fill diffusion coefficients
          old_tracer_collection%tracer(1)%hor_diffusion_coeff => physics_parameters%TracerDiffusion_coeff(:,:,:,1)
          old_tracer_collection%tracer(1)%ver_diffusion_coeff => physics_parameters%a_tracer_v(:,:,:,1)
          DO i = 2, old_tracer_collection%no_of_tracers
            old_tracer_collection%tracer(i)%hor_diffusion_coeff => physics_parameters%TracerDiffusion_coeff(:,:,:,2)
            old_tracer_collection%tracer(i)%ver_diffusion_coeff => physics_parameters%a_tracer_v(:,:,:,2)
          ENDDO
        ENDIF
        !------------------------------------------------------------------------

        DO tracer_index = 1, old_tracer_collection%no_of_tracers
          
          old_tracer => old_tracer_collection%tracer(tracer_index)
          new_tracer => new_tracer_collection%tracer(tracer_index)
          IF ( old_tracer%is_advected) THEN
            
            !-------------------------------------------------------------------------------
            patch_2D        => patch_3d%p_patch_2d(1)
            cells_in_domain => patch_2D%cells%in_domain
            edges_in_domain => patch_2D%edges%in_domain
            delta_t = dtime
            !---------------------------------------------------------------------
         
            ! these are probably not necessary
            div_adv_flux_vert = 0.0_wp
            div_adv_flux_horz = 0.0_wp
            div_diff_flux_horz = 0.0_wp
            
            !---------------------------------------------------------------------
            !-Vertical advection
            !---------------------------------------------------------------------
            IF ( l_with_vert_tracer_advection ) THEN
        
              CALL advect_flux_vertical( patch_3d,&
                & old_tracer%concentration, &
                & transport_state,                           &
                & operators_coefficients,                     &
                & div_adv_flux_vert)
        
            ENDIF  ! l_with_vert_tracer_advection
 
            !---------------------------------------------------------------------
            !-Horizontal  advection
            !---------------------------------------------------------------------
            CALL upwind_zstar_hflux_oce( patch_3d,  &
              & old_tracer%concentration, &
              & eta, &
              & transport_state%mass_flux_e,         &
              & z_adv_flux_h)
            z_adv_low = z_adv_flux_h
 
            CALL map_edges2edges_zstar( patch_3d, transport_state%vn, eta, &
              & old_tracer%concentration, operators_coefficients, z_adv_flux_h)
            z_adv_high = z_adv_flux_h

            CALL limiter_ocean_zalesak_horz_zstar( patch_3d,   &
              & transport_state%w,           &
              & old_tracer%concentration,              &
              & transport_state%mass_flux_e,           &
              & z_adv_low,                             &
              & z_adv_high,                            &
              & div_adv_flux_vert,                     &            
              & eta,                                   &            
              & operators_coefficients,                &
              & z_adv_flux_h)                          

            !Calculate divergence of advective fluxes
            CALL div_oce_3d( z_adv_flux_h, patch_3D, operators_coefficients%div_coeff, &
              & div_adv_flux_horz, subset_range=cells_in_domain )

        !ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, level, &
        !ICON_OMP delta_z, delta_z_new, top_bc) ICON_OMP_DEFAULT_SCHEDULE
            DO jb = cells_in_domain%start_block, cells_in_domain%end_block
              CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
              DO jc = start_cell_index, end_cell_index
                bt_level = patch_3d%p_patch_1d(1)%dolic_c(jc, jb)      

                !! Get height of the cell center at mid point from the bottom
                H_l    = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_level+1, jb)
        
                eta_l  = eta(jc, jb) 

                !! Transform to z from zstar
                coeff_l = (H_l + eta_l)/H_l 

                !! d_z*(coeff*w*C) = coeff*d_z(w*C) since coeff is constant for each column
                div_adv_flux_vert(jc, :, jb) = coeff_l*div_adv_flux_vert(jc, :, jb)

                DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
        
                  new_tracer%concentration(jc,level,jb) =                          &
                    &  old_tracer%concentration(jc,level,jb) -                     &
                    &  (delta_t /  ( coeff_l*patch_3d%p_patch_1D(1)%prism_thick_c(jc,level,jb) ) )    &
                    & * (div_adv_flux_horz(jc,level,jb)  + div_adv_flux_vert(jc,level,jb))

!                  write(*, *) level, ocean_state(jg)%p_diag%w (jc,level,jb), &
!                    & div_adv_flux_horz(jc,level,jb) , div_adv_flux_vert(jc,level,jb)

                ENDDO
        
              END DO
            END DO
        !ICON_OMP_END_PARALLEL_DO
        
            CALL tracer_diffusion_vertical_implicit_zstar( &
                   & patch_3d,                  &
                   & new_tracer,                &
                   & old_tracer_collection%tracer(tracer_index)%ver_diffusion_coeff, &
                   & eta) 
 

            CALL sync_patch_array(sync_c, patch_2D, new_tracer%concentration)
        
          ENDIF
        END DO



        ! One integration cycle finished on the lowest grid level (coarsest
        ! resolution). Set model time.
        model_time_step => newTimedelta('+', 0, 0, 0, 0, 0, NINT(dtime), 0)
        this_datetime = this_datetime + model_time_step
        CALL deallocateTimedelta(model_time_step) 
          
        CALL output_ocean( patch_3d, &
          & ocean_state,             &
          & this_datetime,                &
          & ocean_surface,          &
          & ocean_ice,               &
          & jstep, jstep0)

        ! Shift time indices for the next loop
        ! this HAS to ge into the restart files, because the start with the following loop
        CALL update_time_indices(jg)
        ! update intermediate timestepping variables for the tracers
        ! velocity
        ocean_state(jg)%p_aux%g_nm1 = ocean_state(jg)%p_aux%g_n
        ocean_state(jg)%p_aux%g_n   = 0.0_wp

    END DO
    
    CALL timer_stop(timer_total)
    
  END SUBROUTINE ocean_test_zstar_advection
  !-------------------------------------------------------------------------


END MODULE mo_ocean_testbed_zstar
