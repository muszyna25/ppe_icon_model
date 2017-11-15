!>
!! Provide an implementation of the ocean physics.
!!
!! Provide an implementation of the physical parameters and characteristics
!! for the hydrostatic ocean model.
!!
!! @author Stephan Lorenz, MPI
!! @author Peter Korn, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
!!  Modified by Stephan Lorenz,     MPI-M (2010-07)
!!    adapted to structures discussed in 2010-01.
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
MODULE mo_ocean_pp_scheme
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_ocean_nml,           ONLY: &
    & n_zlev, bottom_drag_coeff,                              &
    & HarmonicViscosity_reference, velocity_VerticalDiffusion_background,                 &
    & Temperature_VerticalDiffusion_background, Salinity_VerticalDiffusion_background, no_tracer,                       &
    & tracer_convection_MixingCoefficient,                                     &
    & BiharmonicViscosity_scaling, HarmonicViscosity_scaling, &
    & VelocityDiffusion_order,                                &
    & BiharmonicViscosity_reference,                          &
    & tracer_RichardsonCoeff, velocity_RichardsonCoeff,                    &
    & PPscheme_type,                                &
    & PPscheme_Constant_type,                       &
    & PPscheme_ICON_type,                        &
    & PPscheme_ICON_Edge_type,                   &
    & PPscheme_ICON_Edge_vnPredict_type,         &
    & PPscheme_MPIOM_type,                       &
    & use_wind_mixing,                                        &
    & HorizontalViscosity_SmoothIterations,                   &
    & convection_InstabilityThreshold,                        &
    & RichardsonDiffusion_threshold,                          &
    & lambda_wind, wma_diff, wma_visc,                        &
    & use_reduced_mixing_under_ice,                           &
    & k_tracer_dianeutral_parameter,                          &
    & k_tracer_isoneutral_parameter, k_tracer_GM_kappa_parameter,    &
    & GMRedi_configuration,GMRedi_combined,                   &
    & GM_only,Redi_only,                                      &
    & laplacian_form,                                         &
    & HorizontalViscosity_SpatialSmoothFactor,                &
    & VerticalViscosity_TimeWeight, OceanReferenceDensity,    &
    & tracer_TopWindMixing, WindMixingDecayDepth,             &
    & velocity_TopWindMixing, TracerHorizontalDiffusion_scaling, &
    &  Temperature_HorizontalDiffusion_Background,            &
    &  Temperature_HorizontalDiffusion_Reference,             &
    &  Salinity_HorizontalDiffusion_Background,               &
    &  Salinity_HorizontalDiffusion_Reference,                &
    &  HarmonicViscosity_background,                          &
    &  BiharmonicViscosity_background,                        &
    &  LeithHarmonicViscosity_background, LeithHarmonicViscosity_reference,    &
    &  LeithHarmonicViscosity_scaling,                                         &
    &  LeithBiharmonicViscosity_background, LeithBiharmonicViscosity_reference,&
    &  LeithBiharmonicViscosity_scaling,                       &
    &  LeithClosure_order,   LeithClosure_form, &
    &  TracerDiffusion_LeithWeight, Salinity_ConvectionRestrict, &
    &  max_turbulenece_TracerDiffusion_amplification, &
    &  ReferencePressureIndbars

  USE mo_ocean_physics_types, ONLY: t_ho_params, v_params, WindMixingDecay, WindMixingLevel
   !, l_convection, l_pp_scheme
  USE mo_parallel_config,     ONLY: nproma
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_impl_constants,      ONLY: success, max_char_length, min_dolic, sea
  USE mo_cdi_constants,       ONLY: grid_cell, grid_edge,                           &
    &                               grid_unstructured_edge, grid_unstructured_cell
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_util_dbg_prnt,       ONLY: dbg_print, debug_print_MaxMinMean
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state, t_onEdges_Pointer_3d_wp, t_onCells_HalfLevels_Pointer_wp, t_operator_coeff
  USE mo_ocean_state,         ONLY: oce_config
  USE mo_physical_constants,  ONLY: grav, sitodbar,sal_ref
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_dynamics_config,     ONLY: nold!, nnew
  USE mo_sea_ice_types,       ONLY: t_sfc_flx
  USE mo_run_config,          ONLY: dtime
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_var_list,            ONLY: add_var,                  &
    & new_var_list,             &
    & delete_var_list,          &
    & default_var_list_settings,&
    & add_ref
  USE mo_var_metadata,        ONLY: groups
  USE mo_cf_convention
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                 ONLY: datatype_pack16, DATATYPE_FLT32, DATATYPE_FLT64, filetype_nc2, &
    &                               GRID_UNSTRUCTURED
  USE mo_zaxis_type,          ONLY: &
    & za_depth_below_sea, za_depth_below_sea_half, za_surface
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, global_max, sync_patch_array_mult
  USE  mo_ocean_thermodyn,    ONLY: calculate_density_onColumn
  USE mo_ocean_math_operators,ONLY: div_oce_3d
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, &
    & timer_extra10, timer_extra11
  USE mo_statistics,          ONLY: global_minmaxmean
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = 'mo_pp_scheme'
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug

  PUBLIC :: update_PP_scheme
  PUBLIC :: ICON_PP_Edge_vnPredict_scheme

  REAL(wp), POINTER :: WindAmplitude_at10m(:,:)  
  REAL(wp), POINTER :: SeaIceConcentration(:,:)
  PUBLIC :: calculate_rho4GMRedi
  

CONTAINS

 !-------------------------------------------------------------------------
  !>
  !! 
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2016-12)
  !<Optimize:inUse:done>
  SUBROUTINE calculate_rho4GMRedi(patch_3d, temperature, salinity, rho_GM) 

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    !TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    REAL(wp), INTENT(IN)                 :: temperature(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(IN)                 :: salinity(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(INOUT)              :: rho_GM(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !TYPE(t_ho_params), INTENT(inout)     :: params_oce

    ! Local variables
    INTEGER :: jc, blockNo, je,jk, tracer_index
    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
    INTEGER :: start_index, end_index
    INTEGER :: levels

    REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev), density(n_zlev)
    REAL(wp) :: pressure(n_zlev), sal(n_zlev)
    !REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D

    !-------------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    !edges_in_domain => patch_2D%edges%in_domain
    !cells_in_domain => patch_2D%cells%in_domain
    all_cells       => patch_2D%cells%ALL
    !z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho
    levels = n_zlev


!     IF (ltimer) CALL timer_start(timer_extra10)

!ICON_OMP_PARALLEL PRIVATE(salinity,z_rho_up, z_rho_down)
    !sal(1:levels) = sal_ref

    z_rho_up(:)=0.0_wp
    z_rho_down(:)=0.0_wp
    rho_GM(:,:,:)=0.0_wp

!ICON_OMP_DO PRIVATE(start_index, end_index, jc, levels, jk, pressure, &
!ICON_OMP  tracer_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index

        levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
        IF (levels < 2) CYCLE

        !--------------------------------------------------------
        sal(1:levels)      = salinity(jc,1:levels,blockNo)

        pressure(2:levels) = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, 2:levels, blockNo) * OceanReferenceDensity * sitodbar

        z_rho_up(1:levels-1)  &
        &= calculate_density_onColumn(temperature(jc,1:levels-1,blockNo),&
                                    & salinity   (jc,1:levels-1,blockNo), pressure(2:levels), levels-1)
        z_rho_down(2:levels)  &
        &= calculate_density_onColumn(temperature(jc, 2:levels, blockNo),&
                                    & salinity   (jc, 2:levels, blockNo), pressure(2:levels), levels-1)


        DO jk = 2, levels
          rho_GM(jc,jk,blockNo)=0.5_wp*(z_rho_up(jk)+z_rho_down(jk))
        
          !z_vert_density_grad_c(jc,jk,blockNo) = (z_rho_down(jk) - z_rho_up(jk-1)) *  &
          !  & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,blockNo)
        END DO ! levels
        !ocean_state%p_diag%grad_rho_PP_vert(jc,2:levels,blockNo)=z_vert_density_grad_c(jc,2:levels,blockNo)
        rho_GM(jc,1,blockNo)=rho_GM(jc,2,blockNo)
 
      END DO ! index


    END DO ! blocks
!ICON_OMP_END_DO

!ICON_OMP_END_PARALLEL
!     IF (ltimer) CALL timer_stop(timer_extra10)
!     IF (ltimer) CALL timer_start(timer_extra11)
! !ICON_OMP_PARALLEL

!     IF (ltimer) CALL timer_stop(timer_extra11)
    DO jk=1,10
      CALL dbg_print('calc_rho4GMRedi: rho_GM',rho_GM(:,jk,:),&
        & module_name, idt_src, in_subset=all_cells) 
    END DO          

  END SUBROUTINE calculate_rho4GMRedi
  !-------------------------------------------------------------------------


!<Optimize:inUse:done>
  SUBROUTINE update_PP_scheme(patch_3d, ocean_state, fu10, concsum, params_oce,op_coeffs) !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    REAL(wp), TARGET                     :: fu10   (:,:) ! t_atmos_for_ocean%fu10
    REAL(wp), TARGET                     :: concsum(:,:) ! t_sea_ice%concsum
    TYPE(t_ho_params), INTENT(inout)     :: params_oce
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs

    INTEGER :: tracer_index
    !-------------------------------------------------------------------------
    WindAmplitude_at10m => fu10
    SeaIceConcentration => concsum

    SELECT CASE (PPscheme_type)
    CASE (PPscheme_Constant_type)
      !nothing to do!In sbr init_ho_params (see above)
      !tracer mixing coefficient params_oce%A_tracer_v(:,:,:, tracer_index) is already
      !initialzed with params_oce%A_tracer_v_back(tracer_index)
      !and velocity diffusion coefficient

      ! prepare independent logicals for PP and convection parametrizations - not yet activated
      ! IF (.NOT. (l_convection .AND. l_pp_scheme)) THEN
      RETURN

    CASE (PPscheme_ICON_type)
      CALL ICON_PP_scheme(patch_3d, ocean_state, params_oce)

    CASE (PPscheme_MPIOM_type)
      CALL MPIOM_PP_scheme(patch_3d, ocean_state, fu10, concsum, params_oce)

    CASE (PPscheme_ICON_Edge_type)
      CALL ICON_PP_Edge_scheme(patch_3d, ocean_state, params_oce)

    CASE (PPscheme_ICON_Edge_vnPredict_type)
!       CALL update_PhysicsParameters_ICON_PP_Tracer(patch_3d, ocean_state)
      CALL ICON_PP_Edge_scheme(patch_3d, ocean_state, params_oce)
      ! the velovity friction will be updated during dynamics

    CASE default
      CALL finish("update_ho_params", "unknown PPscheme_type")
    END SELECT

  END SUBROUTINE update_PP_scheme
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Update of ocean physics: This routine is used used only if time-dependent
  !! changes of physical parametrizations.
  !! Currently vertical mixing coefficients for tracers and vertical diffusivity are updated.
  !! Dependent on the local Richardson number the diffusivity are calculated
  !! (Large & Gent JPO 29, (1999), 449-464).
  !! The formulation follows the MPI-OM implementation as described in Marsland et al. (Ocean
  !! Modelling 5, 2003).
  !! The notational convention is also taken from this paper( cf. eqs (14) and (19)).
  !! What is missing is the fractional ice cover (see eqs. (15-16)).
  !! Eq. (18) is the Redi part that is not implemented, yet
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-02)
  !<Optimize:inUse:done>
  SUBROUTINE ICON_PP_scheme(patch_3d, ocean_state, params_oce) !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_ho_params), INTENT(inout)            :: params_oce

    ! Local variables
    INTEGER :: jc, blockNo, je,jk, tracer_index
    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
    INTEGER :: start_index, end_index
    INTEGER :: levels

    REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev), density(n_zlev)
    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
    REAL(wp) :: z_shear_cell, z_av0
    REAL(wp) :: z_ri_cell               (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)

    !Below is a set of variables and parameters for tracer and velocity
    REAL(wp), PARAMETER :: z_0               = 40.0_wp
    REAL(wp), PARAMETER :: z_c1_t            = 5.0_wp    !  PP diffusivity tuning constant
    REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp    !  PP viscosity tuning constant
    REAL(wp), PARAMETER :: z_threshold       = 5.0E-8_wp
    REAL(wp) :: diffusion_weight
    REAL(wp) :: z_grav_rho, z_inv_OceanReferenceDensity
    REAL(wp) :: density_differ_edge, mean_z_r
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D

    !-------------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    !cells_in_domain => patch_2D%cells%in_domain
    all_cells       => patch_2D%cells%ALL
    z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho
    levels = n_zlev

    !-------------------------------------------------------------------------
    z_av0 = velocity_RichardsonCoeff
    z_grav_rho                   = grav/OceanReferenceDensity
    z_inv_OceanReferenceDensity                = 1.0_wp/OceanReferenceDensity
    !-------------------------------------------------------------------------
!     IF (ltimer) CALL timer_start(timer_extra10)

!ICON_OMP_PARALLEL PRIVATE(salinity,z_rho_up, z_rho_down)
    salinity(1:levels) = sal_ref

    z_rho_up(:)=0.0_wp
    z_rho_down(:)=0.0_wp

!ICON_OMP_DO PRIVATE(start_index, end_index, jc, levels, jk, pressure, &
!ICON_OMP z_shear_cell, tracer_index, diffusion_weight) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      z_ri_cell(:,:, blockNo) = 0.0_wp
      z_vert_density_grad_c(:,:, blockNo) = 0.0_wp
      DO jc = start_index, end_index

        levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
        IF (levels < 2) CYCLE

        IF(no_tracer >= 2) THEN
            salinity(1:levels) = ocean_state%p_prog(nold(1))%tracer(jc,1:levels,blockNo,2)
        ENDIF

        !--------------------------------------------------------
        ! pressure in dbars
        pressure(2:levels) = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, 2:levels, blockNo) * ReferencePressureIndbars
        z_rho_up(1:levels-1)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,1:levels-1,blockNo,1), &
          & salinity(1:levels-1), pressure(2:levels), levels-1)
        z_rho_down(2:levels)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,2:levels,blockNo,1), &
          & salinity(2:levels), pressure(2:levels), levels-1)


        DO jk = 2, levels
        ocean_state%p_diag%rho_GM(jc,jk,blockNo)=0.5_wp*(z_rho_up(jk)+z_rho_down(jk))
        
          z_shear_cell = dbl_eps + &
            & SUM((ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x - ocean_state%p_diag%p_vn(jc,jk,blockNo)%x)**2)
          z_vert_density_grad_c(jc,jk,blockNo) = (z_rho_down(jk) - z_rho_up(jk-1)) *  &
            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,blockNo)
          z_ri_cell(jc, jk, blockNo) = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,blockNo) * z_grav_rho * &
            & (z_rho_down(jk) - z_rho_up(jk-1)) / z_shear_cell, 0.0_wp) ! do not use z_vert_density_grad_c,
                                                                     ! this is canceled out in this formula
        END DO ! levels
        ocean_state%p_diag%grad_rho_PP_vert(jc,2:levels,blockNo)=z_vert_density_grad_c(jc,2:levels,blockNo)
        ocean_state%p_diag%rho_GM(jc,1,blockNo)=ocean_state%p_diag%rho_GM(jc,2,blockNo)
 
      END DO ! index

      DO tracer_index = 1, no_tracer

        params_oce%a_tracer_v(start_index:end_index, 2:n_zlev, blockNo, tracer_index) =   &
          & MERGE(tracer_convection_MixingCoefficient,                    & ! activate convection
          & params_oce%a_tracer_v_back(tracer_index) +   & ! calculate the richardson diffusion
          &   tracer_RichardsonCoeff / ((1.0_wp + z_c1_t *    &
          &   z_ri_cell(start_index:end_index, 2:n_zlev, blockNo))**3), &
          & z_vert_density_grad_c(start_index:end_index, 2:n_zlev,blockNo) < convection_InstabilityThreshold)

        DO jc = start_index, end_index
          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
          DO jk = 2, levels

            IF (z_vert_density_grad_c(jc,jk,blockNo) < RichardsonDiffusion_threshold .AND. &
                z_vert_density_grad_c(jc,jk,blockNo) >= convection_InstabilityThreshold) THEN
              ! interpolate between convection and richardson diffusion
              diffusion_weight =  &
                & (z_vert_density_grad_c(jc,jk,blockNo) - convection_InstabilityThreshold) / &
                & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)
              params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = &
                & tracer_convection_MixingCoefficient * (1.0_wp - diffusion_weight) +&
                & diffusion_weight * params_oce%a_tracer_v(jc,jk,blockNo,tracer_index)
            ENDIF

          ENDDO ! levels
        ENDDO !  block index
      ENDDO ! tracer_index]

    END DO ! blocks
!ICON_OMP_END_DO

! !ICON_OMP_END_PARALLEL
!     IF (ltimer) CALL timer_stop(timer_extra10)
!     IF (ltimer) CALL timer_start(timer_extra11)
! !ICON_OMP_PARALLEL
    !--------------------------------------------
    ! Calculate params_oce%A_veloc_v:
    ! use mean values between the two cells; change to min, max if required
!ICON_OMP_DO PRIVATE(start_index, end_index, je, cell_1_idx, cell_1_block, cell_2_idx, cell_2_block, &
!ICON_OMP jk,  density_differ_edge, mean_z_r, diffusion_weight ) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      DO je = start_index, end_index

        cell_1_idx = patch_2D%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
        cell_2_idx = patch_2D%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)

        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
            ! TODO: the following expect equally sized cells
            ! compute density gradient at edges
!            density_differ_edge = 0.5_wp * &
!              & (z_vert_density_grad_c(cell_1_idx,jk,cell_1_block) + z_vert_density_grad_c(cell_2_idx,jk,cell_2_block))

            !! density gradient smaller then threshold ('semi-stable'): use background value
            ! note if density_differ_edge == z_threshold should be considered
!             IF     (density_differ_edge <  convection_InstabilityThreshold) THEN
!               ! turn on convection
!               params_oce%a_veloc_v(je,jk,blockNo) = max_vert_diff_veloc
!             ELSE
              ! richardson diffusion
              mean_z_r = 0.5_wp * (z_ri_cell(cell_1_idx,jk,cell_1_block) + z_ri_cell(cell_2_idx,jk,cell_2_block))
              params_oce%a_veloc_v(je,jk,blockNo) = &
                & params_oce%a_veloc_v_back +  &
                & z_av0 /                      &
                & ((1.0_wp + z_c1_v * mean_z_r)**2)

!               IF (density_differ_edge < RichardsonDiffusion_threshold) THEN
!                 diffusion_weight =  &
!                   & (density_differ_edge - convection_InstabilityThreshold) / &
!                   & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)
!
!                 params_oce%a_veloc_v(je,jk,blockNo) = &
!                   & max_vert_diff_veloc * (1.0_wp - diffusion_weight) + &
!                   & params_oce%a_veloc_v(je,jk,blockNo) * diffusion_weight
!
!                ENDIF

!             ENDIF

        END DO ! jk = 2, levels
      ENDDO ! je = start_index, end_index
    ENDDO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
!     IF (ltimer) CALL timer_stop(timer_extra11)

  END SUBROUTINE ICON_PP_scheme
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! As in the ICON_PP_scheme, but
  !! velocity gradients for the vertical viscocity are clculated on edges
  !!
  !! @par Revision History
  !! Initial release by Leonidas Linardakis, MPI-M (2011-02)
  !<Optimize:inUse:done>
  SUBROUTINE ICON_PP_Edge_scheme(patch_3d, ocean_state, params_oce) !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_ho_params), INTENT(inout)            :: params_oce

    ! Local variables
    INTEGER :: jc, blockNo, je,jk, tracer_index
    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
    INTEGER :: start_index, end_index
    INTEGER :: levels

    REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev), density(n_zlev)
    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
    REAL(wp) :: z_shear_cell
    REAL(wp) :: z_ri_cell               (nproma, n_zlev)
    REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)

    !Below is a set of variables and parameters for tracer and velocity
    REAL(wp), PARAMETER :: z_0               = 40.0_wp
    REAL(wp), PARAMETER :: z_c1_t            = 5.0_wp    !  PP diffusivity tuning constant
    REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp    !  PP viscosity tuning constant
    REAL(wp), PARAMETER :: z_threshold       = 5.0E-8_wp
    REAL(wp) :: diffusion_weight, instabilitySign
    REAL(wp) :: z_grav_rho, z_inv_OceanReferenceDensity
    REAL(wp) :: density_differ_edge, dz, richardson_edge, z_shear_edge
    REAL(wp), POINTER :: tracer_windMixing(:,:), velocity_windMixing(:,:)
    REAL(wp) :: z_vert_density_grad_e(n_zlev+1)
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D

    !-------------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    !cells_in_domain => patch_2D%cells%in_domain
    all_cells       => patch_2D%cells%ALL
    z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho
    levels = n_zlev

    !-------------------------------------------------------------------------
    z_grav_rho                   = grav/OceanReferenceDensity
    z_inv_OceanReferenceDensity                = 1.0_wp/OceanReferenceDensity
    !-------------------------------------------------------------------------
!     IF (ltimer) CALL timer_start(timer_extra10)

!ICON_OMP_PARALLEL PRIVATE(salinity, z_rho_up, z_rho_down)
    salinity(1:levels) = sal_ref
    z_rho_up(:)=0.0_wp
    z_rho_down(:)=0.0_wp
    
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, levels, jk, pressure, &
!ICON_OMP z_shear_cell, z_ri_cell, tracer_index, diffusion_weight, instabilitySign, tracer_windMixing) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      z_ri_cell(:,:) = 0.0_wp
      z_vert_density_grad_c(:,:, blockNo) = 0.0_wp
      DO jc = start_index, end_index

        levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
        IF (levels < 2) CYCLE

        IF(no_tracer >= 2) THEN
            salinity(1:levels) = ocean_state%p_prog(nold(1))%tracer(jc,1:levels,blockNo,2)
        ENDIF

        !--------------------------------------------------------
        ! pressure in dbars
        pressure(2:levels) = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, 2:levels, blockNo) * ReferencePressureIndbars
        z_rho_up(1:levels-1)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,1:levels-1,blockNo,1), &
          & salinity(1:levels-1), pressure(2:levels), levels-1)
        z_rho_down(2:levels)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,2:levels,blockNo,1), &
          & salinity(2:levels), pressure(2:levels), levels-1)

        DO jk = 2, levels
          ocean_state%p_diag%rho_GM(jc,jk,blockNo)=0.5_wp*(z_rho_up(jk)+z_rho_down(jk))
                    
          z_shear_cell = dbl_eps + &
            & SUM((ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x - ocean_state%p_diag%p_vn(jc,jk,blockNo)%x)**2)
          z_vert_density_grad_c(jc,jk,blockNo) = (z_rho_down(jk) - z_rho_up(jk-1)) *  &
            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,blockNo)
          z_ri_cell(jc, jk) = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,blockNo) * z_grav_rho * &
            & (z_rho_down(jk) - z_rho_up(jk-1)) / z_shear_cell, 0.0_wp) ! do not use z_vert_density_grad_c,
                                                                     ! this is canceled out in this formula
        END DO ! levels
        ocean_state%p_diag%grad_rho_PP_vert(jc,2:levels,blockNo)=z_vert_density_grad_c(jc,2:levels,blockNo)
        ocean_state%p_diag%rho_GM(jc,1,blockNo)=ocean_state%p_diag%rho_GM(jc,2,blockNo)
      END DO ! index
      !-----------------------------------------------------------
      tracer_windMixing => params_oce%tracer_windMixing(:,:,blockNo)
!       tracer_windMixing(:,:) = 0.0_wp
      IF (use_wind_mixing) THEN
        DO jc = start_index, end_index

          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
          ! wind-mixing: Marsland et al., 2003
          ! wind-mixing at surface, eq. (15) of Marsland et al., 2003
          tracer_windMixing(jc,1) = tracer_TopWindMixing  * WindAmplitude_at10m(jc,blockNo)**3 * &
            & (1.0_wp - SeaIceConcentration(jc,blockNo))

          ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
          DO jk = 2, levels
           tracer_windMixing(jc,jk) =  tracer_windMixing(jc,jk-1) * WindMixingDecay(jk) * &
             & WindMixingLevel(jk) / (WindMixingLevel(jk) + MAX(z_vert_density_grad_c(jc,jk,blockNo),0.0_wp))

          END DO! levels

        END DO ! index

      END IF  ! use_wind_mixing
      !-----------------------------------------------------------

      DO tracer_index = 1, no_tracer
        IF (tracer_index == 1) THEN
          instabilitySign = 0.0 ! always enable convection for temperature
        ELSE
          instabilitySign = Salinity_ConvectionRestrict
        ENDIF
        DO jc = start_index, end_index
          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
          DO jk = 2, levels
            ! by default use the richardson formula, no convection
            params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = MIN( &
              & params_oce%a_tracer_v_back(tracer_index) +   &
              & tracer_RichardsonCoeff / ((1.0_wp + z_c1_t * z_ri_cell(jc, jk))**3) + tracer_windMixing(jc,jk), &
              & tracer_convection_MixingCoefficient)
          IF (instabilitySign * &
                 (ocean_state%p_prog(nold(1))%tracer(jc,jk-1,blockNo,tracer_index) -         &
                  ocean_state%p_prog(nold(1))%tracer(jc,jk,blockNo,tracer_index)) >= 0.0_wp) & ! the = is important, do not change!
            THEN
              ! possibly convection
              IF (z_vert_density_grad_c(jc,jk,blockNo) <= convection_InstabilityThreshold) &
              THEN
                ! convection
                params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = tracer_convection_MixingCoefficient
              ELSE
                IF (z_vert_density_grad_c(jc,jk,blockNo) < RichardsonDiffusion_threshold) THEN
                  ! interpolate between convection and richardson diffusion
                  diffusion_weight =  &
                    & (z_vert_density_grad_c(jc,jk,blockNo) - convection_InstabilityThreshold) / &
                    & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)
                  params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = &
                    & tracer_convection_MixingCoefficient * (1.0_wp - diffusion_weight) +&
                    & diffusion_weight * params_oce%a_tracer_v(jc,jk,blockNo,tracer_index)
                ENDIF
              ENDIF
            ENDIF ! possibly convection

          ENDDO ! levels
        ENDDO !  block index
      ENDDO ! tracer_index]

    END DO ! blocks
!ICON_OMP_END_DO

! !ICON_OMP_END_PARALLEL
!     IF (ltimer) CALL timer_stop(timer_extra10)
!     IF (ltimer) CALL timer_start(timer_extra11)
! !ICON_OMP_PARALLEL
    !--------------------------------------------
    ! Calculate params_oce%A_veloc_v:
!ICON_OMP_DO PRIVATE(start_index, end_index, je, cell_1_idx, cell_1_block, cell_2_idx, cell_2_block, &
!ICON_OMP jk, dz, density_differ_edge, z_shear_edge, richardson_edge,z_vert_density_grad_e,velocity_windMixing) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      velocity_windMixing => params_oce%velocity_windMixing(:,:,blockNo)

      DO je = start_index, end_index

        cell_1_idx = patch_2D%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
        cell_2_idx = patch_2D%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)

        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
          z_vert_density_grad_e(jk) =  0.5_wp * &
            & (z_vert_density_grad_c(cell_1_idx,jk,cell_1_block) +   &
            &  z_vert_density_grad_c(cell_2_idx,jk,cell_2_block))
        ENDDO

        IF (use_wind_mixing) THEN
          velocity_windMixing(je, 1) =           &
            & velocity_TopWindMixing &
            & *  (0.5_wp * &
            &    (WindAmplitude_at10m(cell_1_idx,cell_1_block) + WindAmplitude_at10m(cell_2_idx,cell_2_block)))**3 &
            & * (1.0_wp - 0.5_wp *       &
            &    (SeaIceConcentration(cell_1_idx,cell_1_block) + SeaIceConcentration(cell_2_idx,cell_2_block)))

            ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
            DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
              velocity_windMixing(je, jk) =  velocity_windMixing(je, jk-1) * WindMixingDecay(jk) * &
                & WindMixingLevel(jk) / (WindMixingLevel(jk) + MAX(z_vert_density_grad_e(jk),0.0_wp))
            END DO! levels
        END IF  ! use_wind_mixing

        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
          ! TODO: the following expect equally sized cells
          ! compute density gradient at edges
          dz = 0.5_wp * (patch_3d%p_patch_1d(1)%prism_thick_e(je,jk-1,blockNo) + &
            &            patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,blockNo))
          density_differ_edge = z_vert_density_grad_e(jk) * dz
          z_shear_edge = dbl_eps + &
            & (ocean_state%p_prog(nold(1))%vn(je,jk,  blockNo) - &
            &  ocean_state%p_prog(nold(1))%vn(je,jk-1,blockNo)   )**2

          richardson_edge = MAX(dz * z_grav_rho * density_differ_edge / z_shear_edge, 0.0_wp)

          params_oce%a_veloc_v(je,jk,blockNo) =                 &
            & params_oce%a_veloc_v_back * dz +             &
            & velocity_RichardsonCoeff /                           &
            & ((1.0_wp + z_c1_v * richardson_edge)**2) +   &
            & velocity_windMixing(je, jk)

        END DO ! jk = 2, levels
      ENDDO ! je = start_index, end_index
    ENDDO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
!     IF (ltimer) CALL timer_stop(timer_extra11)

  END SUBROUTINE ICON_PP_Edge_scheme
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
  !! As in the ICON_PP_scheme, but
  !! velocity gradients for the vertical viscocity are clculated on edges
  !!
  !! @par Revision History
  !! Initial release by Leonidas Linardakis, MPI-M (2011-02)
  !<Optimize:inUse:done>
  SUBROUTINE ICON_PP_Edge_vnPredict_scheme(patch_3d, &
    & blockNo, start_index, end_index, ocean_state, vn_predict) !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    INTEGER, INTENT(in) :: blockNo, start_index, end_index
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    REAL(wp) :: vn_predict(:,:)

    ! Local variables
    INTEGER :: je,jk
    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
    INTEGER :: levels

    !Below is a set of variables and parameters for tracer and velocity
    REAL(wp), PARAMETER :: z_0               = 40.0_wp
    REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp    !  PP viscosity tuning constant
    REAL(wp) :: diffusion_weight
    REAL(wp) :: z_grav_rho, z_inv_OceanReferenceDensity
    REAL(wp) :: density_differ_edge, dz, richardson_edge, z_shear_edge, vn_diff, new_velocity_friction
    !-------------------------------------------------------------------------
    REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)
    REAL(wp), POINTER :: wind_mixing(:,:)
    REAL(wp) :: z_vert_density_grad_e(1:n_zlev+1)
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_ho_params), POINTER :: params_oce

    !-------------------------------------------------------------------------
    params_oce      => v_params
    patch_2D        => patch_3d%p_patch_2d(1)
    z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho ! already calculated
    wind_mixing     => params_oce%velocity_windMixing(:,:,blockNo) ! already calculated
    levels = n_zlev

    !-------------------------------------------------------------------------
    z_grav_rho                   = grav/OceanReferenceDensity
    z_inv_OceanReferenceDensity  = 1.0_wp/OceanReferenceDensity
    !-------------------------------------------------------------------------

    DO je = start_index, end_index

      cell_1_idx = patch_2D%edges%cell_idx(je,blockNo,1)
      cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
      cell_2_idx = patch_2D%edges%cell_idx(je,blockNo,2)
      cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)

      DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
        z_vert_density_grad_e(jk) =  0.5_wp * &
          & (z_vert_density_grad_c(cell_1_idx,jk,cell_1_block) +   &
          &  z_vert_density_grad_c(cell_2_idx,jk,cell_2_block))
      ENDDO

      DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
        ! TODO: the following expect equally sized cells
        ! compute density gradient at edges
        dz = 0.5_wp * (patch_3d%p_patch_1d(1)%prism_thick_e(je,jk-1,blockNo) + &
          &            patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,blockNo))
        density_differ_edge = z_vert_density_grad_e(jk) * dz
!         vn_diff = MAX( &
!           & ABS(ocean_state%p_prog(nold(1))%vn(je,jk,  blockNo) - &
!           &      ocean_state%p_prog(nold(1))%vn(je,jk-1,blockNo)), &
!           & ABS(vn_predict(je,jk) - &
!           &     vn_predict(je,jk-1)))
        vn_diff =  &
          & ABS(vn_predict(je,jk) - &
          &     vn_predict(je,jk-1))

        z_shear_edge = dbl_eps + vn_diff**2

        richardson_edge = MAX(dz * z_grav_rho * density_differ_edge / z_shear_edge, 0.0_wp)

        new_velocity_friction = &
          & params_oce%a_veloc_v_back * dz +                              &
          & velocity_RichardsonCoeff / ((1.0_wp + z_c1_v * richardson_edge)**2)+  &
          & wind_mixing(je,jk)

        ! the average of the calculated velocity friction based on the old velocity and the predicted one
        params_oce%a_veloc_v(je,jk,blockNo) = &
!           & MAX(params_oce%a_veloc_v(je,jk,blockNo), new_velocity_friction )
           & VerticalViscosity_TimeWeight * params_oce%a_veloc_v(je,jk,blockNo) + &
           & (1.0_wp - VerticalViscosity_TimeWeight) * new_velocity_friction

      END DO ! jk = 2, levels
    ENDDO ! je = start_index, end_index

  END SUBROUTINE ICON_PP_Edge_vnPredict_scheme
  !-------------------------------------------------------------------------

!   !-------------------------------------------------------------------------
!   !>
!   !! As in the ICON_PP_scheme, but
!   !! velocity gradients for the vertical viscocity are calculated on edges
!   !!
!   !! @par Revision History
!   !! Initial release by Leonidas Linardakis, MPI-M (2011-02)
!   !<Optimize:inUse:done>
!   SUBROUTINE update_PhysicsParameters_ICON_PP_Tracer(patch_3d, ocean_state) !, calculate_density_func)
!
!     TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
!     TYPE(t_hydro_ocean_state), TARGET :: ocean_state
!
!     ! Local variables
!     INTEGER :: jc, blockNo, je,jk, tracer_index
!     !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
!     INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
!     INTEGER :: start_index, end_index
!     INTEGER :: levels
!
!     REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev), density(n_zlev)
!     REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
!     REAL(wp) :: z_shear_cell
!     REAL(wp) :: z_ri_cell(nproma, n_zlev+1)
!     REAL(wp) :: wind_mixing(nproma, n_zlev+1)
!     REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)
!
!     !Below is a set of variables and parameters for tracer and velocity
!     REAL(wp), PARAMETER :: z_0               = 40.0_wp
!     REAL(wp), PARAMETER :: z_c1_t            = 5.0_wp    !  PP diffusivity tuning constant
!     REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp    !  PP viscosity tuning constant
!     REAL(wp), PARAMETER :: z_threshold       = 5.0E-8_wp
!     REAL(wp) :: diffusion_weight
!     REAL(wp) :: z_grav_rho, z_inv_OceanReferenceDensity
!     REAL(wp) :: density_differ_edge, dz, richardson_edge, z_shear_edge
!     !-------------------------------------------------------------------------
!     TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
!     TYPE(t_patch), POINTER :: patch_2D
!     TYPE(t_ho_params), POINTER  :: params_oce
!
!     !-------------------------------------------------------------------------
!     params_oce      => v_params
!     patch_2D        => patch_3d%p_patch_2d(1)
!     edges_in_domain => patch_2D%edges%in_domain
!     !cells_in_domain => patch_2D%cells%in_domain
!     all_cells       => patch_2D%cells%ALL
!     z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho
!     levels = n_zlev
!
!     !-------------------------------------------------------------------------
!     z_grav_rho                   = grav/OceanReferenceDensity
!     z_inv_OceanReferenceDensity                = 1.0_wp/OceanReferenceDensity
!     !-------------------------------------------------------------------------
! !     IF (ltimer) CALL timer_start(timer_extra10)
!
! !ICON_OMP_PARALLEL PRIVATE(salinity)
!     salinity(1:levels) = sal_ref
! !ICON_OMP_DO PRIVATE(start_index, end_index, jc, levels, jk, pressure, z_rho_up, z_rho_down, &
! !ICON_OMP z_shear_cell, z_ri_cell, tracer_index, diffusion_weight, wind_mixing) ICON_OMP_DEFAULT_SCHEDULE
!     DO blockNo = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, blockNo, start_index, end_index)
!       z_ri_cell(:,:) = 0.0_wp
!       z_vert_density_grad_c(:,:, blockNo) = 0.0_wp
!       DO jc = start_index, end_index
!
!         levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
!         IF (levels < 2) CYCLE
!
!         IF(no_tracer >= 2) THEN
!             salinity(1:levels) = ocean_state%p_prog(nold(1))%tracer(jc,1:levels,blockNo,2)
!         ENDIF
!
!         !--------------------------------------------------------
!         pressure(2:levels) = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, 2:levels, blockNo) * OceanReferenceDensity * sitodbar
!         z_rho_up(1:levels-1)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,1:levels-1,blockNo,1), &
!           & salinity(1:levels-1), pressure(2:levels), levels-1)
!         z_rho_down(2:levels)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,2:levels,blockNo,1), &
!           & salinity(2:levels), pressure(2:levels), levels-1)
!
!         DO jk = 2, levels
!           z_shear_cell = dbl_eps + &
!             & SUM((ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x - ocean_state%p_diag%p_vn(jc,jk,blockNo)%x)**2)
!           z_vert_density_grad_c(jc,jk,blockNo) = (z_rho_down(jk) - z_rho_up(jk-1)) *  &
!             & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,blockNo)
!           z_ri_cell(jc, jk) = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,blockNo) * z_grav_rho * &
!             & (z_rho_down(jk) - z_rho_up(jk-1)) / z_shear_cell, 0.0_wp) ! do not use z_vert_density_grad_c,
!                                                                         ! this is canceled out in this formula
!         END DO ! levels
!
!       END DO ! index
!       !-----------------------------------------------------------
!       wind_mixing(:,:) = 0.0_wp
!       IF (use_wind_mixing) THEN
!         DO jc = start_index, end_index
!
!           levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
!           ! wind-mixing: Marsland et al., 2003
!           ! wind-mixing at surface, eq. (15) of Marsland et al., 2003
!           wind_mixing(jc,1) = tracer_TopWindMixing  * WindAmplitude_at10m(jc,blockNo)**3 * &
!             & (1.0_wp - SeaIceConcentration(jc,blockNo))
!
!           ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
!           DO jk = 2, levels
!            wind_mixing(jc,jk) =  wind_mixing(jc,jk-1) * WindMixingDecay(jk) * &
!              & WindMixingLevel(jk) / (WindMixingLevel(jk) + MAX(z_vert_density_grad_c(jc,jk,blockNo),0.0_wp))
!
!           END DO! levels
!
!         END DO ! index
!
!       END IF  ! use_wind_mixing
!       !-----------------------------------------------------------
!
!       DO tracer_index = 1, no_tracer
!
!         params_oce%a_tracer_v(start_index:end_index, 2:n_zlev, blockNo, tracer_index) =   &
!           & MERGE(                                       &
!           & tracer_convection_MixingCoefficient,                          & ! activate convection
!           & params_oce%a_tracer_v_back(tracer_index) +   & ! calculate the richardson diffusion
!           &   tracer_RichardsonCoeff / ((1.0_wp + z_c1_t *    &
!           &   z_ri_cell(start_index:end_index, 2:n_zlev))**3) + &
!           &   wind_mixing(start_index:end_index, 2:n_zlev), &
!           & z_vert_density_grad_c(start_index:end_index, 2:n_zlev,blockNo) <= convection_InstabilityThreshold)
!
!         DO jc = start_index, end_index
!           levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
!           DO jk = 2, levels
!
!             IF (z_vert_density_grad_c(jc,jk,blockNo) < RichardsonDiffusion_threshold .AND. &
!                 z_vert_density_grad_c(jc,jk,blockNo) > convection_InstabilityThreshold) THEN
!               ! interpolate between convection and richardson diffusion
!               diffusion_weight =  &
!                 & (z_vert_density_grad_c(jc,jk,blockNo) - convection_InstabilityThreshold) / &
!                 & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)
!               params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = &
!                 & tracer_convection_MixingCoefficient * (1.0_wp - diffusion_weight) +&
!                 & diffusion_weight * params_oce%a_tracer_v(jc,jk,blockNo,tracer_index)
!             ENDIF
!
!           ENDDO ! levels
!         ENDDO !  block index
!       ENDDO ! tracer_index]
!
!     END DO ! blocks
! !ICON_OMP_END_DO NOWAIT
! !ICON_OMP_END_PARALLEL
!
!   END SUBROUTINE update_PhysicsParameters_ICON_PP_Tracer
!   !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Update of ocean physics parameters
  !!
  !! Update of ocean physics: This routine is used used only if time-dependent
  !! changes of physical parametrizations.
  !! Currently vertical mixing coefficients for tracers and vertical diffusivity are updated.
  !! Dependent on the local Richardson number the diffusivity are calculated
  !! (Large & Gent JPO 29, (1999), 449-464).
  !! The formulation follows the MPI-OM implementation as described in Marsland et al. (Ocean
  !! Modelling 5, 2003).
  !! The notational convention is also taken from this paper( cf. eqs (14) and (19)).
  !! What is missing is the fractional ice cover (see eqs. (15-16)).
  !! Eq. (18) is the Redi part that is not implemented, yet
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-02)
!<Optimize:inUse>
  SUBROUTINE MPIOM_PP_scheme(patch_3d, ocean_state, fu10, concsum, params_oce) !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET     :: ocean_state
    REAL(wp),          INTENT(in)         :: fu10   (nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) ! t_atmos_for_ocean%fu10
    REAL(wp),          INTENT(in)         :: concsum(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) ! t_sea_ice%concsum
    TYPE(t_ho_params), INTENT(inout)      :: params_oce
!     INTERFACE !This contains the function version of the actual EOS as chosen in namelist
!       FUNCTION calculate_density_func(tpot, sal, press) result(rho)
!         USE mo_kind, ONLY: wp
!         REAL(wp), INTENT(in) :: tpot
!         REAL(wp), INTENT(in) :: sal
!         REAL(wp), INTENT(in) :: press
!         REAL(wp) :: rho
!       ENDFUNCTION calculate_density_func
!     END INTERFACE

    ! Local variables
    INTEGER :: jc, blockNo, je,jk, tracer_index
    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
    INTEGER :: start_index, end_index
    INTEGER :: levels, jk_max

    REAL(wp) :: rho_up(n_zlev), rho_down(n_zlev)
    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
    REAL(wp) :: vert_velocity_shear
    REAL(wp), POINTER :: vert_density_grad(:,:,:)

    ! Local 3dim variables
    REAL(wp) :: richardson_no(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: dv_wind      (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: av_wind      (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)

    ! Below is a set of variables and parameters for diffusion of tracer and velocity
    ! lambda_diff: relne in MPIOM (0.4), Lambda_D (0.6) in Marsland et al. (2003), same vor Lambda_v
    REAL(wp), PARAMETER :: lambda_diff       = 0.4_wp    !  eddy diffusion relaxation constant
    REAL(wp), PARAMETER :: lambda_visc       = 0.4_wp    !  eddy viscosity relaxation constant
    REAL(wp), PARAMETER :: z0_wind           = 40.0_wp   !  exponential decay of wind mixing with depth
    REAL(wp), PARAMETER :: v10m_ref          = 6.0_wp    !  wind mixing 10m reference windspeed
    REAL(wp), PARAMETER :: crd               = 5.0_wp    !  PP diffusivity tuning constant
    REAL(wp), PARAMETER :: crv               = 5.0_wp    !  PP viscosity tuning constant

    REAL(wp) :: decay_wind_depth, wind_param, vdensgrad_inter, densgrad_k, densgrad_kp1
    REAL(wp) :: v10mexp_3, wma_pv, wma_pd
    REAL(wp) :: diffusion_weight, loc_eps
    REAL(wp) :: dv_old, dv_back, dv_rich
    REAL(wp) :: av_old, av_back, av_rich
    REAL(wp) :: onem1_lambda_d, onem1_lambda_v
    REAL(wp) :: grav_rho
    REAL(wp) :: mean_density_differ_edge, mean_richardson_e, fu10_e, conc_e
    REAL(wp) :: vdgfac_bot(n_zlev)
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
    TYPE(t_patch), POINTER :: p_patch

    !-------------------------------------------------------------------------
    p_patch         => patch_3d%p_patch_2d(1)
    edges_in_domain => p_patch%edges%in_domain
    !cells_in_domain => p_patch%cells%in_domain
    all_cells       => p_patch%cells%ALL
    vert_density_grad => ocean_state%p_diag%zgrad_rho
    levels = n_zlev
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! Attention: with use_constant_mixing=.true. there is no application of
    ! convective mixing parameters in case of instability
    ! max_vert_diff_veloc / tracer_convection_MixingCoefficient
    ! control of convective and constant mixing should be independent

    grav_rho          = grav/OceanReferenceDensity

    onem1_lambda_d    = 1.0_wp-lambda_diff
    onem1_lambda_v    = 1.0_wp-lambda_visc
    dv_rich           = tracer_RichardsonCoeff
    av_rich           = velocity_RichardsonCoeff
    v10mexp_3         = 1.0_wp/v10m_ref**3
    wma_pd            = wma_diff * v10mexp_3   !  scaled wind-mixing amplitude for diffusion
    wma_pv            = wma_visc * v10mexp_3   !  scaled wind-mixing amplitude for viscosity

    dv_wind(:,1:levels,:) = 0.0_wp
    av_wind(:,1:levels,:) = 0.0_wp

    loc_eps = dbl_eps

!ICON_OMP_PARALLEL PRIVATE(salinity, rho_up, rho_down)
    salinity(1:levels) = sal_ref
    rho_up(:)=0.0_wp    
    rho_down(:)=0.0_wp
    
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, levels, jk, pressure, rho_up, rho_down, &
!ICON_OMP vert_velocity_shear, tracer_index, diffusion_weight, decay_wind_depth, wind_param, &
!ICON_OMP jk_max, vdgfac_bot, vdensgrad_inter, dv_old, dv_back) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      richardson_no    (:,:, blockNo) = 0.0_wp
      vert_density_grad(:,:, blockNo) = 0.0_wp
      DO jc = start_index, end_index

        levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
        IF (levels < min_dolic) CYCLE

        IF(no_tracer >= 2) THEN
            salinity(1:levels) = ocean_state%p_prog(nold(1))%tracer(jc,1:levels,blockNo,2)
        ENDIF

        !--------------------------------------------------------
        ! calculate density gradient for stability detection for PP-scheme and convection parameterization:
        !  - pressure at interface between upper and lower layer
        !  - S and T taken from upper and lower level, i.e. 2 times calculation of density per layer

        !  - old formulation without including z in reference interface level
        ! pressure in dbars
        pressure(2:levels) = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, 2:levels, blockNo) * ReferencePressureIndbars
        rho_up(1:levels-1)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,1:levels-1,blockNo,1), &
          & salinity(1:levels-1), pressure(2:levels), levels-1)
        rho_down(2:levels)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,2:levels,blockNo,1), &
          & salinity(2:levels), pressure(2:levels), levels-1)


        DO jk = 2, levels
        
          ocean_state%p_diag%rho_GM(jc,jk,blockNo)=0.5_wp*(rho_up(jk)+rho_down(jk))
        
          ! division by dz**2 is omitted in this calculation of velocity shear: shear = (d_vn)**2
          vert_velocity_shear = loc_eps + &
            & SUM((ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x - ocean_state%p_diag%p_vn(jc,jk,blockNo)%x)**2)
          ! d_rho/dz - full density gradient necessary for wind mixing, stability formula, mixed layer depth calculation
          vert_density_grad(jc,jk,blockNo) = (rho_down(jk) - rho_up(jk-1)) *  &
            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,blockNo)

          ocean_state%p_diag%grad_rho_PP_vert(jc,jk,blockNo)=vert_density_grad(jc,jk,blockNo)


          ! Ri = g/OceanReferenceDensity * dz * d_rho/(d_vn)**2
          richardson_no(jc,jk,blockNo) = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,blockNo) * grav_rho * &
            &                           (rho_down(jk) - rho_up(jk-1)) / vert_velocity_shear, 0.0_wp)
        END DO ! levels
!        ocean_state%p_diag%grad_rho_PP_vert(jc,1,blockNo)=0.0_wp
!        ocean_state%p_diag%grad_rho_PP_vert(jc,levels+1:n_zlev,blockNo)=0.0_wp
        ocean_state%p_diag%rho_GM(jc,1,blockNo)=ocean_state%p_diag%rho_GM(jc,2,blockNo)
!FIXME : this is not done in the other pp schemes
        ocean_state%p_diag%rho_GM(jc,levels+1:n_zlev,blockNo)=ocean_state%p_diag%rho_GM(jc,levels,blockNo)
        
        IF (use_wind_mixing) THEN

          ! wind-mixing at surface, eq. (15) of Marsland et al., 2003

          ! reduced wind-mixing under sea ice, following MPIOM
          IF (use_reduced_mixing_under_ice) THEN
            dv_wind(jc,1,blockNo) = wma_pv * (1.0_wp - concsum(jc,blockNo)) * fu10(jc,blockNo)**3
          ELSE
            dv_wind(jc,1,blockNo) = wma_pv * (1.0_wp - concsum(jc,blockNo))**2 * fu10(jc,blockNo)**3
          ENDIF

          ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
          DO jk = 2, levels

            decay_wind_depth   = EXP(-patch_3d%p_patch_1d(1)%del_zlev_m(jk-1)/z0_wind)

            ! lambda_wind: default changed to 0.05, with 0.03 omip-r2b4 aborted
            !   - in MPIOM it is 0.05, in Marsland et al. it was 0.03,
            !   - for strong winds it might be necessary to increase it
            wind_param         = lambda_wind * patch_3d%p_patch_1d(1)%inv_del_zlev_m(jk)

            ! vertical interpolation of density gradient
            !  - at mid-depth jk the density gradients at interfaces from above (jk) and below (jk+1) are used
            !  - at bottom level the density gradient from below is zero, half density gradient is used
            jk_max             = MIN(jk+1,levels)
            vdgfac_bot(jk)     = 1.0_wp
            vdgfac_bot(levels) = 0.0_wp
            vdensgrad_inter  = 0.5_wp*(vert_density_grad(jc,jk,blockNo)+vdgfac_bot(jk)*vert_density_grad(jc,jk_max,blockNo))

            dv_wind(jc,jk,blockNo) =  dv_wind(jc,jk-1,blockNo)*decay_wind_depth*wind_param / (wind_param + vdensgrad_inter)

            ! cut unphysically negative wind induced mixing
            dv_wind(jc,jk,blockNo) =  MAX(dv_wind(jc,jk,blockNo),0.0_wp)
            ! cut overshoots more than convection maximum - must not be set to low values
          ! IF (tracer_convection_MixingCoefficient .GT. params_oce%a_tracer_v_back(1)) &
          !   &  dv_wind(jc,jk,blockNo) =  MIN(dv_wind(jc,jk,blockNo),0.0_wp)

          END DO

        END IF  ! use_wind_mixing

      ENDDO !  block index

      DO tracer_index = 1, no_tracer
        DO jc = start_index, end_index
          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)

          DO jk = 2, levels

            ! calculate the richardson diffusion using the eddy diffusion relaxation term lambda_diff
            !  - a_tracer_v is relaxed to the last pp-value to avoid relaxing to convection value tracer_convection_MixingCoefficient

              dv_old  = params_oce%a_tracer_v(jc,jk,blockNo,tracer_index)
              dv_back = params_oce%a_tracer_v_back(tracer_index)
              params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = &
                &    onem1_lambda_d*MIN(dv_old, dv_rich + dv_wind(jc,jk,blockNo)) &
                &  + lambda_diff*(dv_rich/((1.0_wp + crd*richardson_no(jc,jk,blockNo))**3) + dv_back + dv_wind(jc,jk,blockNo))


         !  ! clalculate the richardson diffusion - old arithmetic
         !  IF (use_pp_scheme) THEN
         !    params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = &
         !      & params_oce%a_tracer_v_back(tracer_index) + &
         !      & dv_rich / ((1.0_wp + crd *                 &
         !      & richardson_no(jc,jk,blockNo))**3)
         !  ENDIF

                ! #slo# ensure that pp is active for low values of tracer_convection_MixingCoefficient
                dv_old = params_oce%a_tracer_v(jc,jk,blockNo,tracer_index)
                params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = MAX(  &
                  ! #slo# Attention: convection_InstabilityThreshold<0 in old formulation - used with reverted sign
                  &  tracer_convection_MixingCoefficient * &
                  &  (-convection_InstabilityThreshold-vert_density_grad(jc,jk,blockNo)) /     &
                  &  (-convection_InstabilityThreshold+ABS(vert_density_grad(jc,jk,blockNo))), dv_old)

          ENDDO ! levels
        ENDDO !  block index
      ENDDO ! tracer_index

    END DO ! blocks
!ICON_OMP_END_DO


    !--------------------------------------------
    ! Calculate params_oce%A_veloc_v:
    ! use mean values between the two cells; change to min, max if required
!ICON_OMP_DO PRIVATE(start_index, end_index, je, cell_1_idx, cell_1_block, cell_2_idx, cell_2_block, &
!ICON_OMP levels, jk,  mean_density_differ_edge, mean_richardson_e, decay_wind_depth, wind_param, &
!ICON_OMP jk_max, vdgfac_bot, densgrad_k, densgrad_kp1, vdensgrad_inter, av_old, av_back) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      DO je = start_index, end_index

        levels = patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
        IF (levels < min_dolic) CYCLE

        cell_1_idx = p_patch%edges%cell_idx(je,blockNo,1)
        cell_1_block = p_patch%edges%cell_blk(je,blockNo,1)
        cell_2_idx = p_patch%edges%cell_idx(je,blockNo,2)
        cell_2_block = p_patch%edges%cell_blk(je,blockNo,2)

        IF (use_wind_mixing) THEN

          ! TODO: the following expects equally sized cells
          ! TODO: think about boundary values on land
          fu10_e = 0.5_wp * (   fu10(cell_1_idx,cell_1_block) +    fu10(cell_2_idx,cell_2_block))
          !fu10_e = 10.0_wp
          conc_e = 0.5_wp * (concsum(cell_1_idx,cell_1_block) + concsum(cell_2_idx,cell_2_block))

          ! wind-mixing at surface, eq. (15) of Marsland et al., 2003

          ! reduced wind-mixing under sea ice, following MPIOM
          IF (use_reduced_mixing_under_ice) THEN
            av_wind(je,1,blockNo) = wma_pd * (1.0_wp - conc_e)    * fu10_e**3
          ELSE
            av_wind(je,1,blockNo) = wma_pd * (1.0_wp - conc_e)**2 * fu10_e**3
          ENDIF

          ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003, edges
          DO jk = 2, levels
            decay_wind_depth = EXP(-patch_3d%p_patch_1d(1)%del_zlev_m(jk-1)/z0_wind)
            wind_param       = lambda_wind * patch_3d%p_patch_1d(1)%inv_del_zlev_m(jk)

            ! vertical interpolation of density gradient
            !  - at mid-depth jk the density gradients at interfaces from above (jk) and below (jk+1) are used
            !  - at bottom level the density gradient from below is zero, half density gradient is used
            jk_max             = MIN(jk+1,levels)
            vdgfac_bot(jk)     = 1.0_wp
            vdgfac_bot(levels) = 0.0_wp
     !      vdensgrad_inter    = 0.5_wp*(vert_density_grad(jc,jk,blockNo)+vdgfac_bot(jk)*vert_density_grad(jc,jk_max,blockNo))
            densgrad_k       = 0.5_wp * (vert_density_grad(cell_1_idx,jk,cell_1_block) + &
              & vert_density_grad(cell_2_idx,jk,cell_2_block))
            densgrad_kp1     = 0.5_wp*vdgfac_bot(jk) * (vert_density_grad(cell_1_idx,jk_max,cell_1_block) &
              & + vert_density_grad(cell_2_idx,jk_max,cell_2_block))
            vdensgrad_inter  = 0.5_wp*(densgrad_k + densgrad_kp1)

            av_wind(je,jk,blockNo) =  av_wind(je,jk-1,blockNo)*decay_wind_depth*wind_param / (wind_param + vdensgrad_inter)

            ! cut unphysically negative wind induced mixing
            av_wind(je,jk,blockNo) =  MAX(av_wind(je,jk,blockNo),0.0_wp)
            ! cut overshoots more than convection maximum - must not be set to low values
    !       IF (max_vert_diff_veloc .GT. params_oce%a_veloc_v_back) &
    !         &  av_wind(je,jk,blockNo) =  MIN(av_wind(je,jk,blockNo),max_vert_diff_veloc)
          END DO

        END IF  ! use_wind_mixing

        DO jk = 2, levels

          ! Set to zero for land + boundary locations edges
!           IF (patch_3d%lsm_e(je,jk,blockNo) > sea) THEN
!             params_oce%a_veloc_v(je,jk,blockNo) = 0.0_wp
!           ELSE

          ! TODO: the following expects equally sized cells
          ! compute quantities at edges
          mean_density_differ_edge = 0.5_wp * (vert_density_grad(cell_1_idx,jk,cell_1_block) &
            & + vert_density_grad(cell_2_idx,jk,cell_2_block))
          mean_richardson_e   = 0.5_wp * (richardson_no(cell_1_idx,jk,cell_1_block) + &
            & richardson_no(cell_2_idx,jk,cell_2_block))

          ! calculate the richardson viscosity using the eddy viscosity relaxation term lambda_visc
          !  - a_veloc_v is relaxed to the last pp-value to avoid relaxing to convection value max_vert_diff_veloc
!           IF (use_pp_scheme) THEN

            av_old  = params_oce%a_veloc_v(je,jk,blockNo)
            av_back = params_oce%a_veloc_v_back
            params_oce%a_veloc_v(je,jk,blockNo) = &
              &    onem1_lambda_v*MIN(av_old, av_rich + av_wind(je,jk,blockNo)) &
              &  + lambda_visc*(av_rich/((1.0_wp + crv*mean_richardson_e)**2) + av_back + av_wind(je,jk,blockNo))
!           ENDIF

          ! turn on convection - no convection for velocity as in mpiom
!           IF (use_convection) THEN
!
!             ! MPIOM style of convection in PP-scheme: viscosity
!             IF (use_mpiom_pp_form) THEN
!               ! #slo# ensure that pp is active for low values of tracer_convection_MixingCoefficient
!               av_old = params_oce%a_veloc_v(je,jk,blockNo)
!               params_oce%a_veloc_v(je,jk,blockNo) = MAX(  &
!                 ! #slo# Attention: convection_InstabilityThreshold<0 in old formulation - used with reverted sign
!                 &  max_vert_diff_veloc * (-convection_InstabilityThreshold-mean_density_differ_edge) /     &
!                 &                       (-convection_InstabilityThreshold+ABS(mean_density_differ_edge)), av_old)
!
!             ELSE ! do not use_mpiom_pp_form
!
!               ! turn on convection
!               IF (mean_density_differ_edge <  convection_InstabilityThreshold) THEN
!                 ! #slo# ensure that pp is active for low values of tracer_convection_MixingCoefficient
!                 !params_oce%a_veloc_v(je,jk,blockNo) = max_vert_diff_veloc
!                 params_oce%a_veloc_v(je,jk,blockNo) = MAX(max_vert_diff_veloc,params_oce%a_veloc_v(je,jk,blockNo))
!               ELSE

!                 IF (mean_density_differ_edge < RichardsonDiffusion_threshold) THEN
!                   diffusion_weight =  &
!                     & (mean_density_differ_edge - convection_InstabilityThreshold) / &
!                     & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)
!
!                   ! richardson diffusion from above
!                   av_old = params_oce%a_veloc_v(je,jk,blockNo)
!                   params_oce%a_veloc_v(je,jk,blockNo) = &
!                     & max_vert_diff_veloc * (1.0_wp - diffusion_weight) + &
!                     & av_old * diffusion_weight
!
!                 ENDIF  ! grad<RichThreshold
!               ENDIF ! grad<convThreshold
!             ENDIF ! use_mpiom_pp_form
!           ENDIF ! use_convection

        END DO ! jk = 2, levels
      ENDDO ! je = start_index, end_index
    ENDDO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL

    ! Sync the results, the A_tracer_v is only for checking
    ! DO tracer_index = 1, no_tracer
    !   CALL sync_patch_array(sync_c,p_patch,params_oce%a_tracer_v(:,:,:,tracer_index))
    ! END DO
    ! CALL sync_patch_array(sync_e,p_patch,params_oce%a_veloc_v(:,:,:))

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print levels (1-5, fix)
    CALL dbg_print('UpdPar: p_vn%x(1)    ',ocean_state%p_diag%p_vn%x(1),module_name,idt_src,in_subset=p_patch%cells%owned)
  ! CALL dbg_print('UpdPar: p_vn%x(2)    ',ocean_state%p_diag%p_vn%x(2),module_name,idt_src,in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdPar: windMix Diff ',dv_wind                     ,module_name,idt_src,in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdPar: windMix Visc ',av_wind                     ,module_name,idt_src,in_subset=p_patch%edges%owned)
    CALL dbg_print('UpdPar: VertDensGrad ',vert_density_grad           ,module_name,idt_src,in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdPar: Richardson No',richardson_no               ,module_name,idt_src,in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdPar: windsp. fu10 ',fu10                        ,module_name,idt_src,in_subset=p_patch%cells%owned)
    idt_src=5  ! output print levels (1-5, fix)
!     DO tracer_index = 1, no_tracer
!       CALL dbg_print('UpdPar FinalTracerMixing'  ,params_oce%a_tracer_v(:,:,:,tracer_index), module_name,idt_src, &
!         & in_subset=p_patch%cells%owned)
!     ENDDO
!     CALL dbg_print('UpdPar FinalVelocMixing'   ,params_oce%a_veloc_v     ,module_name,idt_src, &
!       & in_subset=p_patch%edges%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE MPIOM_PP_scheme
  !-------------------------------------------------------------------------


END MODULE mo_ocean_pp_scheme
