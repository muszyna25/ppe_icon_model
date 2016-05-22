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
MODULE mo_ocean_physics
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_ocean_nml,           ONLY: n_zlev, bottom_drag_coeff, k_veloc_h, k_veloc_v,        &
    & k_pot_temp_h, k_pot_temp_v, k_sal_h, k_sal_v, no_tracer,&
    & max_vert_diff_veloc, max_vert_diff_trac,                &
    & HorizontalViscosity_type, veloc_diffusion_order,        &
    & n_points_in_munk_layer,                                 &
    & HorizontalViscosityBackground_Biharmonic,               &
    & richardson_tracer, richardson_veloc,                    &
    & physics_parameters_type,                                &
    & physics_parameters_Constant_type,                       &
    & physics_parameters_ICON_PP_type,                        &
    & physics_parameters_ICON_PP_Edge_type,                   &
    & physics_parameters_ICON_PP_Edge_vnPredict_type,         &
    & physics_parameters_MPIOM_PP_type,                       &
    & use_wind_mixing,                                        &
    & HorizontalViscosity_SmoothIterations,                   &
    & convection_InstabilityThreshold,                        &
    & RichardsonDiffusion_threshold,                          &
    & use_convection_parameterization,                        &
    & lambda_wind, wma_diff, wma_visc,                        &
    & use_reduced_mixing_under_ice,                           &
    & k_tracer_dianeutral_parameter,                          &
    & k_tracer_isoneutral_parameter, k_tracer_GM_kappa_parameter,    &
    & GMRedi_configuration,GMRedi_combined,                   &
    & GM_only,Redi_only,                                      &
    & leith_closure, leith_closure_gamma,                     & 
    & veloc_diffusion_form, biharmonic_const,                 &
    & HorizontalViscosity_SpatialSmoothFactor,                &
    & VerticalViscosity_TimeWeight, OceanReferenceDensity,    &
    & HorizontalViscosity_ScaleWeight,                        &
    & tracer_TopWindMixing, WindMixingDecayDepth,             &
    & velocity_TopWindMixing, Salinity_ConvectionRestrict
    
   !, l_convection, l_pp_scheme
  USE mo_parallel_config,     ONLY: nproma
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_impl_constants,      ONLY: success, max_char_length, min_dolic, sea
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_util_dbg_prnt,       ONLY: dbg_print, debug_print_MaxMinMean
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state, t_onEdges_Pointer_3d_wp, t_onCells_HalfLevels_Pointer_wp, t_operator_coeff
  USE mo_ocean_state,         ONLY: oce_config
  USE mo_physical_constants,  ONLY: grav, sitodbar,sal_ref
  USE mo_math_constants,      ONLY: dbl_eps, pi, rad2deg
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
  USE mo_cdi_constants,       ONLY: grid_cell, grid_edge,            &
    & grid_unstructured_edge, grid_unstructured_cell, &
    & za_depth_below_sea, za_depth_below_sea_half
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, global_max
  USE  mo_ocean_thermodyn,      ONLY: calculate_density_onColumn
  USE mo_ocean_math_operators,  ONLY: div_oce_3d
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, timer_upd_phys, &
    & timer_extra10, timer_extra11
  USE mo_statistics,          ONLY: global_minmaxmean
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_ocean_pp_scheme,     ONLY: update_physics_parameters_ICON_PP_scheme,&
    &                               update_physics_parameters_ICON_PP_Edge_vnPredict_scheme,&
    &                               update_physics_parameters_ICON_PP_Edge_scheme,&
    &                               update_physics_parameters_ICON_PP_Tracer,&
    &                               update_physics_parameters_MPIOM_PP_scheme
  USE mo_ocean_physics_types, ONLY:	t_ho_params, v_params, &
   &                                WindAmplitude_at10m, SeaIceConcentration, &
   &                                WindMixingDecay, WindMixingLevel
  

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_physics'
  CHARACTER(LEN=12)           :: str_module    = 'ocePhysics  '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug

  ! Public interface

  !PUBLIC :: init_ho_physics
  PUBLIC :: init_ho_params
  PUBLIC :: update_ho_params
  PRIVATE :: calculate_leith_closure
  PUBLIC  :: calc_characteristic_physical_numbers

  ! variables
  TYPE (t_var_list), PUBLIC :: ocean_params_list

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Initialisation of ocean physics
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07)
!<Optimize:inUse:initOnly>
  SUBROUTINE init_ho_params(  patch_3d, p_phys_param, fu10 )
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE (t_ho_params)                          :: p_phys_param
    REAL(wp), TARGET                     :: fu10   (:,:) ! t_atmos_for_ocean%fu10
 
    ! Local variables
    INTEGER :: i, i_no_trac
    INTEGER :: je, blockNo, jk
    INTEGER :: start_index, end_index
    REAL(wp) :: z_lower_bound_diff, C_MPIOM
    REAL(wp) :: z_largest_edge_length ,z_diff_multfac, z_diff_efdt_ratio
    REAL(wp) :: points_in_munk_layer
    REAL(wp) :: minmaxmean_length(3), reference_scale
    REAL(wp) :: minDualEdgeLength, minEdgeLength, meanDualEdgeLength, meanEdgeLength
    REAL(wp) :: maxDualEdgeLength, maxEdgeLength
    REAL(wp) :: minCellArea, meanCellArea, maxCellArea
    TYPE(t_subset_range), POINTER :: all_edges, owned_edges
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp):: length_scale, dual_length_scale
    !-----------------------------------------------------------------------
    patch_2D   => patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    all_edges => patch_2D%edges%ALL
    owned_edges => patch_2D%edges%owned
    !-------------------------------------------------------------------------
    WindAmplitude_at10m => fu10
    !-------------------------------------------------------------------------
    points_in_munk_layer = REAL(n_points_in_munk_layer,wp)
    !Init from namelist
    p_phys_param%k_veloc_h_back = k_veloc_h
    p_phys_param%k_veloc_h      = k_veloc_h
    p_phys_param%a_veloc_v_back = k_veloc_v
    p_phys_param%a_veloc_v      = k_veloc_v

   IF(GMRedi_configuration==GMRedi_combined&
   &.OR.GMRedi_configuration==GM_only.OR.GMRedi_configuration==Redi_only)THEN

      p_phys_param%k_tracer_isoneutral = k_tracer_isoneutral_parameter
      p_phys_param%k_tracer_dianeutral = k_tracer_dianeutral_parameter

      p_phys_param%k_tracer_GM_kappa = k_tracer_GM_kappa_parameter
    ENDIF

    z_largest_edge_length = global_max(MAXVAL(patch_2D%edges%primal_edge_length))
    minmaxmean_length = global_minmaxmean(patch_2D%edges%dual_edge_length, owned_edges)
    minDualEdgeLength = minmaxmean_length(1)
    meanDualEdgeLength = minmaxmean_length(3)
    maxDualEdgeLength = minmaxmean_length(2)
    minmaxmean_length = global_minmaxmean(patch_2D%edges%primal_edge_length, owned_edges)
    minEdgeLength = minmaxmean_length(1)
    meanEdgeLength = minmaxmean_length(3)
    maxEdgeLength = minmaxmean_length(2)
    minmaxmean_length = global_minmaxmean(patch_2D%cells%area, patch_2D%cells%owned)
    minCellArea = minmaxmean_length(1)
    meanCellArea = minmaxmean_length(3)
    maxCellArea = minmaxmean_length(2)

    !Distinghuish between harmonic and biharmonic laplacian
    !Harmonic laplacian
    IF(veloc_diffusion_order==1)THEN
      p_phys_param%k_veloc_h_back = k_veloc_h
      p_phys_param%k_veloc_h      = k_veloc_h
      SELECT CASE(HorizontalViscosity_type)
      CASE(0)!no friction
        p_phys_param%k_veloc_h(:,:,:) = 0.0_wp

      CASE(1)!use uniform viscosity coefficient from namelist
        CALL calc_lower_bound_veloc_diff(  patch_2D, z_lower_bound_diff )
        IF(z_lower_bound_diff>p_phys_param%k_veloc_h_back)THEN
          ! SX9 cannot handle messages of that size -> split
          CALL message ('init_ho_params','WARNING: Specified diffusivity&
            & does not satisfy Munk criterion.')
          CALL message ('init_ho_params','WARNING: This may lead&
            & to stability problems for experiments with lateral boundaries')
        ENDIF

        p_phys_param%k_veloc_h(:,:,:) = p_phys_param%k_veloc_h_back
        !write(0,*)'lower bound of diffusivity:',z_lower_bound_diff
        WRITE(message_text,'(a,g25.16)') 'Lower bound of diffusivity:',z_lower_bound_diff
        CALL message ('init_ho_params', message_text)

      CASE(2)!calculate uniform viscosity coefficient, according to Munk criterion

        p_phys_param%k_veloc_h(:,:,:) = 3.82E-12_wp&
          & *(points_in_munk_layer*z_largest_edge_length)**3

      CASE(3)! calculate coefficients for each location based on MUNK layer
        DO blockNo = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, blockNo, start_index, end_index)
          DO je = start_index, end_index
            !calculate lower bound for diffusivity
            !The factor cos(lat) is omitted here, because of equatorial reference (cf. Griffies, eq. (18.29))
            p_phys_param%k_veloc_h(je,:,blockNo) = 3.82E-12_wp&
              & *(points_in_munk_layer*patch_2D%edges%primal_edge_length(je,blockNo))**3
          END DO
        END DO

      CASE(4)!calculate coefficients based on Leith closure. Start value is the same as case(2)  

        p_phys_param%k_veloc_h(:,:,:) = 3.82E-12_wp&
          & *(points_in_munk_layer*z_largest_edge_length)**3      

      END SELECT
      CALL dbg_print('horzVelocDiff:',p_phys_param%k_veloc_h ,str_module,0,in_subset=owned_edges)
      
      
      !Biharmonic laplacian
    ELSEIF(veloc_diffusion_order==2)THEN
    
      p_phys_param%k_veloc_h_back = HorizontalViscosityBackground_Biharmonic
!       p_phys_param%k_veloc_h      = HorizontalViscosityBackground_Biharmonic
      SELECT CASE(HorizontalViscosity_type)

      CASE(1)
        p_phys_param%k_veloc_h(:,:,:) = p_phys_param%k_veloc_h_back

      CASE(2)  
        !The number that controls all that the "z_diff_efdt_ratio"
        !is different. Higher z_diff_efdt_ratio decreases the final
        !diffusion coefficient
        z_diff_efdt_ratio = 10000.0_wp * HorizontalViscosityBackground_Biharmonic
        z_diff_multfac = (1._wp/ (z_diff_efdt_ratio*64._wp))/3._wp
        DO blockNo = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, blockNo, start_index, end_index)
          DO je = start_index, end_index
            p_phys_param%k_veloc_h(je,:,blockNo) = &
            & patch_2D%edges%area_edge(je,blockNo)*patch_2D%edges%area_edge(je,blockNo)*z_diff_multfac
          END DO
        END DO
        !          z_diff_multfac = 0.0045_wp*dtime/3600.0_wp
        !         DO blockNo = all_edges%start_block, all_edges%end_block
        !            CALL get_index_range(all_edges, blockNo, start_index, end_index)
        !            DO je = start_index, end_index
        !              p_phys_param%K_veloc_h(je,:,blockNo) = z_diff_multfac*&
        !              &maxval(patch_2D%edges%primal_edge_length)**4
        !            END DO
        !          END DO
      CASE(3)
        ! Simple scaling of the backgound diffusion by the dual edge length^4
        ! This is meant to be used with the non-uniform grids
        !This follows the MPI-OM convention
        C_MPIOM = biharmonic_const*dtime/3600.0_wp
        DO blockNo = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, blockNo, start_index, end_index)
          DO je = start_index, end_index

            length_scale = &
            & sqrt(patch_2D%edges%primal_edge_length(je,blockNo) * patch_2D%edges%dual_edge_length(je,blockNo))

            p_phys_param%k_veloc_h(je,:,blockNo)=C_MPIOM*length_scale**2
          END DO
        END DO

      CASE(4)!calculate coefficients based on Leith closure. Start value is the same as case(2).
      !This is just an initial value that will be overwritten as  soon as the velocities are different from zero.   

        p_phys_param%k_veloc_h(:,:,:) = 3.82E-12_wp&
          & *(points_in_munk_layer*z_largest_edge_length)**3
          
      CASE(5)
        ! Simple scaling of the constant diffusion
        DO blockNo = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, blockNo, start_index, end_index)
          p_phys_param%k_veloc_h(:,:,blockNo) = 0.0_wp
          DO je = start_index, end_index
            
            dual_length_scale = patch_2D%edges%dual_edge_length(je,blockNo) / maxDualEdgeLength
            length_scale = patch_2D%edges%primal_edge_length(je,blockNo) / maxEdgeLength
                        
!             length_scale = 0.5_wp * (length_scale**2 + length_scale**3)
!            length_scale = SQRT(length_scale * dual_length_scale) * dual_length_scale
!            length_scale = length_scale**2
!            length_scale = length_scale**2 * (1.0_wp + length_scale) * 0.5_wp
            length_scale = length_scale**3
            
            DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
              p_phys_param%k_veloc_h(je,jk,blockNo) = &
                & p_phys_param%k_veloc_h_back * &
                & (1.0_wp - HorizontalViscosity_ScaleWeight &
                &  +  HorizontalViscosity_ScaleWeight * length_scale)
            END DO
              
          END DO
        END DO
        
      CASE(6)
        ! Simple scaling of the constant diffusion
        reference_scale = 4.0_wp / (3.0_wp * maxEdgeLength * maxCellArea)
        DO blockNo = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, blockNo, start_index, end_index)
          p_phys_param%k_veloc_h(:,:,blockNo) = 0.0_wp
          DO je = start_index, end_index
            
            dual_length_scale = patch_2D%edges%dual_edge_length(je,blockNo) / maxDualEdgeLength
            length_scale = &
              & patch_2D%edges%primal_edge_length(je,blockNo)**2 * patch_2D%edges%dual_edge_length(je,blockNo) &
              & * reference_scale
                                   
            DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
              p_phys_param%k_veloc_h(je,jk,blockNo) = &
                & p_phys_param%k_veloc_h_back * &
                & (1.0_wp - HorizontalViscosity_ScaleWeight &
                &  +  HorizontalViscosity_ScaleWeight * length_scale)
            END DO
              
          END DO
        END DO
        
      CASE DEFAULT
         CALL finish ('mo_ocean_physics:init_ho_params',  &
          & 'option not supported')

      END SELECT
      CALL dbg_print('horzVelocDiff:',p_phys_param%k_veloc_h ,str_module,0,in_subset=owned_edges)

    ENDIF

    DO i=1, HorizontalViscosity_SmoothIterations
 !      CALL smooth_lapl_diff( patch_2D, patch_3d, p_phys_param%k_veloc_h, HorizontalViscosity_SpatialSmoothFactor )
    ENDDO


    DO i=1,no_tracer

      IF(i==1)THEN!temperature
        p_phys_param%k_tracer_h_back(i) = k_pot_temp_h
        p_phys_param%a_tracer_v_back(i) = k_pot_temp_v

      ELSEIF(i==2)THEN!salinity
        p_phys_param%k_tracer_h_back(2) = k_sal_h
        p_phys_param%a_tracer_v_back(2) = k_sal_v
      ELSE

        p_phys_param%k_tracer_h_back(i) = k_sal_h
        p_phys_param%a_tracer_v_back(i) = k_sal_v
       ! CALL finish ('mo_ocean_physics:init_ho_params',  &
       !   & 'number of tracers exceeds number of background values')
      ENDIF
      p_phys_param%k_tracer_h(:,:,:,i) = p_phys_param%k_tracer_h_back(i)
      p_phys_param%a_tracer_v(:,:,:,i) = p_phys_param%a_tracer_v_back(i)
    END DO

    p_phys_param%bottom_drag_coeff = bottom_drag_coeff

    DO i_no_trac=1, no_tracer
      CALL sync_patch_array(sync_c,patch_2D,p_phys_param%k_tracer_h(:,:,:,i_no_trac))
    END DO
    CALL sync_patch_array(sync_e,patch_2D,p_phys_param%k_veloc_h(:,:,:))
    
    ! precalculate exponential wind mixing decay with depth
    DO jk=2,n_zlev
      WindMixingDecay(jk) = EXP(-patch_3d%p_patch_1d(1)%del_zlev_m(jk-1)/WindMixingDecayDepth)
      WindMixingLevel(jk) = lambda_wind * patch_3d%p_patch_1d(1)%inv_del_zlev_m(jk-1)
    ENDDO 

  END SUBROUTINE init_ho_params
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Calculation of a lower bound for horizontal velocity diffusion of laplacian type ! that is
  !! required to have N (default =1) points in Munk layer. The lower bound is calculated ! with
  !! respect to the equator.
  !! The code is based on  Griffies, Fundamentals of ocean climate modeling, sect 18, p. 413.  !
  !! The lower bound is given in units [m^2/s].
  !!
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-08)
  !
  !
  SUBROUTINE calc_lower_bound_veloc_diff(  patch_2D, lower_bound_diff )
    TYPE(t_patch), TARGET, INTENT(in)  :: patch_2D
    REAL(wp), INTENT(inout)              :: lower_bound_diff

    ! Local variables
    REAL(wp) :: points_in_munk_layer
    REAL(wp)            :: z_largest_edge_length
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-------------------------------------------------------------------------
    edges_in_domain => patch_2D%edges%in_domain

    ! Get the largest edge length globally
    z_largest_edge_length = global_max(MAXVAL(patch_2D%edges%primal_edge_length))

    !calculate lower bound for diffusivity: The factor cos(lat) is omitted here, because of
    !equatorial reference (cf. Griffies, eq.  (18.29))
    points_in_munk_layer = REAL(n_points_in_munk_layer,wp)
    lower_bound_diff = 3.82E-12_wp*(points_in_munk_layer*z_largest_edge_length)**3

  END SUBROUTINE calc_lower_bound_veloc_diff
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-08)
!<Optimize:inUse:initOnly>
  SUBROUTINE smooth_lapl_diff( patch_2D,patch_3d, k_h, smoothFactor )
    TYPE(t_patch), TARGET, INTENT(in)  :: patch_2D
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(inout)    :: k_h(:,:,:)
    REAL(wp), INTENT(in)       ::smoothFactor
    ! Local variables
    INTEGER :: je,jv,blockNo,jk, jev, ile, ibe
    INTEGER :: il_v1,ib_v1, il_v2,ib_v2
    INTEGER :: start_index, end_index
    INTEGER :: i_startidx_v, i_endidx_v
    REAL(wp) :: z_k_ave_v(nproma,n_zlev,patch_2D%nblks_v)
    INTEGER  :: sea_edges_onLevel(n_zlev)
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER ::all_edges, verts_in_domain
    !-------------------------------------------------------------------------
    all_edges => patch_2D%edges%all
    verts_in_domain => patch_2D%verts%in_domain

    z_k_ave_v(:,:,:) = 0.0_wp

    DO blockNo = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, blockNo, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        DO jev = 1, patch_2D%verts%num_edges(jv,blockNo)
          ile = patch_2D%verts%edge_idx(jv,blockNo,jev)
          ibe = patch_2D%verts%edge_blk(jv,blockNo,jev)
          sea_edges_onLevel(:) = 0
          !             write(0,*) jv,blockNo, patch_2D%verts%num_edges(jv,blockNo), ":", ile, ibe
          DO jk = 1, patch_3D%p_patch_1D(1)%dolic_e(ile,ibe)
            z_k_ave_v(jv,jk,blockNo)= z_k_ave_v(jv,jk,blockNo) + k_h(ile,jk,ibe)
            sea_edges_onLevel(jk) = sea_edges_onLevel(jk) + 1
          END DO
            
          !IF(patch_2D%verts%num_edges(jv,blockNo)== 5)THEN
          !  z_K_ave_v(jv,jk,blockNo)=80000_wp!°z_K_max
          !ENDIF
        END DO ! jev = 1, patch_2D%verts%num_edges(jv,blockNo)
        DO jk = 1, patch_3D%p_patch_1D(1)%dolic_e(ile,ibe)
          IF(sea_edges_onLevel(jk) /= 0) & !.and.i_edge_ctr== patch_2D%verts%num_edges(jv,blockNo))THEN
            & z_k_ave_v(jv,jk,blockNo) = z_k_ave_v(jv,jk,blockNo) / REAL(sea_edges_onLevel(jk),wp)
        ENDDO
      ENDDO
    END DO

    ! we do need to sync here
    CALL sync_patch_array(sync_v, patch_2D, z_k_ave_v)    
    
    
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_index, end_index)
      DO je = start_index, end_index

        il_v1 = patch_2D%edges%vertex_idx(je,blockNo,1)
        ib_v1 = patch_2D%edges%vertex_blk(je,blockNo,1)
        il_v2 = patch_2D%edges%vertex_idx(je,blockNo,2)
        ib_v2 = patch_2D%edges%vertex_blk(je,blockNo,2)

        DO jk = 1, patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)
          
          k_h(je,jk,blockNo)= 0.5_wp * smoothFactor * (z_k_ave_v(il_v1,jk,ib_v1) + z_k_ave_v(il_v2,jk,ib_v2)) + &
            & (1.0_wp - smoothFactor) * k_h(je,jk,blockNo)
        END DO
      ENDDO
    END DO

    ! we do not need to sync edge coefficients
    ! CALL sync_patch_array(sync_e, patch_2D, k_h)   
    
    !---------Debug Diagnostics-------------------------------------------
    idt_src=0  ! output print levels - 0: print in any case
    CALL dbg_print('smoothed Laplac Diff.'     ,k_h                     ,str_module,idt_src, &
      & in_subset=patch_2D%edges%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE smooth_lapl_diff
  !-------------------------------------------------------------------------


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
!<Optimize:inUse:done>
  SUBROUTINE update_ho_params(patch_3d, ocean_state, fu10, concsum, params_oce,op_coeffs) !, calculate_density_func)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    REAL(wp), TARGET                     :: fu10   (:,:) ! t_atmos_for_ocean%fu10
    REAL(wp), TARGET                     :: concsum(:,:) ! t_sea_ice%concsum
    TYPE(t_ho_params), INTENT(inout)     :: params_oce
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs

    INTEGER :: tracer_index
    !-------------------------------------------------------------------------
    IF (ltimer) CALL timer_start(timer_upd_phys)
!     WindAmplitude_at10m => fu10
    SeaIceConcentration => concsum


    CALL calc_characteristic_physical_numbers(patch_3d, ocean_state)

    
    SELECT CASE (physics_parameters_type)
    CASE (physics_parameters_Constant_type)
      !nothing to do!In sbr init_ho_params (see above)
      !tracer mixing coefficient params_oce%A_tracer_v(:,:,:, tracer_index) is already
      !initialzed with params_oce%A_tracer_v_back(tracer_index)
      !and velocity diffusion coefficient

      ! prepare independent logicals for PP and convection parametrizations - not yet activated
      ! IF (.NOT. (l_convection .AND. l_pp_scheme)) THEN
      IF (ltimer) CALL timer_stop(timer_upd_phys)
      RETURN

    CASE (physics_parameters_ICON_PP_type)
      CALL update_physics_parameters_ICON_PP_scheme(patch_3d, ocean_state, params_oce)

    CASE (physics_parameters_MPIOM_PP_type)
      CALL update_physics_parameters_MPIOM_PP_scheme(patch_3d, ocean_state, fu10, concsum, params_oce)

    CASE (physics_parameters_ICON_PP_Edge_type)
      CALL update_physics_parameters_ICON_PP_Edge_scheme(patch_3d, ocean_state, params_oce)
      
    CASE (physics_parameters_ICON_PP_Edge_vnPredict_type)
      CALL update_physics_parameters_ICON_PP_Tracer(patch_3d, ocean_state)
!       CALL update_physics_parameters_ICON_PP_Edge_scheme(patch_3d, ocean_state, params_oce)
      ! the velovity friction will be updated during dynamics
      
    CASE default
      CALL finish("update_ho_params", "unknown physics_parameters_type")
    END SELECT

    IF(HorizontalViscosity_type==4)THEN

     IF(veloc_diffusion_order==1)THEN!.AND.veloc_diffusion_form==1)THEN
        CALL calculate_leith_closure(patch_3d, ocean_state, params_oce, op_coeffs)
      ENDIF

      IF(veloc_diffusion_order==2)THEN!.AND.veloc_diffusion_form==2)THEN
        CALL calculate_leith_closure(patch_3d, ocean_state, params_oce, op_coeffs)
      ENDIF

    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    !idt_src=4  ! output print levels (1-5, fix)
    !CALL dbg_print('UpdPar: p_vn%x(1)'         ,ocean_state%p_diag%p_vn%x(1)    ,str_module,idt_src, &
    !  & in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    !CALL dbg_print('UpdPar: p_vn%x(2)'         ,ocean_state%p_diag%p_vn%x(2)    ,str_module,idt_src, &
    !  & in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    idt_src=2  ! output print levels (1-5, fix)
    DO tracer_index = 1, no_tracer
      CALL dbg_print('UpdPar FinalTracerMixing'  ,params_oce%a_tracer_v(:,:,:,tracer_index), str_module,idt_src, &
        & in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    ENDDO
    CALL dbg_print('UpdPar FinalVelocMixing'   ,params_oce%a_veloc_v     ,str_module,idt_src, &
      & in_subset=patch_3d%p_patch_2d(1)%edges%owned)
    !---------------------------------------------------------------------
    IF (ltimer) CALL timer_stop(timer_upd_phys)
  END SUBROUTINE update_ho_params
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the viscosity coefficients following the Leith closure.
  !! Implemented is the pure leith closure and a modified version of it. Both options
  !! are available for harmonic as well as for biharmonic diffusion.
  !!
  !! Liteature: Fox-Kemper, Menemenlis, Can Large Eddy Simulations improve Mesoscale Rich Ocean Models ?
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  SUBROUTINE calculate_leith_closure(patch_3d, ocean_state, param, op_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(in)    :: op_coeff

    !Local variables
    REAL(wp) :: grad_vort_abs(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)

    INTEGER :: jk, blockNo, je, jc
    INTEGER :: start_cell_index, end_cell_index, cell_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level, level,end_level 
    INTEGER :: vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk, vertex3_idx, vertex3_blk
    INTEGER :: cell1_idx, cell1_blk, cell2_idx, cell2_blk    
    INTEGER :: LEITH_EXPONENT
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp):: div_e(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp):: div_c(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp):: vort_c(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp):: vort_e(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)    
    REAL(wp):: laplacian_vort(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp):: laplacian_div(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp):: length_scale(nproma, patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp)::  laplacian_vort_e, laplacian_div_e
    !REAL(wp):: size_grad_S_horz_vec(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1) 
    all_cells       => patch_2D%cells%all
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain

    start_level = 1

    grad_vort_abs(1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e)=0.0_wp
    length_scale (1:nproma, 1:patch_3D%p_patch_2d(1)%nblks_e)=0.0_wp
    div_e        (1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e)=0.0_wp
    div_c        (1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    vort_e        (1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e)=0.0_wp    
    vort_c       (1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    laplacian_vort(1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    laplacian_div(1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
   !-------------------------------------------------------------------------------   


!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, je) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      DO je=start_edge_index, end_edge_index 

!         end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

!         DO level=start_level,end_level

          length_scale(je,blockNo) = &
            & sqrt(patch_2D%edges%primal_edge_length(je,blockNo) * patch_2D%edges%dual_edge_length(je,blockNo))

!         END DO    
      END DO
    END DO ! blocks
!ICON_OMP_END_PARALLEL_DO

   ! Note: this sync is not needed
   CALL sync_patch_array(sync_e, patch_2D, length_scale)    


   !leith closure for Laplacian(harmonic) viscosity
   IF(veloc_diffusion_order==1)THEN

     LEITH_EXPONENT=3


    !1) calculation of horizontal gradient of vertical vorticity
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, je, level, vertex1_idx, &
!ICON_OMP vertex1_blk, vertex2_idx, vertex2_blk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      DO je=start_edge_index, end_edge_index 

!         end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

        vertex1_idx = patch_2d%edges%vertex_idx(je,blockNo,1)
        vertex1_blk = patch_2d%edges%vertex_blk(je,blockNo,1)
        vertex2_idx = patch_2d%edges%vertex_idx(je,blockNo,2)
        vertex2_blk = patch_2d%edges%vertex_blk(je,blockNo,2)

        DO level=start_level, patch_3D%p_patch_1D(1)%dolic_e(je,blockNo) 

          grad_vort_abs(je,level,blockNo)= &
          &sqrt((ocean_state%p_diag%vort(vertex1_idx,level,vertex1_blk)-ocean_state%p_diag%vort(vertex2_idx,level,vertex2_blk)&
          &/patch_2D%edges%primal_edge_length(je,blockNo))**2)

        END DO    
      END DO
    END DO ! blocks
!ICON_OMP_END_PARALLEL_DO

   ! Note: is sync is not needed
   CALL sync_patch_array(sync_e, patch_2D, grad_vort_abs)   


    !1) calculation leith closure or modified leith_closure
    IF(leith_closure==1)THEN
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
        DO je=start_edge_index, end_edge_index   
          end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)
          DO level=start_level,end_level     
            param%k_veloc_h(je,level,blockNo) = &
            &grad_vort_abs(je,level,blockNo) *(leith_closure_gamma* length_scale(je,blockNo))**LEITH_EXPONENT        
          END DO    
        END DO
      END DO ! blocks
    ELSEIF(leith_closure==2)THEN

      CALL div_oce_3d( ocean_state%p_diag%ptp_vn, patch_3D, op_coeff%div_coeff, div_c)

      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
        DO je=start_edge_index, end_edge_index   
          end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

          cell1_idx = patch_2d%edges%cell_idx(je,blockNo,1)
          cell1_blk = patch_2d%edges%cell_blk(je,blockNo,1)
          cell2_idx = patch_2d%edges%cell_idx(je,blockNo,2)
          cell2_blk = patch_2d%edges%cell_blk(je,blockNo,2)

          DO level=start_level,end_level

            div_e(je,level,blockNo)&
            &=0.5_wp*(div_c(cell1_idx,level,cell1_blk) + div_c(cell2_idx,level,cell2_blk))

            param%k_veloc_h(je,level,blockNo) &
            &=(leith_closure_gamma*length_scale(je,blockNo))**LEITH_EXPONENT&
            & *sqrt(grad_vort_abs(je,level,blockNo)**2  + div_e(je,level,blockNo)**2)          
          END DO    
        END DO
      END DO ! blocks

    ENDIF

    !leith closure for Bi-Laplacian(biharmonic) viscosity    
   ELSEIF(veloc_diffusion_order==2)THEN     

     LEITH_EXPONENT=6

    !1) calculation of vertical vorticity at cell centers

!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, jc, level, vertex1_idx, &
!ICON_OMP vertex1_blk, vertex2_idx, vertex2_blk, vertex3_idx, vertex3_blk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      DO jc=start_cell_index, end_cell_index 


        vertex1_idx = patch_2d%cells%vertex_idx(jc,blockNo,1)
        vertex1_blk = patch_2d%cells%vertex_blk(jc,blockNo,1)
        vertex2_idx = patch_2d%cells%vertex_idx(jc,blockNo,2)
        vertex2_blk = patch_2d%cells%vertex_blk(jc,blockNo,2)
        vertex3_idx = patch_2d%cells%vertex_idx(jc,blockNo,3)
        vertex3_blk = patch_2d%cells%vertex_blk(jc,blockNo,3)

        DO level=start_level, patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)

          vort_c(jc,level,blockNo) = &
            &  ocean_state%p_diag%vort(vertex1_idx,level,vertex1_blk) &
            &+ ocean_state%p_diag%vort(vertex2_idx,level,vertex2_blk) &
            &+ ocean_state%p_diag%vort(vertex3_idx,level,vertex3_blk)

        END DO    
      END DO
    END DO ! blocks
!ICON_OMP_END_PARALLEL_DO

   CALL sync_patch_array(sync_c, patch_2D, vort_c)

   !2) calculate laplacian of vertical velocity 
   CALL tracer_diffusion_horz_local(patch_3D, vort_c, ocean_state, vort_e)!,subset_range = edges_in_domain)   
   CALL div_oce_3d( vort_e, patch_3D, op_coeff%div_coeff, laplacian_vort)


   !3a) In case of pure Leith, we have all to calculate the coefficient       
   IF(leith_closure==1)THEN   
     !Now aggregate the final parameter
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, je, level, cell1_idx, &
!ICON_OMP cell1_blk, cell2_idx, cell2_blk, laplacian_vort_e) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
        DO je=start_edge_index, end_edge_index   
!           end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

          cell1_idx = patch_2d%edges%cell_idx(je,blockNo,1)
          cell1_blk = patch_2d%edges%cell_blk(je,blockNo,1)
          cell2_idx = patch_2d%edges%cell_idx(je,blockNo,2)
          cell2_blk = patch_2d%edges%cell_blk(je,blockNo,2)
          DO level=start_level, patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

            laplacian_vort_e&
            &=0.5_wp*(laplacian_vort(cell1_idx,level,cell1_blk) + laplacian_vort(cell2_idx,level,cell2_blk))
!             laplacian_vort_e = laplacian_vort_e*laplacian_vort_e

            !The leith viscosity coefficient
            param%k_veloc_h(je,level,blockNo) = &
              & ((leith_closure_gamma * length_scale(je,blockNo))**LEITH_EXPONENT) &
              & * sqrt(laplacian_vort_e * laplacian_vort_e)
          END DO    
        END DO
      END DO ! blocks
!ICON_OMP_END_PARALLEL_DO

!    DO level= start_level, 5!end_level-1  
!     write(*,*)'SIZES Leith-Viscosity',level,maxval(param%k_veloc_h(:,level,:)),&
!     & minval(param%k_veloc_h(:,level,:)), maxval(laplacian_vort(:,level,:)),minval(laplacian_vort(:,level,:)),&
!      &maxval(length_scale)**LEITH_EXPONENT,minval(length_scale)**LEITH_EXPONENT
!     END DO    

   !3b) In case of modified Leith, we need additionally the Laplacian of divergence    
   ELSEIF(leith_closure==2)THEN

     CALL div_oce_3d(ocean_state%p_diag%vn_time_weighted, patch_3D, op_coeff%div_coeff, div_c)
     !CALL div_oce_3d( ocean_state%p_diag%ptp_vn, patch_3D, op_coeff%div_coeff, div_c)

     !The next two calls calculate the laplacian of the divergence without any diffusion parameter    
     CALL tracer_diffusion_horz_local(patch_3D, div_c, ocean_state, div_e)!,subset_range = edges_in_domain)
     CALL div_oce_3d( div_e, patch_3D, op_coeff%div_coeff, laplacian_div)

     !Now aggregate the final parameter
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, je, level, cell1_idx, &
!ICON_OMP cell1_blk, cell2_idx, cell2_blk, laplacian_div_e, laplacian_vort_e) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
        DO je=start_edge_index, end_edge_index   
!           end_level = patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)

          cell1_idx = patch_2d%edges%cell_idx(je,blockNo,1)
          cell1_blk = patch_2d%edges%cell_blk(je,blockNo,1)
          cell2_idx = patch_2d%edges%cell_idx(je,blockNo,2)
          cell2_blk = patch_2d%edges%cell_blk(je,blockNo,2)

          DO level=start_level,  patch_3D%p_patch_1D(1)%dolic_e(je,blockNo) 

            laplacian_div_e = &
              & 0.5_wp*(laplacian_div(cell1_idx,level,cell1_blk) + laplacian_div(cell2_idx,level,cell2_blk))
            laplacian_div_e = laplacian_div_e * laplacian_div_e

            laplacian_vort_e = &
              & 0.5_wp*(laplacian_vort(cell1_idx,level,cell1_blk) + laplacian_vort(cell2_idx,level,cell2_blk))
            laplacian_vort_e = laplacian_vort_e * laplacian_vort_e

            !The modified leith viscosity coefficient            
            param%k_veloc_h(je,level,blockNo) = &
              & ((leith_closure_gamma*length_scale(je,blockNo))**LEITH_EXPONENT) &
              & * sqrt(laplacian_vort_e+laplacian_div_e)          
          END DO    
        END DO
      END DO ! blocks
!ICON_OMP_END_PARALLEL_DO

!    DO level= start_level, 5!end_level-1  
!     write(*,*)'SIZES Leith-Viscosity',level,maxval(param%k_veloc_h(:,level,:)),&
!     & minval(param%k_veloc_h(:,level,:)), &
!     &maxval(laplacian_vort(:,level,:)),minval(laplacian_vort(:,level,:)),&
!     &maxval(laplacian_div(:,level,:)),minval(laplacian_div(:,level,:)),&    
!     &maxval(length_scale)**LEITH_EXPONENT,minval(length_scale)**LEITH_EXPONENT
!     END DO    

   ENDIF     
 ENDIF  

   !CALL sync_patch_array(sync_e, patch_2D, grad_vort_abs)
   !Note: this sync most probably is nor needed
   CALL sync_patch_array(sync_e, patch_2D, param%k_veloc_h)

!idt_src=1  ! output print level (1-5, fix)
!      CALL debug_print_MaxMinMean('after ocean_gmres: h-new', minmaxmean, str_module, idt_src)

!    !---------DEBUG DIAGNOSTICS-------------------------------------------
     idt_src=2  ! output print level (1-5, fix)
     CALL dbg_print('calculate_leith_closure: viscosity',param%k_veloc_h,&
    &str_module,idt_src, in_subset=edges_in_domain)      
!   !---------------------------------------------------------------------   

!    DO level= start_level, 5!end_level-1  
!    write(*,*)'Leith-Viscosity',level,maxval(param%k_veloc_h(:,level,:)),&
!    & minval(param%k_veloc_h(:,level,:))!, maxval(grad_vort_abs(:,level,:)),minval(grad_vort_abs(:,level,:))!,&
!    !&maxval(grad_vort_abs(:,level,:)),minval(grad_vort_abs(:,level,:))
!    END DO

  END SUBROUTINE calculate_leith_closure
  !-------------------------------------------------------------------------

 !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates importan physical numbers: Richardson number, Buoancy frequency, baroclinic wave speed, Rossby radius.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2016).
  !!
!<Optimize:inUse>
  SUBROUTINE calc_characteristic_physical_numbers(patch_3d, ocean_state)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    !TYPE(t_ho_params),                 INTENT(inout) :: param
    !TYPE(t_operator_coeff),            INTENT(in)    :: op_coeff
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index,level,start_level,end_level, blockNo
    !INTEGER :: start_edge_index, end_edge_index, je   
    TYPE(t_subset_range), POINTER :: cells_in_domain, all_cells
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: lambda
    REAL(wp) :: depth_scale, depth
    REAL(wp) :: Coriolis_abs,z_grav_rho
    REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)
    REAL(wp), PARAMETER :: beta_param=2.29E-11
    REAL(wp), PARAMETER :: latitude_threshold=5.0_wp
    REAL(wp) :: Rossby_radius_equator(nproma)
    REAL(wp) :: Rossby_radius_offequator(nproma)  
    
    
    REAL(wp) :: z_shear_cell  
    REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev), density(n_zlev)
    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
    
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    all_cells       => patch_2D%cells%ALL
    start_level=1
    z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho
    !-------------------------------------------------------------------------------
    z_grav_rho = grav/OceanReferenceDensity

!ICON_OMP_PARALLEL PRIVATE(salinity)
    salinity(1:n_zlev) = sal_ref
!ICON_OMP_DO PRIVATE(start_cell_index, end_cell_index, cell_index, end_level, level, pressure, z_rho_up, z_rho_down, &
!ICON_OMP z_shear_cell) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
      ocean_state%p_diag%Richardson_Number(:,:, blockNo) = 0.0_wp
      z_vert_density_grad_c(:,:, blockNo) = 0.0_wp
      DO cell_index = start_cell_index, end_cell_index

        end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
        IF (end_level < 2) CYCLE

        IF(no_tracer >= 2) THEN
            salinity(1:end_level) = ocean_state%p_prog(nold(1))%tracer(cell_index,1:end_level,blockNo,2)
        ENDIF

        !--------------------------------------------------------
        pressure(2:end_level) = patch_3d%p_patch_1d(1)%depth_CellInterface(cell_index, 2:end_level, blockNo)&
        & * OceanReferenceDensity * sitodbar
        
        z_rho_up(1:end_level-1)  &
        &= calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(cell_index,1:end_level-1,blockNo,1), &
                                    & salinity(1:end_level-1), pressure(2:end_level),end_level-1)
          
        z_rho_down(2:end_level) &
        &= calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(cell_index,2:end_level,blockNo,1), &
                                    & salinity(2:end_level), pressure(2:end_level), end_level-1)

        DO level = 2, end_level
          z_shear_cell = dbl_eps + &
            & SUM((ocean_state%p_diag%p_vn(cell_index,level-1,blockNo)%x - ocean_state%p_diag%p_vn(cell_index,level,blockNo)%x)**2)
            
          z_vert_density_grad_c(cell_index,level,blockNo) = (z_rho_down(level) - z_rho_up(level-1)) *  &
            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(cell_index,level,blockNo)
            
          !adjusted vertical derivative (follows MOM, see Griffies-book, (p. 332, eq. (15.15)) or MOM-5 manual (sect. 23.7.1.1)  
          z_vert_density_grad_c(cell_index,level,blockNo)=min(z_vert_density_grad_c(cell_index,level,blockNo),-dbl_eps)  
            
          ocean_state%p_diag%Richardson_Number(cell_index, level, blockNo) &
          &= MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(cell_index,level,blockNo) * z_grav_rho * &
          & (z_rho_down(level) - z_rho_up(level-1)) / z_shear_cell, 0.0_wp) 
        END DO ! levels

      END DO ! index
    END DO
!ICON_OMP_END_DO
    

!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level, &
!ICON_OMP Rossby_radius_equator, Rossby_radius_offequator) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      
      Rossby_radius_equator   (:)               =0.0_wp
      Rossby_radius_offequator(:)               =0.0_wp
      ocean_state%p_diag%Wavespeed_baroclinic(:,blockNo)=0.0_wp     
       
      DO cell_index = start_cell_index, end_cell_index
            
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

          DO level = start_level, end_level
          
            !calculate buoyancy frequency
            ocean_state%p_diag%Buoyancy_Freq(cell_index,level,blockNo)&
            & = -z_grav_rho*z_vert_density_grad_c(cell_index,level,blockNo)
           
            !calculate baroclinic wave speed
            ocean_state%p_diag%Wavespeed_baroclinic(cell_index,blockNo)&
            &=ocean_state%p_diag%Wavespeed_baroclinic(cell_index,blockNo)&
            &+(1.0_wp/pi)*ocean_state%p_diag%Buoyancy_Freq(cell_index,level,blockNo)&
            &*patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(cell_index,level,blockNo)

          END DO
          
          IF(abs(patch_2d%cells%center(cell_index, blockNo)%lat)*rad2deg>=latitude_threshold)THEN
            Rossby_radius_offequator(cell_index)&
            &=ocean_state%p_diag%Wavespeed_baroclinic(cell_index,blockNo)&
            &/abs(patch_2d%cells%f_c(cell_index,blockNo))
          ENDIF
          Rossby_radius_equator(cell_index)&
            &=ocean_state%p_diag%Wavespeed_baroclinic(cell_index,blockNo)/(2.0_wp*beta_param)  
          
          ocean_state%p_diag%Rossby_Radius(cell_index,blockNo)&
            &=min(Rossby_radius_offequator(cell_index),Rossby_radius_equator(cell_index))
      END DO
    END DO
!ICON_OMP_END_DO


!ICON_OMP_END_PARALLEL

!   Do level=1,n_zlev
!   write(0,*)'max-min',level,&
!    &maxval( ocean_state%p_diag%Richardson_Number(:,level,:)),&
!    &minval( ocean_state%p_diag%Richardson_Number(:,level,:)),& 
!    &maxval( ocean_state%p_diag%Buoyancy_Freq(:,level,:)),&
!    &minval( ocean_state%p_diag%Buoyancy_Freq(:,level,:))
!   End do   
!   write(0,*)'max-min',&     
!    &maxval( ocean_state%p_diag%Wavespeed_baroclinic),&
!    &minval( ocean_state%p_diag%Wavespeed_baroclinic),& 
!    &maxval( ocean_state%p_diag%Rossby_Radius),&
!    &maxval(Rossby_radius_offequator), maxval(Rossby_radius_equator)
  
!  stop
   
  END SUBROUTINE calc_characteristic_physical_numbers
  !-------------------------------------------------------------------------


!  !-------------------------------------------------------------------------
!  !>
!  !! Update of ocean physics: This routine is used used only if time-dependent
!  !! changes of physical parametrizations.
!  !! Currently vertical mixing coefficients for tracers and vertical diffusivity are updated.
!  !! Dependent on the local Richardson number the diffusivity are calculated
!  !! (Large & Gent JPO 29, (1999), 449-464).
!  !! The formulation follows the MPI-OM implementation as described in Marsland et al. (Ocean
!  !! Modelling 5, 2003).
!  !! The notational convention is also taken from this paper( cf. eqs (14) and (19)).
!  !! What is missing is the fractional ice cover (see eqs. (15-16)).
!  !! Eq. (18) is the Redi part that is not implemented, yet
!  !!
!  !! @par Revision History
!  !! Initial release by Peter Korn, MPI-M (2011-02)
!  !<Optimize:inUse:done>
!  SUBROUTINE update_physics_parameters_ICON_PP_scheme(patch_3d, ocean_state, params_oce) !, calculate_density_func)
!
!    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
!    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
!    TYPE(t_ho_params), INTENT(inout)            :: params_oce
!
!    ! Local variables
!    INTEGER :: jc, blockNo, je,jk, tracer_index
!    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
!    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
!    INTEGER :: start_index, end_index
!    INTEGER :: levels
!
!    REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev), density(n_zlev)
!    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
!    REAL(wp) :: z_shear_cell, z_av0
!    REAL(wp) :: z_ri_cell               (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!    REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)
!
!    !Below is a set of variables and parameters for tracer and velocity
!    REAL(wp), PARAMETER :: z_0               = 40.0_wp
!    REAL(wp), PARAMETER :: z_c1_t            = 5.0_wp    !  PP diffusivity tuning constant
!    REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp    !  PP viscosity tuning constant
!    REAL(wp), PARAMETER :: z_threshold       = 5.0E-8_wp
!    REAL(wp) :: diffusion_weight
!    REAL(wp) :: z_grav_rho, z_inv_OceanReferenceDensity
!    REAL(wp) :: density_differ_edge, mean_z_r
!    !-------------------------------------------------------------------------
!    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
!    TYPE(t_patch), POINTER :: patch_2D
!
!    !-------------------------------------------------------------------------
!    patch_2D         => patch_3d%p_patch_2d(1)
!    edges_in_domain => patch_2D%edges%in_domain
!    !cells_in_domain => patch_2D%cells%in_domain
!    all_cells       => patch_2D%cells%ALL
!    z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho
!    levels = n_zlev
!
!    !-------------------------------------------------------------------------
!    z_av0 = richardson_veloc
!    z_grav_rho                   = grav/OceanReferenceDensity
!    z_inv_OceanReferenceDensity                = 1.0_wp/OceanReferenceDensity
!    !-------------------------------------------------------------------------
!!     IF (ltimer) CALL timer_start(timer_extra10)
!
!!ICON_OMP_PARALLEL PRIVATE(salinity)
!    salinity(1:levels) = sal_ref
!!ICON_OMP_DO PRIVATE(start_index, end_index, jc, levels, jk, pressure, z_rho_up, z_rho_down, &
!!ICON_OMP z_shear_cell, tracer_index, diffusion_weight) ICON_OMP_DEFAULT_SCHEDULE
!    DO blockNo = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, blockNo, start_index, end_index)
!      z_ri_cell(:,:, blockNo) = 0.0_wp
!      z_vert_density_grad_c(:,:, blockNo) = 0.0_wp
!      DO jc = start_index, end_index
!
!        levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
!        IF (levels < 2) CYCLE
!
!        IF(no_tracer >= 2) THEN
!            salinity(1:levels) = ocean_state%p_prog(nold(1))%tracer(jc,1:levels,blockNo,2)
!        ENDIF
!
!        !--------------------------------------------------------
!        pressure(2:levels) = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, 2:levels, blockNo) * OceanReferenceDensity * sitodbar
!        z_rho_up(1:levels-1)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,1:levels-1,blockNo,1), &
!          & salinity(1:levels-1), pressure(2:levels), levels-1)
!        z_rho_down(2:levels)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,2:levels,blockNo,1), &
!          & salinity(2:levels), pressure(2:levels), levels-1)
!
!        DO jk = 2, levels
!          z_shear_cell = dbl_eps + &
!            & SUM((ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x - ocean_state%p_diag%p_vn(jc,jk,blockNo)%x)**2)
!          z_vert_density_grad_c(jc,jk,blockNo) = (z_rho_down(jk) - z_rho_up(jk-1)) *  &
!            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,blockNo)
!          z_ri_cell(jc, jk, blockNo) = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,blockNo) * z_grav_rho * &
!            & (z_rho_down(jk) - z_rho_up(jk-1)) / z_shear_cell, 0.0_wp) ! do not use z_vert_density_grad_c,
!                                                                     ! this is canceled out in this formula
!        END DO ! levels
!
!      END DO ! index
!
!      DO tracer_index = 1, no_tracer
!
!        params_oce%a_tracer_v(start_index:end_index, 2:n_zlev, blockNo, tracer_index) =   &
!          & MERGE(max_vert_diff_trac,                    & ! activate convection
!          & params_oce%a_tracer_v_back(tracer_index) +   & ! calculate the richardson diffusion
!          &   richardson_tracer / ((1.0_wp + z_c1_t *    &
!          &   z_ri_cell(start_index:end_index, 2:n_zlev, blockNo))**3), &
!          & z_vert_density_grad_c(start_index:end_index, 2:n_zlev,blockNo) < convection_InstabilityThreshold)
!
!        DO jc = start_index, end_index
!          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
!          DO jk = 2, levels
!
!            IF (z_vert_density_grad_c(jc,jk,blockNo) < RichardsonDiffusion_threshold .AND. &
!                z_vert_density_grad_c(jc,jk,blockNo) >= convection_InstabilityThreshold) THEN
!              ! interpolate between convection and richardson diffusion
!              diffusion_weight =  &
!                & (z_vert_density_grad_c(jc,jk,blockNo) - convection_InstabilityThreshold) / &
!                & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)
!              params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = &
!                & max_vert_diff_trac * (1.0_wp - diffusion_weight) +&
!                & diffusion_weight * params_oce%a_tracer_v(jc,jk,blockNo,tracer_index)
!            ENDIF
!
!          ENDDO ! levels
!        ENDDO !  block index
!      ENDDO ! tracer_index]
!
!    END DO ! blocks
!!ICON_OMP_END_DO
!
!! !ICON_OMP_END_PARALLEL
!!     IF (ltimer) CALL timer_stop(timer_extra10)
!!     IF (ltimer) CALL timer_start(timer_extra11)
!! !ICON_OMP_PARALLEL
!    !--------------------------------------------
!    ! Calculate params_oce%A_veloc_v:
!    ! use mean values between the two cells; change to min, max if required
!!ICON_OMP_DO PRIVATE(start_index, end_index, je, cell_1_idx, cell_1_block, cell_2_idx, cell_2_block, &
!!ICON_OMP jk,  density_differ_edge, mean_z_r, diffusion_weight ) ICON_OMP_DEFAULT_SCHEDULE
!    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
!      DO je = start_index, end_index
!
!        cell_1_idx = patch_2D%edges%cell_idx(je,blockNo,1)
!        cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
!        cell_2_idx = patch_2D%edges%cell_idx(je,blockNo,2)
!        cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)
!
!        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
!            ! TODO: the following expect equally sized cells
!            ! compute density gradient at edges
!!            density_differ_edge = 0.5_wp * &
!!              & (z_vert_density_grad_c(cell_1_idx,jk,cell_1_block) + z_vert_density_grad_c(cell_2_idx,jk,cell_2_block))
!
!            !! density gradient smaller then threshold ('semi-stable'): use background value
!            ! note if density_differ_edge == z_threshold should be considered
!!             IF     (density_differ_edge <  convection_InstabilityThreshold) THEN
!!               ! turn on convection
!!               params_oce%a_veloc_v(je,jk,blockNo) = max_vert_diff_veloc
!!             ELSE
!              ! richardson diffusion
!              mean_z_r = 0.5_wp * (z_ri_cell(cell_1_idx,jk,cell_1_block) + z_ri_cell(cell_2_idx,jk,cell_2_block))
!              params_oce%a_veloc_v(je,jk,blockNo) = &
!                & params_oce%a_veloc_v_back +  &
!                & z_av0 /                      &
!                & ((1.0_wp + z_c1_v * mean_z_r)**2)
!
!!               IF (density_differ_edge < RichardsonDiffusion_threshold) THEN
!!                 diffusion_weight =  &
!!                   & (density_differ_edge - convection_InstabilityThreshold) / &
!!                   & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)
!! 
!!                 params_oce%a_veloc_v(je,jk,blockNo) = &
!!                   & max_vert_diff_veloc * (1.0_wp - diffusion_weight) + &
!!                   & params_oce%a_veloc_v(je,jk,blockNo) * diffusion_weight
!! 
!!                ENDIF
!
!!             ENDIF
!
!        END DO ! jk = 2, levels
!      ENDDO ! je = start_index, end_index
!    ENDDO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!!ICON_OMP_END_DO NOWAIT
!!ICON_OMP_END_PARALLEL
!!     IF (ltimer) CALL timer_stop(timer_extra11)
!
!  END SUBROUTINE update_physics_parameters_ICON_PP_scheme
!  !-------------------------------------------------------------------------
!
!  !-------------------------------------------------------------------------
!  !>
!  !! As in the update_physics_parameters_ICON_PP_scheme, but
!  !! velocity gradients for the vertical viscocity are clculated on edges
!  !!
!  !! @par Revision History
!  !! Initial release by Leonidas Linardakis, MPI-M (2011-02)
!  !<Optimize:inUse:done>
!  SUBROUTINE update_physics_parameters_ICON_PP_Edge_scheme(patch_3d, ocean_state, params_oce) !, calculate_density_func)
!
!    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
!    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
!    TYPE(t_ho_params), INTENT(inout)            :: params_oce
!
!    ! Local variables
!    INTEGER :: jc, blockNo, je,jk, tracer_index
!    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
!    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
!    INTEGER :: start_index, end_index
!    INTEGER :: levels
!
!    REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev), density(n_zlev)
!    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
!    REAL(wp) :: z_shear_cell
!    REAL(wp) :: z_ri_cell               (nproma, n_zlev)
!    REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)
!
!    !Below is a set of variables and parameters for tracer and velocity
!    REAL(wp), PARAMETER :: z_0               = 40.0_wp
!    REAL(wp), PARAMETER :: z_c1_t            = 5.0_wp    !  PP diffusivity tuning constant
!    REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp    !  PP viscosity tuning constant
!    REAL(wp), PARAMETER :: z_threshold       = 5.0E-8_wp
!    REAL(wp) :: diffusion_weight, instabilitySign
!    REAL(wp) :: z_grav_rho, z_inv_OceanReferenceDensity
!    REAL(wp) :: density_differ_edge, dz, richardson_edge, z_shear_edge
!    REAL(wp) :: tracer_windMixing(nproma, n_zlev+1), velocity_windMixing(n_zlev+1), z_vert_density_grad_e(n_zlev+1)
!    !-------------------------------------------------------------------------
!    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
!    TYPE(t_patch), POINTER :: patch_2D
!
!    !-------------------------------------------------------------------------
!    patch_2D         => patch_3d%p_patch_2d(1)
!    edges_in_domain => patch_2D%edges%in_domain
!    !cells_in_domain => patch_2D%cells%in_domain
!    all_cells       => patch_2D%cells%ALL
!    z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho
!    levels = n_zlev
!
!    !-------------------------------------------------------------------------
!    z_grav_rho                   = grav/OceanReferenceDensity
!    z_inv_OceanReferenceDensity                = 1.0_wp/OceanReferenceDensity
!    !-------------------------------------------------------------------------
!!     IF (ltimer) CALL timer_start(timer_extra10)
!
!!ICON_OMP_PARALLEL PRIVATE(salinity)
!    salinity(1:levels) = sal_ref
!!ICON_OMP_DO PRIVATE(start_index, end_index, jc, levels, jk, pressure, z_rho_up, z_rho_down, &
!!ICON_OMP z_shear_cell, z_ri_cell, tracer_index, diffusion_weight, instabilitySign, tracer_windMixing) ICON_OMP_DEFAULT_SCHEDULE
!    DO blockNo = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, blockNo, start_index, end_index)
!      z_ri_cell(:,:) = 0.0_wp
!      z_vert_density_grad_c(:,:, blockNo) = 0.0_wp
!      DO jc = start_index, end_index
!
!        levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
!        IF (levels < 2) CYCLE
!
!        IF(no_tracer >= 2) THEN
!            salinity(1:levels) = ocean_state%p_prog(nold(1))%tracer(jc,1:levels,blockNo,2)
!        ENDIF
!
!        !--------------------------------------------------------
!        pressure(2:levels) = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, 2:levels, blockNo) * OceanReferenceDensity * sitodbar
!        z_rho_up(1:levels-1)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,1:levels-1,blockNo,1), &
!          & salinity(1:levels-1), pressure(2:levels), levels-1)
!        z_rho_down(2:levels)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,2:levels,blockNo,1), &
!          & salinity(2:levels), pressure(2:levels), levels-1)
!
!        DO jk = 2, levels
!          z_shear_cell = dbl_eps + &
!            & SUM((ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x - ocean_state%p_diag%p_vn(jc,jk,blockNo)%x)**2)
!          z_vert_density_grad_c(jc,jk,blockNo) = (z_rho_down(jk) - z_rho_up(jk-1)) *  &
!            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,blockNo)
!          z_ri_cell(jc, jk) = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,blockNo) * z_grav_rho * &
!            & (z_rho_down(jk) - z_rho_up(jk-1)) / z_shear_cell, 0.0_wp) ! do not use z_vert_density_grad_c,
!                                                                     ! this is canceled out in this formula
!        END DO ! levels
!
!      END DO ! index
!      !-----------------------------------------------------------
!      tracer_windMixing(:,:) = 0.0_wp
!      IF (use_wind_mixing) THEN
!        DO jc = start_index, end_index
!
!          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
!          ! wind-mixing: Marsland et al., 2003
!          ! wind-mixing at surface, eq. (15) of Marsland et al., 2003
!          tracer_windMixing(jc,1) = tracer_TopWindMixing  * WindAmplitude_at10m(jc,blockNo)**3 * &
!            & (1.0_wp - SeaIceConcentration(jc,blockNo))
!
!          ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
!          DO jk = 2, levels
!           tracer_windMixing(jc,jk) =  tracer_windMixing(jc,jk-1) * WindMixingDecay(jk) * &
!             & WindMixingLevel(jk) / (WindMixingLevel(jk) + MAX(z_vert_density_grad_c(jc,jk,blockNo),0.0_wp))
!
!          END DO! levels
!
!        END DO ! index
!
!      END IF  ! use_wind_mixing
!      !-----------------------------------------------------------
!
!      DO tracer_index = 1, no_tracer
!        IF (tracer_index == 1) THEN
!          instabilitySign = 0.0 ! always enable convection for temperature
!        ELSE
!          instabilitySign = Salinity_ConvectionRestrict
!        ENDIF
!        DO jc = start_index, end_index
!          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
!          DO jk = 2, levels
!            ! by default use the richardson formula, no convection
!            params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = MIN( &
!              & params_oce%a_tracer_v_back(tracer_index) +   &
!              & richardson_tracer / ((1.0_wp + z_c1_t * z_ri_cell(jc, jk))**3) + tracer_windMixing(jc,jk), &
!              & max_vert_diff_trac)
!          IF (instabilitySign * &
!                 (ocean_state%p_prog(nold(1))%tracer(jc,jk-1,blockNo,tracer_index) -         &
!                  ocean_state%p_prog(nold(1))%tracer(jc,jk,blockNo,tracer_index)) >= 0.0_wp) & ! the = is important, do not change!
!            THEN
!              ! possibly convection
!              IF (z_vert_density_grad_c(jc,jk,blockNo) <= convection_InstabilityThreshold) &
!              THEN
!                ! convection
!                params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = max_vert_diff_trac
!              ELSE
!                IF (z_vert_density_grad_c(jc,jk,blockNo) < RichardsonDiffusion_threshold) THEN
!                  ! interpolate between convection and richardson diffusion
!                  diffusion_weight =  &
!                    & (z_vert_density_grad_c(jc,jk,blockNo) - convection_InstabilityThreshold) / &
!                    & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)
!                  params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = &
!                    & max_vert_diff_trac * (1.0_wp - diffusion_weight) +&
!                    & diffusion_weight * params_oce%a_tracer_v(jc,jk,blockNo,tracer_index)
!                ENDIF
!              ENDIF
!            ENDIF ! possibly convection
!            
!          ENDDO ! levels
!        ENDDO !  block index
!      ENDDO ! tracer_index]
!
!    END DO ! blocks
!!ICON_OMP_END_DO
!
!! !ICON_OMP_END_PARALLEL
!!     IF (ltimer) CALL timer_stop(timer_extra10)
!!     IF (ltimer) CALL timer_start(timer_extra11)
!! !ICON_OMP_PARALLEL
!    !--------------------------------------------
!    ! Calculate params_oce%A_veloc_v:
!!ICON_OMP_DO PRIVATE(start_index, end_index, je, cell_1_idx, cell_1_block, cell_2_idx, cell_2_block, &
!!ICON_OMP jk, dz, density_differ_edge, z_shear_edge, richardson_edge,z_vert_density_grad_e,velocity_windMixing) ICON_OMP_DEFAULT_SCHEDULE
!    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
!      DO je = start_index, end_index
!
!        cell_1_idx = patch_2D%edges%cell_idx(je,blockNo,1)
!        cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
!        cell_2_idx = patch_2D%edges%cell_idx(je,blockNo,2)
!        cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)
!
!        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
!          z_vert_density_grad_e(jk) =  0.5_wp * &
!            & (z_vert_density_grad_c(cell_1_idx,jk,cell_1_block) +   &
!            &  z_vert_density_grad_c(cell_2_idx,jk,cell_2_block))
!        ENDDO
!
!        velocity_windMixing(:) = 0.0_wp
!        IF (use_wind_mixing) THEN
!          velocity_windMixing(1) =           &
!            & velocity_TopWindMixing &
!            & *  (0.5_wp * &
!            &    (WindAmplitude_at10m(cell_1_idx,cell_1_block) + WindAmplitude_at10m(cell_2_idx,cell_2_block)))**3 &
!            & * (1.0_wp - 0.5_wp *       &
!            &    (SeaIceConcentration(cell_1_idx,cell_1_block) + SeaIceConcentration(cell_2_idx,cell_2_block)))
!
!            ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
!            DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
!              velocity_windMixing(jk) =  velocity_windMixing(jk-1) * WindMixingDecay(jk) * &
!                & WindMixingLevel(jk) / (WindMixingLevel(jk) + MAX(z_vert_density_grad_e(jk),0.0_wp))
!            END DO! levels
!        END IF  ! use_wind_mixing
!
!        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
!          ! TODO: the following expect equally sized cells
!          ! compute density gradient at edges
!          dz = 0.5_wp * (patch_3d%p_patch_1d(1)%prism_thick_e(je,jk-1,blockNo) + &
!            &            patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,blockNo))
!          density_differ_edge = z_vert_density_grad_e(jk) * dz
!          z_shear_edge = dbl_eps + &
!            & (ocean_state%p_prog(nold(1))%vn(je,jk,  blockNo) - &
!            &  ocean_state%p_prog(nold(1))%vn(je,jk-1,blockNo)   )**2
!
!          richardson_edge = MAX(dz * z_grav_rho * density_differ_edge / z_shear_edge, 0.0_wp)
!          
!          params_oce%a_veloc_v(je,jk,blockNo) =                 &
!            & params_oce%a_veloc_v_back * dz +             &
!            & richardson_veloc /                           &
!            & ((1.0_wp + z_c1_v * richardson_edge)**2) +   &
!            & velocity_windMixing(jk)
!
!        END DO ! jk = 2, levels
!      ENDDO ! je = start_index, end_index
!    ENDDO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!!ICON_OMP_END_DO NOWAIT
!!ICON_OMP_END_PARALLEL
!!     IF (ltimer) CALL timer_stop(timer_extra11)
!
!  END SUBROUTINE update_physics_parameters_ICON_PP_Edge_scheme
!  !-------------------------------------------------------------------------
!
!
!  !-------------------------------------------------------------------------
!
!  !-------------------------------------------------------------------------
!  !>
!  !! As in the update_physics_parameters_ICON_PP_scheme, but
!  !! velocity gradients for the vertical viscocity are clculated on edges
!  !!
!  !! @par Revision History
!  !! Initial release by Leonidas Linardakis, MPI-M (2011-02)
!  !<Optimize:inUse:done>
!  SUBROUTINE update_physics_parameters_ICON_PP_Edge_vnPredict_scheme(patch_3d, &
!    & blockNo, start_index, end_index, ocean_state, vn_predict) !, calculate_density_func)
!
!    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
!    INTEGER, INTENT(in) :: blockNo, start_index, end_index
!    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
!    REAL(wp) :: vn_predict(:,:)
!
!    ! Local variables
!    INTEGER :: je,jk
!    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
!    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
!    INTEGER :: levels
!
!    !Below is a set of variables and parameters for tracer and velocity
!    REAL(wp), PARAMETER :: z_0               = 40.0_wp
!    REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp    !  PP viscosity tuning constant
!    REAL(wp) :: diffusion_weight
!    REAL(wp) :: z_grav_rho, z_inv_OceanReferenceDensity
!    REAL(wp) :: density_differ_edge, dz, richardson_edge, z_shear_edge, vn_diff, new_velocity_friction
!    !-------------------------------------------------------------------------
!    REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)
!    REAL(wp) :: wind_mixing(1:n_zlev+1)
!    REAL(wp) :: z_vert_density_grad_e(1:n_zlev+1)
!    TYPE(t_patch), POINTER :: patch_2D
!    TYPE(t_ho_params), POINTER :: params_oce
!
!    !-------------------------------------------------------------------------
!    params_oce      => v_params
!    patch_2D        => patch_3d%p_patch_2d(1)
!    z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho ! already calculated
!    levels = n_zlev
!
!    !-------------------------------------------------------------------------
!    z_grav_rho                   = grav/OceanReferenceDensity
!    z_inv_OceanReferenceDensity                = 1.0_wp/OceanReferenceDensity
!    !-------------------------------------------------------------------------
!
!    DO je = start_index, end_index
!
!      cell_1_idx = patch_2D%edges%cell_idx(je,blockNo,1)
!      cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
!      cell_2_idx = patch_2D%edges%cell_idx(je,blockNo,2)
!      cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)
!      
!      DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
!        z_vert_density_grad_e(jk) =  0.5_wp * &
!          & (z_vert_density_grad_c(cell_1_idx,jk,cell_1_block) +   &
!          &  z_vert_density_grad_c(cell_2_idx,jk,cell_2_block))
!      ENDDO
!      
!      wind_mixing(:) = 0.0_wp
!      IF (use_wind_mixing .AND. patch_3d%p_patch_1d(1)%dolic_e(je, blockNo) > 0) THEN
!!         IF (cell_1_idx < 1 .or. cell_1_block < 1 .or. &
!!           & cell_2_idx < 1 .or. cell_2_block < 1)     &
!!           & CALL finish("ICON_PP_Edge_vnPredict wind mixing", "invalid cell pointers")
!        ! wind-mixing: Marsland et al., 2003
!        ! wind-mixing at surface, eq. (15) of Marsland et al., 2003
!        wind_mixing(1) =           &
!          & velocity_TopWindMixing &
!          & *  (0.5_wp *           &
!          &    (WindAmplitude_at10m(cell_1_idx,cell_1_block) + WindAmplitude_at10m(cell_2_idx,cell_2_block)))**3 &
!          & * (1.0_wp - 0.5_wp *       &
!          &    (SeaIceConcentration(cell_1_idx,cell_1_block) + SeaIceConcentration(cell_2_idx,cell_2_block))) 
!
!        ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
!        DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
!          wind_mixing(jk) =  wind_mixing(jk-1) * WindMixingDecay(jk) * &
!            & WindMixingLevel(jk) / (WindMixingLevel(jk) + MAX(z_vert_density_grad_e(jk),0.0_wp)) 
!          ! write(0,*) jk, wind_mixing(jk)
!        END DO! levels
!
!      END IF  ! use_wind_mixing
!
!      DO jk = 2, patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
!        ! TODO: the following expect equally sized cells
!        ! compute density gradient at edges
!        dz = 0.5_wp * (patch_3d%p_patch_1d(1)%prism_thick_e(je,jk-1,blockNo) + &
!          &            patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,blockNo))
!        density_differ_edge = z_vert_density_grad_e(jk) * dz
!!         vn_diff = MAX( &
!!           & ABS(ocean_state%p_prog(nold(1))%vn(je,jk,  blockNo) - &
!!           &      ocean_state%p_prog(nold(1))%vn(je,jk-1,blockNo)), &
!!           & ABS(vn_predict(je,jk) - &
!!           &     vn_predict(je,jk-1)))
!        vn_diff =  &
!          & ABS(vn_predict(je,jk) - &
!          &     vn_predict(je,jk-1))
!
!        z_shear_edge = dbl_eps + vn_diff**2
!
!        richardson_edge = MAX(dz * z_grav_rho * density_differ_edge / z_shear_edge, 0.0_wp)
!        
!        new_velocity_friction = &
!          & params_oce%a_veloc_v_back * dz +                              &
!          & richardson_veloc / ((1.0_wp + z_c1_v * richardson_edge)**2)+  &
!          & wind_mixing(jk)
!
!        ! the average of the calculated velocity friction based on the old velocity and the predicted one
!        params_oce%a_veloc_v(je,jk,blockNo) = &
!!           & MAX(params_oce%a_veloc_v(je,jk,blockNo), new_velocity_friction )
!           & VerticalViscosity_TimeWeight * params_oce%a_veloc_v(je,jk,blockNo) + &
!           & (1.0_wp - VerticalViscosity_TimeWeight) * new_velocity_friction
!
!      END DO ! jk = 2, levels
!    ENDDO ! je = start_index, end_index
!
!  END SUBROUTINE update_physics_parameters_ICON_PP_Edge_vnPredict_scheme
!  !-------------------------------------------------------------------------
!
!  !-------------------------------------------------------------------------
!  !>
!  !! As in the update_physics_parameters_ICON_PP_scheme, but
!  !! velocity gradients for the vertical viscocity are calculated on edges
!  !!
!  !! @par Revision History
!  !! Initial release by Leonidas Linardakis, MPI-M (2011-02)
!  !<Optimize:inUse:done>
!  SUBROUTINE update_physics_parameters_ICON_PP_Tracer(patch_3d, ocean_state) !, calculate_density_func)
!
!    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
!    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
!
!    ! Local variables
!    INTEGER :: jc, blockNo, je,jk, tracer_index
!    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
!    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
!    INTEGER :: start_index, end_index
!    INTEGER :: levels
!
!    REAL(wp) :: z_rho_up(n_zlev), z_rho_down(n_zlev), density(n_zlev)
!    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
!    REAL(wp) :: z_shear_cell
!    REAL(wp) :: z_ri_cell(nproma, n_zlev+1)
!    REAL(wp) :: wind_mixing(nproma, n_zlev+1)
!    REAL(wp), POINTER :: z_vert_density_grad_c(:,:,:)
!
!    !Below is a set of variables and parameters for tracer and velocity
!    REAL(wp), PARAMETER :: z_0               = 40.0_wp
!    REAL(wp), PARAMETER :: z_c1_t            = 5.0_wp    !  PP diffusivity tuning constant
!    REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp    !  PP viscosity tuning constant
!    REAL(wp), PARAMETER :: z_threshold       = 5.0E-8_wp
!    REAL(wp) :: diffusion_weight
!    REAL(wp) :: z_grav_rho, z_inv_OceanReferenceDensity
!    REAL(wp) :: density_differ_edge, dz, richardson_edge, z_shear_edge
!    !-------------------------------------------------------------------------
!    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
!    TYPE(t_patch), POINTER :: patch_2D
!    TYPE(t_ho_params), POINTER  :: params_oce
!
!    !-------------------------------------------------------------------------
!    params_oce      => v_params
!    patch_2D        => patch_3d%p_patch_2d(1)
!    edges_in_domain => patch_2D%edges%in_domain
!    !cells_in_domain => patch_2D%cells%in_domain
!    all_cells       => patch_2D%cells%ALL
!    z_vert_density_grad_c => ocean_state%p_diag%zgrad_rho
!    levels = n_zlev
!
!    !-------------------------------------------------------------------------
!    z_grav_rho                   = grav/OceanReferenceDensity
!    z_inv_OceanReferenceDensity                = 1.0_wp/OceanReferenceDensity
!    !-------------------------------------------------------------------------
!!     IF (ltimer) CALL timer_start(timer_extra10)
!
!!ICON_OMP_PARALLEL PRIVATE(salinity)
!    salinity(1:levels) = sal_ref
!!ICON_OMP_DO PRIVATE(start_index, end_index, jc, levels, jk, pressure, z_rho_up, z_rho_down, &
!!ICON_OMP z_shear_cell, z_ri_cell, tracer_index, diffusion_weight, wind_mixing) ICON_OMP_DEFAULT_SCHEDULE
!    DO blockNo = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, blockNo, start_index, end_index)
!      z_ri_cell(:,:) = 0.0_wp
!      z_vert_density_grad_c(:,:, blockNo) = 0.0_wp
!      DO jc = start_index, end_index
!
!        levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
!        IF (levels < 2) CYCLE
!
!        IF(no_tracer >= 2) THEN
!            salinity(1:levels) = ocean_state%p_prog(nold(1))%tracer(jc,1:levels,blockNo,2)
!        ENDIF
!
!        !--------------------------------------------------------
!        pressure(2:levels) = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, 2:levels, blockNo) * OceanReferenceDensity * sitodbar
!        z_rho_up(1:levels-1)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,1:levels-1,blockNo,1), &
!          & salinity(1:levels-1), pressure(2:levels), levels-1)
!        z_rho_down(2:levels)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,2:levels,blockNo,1), &
!          & salinity(2:levels), pressure(2:levels), levels-1)
!
!        DO jk = 2, levels
!          z_shear_cell = dbl_eps + &
!            & SUM((ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x - ocean_state%p_diag%p_vn(jc,jk,blockNo)%x)**2)
!          z_vert_density_grad_c(jc,jk,blockNo) = (z_rho_down(jk) - z_rho_up(jk-1)) *  &
!            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,blockNo)
!          z_ri_cell(jc, jk) = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,blockNo) * z_grav_rho * &
!            & (z_rho_down(jk) - z_rho_up(jk-1)) / z_shear_cell, 0.0_wp) ! do not use z_vert_density_grad_c,
!                                                                        ! this is canceled out in this formula
!        END DO ! levels
!
!      END DO ! index
!      !-----------------------------------------------------------
!      wind_mixing(:,:) = 0.0_wp
!      IF (use_wind_mixing) THEN
!        DO jc = start_index, end_index
!
!          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
!          ! wind-mixing: Marsland et al., 2003
!          ! wind-mixing at surface, eq. (15) of Marsland et al., 2003
!          wind_mixing(jc,1) = tracer_TopWindMixing  * WindAmplitude_at10m(jc,blockNo)**3 * &
!            & (1.0_wp - SeaIceConcentration(jc,blockNo))
!
!          ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
!          DO jk = 2, levels
!           wind_mixing(jc,jk) =  wind_mixing(jc,jk-1) * WindMixingDecay(jk) * &
!             & WindMixingLevel(jk) / (WindMixingLevel(jk) + MAX(z_vert_density_grad_c(jc,jk,blockNo),0.0_wp))
!
!          END DO! levels
!          
!        END DO ! index
!
!      END IF  ! use_wind_mixing
!      !-----------------------------------------------------------
!
!      DO tracer_index = 1, no_tracer
!
!        params_oce%a_tracer_v(start_index:end_index, 2:n_zlev, blockNo, tracer_index) =   &
!          & MERGE(                                       &
!          & max_vert_diff_trac,                          & ! activate convection
!          & params_oce%a_tracer_v_back(tracer_index) +   & ! calculate the richardson diffusion
!          &   richardson_tracer / ((1.0_wp + z_c1_t *    &
!          &   z_ri_cell(start_index:end_index, 2:n_zlev))**3) + &
!          &   wind_mixing(start_index:end_index, 2:n_zlev), &
!          & z_vert_density_grad_c(start_index:end_index, 2:n_zlev,blockNo) <= convection_InstabilityThreshold)
!
!        DO jc = start_index, end_index
!          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
!          DO jk = 2, levels
!
!            IF (z_vert_density_grad_c(jc,jk,blockNo) < RichardsonDiffusion_threshold .AND. &
!                z_vert_density_grad_c(jc,jk,blockNo) > convection_InstabilityThreshold) THEN
!              ! interpolate between convection and richardson diffusion
!              diffusion_weight =  &
!                & (z_vert_density_grad_c(jc,jk,blockNo) - convection_InstabilityThreshold) / &
!                & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)
!              params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = &
!                & max_vert_diff_trac * (1.0_wp - diffusion_weight) +&
!                & diffusion_weight * params_oce%a_tracer_v(jc,jk,blockNo,tracer_index)
!            ENDIF
!
!          ENDDO ! levels
!        ENDDO !  block index
!      ENDDO ! tracer_index]
!
!    END DO ! blocks
!!ICON_OMP_END_DO NOWAIT
!!ICON_OMP_END_PARALLEL
!
!  END SUBROUTINE update_physics_parameters_ICON_PP_Tracer
!  !-------------------------------------------------------------------------
!
!  !-------------------------------------------------------------------------
!  !>
!  !! Update of ocean physics parameters
!  !!
!  !! Update of ocean physics: This routine is used used only if time-dependent
!  !! changes of physical parametrizations.
!  !! Currently vertical mixing coefficients for tracers and vertical diffusivity are updated.
!  !! Dependent on the local Richardson number the diffusivity are calculated
!  !! (Large & Gent JPO 29, (1999), 449-464).
!  !! The formulation follows the MPI-OM implementation as described in Marsland et al. (Ocean
!  !! Modelling 5, 2003).
!  !! The notational convention is also taken from this paper( cf. eqs (14) and (19)).
!  !! What is missing is the fractional ice cover (see eqs. (15-16)).
!  !! Eq. (18) is the Redi part that is not implemented, yet
!  !!
!  !! @par Revision History
!  !! Initial release by Peter Korn, MPI-M (2011-02)
!!<Optimize:inUse>
!  SUBROUTINE update_physics_parameters_MPIOM_PP_scheme(patch_3d, ocean_state, fu10, concsum, params_oce) !, calculate_density_func)
!
!    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
!    TYPE(t_hydro_ocean_state), TARGET     :: ocean_state
!    REAL(wp),          INTENT(in)         :: fu10   (nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) ! t_atmos_for_ocean%fu10
!    REAL(wp),          INTENT(in)         :: concsum(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) ! t_sea_ice%concsum
!    TYPE(t_ho_params), INTENT(inout)      :: params_oce
!!     INTERFACE !This contains the function version of the actual EOS as chosen in namelist
!!       FUNCTION calculate_density_func(tpot, sal, press) result(rho)
!!         USE mo_kind, ONLY: wp
!!         REAL(wp), INTENT(in) :: tpot
!!         REAL(wp), INTENT(in) :: sal
!!         REAL(wp), INTENT(in) :: press
!!         REAL(wp) :: rho
!!       ENDFUNCTION calculate_density_func
!!     END INTERFACE
!
!    ! Local variables
!    INTEGER :: jc, blockNo, je,jk, tracer_index
!    !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
!    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx,cell_2_block
!    INTEGER :: start_index, end_index
!    INTEGER :: levels, jk_max
!
!    REAL(wp) :: rho_up(n_zlev), rho_down(n_zlev)
!    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)
!    REAL(wp) :: vert_velocity_shear
!    REAL(wp), POINTER :: vert_density_grad(:,:,:)
!
!    ! Local 3dim variables
!    REAL(wp) :: richardson_no(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!    REAL(wp) :: dv_wind      (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!    REAL(wp) :: av_wind      (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
!
!    ! Below is a set of variables and parameters for diffusion of tracer and velocity
!    ! lambda_diff: relne in MPIOM (0.4), Lambda_D (0.6) in Marsland et al. (2003), same vor Lambda_v
!    REAL(wp), PARAMETER :: lambda_diff       = 0.4_wp    !  eddy diffusion relaxation constant
!    REAL(wp), PARAMETER :: lambda_visc       = 0.4_wp    !  eddy viscosity relaxation constant
!    REAL(wp), PARAMETER :: z0_wind           = 40.0_wp   !  exponential decay of wind mixing with depth
!    REAL(wp), PARAMETER :: v10m_ref          = 6.0_wp    !  wind mixing 10m reference windspeed
!    REAL(wp), PARAMETER :: crd               = 5.0_wp    !  PP diffusivity tuning constant
!    REAL(wp), PARAMETER :: crv               = 5.0_wp    !  PP viscosity tuning constant
!
!    REAL(wp) :: decay_wind_depth, wind_param, vdensgrad_inter, densgrad_k, densgrad_kp1
!    REAL(wp) :: v10mexp_3, wma_pv, wma_pd
!    REAL(wp) :: diffusion_weight, loc_eps
!    REAL(wp) :: dv_old, dv_back, dv_rich
!    REAL(wp) :: av_old, av_back, av_rich
!    REAL(wp) :: onem1_lambda_d, onem1_lambda_v
!    REAL(wp) :: grav_rho
!    REAL(wp) :: mean_density_differ_edge, mean_richardson_e, fu10_e, conc_e
!    REAL(wp) :: vdgfac_bot(n_zlev)
!    !-------------------------------------------------------------------------
!    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
!    TYPE(t_patch), POINTER :: p_patch
!
!    !-------------------------------------------------------------------------
!    p_patch         => patch_3d%p_patch_2d(1)
!    edges_in_domain => p_patch%edges%in_domain
!    !cells_in_domain => p_patch%cells%in_domain
!    all_cells       => p_patch%cells%ALL
!    vert_density_grad => ocean_state%p_diag%zgrad_rho
!    levels = n_zlev
!    !-------------------------------------------------------------------------
!
!    !-------------------------------------------------------------------------
!    ! Attention: with use_constant_mixing=.true. there is no application of
!    ! convective mixing parameters in case of instability
!    ! max_vert_diff_veloc / max_vert_diff_trac
!    ! control of convective and constant mixing should be independent
!
!    grav_rho          = grav/OceanReferenceDensity
!
!    onem1_lambda_d    = 1.0_wp-lambda_diff
!    onem1_lambda_v    = 1.0_wp-lambda_visc
!    dv_rich           = richardson_tracer
!    av_rich           = richardson_veloc
!    v10mexp_3         = 1.0_wp/v10m_ref**3
!    wma_pd            = wma_diff * v10mexp_3   !  scaled wind-mixing amplitude for diffusion
!    wma_pv            = wma_visc * v10mexp_3   !  scaled wind-mixing amplitude for viscosity
!
!    dv_wind(:,1:levels,:) = 0.0_wp
!    av_wind(:,1:levels,:) = 0.0_wp
!
!    loc_eps = dbl_eps
!
!!ICON_OMP_PARALLEL PRIVATE(salinity)
!    salinity(1:levels) = sal_ref
!!ICON_OMP_DO PRIVATE(start_index, end_index, jc, levels, jk, pressure, rho_up, rho_down, &
!!ICON_OMP vert_velocity_shear, tracer_index, diffusion_weight, decay_wind_depth, wind_param, &
!!ICON_OMP jk_max, vdgfac_bot, vdensgrad_inter, dv_old, dv_back) ICON_OMP_DEFAULT_SCHEDULE
!    DO blockNo = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, blockNo, start_index, end_index)
!      richardson_no    (:,:, blockNo) = 0.0_wp
!      vert_density_grad(:,:, blockNo) = 0.0_wp
!      DO jc = start_index, end_index
!
!        levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
!        IF (levels < min_dolic) CYCLE
!
!        IF(no_tracer >= 2) THEN
!            salinity(1:levels) = ocean_state%p_prog(nold(1))%tracer(jc,1:levels,blockNo,2)
!        ENDIF
!
!        !--------------------------------------------------------
!        ! calculate density gradient for stability detection for PP-scheme and convection parameterization:
!        !  - pressure at interface between upper and lower layer
!        !  - S and T taken from upper and lower level, i.e. 2 times calculation of density per layer
!
!        !  - old formulation without including z in reference interface level
!        pressure(2:levels) = patch_3d%p_patch_1d(1)%zlev_i(2:levels) * OceanReferenceDensity * sitodbar
!        rho_up(1:levels-1)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,1:levels-1,blockNo,1), &
!          & salinity(1:levels-1), pressure(2:levels), levels-1)
!        rho_down(2:levels)  = calculate_density_onColumn(ocean_state%p_prog(nold(1))%tracer(jc,2:levels,blockNo,1), &
!          & salinity(2:levels), pressure(2:levels), levels-1)
!
!        DO jk = 2, levels
!          ! division by dz**2 is omitted in this calculation of velocity shear: shear = (d_vn)**2
!          vert_velocity_shear = loc_eps + &
!            & SUM((ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x - ocean_state%p_diag%p_vn(jc,jk,blockNo)%x)**2)
!          ! d_rho/dz - full density gradient necessary for wind mixing, stability formula, mixed layer depth calculation
!          vert_density_grad(jc,jk,blockNo) = (rho_down(jk) - rho_up(jk-1)) *  &
!            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,blockNo)
!          ! Ri = g/OceanReferenceDensity * dz * d_rho/(d_vn)**2
!          richardson_no(jc,jk,blockNo) = MAX(patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,blockNo) * grav_rho * &
!            &                           (rho_down(jk) - rho_up(jk-1)) / vert_velocity_shear, 0.0_wp)
!        END DO ! levels
!
!        IF (use_wind_mixing) THEN
!
!          ! wind-mixing at surface, eq. (15) of Marsland et al., 2003
!
!          ! reduced wind-mixing under sea ice, following MPIOM
!          IF (use_reduced_mixing_under_ice) THEN
!            dv_wind(jc,1,blockNo) = wma_pv * (1.0_wp - concsum(jc,blockNo)) * fu10(jc,blockNo)**3
!          ELSE
!            dv_wind(jc,1,blockNo) = wma_pv * (1.0_wp - concsum(jc,blockNo))**2 * fu10(jc,blockNo)**3
!          ENDIF
!
!          ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003
!          DO jk = 2, levels
!
!            decay_wind_depth   = EXP(-patch_3d%p_patch_1d(1)%del_zlev_m(jk-1)/z0_wind)
!
!            ! lambda_wind: default changed to 0.05, with 0.03 omip-r2b4 aborted
!            !   - in MPIOM it is 0.05, in Marsland et al. it was 0.03,
!            !   - for strong winds it might be necessary to increase it
!            wind_param         = lambda_wind * patch_3d%p_patch_1d(1)%inv_del_zlev_m(jk)
!
!            ! vertical interpolation of density gradient
!            !  - at mid-depth jk the density gradients at interfaces from above (jk) and below (jk+1) are used
!            !  - at bottom level the density gradient from below is zero, half density gradient is used
!            jk_max             = MIN(jk+1,levels)
!            vdgfac_bot(jk)     = 1.0_wp
!            vdgfac_bot(levels) = 0.0_wp
!            vdensgrad_inter  = 0.5_wp*(vert_density_grad(jc,jk,blockNo)+vdgfac_bot(jk)*vert_density_grad(jc,jk_max,blockNo))
!
!            dv_wind(jc,jk,blockNo) =  dv_wind(jc,jk-1,blockNo)*decay_wind_depth*wind_param / (wind_param + vdensgrad_inter)
!
!            ! cut unphysically negative wind induced mixing
!            dv_wind(jc,jk,blockNo) =  MAX(dv_wind(jc,jk,blockNo),0.0_wp)
!            ! cut overshoots more than convection maximum - must not be set to low values
!          ! IF (max_vert_diff_trac .GT. params_oce%a_tracer_v_back(1)) &
!          !   &  dv_wind(jc,jk,blockNo) =  MIN(dv_wind(jc,jk,blockNo),0.0_wp)
!
!          END DO
!
!        END IF  ! use_wind_mixing
!
!      ENDDO !  block index
!
!      DO tracer_index = 1, no_tracer
!        DO jc = start_index, end_index
!          levels = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
!
!          DO jk = 2, levels
!
!            ! calculate the richardson diffusion using the eddy diffusion relaxation term lambda_diff
!            !  - a_tracer_v is relaxed to the last pp-value to avoid relaxing to convection value max_vert_diff_trac
!
!              dv_old  = params_oce%a_tracer_v(jc,jk,blockNo,tracer_index)
!              dv_back = params_oce%a_tracer_v_back(tracer_index)
!              params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = &
!                &    onem1_lambda_d*MIN(dv_old, dv_rich + dv_wind(jc,jk,blockNo)) &
!                &  + lambda_diff*(dv_rich/((1.0_wp + crd*richardson_no(jc,jk,blockNo))**3) + dv_back + dv_wind(jc,jk,blockNo))
!
!
!         !  ! clalculate the richardson diffusion - old arithmetic
!         !  IF (use_pp_scheme) THEN
!         !    params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = &
!         !      & params_oce%a_tracer_v_back(tracer_index) + &
!         !      & dv_rich / ((1.0_wp + crd *                 &
!         !      & richardson_no(jc,jk,blockNo))**3)
!         !  ENDIF
!
!                ! #slo# ensure that pp is active for low values of max_vert_diff_trac
!                dv_old = params_oce%a_tracer_v(jc,jk,blockNo,tracer_index)
!                params_oce%a_tracer_v(jc,jk,blockNo,tracer_index) = MAX(  &
!                  ! #slo# Attention: convection_InstabilityThreshold<0 in old formulation - used with reverted sign
!                  &  max_vert_diff_trac * (-convection_InstabilityThreshold-vert_density_grad(jc,jk,blockNo)) /     &
!                  &                       (-convection_InstabilityThreshold+ABS(vert_density_grad(jc,jk,blockNo))), dv_old)
!
!          ENDDO ! levels
!        ENDDO !  block index
!      ENDDO ! tracer_index
!
!    END DO ! blocks
!!ICON_OMP_END_DO
!
!
!    !--------------------------------------------
!    ! Calculate params_oce%A_veloc_v:
!    ! use mean values between the two cells; change to min, max if required
!!ICON_OMP_DO PRIVATE(start_index, end_index, je, cell_1_idx, cell_1_block, cell_2_idx, cell_2_block, &
!!ICON_OMP levels, jk,  mean_density_differ_edge, mean_richardson_e, decay_wind_depth, wind_param, &
!!ICON_OMP jk_max, vdgfac_bot, densgrad_k, densgrad_kp1, vdensgrad_inter, av_old, av_back) ICON_OMP_DEFAULT_SCHEDULE
!    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
!      DO je = start_index, end_index
!
!        levels = patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
!        IF (levels < min_dolic) CYCLE
!
!        cell_1_idx = p_patch%edges%cell_idx(je,blockNo,1)
!        cell_1_block = p_patch%edges%cell_blk(je,blockNo,1)
!        cell_2_idx = p_patch%edges%cell_idx(je,blockNo,2)
!        cell_2_block = p_patch%edges%cell_blk(je,blockNo,2)
!
!        IF (use_wind_mixing) THEN
!
!          ! TODO: the following expects equally sized cells
!          ! TODO: think about boundary values on land
!          fu10_e = 0.5_wp * (   fu10(cell_1_idx,cell_1_block) +    fu10(cell_2_idx,cell_2_block))
!          !fu10_e = 10.0_wp
!          conc_e = 0.5_wp * (concsum(cell_1_idx,cell_1_block) + concsum(cell_2_idx,cell_2_block))
!
!          ! wind-mixing at surface, eq. (15) of Marsland et al., 2003
!
!          ! reduced wind-mixing under sea ice, following MPIOM
!          IF (use_reduced_mixing_under_ice) THEN
!            av_wind(je,1,blockNo) = wma_pd * (1.0_wp - conc_e)    * fu10_e**3
!          ELSE
!            av_wind(je,1,blockNo) = wma_pd * (1.0_wp - conc_e)**2 * fu10_e**3
!          ENDIF
!
!          ! exponential decay of wind-mixing, eq. (16) of Marsland et al., 2003, edges
!          DO jk = 2, levels
!            decay_wind_depth = EXP(-patch_3d%p_patch_1d(1)%del_zlev_m(jk-1)/z0_wind)
!            wind_param       = lambda_wind * patch_3d%p_patch_1d(1)%inv_del_zlev_m(jk)
!
!            ! vertical interpolation of density gradient
!            !  - at mid-depth jk the density gradients at interfaces from above (jk) and below (jk+1) are used
!            !  - at bottom level the density gradient from below is zero, half density gradient is used
!            jk_max             = MIN(jk+1,levels)
!            vdgfac_bot(jk)     = 1.0_wp
!            vdgfac_bot(levels) = 0.0_wp
!     !      vdensgrad_inter    = 0.5_wp*(vert_density_grad(jc,jk,blockNo)+vdgfac_bot(jk)*vert_density_grad(jc,jk_max,blockNo))
!            densgrad_k       = 0.5_wp * (vert_density_grad(cell_1_idx,jk,cell_1_block) + &
!              & vert_density_grad(cell_2_idx,jk,cell_2_block))
!            densgrad_kp1     = 0.5_wp*vdgfac_bot(jk) * (vert_density_grad(cell_1_idx,jk_max,cell_1_block) &
!              & + vert_density_grad(cell_2_idx,jk_max,cell_2_block))
!            vdensgrad_inter  = 0.5_wp*(densgrad_k + densgrad_kp1)
!
!            av_wind(je,jk,blockNo) =  av_wind(je,jk-1,blockNo)*decay_wind_depth*wind_param / (wind_param + vdensgrad_inter)
!
!            ! cut unphysically negative wind induced mixing
!            av_wind(je,jk,blockNo) =  MAX(av_wind(je,jk,blockNo),0.0_wp)
!            ! cut overshoots more than convection maximum - must not be set to low values
!    !       IF (max_vert_diff_veloc .GT. params_oce%a_veloc_v_back) &
!    !         &  av_wind(je,jk,blockNo) =  MIN(av_wind(je,jk,blockNo),max_vert_diff_veloc)
!          END DO
!
!        END IF  ! use_wind_mixing
!
!        DO jk = 2, levels
!
!          ! Set to zero for land + boundary locations edges
!!           IF (patch_3d%lsm_e(je,jk,blockNo) > sea) THEN
!!             params_oce%a_veloc_v(je,jk,blockNo) = 0.0_wp
!!           ELSE
!
!          ! TODO: the following expects equally sized cells
!          ! compute quantities at edges
!          mean_density_differ_edge = 0.5_wp * (vert_density_grad(cell_1_idx,jk,cell_1_block) &
!            & + vert_density_grad(cell_2_idx,jk,cell_2_block))
!          mean_richardson_e   = 0.5_wp * (richardson_no(cell_1_idx,jk,cell_1_block) + &
!            & richardson_no(cell_2_idx,jk,cell_2_block))
!
!          ! calculate the richardson viscosity using the eddy viscosity relaxation term lambda_visc
!          !  - a_veloc_v is relaxed to the last pp-value to avoid relaxing to convection value max_vert_diff_veloc
!!           IF (use_pp_scheme) THEN
!
!            av_old  = params_oce%a_veloc_v(je,jk,blockNo)
!            av_back = params_oce%a_veloc_v_back
!            params_oce%a_veloc_v(je,jk,blockNo) = &
!              &    onem1_lambda_v*MIN(av_old, av_rich + av_wind(je,jk,blockNo)) &
!              &  + lambda_visc*(av_rich/((1.0_wp + crv*mean_richardson_e)**2) + av_back + av_wind(je,jk,blockNo))
!!           ENDIF
!
!          ! turn on convection - no convection for velocity as in mpiom
!!           IF (use_convection) THEN
!! 
!!             ! MPIOM style of convection in PP-scheme: viscosity
!!             IF (use_mpiom_pp_form) THEN
!!               ! #slo# ensure that pp is active for low values of max_vert_diff_trac
!!               av_old = params_oce%a_veloc_v(je,jk,blockNo)
!!               params_oce%a_veloc_v(je,jk,blockNo) = MAX(  &
!!                 ! #slo# Attention: convection_InstabilityThreshold<0 in old formulation - used with reverted sign
!!                 &  max_vert_diff_veloc * (-convection_InstabilityThreshold-mean_density_differ_edge) /     &
!!                 &                       (-convection_InstabilityThreshold+ABS(mean_density_differ_edge)), av_old)
!! 
!!             ELSE ! do not use_mpiom_pp_form
!! 
!!               ! turn on convection
!!               IF (mean_density_differ_edge <  convection_InstabilityThreshold) THEN
!!                 ! #slo# ensure that pp is active for low values of max_vert_diff_trac
!!                 !params_oce%a_veloc_v(je,jk,blockNo) = max_vert_diff_veloc
!!                 params_oce%a_veloc_v(je,jk,blockNo) = MAX(max_vert_diff_veloc,params_oce%a_veloc_v(je,jk,blockNo))
!!               ELSE
!
!!                 IF (mean_density_differ_edge < RichardsonDiffusion_threshold) THEN
!!                   diffusion_weight =  &
!!                     & (mean_density_differ_edge - convection_InstabilityThreshold) / &
!!                     & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)
!! 
!!                   ! richardson diffusion from above
!!                   av_old = params_oce%a_veloc_v(je,jk,blockNo)
!!                   params_oce%a_veloc_v(je,jk,blockNo) = &
!!                     & max_vert_diff_veloc * (1.0_wp - diffusion_weight) + &
!!                     & av_old * diffusion_weight
!! 
!!                 ENDIF  ! grad<RichThreshold
!!               ENDIF ! grad<convThreshold
!!             ENDIF ! use_mpiom_pp_form
!!           ENDIF ! use_convection
!
!        END DO ! jk = 2, levels
!      ENDDO ! je = start_index, end_index
!    ENDDO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!!ICON_OMP_END_DO NOWAIT
!!ICON_OMP_END_PARALLEL
!
!    ! Sync the results, the A_tracer_v is only for checking
!    ! DO tracer_index = 1, no_tracer
!    !   CALL sync_patch_array(sync_c,p_patch,params_oce%a_tracer_v(:,:,:,tracer_index))
!    ! END DO
!    ! CALL sync_patch_array(sync_e,p_patch,params_oce%a_veloc_v(:,:,:))
!
!    !---------DEBUG DIAGNOSTICS-------------------------------------------
!    idt_src=4  ! output print levels (1-5, fix)
!    CALL dbg_print('UpdPar: p_vn%x(1)    ',ocean_state%p_diag%p_vn%x(1),str_module,idt_src,in_subset=p_patch%cells%owned)
!  ! CALL dbg_print('UpdPar: p_vn%x(2)    ',ocean_state%p_diag%p_vn%x(2),str_module,idt_src,in_subset=p_patch%cells%owned)
!    CALL dbg_print('UpdPar: windMix Diff ',dv_wind                     ,str_module,idt_src,in_subset=p_patch%cells%owned)
!    CALL dbg_print('UpdPar: windMix Visc ',av_wind                     ,str_module,idt_src,in_subset=p_patch%edges%owned)
!    CALL dbg_print('UpdPar: VertDensGrad ',vert_density_grad           ,str_module,idt_src,in_subset=p_patch%cells%owned)
!    CALL dbg_print('UpdPar: Richardson No',richardson_no               ,str_module,idt_src,in_subset=p_patch%cells%owned)
!    CALL dbg_print('UpdPar: windsp. fu10 ',fu10                        ,str_module,idt_src,in_subset=p_patch%cells%owned)
!    idt_src=5  ! output print levels (1-5, fix)
!    DO tracer_index = 1, no_tracer
!      CALL dbg_print('UpdPar FinalTracerMixing'  ,params_oce%a_tracer_v(:,:,:,tracer_index), str_module,idt_src, &
!        & in_subset=p_patch%cells%owned)
!    ENDDO
!    CALL dbg_print('UpdPar FinalVelocMixing'   ,params_oce%a_veloc_v     ,str_module,idt_src, &
!      & in_subset=p_patch%edges%owned)
!    !---------------------------------------------------------------------
!
!  END SUBROUTINE update_physics_parameters_MPIOM_PP_scheme
!  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !Subroutine computes the horizontal diffusive flux of an arbitrary tracer.
   SUBROUTINE tracer_diffusion_horz_local(patch_3D, trac_in, p_os, diff_flx, k_t, subset_range)
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)              :: trac_in(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET :: p_os
    REAL(wp), INTENT(inout)           :: diff_flx(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp), OPTIONAL                :: k_t(:,:,:) !mixing coefficient for tracer    
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    !
    !Local variables
    INTEGER :: level, blockNo, edge_index
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: start_edge_index, end_edge_index
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_ocediffusion:tracer_diffusion_horz')
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    !-------------------------------------------------------------------------------

  IF(PRESENT(k_t))THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, edge_index, level, &
!ICON_OMP il_c1, ib_c1, il_c2, ib_c2) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      diff_flx(:,:,blockNo) = 0.0_wp
      DO edge_index = start_edge_index, end_edge_index
        !Get indices of two adjacent triangles
        il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
        ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
        il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
        ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)

        DO level=1,  patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)

          diff_flx(edge_index,level,blockNo) = &
            &   k_t(edge_index,level,blockNo) &
            & * patch_3D%p_patch_1d(1)%prism_thick_e(edge_index,level,blockNo)  &
            & * (trac_in(il_c2,level,ib_c2) - trac_in(il_c1,level,ib_c1))       &
            & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)

        ENDDO

      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF (PRESENT(subset_range)) THEN
      IF (.NOT. subset_range%is_in_domain) &
        & CALL sync_patch_array(sync_e, patch_2D, diff_flx)
    ENDIF    

  ELSEIF(.NOT.PRESENT(k_t))THEN  
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, edge_index, level, &
!ICON_OMP il_c1, ib_c1, il_c2, ib_c2) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      diff_flx(:,:,blockNo) = 0.0_wp
      DO edge_index = start_edge_index, end_edge_index
        !Get indices of two adjacent triangles
        il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
        ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
        il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
        ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)

        DO level=1,  patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)

          diff_flx(edge_index,level,blockNo) = &
            & patch_3D%p_patch_1d(1)%prism_thick_e(edge_index,level,blockNo)  &
            & * (trac_in(il_c2,level,ib_c2) - trac_in(il_c1,level,ib_c1))       &
            & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)

        ENDDO

      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF (PRESENT(subset_range)) THEN
      IF (.NOT. subset_range%is_in_domain) &
        & CALL sync_patch_array(sync_e, patch_2D, diff_flx)
    ENDIF 
  ENDIF  

  END SUBROUTINE tracer_diffusion_horz_local
  !-------------------------------------------------------------------------



END MODULE mo_ocean_physics
