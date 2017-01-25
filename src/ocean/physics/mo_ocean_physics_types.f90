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
MODULE mo_ocean_physics_types
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
    & velocity_TopWindMixing
    
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

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_physics'
  CHARACTER(LEN=12)           :: str_module    = 'ocePhysics  '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug

  ! Public interface
  PUBLIC :: t_ho_params
  PUBLIC :: WindAmplitude_at10m
  PUBLIC :: SeaIceConcentration
  PUBLIC :: WindMixingDecay
  PUBLIC :: WindMixingLevel

  !PUBLIC :: init_ho_physics
  PUBLIC :: construct_ho_params
  PUBLIC :: destruct_ho_params

  ! variables
  TYPE (t_var_list), PUBLIC :: ocean_params_list

  REAL(wp), POINTER :: WindAmplitude_at10m(:,:)  ! can be single precision
  REAL(wp), POINTER :: SeaIceConcentration(:,:)
  REAL(wp), POINTER :: WindMixingDecay(:), WindMixingLevel(:)


  ! Parameters below appear directly in the ocean model/equation. They are eventually
  ! dynamically updated by using the "ocean-physics" structure. #slo# - not yet
  TYPE t_ho_params

    ! diffusion coefficients for horizontal velocity, temp. and salinity, dim=(nproma,n_zlev,nblks_e)
    REAL(wp),POINTER ::     &
      & k_veloc_h(:,:,:),  & ! coefficient of horizontal velocity diffusion
      & k_tracer_h(:,:,:,:)  ! coefficient of horizontal tracer diffusion
    TYPE(t_onEdges_Pointer_3d_wp),ALLOCATABLE :: tracer_h_ptr(:)

    ! diffusion coefficients for vertical velocity, temp. and salinity, dim=(nproma,n_zlev+1,nblks_e)
    REAL(wp),POINTER ::     &
      & a_veloc_v(:,:,:),  & ! coefficient of vertical velocity diffusion
      & a_tracer_v(:,:,:,:)  ! coefficient of vertical tracer diffusion
    TYPE(t_onCells_HalfLevels_Pointer_wp), ALLOCATABLE :: tracer_v_ptr(:)

    !constant background values of coefficients above
    REAL(wp) :: k_veloc_h_back, &! coefficient of horizontal velocity diffusion
      & a_veloc_v_back   ! coefficient of vertical velocity diffusion

    REAL(wp),ALLOCATABLE ::     &
      & k_tracer_h_back(:),    & ! coefficient of horizontal tracer diffusion dim=no_tracer
      & a_tracer_v_back(:)       ! coefficient of vertical tracer diffusion dim=no_tracer

    REAL(wp),POINTER ::     &
      & k_tracer_isoneutral(:,:,:),  & ! coefficient of isoneutral tracer diffusion diffusion at cells
      & k_tracer_dianeutral(:,:,:)  ! coefficient of dianeutral tracer diffusion
   REAL(wp),POINTER ::     &
      & k_tracer_GM_kappa(:,:,:) ! coefficient of Gent-McWilliams mesoscale eddyparametrizations
    REAL(wp) :: bottom_drag_coeff

  END TYPE t_ho_params

  TYPE(t_ho_params),PUBLIC,TARGET :: v_params
  

CONTAINS

 
  !-------------------------------------------------------------------------
  !>
  !! Construction of arrays for ocean physics
  !!
  !! Construction of arrays for ocean physics ...
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07)
  !
  !
!<Optimize:inUse:initOnly>
  SUBROUTINE construct_ho_params(patch_2D, params_oce, ocean_restart_list)

    TYPE(t_patch),      INTENT(IN)    :: patch_2D
    TYPE (t_ho_params), INTENT(INOUT) :: params_oce
    TYPE (t_var_list),  INTENT(INOUT) :: ocean_restart_list

    ! Local variables
    INTEGER :: ist, i,jtrc
    INTEGER :: alloc_cell_blocks, nblks_e
    INTEGER :: datatype_flt

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = this_mod_name//':construct_ho_physics'

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'construct hydro ocean physics')

    CALL new_var_list(ocean_params_list, 'ocean_params_list', patch_id=patch_2D%id)
    CALL default_var_list_settings( ocean_params_list, lrestart=.FALSE.,  model_type='oce' )

    ! determine size of arrays
    alloc_cell_blocks = patch_2D%alloc_cell_blocks
    nblks_e = patch_2D%nblks_e

    CALL add_var(ocean_params_list, 'K_veloc_h', params_oce%k_veloc_h , grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('K_veloc_h', 'kg/kg', 'horizontal velocity diffusion', datatype_flt),&
      & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_physics"))

    CALL add_var(ocean_params_list, 'A_veloc_v', params_oce%a_veloc_v , grid_unstructured_edge,&
      & za_depth_below_sea_half, &
      & t_cf_var('A_veloc_v', 'kg/kg', 'vertical velocity diffusion', datatype_flt),&
      & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
      & ldims=(/nproma,n_zlev+1,nblks_e/),in_group=groups("oce_physics","oce_default"))


    !! Tracers
    IF ( no_tracer > 0 ) THEN
      CALL add_var(ocean_params_list, 'K_tracer_h', params_oce%k_tracer_h , &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('K_tracer_h', '', '1:temperature 2:salinity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e,no_tracer/), &
        & lcontainer=.TRUE., loutput=.FALSE., lrestart=.FALSE.)
      CALL add_var(ocean_params_list, 'A_tracer_v', params_oce%a_tracer_v , &
        & grid_unstructured_cell, za_depth_below_sea_half, &
        & t_cf_var('A_tracer_v', '', '1:temperature 2:salinity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev+1,alloc_cell_blocks,no_tracer/), &
        & lcontainer=.TRUE., loutput=.FALSE., lrestart=.FALSE.)

      ! Reference to individual tracer, for I/O

      ALLOCATE(params_oce%tracer_h_ptr(no_tracer))
      ALLOCATE(params_oce%tracer_v_ptr(no_tracer))
      DO jtrc = 1,no_tracer
        CALL add_ref( ocean_params_list, 'K_tracer_h',&
          & 'K_tracer_h_'//TRIM(oce_config%tracer_names(jtrc)),     &
          & params_oce%tracer_h_ptr(jtrc)%p,                             &
          & grid_unstructured_edge, za_depth_below_sea,               &
          & t_cf_var('K_tracer_h_'//TRIM(oce_config%tracer_names(jtrc)), &
          & 'kg/kg', &
          & TRIM(oce_config%tracer_longnames(jtrc))//'(K_tracer_h_)', &
          & datatype_flt), &
          & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_edge),&
          & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_physics"))
        CALL add_ref( ocean_params_list, 'A_tracer_v',&
          & 'A_tracer_v_'//TRIM(oce_config%tracer_names(jtrc)),     &
          & params_oce%tracer_h_ptr(jtrc)%p,                             &
          & grid_unstructured_cell, za_depth_below_sea_half,            &
          & t_cf_var('A_tracer_v_'//TRIM(oce_config%tracer_names(jtrc)), &
          & 'kg/kg', &
          & TRIM(oce_config%tracer_longnames(jtrc))//'(A_tracer_v)', &
          & datatype_flt), &
          & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
          & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_physics"))

      END DO
      !TODO     use the following code, if add_var support 1d arrays:
      !TODO     CALL add_var(ocean_params_list, 'K_tracer_h_back', params_oce%K_tracer_h_back , &
      !TODO     &            GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, &
      !TODO     &            t_cf_var('K_tracer_h_back', '', '1:temperature 2:salinity', datatype_flt),&
      !TODO     &            grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_EDGE),&
      !TODO     &            ldims=(/ no_tracer /))
      !TODO     CALL add_var(ocean_params_list, 'A_tracer_v_back', params_oce%A_tracer_v_back , &
      !TODO     &            GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      !TODO     &            t_cf_var('A_tracer_v_back', '', '1:temperature 2:salinity', datatype_flt),&
      !TODO     &            grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      !TODO     &            ldims=(/no_tracer/))
    ENDIF ! no_tracer > 0


    ALLOCATE(params_oce%k_tracer_h_back(no_tracer), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for horizontal background tracer diffusion failed')
    END IF

    ALLOCATE(params_oce%a_tracer_v_back(no_tracer), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation for vertical tracer background diffusion failed')
    END IF

   !IF(GMRedi_configuration==GMRedi_combined&
   !&.OR.GMRedi_configuration==GM_only.OR.GMRedi_configuration==Redi_only)THEN


     CALL add_var(ocean_params_list, 'k_tracer_isoneutral', params_oce%k_tracer_isoneutral, &
        & grid_unstructured_cell, za_depth_below_sea, &
        & t_cf_var('k_tracer_isoneutral at cells', '', '1:temperature 2:salinity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/), &
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.TRUE.)


      CALL add_var(ocean_params_list, 'k_tracer_dianeutral', params_oce%k_tracer_dianeutral, &
        & grid_unstructured_cell, za_depth_below_sea_half, &
        & t_cf_var('A_tracer_v', '', '1:temperature 2:salinity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/), &
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


     CALL add_var(ocean_params_list, 'k_tracer_GM_kappa', params_oce%k_tracer_GM_kappa, &
        & grid_unstructured_cell, za_depth_below_sea, &
        & t_cf_var('k_tracer_GM_kappa at cells', '', '1:temperature 2:salinity', datatype_flt),&
        & grib2_var(255, 255, 255, datatype_pack16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/), &
        & lcontainer=.TRUE., lrestart=.FALSE., loutput=.TRUE.)

    !ENDIF


    DO i=1,no_tracer
      params_oce%k_tracer_h_back(i)  = 0.0_wp
      params_oce%a_tracer_v_back(i)  = 0.0_wp
    END DO

    ALLOCATE(WindMixingDecay(1:n_zlev+1), WindMixingLevel(1:n_zlev+1))
    
  END SUBROUTINE construct_ho_params
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Destruction of arrays for ocean physics
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07)
  !
!<Optimize:inUse:initOnly>
  SUBROUTINE destruct_ho_params(params_oce)

    TYPE (t_ho_params), INTENT(inout) :: params_oce

    ! Local variables
    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = this_mod_name//':destruct_ho_physics'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'destruct hydro ocean physics')

    CALL delete_var_list(ocean_params_list)

    DEALLOCATE(params_oce%k_tracer_h_back, stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'deallocation for horizontal tracer &
        & background iffusion failed')
    END IF

    DEALLOCATE(params_oce%a_tracer_v_back, stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'deallocation for vertical background &
        & temperaure diffusion failed')
    END IF
    
    DEALLOCATE(WindMixingDecay, WindMixingLevel)
    
  END SUBROUTINE destruct_ho_params
  !-------------------------------------------------------------------------



END MODULE mo_ocean_physics_types
