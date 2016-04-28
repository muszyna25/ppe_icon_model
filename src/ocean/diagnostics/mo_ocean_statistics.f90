!>
!! Contains the main stepping routine the 3-dim hydrostatic ocean model.
!!
!! @author Peter Korn, Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Initial version by Stephan Lorenz (MPI-M), (2010-04).
!!   - renaming and adjustment of hydrostatic ocean model V1.0.3 to ocean domain and patch_oce
!!  Modification by Stephan Lorenz, MPI-M, 2010-10
!!   - new module mo_ocean_statistics including updated reconstructions
!
!
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
MODULE mo_ocean_statistics
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: max_char_length
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d, t_subset_range, t_patch_vert
  USE mo_grid_config,            ONLY: n_dom
  USE mo_grid_subset,            ONLY: get_index_range
  USE mo_sync,                   ONLY: sync_patch_array, sync_e, sync_c !, sync_v
  USE mo_ocean_nml,              ONLY: iforc_oce, No_Forcing, iswm_oce, n_zlev, no_tracer, &
    & diagnostics_level, &
    & eos_type, i_sea_ice, gibraltar, &
    & iforc_oce, No_Forcing
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_io_config,              ONLY: n_checkpoints
  USE mo_run_config,             ONLY: nsteps, dtime, ltimer, output_mode
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ext_data_types,         ONLY: t_external_data
  !USE mo_io_units,               ONLY: filename_max
  USE mo_datetime,               ONLY: t_datetime, print_datetime, add_time, datetime_to_string
  USE mo_ocean_types,              ONLY: t_hydro_ocean_state, t_hydro_ocean_acc, t_hydro_ocean_diag, &
    & t_hydro_ocean_prog
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_sea_ice,                ONLY: compute_mean_ice_statistics, reset_ice_statistics
  USE mo_sea_ice_types,          ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
    & t_sea_ice
  USE mo_ocean_physics_types,    ONLY: t_ho_params
  USE mo_name_list_output,       ONLY: write_name_list_output, istime4name_list_output
  USE mo_linked_list,            ONLY: t_list_element, find_list_element
  USE mo_var_list,               ONLY: print_var_list
  USE mo_mpi,                    ONLY: my_process_is_stdio
  USE mo_time_config,            ONLY: time_config
  USE mo_statistics
  USE mo_sea_ice_nml,            ONLY: i_ice_dyn
  USE mo_util_dbg_prnt,          ONLY: dbg_print
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: update_ocean_statistics
  PUBLIC :: compute_mean_ocean_statistics
  PUBLIC :: reset_ocean_statistics
  PUBLIC :: new_ocean_statistics  
  !-------------------------------------------------------------------------
  
CONTAINS
    
  !---------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE update_ocean_statistics(ocean_state,p_sfc_flx,cells,edges,verts,max_zlev,p_phys_param)
    TYPE(t_hydro_ocean_state), INTENT(inout) :: ocean_state
    TYPE(t_sfc_flx),           INTENT(inout) :: p_sfc_flx
    TYPE(t_ho_params),OPTIONAL,INTENT(IN)    :: p_phys_param
    TYPE(t_subset_range),      INTENT(IN)    :: cells,edges,verts
    INTEGER, INTENT(in)                      :: max_zlev
    
    INTEGER :: jtrc,i
    
    ! update ocean state accumulated values
    CALL add_fields(ocean_state%p_acc%h     , ocean_state%p_prog(nnew(1))%h, cells)
    CALL add_sqr_fields(ocean_state%p_acc%h_sqr , ocean_state%p_prog(nnew(1))%h, cells)
    CALL add_fields(ocean_state%p_acc%u     , ocean_state%p_diag%u         , cells)
    CALL add_fields(ocean_state%p_acc%v     , ocean_state%p_diag%v         , cells)
    CALL add_fields(ocean_state%p_acc%rhopot, ocean_state%p_diag%rhopot    , cells)
    DO jtrc=1,no_tracer
      CALL add_fields(ocean_state%p_acc%tracer(:,:,:,jtrc),ocean_state%p_prog(nnew(1))%tracer(:,:,:,jtrc),cells)
    END DO
    CALL add_fields(ocean_state%p_acc%u_vint        , ocean_state%p_diag%u_vint        , cells)
    CALL add_fields(ocean_state%p_acc%v_vint        , ocean_state%p_diag%v_vint        , cells)
    CALL add_fields(ocean_state%p_acc%w             , ocean_state%p_diag%w             , cells,levels=max_zlev+1)
    CALL add_fields(ocean_state%p_acc%div_mass_flx_c, ocean_state%p_diag%div_mass_flx_c, cells)
    CALL add_fields(ocean_state%p_acc%rho           , ocean_state%p_diag%rho           , cells)
    CALL add_fields(ocean_state%p_acc%mass_flx_e    , ocean_state%p_diag%mass_flx_e    , edges)
    CALL add_fields(ocean_state%p_acc%vort          , ocean_state%p_diag%vort          , verts,levels=max_zlev)
    CALL add_fields(ocean_state%p_acc%kin           , ocean_state%p_diag%kin           , cells)
    CALL add_verticalSum_field(ocean_state%p_acc%edgeFlux_total, ocean_state%p_diag%mass_flx_e, edges)
    ! CALL dbg_print('nnew(1))%vn', ocean_state%p_prog(nnew(1))%vn ,"statistics",1,in_subset=edges)
    ! CALL dbg_print('edgeFlux_total_acc', ocean_state%p_acc%edgeFlux_total, "statistics",1,in_subset=edges)


    IF (PRESENT(p_phys_param)) THEN
      ! physics
      DO jtrc=1,no_tracer
        CALL add_fields(ocean_state%p_acc%k_tracer_h(:,:,:,jtrc),p_phys_param%k_tracer_h(:,:,:,jtrc),edges)
        CALL add_fields(ocean_state%p_acc%a_tracer_v(:,:,:,jtrc),p_phys_param%a_tracer_v(:,:,:,jtrc),cells)
      END DO
      CALL add_fields(ocean_state%p_acc%k_veloc_h(:,:,:),p_phys_param%k_veloc_h(:,:,:),edges)
      CALL add_fields(ocean_state%p_acc%a_veloc_v(:,:,:),p_phys_param%a_veloc_v(:,:,:),edges)
    END IF



    IF (iforc_oce > No_Forcing) THEN
      ! update forcing accumulated values
      CALL add_fields(p_sfc_flx%topBoundCond_windStress_u_acc   , p_sfc_flx%topBoundCond_windStress_u   , cells)
      CALL add_fields(p_sfc_flx%topBoundCond_windStress_v_acc   , p_sfc_flx%topBoundCond_windStress_v   , cells)
      CALL add_fields(p_sfc_flx%HeatFlux_ShortWave_acc          , p_sfc_flx%HeatFlux_ShortWave          , cells)
      CALL add_fields(p_sfc_flx%HeatFlux_LongWave_acc           , p_sfc_flx%HeatFlux_LongWave           , cells)
      CALL add_fields(p_sfc_flx%HeatFlux_Sensible_acc           , p_sfc_flx%HeatFlux_Sensible           , cells)
      CALL add_fields(p_sfc_flx%HeatFlux_Latent_acc             , p_sfc_flx%HeatFlux_Latent             , cells)
      CALL add_fields(p_sfc_flx%HeatFlux_Total_acc              , p_sfc_flx%HeatFlux_Total              , cells)
      CALL add_fields(p_sfc_flx%FrshFlux_Precipitation_acc      , p_sfc_flx%FrshFlux_Precipitation      , cells)
      CALL add_fields(p_sfc_flx%FrshFlux_SnowFall_acc           , p_sfc_flx%FrshFlux_SnowFall           , cells)
      CALL add_fields(p_sfc_flx%FrshFlux_Evaporation_acc        , p_sfc_flx%FrshFlux_Evaporation        , cells)
      CALL add_fields(p_sfc_flx%FrshFlux_Runoff_acc             , p_sfc_flx%FrshFlux_Runoff             , cells)
      CALL add_fields(p_sfc_flx%FrshFlux_TotalSalt_acc          , p_sfc_flx%FrshFlux_TotalSalt          , cells)
      CALL add_fields(p_sfc_flx%FrshFlux_TotalOcean_acc         , p_sfc_flx%FrshFlux_TotalOcean         , cells)
      CALL add_fields(p_sfc_flx%FrshFlux_TotalIce_acc           , p_sfc_flx%FrshFlux_TotalIce           , cells)
      CALL add_fields(p_sfc_flx%FrshFlux_VolumeIce_acc          , p_sfc_flx%FrshFlux_VolumeIce          , cells)
      CALL add_fields(p_sfc_flx%FrshFlux_VolumeTotal_acc        , p_sfc_flx%FrshFlux_VolumeTotal        , cells)
      CALL add_fields(p_sfc_flx%HeatFlux_Relax_acc              , p_sfc_flx%HeatFlux_Relax              , cells)
      CALL add_fields(p_sfc_flx%FrshFlux_Relax_acc              , p_sfc_flx%FrshFlux_Relax              , cells)
      CALL add_fields(p_sfc_flx%data_surfRelax_Temp_acc(:,:)    , p_sfc_flx%data_surfRelax_Temp(:,:)    , cells)
      CALL add_fields(p_sfc_flx%data_surfRelax_Salt_acc(:,:)    , p_sfc_flx%data_surfRelax_Salt(:,:)    , cells)
      CALL add_fields(p_sfc_flx%topBoundCond_Salt_vdiff_acc(:,:), p_sfc_flx%topBoundCond_Temp_vdiff(:,:), cells)
      CALL add_fields(p_sfc_flx%topBoundCond_Salt_vdiff_acc(:,:), p_sfc_flx%topBoundCond_Temp_vdiff(:,:), cells)
    ENDIF
    
  END SUBROUTINE update_ocean_statistics
  !---------------------------------------------------------------------
  
  !---------------------------------------------------------------------
  SUBROUTINE compute_mean_ocean_statistics(p_acc,p_sfc_flx,nsteps_since_last_output)
    TYPE(t_hydro_ocean_acc), INTENT(inout) :: p_acc
    TYPE(t_sfc_flx),         INTENT(inout) :: p_sfc_flx
    INTEGER,INTENT(in)                     :: nsteps_since_last_output
    
    
    !TODO [ram] trigger aoutput wrt to multiple output intervals
    !TODO [ram] CALL collect_group(TRIM(grp_name), grp_vars_fg_sfc, ngrp_vars_fg_sfc,    &
    !TODO [ram]       &                loutputvars_only=.FALSE.,lremap_lonlat=.FALSE.)
    
    !TODO [ram]   IF (ALLOCATED(output_file)) THEN
    !TODO [ram]      DO i = 1, nvar_lists
    !TODO [ram]         element => var_lists(i)%p%first_list_element
    !TODO [ram]         DO
    !TODO [ram]           IF(.NOT. ASSOCIATED(element)) EXIT
    !TODO [ram]           if (.not. element%isteptype==TSTEP_ACCUM) cycle
    !TODO [ram]           DO i = 1, SIZE(output_file)
    !TODO [ram]              if (istime4output*output_file(i)%out_event) then
    !TODO [ram]                 if (one_of(element%name, output_file(i)%namelist%ml_varlist(1:output_file(i)%num_vars)) then
    !TODO [ram]               ! do accumulation for variable varlist(j)
    !TODO [ram]                   element%field%rptr = element%field%rptr /
    !TODO [ram]            end if
    !TODO [ram]         end if
    !TODO [ram]      end DO
    !TODO [ram]   end IF
    

      
!ICON_OMP_PARALLEL
!ICON_OMP_WORKSHARE
    p_acc%h                         = p_acc%h                        /REAL(nsteps_since_last_output,wp)
    p_acc%h_sqr                     = p_acc%h_sqr                    /REAL(nsteps_since_last_output,wp)
    p_acc%u                         = p_acc%u                        /REAL(nsteps_since_last_output,wp)
    p_acc%v                         = p_acc%v                        /REAL(nsteps_since_last_output,wp)
    p_acc%rhopot                    = p_acc%rhopot                   /REAL(nsteps_since_last_output,wp)
    p_acc%u_vint                    = p_acc%u_vint                   /REAL(nsteps_since_last_output,wp)
    p_acc%v_vint                    = p_acc%v_vint                   /REAL(nsteps_since_last_output,wp)
    p_acc%edgeFlux_total                   = p_acc%edgeFlux_total                  /REAL(nsteps_since_last_output,wp)
    p_acc%w                         = p_acc%w                        /REAL(nsteps_since_last_output,wp)
    p_acc%div_mass_flx_c            = p_acc%div_mass_flx_c           /REAL(nsteps_since_last_output,wp)
    p_acc%rho                       = p_acc%rho                      /REAL(nsteps_since_last_output,wp)
    p_acc%mass_flx_e                = p_acc%mass_flx_e               /REAL(nsteps_since_last_output,wp)
    p_acc%vort                      = p_acc%vort                     /REAL(nsteps_since_last_output,wp)
    p_acc%kin                       = p_acc%kin                      /REAL(nsteps_since_last_output,wp)
    p_acc%k_veloc_h                 = p_acc%k_veloc_h                /REAL(nsteps_since_last_output)
    p_acc%a_veloc_v                 = p_acc%a_veloc_v                /REAL(nsteps_since_last_output)
!ICON_OMP_END_WORKSHARE
    
    IF(no_tracer>0)THEN
!ICON_OMP_WORKSHARE
      p_acc%tracer                    = p_acc%tracer                   /REAL(nsteps_since_last_output,wp)
      p_acc%k_tracer_h                = p_acc%k_tracer_h               /REAL(nsteps_since_last_output)
      p_acc%a_tracer_v                = p_acc%a_tracer_v               /REAL(nsteps_since_last_output)
!ICON_OMP_END_WORKSHARE
    ENDIF
      
    IF (iforc_oce > No_Forcing) THEN
!ICON_OMP_WORKSHARE
      p_sfc_flx%topBoundCond_windStress_u_acc = p_sfc_flx%topBoundCond_windStress_u_acc/REAL(nsteps_since_last_output,wp)
      p_sfc_flx%topBoundCond_windStress_v_acc = p_sfc_flx%topBoundCond_windStress_v_acc/REAL(nsteps_since_last_output,wp)
      p_sfc_flx%HeatFlux_ShortWave_acc        = p_sfc_flx%HeatFlux_ShortWave_acc       /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%HeatFlux_LongWave_acc         = p_sfc_flx%HeatFlux_LongWave_acc        /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%HeatFlux_Sensible_acc         = p_sfc_flx%HeatFlux_Sensible_acc        /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%HeatFlux_Latent_acc           = p_sfc_flx%HeatFlux_Latent_acc          /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%HeatFlux_Total_acc            = p_sfc_flx%HeatFlux_Total_acc           /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%FrshFlux_Precipitation_acc    = p_sfc_flx%FrshFlux_Precipitation_acc   /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%FrshFlux_SnowFall_acc         = p_sfc_flx%FrshFlux_SnowFall_acc        /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%FrshFlux_Evaporation_acc      = p_sfc_flx%FrshFlux_Evaporation_acc     /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%FrshFlux_Runoff_acc           = p_sfc_flx%FrshFlux_Runoff_acc          /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%FrshFlux_TotalSalt_acc        = p_sfc_flx%FrshFlux_TotalSalt_acc       /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%FrshFlux_TotalOcean_acc       = p_sfc_flx%FrshFlux_TotalOcean_acc      /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%FrshFlux_TotalIce_acc         = p_sfc_flx%FrshFlux_TotalIce_acc        /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%FrshFlux_Relax_acc            = p_sfc_flx%FrshFlux_Relax_acc           /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%FrshFlux_VolumeIce_acc        = p_sfc_flx%FrshFlux_VolumeIce_acc       /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%FrshFlux_VolumeTotal_acc      = p_sfc_flx%FrshFlux_VolumeTotal_acc     /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%HeatFlux_Relax_acc            = p_sfc_flx%HeatFlux_Relax_acc           /REAL(nsteps_since_last_output,wp)
!ICON_OMP_END_WORKSHARE

      IF(no_tracer>0)THEN
!ICON_OMP_WORKSHARE
        p_sfc_flx%data_surfRelax_Temp_acc     = p_sfc_flx%data_surfRelax_Temp_acc        /REAL(nsteps_since_last_output,wp)
        p_sfc_flx%topBoundCond_Temp_vdiff_acc = p_sfc_flx%topBoundCond_Temp_vdiff_acc    /REAL(nsteps_since_last_output,wp)
!ICON_OMP_END_WORKSHARE
      ENDIF
      IF(no_tracer>1)THEN
!ICON_OMP_WORKSHARE
        p_sfc_flx%topBoundCond_Salt_vdiff_acc = p_sfc_flx%topBoundCond_Salt_vdiff_acc    /REAL(nsteps_since_last_output,wp)
!ICON_OMP_END_WORKSHARE
      ENDIF
    ENDIF ! iforc_oce > No_Forcing
!ICON_OMP_END_PARALLEL
      
  END SUBROUTINE compute_mean_ocean_statistics
  !---------------------------------------------------------------------
  
  !---------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE reset_ocean_statistics(p_acc,p_diag,p_sfc_flx,nsteps_since_last_output)
    TYPE(t_hydro_ocean_acc), INTENT(inout)  :: p_acc
    TYPE(t_hydro_ocean_diag), INTENT(inout) :: p_diag
    TYPE(t_sfc_flx),         INTENT(inout)  :: p_sfc_flx
    INTEGER,OPTIONAL,        INTENT(inout)  :: nsteps_since_last_output
    
    IF (PRESENT(nsteps_since_last_output)) nsteps_since_last_output        = 0
    
    p_acc%h                         = 0.0_wp
    p_acc%h_sqr                     = 0.0_wp
    p_acc%u                         = 0.0_wp
    p_acc%v                         = 0.0_wp
    p_acc%rhopot                    = 0.0_wp
    p_acc%u_vint                    = 0.0_wp
    p_acc%v_vint                    = 0.0_wp
    p_acc%edgeFlux_total                   = 0.0_wp
    p_acc%w                         = 0.0_wp
    p_acc%div_mass_flx_c            = 0.0_wp
    p_acc%rho                       = 0.0_wp
    p_acc%mass_flx_e                = 0.0_wp
    p_acc%vort                      = 0.0_wp
    p_acc%kin                       = 0.0_wp
    p_acc%k_veloc_h                 = 0.0_wp
    p_acc%a_veloc_v                 = 0.0_wp
    p_sfc_flx%topBoundCond_windStress_u_acc = 0.0_wp
    p_sfc_flx%topBoundCond_windStress_v_acc = 0.0_wp
    
    IF (no_tracer>0) THEN
      p_acc%tracer                    = 0.0_wp
      p_acc%k_tracer_h                = 0.0_wp
      p_acc%a_tracer_v                = 0.0_wp
      
      p_sfc_flx%HeatFlux_ShortWave_acc          = 0.0_wp
      p_sfc_flx%HeatFlux_LongWave_acc           = 0.0_wp
      p_sfc_flx%HeatFlux_Sensible_acc           = 0.0_wp
      p_sfc_flx%HeatFlux_Latent_acc             = 0.0_wp
      p_sfc_flx%HeatFlux_Total_acc              = 0.0_wp
      p_sfc_flx%topBoundCond_Temp_vdiff_acc     = 0.0_wp
      p_sfc_flx%HeatFlux_Relax_acc              = 0.0_wp
      IF (no_tracer>1) THEN
        p_sfc_flx%FrshFlux_Precipitation_acc    = 0.0_wp
        p_sfc_flx%FrshFlux_SnowFall_acc         = 0.0_wp
        p_sfc_flx%FrshFlux_Evaporation_acc      = 0.0_wp
        p_sfc_flx%FrshFlux_Runoff_acc           = 0.0_wp
        p_sfc_flx%FrshFlux_TotalSalt_acc        = 0.0_wp
        p_sfc_flx%FrshFlux_TotalOcean_acc       = 0.0_wp
        p_sfc_flx%FrshFlux_TotalIce_acc         = 0.0_wp
        p_sfc_flx%FrshFlux_Relax_acc            = 0.0_wp
        p_sfc_flx%FrshFlux_VolumeIce_acc        = 0.0_wp
        p_sfc_flx%FrshFlux_VolumeTotal_acc      = 0.0_wp
        p_sfc_flx%topBoundCond_Salt_vdiff_acc   = 0.0_wp
        p_sfc_flx%data_surfRelax_Temp_acc = 0.0_wp
        p_sfc_flx%data_surfRelax_Salt_acc = 0.0_wp
      ENDIF
    ENDIF

    ! reset mixed layer depth to zero
    p_diag%mld = 0.0_wp
  END SUBROUTINE reset_ocean_statistics
  !---------------------------------------------------------------------
  
  !---------------------------------------------------------------------
  SUBROUTINE new_ocean_statistics()
  END SUBROUTINE new_ocean_statistics
  !---------------------------------------------------------------------
  
END MODULE mo_ocean_statistics
