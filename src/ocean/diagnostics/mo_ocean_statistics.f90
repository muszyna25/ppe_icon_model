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
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
!!
MODULE mo_ocean_statistics
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: max_char_length
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d, t_subset_range, t_patch_vert
  USE mo_grid_config,            ONLY: n_dom
  USE mo_grid_subset,            ONLY: get_index_range
  USE mo_sync,                   ONLY: sync_patch_array, sync_e, sync_c !, sync_v
  USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer, &
    & diagnostics_level, &
    & eos_type, i_sea_ice, l_staggered_timestep, gibraltar
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_io_config,              ONLY: n_checkpoints
  USE mo_run_config,             ONLY: nsteps, dtime, ltimer, output_mode
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ext_data_types,         ONLY: t_external_data
  !USE mo_io_units,               ONLY: filename_max
  USE mo_datetime,               ONLY: t_datetime, print_datetime, add_time, datetime_to_string
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_total, timer_solve_ab,  &
    & timer_tracer_ab, timer_vert_veloc, timer_normal_veloc, &
    & timer_upd_phys, timer_upd_flx
!  USE mo_oce_ab_timestepping,    ONLY: solve_free_surface_eq_ab, &
!    & calc_normal_velocity_ab,  &
!    & calc_vert_velocity,       &
!    & update_time_indices
  USE mo_oce_types,              ONLY: t_hydro_ocean_state, t_hydro_ocean_acc, t_hydro_ocean_diag, &
    & t_hydro_ocean_prog
  USE mo_oce_state,              ONLY: ocean_restart_list
 ! USE mo_ocean_initialization,   ONLY: set_lateral_boundary_values
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_scalar_product,         ONLY: calc_scalar_product_veloc_3d
  USE mo_oce_tracer,             ONLY: advect_tracer_ab
  USE mo_io_restart,             ONLY: write_restart_info_file, create_restart_file
  USE mo_sea_ice,                ONLY: update_ice_statistic, compute_mean_ice_statistics, reset_ice_statistics
  USE mo_sea_ice_types,          ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
    & t_sea_ice
  USE mo_oce_physics,            ONLY: t_ho_params, update_ho_params
  USE mo_oce_thermodyn,          ONLY: calc_density_mpiom_func, calc_density_lin_eos_func,&
    & calc_density_jmdwfg06_eos_func, calc_potential_density, &
    & calc_density
  USE mo_name_list_output,       ONLY: write_name_list_output, istime4name_list_output
  USE mo_oce_diagnostics,        ONLY: calc_slow_oce_diagnostics, calc_fast_oce_diagnostics, &
    & construct_oce_diagnostics,&
    & destruct_oce_diagnostics, t_oce_timeseries, &
    & calc_moc, calc_psi
  USE mo_oce_ab_timestepping_mimetic, ONLY: init_ho_lhs_fields_mimetic
  USE mo_linked_list,            ONLY: t_list_element, find_list_element
  USE mo_var_list,               ONLY: print_var_list
  USE mo_io_restart_attributes,  ONLY: get_restart_attribute
  USE mo_mpi,                    ONLY: my_process_is_stdio
  USE mo_time_config,            ONLY: time_config
  USE mo_master_control,         ONLY: is_restart_run
  USE mo_statistics
  USE mo_sea_ice_nml,            ONLY: i_ice_dyn
  USE mo_util_dbg_prnt,          ONLY: dbg_print
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: update_ocean_statistics
  PUBLIC :: compute_mean_ocean_statistics
  PUBLIC :: reset_ocean_statistics
  PUBLIC :: new_ocean_statistics

  
  !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  

  !-------------------------------------------------------------------------
  
CONTAINS

    
  !---------------------------------------------------------------------
  SUBROUTINE update_ocean_statistics(ocean_state,p_sfc_flx,cells,edges,verts,max_zlev)
    TYPE(t_hydro_ocean_state), INTENT(inout) :: ocean_state
    TYPE(t_sfc_flx),           INTENT(inout) :: p_sfc_flx
    TYPE(t_subset_range),      INTENT(in)    :: cells,edges,verts
    INTEGER, INTENT(in)                      :: max_zlev
    
    INTEGER :: jtrc,i
    
    
    ! update ocean state accumulated values
    CALL add_fields(ocean_state%p_acc%h     , ocean_state%p_prog(nnew(1))%h, cells)
    CALL add_fields(ocean_state%p_acc%u     , ocean_state%p_diag%u         , cells)
    CALL add_fields(ocean_state%p_acc%v     , ocean_state%p_diag%v         , cells)
    CALL add_fields(ocean_state%p_acc%rhopot, ocean_state%p_diag%rhopot    , cells)
    DO jtrc=1,no_tracer
      CALL add_fields(ocean_state%p_acc%tracer(:,:,:,jtrc),           &
        & ocean_state%p_prog(nnew(1))%tracer(:,:,:,jtrc), &
        & cells)
    END DO
    CALL add_fields(ocean_state%p_acc%u_vint        , ocean_state%p_diag%u_vint        , cells)
    CALL add_fields(ocean_state%p_acc%w             , ocean_state%p_diag%w             , cells,max_zlev+1)
    CALL add_fields(ocean_state%p_acc%div_mass_flx_c, ocean_state%p_diag%div_mass_flx_c, cells)
    CALL add_fields(ocean_state%p_acc%rho           , ocean_state%p_diag%rho           , cells)
    CALL add_fields(ocean_state%p_acc%mass_flx_e    , ocean_state%p_diag%mass_flx_e    , edges)
    CALL add_fields(ocean_state%p_acc%vort          , ocean_state%p_diag%vort          , verts,max_zlev)
    CALL add_fields(ocean_state%p_acc%kin           , ocean_state%p_diag%kin           , cells)
    
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
    CALL add_fields(p_sfc_flx%forc_hfrelax_acc                , p_sfc_flx%forc_hfrelax                , cells)
    CALL add_fields(p_sfc_flx%forc_fwrelax_acc                , p_sfc_flx%forc_fwrelax                , cells)
    CALL add_fields(p_sfc_flx%data_surfRelax_Temp_acc(:,:)    , p_sfc_flx%data_surfRelax_Temp(:,:)    , cells)
    CALL add_fields(p_sfc_flx%data_surfRelax_Salt_acc(:,:)    , p_sfc_flx%data_surfRelax_Salt(:,:)    , cells)
    CALL add_fields(p_sfc_flx%topBoundCond_Salt_vdiff_acc(:,:), p_sfc_flx%topBoundCond_Temp_vdiff(:,:), cells)
    CALL add_fields(p_sfc_flx%topBoundCond_Salt_vdiff_acc(:,:), p_sfc_flx%topBoundCond_Temp_vdiff(:,:), cells)
    
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
    
    
    p_acc%tracer                    = p_acc%tracer                   /REAL(nsteps_since_last_output,wp)
    p_acc%h                         = p_acc%h                        /REAL(nsteps_since_last_output,wp)
    p_acc%u                         = p_acc%u                        /REAL(nsteps_since_last_output,wp)
    p_acc%v                         = p_acc%v                        /REAL(nsteps_since_last_output,wp)
    p_acc%rhopot                    = p_acc%rhopot                   /REAL(nsteps_since_last_output,wp)
    p_acc%u_vint                    = p_acc%u_vint                   /REAL(nsteps_since_last_output,wp)
    p_acc%w                         = p_acc%w                        /REAL(nsteps_since_last_output,wp)
    p_acc%div_mass_flx_c            = p_acc%div_mass_flx_c           /REAL(nsteps_since_last_output,wp)
    p_acc%rho                       = p_acc%rho                      /REAL(nsteps_since_last_output,wp)
    p_acc%mass_flx_e                = p_acc%mass_flx_e               /REAL(nsteps_since_last_output,wp)
    p_acc%vort                      = p_acc%vort                     /REAL(nsteps_since_last_output,wp)
    p_acc%kin                       = p_acc%kin                      /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%topBoundCond_windStress_u_acc = p_sfc_flx%topBoundCond_windStress_u_acc/REAL(nsteps_since_last_output,wp)
    p_sfc_flx%topBoundCond_windStress_v_acc = p_sfc_flx%topBoundCond_windStress_v_acc/REAL(nsteps_since_last_output,wp)
    p_sfc_flx%HeatFlux_ShortWave_acc        = p_sfc_flx%HeatFlux_ShortWave           /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%HeatFlux_LongWave_acc         = p_sfc_flx%HeatFlux_LongWave            /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%HeatFlux_Sensible_acc         = p_sfc_flx%HeatFlux_Sensible            /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%HeatFlux_Latent_acc           = p_sfc_flx%HeatFlux_Latent              /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%HeatFlux_Total_acc            = p_sfc_flx%HeatFlux_Total               /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%FrshFlux_Precipitation_acc    = p_sfc_flx%FrshFlux_Precipitation_acc   /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%FrshFlux_SnowFall_acc         = p_sfc_flx%FrshFlux_SnowFall_acc        /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%FrshFlux_Evaporation_acc      = p_sfc_flx%FrshFlux_Evaporation_acc     /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%FrshFlux_Runoff_acc           = p_sfc_flx%FrshFlux_Runoff_acc          /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%FrshFlux_TotalSalt_acc        = p_sfc_flx%FrshFlux_TotalSalt_acc       /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%FrshFlux_TotalOcean_acc       = p_sfc_flx%FrshFlux_TotalOcean_acc      /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%FrshFlux_TotalIce_acc         = p_sfc_flx%FrshFlux_TotalIce_acc        /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_fwrelax_acc              = p_sfc_flx%forc_fwrelax_acc             /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%FrshFlux_VolumeIce_acc        = p_sfc_flx%FrshFlux_VolumeIce_acc       /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%FrshFlux_VolumeTotal_acc      = p_sfc_flx%FrshFlux_VolumeTotal_acc     /REAL(nsteps_since_last_output,wp)
    p_sfc_flx%forc_hfrelax_acc              = p_sfc_flx%forc_hfrelax_acc             /REAL(nsteps_since_last_output,wp)
    IF(no_tracer>0)THEN
      p_sfc_flx%data_surfRelax_Temp_acc     = p_sfc_flx%data_surfRelax_Temp_acc        /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%topBoundCond_Temp_vdiff_acc = p_sfc_flx%topBoundCond_Temp_vdiff_acc    /REAL(nsteps_since_last_output,wp)
      p_sfc_flx%topBoundCond_Salt_vdiff_acc = p_sfc_flx%topBoundCond_Salt_vdiff_acc    /REAL(nsteps_since_last_output,wp)
    ENDIF  
  END SUBROUTINE compute_mean_ocean_statistics
  !---------------------------------------------------------------------
  
  !---------------------------------------------------------------------
  SUBROUTINE reset_ocean_statistics(p_acc,p_sfc_flx,nsteps_since_last_output)
    TYPE(t_hydro_ocean_acc), INTENT(inout) :: p_acc
    TYPE(t_sfc_flx),         INTENT(inout) :: p_sfc_flx
    INTEGER,OPTIONAL,        INTENT(inout) :: nsteps_since_last_output
    
    IF (PRESENT(nsteps_since_last_output)) nsteps_since_last_output        = 0
    
    p_acc%tracer                    = 0.0_wp
    p_acc%h                         = 0.0_wp
    p_acc%u                         = 0.0_wp
    p_acc%v                         = 0.0_wp
    p_acc%rhopot                    = 0.0_wp
    p_acc%u_vint                    = 0.0_wp
    p_acc%w                         = 0.0_wp
    p_acc%div_mass_flx_c            = 0.0_wp
    p_acc%rho                       = 0.0_wp
    p_acc%mass_flx_e                = 0.0_wp
    p_acc%vort                      = 0.0_wp
    p_acc%kin                       = 0.0_wp
    p_sfc_flx%topBoundCond_windStress_u_acc = 0.0_wp
    p_sfc_flx%topBoundCond_windStress_v_acc = 0.0_wp
    IF (no_tracer>0) THEN
      p_sfc_flx%HeatFlux_ShortWave_acc          = 0.0_wp
      p_sfc_flx%HeatFlux_LongWave_acc           = 0.0_wp
      p_sfc_flx%HeatFlux_Sensible_acc           = 0.0_wp
      p_sfc_flx%HeatFlux_Latent_acc             = 0.0_wp
      p_sfc_flx%HeatFlux_Total_acc              = 0.0_wp
      p_sfc_flx%topBoundCond_Temp_vdiff_acc     = 0.0_wp
      p_sfc_flx%forc_hfrelax_acc                = 0.0_wp
      IF (no_tracer>1) THEN
        p_sfc_flx%FrshFlux_Precipitation_acc    = 0.0_wp
        p_sfc_flx%FrshFlux_SnowFall_acc         = 0.0_wp
        p_sfc_flx%FrshFlux_Evaporation_acc      = 0.0_wp
        p_sfc_flx%FrshFlux_Runoff_acc           = 0.0_wp
        p_sfc_flx%FrshFlux_TotalSalt_acc        = 0.0_wp
        p_sfc_flx%FrshFlux_TotalOcean_acc       = 0.0_wp
        p_sfc_flx%FrshFlux_TotalIce_acc         = 0.0_wp
        p_sfc_flx%forc_fwrelax_acc              = 0.0_wp
        p_sfc_flx%FrshFlux_VolumeIce_acc        = 0.0_wp
        p_sfc_flx%FrshFlux_VolumeTotal_acc      = 0.0_wp
        p_sfc_flx%topBoundCond_Salt_vdiff_acc   = 0.0_wp
        p_sfc_flx%data_surfRelax_Temp_acc = 0.0_wp
        p_sfc_flx%data_surfRelax_Salt_acc = 0.0_wp
      ENDIF
    ENDIF
  END SUBROUTINE reset_ocean_statistics
  !---------------------------------------------------------------------
  
  !---------------------------------------------------------------------
  SUBROUTINE new_ocean_statistics()
  END SUBROUTINE new_ocean_statistics
  !---------------------------------------------------------------------
  
END MODULE mo_ocean_statistics
