!>
!! Provide an implementation of the ocean forcing.
!!
!! Provide an implementation of the parameters used for surface forcing
!! of the hydrostatic ocean model.
!!
!! @author Peter Korn, MPI
!! @author Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
!!  Modification by Stephan Lorenz, MPI-M (2010-06):
!!   - renaming and adjustment to ocean domain and patch_oce
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
MODULE mo_ocean_forcing
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: filename_max
  USE mo_grid_config,         ONLY: nroot
  USE mo_parallel_config,     ONLY: nproma
  USE mo_ocean_nml,           ONLY: no_tracer,                                &
    & basin_height_deg, basin_width_deg, basin_center_lat, basin_center_lon,  &
    & forcing_windstress_zonal_waveno, forcing_windstress_merid_waveno,       &
    & init_oce_relax, type_3dimRelax_Salt, type_3dimRelax_Temp,               &
    & type_surfRelax_Salt, type_surfRelax_Temp,                               &
    & forcing_windStress_u_amplitude, forcing_windStress_v_amplitude,         &
    & forcing_windstress_u_type, forcing_windstress_v_type,                   &
    & forcing_windspeed_type, forcing_windspeed_amplitude,                    &
    & relax_temperature_min, relax_temperature_max, initial_temperature_type, &
	& initial_temperature_bottom,initial_temperature_north,initial_temperature_south,&
#ifdef __SX__                                                                
    & forcing_windstress_zonalWavePhas,                                       &
#else                                                                        
    & forcing_windstress_zonalWavePhase,                                      &
#endif
    & forcing_temperature_poleLat                                         
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_math_constants,      ONLY: pi, deg2rad, pi_2
  USE mo_impl_constants,      ONLY: max_char_length, sea_boundary, success
  USE mo_sea_ice_types,       ONLY: t_sfc_flx, t_atmos_fluxes
  USE mo_ocean_state,           ONLY: set_oce_tracer_info
  USE mo_ocean_types,           ONLY: t_hydro_ocean_state
  USE mo_dynamics_config,     ONLY: nold
!   USE mo_ocean_initial_conditions, ONLY: SST_LinearMeridional, increaseTracerLevelsLinearly

  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_var_list,            ONLY: add_var, add_ref
  USE mo_var_metadata,        ONLY: groups
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_fortran_tools,       ONLY: assign_if_present
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_read_interface,     ONLY: openInputFile, closeFile, t_stream_id, &
    &                              onCells, read_2D_1lev_1time

  IMPLICIT NONE
  PRIVATE
  INCLUDE 'netcdf.inc'

  CHARACTER(LEN=12)           :: str_module    = 'oceForcing  '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug

  ! Public interface
  PUBLIC :: construct_ocean_forcing, destruct_ocean_forcing
  PUBLIC :: init_ocean_forcing

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Constructor of surface fluxes for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
!<Optimize:inUse>
  SUBROUTINE construct_ocean_forcing(p_patch, p_sfc_flx, var_list)
    !
    TYPE(t_patch),   INTENT(IN)    :: p_patch
    TYPE(t_sfc_flx), INTENT(INOUT) :: p_sfc_flx
    TYPE(t_var_list),INTENT(INOUT) :: var_list

    ! Local variables
    INTEGER                        :: alloc_cell_blocks, ist
    ! CHARACTER(len=max_char_length) :: oce_tracer_names(no_tracer),&
    ! &                                 oce_tracer_units(no_tracer),&
    ! &                                 oce_tracer_longnames(no_tracer)
    ! INTEGER                        :: oce_tracer_codes(no_tracer)

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_ocean_forcing'

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    alloc_cell_blocks = p_patch%alloc_cell_blocks

    CALL add_var(var_list, 'topBoundCond_windStress_u', p_sfc_flx%topBoundCond_windStress_u , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('topBoundCond_windStress_u', 'Pa', 'topBoundCond_windStress_u', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'topBoundCond_windStress_v', p_sfc_flx%topBoundCond_windStress_v , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('topBoundCond_windStress_v', 'Pa', 'topBoundCond_windStress_v', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'HeatFlux_Total', p_sfc_flx%HeatFlux_Total , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Total', 'W/m2', 'HeatFlux_Total', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'FrshFlux_VolumeTotal', p_sfc_flx%FrshFlux_VolumeTotal , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_VolumeTotal', 'm/s', 'FrshFlux_VolumeTotal', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'HeatFlux_ShortWave', p_sfc_flx%HeatFlux_ShortWave , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_ShortWave', 'W/m2', 'HeatFlux_ShortWave', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'HeatFlux_LongWave', p_sfc_flx%HeatFlux_LongWave , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_LongWave', 'W/m2', 'HeatFlux_LongWave', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'HeatFlux_Sensible', p_sfc_flx%HeatFlux_Sensible , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Sensible', 'W/m2', 'HeatFlux_Sensible', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'HeatFlux_Latent', p_sfc_flx%HeatFlux_Latent , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Latent', 'W/m2', 'HeatFlux_Latent', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'FrshFlux_Precipitation', p_sfc_flx%FrshFlux_Precipitation , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Precipitation', 'm/s', 'FrshFlux_Precipitation', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'FrshFlux_SnowFall', p_sfc_flx%FrshFlux_SnowFall , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_SnowFall', 'm/s', 'FrshFlux_SnowFall', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'FrshFlux_Evaporation', p_sfc_flx%FrshFlux_Evaporation , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Evaporation', 'm/s', 'FrshFlux_Evaporation', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'FrshFlux_Runoff', p_sfc_flx%FrshFlux_Runoff , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Runoff', 'm/s', 'FrshFlux_Runoff', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'FrshFlux_TotalSalt', p_sfc_flx%FrshFlux_TotalSalt , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_TotalSalt', 'm/s', 'FrshFlux_TotalSalt', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'FrshFlux_TotalOcean', p_sfc_flx%FrshFlux_TotalOcean , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_TotalOcean', 'm/s', 'FrshFlux_TotalOcean', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

    CALL add_var(var_list, 'FrshFlux_TotalIce', p_sfc_flx%FrshFlux_TotalIce , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_TotalIce', 'm/s', 'FrshFlux_TotalIce', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'FrshFlux_VolumeIce', p_sfc_flx%FrshFlux_VolumeIce, &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_VolumeIce', 'm/s', 'FrshFlux_VolumeIce', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'HeatFlux_Relax', p_sfc_flx%HeatFlux_Relax , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Relax', 'W/m2', 'HeatFlux_Relax', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'FrshFlux_Relax', p_sfc_flx%FrshFlux_Relax , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Relax', 'm/s', 'FrshFlux_Relax', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'TempFlux_Relax', p_sfc_flx%TempFlux_Relax , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('TempFlux_Relax', 'K/s', 'TempFlux_Relax', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'SaltFlux_Relax', p_sfc_flx%SaltFlux_Relax , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('SaltFlux_Relax', 'psu/s', 'SaltFlux_Relax', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    IF(no_tracer>=1) THEN
      ! there are four tracer related fields: tracer focing, tracer relaxation
      ! and both accumulated
      CALL add_var(var_list, 'topBoundCond_Temp_vdiff', p_sfc_flx%topBoundCond_Temp_vdiff, &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('topBoundCond_Temp_vdiff', 'K*m/s', 'topBoundCond_Temp_vdiff', DATATYPE_FLT32),&
        &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      CALL add_var(var_list, 'topBoundCond_Salt_vdiff', p_sfc_flx%topBoundCond_Salt_vdiff, &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('topBoundCond_Salt_vdiff', 'psu*m/s', 'topBoundCond_Salt_vdiff', DATATYPE_FLT32),&
        &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      CALL add_var(var_list, 'data_surfRelax_Temp', p_sfc_flx%data_surfRelax_Temp, &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('data_surfRelax_Temp', 'C', 'data_surfRelax_Temp', DATATYPE_FLT32),&
        &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lrestart=.FALSE., loutput=.TRUE.)		
!         &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.TRUE.)

      CALL add_var(var_list, 'data_surfRelax_Salt', p_sfc_flx%data_surfRelax_Salt, &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('data_surfRelax_Salt', 'psu', 'data_surfRelax_Salt', DATATYPE_FLT32),&
        &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      CALL add_var(var_list, 'topBoundCond_Temp_vdiff_acc', p_sfc_flx%topBoundCond_Temp_vdiff_acc , &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('topBoundCond_Temp_vdiff_acc', 'K*m/s', 'topBoundCond_Temp_vdiff_acc', DATATYPE_FLT32),&
        &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      CALL add_var(var_list, 'topBoundCond_Salt_vdiff_acc', p_sfc_flx%topBoundCond_Salt_vdiff_acc , &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('topBoundCond_Salt_vdiff_acc', 'psu*m/s', 'topBoundCond_Salt_vdiff_acc', DATATYPE_FLT32),&
        &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      CALL add_var(var_list, 'data_surfRelax_Temp_acc', p_sfc_flx%data_surfRelax_Temp_acc , &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('data_surfRelax_Temp_acc', 'C', 'data_surfRelax_Temp_acc', DATATYPE_FLT32),&
        &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      CALL add_var(var_list, 'data_surfRelax_Salt_acc', p_sfc_flx%data_surfRelax_Salt_acc , &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('data_surfRelax_Salt_acc', 'psu', 'data_surfRelax_Salt_acc', DATATYPE_FLT32),&
        &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
    ENDIF

    CALL add_var(var_list, 'topBoundCond_windStress_u_acc', p_sfc_flx%topBoundCond_windStress_u_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('topBoundCond_windStress_u_acc', 'Pa', 'topBoundCond_windStress_u_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))

    CALL add_var(var_list, 'topBoundCond_windStress_v_acc', p_sfc_flx%topBoundCond_windStress_v_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('topBoundCond_windStress_v_acc', 'Pa', 'topBoundCond_windStress_v_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))

    CALL add_var(var_list, 'HeatFlux_Total_acc', p_sfc_flx%HeatFlux_Total_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Total_acc', 'W/m2', 'HeatFlux_Total_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))

    CALL add_var(var_list, 'FrshFlux_VolumeTotal_acc', p_sfc_flx%FrshFlux_VolumeTotal_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_VolumeTotal_acc', 'm/s', 'FrshFlux_VolumeTotal_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))

    CALL add_var(var_list, 'HeatFlux_ShortWave_acc', p_sfc_flx%HeatFlux_ShortWave_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_ShortWave_acc', 'W/m2', 'HeatFlux_ShortWave_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))

    CALL add_var(var_list, 'HeatFlux_LongWave_acc', p_sfc_flx%HeatFlux_LongWave_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_LongWave_acc', 'W/m2', 'HeatFlux_LongWave_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))

    CALL add_var(var_list, 'HeatFlux_Sensible_acc', p_sfc_flx%HeatFlux_Sensible_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Sensible_acc', 'W/m2', 'HeatFlux_Sensible_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))

    CALL add_var(var_list, 'HeatFlux_Latent_acc', p_sfc_flx%HeatFlux_Latent_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Latent_acc', 'W/m2', 'HeatFlux_Latent_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))

    CALL add_var(var_list, 'FrshFlux_Precipitation_acc', p_sfc_flx%FrshFlux_Precipitation_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Precipitation_acc', 'm/s', 'FrshFlux_Precipitation_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))

    CALL add_var(var_list, 'FrshFlux_SnowFall_acc', p_sfc_flx%FrshFlux_SnowFall_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_SnowFall_acc', 'm/s', 'FrshFlux_SnowFall_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))

    CALL add_var(var_list, 'FrshFlux_Evaporation_acc', p_sfc_flx%FrshFlux_Evaporation_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Evaporation_acc', 'm/s', 'FrshFlux_Evaporation_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))

    CALL add_var(var_list, 'FrshFlux_Runoff_acc', p_sfc_flx%FrshFlux_Runoff_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Runoff_acc', 'm/s', 'FrshFlux_Runoff_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))

    CALL add_var(var_list, 'FrshFlux_TotalSalt_acc', p_sfc_flx%FrshFlux_TotalSalt_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_TotalSalt_acc', 'm/s', 'FrshFlux_TotalSalt_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))

    CALL add_var(var_list, 'FrshFlux_TotalOcean_acc', p_sfc_flx%FrshFlux_TotalOcean_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_TotalOcean_acc', 'm/s', 'FrshFlux_TotalOcean_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))

    CALL add_var(var_list, 'FrshFlux_TotalIce_acc', p_sfc_flx%FrshFlux_TotalIce_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_TotalIce_acc', 'm/s', 'FrshFlux_TotalIce_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'FrshFlux_VolumeIce_acc', p_sfc_flx%FrshFlux_VolumeIce_acc, &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_VolumeIce_acc', 'm/s', 'FrshFlux_VolumeIce_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'HeatFlux_Relax_acc', p_sfc_flx%HeatFlux_Relax_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Relax_acc', 'W/m2', 'HeatFlux_Relax_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'FrshFlux_Relax_acc', p_sfc_flx%FrshFlux_Relax_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Relax_acc', 'm/s', 'FrshFlux_Relax_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'TempFlux_Relax_acc', p_sfc_flx%TempFlux_Relax_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('TempFlux_Relax_acc', 'K/s', 'TempFlux_Relax_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'SaltFlux_Relax_acc', p_sfc_flx%SaltFlux_Relax_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('SaltFlux_Relax_acc', 'psu/s', 'SaltFlux_Relax_acc', DATATYPE_FLT32),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))

    ! cartesians
    ALLOCATE(p_sfc_flx%topBoundCond_windStress_cc(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for forcing wind stress cc  failed')
    END IF

    ! init of cartesian coordinates:
    p_sfc_flx%topBoundCond_windStress_cc(:,:)%x(1) = 0.0_wp
    p_sfc_flx%topBoundCond_windStress_cc(:,:)%x(2) = 0.0_wp
    p_sfc_flx%topBoundCond_windStress_cc(:,:)%x(3) = 0.0_wp

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE construct_ocean_forcing
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !>
  !! Destructor surface flux forcing for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
!<Optimize:inUse>
  SUBROUTINE destruct_ocean_forcing(p_sfc_flx)
    TYPE(t_sfc_flx), INTENT(INOUT) :: p_sfc_flx
    !
    ! Local variables

    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_ocean_forcing'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    ! forcing fields are handled by the ocean_default_list

    DEALLOCATE(p_sfc_flx%topBoundCond_windStress_cc, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for forcing wind stress cc failed')
    END IF
    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE destruct_ocean_forcing
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE init_ocean_forcing(patch_2d, patch_3d, ocean_state, atmos_fluxes, fu10)
    !
    TYPE(t_patch),TARGET, INTENT(in)        :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_atmos_fluxes)                    :: atmos_fluxes
    REAL(wp), INTENT(INOUT)                 :: fu10(:,:)      !  windspeed: p_as%fu10

    TYPE(t_subset_range), POINTER :: all_cells, owned_cells

    
    CALL init_ocean_WindForcing(patch_3d, ocean_state, atmos_fluxes, fu10)

    IF (init_oce_relax > 0) THEN
      CALL init_ho_relaxation(patch_3d, ocean_state, atmos_fluxes)
    END IF
  END SUBROUTINE init_ocean_forcing
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Initialization of temperature and salinity relaxation for the hydrostatic ocean model.
  !! Temperature and salinity relaxation data are read from external data
  !
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M, 2011-11
  !
  !-------------------------------------------------------------------------
  !
!<Optimize:inUse>
  SUBROUTINE init_ho_relaxation(patch_3d, ocean_state, atmos_fluxes)

    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_atmos_fluxes)                    :: atmos_fluxes

    ! Local Variables

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_initial_conditions:init_ho_relaxation'
    CHARACTER(filename_max) :: relax_init_file   !< file name for reading in

    LOGICAL :: l_exist
    INTEGER :: i_lev, no_cells, no_levels, jb, jc, jk
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: ncid, dimid
    TYPE(t_stream_id) :: stream_id

    REAL(wp):: z_c(nproma,1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp):: z_surfRelax(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: temperature_difference, poleLat, waveNo
	
    REAL(wp):: distan
    REAL(wp):: perturbation_lat, perturbation_lon,  max_perturbation, perturbation_width
    REAL(wp) :: basin_NorthBoundary, basin_SouthBoundary, lat_diff
 	

    TYPE(t_patch),POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    !-------------------------------------------------------------------------
    patch_2d  => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    basin_NorthBoundary    = (basin_center_lat + 0.5_wp*basin_height_deg) * deg2rad
    basin_SouthBoundary    = (basin_center_lat - 0.5_wp*basin_height_deg) * deg2rad

    ! Read relaxation data from file
    IF (init_oce_relax == 1) THEN

       ! sphere_radius = grid_sphere_radius
      ! u0 =(2.0_wp*pi*sphere_radius)/(12.0_wp*24.0_wp*3600.0_wp)

      CALL message (TRIM(routine), 'start')

      i_lev        = patch_2d%level

      ! Relaxation variables are read from relax_init_file
!       WRITE (relax_init_file,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',i_lev, '-relax.nc'
      relax_init_file='ocean-relax.nc'
      
      IF (my_process_is_stdio()) THEN

        INQUIRE (FILE=relax_init_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          WRITE(message_text,'(3a)') 'netcdf file named ', TRIM(relax_init_file),' not found!'
          CALL message(TRIM(routine),TRIM(message_text))
          CALL finish(TRIM(routine),'netcdf file for reading T/S relax. input not found - ABORT')
        ENDIF

        WRITE(message_text,'(3a)') 'netcdf file named ', TRIM(relax_init_file), &
          & ' opened for reading'
        CALL message(TRIM(routine),TRIM(message_text))

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(relax_init_file), nf_nowrite, ncid))

        !
        ! get number of cells
        !
        CALL nf(nf_inq_dimid(ncid, 'ncells', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))

        !
        ! check the number of cells
        !
        WRITE(message_text,'(a,i6)') 'No of cells =', no_cells
        CALL message(TRIM(routine),TRIM(message_text))
        IF (patch_2d%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(routine),&
            & 'Number of patch cells and cells in T/S relaxation input file do not match - ABORT')
        ENDIF
        !
        ! get number of levels
        !
        CALL nf(nf_inq_dimid(ncid, 'level', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, no_levels))

        !
        ! check the number of cells
        !
        WRITE(message_text,'(a,i6)') 'No of vertical levels =', no_levels
        CALL message(TRIM(routine),TRIM(message_text))
        IF (no_levels /= 1) THEN
          CALL finish(TRIM(routine),'Number of vertical levels is not equal 1 - ABORT')
        ENDIF

        CALL nf(nf_close(ncid))

      ENDIF  !  stdio

      !-------------------------------------------------------
      !
      ! Read ocean relaxation data at cells
      !
      !-------------------------------------------------------

      stream_id = openInputFile(relax_init_file, patch_2d)

      ! triangle center and edges

      ! read temperature
      !  - read one data set, annual mean only
      !  - "T": annual mean temperature
      CALL read_2D_1lev_1time(stream_id, onCells, 'T', z_surfRelax)

      IF (no_tracer>=1) THEN
        atmos_fluxes%data_surfRelax_Temp(:,:) = z_surfRelax(:,:)
      ELSE
        CALL message( TRIM(routine),'WARNING: no tracer used, but init relaxation attempted')
      END IF

      ! read salinity
      !  - "S": annual mean salinity
      IF (no_tracer > 1) THEN
        CALL read_2D_1lev_1time(stream_id, onCells, 'S', z_surfRelax)
        atmos_fluxes%data_surfRelax_Salt(:,:) = z_surfRelax(:,:)
      END IF

      ! close file
      CALL closeFile(stream_id)

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          IF ( patch_3d%lsm_c(jc,1,jb) > sea_boundary ) THEN
            atmos_fluxes%data_surfRelax_Temp(jc,jb) = 0.0_wp
            IF (no_tracer>1) atmos_fluxes%data_surfRelax_Salt(jc,jb) = 0.0_wp
          ENDIF
        END DO
      END DO

      CALL message( TRIM(routine),'Ocean T/S relaxation reading finished' )

    ENDIF  !  read relaxation data from file

    !-------------------------------------------------------
    !
    ! use initialized temperature/salinity, assigned to tracer, for 2-dim/3-dim relaxation
    !  - relaxation switch equals 3
    !
    !-------------------------------------------------------

    SELECT CASE(type_surfRelax_Temp)
    CASE(3)
      atmos_fluxes%data_surfRelax_Temp(:,:) = ocean_state%p_prog(nold(1))%tracer(:,1,:,1)
      
      

    CASE(4)
      ! smooth ape relaxation, as in temperature_smoothAPE in mo_cean_initial_conditions
      temperature_difference = relax_temperature_max - relax_temperature_min
      poleLat = ABS(forcing_temperature_poleLat * deg2rad)
      waveNo = pi_2 / poleLat

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          DO jk=1, MIN(1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb))

            atmos_fluxes%data_surfRelax_Temp(jc,jb) = relax_temperature_min    + &
              & (COS(waveNo * MIN(ABS(patch_2d%cells%center(jc,jb)%lat), poleLat))**2) * temperature_difference

          END DO
        END DO
      END DO
      CASE(5)
      !Only valid for this specific testcase: temperature is resored to initial tempertaure without perturbation
      IF(initial_temperature_type==215.OR.initial_temperature_type==214)THEN

        atmos_fluxes%data_surfRelax_Temp(:,:) = 0.0_wp!initial_temperature_south

       temperature_difference = initial_temperature_north - initial_temperature_south
       lat_diff               = basin_NorthBoundary - basin_SouthBoundary  !  basin_height_deg*deg2rad

        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index
            jk=1
            atmos_fluxes%data_surfRelax_Temp(jc,jb) = initial_temperature_north &
            &- temperature_difference*((basin_NorthBoundary-patch_2d%cells%center(jc,jb)%lat)/lat_diff)
  
            atmos_fluxes%data_surfRelax_Temp(jc,jb) = MERGE(ocean_state%p_prog(nold(1))%tracer(jc,1,jb,1), &
            &initial_temperature_north, patch_2d%cells%center(jc,jb)%lat>basin_SouthBoundary)
            
			
            atmos_fluxes%data_surfRelax_Temp(jc,jb) = MERGE(&
            &initial_temperature_north-1.0_wp, &			
			      &ocean_state%p_prog(nold(1))%tracer(jc,1,jb,1), &
			      &patch_2d%cells%center(jc,jb)%lat>(basin_center_lat + 1.0_wp*basin_height_deg) * deg2rad)!&
			!&.AND.lat(jc,jb)<(basin_center_lat + 1.25_wp*basin_height_deg) * deg2rad)                    
          END DO
        END DO

        !CALL SST_LinearMeridional(patch_3d, ocean_state%p_prog(nold(1))%tracer(:,:,:,1))
        !atmos_fluxes%data_surfRelax_Temp(:,:)=ocean_state%p_prog(nold(1))%tracer(:,1,:,1)
      !ENDIF
      ELSEIF(initial_temperature_type==203)THEN
        perturbation_lat = basin_center_lat + 0.1_wp * basin_height_deg
        perturbation_lon = basin_center_lon + 0.1_wp * basin_width_deg
        max_perturbation  = 0.1_wp!20.1_wp
        perturbation_width  = 10.0_wp!1.5_wp

        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index
            distan = SQRT((patch_2d%cells%center(jc,jb)%lat - perturbation_lat * deg2rad)**2 + &
            & (patch_2d%cells%center(jc,jb)%lon - perturbation_lon * deg2rad)**2)
  
            atmos_fluxes%data_surfRelax_Temp(jc,jb)=ocean_state%p_prog(nold(1))%tracer(jc,1,jb,1)
            !Local cold perturbation
            IF(distan<=5.0_wp*deg2rad)THEN
              atmos_fluxes%data_surfRelax_Temp(jc,jb) =          &
              & ocean_state%p_prog(nold(1))%tracer(jc,1,jb,1)        &
              & - 2.0_wp*max_perturbation*EXP(-(distan/(perturbation_width*deg2rad))**2) &
              !                &   * sin(pi*v_base%zlev_m(jk)/4000.0_wp)!&
              & * SIN(pi*patch_3d%p_patch_1d(1)%zlev_m(1) / patch_3d%p_patch_1d(1)%zlev_i(2))
             ENDIF !Local cold perturbation
           END DO
         END DO
      ENDIF
  
    END SELECT

    IF (type_surfRelax_Salt == 3) THEN
      IF (no_tracer > 1) THEN
        atmos_fluxes%data_surfRelax_Salt(:,:) = ocean_state%p_prog(nold(1))%tracer(:,1,:,2)
      END IF
    END IF

    IF (type_3dimRelax_Temp >= 3) THEN
      ocean_state%p_aux%data_3dimRelax_Temp(:,:,:) = ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
      IF (type_3dimRelax_Temp ==4) THEN
        CALL init_3Drelax_coefficient_NS( &
          & patch_3D = patch_3d, &
          & relax_coefficient=ocean_state%p_aux%relax_3dim_coefficient, &
          & SouthBoundary=basin_SouthBoundary, &
          & NorthBoundary=basin_NorthBoundary, &
          & relaxWidth=1.5_wp * deg2rad)
      END IF
    END IF
    IF (type_3dimRelax_Salt == 3) THEN
      IF (no_tracer > 1) THEN
        ocean_state%p_aux%data_3dimRelax_Salt(:,:,:) = ocean_state%p_prog(nold(1))%tracer(:,:,:,2)
      ELSE
        CALL finish(TRIM(routine),' type_3dimRelax_Salt=3 and no_tracer<2 - ABORT')
      END IF
    END IF

    !---------Debug Diagnostics-------------------------------------------
    IF (type_surfRelax_Temp > 0) THEN
      idt_src=0  ! output print level - 0: print in any case
      ! z_c(:,1,:) = atmos_fluxes%data_surfRelax_Temp(:,:)
      CALL dbg_print('init relaxation - T'       ,atmos_fluxes%data_surfRelax_Temp(:,:), &
        & str_module,idt_src, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    END IF
    IF (type_surfRelax_Salt > 0) THEN
      IF (no_tracer > 1) THEN
        z_c(:,1,:) = atmos_fluxes%data_surfRelax_Salt(:,:)
        CALL dbg_print('init relaxation - S'       ,z_c   ,str_module,idt_src, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
      ELSE
        CALL finish(TRIM(routine),' type_surfRelax_Salt>0 and no_tracer<2 - ABORT')
      END IF
    END IF

    CALL message( TRIM(routine),'end' )

  END SUBROUTINE init_ho_relaxation
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE  init_3Drelax_coefficient_NS(patch_3D, relax_coefficient, SouthBoundary, NorthBoundary, relaxWidth)

    TYPE(t_patch_3d ),TARGET :: patch_3D
    REAL(wp) :: relax_coefficient(:,:,:)
    REAL(wp) :: SouthBoundary, NorthBoundary, relaxWidth

    INTEGER :: jb, jc, start_cell_index, end_cell_index
    REAL(wp) :: lat_diff
    TYPE(t_patch),POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    !-------------------------------------------------------------------------
    patch_2d  => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    relax_coefficient(:,:,:) = 0.0_wp
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        lat_diff = patch_2d%cells%center(jc,jb)%lat - SouthBoundary
        IF (lat_diff > relaxWidth) & ! check the north boundary
          lat_diff = NorthBoundary - patch_2d%cells%center(jc,jb)%lat
          
        IF (lat_diff >= 0.0_wp .AND. lat_diff <  relaxWidth) THEN
          relax_coefficient(jc,:,jb) = (relaxWidth - lat_diff) / relaxWidth
        ENDIF
        
      END DO
    END DO

  END SUBROUTINE  init_3Drelax_coefficient_NS
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE init_ocean_WindForcing(patch_3d, ocean_state, atmos_fluxes, fu10)
    
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_atmos_fluxes)                    :: atmos_fluxes
    REAL(wp), INTENT(INOUT)                 :: fu10(:,:)      !  windspeed: p_as%fu10

    TYPE(t_patch), POINTER        :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges, owned_edges
    REAL(wp):: windStress_u(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp):: windStress_v(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    INTEGER :: start_edge_index, end_edge_index !,i_startblk_c, i_endblk_c,
    INTEGER :: je,blockNo, level

    patch_2d    => patch_3d%p_patch_2d(1)
    all_edges   => patch_2d%edges%all
    owned_edges => patch_2d%edges%owned

    CALL set_windstress(patch_2d, windStress_u, forcing_windstress_u_type, &
      & forcing_windStress_u_amplitude, forcing_windstress_zonal_waveno, forcing_windstress_merid_waveno)

    CALL set_windstress(patch_2d, windStress_v, forcing_windstress_v_type, &
      & forcing_windStress_u_amplitude, forcing_windstress_zonal_waveno, forcing_windstress_merid_waveno)

    ! calculate normal edge stress from u,v
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, je, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_edges%start_block,all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)
      ocean_state%p_aux%bc_top_WindStress(:,blockNo) = 0.0_wp
      DO je =  start_edge_index, end_edge_index
        DO level=1,MAX(patch_3D%p_patch_1D(1)%dolic_e(je,blockNo),1)
        
          ocean_state%p_aux%bc_top_WindStress(je,blockNo) = &
            & windStress_u(je,blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v1 + &
            & windStress_v(je,blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v2
            
          ENDDO
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
    
    
!     field_2d(:,:) = MERGE(amplitude * COS(zonal_waveno*pi*(lat(:,:)-center)/length),0.0_wp,mask(:,:) <= threshold)
!     CALL set_windspeed(all_cells, patch_3d%lsm_c(:,1,:), sea_boundary, fu10,                                      &
!       & atmos_fluxes%topBoundCond_windStress_u, atmos_fluxes%topBoundCond_windStress_v, forcing_windspeed_type,   &
!       & forcing_windspeed_amplitude)

    idt_src=0  ! output print level - 0: print in any case
    CALL dbg_print('init zonal wind stress'    ,windStress_u,str_module,idt_src,in_subset=owned_edges)
    CALL dbg_print('init merid. wind stress'   ,windStress_v,str_module,idt_src,in_subset=owned_edges)
!     CALL dbg_print('init wind speed'           ,fu10                                  ,str_module,idt_src,in_subset=owned_cells)

  END SUBROUTINE init_ocean_WindForcing

!<Optimize:inUse>
  SUBROUTINE set_windstress(patch_2D, windstress, &
      &                     control, amplitude, zonal_waveno, meridional_waveno,center,length)
    TYPE(t_patch), POINTER           :: patch_2d
    REAL(wp), INTENT(INOUT)          :: windstress(:,:)
    INTEGER,  INTENT(IN)             :: control
    REAL(wp), INTENT(IN)             :: amplitude,zonal_waveno,meridional_waveno
    REAL(wp), INTENT(IN),OPTIONAL    :: center, length

!     REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp), TARGET                 :: lon(nproma,patch_2D%nblks_e), lat(nproma,patch_2D%nblks_e)

    lat(:,:) = patch_2D%edges%center(:,:)%lat
    lon(:,:) = patch_2D%edges%center(:,:)%lon
    
    SELECT CASE (control)
    CASE (0) ! NO FORCING, SET TO ZERO ========================================
      windstress = 0.0_wp
    CASE (1:100)      ! FILE INPUT, DONE ELSEWHERE ============================
      windstress = 0.0_wp
      CALL message('windstress forcing','file input')
    CASE (101:200)    ! ANALYTIC SETUP ========================================
      SELECT CASE (control)
      CASE(101)       ! constant amplitude
        windstress = amplitude
      CASE(102)       ! basin setup, zonally changed
        CALL basin_zonal(lon, lat, windstress,amplitude,length)
      CASE(103)       ! basin setup, meridionally changed
        CALL basin_meridional(lon, lat, windstress,amplitude,length)
      CASE(104)       ! zonally periodic, nonzero at pols, meridionally constant
        CALL zonal_periodic_nonzero_around_center_zero_at_pols(lon, lat,  windstress, amplitude)
      CASE(105)
        CALL meridional_periodic_around_center_zero_at_pols(lon, lat, windstress, amplitude)
      CASE(106)       ! zonally periodic around a given center, zero at pols, meridionally constant
        CALL zonal_periodic_zero_at_pols(lon, lat, windstress,amplitude,zonal_waveno)
      CASE(107)       ! latteral cells, zonal period only
        CALL cells_zonal_periodic(lon, lat, windstress,amplitude,zonal_waveno)
      CASE(108)       ! latteral cells, zonally and meridionally periodic
        CALL cells_zonal_and_meridional_periodic(lon, lat, windstress,amplitude,zonal_waveno,meridional_waveno)
      CASE(109)
        CALL cells_zonal_and_meridional_periodic_constant_amplitude_sin(lon, lat, windstress, amplitude)
      CASE(110)
        CALL cells_zonal_and_meridional_periodic_constant_amplitude_cosin(lon, lat, windstress, amplitude)
      CASE(111)
        CALL Wolfe_Cessi_TestCase(lon, lat, windstress, amplitude)
      CASE(112)       ! basin setup, zonally changed for Abernathey test case
        CALL basin_zonal_zeroOutside(lon, lat, windstress,amplitude,length)     
      CASE(113)
        CALL zentral_jet(lon, lat, windstress, amplitude)
        CASE(114)
          CALL wind_shear_u(patch_2D, lon, lat, windstress, amplitude)		        
          CASE(115)
            CALL wind_shear_v(patch_2D, lon, lat, windstress, amplitude)		        
				        
      END SELECT
    END SELECT
    

  END SUBROUTINE set_windstress

  SUBROUTINE basin_zonal(lon, lat, field_2d, amplitude,length_opt, zonal_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN) , OPTIONAL  :: length_opt,zonal_waveno_opt

    REAL(wp) :: length, zonal_waveno

    length = basin_height_deg * deg2rad
    zonal_waveno = forcing_windstress_zonal_waveno

    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(zonal_waveno,zonal_waveno_opt)

    field_2d(:,:) = amplitude * COS(zonal_waveno*pi*(lat(:,:)-length)/length)

  END SUBROUTINE basin_zonal

  SUBROUTINE basin_zonal_zeroOutside(lon, lat, field_2d, amplitude, center_opt, length_opt, zonal_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN) , OPTIONAL  :: center_opt,length_opt,zonal_waveno_opt

    REAL(wp) :: center, length, zonal_waveno

    center = basin_center_lat * deg2rad
    length = basin_height_deg * deg2rad
    zonal_waveno = forcing_windstress_zonal_waveno

    CALL assign_if_present(center,center_opt)
    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(zonal_waveno,zonal_waveno_opt)

    field_2d(:,:) = amplitude * COS(zonal_waveno*pi*(lat(:,:)-center)/length)
    field_2d(:,:) = MERGE(field_2d(:,:),0.0_wp, lat(:,:) > center-0.5_wp*length)
    field_2d(:,:) = MERGE(field_2d(:,:),0.0_wp, lat(:,:) < center+0.5_wp*length)

  END SUBROUTINE basin_zonal_zeroOutside

  SUBROUTINE basin_meridional(lon, lat, field_2d, amplitude,length_opt,meridional_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN) , OPTIONAL  :: length_opt,meridional_waveno_opt

    REAL(wp) :: length, meridional_waveno

    length            = basin_width_deg * deg2rad
    meridional_waveno = forcing_windstress_merid_waveno

    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(meridional_waveno,meridional_waveno_opt)

    field_2d(:,:) = amplitude * COS(meridional_waveno*pi*(lon(:,:)-length)/length)
  END SUBROUTINE basin_meridional

  SUBROUTINE zonal_periodic_nonzero_around_center_zero_at_pols(lon, lat, field_2d, amplitude,&
      & center_opt,length_opt,zonal_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: center_opt, length_opt, zonal_waveno_opt

    REAL(wp) :: center, length, zonal_waveno

    length       = 180.0_wp * deg2rad
    center       = -60.0_wp * deg2rad
    zonal_waveno = forcing_windstress_zonal_waveno

    CALL assign_if_present(center,center_opt)
    CALL assign_if_present(length,length_opt)
    CaLL assign_if_present(zonal_waveno, zonal_waveno_opt)

    field_2d(:,:) = amplitude * COS(zonal_waveno*pi*(lat(:,:)-center)/length)
  END SUBROUTINE zonal_periodic_nonzero_around_center_zero_at_pols

  SUBROUTINE meridional_periodic_around_center_zero_at_pols(lon, lat, field_2d, amplitude, &
      & center_opt,length_opt,meridional_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: center_opt, length_opt, meridional_waveno_opt

    REAL(wp) :: center, length, meridional_waveno

    length                        =  90.0_wp * deg2rad
    center                        = -20.0_wp * deg2rad
    meridional_waveno             = forcing_windstress_merid_waveno

    CALL assign_if_present(center,center_opt)
    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(meridional_waveno, meridional_waveno_opt)

    field_2d(:,:) = amplitude * COS(meridional_waveno*pi*(lon(:,:)-center)/length)
  END SUBROUTINE meridional_periodic_around_center_zero_at_pols

  SUBROUTINE zonal_periodic_zero_at_pols(lon, lat, field_2d, amplitude, zonal_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: zonal_waveno_opt

    REAL(wp) :: zonal_waveno

    zonal_waveno = forcing_windstress_zonal_waveno

    CaLL assign_if_present(zonal_waveno, zonal_waveno_opt)

    field_2d(:,:) = amplitude * COS(lat(:,:)) * &
      & COS(forcing_windstress_zonalWavePhase + zonal_waveno * ABS(lat(:,:)))

  END SUBROUTINE zonal_periodic_zero_at_pols

  SUBROUTINE Wolfe_Cessi_TestCase(lon, lat, field_2d, amplitude)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude

    field_2d(:,:) = amplitude *                                   &
      & (0.8_wp*EXP(-1*lat(:,:)*lat(:,:)/0.01) - COS(6.0_wp*lat(:,:)))

  END SUBROUTINE Wolfe_Cessi_TestCase

  SUBROUTINE zentral_jet(lon, lat, field_2d, amplitude)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude

    field_2d(:,:) = amplitude *                                   &
      & (0.8_wp*EXP(-1*lat(:,:)*lat(:,:)/0.5) )

  END SUBROUTINE zentral_jet

  SUBROUTINE wind_shear_u(patch_2d, lon, lat, field_2d, amplitude)
    TYPE(t_patch), POINTER           :: patch_2d	  
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude

	INTEGER :: edge_index,edge_block, start_edges_index, end_edges_index
    TYPE(t_subset_range), POINTER :: all_edges	
	REAL(wp) :: point_lon, point_lat, uu, vv	
		
		
    all_edges => patch_2d%edges%ALL		
    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index
      
         point_lon = patch_2d%edges%center(edge_index,edge_block)%lon
         point_lat = patch_2d%edges%center(edge_index,edge_block)%lat
		 
      	 IF(point_lat>=basin_center_lat)THEN		 
  		   !uu=tanh((point_lat-0.025)*30.0_wp)	 
 		   uu=tanh((point_lat)*30.0_wp)	 		   
  		 ELSEIF(point_lat<basin_center_lat)THEN
 		   !uu=tanh((0.75_wp-point_lat)*30.0_wp)	 			 
  		   uu=tanh((-point_lat)*30.0_wp)	 
  		 ENDIF
  		 vv=0.1_wp*sin(2.0_wp*pi*point_lon)	   	 
		   
!      	 field_2d(edge_index,edge_block) =amplitude*(uu * patch_2d%edges%primal_normal(edge_index,edge_block)%v1 &
!    		         & + vv * patch_2d%edges%primal_normal(edge_index,edge_block)%v2) 
        field_2d(edge_index,edge_block) =amplitude*uu * patch_2d%edges%primal_normal(edge_index,edge_block)%v1 

		 
        ENDDO
      ENDDO

  END SUBROUTINE wind_shear_u
  
  SUBROUTINE wind_shear_v(patch_2d, lon, lat, field_2d, amplitude)
    TYPE(t_patch), POINTER           :: patch_2d	  
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude

	INTEGER :: edge_index,edge_block, start_edges_index, end_edges_index
    TYPE(t_subset_range), POINTER :: all_edges	
	REAL(wp) :: point_lon, point_lat, uu, vv	
		
		
    all_edges => patch_2d%edges%ALL		
    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index
      
         point_lon = patch_2d%edges%center(edge_index,edge_block)%lon
         point_lat = patch_2d%edges%center(edge_index,edge_block)%lat
		 
      	 !IF(point_lat>=basin_center_lat)THEN		 
  		 !  uu=tanh((point_lat-0.025)*30.0_wp)	 
  		 !ELSEIF(point_lat<basin_center_lat)THEN
  		 !  uu=tanh((0.75_wp-point_lat)*30.0_wp)	 
  		 !ENDIF
  		 vv=0.1_wp*sin(2.0_wp*pi*point_lon)	   	 
		   
      	 field_2d(edge_index,edge_block) =amplitude*vv * patch_2d%edges%primal_normal(edge_index,edge_block)%v2
		 
        ENDDO
      ENDDO



  END SUBROUTINE wind_shear_v


  SUBROUTINE cells_zonal_periodic(lon, lat, field_2d, amplitude, zonal_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: zonal_waveno_opt

    REAL(wp) :: zonal_waveno

    zonal_waveno = forcing_windstress_zonal_waveno
    CALL assign_if_present(zonal_waveno, zonal_waveno_opt)

    field_2d(:,:) = amplitude &
                & * COS(lat) &
                & * COS(zonal_waveno * lat) &
                & * COS(lon)
  END SUBROUTINE cells_zonal_periodic

  SUBROUTINE cells_zonal_and_meridional_periodic(lon, lat, field_2d, amplitude, &
      & zonal_waveno_opt,meridional_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: zonal_waveno_opt, meridional_waveno_opt

    REAL(wp) :: zonal_waveno, meridional_waveno

    zonal_waveno      = forcing_windstress_zonal_waveno
    meridional_waveno = forcing_windstress_merid_waveno

    CALL assign_if_present(zonal_waveno, zonal_waveno_opt)
    CALL assign_if_present(meridional_waveno, meridional_waveno_opt)

    field_2d(:,:) = amplitude &
      & * COS(lat(:,:)) &
      & * COS(zonal_waveno * lat(:,:)) &
      & * SIN(meridional_waveno * lon(:,:))

  END SUBROUTINE cells_zonal_and_meridional_periodic

  SUBROUTINE cells_zonal_and_meridional_periodic_constant_amplitude_sin(lon, lat, field_2d, &
      & amplitude,length_opt, zonal_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp)                         :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN) , OPTIONAL  :: length_opt,zonal_waveno_opt

    REAL(wp) :: length, zonal_waveno

    length       = basin_height_deg * deg2rad
    zonal_waveno = forcing_windstress_zonal_waveno

    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(zonal_waveno,zonal_waveno_opt)

    ! -length/length ???
!     strength = amplitude*cos(zonal_waveno*pi*lat(:,:)-length/length)

    field_2d(:,:) = amplitude*amplitude*cos(zonal_waveno*pi*lat(:,:)-length/length)*sin(lon(:,:))

  END SUBROUTINE cells_zonal_and_meridional_periodic_constant_amplitude_sin

  SUBROUTINE cells_zonal_and_meridional_periodic_constant_amplitude_cosin(lon, lat, field_2d, &
      & amplitude,length_opt, zonal_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp)                         :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN) , OPTIONAL  :: length_opt,zonal_waveno_opt

    REAL(wp) :: length, zonal_waveno

    length       = basin_height_deg * deg2rad
    zonal_waveno = forcing_windstress_zonal_waveno

    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(zonal_waveno,zonal_waveno_opt)

    ! -length/length ???
!     strength = amplitude*cos(zonal_waveno*pi*lat(:,:)-length/length)

    field_2d(:,:) = amplitude*amplitude*cos(zonal_waveno*pi*lat(:,:)-length/length)*cos(lon(:,:))

  END SUBROUTINE cells_zonal_and_meridional_periodic_constant_amplitude_cosin
  !-------------------------------------------------------------------------

!<Optimize:inUse>
  SUBROUTINE set_windspeed( windspeed, windstress_u, windstress_v, control, amplitude)
    REAL(wp), INTENT(INOUT)          :: windspeed(:,:)
    REAL(wp), INTENT(INOUT)          :: windstress_u(:,:)
    REAL(wp), INTENT(INOUT)          :: windstress_v(:,:)
    INTEGER,  INTENT(IN)             :: control
    REAL(wp), INTENT(IN)             :: amplitude

    SELECT CASE (control)
    CASE (0) ! NO FORCING, SET TO ZERO ========================================
      windspeed = 0.0_wp
    CASE (1:100)      ! FILE INPUT, DONE ELSEWHERE ============================
      CALL message('windspeed forcing','file input')
    CASE (101:200)    ! ANALYTIC SETUP ========================================
      SELECT CASE (control)
      CASE(101)       ! constant amplitude
        windspeed(:,:) = amplitude
      CASE(102)       ! windspeed backward calculated from windstress amplitude
        CALL calc_windspeed_fromwindstress( windspeed, windstress_u, windstress_v, amplitude)
   !  CASE(103)       ! basin setup, zonally changed
   !    CALL basin_zonal(subset,mask,threshold,windstress,amplitude,length)
      END SELECT
    END SELECT

  END SUBROUTINE set_windspeed
  !-------------------------------------------------------------------------
  !>
  !! Initialization of wind speed if wind stress is available
  !!  Following common relation of wind stress to wind speed (e.g. Smith 1980)
  !!  
  !
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M, 2014-07
  !
  !-------------------------------------------------------------------------
  !
  SUBROUTINE calc_windspeed_fromwindstress( windspeed, windstress_u, windstress_v, amplitude)
  !!
  !!  
    REAL(wp)                         :: windspeed(:,:)
    REAL(wp)                         :: windstress_u(:,:)
    REAL(wp)                         :: windstress_v(:,:)
    REAL(wp), INTENT(IN)             :: amplitude

    ! |tau|  = Rho_air * C_D * U_10**2 
    ! |U_10| = ( |tau| / (Rho_air*C_D) )**(1/2) or
    ! |U_10| = (Rho_air*C_D)**(-1/2) * |tau|**(1/2)
    !  with Rho_air=1.22 kg/m3, drag coefficient C_D ~ 1.3e-3, U_10 the wind speed in 10m height
    !  (Rho_air * C_D)**(-1/2) ~ (1.5 e-3 )**(-1/2) ~ 26
    !  or: |U_10| ~ 26 m/s for tau=1 Pa, or 8.2 m/s for tau=0.1 Pa
    !  |tau| = (tau_x**2 + tau_y**2)**(1/2)

    windspeed(:,:) = amplitude*26.0_wp*                                                                 &
      &              SQRT(SQRT(windstress_u(:,:)*windstress_u(:,:) + windstress_v(:,:)*windstress_v(:,:)))

  END SUBROUTINE calc_windspeed_fromwindstress

  SUBROUTINE nf(STATUS)

    INTEGER, INTENT(in) :: STATUS

    IF (STATUS /= nf_noerr) THEN
      CALL finish('mo_ext_data netCDF error', nf_strerror(STATUS))
    ENDIF

  END SUBROUTINE nf
  !-------------------------------------------------------------------------
END MODULE mo_ocean_forcing
