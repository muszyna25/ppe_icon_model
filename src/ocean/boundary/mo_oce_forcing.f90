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
MODULE mo_oce_forcing
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: filename_max
  USE mo_grid_config,         ONLY: nroot
  USE mo_parallel_config,     ONLY: nproma
  USE mo_ocean_nml,           ONLY: basin_height_deg, basin_width_deg, no_tracer,   &
    & forcing_windstress_zonal_waveno, forcing_windstress_merid_waveno,             &
    & init_oce_relax, type_3dimRelax_Salt, type_3dimRelax_Temp,                     &
    & type_surfRelax_Salt, type_surfRelax_Temp,                                     &
    & forcing_windStress_u_amplitude, forcing_windStress_v_amplitude,               &
    & forcing_windstress_u_type, forcing_windstress_v_type,    &
#ifdef __SX__
    & forcing_windstress_zonalWavePhas,                        &
#else
    & forcing_windstress_zonalWavePhase,                       &
#endif
    & relax_temperature_min, relax_temperature_max,            &
    & forcing_temperature_poleLat
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_math_constants,      ONLY: pi, deg2rad, pi_2
  USE mo_impl_constants,      ONLY: max_char_length, sea_boundary, success
  USE mo_math_utilities,      ONLY: gvec2cvec, cvec2gvec, t_cartesian_coordinates
  USE mo_sea_ice_types,       ONLY: t_sfc_flx
  USE mo_oce_state,           ONLY: set_oce_tracer_info
  USE mo_oce_types,           ONLY: t_hydro_ocean_state
  USE mo_dynamics_config,     ONLY: nold

  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_var_list,            ONLY: add_var, add_ref
  USE mo_var_metadata,        ONLY: groups
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_fortran_tools,       ONLY: assign_if_present
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_netcdf_read,        ONLY: read_netcdf_data

  IMPLICIT NONE
  PRIVATE
  INCLUDE 'netcdf.inc'

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
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
!<Optimize_Used>
  SUBROUTINE construct_ocean_forcing(p_patch, p_sfc_flx, var_list)
    !
    TYPE(t_patch),   INTENT(IN)    :: p_patch
    TYPE(t_sfc_flx), INTENT(INOUT) :: p_sfc_flx
    TYPE(t_var_list),INTENT(INOUT) :: var_list

    ! Local variables
    INTEGER                        :: alloc_cell_blocks, ist, jtrc, i
    CHARACTER(len=max_char_length) :: oce_tracer_names(no_tracer),&
    &                                 oce_tracer_units(no_tracer),&
    &                                 oce_tracer_longnames(no_tracer)
    INTEGER                        :: oce_tracer_codes(no_tracer)

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_ocean_forcing'

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    alloc_cell_blocks = p_patch%alloc_cell_blocks

    CALL add_var(var_list, 'topBoundCond_windStress_u', p_sfc_flx%topBoundCond_windStress_u , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('topBoundCond_windStress_u', 'Pa', 'topBoundCond_windStress_u', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'topBoundCond_windStress_v', p_sfc_flx%topBoundCond_windStress_v , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('topBoundCond_windStress_v', 'Pa', 'topBoundCond_windStress_v', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'HeatFlux_Total', p_sfc_flx%HeatFlux_Total , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Total', 'W/m2', 'HeatFlux_Total', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'FrshFlux_VolumeTotal', p_sfc_flx%FrshFlux_VolumeTotal , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_VolumeTotal', 'm/s', 'FrshFlux_VolumeTotal', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'HeatFlux_ShortWave', p_sfc_flx%HeatFlux_ShortWave , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_ShortWave', 'W/m2', 'HeatFlux_ShortWave', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'HeatFlux_LongWave', p_sfc_flx%HeatFlux_LongWave , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_LongWave', 'W/m2', 'HeatFlux_LongWave', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'HeatFlux_Sensible', p_sfc_flx%HeatFlux_Sensible , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Sensible', 'W/m2', 'HeatFlux_Sensible', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'HeatFlux_Latent', p_sfc_flx%HeatFlux_Latent , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Latent', 'W/m2', 'HeatFlux_Latent', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'FrshFlux_Precipitation', p_sfc_flx%FrshFlux_Precipitation , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Precipitation', 'm/s', 'FrshFlux_Precipitation', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'FrshFlux_SnowFall', p_sfc_flx%FrshFlux_SnowFall , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_SnowFall', 'm/s', 'FrshFlux_SnowFall', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'FrshFlux_Evaporation', p_sfc_flx%FrshFlux_Evaporation , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Evaporation', 'm/s', 'FrshFlux_Evaporation', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'FrshFlux_Runoff', p_sfc_flx%FrshFlux_Runoff , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Runoff', 'm/s', 'FrshFlux_Runoff', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'FrshFlux_TotalSalt', p_sfc_flx%FrshFlux_TotalSalt , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_TotalSalt', 'm/s', 'FrshFlux_TotalSalt', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'FrshFlux_TotalOcean', p_sfc_flx%FrshFlux_TotalOcean , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_TotalOcean', 'm/s', 'FrshFlux_TotalOcean', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'FrshFlux_TotalIce', p_sfc_flx%FrshFlux_TotalIce , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_TotalIce', 'm/s', 'FrshFlux_TotalIce', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_hfrelax', p_sfc_flx%forc_hfrelax , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_hfrelax', 'm/s', 'forc_hfrelax', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'forc_fwrelax', p_sfc_flx%forc_fwrelax , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fwrelax', 'm/s', 'forc_fwrelax', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    CALL add_var(var_list, 'FrshFlux_VolumeIce', p_sfc_flx%FrshFlux_VolumeIce, &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_VolumeIce', 'm/s', 'FrshFlux_VolumeIce', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))
    IF(no_tracer>=1) THEN
      ! there are four tracer related fields: tracer focing, tracer relaxation
      ! and both accumulated
      CALL add_var(var_list, 'topBoundCond_Temp_vdiff', p_sfc_flx%topBoundCond_Temp_vdiff, &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('topBoundCond_Temp_vdiff', 'K*m/s', 'topBoundCond_Temp_vdiff', DATATYPE_FLT32),&
        &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      CALL add_var(var_list, 'topBoundCond_Salt_vdiff', p_sfc_flx%topBoundCond_Salt_vdiff, &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('topBoundCond_Salt_vdiff', 'psu*m/s', 'topBoundCond_Salt_vdiff', DATATYPE_FLT32),&
        &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      CALL add_var(var_list, 'data_surfRelax_Temp', p_sfc_flx%data_surfRelax_Temp, &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('data_surfRelax_Temp', 'C', 'data_surfRelax_Temp', DATATYPE_FLT32),&
        &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      CALL add_var(var_list, 'data_surfRelax_Salt', p_sfc_flx%data_surfRelax_Salt, &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('data_surfRelax_Salt', 'psu', 'data_surfRelax_Salt', DATATYPE_FLT32),&
        &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      CALL add_var(var_list, 'topBoundCond_Temp_vdiff_acc', p_sfc_flx%topBoundCond_Temp_vdiff_acc , &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('topBoundCond_Temp_vdiff_acc', 'K*m/s', 'topBoundCond_Temp_vdiff_acc', DATATYPE_FLT32),&
        &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      CALL add_var(var_list, 'topBoundCond_Salt_vdiff_acc', p_sfc_flx%topBoundCond_Salt_vdiff_acc , &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('topBoundCond_Salt_vdiff_acc', 'psu*m/s', 'topBoundCond_Salt_vdiff_acc', DATATYPE_FLT32),&
        &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      CALL add_var(var_list, 'data_surfRelax_Temp_acc', p_sfc_flx%data_surfRelax_Temp_acc , &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('data_surfRelax_Temp_acc', 'C', 'data_surfRelax_Temp_acc', DATATYPE_FLT32),&
        &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      CALL add_var(var_list, 'data_surfRelax_Salt_acc', p_sfc_flx%data_surfRelax_Salt_acc , &
        &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &          t_cf_var('data_surfRelax_Salt_acc', 'psu', 'data_surfRelax_Salt_acc', DATATYPE_FLT32),&
        &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,alloc_cell_blocks/), &
        &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

   !  ALLOCATE(p_sfc_flx%tracer_ptr(no_tracer*4))

   !  CALL set_oce_tracer_info(no_tracer           , &
   !    &                      oce_tracer_names    , &
   !    &                      oce_tracer_longnames, &
   !    &                      oce_tracer_codes    , &
   !    &                      oce_tracer_units)
   !  DO jtrc = 1,no_tracer
   !    CALL add_ref( var_list, 'forc_tracer', &
   !      &           'forc_tracer_'//TRIM(oce_tracer_names(jtrc)),          &
   !      &           p_sfc_flx%tracer_ptr(jtrc)%p,    &
   !      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,&
   !      &           t_cf_var('forc_tracer'//TRIM(oce_tracer_names(jtrc)), &
   !      &                    oce_tracer_units(jtrc), &
   !      &                    'forcing: '//TRIM(oce_tracer_longnames(jtrc)), DATATYPE_FLT32), &
   !      &           t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
   !      &           ldims=(/nproma,alloc_cell_blocks/))
   !  END DO
   !  DO jtrc = no_tracer+1,2*no_tracer
   !    i = jtrc - no_tracer
   !    CALL add_ref( var_list, 'forc_tracer_relax', &
   !      &           'forc_tracer_relax_'//TRIM(oce_tracer_names(i)),          &
   !      &           p_sfc_flx%tracer_ptr(jtrc)%p,    &
   !      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,&
   !      &           t_cf_var('forc_tracer_relax'//TRIM(oce_tracer_names(i)), &
   !      &                    oce_tracer_units(i), &
   !      &                    'forcing relaxation accumulated: '//TRIM(oce_tracer_longnames(i)), DATATYPE_FLT32), &
   !      &           t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
   !      &           ldims=(/nproma,alloc_cell_blocks/))
   !  END DO
   !  DO jtrc = (2*no_tracer)+1,3*no_tracer
   !    i = jtrc - 2*no_tracer
   !    CALL add_ref( var_list, 'forc_tracer_acc', &
   !      &           'forc_tracer_acc_'//TRIM(oce_tracer_names(i)),          &
   !      &           p_sfc_flx%tracer_ptr(jtrc)%p,    &
   !      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,&
   !      &           t_cf_var('forc_tracer_acc'//TRIM(oce_tracer_names(i)), &
   !      &                    oce_tracer_units(i), &
   !      &                    'forcing accumulated: '//TRIM(oce_tracer_longnames(i)), DATATYPE_FLT32), &
   !      &           t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
   !      &           ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
   !  END DO
   !  DO jtrc = (3*no_tracer)+1,4*no_tracer
   !    i = jtrc - 3*no_tracer
   !    CALL add_ref( var_list, 'forc_tracer_relax_acc', &
   !      &           'forc_tracer_relax_acc_'//TRIM(oce_tracer_names(i)),          &
   !      &           p_sfc_flx%tracer_ptr(jtrc)%p,    &
   !      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,&
   !      &           t_cf_var('forc_tracer_relax_acc'//TRIM(oce_tracer_names(i)), &
   !      &                    oce_tracer_units(i), &
   !      &                    'forcing relaxation accumulated: '//TRIM(oce_tracer_longnames(i)), DATATYPE_FLT32), &
   !      &           t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
   !      &           ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
   !  END DO
    ENDIF
    CALL add_var(var_list, 'topBoundCond_windStress_u_acc', p_sfc_flx%topBoundCond_windStress_u_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('topBoundCond_windStress_u_acc', 'Pa', 'topBoundCond_windStress_u_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'topBoundCond_windStress_v_acc', p_sfc_flx%topBoundCond_windStress_v_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('topBoundCond_windStress_v_acc', 'Pa', 'topBoundCond_windStress_v_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'HeatFlux_Total_acc', p_sfc_flx%HeatFlux_Total_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Total_acc', 'W/m2', 'HeatFlux_Total_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'FrshFlux_VolumeTotal_acc', p_sfc_flx%FrshFlux_VolumeTotal_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_VolumeTotal_acc', 'm/s', 'FrshFlux_VolumeTotal_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'HeatFlux_ShortWave_acc', p_sfc_flx%HeatFlux_ShortWave_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_ShortWave_acc', 'W/m2', 'HeatFlux_ShortWave_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'HeatFlux_LongWave_acc', p_sfc_flx%HeatFlux_LongWave_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_LongWave_acc', 'W/m2', 'HeatFlux_LongWave_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'HeatFlux_Sensible_acc', p_sfc_flx%HeatFlux_Sensible_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Sensible_acc', 'W/m2', 'HeatFlux_Sensible_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'HeatFlux_Latent_acc', p_sfc_flx%HeatFlux_Latent_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('HeatFlux_Latent_acc', 'W/m2', 'HeatFlux_Latent_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'FrshFlux_Precipitation_acc', p_sfc_flx%FrshFlux_Precipitation_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Precipitation_acc', 'm/s', 'FrshFlux_Precipitation_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'FrshFlux_SnowFall_acc', p_sfc_flx%FrshFlux_SnowFall_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_SnowFall_acc', 'm/s', 'FrshFlux_SnowFall_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'FrshFlux_Evaporation_acc', p_sfc_flx%FrshFlux_Evaporation_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Evaporation_acc', 'm/s', 'FrshFlux_Evaporation_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default", "oce_force_essentials"))
    CALL add_var(var_list, 'FrshFlux_Runoff_acc', p_sfc_flx%FrshFlux_Runoff_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_Runoff_acc', 'm/s', 'FrshFlux_Runoff_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'FrshFlux_TotalSalt_acc', p_sfc_flx%FrshFlux_TotalSalt_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_TotalSalt_acc', 'm/s', 'FrshFlux_TotalSalt_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'FrshFlux_TotalOcean_acc', p_sfc_flx%FrshFlux_TotalOcean_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_TotalOcean_acc', 'm/s', 'FrshFlux_TotalOcean_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'FrshFlux_TotalIce_acc', p_sfc_flx%FrshFlux_TotalIce_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_TotalIce_acc', 'm/s', 'FrshFlux_TotalIce_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'forc_hfrelax_acc', p_sfc_flx%forc_hfrelax_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_hfrelax_acc', 'm/s', 'forc_hfrelax_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'forc_fwrelax_acc', p_sfc_flx%forc_fwrelax_acc , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('forc_fwrelax_acc', 'm/s', 'forc_fwrelax_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups("oce_default"))
    CALL add_var(var_list, 'FrshFlux_VolumeIce_acc', p_sfc_flx%FrshFlux_VolumeIce_acc, &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_VolumeIce_acc', 'm/s', 'FrshFlux_VolumeIce_acc', DATATYPE_FLT32),&
    &          t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/))

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
!<Optimize_Used>
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
  !>
  !! Initialization of temperature and salinity relaxation for the hydrostatic ocean model.
  !! Temperature and salinity relaxation data are read from external data
  !
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M, 2011-11
  !
  !-------------------------------------------------------------------------
  !
  SUBROUTINE init_ho_relaxation(patch_2d, patch_3d, ocean_state, p_sfc_flx)

    TYPE(t_patch),TARGET, INTENT(in)  :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_sfc_flx)                   :: p_sfc_flx

    ! Local Variables

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_initial_conditions:init_ho_relaxation'
    CHARACTER(filename_max) :: relax_init_file   !< file name for reading in

    LOGICAL :: l_exist
    INTEGER :: i_lev, no_cells, no_levels, jb, jc, jk
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: ncid, dimid

    REAL(wp):: z_c(nproma,1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp):: z_surfRelax(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: temperature_difference, poleLat, waveNo

    TYPE(t_subset_range), POINTER :: all_cells
    !-------------------------------------------------------------------------
    all_cells => patch_2d%cells%ALL

    ! Read relaxation data from file
    IF (init_oce_relax == 1) THEN

      ! sphere_radius = grid_sphere_radius
      ! u0 =(2.0_wp*pi*sphere_radius)/(12.0_wp*24.0_wp*3600.0_wp)

      CALL message (TRIM(routine), 'start')

      i_lev        = patch_2d%level

      IF (my_process_is_stdio()) THEN
        !
        ! Relaxation variables are read from relax_init_file
        WRITE (relax_init_file,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',i_lev, '-relax.nc'

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

      ENDIF  !  stdio


      !-------------------------------------------------------
      !
      ! Read ocean relaxation data at cells
      !
      !-------------------------------------------------------

      ! triangle center and edges

      ! read temperature
      !  - read one data set, annual mean only
      !  - "T": annual mean temperature
      CALL read_netcdf_data (ncid, 'T', patch_2d%n_patch_cells_g, patch_2d%n_patch_cells, &
        & patch_2d%cells%decomp_info%glb_index, z_surfRelax)

      IF (no_tracer>=1) THEN
        p_sfc_flx%data_surfRelax_Temp(:,:) = z_surfRelax(:,:)
      ELSE
        CALL message( TRIM(routine),'WARNING: no tracer used, but init relaxation attempted')
      END IF

      ! read salinity
      !  - "S": annual mean salinity
      IF (no_tracer > 1) THEN
        CALL read_netcdf_data (ncid, 'S', patch_2d%n_patch_cells_g, patch_2d%n_patch_cells, &
          & patch_2d%cells%decomp_info%glb_index, z_surfRelax)
        p_sfc_flx%data_surfRelax_Salt(:,:) = z_surfRelax(:,:)
      END IF

      ! close file
      IF(my_process_is_stdio()) CALL nf(nf_close(ncid))

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          IF ( patch_3d%lsm_c(jc,1,jb) > sea_boundary ) THEN
            p_sfc_flx%data_surfRelax_Temp(jc,jb) = 0.0_wp
            IF (no_tracer>1) p_sfc_flx%data_surfRelax_Salt(jc,jb) = 0.0_wp
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
      p_sfc_flx%data_surfRelax_Temp(:,:) = ocean_state%p_prog(nold(1))%tracer(:,1,:,1)

    CASE(4)
      ! smooth ape relaxation, as in temperature_smoothAPE in mo_cean_initial_conditions
      temperature_difference = relax_temperature_max - relax_temperature_min
      poleLat = ABS(forcing_temperature_poleLat * deg2rad)
      waveNo = pi_2 / poleLat

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          DO jk=1, MIN(1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb))

            p_sfc_flx%data_surfRelax_Temp(jc,jb) = relax_temperature_min    + &
              & (COS(waveNo * MIN(ABS(patch_2d%cells%center(jc,jb)%lat), poleLat))**2) * temperature_difference

          END DO
        END DO
      END DO

    END SELECT

    IF (type_surfRelax_Salt == 3) THEN
      IF (no_tracer > 1) THEN
        p_sfc_flx%data_surfRelax_Salt(:,:) = ocean_state%p_prog(nold(1))%tracer(:,1,:,2)
      END IF
    END IF

    IF (type_3dimRelax_Temp == 3) THEN
      ocean_state%p_aux%data_3dimRelax_Temp(:,:,:) = ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
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
      z_c(:,1,:) = p_sfc_flx%data_surfRelax_Temp(:,:)
      CALL dbg_print('init relaxation - T'       ,z_c      ,str_module,idt_src, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    END IF
    IF (type_surfRelax_Salt > 0) THEN
      IF (no_tracer > 1) THEN
        z_c(:,1,:) = p_sfc_flx%data_surfRelax_Salt(:,:)
        CALL dbg_print('init relaxation - S'       ,z_c   ,str_module,idt_src, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
      ELSE
        CALL finish(TRIM(routine),' type_surfRelax_Salt>0 and no_tracer<2 - ABORT')
      END IF
    END IF

    CALL message( TRIM(routine),'end' )

  END SUBROUTINE init_ho_relaxation
  !-------------------------------------------------------------------------
!<Optimize_Used>
  SUBROUTINE init_ocean_forcing(patch_2d, patch_3d, ocean_state, p_sfc_flx)
    !
    TYPE(t_patch),TARGET, INTENT(in)        :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_sfc_flx)                         :: p_sfc_flx

    TYPE(t_subset_range), POINTER :: all_cells

    all_cells => patch_3d%p_patch_2d(1)%cells%All

    CALL set_windstress_u(all_cells, patch_3d%lsm_c(:,1,:), sea_boundary, p_sfc_flx%topBoundCond_windStress_u,&
      & forcing_windStress_u_amplitude, forcing_windstress_zonal_waveno, forcing_windstress_merid_waveno)

    CALL set_windstress_v(all_cells, patch_3d%lsm_c(:,1,:), sea_boundary, p_sfc_flx%topBoundCond_windStress_v,&
      & forcing_windStress_v_amplitude, forcing_windstress_zonal_waveno, forcing_windstress_merid_waveno)

    IF (init_oce_relax > 0) THEN
      CALL init_ho_relaxation(patch_2d, patch_3d, ocean_state, p_sfc_flx)
    END IF
  END SUBROUTINE init_ocean_forcing

!<Optimize_Used>
  SUBROUTINE set_windstress_u(subset, mask, threshold, windstress, &
      &                       amplitude, zonal_waveno, meridional_waveno, center, length)
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER, INTENT(IN)              :: mask(:,:)
    INTEGER, INTENT(IN)              :: threshold
    REAL(wp), INTENT(INOUT)          :: windstress(:,:)
    REAL(wp), INTENT(IN)             :: amplitude, zonal_waveno, meridional_waveno
    REAL(wp), INTENT(IN), OPTIONAL   :: center, length

    CALL set_windstress(subset, mask, threshold, windstress, &
      & forcing_windstress_u_type, amplitude, zonal_waveno, meridional_waveno, center, length)
  END SUBROUTINE set_windstress_u

!<Optimize_Used>
  SUBROUTINE set_windstress_v(subset, mask, threshold, windstress, &
      &                       amplitude, zonal_waveno, meridional_waveno, center, length)
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER, INTENT(IN)              :: mask(:,:)
    INTEGER, INTENT(IN)              :: threshold
    REAL(wp),INTENT(INOUT)           :: windstress(:,:)
    REAL(wp), INTENT(IN)             :: amplitude, zonal_waveno, meridional_waveno
    REAL(wp), INTENT(IN), OPTIONAL   :: center, length

    CALL set_windstress(subset, mask, threshold, windstress, &
      & forcing_windstress_v_type, amplitude, zonal_waveno, meridional_waveno, center, length)
  END SUBROUTINE set_windstress_v

!TODO disabled because of compiler error on blizzard
! SUBROUTINE update_cartesian_coords_for_cells(cell_subset,field_x,field_y, field_cc)
!   TYPE(t_subset_range),INTENT(IN)              :: cell_subset
!   REAL(wp), INTENT(IN)                         :: field_x(:,:), field_y(:,:)
!   TYPE(t_cartesian_coordinates), INTENT(OUT)   :: field_cc(:,:)
!   CALL gvec2cvec(  field_x(:,:),&
!     &              field_y(:,:),&
!     &              cell_subset%patch%cells%center(:,:)%lon,&
!     &              cell_subset%patch%cells%center(:,:)%lat,&
!     &              field_cc(:,:)%x(1),&
!     &              field_cc(:,:)%x(2),&
!     &              field_cc(:,:)%x(3))
! END SUBROUTINE update_cartesian_coords_for_cells
! SUBROUTINE update_from_cartesian_coords_for_cells(cell_subset,field_x,field_y, field_cc)
!   TYPE(t_subset_range),INTENT(IN)           :: cell_subset
!   REAL(wp), INTENT(INOUT)                   :: field_x(:,:), field_y(:,:)
!   TYPE(t_cartesian_coordinates), INTENT(IN) :: field_cc(:,:)
!   CALL cvec2gvec(field_cc(:,:)%x(1),&
!     &            field_cc(:,:)%x(2),&
!     &            field_cc(:,:)%x(3),&
!     &            cell_subset%patch%cells%center(:,:)%lon,&
!     &            cell_subset%patch%cells%center(:,:)%lat,&
!     &            field_x(:,:),        &
!     &            field_y(:,:))
! END SUBROUTINE update_from_cartesian_coords_for_cells

!<Optimize_Used>
  SUBROUTINE set_windstress(subset, mask, threshold, windstress, &
      &                     control, amplitude, zonal_waveno, meridional_waveno,center,length)
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER, INTENT(IN)              :: mask(:,:)
    INTEGER, INTENT(IN)              :: threshold
    REAL(wp), INTENT(INOUT)          :: windstress(:,:)
    INTEGER,  INTENT(IN)             :: control
    REAL(wp), INTENT(IN)             :: amplitude,zonal_waveno,meridional_waveno
    REAL(wp), INTENT(IN),OPTIONAL    :: center, length

    SELECT CASE (control)
    CASE (0) ! NO FORCING, SET TO ZERO ========================================
      windstress = 0.0_wp
    CASE (1:100)      ! FILE INPUT, DONE ELSEWHERE ============================
      CALL message('windstress forcing','file input')
    CASE (101:200)    ! ANALYTIC SETUP ========================================
      SELECT CASE (control)
      CASE(101)       ! constant amplitude
        windstress = amplitude
      CASE(102)       ! basin setup, zonally changed
        CALL basin_zonal(subset,mask,threshold,windstress,amplitude,length)
      CASE(103)       ! basin setup, meridionally changed
        CALL basin_meridional(subset,mask,threshold,windstress,amplitude,length)
      CASE(104)       ! zonally periodic, nonzero at pols, meridionally constant
        CALL zonal_periodic_nonzero_around_center_zero_at_pols(subset, mask, threshold, windstress, amplitude)
      CASE(105)
        CALL meridional_periodic_around_center_zero_at_pols(subset,mask,threshold,windstress, amplitude)
      CASE(106)       ! zonally periodic around a given center, zero at pols, meridionally constant
        CALL zonal_periodic_zero_at_pols(subset,mask,threshold,windstress,amplitude,zonal_waveno)
      CASE(107)       ! latteral cells, zonal period only
        CALL cells_zonal_periodic(subset,mask,threshold,windstress,amplitude,zonal_waveno)
      CASE(108)       ! latteral cells, zonally and meridionally periodic
        CALL cells_zonal_and_meridional_periodic(subset,mask,threshold,windstress,amplitude,zonal_waveno,meridional_waveno)
      CASE(109)
        CALL cells_zonal_and_meridional_periodic_constant_amplitude_sin(subset, mask, threshold, windstress, amplitude)
      CASE(110)
        CALL cells_zonal_and_meridional_periodic_constant_amplitude_cosin(subset, mask, threshold, windstress, amplitude)
      CASE(111)
        CALL Wolfe_Cessi_TestCase(subset, mask, threshold, windstress, amplitude)
      END SELECT
    END SELECT

  END SUBROUTINE set_windstress

  SUBROUTINE basin_zonal(subset, mask, threshold, field_2d, amplitude,length_opt, zonal_waveno_opt)
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER,  INTENT(IN)             :: mask(:,:)
    INTEGER, INTENT(IN)              :: threshold
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN) , OPTIONAL  :: length_opt,zonal_waveno_opt

    REAL(wp) :: length, zonal_waveno
    REAL(wp) :: lat(nproma,subset%patch%alloc_cell_blocks), lon(nproma,subset%patch%alloc_cell_blocks)

    length = basin_height_deg * deg2rad
    zonal_waveno = forcing_windstress_zonal_waveno

    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(zonal_waveno,zonal_waveno_opt)

    lat(:,:) = subset%patch%cells%center(:,:)%lat
    lon(:,:) = subset%patch%cells%center(:,:)%lon

    field_2d(:,:) = MERGE(amplitude * COS(zonal_waveno*pi*(lat(:,:)-length)/length),0.0_wp,mask(:,:) <= threshold)

  END SUBROUTINE basin_zonal

  SUBROUTINE basin_meridional(subset, mask, threshold, field_2d, amplitude,length_opt,meridional_waveno_opt)
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER, INTENT(IN)              :: mask(:,:)
    INTEGER, INTENT(IN)              :: threshold
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN) , OPTIONAL  :: length_opt,meridional_waveno_opt

    REAL(wp) :: length, meridional_waveno
    REAL(wp) :: lat(nproma,subset%patch%alloc_cell_blocks), lon(nproma,subset%patch%alloc_cell_blocks)

    length            = basin_width_deg * deg2rad
    meridional_waveno = forcing_windstress_merid_waveno

    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(meridional_waveno,meridional_waveno_opt)

    lat(:,:) = subset%patch%cells%center(:,:)%lat
    lon(:,:) = subset%patch%cells%center(:,:)%lon

    field_2d(:,:) = MERGE(amplitude * COS(meridional_waveno*pi*(lon(:,:)-length)/length),0.0_wp,mask(:,:) <= threshold)
  END SUBROUTINE basin_meridional

  SUBROUTINE zonal_periodic_nonzero_around_center_zero_at_pols(subset, mask, threshold, field_2d, amplitude,&
      & center_opt,length_opt,zonal_waveno_opt)
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER, INTENT(IN)              :: mask(:,:)
    INTEGER, INTENT(IN)              :: threshold
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: center_opt, length_opt, zonal_waveno_opt

    REAL(wp) :: center, length, zonal_waveno
    REAL(wp) :: lat(nproma,subset%patch%alloc_cell_blocks), lon(nproma,subset%patch%alloc_cell_blocks)

    length       = 180.0_wp * deg2rad
    center       = -60.0_wp * deg2rad
    zonal_waveno = forcing_windstress_zonal_waveno

    CALL assign_if_present(center,center_opt)
    CALL assign_if_present(length,length_opt)
    CaLL assign_if_present(zonal_waveno, zonal_waveno_opt)

    lat(:,:) = subset%patch%cells%center(:,:)%lat
    lon(:,:) = subset%patch%cells%center(:,:)%lon

    field_2d(:,:) = MERGE(amplitude * COS(zonal_waveno*pi*(lat(:,:)-center)/length),0.0_wp,mask(:,:) <= threshold)
  END SUBROUTINE zonal_periodic_nonzero_around_center_zero_at_pols

  SUBROUTINE meridional_periodic_around_center_zero_at_pols(subset, mask, threshold, field_2d, amplitude, &
      & center_opt,length_opt,meridional_waveno_opt)
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER, INTENT(IN)              :: mask(:,:)
    INTEGER, INTENT(IN)              :: threshold
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: center_opt, length_opt, meridional_waveno_opt

    REAL(wp) :: center, length, meridional_waveno
    REAL(wp) :: lat(nproma,subset%patch%alloc_cell_blocks), lon(nproma,subset%patch%alloc_cell_blocks)

    length                        =  90.0_wp * deg2rad
    center                        = -20.0_wp * deg2rad
    meridional_waveno             = forcing_windstress_merid_waveno

    CALL assign_if_present(center,center_opt)
    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(meridional_waveno, meridional_waveno_opt)

    lat(:,:) = subset%patch%cells%center(:,:)%lat
    lon(:,:) = subset%patch%cells%center(:,:)%lon

    field_2d(:,:) = MERGE(amplitude * COS(meridional_waveno*pi*(lon(:,:)-center)/length),0.0_wp,mask(:,:) <= threshold)
  END SUBROUTINE meridional_periodic_around_center_zero_at_pols

!<Optimize_Used>
  SUBROUTINE zonal_periodic_zero_at_pols(subset, mask, threshold, field_2d, amplitude, zonal_waveno_opt)
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER, INTENT(IN)              :: mask(:,:)
    INTEGER, INTENT(IN)              :: threshold
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: zonal_waveno_opt

    REAL(wp) :: zonal_waveno
    REAL(wp) :: lat(nproma,subset%patch%alloc_cell_blocks), lon(nproma,subset%patch%alloc_cell_blocks)

    zonal_waveno = forcing_windstress_zonal_waveno

    CaLL assign_if_present(zonal_waveno, zonal_waveno_opt)

    lat(:,:) = subset%patch%cells%center(:,:)%lat
    lon(:,:) = subset%patch%cells%center(:,:)%lon

#ifdef __SX__
    field_2d(:,:) = MERGE(amplitude * COS(lat(:,:)) * &
      & COS(forcing_windstress_zonalWavePhas + zonal_waveno * ABS(lat(:,:))), 0.0_wp, mask(:,:) <= threshold)
#else
    field_2d(:,:) = MERGE(amplitude * COS(lat(:,:)) * &
      & COS(forcing_windstress_zonalWavePhase + zonal_waveno * ABS(lat(:,:))), 0.0_wp, mask(:,:) <= threshold)
#endif

  END SUBROUTINE zonal_periodic_zero_at_pols

  SUBROUTINE Wolfe_Cessi_TestCase(subset, mask, threshold, field_2d, amplitude)
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER, INTENT(IN)              :: mask(:,:)
    INTEGER, INTENT(IN)              :: threshold
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude

    ! REAL(wp) :: zonal_waveno
    ! REAL(wp) :: reference_value
    REAL(wp) :: lat(nproma,subset%patch%alloc_cell_blocks), lon(nproma,subset%patch%alloc_cell_blocks)

    lat(:,:) = subset%patch%cells%center(:,:)%lat
    lon(:,:) = subset%patch%cells%center(:,:)%lon

    ! LL init:
  ! reference_value = COS(3.5_wp * pi)
  ! field_2d(:,:) = MERGE(amplitude *      &
  !   & (- ABS(SIN(4.0_wp * lat(:,:))) -   &
  !   & COS(7.0_wp * lat(:,:)) +           &
  !   & reference_value), 0.0_wp, mask(:,:) <= threshold)

  ! SLO init:
    field_2d(:,:) = MERGE(amplitude *                                   &
      & (0.8_wp*EXP(-1*lat(:,:)*lat(:,:)/0.01) - COS(6.0_wp*lat(:,:))), &
      & 0.0_wp, mask(:,:) <= threshold)

  END SUBROUTINE Wolfe_Cessi_TestCase

  SUBROUTINE cells_zonal_periodic(subset, mask, threshold, field_2d, amplitude, zonal_waveno_opt)
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER, INTENT(IN)              :: mask(:,:)
    INTEGER, INTENT(IN)              :: threshold
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: zonal_waveno_opt

    REAL(wp) :: zonal_waveno
    REAL(wp) :: lat(nproma,subset%patch%alloc_cell_blocks), lon(nproma,subset%patch%alloc_cell_blocks)

    zonal_waveno = forcing_windstress_zonal_waveno
    CALL assign_if_present(zonal_waveno, zonal_waveno_opt)

    lat(:,:) = subset%patch%cells%center(:,:)%lat
    lon(:,:) = subset%patch%cells%center(:,:)%lon

    field_2d(:,:) = MERGE(amplitude &
                & * COS(lat) &
                & * COS(zonal_waveno * lat) &
                & * COS(lon),0.0_wp,mask(:,:) <= threshold)
  END SUBROUTINE cells_zonal_periodic

  SUBROUTINE cells_zonal_and_meridional_periodic(subset, mask, threshold, field_2d, amplitude, &
      & zonal_waveno_opt,meridional_waveno_opt)
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER, INTENT(IN)              :: mask(:,:)
    INTEGER, INTENT(IN)              :: threshold
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: zonal_waveno_opt, meridional_waveno_opt

    REAL(wp) :: zonal_waveno, meridional_waveno
    REAL(wp) :: lat(nproma,subset%patch%alloc_cell_blocks), lon(nproma,subset%patch%alloc_cell_blocks)

    zonal_waveno      = forcing_windstress_zonal_waveno
    meridional_waveno = forcing_windstress_merid_waveno

    CALL assign_if_present(zonal_waveno, zonal_waveno_opt)
    CALL assign_if_present(meridional_waveno, meridional_waveno_opt)

    lat(:,:) = subset%patch%cells%center(:,:)%lat
    lon(:,:) = subset%patch%cells%center(:,:)%lon

    field_2d(:,:) = MERGE(amplitude &
      & * COS(lat(:,:)) &
      & * COS(zonal_waveno * lat(:,:)) &
      & * SIN(meridional_waveno * lon(:,:)),0.0_wp,mask(:,:) <= threshold)

  END SUBROUTINE cells_zonal_and_meridional_periodic

  SUBROUTINE cells_zonal_and_meridional_periodic_constant_amplitude_sin(subset, mask, threshold, field_2d, &
      & amplitude,length_opt, zonal_waveno_opt)
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER,  INTENT(IN)             :: mask(:,:)
    INTEGER, INTENT(IN)              :: threshold
    REAL(wp)                         :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN) , OPTIONAL  :: length_opt,zonal_waveno_opt

    REAL(wp) :: length, zonal_waveno
    REAL(wp) :: lat(nproma,subset%patch%alloc_cell_blocks), &
      &         lon(nproma,subset%patch%alloc_cell_blocks), &
      &         strength(nproma,subset%patch%alloc_cell_blocks)

    length       = basin_height_deg * deg2rad
    zonal_waveno = forcing_windstress_zonal_waveno

    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(zonal_waveno,zonal_waveno_opt)

    lat(:,:) = subset%patch%cells%center(:,:)%lat
    lon(:,:) = subset%patch%cells%center(:,:)%lon

    strength = MERGE(amplitude*cos(zonal_waveno*pi*lat(:,:)-length/length),0.0_wp, mask(:,:) <= threshold)

    field_2d(:,:) = amplitude*strength*sin(lon(:,:))

  END SUBROUTINE cells_zonal_and_meridional_periodic_constant_amplitude_sin
  SUBROUTINE cells_zonal_and_meridional_periodic_constant_amplitude_cosin(subset, mask, threshold, field_2d, &
      & amplitude,length_opt, zonal_waveno_opt)
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER,  INTENT(IN)             :: mask(:,:)
    INTEGER, INTENT(IN)              :: threshold
    REAL(wp)                         :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN) , OPTIONAL  :: length_opt,zonal_waveno_opt

    REAL(wp) :: length, zonal_waveno
    REAL(wp) :: lat(nproma,subset%patch%alloc_cell_blocks), &
      &         lon(nproma,subset%patch%alloc_cell_blocks), &
      &         strength(nproma,subset%patch%alloc_cell_blocks)

    length       = basin_height_deg * deg2rad
    zonal_waveno = forcing_windstress_zonal_waveno

    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(zonal_waveno,zonal_waveno_opt)

    lat(:,:) = subset%patch%cells%center(:,:)%lat
    lon(:,:) = subset%patch%cells%center(:,:)%lon

    strength = MERGE(amplitude*cos(zonal_waveno*pi*lat(:,:)-length/length),0.0_wp, mask(:,:) <= threshold)

    field_2d(:,:) = amplitude*strength*cos(lon(:,:))

  END SUBROUTINE cells_zonal_and_meridional_periodic_constant_amplitude_cosin

  SUBROUTINE nf(STATUS)

    INTEGER, INTENT(in) :: STATUS

    IF (STATUS /= nf_noerr) THEN
      CALL finish('mo_ext_data netCDF error', nf_strerror(STATUS))
    ENDIF

  END SUBROUTINE nf
  !-------------------------------------------------------------------------
END MODULE mo_oce_forcing
