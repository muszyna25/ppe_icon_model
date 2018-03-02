!>
!! Provide an implementation of the ocean surface module.
!!
!! Provide an implementation of the parameters used for surface forcing
!! of the hydrostatic ocean model.
!!
!! @author Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
!!  Restructured code by Stephan Lorenz, MPI-M: (2015-04)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ocean_surface
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2007
!
!-------------------------------------------------------------------------
!
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_sync,                ONLY: global_sum_array
  USE mo_io_units,            ONLY: filename_max
  USE mo_mpi,                 ONLY: my_process_is_stdio, p_io, p_bcast, p_comm_work_test, p_comm_work
  USE mo_parallel_config,     ONLY: p_test_run
  USE mo_read_interface,      ONLY: openInputFile, closeFile, t_stream_id, &
    &                               on_cells, read_2D_time  !, read_3D
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_ocean_ext_data,      ONLY: ext_data
  USE mo_ocean_nml,           ONLY: iforc_oce, forcing_timescale, relax_analytical_type,  &
    &  no_tracer, n_zlev, basin_center_lat,                  &
    &  basin_center_lon, basin_width_deg, basin_height_deg,  &
    &  para_surfRelax_Temp, type_surfRelax_Temp,             &
    &  para_surfRelax_Salt, type_surfRelax_Salt,             &
    &  No_Forcing, Analytical_Forcing, OMIP_FluxFromFile,    &
    &  Coupled_FluxFromAtmo,                                 &
    &  i_sea_ice, forcing_enable_freshwater, zero_freshwater_flux,    &
    &  forcing_set_runoff_to_zero,           &
    &  atmos_flux_analytical_type, atmos_precip_const, &  ! atmos_evap_constant
    &  atmos_SWnet_const, atmos_LWnet_const, atmos_lat_const, atmos_sens_const, &
    &  atmos_SWnetw_const, atmos_LWnetw_const, atmos_latw_const, atmos_sensw_const, &
    &  limit_elevation, l_relaxsal_ice, &
    &  relax_width, forcing_HeatFlux_amplitude, forcing_HeatFlux_base, &
    &  OceanReferenceDensity, &
    &  use_wind_mixing, nbgctra,lhamocc
  USE mo_dynamics_config,     ONLY: nold
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_dbg_nml,             ONLY: idbg_mxmn
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_math_constants,      ONLY: pi, deg2rad, rad2deg
  USE mo_physical_constants,  ONLY: rho_ref, alv, tmelt, tf, clw, albedoW_sim, stbo, zemiss_def
  USE mo_physical_constants,  ONLY: rd, cpd, fr_fac, alf  ! cd_ia, used for omip bulk formula
  USE mo_impl_constants,      ONLY: max_char_length, sea_boundary, MIN_DOLIC
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL
  USE mo_math_utilities,      ONLY: gvec2cvec
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_ocean_surface_types, ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_sea_ice,             ONLY: ice_fast, ice_clean_up_dyn, ice_ocean_stress
  USE mo_sea_ice_refactor,    ONLY: ice_slow_slo
  USE mo_ice_fem_interface,       ONLY: ice_fem_interface
  USE mo_ice_advection,       ONLY: ice_advection_upwind !, ice_advection_upwind_einar
  USE mo_sea_ice_nml,         ONLY: use_calculated_ocean_stress, i_ice_dyn
  USE mo_timer,               ONLY: timers_level, timer_start, timer_stop, timer_extra40
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mtime,                  ONLY: datetime, &
       &                            getDayOfYearFromDateTime, &
       &                            getNoOfDaysInYearDateTime, &
       &                            getNoOfSecondsElapsedInDayDateTime

  
  IMPLICIT NONE
  
  ! required for reading netcdf files
  INCLUDE 'netcdf.inc'

  PRIVATE

  ! Public interface
  PUBLIC  :: update_ocean_surface

  ! private routines
  PRIVATE :: update_flux_analytical
  PRIVATE :: update_flux_fromFile
  PRIVATE :: update_surface_relaxation
  PRIVATE :: apply_surface_relaxation
  PRIVATE :: calc_omip_budgets_ice
  PRIVATE :: calc_omip_budgets_oce
  PRIVATE :: read_forc_data_oce
  PRIVATE :: balance_elevation

  CHARACTER(len=12)           :: str_module    = 'oceanSurface'  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug
  REAL(wp), PARAMETER         :: seconds_per_month = 2.592e6_wp  ! TODO: use real month length

CONTAINS

  !-------------------------------------------------------------------------
  !
  !>
  !! Update ocean surface by applying flux forcing for hydrostatic ocean
  !!
  !!
  !! @par Revision History
  !! Initial release (mo_oce_bulk)        by Stephan Lorenz, MPI-M (2010-07)
  !! restructured code (mo_ocean_surface) by Stephan Lorenz, MPI-M (2015-04)
  !
!<Optimize_Used>
  SUBROUTINE update_ocean_surface(p_patch_3D, p_os, p_as, p_ice, atmos_fluxes, &
       &                          p_oce_sfc, jstep, this_datetime, p_op_coeff)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_atmos_for_ocean)                     :: p_as
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes
    TYPE(t_ocean_surface)                       :: p_oce_sfc      !  new forcing for ocean surface, replaces p_oce_sfc
    INTEGER, INTENT(IN)                         :: jstep
    TYPE(datetime), POINTER                     :: this_datetime
    TYPE(t_operator_coeff),   INTENT(IN)        :: p_op_coeff
    !
    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_surface:update_ocean_surface'
    INTEGER               :: jc, jb, trac_no
    INTEGER               :: i_startidx_c, i_endidx_c
    REAL(wp)              :: dsec
    REAL(wp)              :: Tfw(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: sss_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zUnderIceOld(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zUnderIceIni(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zUnderIceArt(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    
    REAL(wp) :: h_old_test

    TYPE(t_patch), POINTER:: p_patch 
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER:: i_bgc_tra
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------

    IF (no_tracer>=1) p_oce_sfc%sst => p_os%p_prog(nold(1))%tracer(:,1,:,1)
    IF (no_tracer>=2) p_oce_sfc%sss => p_os%p_prog(nold(1))%tracer(:,1,:,2)

    IF (iforc_oce == No_Forcing) RETURN  !  forcing for ocean not defined

    sss_inter(:,:)  = p_oce_sfc%sss(:,:)
    zUnderIceOld(:,:) = 0.0_wp
    zUnderIceIni(:,:) = 0.0_wp
    zUnderIceArt(:,:) = 0.0_wp

    


    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('on entry: hi     ',p_ice%hi       ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('on entry: hs     ',p_ice%hs       ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('on entry: concSum',p_ice%concSum  ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('on entry: SST    ',p_oce_sfc%sst  ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('on entry: TotalHeat', atmos_fluxes%HeatFlux_Total,     str_module,4,in_subset=p_patch%cells%owned)
    CALL dbg_print('ocesfc%windStr-u ',p_oce_sfc%topBC_windStress_u,       str_module,3,in_subset=p_patch%cells%owned)
    CALL dbg_print('sfcflx%windStr-u ',p_oce_sfc%TopBC_WindStress_u,str_module,3,in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  (1) Apply relaxation to surface temperature and salinity
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    IF (type_surfRelax_Temp >= 1) THEN
      trac_no = 1   !  tracer no 1: temperature
      CALL update_surface_relaxation(p_patch_3D, p_os, p_ice, p_oce_sfc, trac_no)

      !  apply restoring to surface temperature directly
      CALL apply_surface_relaxation(p_patch_3D, p_os, p_oce_sfc, trac_no)

    END IF

    IF (type_surfRelax_Salt >= 1 .AND. no_tracer >1) THEN
      trac_no = 2   !  tracer no 2: salinity
      CALL update_surface_relaxation(p_patch_3D, p_os, p_ice, p_oce_sfc, trac_no)

      !  apply restoring to surface salinity directly
      CALL apply_surface_relaxation(p_patch_3D, p_os, p_oce_sfc, trac_no)

    ENDIF

    ! assign freeboard before sea ice model
    IF (i_sea_ice==0) THEN
      p_ice%zUnderIce(:,:) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,1,:) + p_os%p_prog(nold(1))%h(:,:)
      p_oce_sfc%cellThicknessUnderIce(:,:) = p_ice%zUnderIce(:,:)
    ENDIF

    ! freeboard used for thermal boundary condition (Eq.1)
    zUnderIceIni(:,:) = p_ice%zUnderIce(:,:)

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  (2) Receive surface fluxes for ocean forcing
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    !   fluxes from analytical / OMIP and bulk / coupling
    SELECT CASE (iforc_oce)

    CASE (Analytical_Forcing)        !  11

      !  Driving the ocean with analytically calculated fluxes
      CALL update_flux_analytical(p_patch_3D, p_os, atmos_fluxes, p_oce_sfc)

      !  needed for old heat flux BC
      atmos_fluxes%HeatFlux_Total=p_oce_sfc%TopBC_Temp_vdiff

      ! provide dLWdt for ice_fast as for OMIP
      atmos_fluxes%dLWdT (:,:,:)  = -4._wp*zemiss_def*stbo*(p_ice%tsurf(:,:,:)+tmelt)**3

      ! provide constant water fluxes for special analytical cases
      IF (atmos_flux_analytical_type >= 101) THEN
        p_as%FrshFlux_Precipitation        (:,:) = atmos_precip_const
        atmos_fluxes%FrshFlux_Precipitation(:,:) = p_as%FrshFlux_Precipitation(:,:)
        atmos_fluxes%latw                  (:,:) = atmos_latw_const
      ENDIF

    CASE (OMIP_FluxFromFile)         !  12

      !  Driving the ocean with OMIP fluxes
      !   a) read OMIP data (read relaxation data)
      CALL update_flux_fromFile(p_patch_3D, p_as, jstep, this_datetime)
      !   b) calculate OMIP flux data

      ! bulk formula for heat flux are calculated globally using specific OMIP fluxes
      ! TODO: include calculation of icefree part of cells in calc_omip_budgets_ice

      CALL calc_omip_budgets_oce(p_patch, p_as, p_os, atmos_fluxes)

      ! pass parameters read from OMIP into budget calculation directly
      CALL calc_omip_budgets_ice(                            &
        &                        p_patch%cells%center(:,:)%lat,                  &  !  input parameter
        &                        p_as%tafo(:,:), p_as%ftdew(:,:)-tmelt,          &  !  input parameter
        &                        p_as%fu10(:,:), p_as%fclou(:,:), p_as%pao(:,:), &  !  input parameter
        &                        p_as%fswr(:,:), p_ice%kice, p_ice%Tsurf(:,:,:), &  !  input parameter
        &                        p_ice%hi(:,:,:),                                &  !  input parameter
        &                        atmos_fluxes%albvisdir(:,:,:),                  &  !  input parameter
        &                        atmos_fluxes%albvisdif(:,:,:),                  &  !  input parameter
        &                        atmos_fluxes%albnirdir(:,:,:),                  &  !  input parameter
        &                        atmos_fluxes%albnirdif(:,:,:),                  &  !  input parameter
        &                        atmos_fluxes%LWnet    (:,:,:),                  &  !  output parameter
        &                        atmos_fluxes%SWnet    (:,:,:),                  &  !  output parameter
        &                        atmos_fluxes%sens     (:,:,:),                  &  !  output parameter
        &                        atmos_fluxes%lat      (:,:,:),                  &  !  output parameter
        &                        atmos_fluxes%dLWdT    (:,:,:),                  &  !  output parameter
        &                        atmos_fluxes%dsensdT  (:,:,:),                  &  !  output parameter
        &                        atmos_fluxes%dlatdT   (:,:,:) )                    !  output parameter

      ! wind stress over ice is provided by OMIP data
      atmos_fluxes%stress_x(:,:) = p_as%topBoundCond_windStress_u(:,:)
      atmos_fluxes%stress_y(:,:) = p_as%topBoundCond_windStress_v(:,:)

      ! wind stress over water (stress_xw, stress_yw) is the same and read from OMIP, see calc_omip_budgets_oce

    CASE (Coupled_FluxFromAtmo)                                       !  14

      !  Driving the ocean in a coupled mode:
      !  nothing to be done, atmospheric fluxes are provided at the end of time stepping
      !  atmospheric fluxes drive the ocean; fluxes are calculated by atmospheric model
      !  use atmospheric fluxes directly, i.e. no bulk formula as for OMIP is applied
       
       ! HAMOCC uses p_as to get SW radiation and wind, so we need to copy
       ! the SW radiation onto it in the coupled case 
       if(lhamocc) p_as%fswr(:,:) = atmos_fluxes%HeatFlux_ShortWave(:,:)
      CONTINUE

    CASE DEFAULT

      CALL message(TRIM(routine), 'STOP: Ocean Forcing option not implemented' )
      CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION DOES NOT EXIST - TERMINATE')

    END SELECT

    ! copy atmospheric wind speed from p_as%fu10 into new forcing variable for output purpose - not accumulated yet
    p_oce_sfc%Wind_Speed_10m(:,:) = p_as%fu10(:,:)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('bef.fast: Tsurf  ',        p_ice%tsurf          ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.fast:atmflx%LWnetIce', atmos_fluxes%LWnet   ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.fast:atmflx%SensIce',  atmos_fluxes%sens    ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.fast:atmflx%LatentIce',atmos_fluxes%lat     ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.fast:atmflx%dsensdT'  ,atmos_fluxes%dsensdT ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.fast:atmflx%dlatdT'   ,atmos_fluxes%dlatdT  ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.fast:atmflx%dLWdT'    ,atmos_fluxes%dLWdt   ,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.fast:stress_x'        ,atmos_fluxes%stress_x,str_module,idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.fast:TotalHeat', atmos_fluxes%HeatFlux_Total,str_module,idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  (3) Calculate fast sea ice thermodynamics
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    !  for analytical and OMIP/bulk calculated fluxes only
    !  in coupled case ice_fast is called within atmosphere model
    IF (iforc_oce == Analytical_Forcing .OR. iforc_oce == OMIP_FluxFromFile)  THEN  !  11 or 12

      IF (i_sea_ice >0 ) THEN

      ! Calculate the sea surface freezing temperature
      !  2015-06: array used in ice_fast, set to constant Tf
      !  if Tfw is variable it should be included in ice-variables and initialized in ice_init
      Tfw(:,:) = Tf

        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
          CALL ice_fast(i_startidx_c, i_endidx_c, nproma, p_ice%kice, dtime, &
            &   p_ice% Tsurf(:,:,jb),   &          !  intent(inout)
            &   p_ice% T1   (:,:,jb),   &          !  intent(out)   dummy for zerolayer model
            &   p_ice% T2   (:,:,jb),   &          !  intent(out)   dummy for zerolayer model
            &   p_ice% hi   (:,:,jb),   &          !  intent(in)
            &   p_ice% hs   (:,:,jb),   &          !  intent(in)
            &   p_ice% Qtop (:,:,jb),   &          !  intent(out)
            &   p_ice% Qbot (:,:,jb),   &          !  intent(out)
            &   atmos_fluxes%SWnet  (:,:,jb),   &  !  following: intent(in)
            &   atmos_fluxes%lat(:,:,jb) + atmos_fluxes%sens(:,:,jb) + atmos_fluxes%LWnet(:,:,jb),   & 
            &   atmos_fluxes%dlatdT(:,:,jb) + atmos_fluxes%dsensdT(:,:,jb) + atmos_fluxes%dLWdT(:,:,jb),   & 
            &   Tfw         (:,  jb),   &
            &   atmos_fluxes%albvisdir(:,:,jb), &  !  albedos: intent(out)
            &   atmos_fluxes%albvisdif(:,:,jb), &
            &   atmos_fluxes%albnirdir(:,:,jb), &
            &   atmos_fluxes%albnirdif(:,:,jb), &
            &   doy = getDayOfYearFromDateTime(this_datetime))
        ENDDO
       
        ! Unique albedo for analytical and OMIP cases (i_ice_albedo=1)
        !  the near infrared and diffuse albedos may be used for calculation of tsurf
        atmos_fluxes%albvisdirw = albedoW_sim
        atmos_fluxes%albvisdifw = albedoW_sim
        atmos_fluxes%albnirdirw = albedoW_sim
        atmos_fluxes%albnirdifw = albedoW_sim
       
        ! provide constant heat fluxes for special analytical cases
        IF (atmos_flux_analytical_type == 102) THEN
          p_ice%Qtop(:,1,:) = atmos_SWnet_const
          p_ice%Qbot(:,1,:) = 0.0_wp
        ENDIF
        IF (atmos_flux_analytical_type == 103) THEN
          p_ice%Qtop(:,1,:) = 0.0_wp
          p_ice%Qbot(:,1,:) = atmos_sens_const
        ENDIF

      ENDIF  !  sea ice

    ENDIF  !  iforc_oce = analytical or OMIP, 11 or 12

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('aft.fast: Tsurf  ',p_ice%tsurf    ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.fast: Qtop   ',p_ice%Qtop     ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.fast: Qbot   ',p_ice%Qbot     ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.fast: dLWdT  ',atmos_fluxes%dLWdT,     str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.fast: albvdir',atmos_fluxes%albvisdir ,str_module,3, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  (4a) Provide fluxes for slow sea ice thermodynamics
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    IF (iforc_oce == Analytical_Forcing .OR. iforc_oce == OMIP_FluxFromFile)  THEN  !  11 or 12

      ! provide evaporation from latent heat flux for OMIP case
      ! under sea ice evaporation is neglected, atmos_fluxes%latw is flux in the absence of sea ice
      atmos_fluxes%FrshFlux_Evaporation(:,:) = atmos_fluxes%latw(:,:) / (alv*rho_ref)

      !  copy variables into atmos_fluxes
      atmos_fluxes%FrshFlux_Runoff(:,:)      = p_as%FrshFlux_Runoff(:,:)
    
      ! Precipitation on ice is snow when tsurf is below the freezing point
      !  - no snowfall from OMIP data
      !  - rprecw, rpreci are water equivalent over whole grid-area
      WHERE ( ALL( p_ice%Tsurf(:,:,:) < 0._wp, 2 ) )  !  Tsurf is -1.8 over open water, incorrect specification
        atmos_fluxes%rpreci(:,:) = p_as%FrshFlux_Precipitation(:,:)
        atmos_fluxes%rprecw(:,:) = 0._wp
      ELSEWHERE
        ! not considered in ice_growth_zero
        atmos_fluxes%rpreci(:,:) = 0._wp
        atmos_fluxes%rprecw(:,:) = p_as%FrshFlux_Precipitation(:,:)
      ENDWHERE

      ! evaporation and runoff not used in sea ice but in VolumeTotal, evaporation used for TotalOcean only
      atmos_fluxes%FrshFlux_TotalOcean(:,:) = p_patch_3d%wet_c(:,1,:)*( 1.0_wp-p_ice%concSum(:,:) ) * &
        &  (p_as%FrshFlux_Precipitation(:,:) + atmos_fluxes%FrshFlux_Evaporation(:,:))

    ELSEIF (iforc_oce == Coupled_FluxFromAtmo)  THEN

      ! these 4 fluxes over open ocean are used in sea ice thermodynamics
      atmos_fluxes%SWnetw (:,:)   = atmos_fluxes%HeatFlux_ShortWave(:,:)
      atmos_fluxes%LWnetw (:,:)   = atmos_fluxes%HeatFlux_LongWave (:,:)
      atmos_fluxes%sensw  (:,:)   = atmos_fluxes%HeatFlux_Sensible (:,:)
      atmos_fluxes%latw   (:,:)   = atmos_fluxes%HeatFlux_Latent   (:,:)
     
      WHERE ( p_ice%concSum(:,:) > 0._wp) !  corresponding to (1-concSum)*Precip in TotalOcean
   !  WHERE ( ALL( p_ice%hi   (:,:,:) > 0._wp, 2 ) )  !  corresponding to hi>0 in ice_growth_zero
   !  WHERE ( ALL( p_ice%Tsurf(:,:,:) < 0._wp, 2 ) )  !  Tsurf is -1.8 over open water, incorrect specification
        ! SnowFall and liquid rain over ice-covered part of ocean are taken from the atmosphere model
        atmos_fluxes%rpreci(:,:) = atmos_fluxes%FrshFlux_SnowFall(:,:)
        atmos_fluxes%rprecw(:,:) = atmos_fluxes%FrshFlux_Precipitation(:,:) - atmos_fluxes%FrshFlux_SnowFall(:,:)
      ELSEWHERE
        ! not considered in ice_growth_zero
        atmos_fluxes%rpreci(:,:) = 0._wp
        atmos_fluxes%rprecw(:,:) = atmos_fluxes%FrshFlux_Precipitation(:,:)
      ENDWHERE

      ! copy flux for use in TotalOcean, since analytical/omip use p_as:
      !p_as%FrshFlux_Precipitation      = atmos_fluxes%FrshFlux_Precipitation

      ! total water flux over ice-free ocean water: P*(1-C)+E
      !  - whole evaporation over grid-box enters open ocean, this includes evaporation over sea ice covered part
      !  - snowfall is included as (melted) water equivalent
      !  - runoff is added to VolumeTotal below
      atmos_fluxes%FrshFlux_TotalOcean(:,:) = p_patch_3d%wet_c(:,1,:)* &
        &  (( 1.0_wp-p_ice%concSum(:,:) ) * atmos_fluxes%FrshFlux_Precipitation(:,:) + atmos_fluxes%FrshFlux_Evaporation(:,:))

    ENDIF  ! iforc_oce

    IF (zero_freshwater_flux) THEN
      ! since latw<>0. we must set evap and TotalOcean again to zero:
      atmos_fluxes%FrshFlux_Evaporation  (:,:) = 0.0_wp
      atmos_fluxes%FrshFlux_TotalOcean   (:,:) = 0.0_wp
      atmos_fluxes%FrshFlux_Precipitation(:,:) = 0.0_wp
      atmos_fluxes%FrshFlux_SnowFall     (:,:) = 0.0_wp
      atmos_fluxes%FrshFlux_Evaporation  (:,:) = 0.0_wp
      atmos_fluxes%FrshFlux_Runoff       (:,:) = 0.0_wp
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('aftAtmFB: Precipitation', atmos_fluxes%FrshFlux_Precipitation,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFB: Evaporation'  , atmos_fluxes%FrshFlux_Evaporation  ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFB: SnowFall'     , atmos_fluxes%FrshFlux_SnowFall     ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFB: Runoff'       , atmos_fluxes%FrshFlux_Runoff       ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFB: TotalOcean'   , atmos_fluxes%FrshFlux_TotalOcean   ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFB: rprecw'       , atmos_fluxes%rprecw                ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFB: rpreci'       , atmos_fluxes%rpreci                ,str_module, 3, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  (4b) Call sea ice dynamics
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    ! ocean stress calculated independent of ice dynamics
    CALL ice_ocean_stress( p_patch, atmos_fluxes, p_ice, p_os )

    CALL dbg_print('bef.icedyn: hi   ',p_ice%hi,       str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.icedyn: hs   ',p_ice%hs,       str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.icedyn: Conc.',p_ice%conc     ,str_module, 3, in_subset=p_patch%cells%owned)

    IF ( i_ice_dyn >= 1 ) THEN
      ! AWI FEM model wrapper
      IF (timers_level > 1) CALL timer_start(timer_extra40)
      CALL ice_fem_interface ( p_patch_3D, p_ice, p_os, atmos_fluxes, p_op_coeff)
!      CALL ice_advection_upwind_einar( p_patch_3D, p_op_coeff, p_ice ) ! messy advection routine, bugs fixed; renamed as ice_advection_upwind
      CALL ice_advection_upwind( p_patch_3D, p_op_coeff, p_ice )

      ! the original clean up routine has been split into two: ice_clean_up_dyn, ice_clean_up_thd
      ! here we fix possible overshoots in conc afther the advection step
      CALL ice_clean_up_dyn( p_patch_3D, p_ice )
      IF (timers_level > 1) CALL timer_stop(timer_extra40)

    ELSE
      p_ice%u = 0._wp
      p_ice%v = 0._wp
    ENDIF

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  (4c) Call slow sea ice thermodynamics
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('bef.slow: hi     ',p_ice%hi       ,str_module,4, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.slow: Tsurf  ',p_ice%tsurf    ,str_module,4, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.slow: SST    ',p_os%p_prog(nold(1))%tracer(:,1,:,1), str_module, 4, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    IF (i_sea_ice >= 1) THEN

      ! call to refactored ice thermodynamics
      CALL ice_slow_slo(p_patch_3D, p_os, p_ice, atmos_fluxes, p_op_coeff)

      ! ice_clean_up_thd routine is called inside ice_slow_slo
      ! it fixes undershoots in concentation;
      ! limits sea ice thickness to seaice_limit of surface layer depth after changes due to the thermodynamic growth/melt;
      ! calculates the new freeboard (used below at step (6))

    ELSE   !  no sea ice
     
      ! apply wind stress to forcing variable since no ice_ocean_stress routine is called
      atmos_fluxes%topBoundCond_windStress_u(:,:) = atmos_fluxes%stress_xw(:,:)
      atmos_fluxes%topBoundCond_windStress_v(:,:) = atmos_fluxes%stress_yw(:,:)

      ! apply net surface heat flux in W/m2 for OMIP case, since these fluxes are calculated in calc_omip_budgets_oce
      IF (iforc_oce == OMIP_FluxFromFile) THEN
        WHERE (p_patch_3D%lsm_c(:,1,:) <= sea_boundary)
          atmos_fluxes%HeatFlux_ShortWave(:,:) = atmos_fluxes%SWnetw(:,:) ! net SW radiation flux over water
          atmos_fluxes%HeatFlux_LongWave (:,:) = atmos_fluxes%LWnetw(:,:) ! net LW radiation flux over water
          atmos_fluxes%HeatFlux_Sensible (:,:) = atmos_fluxes%sensw (:,:) ! Sensible heat flux over water
          atmos_fluxes%HeatFlux_Latent   (:,:) = atmos_fluxes%latw  (:,:) ! Latent heat flux over water
        ELSEWHERE
          atmos_fluxes%HeatFlux_ShortWave(:,:) = 0.0_wp
          atmos_fluxes%HeatFlux_LongWave (:,:) = 0.0_wp
          atmos_fluxes%HeatFlux_Sensible (:,:) = 0.0_wp
          atmos_fluxes%HeatFlux_Latent   (:,:) = 0.0_wp
        ENDWHERE
      ENDIF  ! OMIP case
     
      ! sum of ocean heat fluxes for ocean boundary condition without ice, generally aggregated in ice thermodynamics
      atmos_fluxes%HeatFlux_Total(:,:) = atmos_fluxes%HeatFlux_ShortWave(:,:) + atmos_fluxes%HeatFlux_LongWave(:,:) &
        &                              + atmos_fluxes%HeatFlux_Sensible(:,:)  + atmos_fluxes%HeatFlux_Latent(:,:)
     
      ! for the setup without sea ice the SST is set to freezing temperature Tf
      WHERE (p_oce_sfc%SST(:,:) .LT. Tf)
        p_oce_sfc%SST(:,:) = Tf
      ENDWHERE

    ENDIF  !  sea ice

    ! provide total salinity forcing flux for diagnostics only - salt_content_in_surface
    atmos_fluxes%FrshFlux_TotalSalt(:,:) = atmos_fluxes%FrshFlux_Runoff(:,:)     &
      &                                  + atmos_fluxes%FrshFlux_TotalIce(:,:)   &
      &                                  + atmos_fluxes%FrshFlux_TotalOcean(:,:)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('aft.slow: hi     ',   p_ice%hi       ,str_module,                    4, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.slow: concSum',   p_ice%concSum  ,str_module,                    4, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.slow: TotalSalt', atmos_fluxes%FrshFlux_TotalSalt,   str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.slow: TotalOcean',atmos_fluxes%FrshFlux_TotalOcean,  str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.slow: TotalHeat', atmos_fluxes%HeatFlux_Total,       str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.slow: tracer1',p_os%p_prog(nold(1))%tracer(:,1,:,1), str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.slow: sfc%SST',   p_oce_sfc%SST,                     str_module, 3, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  (5) Set wind stress boundary condition
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    ! windstress
    p_oce_sfc%TopBC_WindStress_u(:,:) = atmos_fluxes%topBoundCond_windStress_u(:,:)
    p_oce_sfc%TopBC_WindStress_v(:,:) = atmos_fluxes%topBoundCond_windStress_v(:,:)

    ! After final updating of zonal and merdional components (from file, bulk formula, or coupling)
    ! cartesian coordinates are calculated
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
          CALL gvec2cvec(  p_oce_sfc%TopBC_WindStress_u(jc,jb),&
                         & p_oce_sfc%TopBC_WindStress_v(jc,jb),&
                         & p_patch%cells%center(jc,jb)%lon,&
                         & p_patch%cells%center(jc,jb)%lat,&
                         & p_oce_sfc%TopBC_WindStress_cc(jc,jb)%x(1),&
                         & p_oce_sfc%TopBC_WindStress_cc(jc,jb)%x(2),&
                         & p_oce_sfc%TopBC_WindStress_cc(jc,jb)%x(3))
        ELSE
          p_oce_sfc%TopBC_WindStress_u(jc,jb)         = 0.0_wp
          p_oce_sfc%TopBC_WindStress_v(jc,jb)         = 0.0_wp
          p_oce_sfc%TopBC_WindStress_cc(jc,jb)%x      = 0.0_wp
        ENDIF
      END DO
    END DO

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('UpdSfc: aft.Bulk/Ice: hi' , p_ice%hi                 , str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: aft.Bulk/Ice: hs' , p_ice%hs                 , str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: aft.Bulk/Ice:conc', p_ice%conc               , str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('sfc_flx: windStress_u',p_oce_sfc%TopBC_WindStress_u, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('sfc_flx: windStress_v',p_oce_sfc%TopBC_WindStress_v, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('sfc_flx: windStress_cc1',p_oce_sfc%TopBC_WindStress_cc%x(1), &
      &             str_module,3, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

 !  ! Prepare windstress boundary condition
 !  !  - on new surface type - not yet, since p_oce_sfc argument must be passed through some routines
 !  !  - needs allocation of cartesian coordinate variable
 !  p_oce_sfc%TopBC_windStress_u(:,:) = atmos_fluxes%topBoundCond_windStress_u(:,:)
 !  p_oce_sfc%TopBC_windStress_v(:,:) = atmos_fluxes%topBoundCond_windStress_v(:,:)

 !  !
 !  ! After final updating of zonal and merdional components (from file, bulk formula, or coupling)
 !  ! cartesian coordinates are calculated
 !  !
 !  DO jb = all_cells%start_block, all_cells%end_block
 !    CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
 !    DO jc = i_startidx_c, i_endidx_c
 !      IF(p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
 !        CALL gvec2cvec(  p_oce_sfc%TopBC_windStress_u(jc,jb),&
 !                       & p_oce_sfc%TopBC_windStress_v(jc,jb),&
 !                       & p_patch%cells%center(jc,jb)%lon,&
 !                       & p_patch%cells%center(jc,jb)%lat,&
 !                       & p_oce_sfc%TopBC_windStress_cc(jc,jb)%x(1),&
 !                       & p_oce_sfc%TopBC_windStress_cc(jc,jb)%x(2),&
 !                       & p_oce_sfc%TopBC_windStress_cc(jc,jb)%x(3))
 !      ENDIF
 !    END DO
 !  END DO


    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  (6) Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    !!  Provide total ocean forcing:
    !    - total heat fluxes are aggregated for ice/ocean in ice thermodynamics
    p_oce_sfc%HeatFlux_Total(:,:)       = atmos_fluxes%HeatFlux_Total(:,:)
    !    - shortwave heat flux for calculation of penetration depth
    p_oce_sfc%HeatFlux_Shortwave(:,:)   = atmos_fluxes%HeatFlux_Shortwave(:,:)
    !    - total internal salt flux atmos_fluxes%FrshFlux_TotalIce is calculated in sea ice model
    p_oce_sfc%FrshFlux_TotalIce(:,:)    = atmos_fluxes%FrshFlux_TotalIce(:,:)
    !    - total freshwater volume forcing
    p_oce_sfc%FrshFlux_VolumeTotal(:,:) = atmos_fluxes%FrshFlux_Runoff    (:,:) &
      &                                 + atmos_fluxes%FrshFlux_VolumeIce (:,:) &
      &                                 + atmos_fluxes%FrshFlux_TotalOcean(:,:) &
      &                                 + p_oce_sfc%FrshFlux_Relax     (:,:)

    !  ******  (Thermodynamic Eq. 1)  ******
    ! Apply net surface heat flux to ocean surface (new p_oce_flx%SST)
    IF (no_tracer > 0) THEN

      ! total heat flux calculated in sea ice module
      p_oce_sfc%HeatFlux_Total(:,:) = atmos_fluxes%HeatFlux_Total(:,:)

      ! sst-change in surface module after sea-ice thermodynamics using HeatFlux_Total and old freeboard zUnderIceIni
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN
            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    p_oce_sfc%HeatFlux_Total(jc,jb)*dtime/(clw*OceanReferenceDensity*zUnderIceIni(jc,jb))
          ENDIF
        ENDDO
      ENDDO

    END IF
    
    CALL dbg_print('UpdSfc: h-old',p_os%p_prog(nold(1))%h,         str_module, 1, in_subset=p_patch%cells%owned)
    ! apply volume flux to surface elevation
    !  - add to h_old before explicit term
    !  - change in salt concentration applied here
    !    i.e. for salinity relaxation only, no volume flux is applied
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      
      DO jc = i_startidx_c, i_endidx_c
        IF (p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0) THEN

          !******  (Thermodynamic Eq. 2)  ******
          !! Calculate the new freeboard caused by changes in ice thermodynamics
          !!  zUnderIce = z_surf + h_old - (z_draft - z_snowfall)
          !!  - new zUnderIce is calculated in thermodynamics in routine ice_clean_up_thd

          !******  (Thermodynamic Eq. 3)  ******
          !! First, calculate internal salinity change caused by melting of snow and melt or growth of ice:
          !!   SSS_new * zUnderIce = SSS_old * zUnderIceArt
          !!   artificial freeboard zUnderIceArt is used for internal Salinity change only:
          !!   - melt/growth of ice and snow to ice conversion imply a reduced water flux compared to saltfree water
          !!   - reduced water flux is calculated in FrshFlux_TotalIce by the term  (1-Sice/SSS)
          !!   - respective zUnderIceArt for calculating salt change is derived from these fluxes
          !!     which are calculated in sea ice thermodynamics (upper_ocean_TS)
          !    - for i_sea_ice=0 it is FrshFlux_TotalIce=0 and no change here
          zUnderIceArt(jc,jb)= p_ice%zUnderIce(jc,jb) - p_oce_sfc%FrshFlux_TotalIce(jc,jb)*dtime
          sss_inter(jc,jb)   = p_oce_sfc%sss(jc,jb) * zUnderIceArt(jc,jb) / p_ice%zUnderIce(jc,jb)
        
          

          !******  (Thermodynamic Eq. 4)  ******
          !! Next, calculate salinity change caused by rain and runoff without snowfall by adding their freshwater to zUnderIce
          zUnderIceOld(jc,jb)    = p_ice%zUnderIce(jc,jb)
          p_ice%zUnderIce(jc,jb) = zUnderIceOld(jc,jb) + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb) * dtime
          p_oce_sfc%SSS(jc,jb)   = sss_inter(jc,jb) * zUnderIceOld(jc,jb) / p_ice%zUnderIce(jc,jb)

          h_old_test=  (p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,1,jb)+p_os%p_prog(nold(1))%h(jc,jb))


          !******  (Thermodynamic Eq. 5)  ******
          !! Finally, let sea-level change from P-E+RO plus snow fall on ice, net total volume forcing to ocean surface
          p_os%p_prog(nold(1))%h(jc,jb) = p_os%p_prog(nold(1))%h(jc,jb)               &
            &                           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
            &                           + p_ice%totalsnowfall(jc,jb)
         
          !! update zunderice
          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) + p_os%p_prog(nold(1))%h(jc,jb) &
            &                    - p_ice%draftave(jc,jb)
    
 
          if(lhamocc.and. (p_os%p_prog(nold(1))%h(jc,jb)+ p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,1,jb)) > 0._wp)then 
          DO i_bgc_tra = no_tracer+1, no_tracer+nbgctra
           ! for HAMOCC tracer dilution
             p_os%p_prog(nold(1))%tracer(jc,1,jb,i_bgc_tra)  = p_os%p_prog(nold(1))%tracer(jc,1,jb,i_bgc_tra) &
           &        * h_old_test/(p_os%p_prog(nold(1))%h(jc,jb) + p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,1,jb))
          ENDDO
          endif

          

        ENDIF  !  dolic>0
      END DO
    END DO
         
    !! set correct cell thickness under ice
    p_oce_sfc%cellThicknessUnderIce   (:,:) = p_ice%zUnderIce(:,:)
!    atmos_fluxes%cellThicknessUnderIce(:,:) = p_ice%zUnderIce(:,:)  ! for diagnosis only
      
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('UpdSfc: oce_sfc%HFTot ', p_oce_sfc%HeatFlux_Total,       str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: oce_sfc%VolTot', p_oce_sfc%FrshFlux_VolumeTotal, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: oce_sfc%TotIce', p_oce_sfc%FrshFlux_TotalIce,    str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIceIni',   zUnderIceIni,                   str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIceArt',   zUnderIceArt,                   str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIceOld',   zUnderIceOld,                   str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIce   ',   p_ice%zUnderIce,                str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: sss_inter   ',   sss_inter,                      str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST,                  str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SSS ',p_oce_sfc%SSS,                  str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEnd: h-old+fwfVol',p_os%p_prog(nold(1))%h,         str_module, 2, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    ! copy fluxes to bulk-type variables for output and average statistic purposes only:
    p_oce_sfc%FrshFlux_Precipitation = atmos_fluxes%FrshFlux_Precipitation
    p_oce_sfc%FrshFlux_Evaporation   = atmos_fluxes%FrshFlux_Evaporation
    p_oce_sfc%FrshFlux_SnowFall      = atmos_fluxes%FrshFlux_SnowFall
    p_oce_sfc%FrshFlux_Runoff        = atmos_fluxes%FrshFlux_Runoff
    p_oce_sfc%HeatFlux_Total         = atmos_fluxes%HeatFlux_Total
    p_oce_sfc%HeatFlux_ShortWave     = atmos_fluxes%HeatFlux_ShortWave
    p_oce_sfc%HeatFlux_Longwave      = atmos_fluxes%HeatFlux_Longwave
    p_oce_sfc%HeatFlux_Sensible      = atmos_fluxes%HeatFlux_Sensible
    p_oce_sfc%HeatFlux_Latent        = atmos_fluxes%HeatFlux_Latent

    p_oce_sfc%FrshFlux_TotalOcean    = atmos_fluxes%FrshFlux_TotalOcean
    p_oce_sfc%FrshFlux_TotalSalt     = atmos_fluxes%FrshFlux_TotalSalt
!    p_oce_sfc%FrshFlux_TotalIce      = p_oce_sfc%FrshFlux_TotalIce
!    p_oce_sfc%FrshFlux_VolumeTotal   = p_oce_sfc%FrshFlux_VolumeTotal

    
    ! apply volume flux correction: 
    !  - sea level is balanced to zero over ocean surface
    !  - correction applied daily
    !  calculate time
    dsec  = REAL(getNoOfSecondsElapsedInDayDateTime(this_datetime), wp)
    ! event at end of first timestep of day - tbd: use mtime
    IF (limit_elevation .AND. (dsec-dtime)<0.1 ) THEN
      CALL balance_elevation(p_patch_3D, p_os%p_prog(nold(1))%h)
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('UpdSfc: h-old+BalElev',p_os%p_prog(nold(1))%h  ,str_module, 2, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------
    END IF

  END SUBROUTINE update_ocean_surface

  !-------------------------------------------------------------------------
  !
  !>
  !! Update surface flux forcing from file
  !!
  !! Provides surface forcing fluxes for ocean model from file.
  !!  Reads OMIP/NCEP fluxes via netcdf for bulk formula
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010/2014)
  !
  SUBROUTINE update_flux_fromFile(p_patch_3D, p_as, jstep, this_datetime)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_atmos_for_ocean)                     :: p_as
    INTEGER, INTENT(IN)                         :: jstep
    TYPE(datetime), POINTER                     :: this_datetime
    !
    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_surface:update_flux_fromFile'
    INTEGER  :: jmon, jdmon, jmon1, jmon2, ylen, yday
    REAL(wp) :: rday1, rday2
    REAL(wp) ::  z_c2(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    TYPE(t_patch), POINTER:: p_patch 
    !TYPE(t_subset_range), POINTER :: all_cells

    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------

    !all_cells       => p_patch%cells%all

    !  calculate day and month
    jmon  = this_datetime%date%month
    jdmon = this_datetime%date%day
    yday  = getDayOfYearFromDateTime(this_datetime)
    ylen  = getNoOfDaysInYearDateTime(this_datetime)

    !
    ! use annual forcing-data:
    !
    IF (forcing_timescale == 1)  THEN

      jmon1=1
      jmon2=1
      rday1=0.5_wp
      rday2=0.5_wp

    !
    ! interpolate monthly forcing-data daily:
    !
    ELSE IF (forcing_timescale == 12)  THEN

      jmon1=jmon-1
      jmon2=jmon
      rday1=REAL(15-jdmon,wp)/30.0_wp
      rday2=REAL(15+jdmon,wp)/30.0_wp
      IF (jdmon > 15)  THEN
        jmon1=jmon
        jmon2=jmon+1
        rday1=REAL(45-jdmon,wp)/30.0_wp
        rday2=REAL(jdmon-15,wp)/30.0_wp
      END IF

      IF (jmon1 ==  0) jmon1=12
      IF (jmon1 == 13) jmon1=1
      IF (jmon2 ==  0) jmon2=12
      IF (jmon2 == 13) jmon2=1

    !
    ! apply daily forcing-data directly:
    !
    ELSE

      ! - now daily data sets are read in mo_ext_data
      ! - use rday1, rday2, jmon1 = jmon2 = yday for controling correct day in year
      ! - no interpolation applied, 
      jmon1 = yday
      jmon2 = jmon1
      rday1 = 1.0_wp
      rday2 = 0.0_wp

      ! Leap year: read Feb, 28 twice since only 365 data-sets are available
      IF (ylen == 366) then
        IF (yday>59) jmon1=yday-1
        jmon2=jmon1
      ENDIF

    END IF

    !
    ! OMIP data read in mo_ext_data into variable ext_data
    !

    ! file based wind forcing:
    ! provide OMIP fluxes for wind stress forcing
    ! data set 1:  wind_u(:,:)   !  'stress_x': zonal wind stress       [Pa]
    ! data set 2:  wind_v(:,:)   !  'stress_y': meridional wind stress  [Pa]
    !  - forcing_windstress_u_type and v_type not used anymore
    !  - full OMIP data read if iforc_oce=OMIP_FluxFromFile (=11)

    ! ext_data has rank n_dom due to grid refinement in the atmosphere but not in the ocean
    !IF (forcing_windstress_u_type == 1)
    p_as%topBoundCond_windStress_u(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,1) + &
      &                                   rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,1)

    !IF (forcing_windstress_v_type == 1) THEN
    p_as%topBoundCond_windStress_v(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,2) + &
      &                                   rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,2)


    !-------------------------------------------------------------------------
    ! provide OMIP fluxes for sea ice (interface to ocean)
    ! data set 4:  tafo(:,:),   &  ! 2 m air temperature                              [C]
    ! data set 5:  ftdew(:,:),  &  ! 2 m dew-point temperature                        [K]
    ! data set 6:  fu10(:,:) ,  &  ! 10 m wind speed                                  [m/s]
    ! data set 7:  fclou(:,:),  &  ! Fractional cloud cover
    ! data set 8:  pao(:,:),    &  ! Surface atmospheric pressure                     [hPa]
    ! data set 9:  fswr(:,:),   &  ! Incoming surface solar radiation                 [W/m]
    ! data set 10:  precip(:,:), &  ! precipitation rate                              [m/s]
    ! data set 11:  evap  (:,:), &  ! evaporation   rate                              [m/s]
    ! data set 12:  runoff(:,:)     ! river runoff  rate                              [m/s]
    ! data set 13: u(:,:),      &  ! 10m zonal wind speed                             [m/s]
    ! data set 14: v(:,:),      &  ! 10m meridional wind speed                        [m/s]

    !IF (iforc_type == 2 .OR. iforc_type == 5) THEN
    !IF (forcing_fluxes_type > 0 .AND. forcing_fluxes_type < 101 ) THEN
    !  - forcing_fluxes_type = 1 not used anymore,
    !  - full OMIP data read if iforc_oce=OMIP_FluxFromFile (=11)

    p_as%tafo(:,:)  = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,4) + &
      &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,4)
    !  - change units to deg C, subtract tmelt (0 deg C, 273.15)
    p_as%tafo(:,:)  = p_as%tafo(:,:) - tmelt
    p_as%ftdew(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,5) + &
      &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,5)
    p_as%fu10(:,:)  = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,6) + &
      &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,6)
    p_as%fclou(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,7) + &
      &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,7)
    p_as%pao(:,:)   = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,8) + &
      &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,8)
    !  don't - change units to mb/hPa
    !p_as%pao(:,:)   = p_as%pao(:,:) !* 0.01
    p_as%fswr(:,:)  = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,9) + &
      &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,9)
      
    p_as%u(:,:)     = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,13) + &
      &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,13)
    p_as%v(:,:)     = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,14) + &
      &               rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,14)

    ! provide precipitation, evaporation, runoff flux data for freshwater forcing of ocean 
    !  - not changed via bulk formula, stored in surface flux data
    !  - Attention: as in MPIOM evaporation is calculated from latent heat flux (which is depentent on current SST)
    !               therefore not applied here
    p_as%FrshFlux_Precipitation(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,10) + &
      &                                     rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,10)
    !p_as%FrshFlux_Evaporation  (:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,11) + &
    !  &                                     rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,11)
    IF (forcing_set_runoff_to_zero) THEN
      p_as%FrshFlux_Runoff(:,:) = 0.0_wp
    ELSE
      p_as%FrshFlux_Runoff(:,:) = rday1*ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,12) + &
        &                              rday2*ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,12)
    ENDIF

 !  ! for test only - introduced temporarily
 !  p_as%tafo(:,:)  = 292.9_wp
 !  !  - change units to deg C, subtract tmelt (0 deg C, 273.15)
 !  p_as%tafo(:,:)  = p_as%tafo(:,:) - 273.15
 !  p_as%ftdew(:,:) = 289.877
 !  p_as%fu10(:,:)  = 7.84831
 !  p_as%fclou(:,:) = 0.897972
 !  p_as%fswr(:,:)  = 289.489
 !  p_as%u(:,:)     = 0.0_wp
 !  p_as%v(:,:)     = 0.0_wp
 !  p_as%topBoundCond_windStress_u(:,:) = 0.0_wp
 !  p_as%topBoundCond_windStress_v(:,:) = 0.0_wp
 !  p_as%FrshFlux_Precipitation(:,:) = 1.04634e-8
 !  p_as%FrshFlux_Runoff(:,:) = 0.0_wp
 !  p_as%pao(:,:)   = 101300.0_wp

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,4)
    CALL dbg_print('FlxFil: Ext data4-ta/mon1' ,z_c2        ,str_module,3, in_subset=p_patch%cells%owned)
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,4)
    CALL dbg_print('FlxFil: Ext data4-ta/mon2' ,z_c2        ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: p_as%tafo'         ,p_as%tafo   ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: p_as%windStr-u',p_as%topBoundCond_windStress_u, str_module,3,in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: p_as%windStr-v',p_as%topBoundCond_windStress_v, str_module,4,in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: Precipitation' ,p_as%FrshFlux_Precipitation,str_module,3,in_subset=p_patch%cells%owned)
    CALL dbg_print('FlxFil: Runoff'        ,p_as%FrshFlux_Runoff       ,str_module,3,in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    IF (type_surfRelax_Temp == 2)  THEN

      !-------------------------------------------------------------------------
      ! Apply temperature relaxation data (record 3) from stationary forcing
      !  - change units to deg C, subtract tmelt (0 deg C, 273.15)
      !  - this is not done for type_surfRelax_Temp=3, since init-data is in Celsius

       p_as%data_surfRelax_Temp(:,:) = &
         &  rday1*(ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,3)-tmelt) + &
         &  rday2*(ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,3)-tmelt)

    END IF

    IF (type_surfRelax_Salt == 2 .AND. no_tracer >1) THEN

      !-------------------------------------------------------------------------
      ! Apply salinity relaxation data (record ??) from stationary forcing
      CALL finish(TRIM(ROUTINE),' type_surfRelax_Salt=2 (reading from flux file) not yet implemented')

    END IF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    IF (idbg_mxmn >= idt_src) THEN
      WRITE(message_text,'(a,i6,2(a,i4),2(a,f12.8))') 'FLUX time interpolation: jt=',jstep, &
        &  ' mon1=',jmon1,' mon2=',jmon2,' day1=',rday1,' day2=',rday2
      CALL message (' ', message_text)
    END IF
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,1)
    CALL dbg_print('FlxFil: Ext data1-u/mon1'  ,z_c2 ,str_module,idt_src, in_subset=p_patch%cells%owned)
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,1)
    CALL dbg_print('FlxFil: Ext data1-u/mon2'  ,z_c2 ,str_module,idt_src, in_subset=p_patch%cells%owned)
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,2)
    CALL dbg_print('FlxFil: Ext data2-v/mon1'  ,z_c2 ,str_module,idt_src, in_subset=p_patch%cells%owned)
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,2)
    CALL dbg_print('FlxFil: Ext data2-v/mon2'  ,z_c2 ,str_module,idt_src, in_subset=p_patch%cells%owned)
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon1,:,3)
    CALL dbg_print('FlxFil: Ext data3-t/mon1'  ,z_c2 ,str_module,idt_src, in_subset=p_patch%cells%owned)
    z_c2(:,:)=ext_data(1)%oce%flux_forc_mon_c(:,jmon2,:,3)
    CALL dbg_print('FlxFil: Ext data3-t/mon2'  ,z_c2 ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE update_flux_fromFile

  !-------------------------------------------------------------------------
  !
  !> Calculates surface temperature and salinity tracer relaxation
  !!   relaxation terms for tracer equation and surface fluxes are calculated
  !!   in addition to other surface tracer fluxes 
  !!   surface tracer restoring is applied either in apply_surface_relaxation
  !!   or in adding surface relaxation fluxes to total forcing fluxes
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2014)
  !
  SUBROUTINE update_surface_relaxation(p_patch_3D, p_os, p_ice, p_oce_sfc, tracer_no)

    TYPE (t_patch_3D ),    TARGET, INTENT(IN) :: p_patch_3D
    TYPE (t_hydro_ocean_state), INTENT(INOUT) :: p_os
    TYPE (t_sea_ice),              INTENT(IN) :: p_ice
    TYPE (t_ocean_surface)                    :: p_oce_sfc
    INTEGER,                       INTENT(IN) :: tracer_no       !  no of tracer: 1=temperature, 2=salinity

    !Local variables 
    INTEGER                       :: jc, jb
    INTEGER                       :: i_startidx_c, i_endidx_c
    REAL(wp)                      :: relax_strength, thick
    TYPE(t_patch), POINTER        :: p_patch
    REAL(wp),      POINTER        :: t_top(:,:), s_top(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    !-----------------------------------------------------------------------  
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    !-------------------------------------------------------------------------

    t_top => p_os%p_prog(nold(1))%tracer(:,1,:,1)


    IF (tracer_no == 1) THEN  ! surface temperature relaxation
      !
      ! Temperature relaxation activated as additonal forcing term in the tracer equation
      ! implemented as simple time-dependent relaxation (time needed to restore tracer completely back to T*)
      !    F_T  = Q_T/dz = -1/tau * (T-T*) [ K/s ]  (where Q_T is boundary condition for vertical diffusion in [K*m/s])
      ! when using the sign convention
      !   dT/dt = Operators + F_T
      ! i.e. F_T <0 for  T-T* >0 (i.e. decreasing temperature T if T is warmer than relaxation data T*) 
    
      ! EFFECTIVE RESTORING PARAMETER: 1.0_wp/(para_surfRelax_Temp*seconds_per_month)
    
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN

            !relax_strength = (p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) + p_os%p_prog(nold(1))%h(jc,jb)) / &
            !  &       (para_surfRelax_Temp*seconds_per_month)
!           p_oce_sfc%TopBC_Temp_vdiff(jc,jb) = -relax_strength*(t_top(jc,jb)-p_oce_sfc%data_surfRelax_Temp(jc,jb))
            relax_strength = 1.0_wp / (para_surfRelax_Temp*seconds_per_month)

            ! calculate additional temperature restoring rate F_T due to relaxation [K/s]
            p_oce_sfc%TempFlux_Relax(jc,jb) = -relax_strength*(t_top(jc,jb)-p_oce_sfc%data_surfRelax_Temp(jc,jb))

            ! Diagnosed heat flux Q_surf due to relaxation
            !  Q_surf = F_T*dz * (rho*Cp) = -dz/tau*(T-T*) * (rho*Cp)  [W/m2]
            !  HeatFlux_Relax = thick * TempFlux_Relax * (OceanReferenceDensity*clw)
            ! this heat flux is negative if relaxation flux is negative, i.e. heat is released if temperature decreases
            ! this flux is for diagnosis only and not added to tracer forcing

            thick = (p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)+p_os%p_prog(nold(1))%h(jc,jb))
            p_oce_sfc%HeatFlux_Relax(jc,jb) = p_oce_sfc%TempFlux_Relax(jc,jb) * thick * OceanReferenceDensity*clw

          ENDIF
    
        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('UpdSfcRlx:HeatFlx_Rlx[W/m2]',p_oce_sfc%HeatFlux_Relax     ,str_module,2, in_subset=p_patch%cells%owned)
      CALL dbg_print('UpdSfcRlx: T* to relax to'  ,p_oce_sfc%data_surfRelax_Temp,str_module,4, in_subset=p_patch%cells%owned)
      CALL dbg_print('UpdSfcRlx: 1/tau*(T*-T)'    ,p_oce_sfc%TempFlux_Relax     ,str_module,3, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

    ELSE IF (tracer_no == 2) THEN  ! surface salinity relaxation
      !
      ! Salinity relaxation activated as additonal forcing term in the tracer equation
      ! implemented as simple time-dependent relaxation (time needed to restore tracer completely back to S*)
      !    F_S  = -1/tau * (S-S*) [ psu/s ]
      ! when using the sign convention
      !   dS/dt = Operators + F_S
      ! i.e. F_S <0 for  S-S* >0 (i.e. decreasing salinity S if S is saltier than relaxation data S*)
      ! note that the freshwater flux is opposite in sign to F_S, see below,
      ! i.e. fwf >0 for  S-S* >0 (i.e. increasing freshwater flux to decrease salinity)

      s_top => p_os%p_prog(nold(1))%tracer(:,1,:,2)

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN

            relax_strength = 1.0_wp / (para_surfRelax_Salt*seconds_per_month)
            ! 
            ! If sea ice is present (and l_relaxsal_ice), salinity relaxation is proportional to open water,
            !   under sea ice, no relaxation is applied, according to the procedure in MPIOM
            IF (l_relaxsal_ice .AND. i_sea_ice >=1) relax_strength = (1.0_wp-p_ice%concsum(jc,jb))*relax_strength

            ! calculate additional salt restoring rate F_S due to relaxation [psu/s]
            p_oce_sfc%SaltFlux_Relax(jc,jb) = -relax_strength*(s_top(jc,jb)-p_oce_sfc%data_surfRelax_Salt(jc,jb))

            ! Diagnosed freshwater flux due to relaxation (equivalent to heat flux Q)
            !  Fw_S = F_S*dz/S = dz/tau * (S-S*)/S  [m/s]
            ! this flux is applied as volume forcing in surface equation in fill_rhs4surface_eq_ab
            thick = (p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)+p_os%p_prog(nold(1))%h(jc,jb))
            p_oce_sfc%FrshFlux_Relax(jc,jb) = -p_oce_sfc%SaltFlux_Relax(jc,jb) * thick / s_top(jc,jb)

          ENDIF
        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('UpdSfcRlx:s_top ',s_top  ,str_module,4, in_subset=p_patch%cells%owned)
      CALL dbg_print('UpdSfcRlx:FrshFlxRelax[m/s]',p_oce_sfc%FrshFlux_Relax     ,str_module,4, in_subset=p_patch%cells%owned)
      CALL dbg_print('UpdSfcRlx: S* to relax to'  ,p_oce_sfc%data_surfRelax_Salt,str_module,4, in_subset=p_patch%cells%owned)
      CALL dbg_print('UpdSfcRlx: 1/tau*(S*-S)'    ,p_oce_sfc%SaltFlux_Relax     ,str_module,4, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

    END IF  ! tracer_no

  END SUBROUTINE update_surface_relaxation

  !-------------------------------------------------------------------------
  !
  !> Calculates surface temperature and salinity tracer relaxation
  !!   relaxation terms for tracer equation and surface fluxes are calculated
  !!   in addition to other surface tracer fluxes 
  !!   surface tracer restoring is applied either in apply_surface_relaxation
  !!   or in adding surface relaxation fluxes to total forcing fluxes
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2014)
  !
  SUBROUTINE apply_surface_relaxation(p_patch_3D, p_os, p_oce_sfc, tracer_no)

    TYPE (t_patch_3D ),    TARGET, INTENT(IN)    :: p_patch_3D
    TYPE (t_hydro_ocean_state),    INTENT(INOUT) :: p_os
    TYPE (t_ocean_surface), INTENT(IN)           :: p_oce_sfc
    INTEGER,                      INTENT(IN)     :: tracer_no       !  no of tracer: 1=temperature, 2=salinity

    !Local variables 
    INTEGER :: jc, jb
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp) :: t_top_old  (nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: s_top_old  (nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_patch), POINTER :: p_patch
    REAL(wp),      POINTER :: t_top(:,:), s_top(:,:)
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------  
    p_patch         => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------

    all_cells => p_patch%cells%all

    t_top =>p_os%p_prog(nold(1))%tracer(:,1,:,1)
    t_top_old(:,:) = t_top(:,:)


    ! add relaxation term to temperature tracer
    IF (tracer_no == 1) THEN
    
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
    
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN
            t_top_old(jc,jb) = t_top(jc,jb)
            t_top(jc,jb)     = t_top_old(jc,jb) + p_oce_sfc%TempFlux_Relax(jc,jb)*dtime
          ENDIF
    
        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('AppTrcRlx: TempFluxRelax'  , p_oce_sfc%TempFlux_Relax, str_module, 3, in_subset=p_patch%cells%owned)
      CALL dbg_print('AppTrcRlx: Old Temperature', t_top_old                  , str_module, 3, in_subset=p_patch%cells%owned)
      CALL dbg_print('AppTrcRlx: New Temperature', t_top                      , str_module, 2, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

    ! add relaxation term to salinity tracer
    ELSE IF (tracer_no == 2) THEN

      s_top =>p_os%p_prog(nold(1))%tracer(:,1,:,2)
      s_top_old(:,:) = s_top(:,:)

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN
            s_top(jc,jb)     = s_top_old(jc,jb) + p_oce_sfc%SaltFlux_Relax(jc,jb)*dtime
          ENDIF
        END DO
      END DO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('AppTrcRlx: SaltFluxRelax', p_oce_sfc%SaltFlux_Relax, str_module, 3, in_subset=p_patch%cells%owned)
      CALL dbg_print('AppTrcRlx: Old Salt'     , s_top_old                  , str_module, 3, in_subset=p_patch%cells%owned)
      CALL dbg_print('AppTrcRlx: New Salt'     , s_top                      , str_module, 2, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

    END IF  ! tracer_no

  END SUBROUTINE apply_surface_relaxation

  !-------------------------------------------------------------------------
  !
  !>
  !! Update surface flux forcing for hydrostatic ocean
  !!
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
!<Optimize_Used>
  SUBROUTINE update_flux_analytical(p_patch_3D, p_os, atmos_fluxes, p_oce_sfc)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)    :: p_patch_3D
    TYPE(t_hydro_ocean_state), INTENT(IN)   :: p_os
    TYPE(t_atmos_fluxes)                    :: atmos_fluxes
    TYPE(t_ocean_surface)                   :: p_oce_sfc
    !
    ! local variables
    INTEGER :: jc, jb
    INTEGER :: i_startblk_c, i_endblk_c, start_cell_index, end_cell_index
    !INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
    !INTEGER :: rl_start_c, rl_end_c

    REAL(wp) :: z_lat, z_lon, z_lat_deg
    !REAL(wp) :: y_length               !basin extension in y direction in degrees
    REAL(wp) :: z_T_init(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: z_perlat, z_perlon, z_permax, z_perwid, z_relax, z_dst
    INTEGER  :: z_dolic
    REAL(wp) :: z_temp_max, z_temp_min, z_temp_incr
    REAL(wp) :: center, length, zonal_waveno,amplitude
    REAL(wp) :: no_flux_length, south_bound, max_flux_y

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_surface:update_flux_analytical'
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------  
    p_patch         => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------
    all_cells => p_patch%cells%all

    ! atmosphere fluxes for analytical testcased similar to mo_ocean_initial_conditions:
    SELECT CASE (atmos_flux_analytical_type)

    CASE(0)  !  all fluxes are unchanged, zero as default
      CONTINUE

    CASE(101,102,103)  !  constant fluxes for test of sea-ice processes
      ! set LW/SW/sensible/latent heat fluxes over ice to constant
      atmos_fluxes%SWnet(:,1,:) = atmos_SWnet_const
      atmos_fluxes%LWnet(:,1,:) = atmos_LWnet_const
      atmos_fluxes%sens (:,1,:) = atmos_sens_const
      atmos_fluxes%lat  (:,1,:) = atmos_lat_const
      ! set LW/SW/sensible/latent heat fluxes over water to constant
      atmos_fluxes%SWnetw(:,:)  = atmos_SWnetw_const
      atmos_fluxes%LWnetw(:,:)  = atmos_LWnetw_const
      atmos_fluxes%sensw(:,:)   = atmos_sensw_const
      atmos_fluxes%latw(:,:)    = atmos_lat_const
      ! set water fluxes over water to constant
      atmos_fluxes%FrshFlux_Precipitation(:,:) = atmos_precip_const
  
      CASE(200) 

        IF(no_tracer>=1.AND.type_surfRelax_Temp==0)THEN

          center = basin_center_lat * deg2rad
          length = basin_height_deg * deg2rad
          zonal_waveno = 2.0_wp
          amplitude    =10.0_wp

          DO jb = all_cells%start_block, all_cells%end_block
            CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
            DO jc = start_cell_index, end_cell_index

              IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN

                z_lat    = p_patch%cells%center(jc,jb)%lat
                z_lon    = p_patch%cells%center(jc,jb)%lon
  
                p_oce_sfc%data_surfRelax_Temp(jc,jb) =0.0_wp

                p_oce_sfc%TopBC_Temp_vdiff(jc,jb) &
                &= amplitude * COS(zonal_waveno*pi*(z_lat-center)/length)
                
                !IF(z_lat<= center-0.5_wp*length)p_oce_sfc%TopBC_Temp_vdiff(jc,jb)=0.0_wp
                !IF(z_lat>= center+0.0_wp*length)p_oce_sfc%TopBC_Temp_vdiff(jc,jb)=10.0_wp
                !IF(z_lat>= center+0.75_wp*length)p_oce_sfc%TopBC_Temp_vdiff(jc,jb)=&
                !& -10.0_wp!amplitude * (COS(zonal_waveno*pi*(z_lat-center)/length))  
              ENDIF 
            END DO
          END DO
        ENDIF

      CASE(201) ! Abernathey 2011

        IF(no_tracer>=1.AND.type_surfRelax_Temp==0)THEN

          no_flux_length = 3.0_wp * relax_width * deg2rad
          south_bound = (basin_center_lat - 0.5_wp * basin_height_deg) * deg2rad
          length      = basin_height_deg * deg2rad
          max_flux_y  = length - no_flux_length
          zonal_waveno =  2.5_wp
          amplitude    = forcing_HeatFlux_amplitude

          DO jb = all_cells%start_block, all_cells%end_block
            CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
            DO jc = start_cell_index, end_cell_index
              p_oce_sfc%data_surfRelax_Temp(jc,jb)     = 0.0_wp
              p_oce_sfc%TopBC_Temp_vdiff(jc,jb) = 0.0_wp
              
              z_lat = p_patch%cells%center(jc,jb)%lat - south_bound
              
              IF(p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary .AND. z_lat < max_flux_y) THEN

                p_oce_sfc%TopBC_Temp_vdiff(jc,jb) = &
                  & - amplitude * COS(zonal_waveno * pi * z_lat/max_flux_y) &
                  & + forcing_HeatFlux_base

              ENDIF
            END DO
          END DO
        ENDIF
        
      CASE default

        CALL finish(routine, "unknown atmos_flux_analytical_type")

    END SELECT

    SELECT CASE (relax_analytical_type)

    CASE(27,30,32)

     IF(no_tracer>=1.AND.type_surfRelax_Temp/=0)THEN

       !y_length = basin_height_deg * deg2rad
       DO jb = all_cells%start_block, all_cells%end_block
         CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
         DO jc = start_cell_index, end_cell_index

           IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN

             z_T_init(jc,jb) = 20.0_wp- p_patch_3D%p_patch_1D(1)%zlev_m(1)*15.0_wp/4000.0_wp

             z_lat    = p_patch%cells%center(jc,jb)%lat
             z_lon    = p_patch%cells%center(jc,jb)%lon

             ! Add temperature perturbation at new values
             z_perlat = basin_center_lat + 0.1_wp*basin_height_deg
             z_perlon = basin_center_lon + 0.1_wp*basin_width_deg
             z_permax = 0.1_wp
             z_perwid = 10.0_wp

             z_relax  = para_surfRelax_Temp/(30.0_wp*24.0_wp*3600.0_wp)

             z_dolic  = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
             IF (z_dolic > MIN_DOLIC) THEN

               z_dst = sqrt((z_lat-z_perlat*deg2rad)**2+(z_lon-z_perlon*deg2rad)**2)

               IF(z_dst <= 5.0_wp*deg2rad)THEN
                 z_T_init = z_T_init &
                 &        + z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2) &
                 &        * sin(pi*p_patch_3D%p_patch_1D(1)%zlev_m(1)/4000.0_wp)
               ENDIF
               ! up to here z_init is identically initialized than temperature

               !add local cold perturbation 
               IF(z_dst <= 10.5_wp*deg2rad)THEN
                 z_T_init(jc,jb) = z_T_init(jc,jb) - exp(-(z_dst/(z_perwid*deg2rad))**2)
               ENDIF

               p_oce_sfc%data_surfRelax_Temp(jc,jb)     = z_T_init(jc,jb)

               p_oce_sfc%TopBC_Temp_vdiff(jc,jb) = z_relax * &          
                 &  ( p_oce_sfc%data_surfRelax_Temp(jc,jb)-p_os%p_prog(nold(1))%tracer(jc,1,jb,1) )

             END IF
           ELSE
             atmos_fluxes%topBoundCond_windStress_cc(jc,jb)%x(:) = 0.0_wp
!            atmos_fluxes%topBoundCond_windStress_u(jc,jb)       = 0.0_wp
!            atmos_fluxes%topBoundCond_windStress_v(jc,jb)       = 0.0_wp
           ENDIF 
       END DO
      END DO

    ENDIF

    CASE (33)
      IF(type_surfRelax_Temp>=1)THEN
        z_relax = para_surfRelax_Temp/(30.0_wp*24.0_wp*3600.0_wp)

        p_oce_sfc%TopBC_Temp_vdiff(:,:) = z_relax*( p_oce_sfc%data_surfRelax_Temp(:,:) &
          &                                               -p_os%p_prog(nold(1))%tracer(:,1,:,1) )

      END IF
  
    CASE(51)

      IF(type_surfRelax_Temp>=1)THEN

        z_relax = para_surfRelax_Temp/(30.0_wp*24.0_wp*3600.0_wp)

        z_temp_max  = 30.5_wp
        z_temp_min  = 0.5_wp
        z_temp_incr = (z_temp_max-z_temp_min)/(n_zlev-1.0_wp)

      !Add horizontal variation
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          z_lat = p_patch%cells%center(jc,jb)%lat
          z_lat_deg = z_lat*rad2deg

            IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN

              z_temp_max     =0.01_wp*(z_lat_deg-basin_center_lat)*(z_lat_deg-basin_center_lat)
              z_T_init(jc,jb)=30.5_wp

              z_T_init(jc,jb)&
              &=z_T_init(jc,jb)*exp(-z_temp_max/basin_height_deg)
            ELSE
              z_T_init(jc,jb)=0.0_wp
            ENDIF
        END DO
      END DO
      p_oce_sfc%data_surfRelax_Temp(:,:)=z_T_init(:,:)

      p_oce_sfc%TopBC_Temp_vdiff(:,:) = z_relax*( p_oce_sfc%data_surfRelax_Temp(:,:) &
        &                                               -p_os%p_prog(nold(1))%tracer(:,1,:,1) )

      END IF

    END SELECT

  END SUBROUTINE update_flux_analytical

  !-------------------------------------------------------------------------
  !
  !> Calc_omip_budgets_ice equals sbr "Budget" in MPIOM.
  !! Sets the atmospheric fluxes over *SEA ICE ONLY* for the update of the ice
  !! temperature and ice growth rates for OMIP forcing
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-07). Originally code written by
  !! Dirk Notz, following MPIOM. Code transfered to ICON.
  !! Einar Olason, split calc_atm_fluxes_from_bulk into calc_bulk_flux_ice and calc_bulk_flux_oce
  !! so that the ocean model can be run without the ice model, but with OMIP fluxes.
  !!
  !! Rewritten by Stephan Lorenz, MPI-M (2015-06).
  !!  Using interface with parameters in order to call budget routine independent of ocean model

  SUBROUTINE calc_omip_budgets_ice(geolat, tafoC, ftdewC, fu10, fclou, pao, fswr,                &
    &                              kice, tice, hice, albvisdir, albvisdif, albnirdir, albnirdif, &
    &                              LWnetIce, SWnetIce, sensIce, latentIce,                       &
    &                              dLWdTIce, dsensdTIce, dlatdTIce)                              

 !  INPUT variables for OMIP via parameter:
    REAL(wp), INTENT(in)    :: geolat(:,:)      ! latitude                             [rad]
    REAL(wp), INTENT(in)    :: tafoC(:,:)       ! 2 m air temperature in Celsius       [C]
    REAL(wp), INTENT(in)    :: ftdewC(:,:)      ! 2 m dew point temperature in Celsius [C]
    REAL(wp), INTENT(in)    :: fu10(:,:)        ! 10 m wind speed                      [m/s]
    REAL(wp), INTENT(in)    :: fclou(:,:)       ! Fractional cloud cover               [frac]
    REAL(wp), INTENT(in)    :: pao(:,:)         ! Surface atmospheric pressure         [hPa]
    REAL(wp), INTENT(in)    :: fswr(:,:)        ! Incoming surface solar radiation     [W/m2]
    INTEGER,  INTENT(in)    :: kice             ! number of ice classes (currently 1)
    REAL(wp), INTENT(in)    :: tice(:,:,:)      ! surface ice temperature per class    [C]
    REAL(wp), INTENT(in)    :: hice(:,:,:)      ! ice thickness per class              [m]
    REAL(wp), INTENT(in)    :: albvisdir(:,:,:) ! direct ice albedo per class
    REAL(wp), INTENT(in)    :: albvisdif(:,:,:) ! diffuse ice albedo per class
    REAL(wp), INTENT(in)    :: albnirdir(:,:,:) ! direct near infrared ice albedo per class
    REAL(wp), INTENT(in)    :: albnirdif(:,:,:) ! diffuse near infrared ice albedo per class

 !  OUTPUT variables for sea ice model via parameter (inout since icefree part is not touched)
    REAL(wp), INTENT(inout) :: LWnetIce (:,:,:) ! net longwave heat flux over ice      [W/m2]
    REAL(wp), INTENT(inout) :: SWnetIce (:,:,:) ! net shortwave heat flux over ice     [W/m2]
    REAL(wp), INTENT(inout) :: sensIce  (:,:,:) ! sensible heat flux over ice          [W/m2]
    REAL(wp), INTENT(inout) :: latentIce(:,:,:) ! latent heat flux over ice            [W/m2]
    REAL(wp), INTENT(inout) :: dLWdTIce (:,:,:) ! derivitave of LWnetIce w.r.t temperature
    REAL(wp), INTENT(inout) :: dsensdTIce(:,:,:)! derivitave of sensIce w.r.t temperature
    REAL(wp), INTENT(inout) :: dlatdTIce(:,:,:) ! derivitave of latentIce w.r.t temperature

 !  Local variables
 !  REAL(wp), DIMENSION (nproma,p_patch%alloc_cell_blocks) :: &
    REAL(wp), DIMENSION (SIZE(tafoC,1), SIZE(tafoC,2)) ::  &
      & Tsurf,          &  ! Surface temperature in Celsius                  [C]
      & tafoK,          &  ! Air temperature at 2 m in Kelvin                [K]
      & fu10lim,        &  ! wind speed at 10 m height in range 2.5...32     [m/s]
      & esta,           &  ! water vapor pressure at 2 m height              [Pa]
      & esti,           &  ! water vapor pressure at ice surface             [Pa]
      & sphumida,       &  ! Specific humididty at 2 m height
      & sphumidi,       &  ! Specific humididty at ice surface
      & rhoair,         &  ! air density                                     [kg/m^3]
      & dragl0,         &  ! part of dragl
      & dragl1,         &  ! part of dragl
      & dragl,          &  ! Drag coefficient for latent   heat flux
      & drags,          &  ! Drag coefficient for sensible heat flux (=0.95 dragl)
      & fakts,          &  ! Effect of cloudiness on LW radiation
      & humi,           &  ! Effect of air humidity on LW radiation
      & fa, fi,         &  ! Enhancment factor for vapor pressure
      & dsphumididesti, &  ! Derivative of sphumidi w.r.t. esti
      & destidT,        &  ! Derivative of esti w.r.t. T
      & dfdT               ! Derivative of f w.r.t. T
 !    & wspeed             ! Wind speed                                      [m/s]

    INTEGER :: i
    REAL(wp) :: aw,bw,cw,dw,ai,bi,ci,di,AAw,BBw,CCw,AAi,BBi,CCi,alpha,beta
    REAL(wp) :: fvisdir, fvisdif, fnirdir, fnirdif, local_rad2deg

    tafoK(:,:)  = tafoC(:,:) + tmelt  ! Change units of tafo  to Kelvin

    ! set to zero for NAG
    sphumida(:,:)  = 0.0_wp
    fa      (:,:)  = 0.0_wp
    esta    (:,:)  = 0.0_wp
    rhoair  (:,:)  = 0.0_wp

    !-----------------------------------------------------------------------
    ! Compute water vapor pressure and specific humididty in 2m height (esta)
    ! and at water surface (estw) according to "Buck Research Manual (1996)
    ! (see manuals for instruments at http://www.buck-research.com/);
    ! updated from Buck, A. L., New equations for computing vapor pressure and
    ! enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981"
    !-----------------------------------------------------------------------
    ! #slo# 2015-03: the comment above is now valid 
    ! the values for ice are not changed in Buck (1996) in comparison to Buck (1981)

    ! the following commented values are from Buck (1981)
    ! aw=611.21_wp; bw=18.729_wp; cw=257.87_wp; dw=227.3_wp

    ! here are the updated values for open water according to Buck (1996)
    ai=611.15_wp; bi=23.036_wp; ci=279.82_wp; di=333.7_wp
    aw=611.21_wp; bw=18.678_wp; cw=257.14_wp; dw=234.5_wp

    AAw=7.2e-4_wp; BBw=3.20e-6_wp; CCw=5.9e-10_wp
    AAi=2.2e-4_wp; BBi=3.83e-6_wp; CCi=6.4e-10_wp

    alpha=0.62197_wp; beta=0.37803_wp

    ! #slo# correction: pressure in enhancement formula is in mb (hPa) according to Buck 1981 and 1996
    fa(:,:)        = 1.0_wp+AAw+pao*0.01_wp*(BBw+CCw*ftdewC**2)
    esta(:,:)      = fa * aw*EXP((bw-ftdewC/dw)*ftdewC/(ftdewC+cw))
    sphumida(:,:)  = alpha * esta/(pao-beta*esta)

    !-----------------------------------------------------------------------
    !  Compute longwave radiation according to
    !         Berliand, M. E., and T. G. Berliand, 1952: Determining the net
    !         long-wave radiation of the Earth with consideration of the effect
    !         of cloudiness. Izv. Akad. Nauk SSSR, Ser. Geofiz., 1, 6478.
    !         cited by: Budyko, Climate and Life, 1974.
    !         Note that for humi, esta is given in [mmHg] in the original
    !         publication. Therefore, 0.05*sqrt(esta/100) is used rather than
    !         0.058*sqrt(esta)
    !  This is the formula used in MPI-OM when using the QLOBERL preprocessing option (currently
    !  the default usage).
    !-----------------------------------------------------------------------

    ! Berliand & Berliand ('52) calculate only LWnet
    humi    = 0.39_wp - 0.05_wp*SQRT(esta/100._wp)

    ! icon-identical calculation of rad2deg
    local_rad2deg = 180.0_wp / 3.14159265358979323846264338327950288_wp
    fakts   =  1.0_wp - ( 0.5_wp + 0.4_wp/90._wp &
      &         *MIN(ABS(local_rad2deg*geolat(:,:)),60._wp) ) * fclou(:,:)**2
  !   &         *MIN(ABS(rad2deg*p_patch%cells%center(:,:)%lat),60._wp) ) * fclou(:,:)**2

    !-----------------------------------------------------------------------
    !  Calculate bulk equations according to
    !      Kara, B. A., P. A. Rochford, and H. E. Hurlburt, 2002:
    !      Air-Sea Flux Estimates And The 19971998 Enso Event,  Bound.-Lay.
    !      Met., 103(3), 439-458, doi: 10.1023/A:1014945408605.
    !-----------------------------------------------------------------------

    ! with nag there is floating invalid operation on rest of last nproma-block only due to pao=nan
    ! rhoair(:,:) = pao(:,:) / (rd*tafoK(:,:)*(1.0_wp+0.61_wp*sphumida(:,:)) ) !  error with nag
    WHERE (pao(:,:)>0.0_wp) rhoair(:,:) = pao(:,:) / (rd*tafoK(:,:)*(1.0_wp+0.61_wp*sphumida(:,:)) )

    fu10lim(:,:)    = MAX (2.5_wp, MIN(32.5_wp,fu10(:,:)) )
    dragl1(:,:)     = 1e-3_wp*(-0.0154_wp + 0.5698_wp/fu10lim(:,:) &
      &               - 0.6743_wp/(fu10lim(:,:) * fu10lim(:,:)))
    dragl0(:,:)     = 1e-3_wp*(0.8195_wp+0.0506_wp*fu10lim(:,:) &
      &               - 0.0009_wp*fu10lim(:,:)*fu10lim(:,:))

    ! Fractions of SWin in each band (from cice)
    fvisdir=0.28_wp; fvisdif=0.24_wp; fnirdir=0.31_wp; fnirdif=0.17_wp
    Tsurf(:,:) = 0.0_wp ! For debug output

    ! Over sea ice area only
    !  TODO: in case of no ice model, ice variables cannot be used here
    !  ice classes: currently one class (kice=1) is used, therefore formulation can be simplified to 2-dim variables as in mpiom
    DO i = 1, kice
      WHERE (hice(:,i,:)>0._wp)
        
        !  albedo model: atmos_fluxes%albvisdir, albvisdif, albnirdir, albnirdif
        !   - all 4 albedos are the same (i_ice_albedo = 1), they are calculated in ice_fast and should be stored in p_ice
        SWnetIce(:,i,:)  = ( 1._wp-albvisdir(:,i,:) )*fvisdir*fswr(:,:) +   &
          &                ( 1._wp-albvisdif(:,i,:) )*fvisdif*fswr(:,:) +   &
          &                ( 1._wp-albnirdir(:,i,:) )*fnirdir*fswr(:,:) +   &
          &                ( 1._wp-albnirdif(:,i,:) )*fnirdif*fswr(:,:)
      ! Tsurf(:,:)       = p_ice%Tsurf(:,i,:)
        Tsurf(:,:)       = tice(:,i,:)
        ! pressure in enhancement formula is in mb (hPa) according to Buck 1981 and 1996
        fi(:,:)          = 1.0_wp+AAi+pao(:,:)*0.01_wp*(BBi+CCi*Tsurf(:,:) **2)
        esti(:,:)        = fi(:,:)*ai*EXP((bi-Tsurf(:,:) /di)*Tsurf(:,:) /(Tsurf(:,:) +ci))
        sphumidi(:,:)    = alpha*esti(:,:)/(pao(:,:)-beta*esti(:,:))
        ! This may not be the best drag parametrisation to use over ice
        dragl(:,:)       = dragl0(:,:) + dragl1(:,:) * (Tsurf(:,:)-tafoC(:,:))
        ! A reasonableee maximum and minimum is needed for dragl in case there's a large difference
        ! between the 2-m and surface temperatures.
        dragl(:,:)       = MAX(0.5e-3_wp, MIN(3.0e-3_wp,dragl(:,:)))
        drags(:,:)       = 0.95_wp * dragl(:,:)

        LWnetIce(:,i,:)  = -fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4 &
           &               -4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - tafoC(:,:))
        ! same form as MPIOM:
        !atmos_fluxes%LWnet (:,i,:)  = - (fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4 &
        !  &         + 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:)))
        dLWdTIce(:,i,:)  = -4._wp*zemiss_def*stbo*tafoK(:,:)**3
        sensIce(:,i,:)   = drags(:,:) * rhoair(:,:)*cpd*fu10(:,:) * fr_fac * (tafoC(:,:) -Tsurf(:,:))
        latentIce(:,i,:) = dragl(:,:) * rhoair(:,:)* alf *fu10(:,:) * fr_fac &
          &                   * (sphumida(:,:)-sphumidi(:,:))

        dsensdTIce(:,i,:)   = 0.95_wp*cpd*rhoair(:,:)*fu10(:,:)&
          &                   *(dragl0(:,:) - 2.0_wp*dragl(:,:))
        dsphumididesti(:,:) = alpha/(pao(:,:)-beta*esti(:,:)) &
          &                   * (1.0_wp + beta*esti(:,:)/(pao(:,:)-beta*esti(:,:)))
        destidT(:,:)        = (bi*ci*di-Tsurf(:,:)*(2.0_wp*ci+Tsurf(:,:)))&
          &                   /(di*(ci+Tsurf(:,:))**2) * esti(:,:)
        dfdT(:,:)           = 2.0_wp*CCi*BBi*Tsurf(:,:)
        dlatdTIce(:,i,:)    = alf*rhoair(:,:)*fu10(:,:)* &
          &                  ( (sphumida(:,:)-sphumidi(:,:))*dragl1(:,:) &
          &                    - dragl(:,:)*dsphumididesti(:,:)*(fi(:,:)*destidT(:,:) &
          &                    + esti(:,:)*dfdT(:,:)) )
      ENDWHERE
    ENDDO

  ! IF (use_calculated_ocean_stress) THEN
  !   !-----------------------------------------------------------------------
  !   !  Calculate wind stress over ice covered part
  !   !  TODO: should be moved to ice_ocean_stress
  !   !-----------------------------------------------------------------------
  !   WHERE (hice(:,i,:)>0._wp)
  !     wspeed(:,:) = SQRT( p_as%u**2 + p_as%v**2 )
  !     atmos_fluxes%stress_x(:,:) = Cd_ia*rhoair(:,:)*wspeed(:,:)*p_as%u(:,:)
  !     atmos_fluxes%stress_y(:,:) = Cd_ia*rhoair(:,:)*wspeed(:,:)*p_as%v(:,:)
  !   ENDWHERE
  ! ELSE
  !   ! use wind stress provided by OMIP data
  !   atmos_fluxes%stress_x(:,:) = p_as%topBoundCond_windStress_u(:,:)
  !   atmos_fluxes%stress_y(:,:) = p_as%topBoundCond_windStress_v(:,:)
  ! ENDIF

  END SUBROUTINE calc_omip_budgets_ice

  !-------------------------------------------------------------------------
  !
  !> Forcing_from_bulk equals sbr "Budget_omip" in MPIOM.
  !! Sets the atmospheric fluxes for the update of
  !! temperature of *OPEN WATER* for OMIP forcing.
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2012-08). Originally code written by
  !! Dirk Notz, following MPIOM. Code transfered to ICON.
   
  SUBROUTINE calc_omip_budgets_oce(p_patch, p_as, p_os, atmos_fluxes)
    TYPE(t_patch),            INTENT(IN), TARGET    :: p_patch
    TYPE(t_atmos_for_ocean),  INTENT(IN)    :: p_as
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os
    TYPE(t_atmos_fluxes),     INTENT(INOUT) :: atmos_fluxes

 !  INPUT variables:
 !  p_as%tafo(:,:),      : 2 m air temperature                              [C]
 !  p_as%ftdew(:,:),     : 2 m dew-point temperature                        [K]
 !  p_as%fu10(:,:) ,     : 10 m wind speed                                  [m/s]
 !  p_as%fclou(:,:),     : Fractional cloud cover
 !  p_as%pao(:,:),       : Surface atmospheric pressure                     [hPa]
 !  p_as%fswr(:,:),      : Incoming surface solar radiation                 [W/m]
 !  p_os%tracer(:,1,:,1) : SST
 !  atmos_fluxes%albvisdirw, albvisdifw, albnirdirw, albnirdifw
 !
 !  OUTPUT variables:  atmos_fluxes - heat fluxes and wind stress over open ocean
 !  atmos_fluxes%LWnetw   : long wave
 !  atmos_fluxes%SWnetw   : short wave
 !  atmos_fluxes%sensw    : sensible
 !  atmos_fluxes%latw     : latent
 !  atmos_fluxes%stress_xw: zonal stress
 !  atmos_fluxes%stress_yw: meridional stress

 !  Local variables
    REAL(wp), DIMENSION (nproma,p_patch%alloc_cell_blocks) ::           &
      & Tsurf,          &  ! Surface temperature                             [C]
      & tafoK,          &  ! Air temperature at 2 m in Kelvin                [K]
      & fu10lim,        &  ! wind speed at 10 m height in range 2.5...32     [m/s]
      & esta,           &  ! water vapor pressure at 2 m height              [Pa]
      & estw,           &  ! water vapor pressure at water surface           [Pa]
      & sphumida,       &  ! Specific humididty at 2 m height
      & sphumidw,       &  ! Specific humididty at water surface
      & ftdewC,         &  ! Dew point temperature in Celsius                [C]
      & rhoair,         &  ! air density                                     [kg/m^3]
      & dragl0,         &  ! part of dragl
      & dragl1,         &  ! part of dragl
      & dragl,          &  ! Drag coefficient for latent   heat flux
      & drags,          &  ! Drag coefficient for sensible heat flux (=0.95 dragl)
      & fakts,          &  ! Effect of cloudiness on LW radiation
      & humi,           &  ! Effect of air humidity on LW radiation
      & fa, fw,         &  ! Enhancment factor for vapor pressure
      & wspeed,         &  ! Wind speed                                      [m/s]
      & C_ao               ! Drag coefficient for atm-ocean stress           [m/s]

    INTEGER :: jb, jc, i_startidx_c, i_endidx_c
    REAL(wp) :: aw,bw,cw,dw,AAw,BBw,CCw,alpha,beta
    REAL(wp) :: fvisdir, fvisdif, fnirdir, fnirdif

    TYPE(t_subset_range), POINTER :: all_cells

    Tsurf(:,:)  = p_os%p_prog(nold(1))%tracer(:,1,:,1)  ! set surface temp = mixed layer temp
    tafoK(:,:)  = p_as%tafo(:,:)  + tmelt               ! Change units of tafo  to Kelvin
    ftdewC(:,:) = p_as%ftdew(:,:) - tmelt               ! Change units of ftdew to Celsius

    ! subset range pointer
    all_cells => p_patch%cells%all

    !-----------------------------------------------------------------------
    ! Compute water vapor pressure and specific humididty in 2m height (esta)
    ! and at water surface (estw) according to "Buck Research Manual (1996)
    ! (see manuals for instruments at http://www.buck-research.com/);
    ! updated from Buck, A. L., New equations for computing vapor pressure and
    ! enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981"
    !-----------------------------------------------------------------------
    AAw   = 7.2e-4_wp;  BBw  = 3.20e-6_wp; CCw = 5.9e-10_wp
    alpha = 0.62197_wp; beta = 0.37803_wp
    ! #slo# 2015-03: the comment above is now valid - the following commented values are from Buck (1981)
    ! aw    = 611.21_wp; bw    = 18.729_wp;  cw  = 257.87_wp; dw = 227.3_wp
    ! these are the updated values according to Buck (1996)
    aw    = 611.21_wp; bw    = 18.678_wp;  cw  = 257.14_wp; dw = 234.5_wp

    ! #slo# correction: pressure in enhancement formula is in mb (hPa) according to Buck 1981 and 1996
   !fa(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*(BBw+CCw*ftdewC(:,:)**2)
    fa(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*0.01_wp*(BBw+CCw*ftdewC(:,:)**2)
    esta(:,:) = fa(:,:) * aw*EXP((bw-ftdewC(:,:)/dw)*ftdewC(:,:)/(ftdewC(:,:)+cw))
   !esta(:,:) =           aw*EXP((bw-ftdewC(:,:)/dw)*ftdewC(:,:)/(ftdewC(:,:)+cw))
   !fw(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*(BBw+CCw*Tsurf(:,:) **2)
    fw(:,:)   = 1.0_wp+AAw+p_as%pao(:,:)*0.01_wp*(BBw+CCw*Tsurf(:,:) **2)
   !estw(:,:) = fw(:,:) *aw*EXP((bw-Tsurf(:,:) /dw)*Tsurf(:,:) /(Tsurf(:,:) +cw))
    ! For a given surface salinity we should multiply estw with  1 - 0.000537*S
    ! #slo# correction according to MPIOM: lowering of saturation vapor pressure over saline water
    !       is taken constant to 0.9815
    estw(:,:) = 0.9815_wp*fw(:,:)*aw*EXP((bw-Tsurf(:,:) /dw)*Tsurf(:,:) /(Tsurf(:,:) +cw))

    sphumida(:,:)  = alpha * esta(:,:)/(p_as%pao(:,:)-beta*esta(:,:))
    sphumidw(:,:)  = alpha * estw(:,:)/(p_as%pao(:,:)-beta*estw(:,:))

    !-----------------------------------------------------------------------
    !  Compute longwave radiation according to
    !         Berliand, M. E., and T. G. Berliand, 1952: Determining the net
    !         long-wave radiation of the Earth with consideration of the effect
    !         of cloudiness. Izv. Akad. Nauk SSSR, Ser. Geofiz., 1, 6478.
    !         cited by: Budyko, Climate and Life, 1974.
    !         Note that for humi, esta is given in [mmHg] in the original
    !         publication. Therefore, 0.05*sqrt(esta/100) is used rather than
    !         0.058*sqrt(esta)
    !  This is the formula used in MPI-OM when using the QLOBERL preprocessing option (currently
    !  the default usage).
    !-----------------------------------------------------------------------

    humi(:,:)    = 0.39_wp - 0.05_wp*SQRT(esta(:,:)/100._wp)
    fakts(:,:)   =  1.0_wp - ( 0.5_wp + 0.4_wp/90._wp &
      &         *MIN(ABS(rad2deg*p_patch%cells%center(:,:)%lat),60._wp) ) * p_as%fclou(:,:)**2
    ! Berliand & Berliand ('52) calculate only LWnetw

    ! #eoo# 2012-12-14: another bugfix
    ! #slo# #hha# 2012-12-13: bugfix, corrected form
    atmos_fluxes%LWnetw(:,:) = - fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4  &
      &                - 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:))
    ! same form as MPIOM:
    !atmos_fluxes%LWnetw(:,:) = - (fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4  &
    !  &         + 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:)))
    ! bug
    !atmos_fluxes%LWnetw(:,:) = fakts(:,:) * humi(:,:) * zemiss_def*stbo * tafoK(:,:)**4  &
    !  &         - 4._wp*zemiss_def*stbo*tafoK(:,:)**3 * (Tsurf(:,:) - p_as%tafo(:,:))

    ! Fractions of SWin in each band (from cice)
    fvisdir=0.28_wp; fvisdif=0.24_wp; fnirdir=0.31_wp; fnirdif=0.17_wp
    atmos_fluxes%SWnetw(:,:) = ( 1._wp-atmos_fluxes%albvisdirw(:,:) )*fvisdir*p_as%fswr(:,:) +   &
      &                ( 1._wp-atmos_fluxes%albvisdifw(:,:) )*fvisdif*p_as%fswr(:,:) +   &
      &                ( 1._wp-atmos_fluxes%albnirdirw(:,:) )*fnirdir*p_as%fswr(:,:) +   &
      &                ( 1._wp-atmos_fluxes%albnirdifw(:,:) )*fnirdif*p_as%fswr(:,:)

    !-----------------------------------------------------------------------
    !  Calculate bulk equations according to
    !      Kara, B. A., P. A. Rochford, and H. E. Hurlburt, 2002:
    !      Air-Sea Flux Estimates And The 19971998 Enso Event,  Bound.-Lay.
    !      Met., 103(3), 439-458, doi: 10.1023/A:1014945408605.
    !-----------------------------------------------------------------------

    rhoair(:,:) = 0._wp
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c,i_endidx_c

        rhoair(jc,jb) = p_as%pao(jc,jb)                &
          &            /(rd*tafoK(jc,jb)*(1.0_wp+0.61_wp*sphumida(jc,jb)) )

      END DO
    END DO

    fu10lim(:,:)    = MAX (2.5_wp, MIN(32.5_wp,p_as%fu10(:,:)) )
    dragl1(:,:)     = 1e-3_wp*(-0.0154_wp + 0.5698_wp/fu10lim(:,:) &
      &               - 0.6743_wp/(fu10lim(:,:) * fu10lim(:,:)))
    dragl0(:,:)     = 1e-3_wp*(0.8195_wp+0.0506_wp*fu10lim(:,:) &
      &               - 0.0009_wp*fu10lim(:,:)*fu10lim(:,:))
    dragl(:,:)      = dragl0(:,:) + dragl1(:,:) * (Tsurf(:,:)-p_as%tafo(:,:))
    ! A reasonable maximum and minimum is needed for dragl in case there's a large difference
    ! between the 2-m and surface temperatures.
    dragl(:,:)      = MAX(0.5e-3_wp, MIN(3.0e-3_wp,dragl(:,:)))
    drags(:,:)      = 0.95_wp * dragl(:,:)
    atmos_fluxes%sensw(:,:) = drags(:,:)*rhoair(:,:)*cpd*p_as%fu10(:,:) * fr_fac &
      &               * (p_as%tafo(:,:) -Tsurf(:,:))
    atmos_fluxes%latw(:,:)  = dragl(:,:)*rhoair(:,:)*alv*p_as%fu10(:,:) * fr_fac &
      &               * (sphumida(:,:)-sphumidw(:,:))

    IF (use_calculated_ocean_stress) THEN
      !-----------------------------------------------------------------------
      !  Calculate oceanic wind stress according to:
      !   Gill (Atmosphere-Ocean Dynamics, 1982, Academic Press) (see also Smith, 1980, J. Phys
      !   Oceanogr., 10, 709-726)
      !-----------------------------------------------------------------------

      wspeed(:,:) = SQRT( p_as%u**2 + p_as%v**2 )
      C_ao(:,:)   = MIN( 2._wp, MAX(1.1_wp, 0.61_wp+0.063_wp*wspeed ) )*1e-3_wp
      atmos_fluxes%stress_xw(:,:) = C_ao(:,:)*rhoair*wspeed(:,:)*p_as%u(:,:)
      atmos_fluxes%stress_yw(:,:) = C_ao(:,:)*rhoair*wspeed(:,:)*p_as%v(:,:)
    ELSE
      ! use wind stress provided by OMIP data
      atmos_fluxes%stress_xw(:,:) = p_as%topBoundCond_windStress_u(:,:)
      atmos_fluxes%stress_yw(:,:) = p_as%topBoundCond_windStress_v(:,:)
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5          , fix)
    CALL dbg_print('omipBudOce:tafoK'              , tafoK                 , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:tafo'               , p_as%tafo             , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:ftdew'              , p_as%ftdew            , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:ftdewC'             , ftdewC                , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:pao'                , p_as%pao              , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:fa'                 , fa                    , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:fw'                 , fw                    , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:esta'               , esta                  , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:estw'               , estw                  , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:sphumida'           , sphumida              , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:sphumidw'           , sphumidw              , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:rhoair'             , rhoair                , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:dragl'              , dragl                 , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:drags'              , drags                 , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:fu10'               , p_as%fu10             , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:fu10lim'            , fu10lim               , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:stress_xw'          , atmos_fluxes%stress_xw, str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:stress_yw'          , atmos_fluxes%stress_yw, str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:p_as%windStr-u',p_as%topBoundCond_windStress_u,str_module,idt_src, in_subset=p_patch%cells%owned)
    idt_src=3  ! output print level (1-5          , fix)
    CALL dbg_print('omipBudOce:Tsurf ocean'        , Tsurf                 , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:atmflx%SWnetw'      , atmos_fluxes%SWnetw   , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:atmflx%LWnetw'      , atmos_fluxes%LWnetw   , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:atmflx%sensw'       , atmos_fluxes%sensw    , str_module, idt_src, in_subset=p_patch%cells%owned)
    CALL dbg_print('omipBudOce:atmflx%latw'        , atmos_fluxes%latw     , str_module, idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE calc_omip_budgets_oce

  !-------------------------------------------------------------------------
  !>
  !! Balance sea level to zero over global ocean
  !!
  !! Balance sea level to zero over global ocean
  !! This routine uses parts of mo_ocean_diagnostics
  !!
  !! @par Revision History
  !! Initial revision by Stephan Lorenz, MPI (2013-04)
  !!
  !!
  SUBROUTINE balance_elevation (p_patch_3D, h_old)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)    :: p_patch_3D
    REAL(wp), INTENT(INOUT)                 :: h_old(1:nproma,1:p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    TYPE(t_patch), POINTER                  :: p_patch
    TYPE(t_subset_range), POINTER           :: all_cells

    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: jc, jb
    REAL(wp) :: ocean_are, glob_slev, corr_slev

    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
 
    ! parallelize correctly
    ocean_are = p_patch_3D%p_patch_1D(1)%ocean_area(1)
    glob_slev = global_sum_array(p_patch%cells%area(:,:)*h_old(:,:)*p_patch_3D%wet_halo_zero_c(:,1,:))
    corr_slev = glob_slev/ocean_are

    idt_src=2
    IF ((my_process_is_stdio()) .AND. (idbg_mxmn >= idt_src)) &
      & write(0,*)' BALANCE_ELEVATION(Dom): ocean_are, glob_slev, corr_slev =',ocean_are, glob_slev, glob_slev/ocean_are

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc =  i_startidx_c, i_endidx_c
        IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN
          ! subtract or scale?
          h_old(jc,jb) = h_old(jc,jb) - corr_slev
          !h_old(jc,jb) = h_old(jc,jb) * (1.0_wp - corr_slev)
          !h_old(jc,jb) = h_old(jc,jb) - h_old(jc,jb)*corr_slev
        END IF
      END DO
    END DO

  END SUBROUTINE balance_elevation
 
 
  !-------------------------------------------------------------------------
  !>
  !! Read ocean forcing data from netcdf
  !!
  !! Read ocean forcing data for NCEP or other forcing
  !! This routine reads annual data sets of length forcing_timescale
  !!
  !! @par Revision History
  !! Initial revision by Stephan Lorenz, MPI (2012-02-17)
  !!
  !!
  SUBROUTINE read_forc_data_oce (p_patch, ext_data, no_set)

    TYPE(t_patch), INTENT(IN)            :: p_patch
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    INTEGER,       INTENT(IN)            :: no_set          !  no of set in file to be read

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ocean_surface:read_forc_data_oce'

    CHARACTER(filename_max) :: ncep_file   !< file name for reading in

    LOGICAL :: l_exist
    INTEGER :: jg, no_cells, no_tst, jtime, jt !, jc, jb
    INTEGER :: ncid, dimid,mpi_comm
    TYPE(t_stream_id) :: stream_id
    INTEGER :: i_start(2),i_count(2), jcells

    REAL(wp):: z_flux(nproma,p_patch%alloc_cell_blocks,forcing_timescale)  ! set length is forcing_timescale, 3rd dimension
    REAL(wp):: z_c   (nproma,forcing_timescale,p_patch%alloc_cell_blocks)  ! 2nd dimension is forcing_timescale
    !TYPE (t_keyword_list), POINTER :: keywords => NULL()

    !-------------------------------------------------------------------------

    !  READ NCEP FORCING

    !-------------------------------------------------------------------------

    !CALL message (TRIM(routine), 'start')

    IF (iforc_oce == OMIP_FluxFromFile) THEN

    !DO jg = 1,n_dom
      jg = 1

      ! i_lev       = p_patch%level
      ! WRITE (ncep_file,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',i_lev, '-flux.nc'
      ncep_file='ocean-flux.nc'

      !ncep_file=TRIM('/pool/data/ICON/external/iconR2B04-flux.nc')

      IF(my_process_is_stdio()) THEN
        !
        CALL message( TRIM(routine),'Ocean NCEP forcing flux file is: '//TRIM(ncep_file) )
        INQUIRE (FILE=ncep_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'NCEP forcing flux file is not found.')
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(ncep_file), NF_NOWRITE, ncid))
        !CALL message( TRIM(routine),'Ocean NCEP flux file opened for read' )

        !
        ! get and check number of cells in ncep data
        !
        CALL nf(nf_inq_dimid(ncid, 'ncells', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))

        IF(p_patch%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in NCEP flux file do not match.')
        ENDIF

        !
        ! get number of timesteps
        !
        CALL nf(nf_inq_dimid(ncid, 'time', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, no_tst))
        !
        ! check - s.b.

      ENDIF

      IF(my_process_is_stdio()) CALL nf(nf_close(ncid))
      stream_id = openInputFile(ncep_file, p_patch)

      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF
      CALL p_bcast(no_tst, p_io, mpi_comm)

      !-------------------------------------------------------
      !
      ! Read 12 monthly NCEP data sets for triangle centers using 4-dim routine
      !
      !-------------------------------------------------------

      jcells = p_patch%n_patch_cells  !  global dimension
      jtime  = forcing_timescale              !  time period to read (not yet)


      ! provide NCEP fluxes for sea ice (interface to ocean)
      ! 1:  'stress_x': zonal wind stress       [Pa]
      ! 2:  'stress_y': meridional wind stress  [Pa]
      ! 3:  'SST"     : sea surface temperature [K]

      ! zonal wind stress
      !write(0,*) ' ncep set 1: dimensions:',p_patch%n_patch_cells_g, p_patch%n_patch_cells, &
      ! &  forcing_timescale, nproma, p_patch%nblks_c
      !CALL read_3D(stream_id, on_cells, 'stress_x', z_flx2(:,:,:))
      !write(0,*) ' READ_FORC, READ 1: first data sets: stress-x, block=5, index=1,5:'
      !do jt=1,jtime
      !  write(0,*) 'jt=',jt,' val:',(z_flx2(jc,jt,5),jc=1,5)
      !enddo

      ! start-pointer and length of pointer for reading data:
      ! start: first set (1,1); second year (1,jtime+1)
      i_start(1) = 1
      i_start(2) = jtime*(no_set-1) + 1  ! position pointer to set no_set
      i_count(1) = jcells                ! length of pointer, dim 1 of z_dummy_array
      i_count(2) = jtime                 ! length of pointer, dim 2 of z_dummy_array

      idt_src=2  ! output print level (1-5, fix)
      IF (idbg_mxmn >= idt_src) THEN
        !
        WRITE(message_text,'(A,I6,A)')  'Ocean NCEP flux file contains',no_tst,' data sets'
        CALL message( TRIM(routine), TRIM(message_text) )

        WRITE(message_text,'(4(A,I4))')  'NCEP data set: length = ',jtime, &
          &   '; no. of set =',no_set,                                     &
          &   '; pos. of ptr =', i_start(2)
        CALL message( TRIM(routine), TRIM(message_text) )
      END IF

      CALL read_2D_time(stream_id, on_cells, 'stress_x', &
        &               fill_array=z_flux(:,:,:), start_timestep=i_start(2), &
        &               end_timestep=i_count(2) + i_start(2))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,1) = z_flux(:,:,jt)
      END DO

      ! meridional wind stress
      CALL read_2D_time(stream_id, on_cells, 'stress_y', &
        &               fill_array=z_flux(:,:,:), start_timestep=i_start(2), &
        &               end_timestep=i_count(2) + i_start(2))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,2) = z_flux(:,:,jt)
      END DO

      ! SST
      CALL read_2D_time(stream_id, on_cells, 'SST', &
        &               fill_array=z_flux(:,:,:), start_timestep=i_start(2), &
        &               end_timestep=i_count(2) + i_start(2))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,3) = z_flux(:,:,jt)
      END DO

 !    ! Read complete NCEP data sets for focing ocean model (iforc_type=5)
 !    ! 4:  tafo(:,:),   &  ! 2 m air temperature                              [C]
 !    ! 5:  ftdew(:,:),  &  ! 2 m dew-point temperature                        [K]
 !    ! 6:  fu10(:,:) ,  &  ! 10 m wind speed                                  [m/s]
 !    ! 7:  fclou(:,:),  &  ! Fractional cloud cover
 !    ! 8:  pao(:,:),    &  ! Surface atmospheric pressure                     [hPa]
 !    ! 9:  fswr(:,:),   &  ! Incoming surface solar radiation                 [W/m]

      ! 2m-temperature
      CALL read_2D_time(stream_id, on_cells, 'temp_2m', &
        &               fill_array=z_flux(:,:,:), start_timestep=i_start(2), &
        &               end_timestep=i_count(2) + i_start(2))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,4) = z_flux(:,:,jt)
      END DO

      ! 2m dewpoint temperature
      CALL read_2D_time(stream_id, on_cells, 'dpt_temp_2m', &
        &               fill_array=z_flux(:,:,:), start_timestep=i_start(2), &
        &               end_timestep=i_count(2) + i_start(2))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,5) = z_flux(:,:,jt)
      END DO

      ! Scalar wind
      CALL read_2D_time(stream_id, on_cells, 'scalar_wind', &
        &               fill_array=z_flux(:,:,:), start_timestep=i_start(2), &
        &               end_timestep=i_count(2) + i_start(2))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,6) = z_flux(:,:,jt)
      END DO

      ! cloud cover
      CALL read_2D_time(stream_id, on_cells, 'cloud', &
        &               fill_array=z_flux(:,:,:), start_timestep=i_start(2), &
        &               end_timestep=i_count(2) + i_start(2))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,7) = z_flux(:,:,jt)
      END DO

      ! sea level pressure
      CALL read_2D_time(stream_id, on_cells, 'pressure', &
        &               fill_array=z_flux(:,:,:), start_timestep=i_start(2), &
        &               end_timestep=i_count(2) + i_start(2))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,8) = z_flux(:,:,jt)
      END DO

      ! total solar radiation
      CALL read_2D_time(stream_id, on_cells, 'tot_solar', &
        &               fill_array=z_flux(:,:,:), start_timestep=i_start(2), &
        &               end_timestep=i_count(2) + i_start(2))
      DO jt = 1, jtime
        ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,9) = z_flux(:,:,jt)
      END DO

      ! precipitation
  !   CALL read_2D_time(stream_id, on_cells, 'precip', &
  !     &               fill_array=z_flux(:,:,:), start_timestep=i_start(2), &
  !     &               end_timestep=i_count(2) + i_start(2))
  !   DO jt = 1, jtime
  !     ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,10) = z_flux(:,:,jt)
  !   END DO

      ! evaporation or downward surface LW flux
  !   CALL read_2D_time(stream_id, on_cells, 'evap', &
  !     &               fill_array=z_flux(:,:,:), start_timestep=i_start(2), &
  !     &               end_timestep=i_count(2) + i_start(2))
  !   DO jt = 1, jtime
  !     ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,11) = z_flux(:,:,jt)
  !   END DO
  !   CALL read_2D_time(stream_id, on_cells, 'dlwrf', &
  !     &               fill_array=z_flux(:,:,:), start_timestep=i_start(2), &
  !     &               end_timestep=i_count(2) + i_start(2))
  !   DO jt = 1, jtime
  !     ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,11) = z_flux(:,:,jt)
  !   END DO

      ! runoff
  !   CALL read_2D_time(stream_id, on_cells, 'runoff', &
  !     &               fill_array=z_flux(:,:,:), start_timestep=i_start(2), &
  !     &               end_timestep=i_count(2) + i_start(2))
  !   DO jt = 1, jtime
  !     ext_data(jg)%oce%flux_forc_mon_c(:,jt,:,12) = z_flux(:,:,jt)
  !   END DO


      !
      ! close file
      !
      CALL closeFile(stream_id)

    !ENDDO

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,1)
      CALL dbg_print('ReadFc: NCEP: stress-x'    ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,2)
      CALL dbg_print('ReadFc: NCEP: stress-y'    ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,3)
      CALL dbg_print('ReadFc: NCEP: SST'         ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      idt_src=4  ! output print level (1-5, fix)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,4)
      CALL dbg_print('ReadFc: NCEP: temp_2m'     ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,5)
      CALL dbg_print('ReadFc: NCEP: dpt_temp_2m' ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,6)
      CALL dbg_print('ReadFc: NCEP: scalar_wind' ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,7)
      CALL dbg_print('ReadFc: NCEP: cloudiness'  ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,8)
      CALL dbg_print('ReadFc: NCEP: pressure'    ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,9)
      CALL dbg_print('ReadFc: NCEP: total solar' ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
    ! z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,10)
    ! CALL dbg_print('ReadFc: NCEP: precip.'     ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
    ! z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,11)
    ! CALL dbg_print('ReadFc: NCEP: evaporation' ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
    ! z_c(:,:,:) = ext_data(jg)%oce%flux_forc_mon_c(:,:,:,12)
    ! CALL dbg_print('ReadFc: NCEP: runoff'      ,z_c         ,str_module,idt_src, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------

      idt_src=2  ! output print level (1-5, fix)
      IF (idbg_mxmn >= idt_src) &
        & CALL message( TRIM(routine),'Ocean NCEP fluxes for external data read' )

    END IF ! iforc_oce=OMIP_FluxFromFile

  END SUBROUTINE read_forc_data_oce
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE nf(status)

    INTEGER, INTENT(in) :: status

    IF (status /= nf_noerr) THEN
      CALL finish('mo_ocean_surface netCDF error', nf_strerror(status))
    ENDIF

  END SUBROUTINE nf
  !-------------------------------------------------------------------------

END MODULE mo_ocean_surface
