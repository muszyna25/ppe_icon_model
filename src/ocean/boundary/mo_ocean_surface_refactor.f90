!>
!! Provide an implementation of the ocean surface module.
!!
!! Provide an implementation of the parameters used for surface forcing
!! of the hydrostatic ocean model.
!!
!! @author Vladimir Lapin, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
!!  Restructured code by Stephan Lorenz, MPI-M: (2015-04)
!!  Restructured code by Vladimir Lapin, MPI-M: (2017-01)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ocean_surface_refactor
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
  USE mo_dynamics_config,     ONLY: nold
  USE mo_exception,           ONLY: finish, message
  USE mo_util_dbg_prnt,       ONLY: dbg_print

  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

  USE mo_ocean_nml,           ONLY: iforc_oce, no_tracer, type_surfRelax_Temp, type_surfRelax_Salt, &
    &  No_Forcing, Analytical_Forcing, OMIP_FluxFromFile, Coupled_FluxFromAtmo,                     &
    &  i_sea_ice, zero_freshwater_flux, atmos_flux_analytical_type, atmos_precip_const, &  ! atmos_evap_constant
    &  limit_elevation, lhamocc, lswr_jerlov, lhamocc, lfb_bgc_oce, lcheck_salt_content, &
    &  lfix_salt_content, surface_flux_type, ice_flux_type

  USE mo_sea_ice_nml,        ONLY: sice

  USE mo_ocean_nml,           ONLY: atmos_flux_analytical_type, relax_analytical_type, &
    &  n_zlev, para_surfRelax_Salt, para_surfRelax_Temp, atmos_precip_const, &  ! atmos_evap_constant
    &  atmos_SWnet_const, atmos_LWnet_const, atmos_lat_const, atmos_sens_const, &
    &  atmos_SWnetw_const, atmos_LWnetw_const, atmos_latw_const, atmos_sensw_const, &
    &  relax_width, forcing_HeatFlux_amplitude, forcing_HeatFlux_base,              &
    &  basin_center_lat, basin_center_lon, basin_width_deg, basin_height_deg

  USE mo_math_constants,      ONLY: pi, deg2rad, rad2deg
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref, alv, tmelt, tf, clw, stbo, zemiss_def
  USE mo_impl_constants,      ONLY: max_char_length, sea_boundary, MIN_DOLIC

  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_ocean_surface_types, ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff

  USE mtime,                  ONLY: datetime, getNoOfSecondsElapsedInDayDateTime
  USE mo_ice_interface,       ONLY: ice_fast_interface, ice_thermodynamics, ice_dynamics
  USE mo_ocean_bulk_forcing,  ONLY: update_surface_relaxation, apply_surface_relaxation, &
                                &   update_flux_fromFile, calc_omip_budgets_ice, calc_omip_budgets_oce, &
                                &   update_ocean_surface_stress, balance_elevation, &
                                & update_surface_relaxation_zstar, balance_elevation_zstar

  USE mo_ocean_diagnostics, ONLY : diag_heat_salt_tendency
  USE mo_name_list_output_init, ONLY: isRegistered

  USE mo_swr_absorption
  USE mo_ocean_check_total_content,       ONLY: check_total_salt_content, check_total_si_volume

  USE mo_mpi, only: get_my_mpi_work_id 

  IMPLICIT NONE
  
  PRIVATE

  ! public interface
  PUBLIC  :: update_ocean_surface_refactor
  ! private routine
  PRIVATE :: apply_surface_fluxes
  PUBLIC  :: apply_surface_fluxes_slo
  PUBLIC  :: update_atmos_fluxes
  PRIVATE :: update_atmos_fluxes_analytical
  PUBLIC  :: update_ocean_surface_refactor_zstar

  CHARACTER(len=12)           :: str_module    = 'OceanSurfaceRefactor'  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Update ocean surface by applying flux forcing for hydrostatic ocean
  !!
  !! This function changes:
  !! p_ice      thermodynamical and dynamical fields of sea ice
  !! p_oce_sfc  surface fluxes and stress, passed to the ocean
  !! p_os       SSH, SST, SSS and HAMMOC tracers (dilution)
  !!
  !! @par Revision History
  !! Initial release (mo_oce_bulk)          by Stephan Lorenz, MPI-M (2010-07)
  !! Restructuring 1 (mo_ocean_surface)     by Stephan Lorenz, MPI-M (2015-04)
  !! Restructuring 2 (mo_ocean_surface_refactor) by Vladimir Lapin, MPI-M (2017-02)
  !
  SUBROUTINE update_ocean_surface_refactor(p_patch_3D, p_os, p_as, p_ice, atmos_fluxes, p_oce_sfc, this_datetime, p_op_coeff)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_atmos_for_ocean)                     :: p_as
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    TYPE(datetime), POINTER                     :: this_datetime
    TYPE(t_operator_coeff),   INTENT(IN)        :: p_op_coeff
    !
    ! local variables
    TYPE(t_patch), POINTER                      :: p_patch
    INTEGER                                     :: trac_no
    REAL(wp)                                    :: dsec

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_surface_refactor:update_ocean_surface'

    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    
    p_oce_sfc%top_dilution_coeff(:,:) = 1.0_wp ! Initilize dilution factor to one

    IF (no_tracer>=1) p_oce_sfc%sst => p_os%p_prog(nold(1))%tracer(:,1,:,1)
    IF (no_tracer>=2) p_oce_sfc%sss => p_os%p_prog(nold(1))%tracer(:,1,:,2)

    IF (iforc_oce == No_Forcing) RETURN  !  forcing for ocean not defined

    ! save values of ice, snow and temperature for the tendendy and hfbasin diagnostic
    IF ( isRegistered('delta_ice') .OR. isRegistered('delta_snow') .OR. &
           isRegistered('delta_thetao') .OR. &
           isRegistered('delta_so') .OR. &
           isRegistered('global_sltbasin') .OR. isRegistered('atlant_sltbasin') .OR. &
           isRegistered('pacind_sltbasin') .OR. &
           isRegistered('global_hfbasin') .OR. isRegistered('atlant_hfbasin') .OR. &
           isRegistered('pacind_hfbasin') ) THEN

      CALL diag_heat_salt_tendency(p_patch_3d, 1, p_ice,            &
         p_os%p_prog(nold(1))%tracer(:,:,:,1),               &
         p_os%p_prog(nold(1))%tracer(:,:,:,2),               &
         p_os%p_diag%delta_ice,                              &
         p_os%p_diag%delta_snow, p_os%p_diag%delta_thetao, p_os%p_diag%delta_so)

     END IF

    !---------------------------------------------------------------------
    ! (1) Apply relaxation to surface temperature and salinity
    !---------------------------------------------------------------------
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

!    IF (lcheck_salt_content) CALL check_total_salt_content(106,p_os%p_prog(nold(1))%tracer(:,:,:,2), p_patch, &
!         p_os%p_prog(nold(1))%h(:,:), p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
!         p_ice, 0)
    !---------------------------------------------------------------------
    ! (2) Receive/calculate surface fluxes and wind stress
    !---------------------------------------------------------------------
    CALL update_atmos_fluxes(p_patch_3D, p_as, atmos_fluxes, p_oce_sfc, p_os, p_ice, this_datetime)

    ! copy atmospheric wind speed from p_as%fu10 into new forcing variable for output purpose - not accumulated yet
    p_oce_sfc%Wind_Speed_10m(:,:) = p_as%fu10(:,:)

    !---------------------------------------------------------------------
    ! (3) Sea ice thermodynamics & dynamics (at ocean time-step)
    !---------------------------------------------------------------------
    p_oce_sfc%cellThicknessUnderIce(:,:) = p_ice%zUnderIce(:,:) ! neccessary, because is not yet in restart
   
    p_ice%draftave_old(:,:) = p_ice%draftave(:,:)

    IF ( i_sea_ice > 0 ) THEN ! sea ice is on

        !  (3a) Fast sea ice thermodynamics (Analytical or OMIP cases only. Otherwise, done in the atmosphere)
        IF (iforc_oce == Analytical_Forcing .OR. iforc_oce == OMIP_FluxFromFile)  THEN
            CALL ice_fast_interface(p_patch, p_ice, atmos_fluxes, this_datetime)
        ENDIF

!    IF (lcheck_salt_content)  CALL check_total_salt_content(1066,p_os%p_prog(nold(1))%tracer(:,:,:,2), p_patch, &
!           p_os%p_prog(nold(1))%h(:,:), p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
!           p_ice, 0)
 
!       CALL check_total_si_volume(1066,  p_patch, p_ice, p_os%p_prog(nold(1))%h(:,:))

        !  (3b) Slow sea ice dynamics and thermodynamics
        CALL ice_dynamics(p_patch_3D, p_ice, p_oce_sfc, atmos_fluxes, p_os, p_as, p_op_coeff)

!       CALL check_total_si_volume(1067,  p_patch, p_ice, p_os%p_prog(nold(1))%h(:,:))

!        IF (lcheck_salt_content) CALL check_total_salt_content(1067,p_os%p_prog(nold(1))%tracer(:,:,:,2), p_patch, &
!               p_os%p_prog(nold(1))%h(:,:), p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
!               p_ice, 1)

       IF (lfix_salt_content) CALL close_salt_budget(1,  p_patch_3D, p_os, p_ice, p_oce_sfc)


        CALL ice_thermodynamics(p_patch_3D, p_ice, p_oce_sfc, atmos_fluxes, p_os, p_as, p_op_coeff)

    ELSE !  sea ice is off

      ! for the setup without sea ice the SST is set to freezing temperature Tf
      ! should not be done here! Move to apply_surface_fluxes
      WHERE (p_oce_sfc%SST(:,:) .LT. Tf)
        p_oce_sfc%SST(:,:) = Tf
      ENDWHERE

    ENDIF
!    IF (lcheck_salt_content) CALL check_total_salt_content(107,p_os%p_prog(nold(1))%tracer(:,:,:,2), p_patch, &
!         p_os%p_prog(nold(1))%h(:,:), p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
!         p_ice, 1)

    !---------------------------------------------------------------------
    ! (4) Ocean surface stress boundary condition (atm-ocean + ice-ocean)
    !---------------------------------------------------------------------
    ! atm-ocean stress is either from file, bulk formula, or coupling
    ! if ice dynamics is off, ocean-ice stress is set to zero (no friction with ice)
    CALL update_ocean_surface_stress(p_patch_3D, p_ice, p_os, atmos_fluxes, p_oce_sfc)

    !---------------------------------------------------------------------
    ! (5) Apply thermal and haline fluxes to the ocean surface layer
    !---------------------------------------------------------------------
!    CALL apply_surface_fluxes(p_patch_3D, p_os, p_ice, p_oce_sfc)

!   calculate the sw flux used for subsurface heating

!   include hamoccs chlorophylls effect sw absorption 
    IF ( lhamocc .AND. lfb_bgc_oce ) CALL dynamic_swr_absorption(p_patch_3d, p_os)


    IF ( lswr_jerlov ) THEN
     p_os%p_diag%heatabs(:,:)=(p_os%p_diag%swsum(:,:)  &
             *p_oce_sfc%HeatFlux_ShortWave(:,:)*(1.0_wp-p_ice%concsum(:,:)))

    ELSE
      p_os%p_diag%heatabs(:,:)=0.0_wp
    ENDIF

!    IF (lcheck_salt_content) CALL check_total_salt_content(108,p_os%p_prog(nold(1))%tracer(:,:,:,2), p_patch, &
!         p_os%p_prog(nold(1))%h(:,:), p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
!         p_ice, 0)

    IF (surface_flux_type .EQ. 2 ) CALL apply_surface_fluxes(p_patch_3D, p_os, p_ice, p_oce_sfc)
    IF (surface_flux_type .EQ. 1 ) CALL apply_surface_fluxes_slo(p_patch_3D, p_os, p_ice, p_oce_sfc)
    
    IF (surface_flux_type .EQ. 3 ) CALL apply_surface_fluxes_slo_3(p_patch_3D, p_os, p_ice, p_oce_sfc)

    IF (surface_flux_type .EQ. 10) CALL apply_surface_fluxes_v10(p_patch_3D, p_os, p_ice, p_oce_sfc)

!   IF (lcheck_salt_content) CALL check_total_salt_content(1080,p_os%p_prog(nold(1))%tracer(:,:,:,2), p_patch, &
!         p_os%p_prog(nold(1))%h(:,:), p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
!         p_ice, 0)

    IF (lfix_salt_content) CALL close_salt_budget(2, p_patch_3D, p_os, p_ice, p_oce_sfc)

    IF (lcheck_salt_content) CALL check_total_salt_content(109,p_os%p_prog(nold(1))%tracer(:,:,:,2), p_patch, &
        p_os%p_prog(nold(1))%h(:,:), p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
         p_ice, 0)

!   apply subsurface heating
    IF ( lswr_jerlov ) THEN

      CALL subsurface_swr_absorption(p_patch_3d, p_os)

    ENDIF


    !---------------------------------------------------------------------
    ! (6) Apply volume flux correction
    !---------------------------------------------------------------------
    !  - sea level is balanced to zero over ocean surface
    !  - correction applied daily
    !  calculate time
    dsec  = REAL(getNoOfSecondsElapsedInDayDateTime(this_datetime), wp)
    ! event at end of first timestep of day - tbd: use mtime
!     IF (limit_elevation .AND. (dsec-dtime)<0.1 ) THEN
!     IF (limit_elevation .AND. dsec < dtime ) THEN
    IF (limit_elevation) THEN

      CALL balance_elevation(p_patch_3D, p_os%p_prog(nold(1))%h,p_oce_sfc,p_ice)
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('UpdSfc: h-old+BalElev',p_os%p_prog(nold(1))%h  ,str_module, 3, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------
    END IF

  END SUBROUTINE update_ocean_surface_refactor

  !-------------------------------------------------------------------------
  !>
  !! Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
  !!
  !! This function changes:
  !! SSH, SST and SSS, zUnderIce (in p_os, p_oce_sfc and p_ice)
  !! HAMMOC tracers by dilution (in p_os)
  !!
  !! @par Revision History
  !! Initial release by Vladimir Lapin, MPI-M (2017-02)
  !
  SUBROUTINE apply_surface_fluxes(p_patch_3D, p_os, p_ice, p_oce_sfc)

    TYPE(t_patch_3D ),TARGET,   INTENT(IN)      :: p_patch_3D
    TYPE(t_hydro_ocean_state),  INTENT(INOUT)   :: p_os
    TYPE(t_sea_ice),            INTENT(INOUT)   :: p_ice
    TYPE(t_ocean_surface),      INTENT(INOUT)   :: p_oce_sfc
    !
    ! local variables
    TYPE(t_patch), POINTER          :: p_patch
    TYPE(t_subset_range), POINTER   :: all_cells
    INTEGER                         :: jc, jb, i_startidx_c, i_endidx_c, i_bgc_tra

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_surface_refactor:apply_surface_fluxes'

    REAL(wp)  :: heatflux_surface_layer ! heatflux into the surface layer
    REAL(wp)  :: zunderice_old

    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------

    ! Provide total freshwater volume forcing:
    p_oce_sfc%FrshFlux_VolumeTotal(:,:) = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_VolumeIce (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)

    ! Note: old freeboard is stored in p_oce_sfc%cellThicknessUnderIce (equiv. to zUnderIceIni in apply_surface_fluxes_slo)

    IF (no_tracer > 0) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, jc, i_bgc_tra, heatflux_surface_layer) SCHEDULE(dynamic)
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN
            ! (5a) Net heat flux changes sst using old freeboard (Thermodynamic Eq. 1)

            ! substract the fraction of heatflux used for subsurface heating
            heatflux_surface_layer=p_oce_sfc%HeatFlux_Total(jc,jb)-p_os%p_diag%heatabs(jc,jb)
            zunderice_old =  p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) &
              &                    + p_os%p_prog(nold(1))%h(jc,jb) - p_ice%draftave_old(jc,jb)

            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    heatflux_surface_layer*dtime/(clw*rho_ref*zunderice_old)

            ! (5b) Net volume flux (plus snow fall on ice) changes ssh (Thermodynamic Eq. 5)
            p_os%p_prog(nold(1))%h(jc,jb) = p_os%p_prog(nold(1))%h(jc,jb)               &
              &                           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
              &                           + p_ice%totalsnowfall(jc,jb)

            ! (5c) New zUnderIce calculated with NEW h and draftave
            p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) &
              &                    + p_os%p_prog(nold(1))%h(jc,jb) - p_ice%draftave(jc,jb)

            ! (5d) New salinity is calculated from conservation formula:
            !      SSS_new * zUnderIce_new = ( SSS_old * zUnderIceOld + SaltFluxFromIce * dtime )
            p_oce_sfc%sss(jc,jb)   = ( p_oce_sfc%sss(jc,jb) * zunderice_old + &
             &                         p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime ) / p_ice%zUnderIce(jc,jb)


            !! update cell thickness under ice in p_oce_sfc
            p_oce_sfc%cellThicknessUnderIce(jc,jb) = p_ice%zUnderIce(jc,jb)

          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
    END IF

    ! provide total salinity forcing flux for diagnostics only
    p_oce_sfc%FrshFlux_TotalSalt(:,:)   = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalIce  (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('UpdSfc: oce_sfc%VolTot', p_oce_sfc%FrshFlux_VolumeTotal, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: oce_sfc%TotIce', p_oce_sfc%FrshFlux_TotalIce,    str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIce   ',   p_ice%zUnderIce,                str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST,                  str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SSS ',p_oce_sfc%SSS,                  str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEnd: h-old+fwfVol',p_os%p_prog(nold(1))%h,         str_module, 2, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE apply_surface_fluxes


  !-------------------------------------------------------------------------
  !>
  !! Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
  !!
  !! @par Revision History
  !! Initial release (mo_oce_bulk)          by Stephan Lorenz, MPI-M (2010-07)
  !! Restructuring (mo_ocean_surface)       by Stephan Lorenz, MPI-M (2015-04)
  !
  SUBROUTINE apply_surface_fluxes_slo(p_patch_3D, p_os, p_ice, p_oce_sfc)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    !
    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_surface_refactor:apply_surface_fluxes_slo'
    INTEGER               :: jc, jb, trac_no
    INTEGER               :: i_startidx_c, i_endidx_c
    REAL(wp)              :: sss_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zUnderIceOld(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zUnderIceIni(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zUnderIceArt(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    REAL(wp)  :: heatflux_surface_layer ! heatflux into the surface layer
    REAL(wp)  :: zunderice_ini


    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells

    
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------
    sss_inter(:,:)  = p_oce_sfc%sss(:,:)
    zUnderIceOld(:,:) = 0.0_wp
    zUnderIceArt(:,:) = 0.0_wp
    ! freeboard before sea ice model (used for thermal boundary condition (Eq.1))
    ! by construction, is stored in p_oce_sfc%cellThicknessUnderIce
    zUnderIceIni(:,:) = p_oce_sfc%cellThicknessUnderIce (:,:)


    !!  Provide total ocean forcing:
    !    - total heat fluxes are aggregated for ice/ocean in ice thermodynamics
    !    - total internal salt flux p_oce_sfc%FrshFlux_TotalIce is calculated in sea ice model
    !    - total freshwater volume forcing
    p_oce_sfc%FrshFlux_VolumeTotal(:,:) = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_VolumeIce (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)
    ! provide total salinity forcing flux for diagnostics only
    p_oce_sfc%FrshFlux_TotalSalt(:,:)   = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalIce  (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)

    !  ******  (Thermodynamic Eq. 1)  ******
    ! Apply net surface heat flux to ocean surface (new p_oce_flx%SST)
    IF (no_tracer > 0) THEN

      ! sst-change in surface module after sea-ice thermodynamics using HeatFlux_Total and old freeboard zUnderIceIni
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN

            ! substract the fraction of heatflux used for subsurface heating
            heatflux_surface_layer=p_oce_sfc%HeatFlux_Total(jc,jb)-p_os%p_diag%heatabs(jc,jb)
            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    heatflux_surface_layer*dtime/(clw*rho_ref*zUnderIceIni(jc,jb))

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
          !  #slo# 2015-01: totalsnowfall is needed for correct salt update (in surface module)
          !                 since draft was increased by snowfall but water below ice is not affected by snowfall
          !                 snow to ice conversion does not effect draft
          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) + &
          &                      p_os%p_prog(nold(1))%h(jc,jb) - p_ice%draftave(jc,jb) + p_ice%totalsnowfall(jc,jb)

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

         zUnderIce_ini=  p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) &
              &                    + p_os%p_prog(nold(1))%h(jc,jb) - p_ice%draftave_old(jc,jb)
          !******  (Thermodynamic Eq. 5)  ******
          !! Finally, let sea-level change from P-E+RO plus snow fall on ice, net total volume forcing to ocean surface
          p_os%p_prog(nold(1))%h(jc,jb) = p_os%p_prog(nold(1))%h(jc,jb)               &
            &                           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
            &                           + p_ice%totalsnowfall(jc,jb)

          !! update zunderice
          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) + p_os%p_prog(nold(1))%h(jc,jb) &
            &                    - p_ice%draftave(jc,jb)
          
          p_oce_sfc%top_dilution_coeff(jc,jb) = zUnderIce_ini / p_ice%zUnderIce(jc,jb)
        ENDIF  !  dolic>0
      END DO
    END DO

    !! set correct cell thickness under ice
    p_oce_sfc%cellThicknessUnderIce   (:,:) = p_ice%zUnderIce(:,:)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('UpdSfc: oce_sfc%HFTot ', p_oce_sfc%HeatFlux_Total,       str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: oce_sfc%VolTot', p_oce_sfc%FrshFlux_VolumeTotal, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: oce_sfc%TotIce', p_oce_sfc%FrshFlux_TotalIce,    str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: ice%totalsnowf', p_ice%totalsnowfall,            str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIceIni',   zUnderIceIni,                   str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIceArt',   zUnderIceArt,                   str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIceOld',   zUnderIceOld,                   str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIce   ',   p_ice%zUnderIce,                str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: sss_inter   ',   sss_inter,                      str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST,                  str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SSS ',p_oce_sfc%SSS,                  str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEnd: h-old+fwfVol',p_os%p_prog(nold(1))%h,         str_module, 2, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------


  END SUBROUTINE apply_surface_fluxes_slo


  !-------------------------------------------------------------------------
  !>
  !! Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
  !!
  !! @par Revision History
  !! Initial release (mo_oce_bulk)          by Stephan Lorenz, MPI-M (2010-07)
  !! Restructuring (mo_ocean_surface)       by Stephan Lorenz, MPI-M (2015-04)
  !
  SUBROUTINE apply_surface_fluxes_slo_3(p_patch_3D, p_os, p_ice, p_oce_sfc)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    !
    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_surface_refactor:apply_surface_fluxes_slo'
    INTEGER               :: jc, jb, trac_no
    INTEGER               :: i_startidx_c, i_endidx_c
    REAL(wp)              :: sss_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zUnderIce_old
    REAL(wp)              :: zUnderIceIni(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zUnderIceArt(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    REAL(wp)  :: heatflux_surface_layer ! heatflux into the surface layer
    REAL(wp)  :: zunderice_ini


    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells

    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------

    ! Provide total freshwater volume forcing:
    p_oce_sfc%FrshFlux_VolumeTotal(:,:) = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_VolumeIce (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)

    ! Note: old freeboard is stored in p_oce_sfc%cellThicknessUnderIce (equiv. to zUnderIceIni in apply_surface_fluxes_slo)

    IF (no_tracer > 0) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, jc, i_bgc_tra, heatflux_surface_layer) SCHEDULE(dynamic)
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN
            ! (5a) Net heat flux changes sst using old freeboard (Thermodynamic Eq. 1)

            ! substract the fraction of heatflux used for subsurface heating
            heatflux_surface_layer=p_oce_sfc%HeatFlux_Total(jc,jb)-p_os%p_diag%heatabs(jc,jb)
            zunderice_old =  p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) &
              &                    + p_os%p_prog(nold(1))%h(jc,jb) 

            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    heatflux_surface_layer*dtime/(clw*rho_ref*zunderice_old)

            ! (5b) Net volume flux (plus snow fall on ice) changes ssh (Thermodynamic Eq. 5)
            p_os%p_prog(nold(1))%h(jc,jb) = p_os%p_prog(nold(1))%h(jc,jb)               &
              &                           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime & 
              &                           + p_oce_sfc%FrshFlux_TotalIce(jc, jb)*dtime  

            ! (5c) New zUnderIce calculated with NEW h and draftave
            p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) &
              &                    + p_os%p_prog(nold(1))%h(jc,jb) 

            ! (5d) New salinity is calculated from conservation formula:
            !      SSS_new * zUnderIce_new = ( SSS_old * zUnderIceOld + SaltFluxFromIce * dtime )
            p_oce_sfc%sss(jc,jb)   = ( p_oce_sfc%sss(jc,jb) * zunderice_old + &
             &                         p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime ) / p_ice%zUnderIce(jc,jb)


            !! update cell thickness under ice in p_oce_sfc
            p_oce_sfc%cellThicknessUnderIce(jc,jb) = p_ice%zUnderIce(jc,jb)

          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
    END IF

    ! provide total salinity forcing flux for diagnostics only
    p_oce_sfc%FrshFlux_TotalSalt(:,:)   = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalIce  (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)
   

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('UpdSfc: oce_sfc%HFTot ', p_oce_sfc%HeatFlux_Total,       str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: oce_sfc%VolTot', p_oce_sfc%FrshFlux_VolumeTotal, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: oce_sfc%TotIce', p_oce_sfc%FrshFlux_TotalIce,    str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: ice%totalsnowf', p_ice%totalsnowfall,            str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIce   ',   p_ice%zUnderIce,                str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: sss_inter   ',   sss_inter,                      str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST,                  str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SSS ',p_oce_sfc%SSS,                  str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEnd: h-old+fwfVol',p_os%p_prog(nold(1))%h,         str_module, 2, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------


  END SUBROUTINE apply_surface_fluxes_slo_3


  !-------------------------------------------------------------------------
  !>
  !! Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
  !!
  !
  SUBROUTINE apply_surface_fluxes_v10(p_patch_3D, p_os, p_ice, p_oce_sfc)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    !
    ! local variables
    INTEGER               :: jc, jb
    INTEGER               :: i_startidx_c, i_endidx_c
    REAL(wp)              :: sss_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: sst_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zunderice_old
    REAL(wp)              :: zUnderIceIni(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    REAL(wp) :: h_old_test, h_new_test
    
    REAL(wp) :: temp_eta, min_h
    REAL(wp) :: temp_stretch(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    
    INTEGER  :: bt_lev, jk, adj_lev
    REAL(wp) :: d_c, dz_old, dz_new, extra_salt_in_col
    REAL(wp) :: dz_ratio, adj_lev_ht
    
    REAL(wp) :: old_sss, new_s1, new_s2
    INTEGER  :: flag 

    REAL(wp)  :: heatflux_surface_layer ! heatflux into the surface layer
    
    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    
    CHARACTER(LEN=max_char_length), PARAMETER :: str_module = 'ocean_sfc_refactor_zstar'
    
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------
    sst_inter(:,:)    = p_oce_sfc%sst(:,:)
    sss_inter(:,:)    = p_oce_sfc%sss(:,:)
    zUnderIceIni(:,:) = p_oce_sfc%cellThicknessUnderIce (:,:)

    CALL dbg_print('UpdSfcSTART: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcSTART: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)

    !!  Provide total ocean forcing:
    !    - total heat fluxes are aggregated for ice/ocean in ice thermodynamics
    !    - total internal salt flux p_oce_sfc%FrshFlux_TotalIce is calculated in sea ice model
    !    - total freshwater volume forcing
    p_oce_sfc%FrshFlux_VolumeTotal(:,:) = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_VolumeIce (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:) &
      &                                 + p_oce_sfc%FrshFlux_Relax     (:,:)
    ! provide total salinity forcing flux for diagnostics only
    p_oce_sfc%FrshFlux_TotalSalt(:,:)   = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalIce  (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)

    
    ! Apply net surface heat flux to ocean surface (new p_oce_flx%SST)
    IF (no_tracer > 0) THEN

      ! sst-change in surface module after sea-ice thermodynamics using HeatFlux_Total and old freeboard zUnderIceIni
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN

            ! substract the fraction of heatflux used for subsurface heating
            heatflux_surface_layer=p_oce_sfc%HeatFlux_Total(jc,jb)-p_os%p_diag%heatabs(jc,jb)
            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    heatflux_surface_layer*dtime/(clw*rho_ref*zUnderIceIni(jc,jb))

          ENDIF
        ENDDO
      ENDDO

    END IF

    ! apply volume flux to surface elevation
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

      DO jc = i_startidx_c, i_endidx_c
        IF (p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0) THEN

          zunderice_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) &
            &                    + p_os%p_prog(nold(1))%h(jc,jb) 

          p_os%p_prog(nold(1))%h(jc,jb) = p_os%p_prog(nold(1))%h(jc,jb)               &
            &                           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
            &                           + p_oce_sfc%FrshFlux_TotalIce(jc, jb)*dtime  

          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) &
            &                    + p_os%p_prog(nold(1))%h(jc,jb) 

          new_s1   = ( p_oce_sfc%sss(jc,jb) * zunderice_old + &
           &                         p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime ) / p_ice%zUnderIce(jc,jb)

          p_oce_sfc%sss(jc,jb) = new_s1

          p_oce_sfc%cellThicknessUnderIce(jc,jb) = p_ice%zUnderIce(jc,jb)

        ENDIF  !  dolic>0
      END DO
    END DO
          
    !! set correct cell thickness under ice
    p_oce_sfc%cellThicknessUnderIce   (:,:) = p_ice%zUnderIce(:,:)

    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND v10: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 1, in_subset=p_patch%cells%owned)


  END SUBROUTINE apply_surface_fluxes_v10


  SUBROUTINE close_salt_budget(step, p_patch_3D, p_os, p_ice, p_oce_sfc)

    TYPE(t_patch_3D ),TARGET,   INTENT(IN)      :: p_patch_3D
    TYPE(t_hydro_ocean_state),  INTENT(INOUT)   :: p_os
    TYPE(t_sea_ice),            INTENT(IN)   :: p_ice
    TYPE(t_ocean_surface),      INTENT(INOUT)   :: p_oce_sfc

    REAL(wp)     :: zUnderIce_old
    !
    ! local variables
    TYPE(t_patch), POINTER          :: p_patch
    TYPE(t_subset_range), POINTER   :: all_cells
    INTEGER                         :: step, jc, jb, i_startidx_c, i_endidx_c



    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------


    IF (no_tracer > 0) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, jc,zunderice_old) SCHEDULE(dynamic)
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN
            IF (step .EQ. 1 ) THEN

              zUnderIce_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) &
                   + p_os%p_prog(nold(1))%h(jc,jb) - p_ice%draftave_old(jc,jb)

              ! new salt from ice after the ice advection ;  liquid water part  before ice advection
              p_oce_sfc%surface_salt_content(jc,jb) =sice*rhoi/rho_ref*SUM(p_ice%hi(jc,:,jb)*p_ice%conc(jc,:,jb)) &
                   +p_os%p_prog(nold(1))%tracer(jc,1,jb,2)*zUnderIce_old

            ENDIF

            IF (step .EQ. 2 ) THEN

              p_ice%draft(jc,:,jb) = (rhos * p_ice%hs(jc,:,jb) + rhoi * p_ice%hi(jc,:,jb))/rho_ref
              p_ice%draftave(jc,jb) = SUM(p_ice%draft(jc,:,jb) * p_ice%conc(jc,:,jb))
              p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) &
                   + p_os%p_prog(nold(1))%h(jc,jb) - p_ice%draftave(jc,jb)

              p_oce_sfc%sss(jc,jb)   = ( p_oce_sfc%surface_salt_content(jc,jb) - sice * rhoi/rho_ref * &
                 & SUM(p_ice%hi(jc,:,jb)*p_ice%conc(jc,:,jb))) &
                 &  /p_ice%zUnderIce(jc,jb)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
    END IF


  END SUBROUTINE close_salt_budget

  
  !-------------------------------------------------------------------------
  !>
  !! Update ocean surface by applying flux forcing for hydrostatic ocean
  !!
  !! This function changes:
  !! p_ice      thermodynamical and dynamical fields of sea ice
  !! p_oce_sfc  surface fluxes and stress, passed to the ocean
  !! p_os       SSH, SST, SSS and HAMMOC tracers (dilution)
  !!
  !  Adapted for zstar
  !
  SUBROUTINE update_ocean_surface_refactor_zstar(p_patch_3D, p_os, p_as, p_ice, &
      & atmos_fluxes, p_oce_sfc, this_datetime, p_op_coeff, eta_c, stretch_c)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_atmos_for_ocean)                     :: p_as
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    TYPE(datetime), POINTER                     :: this_datetime
    TYPE(t_operator_coeff),   INTENT(IN)        :: p_op_coeff
    REAL(wp), INTENT(INOUT) :: eta_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(IN   ) :: stretch_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    !
    ! local variables
    TYPE(t_patch), POINTER                      :: p_patch
    INTEGER                                     :: trac_no
    REAL(wp)                                    :: dsec

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_surface_refactor:update_ocean_surface'

    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------

    IF (no_tracer>=1) p_oce_sfc%sst => p_os%p_prog(nold(1))%tracer(:,1,:,1)
    IF (no_tracer>=2) p_oce_sfc%sss => p_os%p_prog(nold(1))%tracer(:,1,:,2)

    IF (iforc_oce == No_Forcing) RETURN  !  forcing for ocean not defined

    ! save values of ice, snow and temperature for the tendendy and hfbasin diagnostic
    IF ( isRegistered('delta_ice') .OR. isRegistered('delta_snow') .OR. &
           isRegistered('delta_thetao') .OR. &
           isRegistered('delta_so') .OR. &
           isRegistered('global_sltbasin') .OR. isRegistered('atlant_sltbasin') .OR. &
           isRegistered('pacind_sltbasin') .OR. &
           isRegistered('global_hfbasin') .OR. isRegistered('atlant_hfbasin') .OR. &
           isRegistered('pacind_hfbasin') ) THEN

      CALL diag_heat_salt_tendency(p_patch_3d, 1, p_ice,            &
         p_os%p_prog(nold(1))%tracer(:,:,:,1),               &
         p_os%p_prog(nold(1))%tracer(:,:,:,2),               &
         p_os%p_diag%delta_ice,                              &
         p_os%p_diag%delta_snow, p_os%p_diag%delta_thetao, p_os%p_diag%delta_so)

     END IF

    !---------------------------------------------------------------------
    ! (1) Apply relaxation to surface temperature and salinity
    !---------------------------------------------------------------------
    IF (type_surfRelax_Temp >= 1) THEN
      trac_no = 1   !  tracer no 1: temperature
      CALL update_surface_relaxation_zstar(p_patch_3D, p_os, p_ice, p_oce_sfc, trac_no, stretch_c)

      !  apply restoring to surface temperature directly
      CALL apply_surface_relaxation(p_patch_3D, p_os, p_oce_sfc, trac_no)

    END IF

    IF (type_surfRelax_Salt >= 1 .AND. no_tracer >1) THEN
      trac_no = 2   !  tracer no 2: salinity
      CALL update_surface_relaxation_zstar(p_patch_3D, p_os, p_ice, p_oce_sfc, trac_no, stretch_c)

      !  apply restoring to surface salinity directly
      CALL apply_surface_relaxation(p_patch_3D, p_os, p_oce_sfc, trac_no)

    ENDIF

    !---------------------------------------------------------------------
    ! (2) Receive/calculate surface fluxes and wind stress
    !---------------------------------------------------------------------
    CALL update_atmos_fluxes(p_patch_3D, p_as, atmos_fluxes, p_oce_sfc, p_os, p_ice, this_datetime)

    ! copy atmospheric wind speed from p_as%fu10 into new forcing variable for output purpose - not accumulated yet
    p_oce_sfc%Wind_Speed_10m(:,:) = p_as%fu10(:,:)

    !---------------------------------------------------------------------
    ! (3) Sea ice thermodynamics & dynamics (at ocean time-step)
    !---------------------------------------------------------------------
    p_oce_sfc%cellThicknessUnderIce(:,:) = p_ice%zUnderIce(:,:) ! neccessary, because is not yet in restart

    IF ( i_sea_ice > 0 ) THEN ! sea ice is on

        !  (3a) Fast sea ice thermodynamics (Analytical or OMIP cases only. Otherwise, done in the atmosphere)
        IF (iforc_oce == Analytical_Forcing .OR. iforc_oce == OMIP_FluxFromFile)  THEN
            CALL ice_fast_interface(p_patch, p_ice, atmos_fluxes, this_datetime)
        ENDIF

        CALL ice_dynamics(p_patch_3D, p_ice, p_oce_sfc, atmos_fluxes, p_os, p_as, p_op_coeff)

        CALL ice_thermodynamics(p_patch_3D, p_ice, p_oce_sfc, atmos_fluxes, p_os, p_as, p_op_coeff)

    ELSE !  sea ice is off

      ! for the setup without sea ice the SST is set to freezing temperature Tf
      ! should not be done here! Move to apply_surface_fluxes
      WHERE (p_oce_sfc%SST(:,:) .LT. Tf)
        p_oce_sfc%SST(:,:) = Tf
      ENDWHERE

    ENDIF

    !---------------------------------------------------------------------
    ! (4) Ocean surface stress boundary condition (atm-ocean + ice-ocean)
    !---------------------------------------------------------------------
    ! atm-ocean stress is either from file, bulk formula, or coupling
    ! if ice dynamics is off, ocean-ice stress is set to zero (no friction with ice)
    CALL update_ocean_surface_stress(p_patch_3D, p_ice, p_os, atmos_fluxes, p_oce_sfc)

    !---------------------------------------------------------------------
    ! (5) Apply thermal and haline fluxes to the ocean surface layer
    !---------------------------------------------------------------------
!   calculate the sw flux used for subsurface heating
!   include hamoccs chlorophylls effect sw absorption 
!   FIXME zstar: Haven't checked for zstar requirements in hamocc 
    IF ( lhamocc .AND. lfb_bgc_oce ) CALL dynamic_swr_absorption(p_patch_3d, p_os)


    IF ( lswr_jerlov ) THEN
     p_os%p_diag%heatabs(:,:)=(p_os%p_diag%swsum(:,:)  &
             *p_oce_sfc%HeatFlux_ShortWave(:,:)*(1.0_wp-p_ice%concsum(:,:)))

    ELSE
      p_os%p_diag%heatabs(:,:)=0.0_wp
    ENDIF

    IF (surface_flux_type .EQ. 2 ) CALL apply_surface_fluxes_zstar(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)

    IF (surface_flux_type .EQ. 1 ) CALL apply_surface_fluxes_slo_zstar(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)
    
    IF (surface_flux_type .EQ. 10) CALL apply_surface_fluxes_zstar_v10(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)
    IF (surface_flux_type .EQ. 11) CALL apply_surface_fluxes_zstar_v11(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)
    IF (surface_flux_type .EQ. 12) CALL apply_surface_fluxes_zstar_v12(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)
    IF (surface_flux_type .EQ. 14) CALL apply_surface_fluxes_zstar_v14(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)
    IF (surface_flux_type .EQ. 15) CALL apply_surface_fluxes_zstar_v15(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)

!   apply subsurface heating
    IF ( lswr_jerlov ) THEN

      CALL subsurface_swr_absorption_zstar(p_patch_3d, p_os, stretch_c)

    ENDIF


    !---------------------------------------------------------------------
    ! (6) Apply volume flux correction
    !---------------------------------------------------------------------
    !  - sea level is balanced to zero over ocean surface
    !  - correction applied daily
    !  calculate time
    dsec  = REAL(getNoOfSecondsElapsedInDayDateTime(this_datetime), wp)
    ! event at end of first timestep of day - tbd: use mtime
    IF (limit_elevation .AND. (dsec-dtime)<0.1 ) THEN
      CALL balance_elevation_zstar(p_patch_3D, eta_c)
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('UpdSfc: h-old+BalElev', eta_c, routine, 2, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------
    END IF

  END SUBROUTINE update_ocean_surface_refactor_zstar



  !-------------------------------------------------------------------------
  !>
  !! Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
  !!
  !! Adapted for zstar
  !
  SUBROUTINE apply_surface_fluxes_slo_zstar(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    !
    REAL(wp), INTENT(INOUT) :: eta_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(IN   ) :: stretch_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
 
    ! local variables
    INTEGER               :: jc, jb
    INTEGER               :: i_startidx_c, i_endidx_c
    REAL(wp)              :: sss_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zUnderIceOld(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zUnderIceIni(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zUnderIceArt(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    REAL(wp) :: h_old_test, h_new_test
    
    REAL(wp) :: temp_eta, min_h
    REAL(wp) :: temp_stretch(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    
    INTEGER  :: bt_lev
    REAL(wp) :: d_c 

    REAL(wp)  :: heatflux_surface_layer ! heatflux into the surface layer
    
    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    
    CHARACTER(LEN=max_char_length), PARAMETER :: str_module = 'ocean_sfc_refactor_zstar'
    
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------
    sss_inter(:,:)  = p_oce_sfc%sss(:,:)
    zUnderIceOld(:,:) = 0.0_wp
    zUnderIceArt(:,:) = 0.0_wp
    ! freeboard before sea ice model (used for thermal boundary condition (Eq.1))
    ! by construction, is stored in p_oce_sfc%cellThicknessUnderIce
    zUnderIceIni(:,:) = p_oce_sfc%cellThicknessUnderIce (:,:)


    !!  Provide total ocean forcing:
    !    - total heat fluxes are aggregated for ice/ocean in ice thermodynamics
    !    - total internal salt flux p_oce_sfc%FrshFlux_TotalIce is calculated in sea ice model
    !    - total freshwater volume forcing
    p_oce_sfc%FrshFlux_VolumeTotal(:,:) = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_VolumeIce (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:) &
      &                                 + p_oce_sfc%FrshFlux_Relax     (:,:)
    ! provide total salinity forcing flux for diagnostics only
    p_oce_sfc%FrshFlux_TotalSalt(:,:)   = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalIce  (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)

    !  ******  (Thermodynamic Eq. 1)  ******
    ! Apply net surface heat flux to ocean surface (new p_oce_flx%SST)
    IF (no_tracer > 0) THEN

      ! sst-change in surface module after sea-ice thermodynamics using HeatFlux_Total and old freeboard zUnderIceIni
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN

            ! substract the fraction of heatflux used for subsurface heating
            heatflux_surface_layer=p_oce_sfc%HeatFlux_Total(jc,jb)-p_os%p_diag%heatabs(jc,jb)
            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    heatflux_surface_layer*dtime/(clw*rho_ref*zUnderIceIni(jc,jb))

          ENDIF
        ENDDO
      ENDDO

    END IF

    CALL dbg_print('UpdSfc: eta-old', eta_c,    str_module, 1, in_subset=p_patch%cells%owned)
    ! apply volume flux to surface elevation
    !  - add to h_old before explicit term
    !  - change in salt concentration applied here
    !    i.e. for salinity relaxation only, no volume flux is applied
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

      DO jc = i_startidx_c, i_endidx_c
        IF (p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0) THEN

          bt_lev = p_patch_3d%p_patch_1d(1)%dolic_c(jc, jb)      
          d_c    = p_patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, jb)


          !******  (Thermodynamic Eq. 2)  ******
          !! Calculate the new freeboard caused by changes in ice thermodynamics
          !!  zUnderIce = z_surf + h_old - (z_draft - z_snowfall)
          !  #slo# 2015-01: totalsnowfall is needed for correct salt update (in surface module)
          !                 since draft was increased by snowfall but water below ice is not affected by snowfall
          !                 snow to ice conversion does not effect draft
          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * &
          &                      stretch_c(jc,jb) - p_ice%draftave(jc,jb) + p_ice%totalsnowfall(jc,jb)


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

          h_old_test =  p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * stretch_c(jc,jb)

          !******  (Thermodynamic Eq. 5)  ******
          !! Finally, let sea-level change from P-E+RO plus snow fall on ice, net total volume forcing to ocean surface
          temp_eta     = eta_c(jc,jb)              

          eta_c(jc,jb) = eta_c(jc,jb)               &
            &           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
            &           + p_ice%totalsnowfall(jc,jb)

          !! Only change the stretching parameter if it is above a certain threshold
          !! This avoids divide by 0 
          temp_stretch(jc, jb) = stretch_c(jc, jb)
          min_h                = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)
          
          !! Update only if height is atleast dz
          if ( d_c  .GT.  min_h ) &
            & temp_stretch(jc, jb) = ( eta_c(jc, jb) + d_c)/( d_c )
 
          !! update zunderice
          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * temp_stretch(jc, jb) &
            &                    - p_ice%draftave(jc,jb)
    
          h_new_test =  p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)*temp_stretch(jc, jb)
          p_oce_sfc%top_dilution_coeff(jc,jb) = h_old_test/h_new_test
          
        ENDIF  !  dolic>0
      END DO
    END DO
          
    !! set correct cell thickness under ice
    p_oce_sfc%cellThicknessUnderIce   (:,:) = p_ice%zUnderIce(:,:)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('UpdSfc: oce_sfc%HFTot ', p_oce_sfc%HeatFlux_Total,       str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: oce_sfc%VolTot', p_oce_sfc%FrshFlux_VolumeTotal, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: oce_sfc%TotIce', p_oce_sfc%FrshFlux_TotalIce,    str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: ice%totalsnowf', p_ice%totalsnowfall,            str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIceIni',   zUnderIceIni,                   str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIceArt',   zUnderIceArt,                   str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIceOld',   zUnderIceOld,                   str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIce   ',   p_ice%zUnderIce,                str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: sss_inter   ',   sss_inter,                      str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST,                  str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SSS ',p_oce_sfc%SSS,                  str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEnd: h-old+fwfVol',p_os%p_prog(nold(1))%h,         str_module, 2, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------


  END SUBROUTINE apply_surface_fluxes_slo_zstar

 
  !-------------------------------------------------------------------------
  !>
  !! Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
  !!
  !! Adapted for zstar
  !
  SUBROUTINE apply_surface_fluxes_zstar(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    !
    REAL(wp), INTENT(INOUT) :: eta_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(IN   ) :: stretch_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
 
    ! local variables
    INTEGER               :: jc, jb
    INTEGER               :: i_startidx_c, i_endidx_c
    REAL(wp)              :: sss_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: sst_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zunderice_old
    REAL(wp)              :: zUnderIceIni(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    REAL(wp) :: h_old_test, h_new_test
    
    REAL(wp) :: temp_eta, min_h
    REAL(wp) :: temp_stretch(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    
    INTEGER  :: bt_lev, jk, adj_lev
    REAL(wp) :: d_c, dz_old, dz_new, extra_salt_in_col
    REAL(wp) :: dz_ratio, adj_lev_ht
    
    REAL(wp) :: old_sss, new_s1, new_s2
    INTEGER  :: flag 

    REAL(wp)  :: heatflux_surface_layer ! heatflux into the surface layer
    
    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    
    CHARACTER(LEN=max_char_length), PARAMETER :: str_module = 'ocean_sfc_refactor_zstar'
    
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------
    sst_inter(:,:)    = p_oce_sfc%sst(:,:)
    sss_inter(:,:)    = p_oce_sfc%sss(:,:)
    zUnderIceIni(:,:) = p_oce_sfc%cellThicknessUnderIce (:,:)

    CALL dbg_print('UpdSfcSTART: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcSTART: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)

    !! For thin levels we need to adjust surface flux over multiple levels
    !! below we compute number of levels based on arbitrary requirement
    !! of the level that includes 11m depth
    adj_lev    = 1
    adj_lev_ht = p_patch_3d%p_patch_1d(1)%zlev_i(2)
    DO jk = 1, n_zlev 
      IF ( ( p_patch_3d%p_patch_1d(1)%zlev_i(jk) .LE. 11) .AND. &
        & ( p_patch_3d%p_patch_1d(1)%zlev_i(jk + 1) .GT. 11) ) THEN 
        
         adj_lev    = jk
         adj_lev_ht = p_patch_3d%p_patch_1d(1)%zlev_i(jk + 1)

      END IF
    END DO

    !!  Provide total ocean forcing:
    !    - total heat fluxes are aggregated for ice/ocean in ice thermodynamics
    !    - total internal salt flux p_oce_sfc%FrshFlux_TotalIce is calculated in sea ice model
    !    - total freshwater volume forcing
    p_oce_sfc%FrshFlux_VolumeTotal(:,:) = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_VolumeIce (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:) &
      &                                 + p_oce_sfc%FrshFlux_Relax     (:,:)
    ! provide total salinity forcing flux for diagnostics only
    p_oce_sfc%FrshFlux_TotalSalt(:,:)   = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalIce  (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)

    
!    IF (get_my_mpi_work_id() .EQ. 359 ) THEN
!      write(0, *) "====================================================" 
!      write(0, *) "THIN LEVEL ALERT" 
!      write(0, *) p_ice%draftave(4, 1), p_oce_sfc%FrshFlux_VolumeTotal(4,1)*dtime
!      write(0, *) sst_inter(4, 1), p_oce_sfc%sst(4, 1)
!      write(0, *) stretch_c(4, 1) 
!      write(0, *) eta_c(4, 1), p_oce_sfc%FrshFlux_TotalIce(4, 1)*dtime  
!      write(0, *) p_os%p_diag%kin(4, 1, 1), p_os%p_diag%kin(4, 2, 1), p_os%p_diag%kin(4, 4, 1)
!      write(0, *) "====================================================" 
!    END IF


    !  ******  (Thermodynamic Eq. 1)  ******
    ! Apply net surface heat flux to ocean surface (new p_oce_flx%SST)
    IF (no_tracer > 0) THEN

      ! sst-change in surface module after sea-ice thermodynamics using HeatFlux_Total and old freeboard zUnderIceIni
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN

            ! substract the fraction of heatflux used for subsurface heating
            heatflux_surface_layer=p_oce_sfc%HeatFlux_Total(jc,jb)-p_os%p_diag%heatabs(jc,jb)
            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    heatflux_surface_layer*dtime/(clw*rho_ref*zUnderIceIni(jc,jb))

!            IF (p_oce_sfc%sst(jc,jb) > 40 ) THEN
!              write(0, *) "====================================================" 
!              write(0, *) "Hoch temperature" 
!              write(0, *) jc, jb,  p_patch_3d%p_patch_1d(1)%dolic_c(jc, jb), get_my_mpi_work_id()
!              write(0, *) sst_inter(jc,jb), p_oce_sfc%sst(jc,jb)
!              write(0, *) p_oce_sfc%HeatFlux_Total(jc,jb), p_os%p_diag%heatabs(jc,jb)
!              write(0, *) heatflux_surface_layer, clw, rho_ref, zUnderIceIni(jc,jb)
!              write(0, *) p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb), stretch_c(jc,jb) 
!            END IF

          ENDIF
        ENDDO
      ENDDO

    END IF

    CALL dbg_print('UpdSfc: eta-old', eta_c,    str_module, 1, in_subset=p_patch%cells%owned)
    ! apply volume flux to surface elevation
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

      DO jc = i_startidx_c, i_endidx_c
        IF (p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0) THEN

          bt_lev = p_patch_3d%p_patch_1d(1)%dolic_c(jc, jb)      
          d_c    = p_patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, jb)

          !! We are using variable name zUnderIce but here it only refers
          !! to top level thickness
          zunderice_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * &
            &                      stretch_c(jc,jb) 


          h_old_test =  p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * stretch_c(jc,jb)

          temp_eta     = eta_c(jc,jb)              

          eta_c(jc,jb) = eta_c(jc,jb)               &
            &           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
            &           + p_oce_sfc%FrshFlux_TotalIce(jc, jb)*dtime  

          !! Only change the stretching parameter if it is above a certain threshold
          !! This avoids divide by 0 
          temp_stretch(jc, jb) = stretch_c(jc, jb)
          min_h                = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)
          
          !! Update only if height is atleast dz
          if ( d_c  .GT.  min_h ) &
            & temp_stretch(jc, jb) = ( eta_c(jc, jb) + d_c)/( d_c )
 
          !! update zunderice
          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * temp_stretch(jc, jb) 
    
!          p_oce_sfc%sss(jc,jb)   = ( sss_inter(jc,jb) * zunderice_old + &
!             &          p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime ) / p_ice%zUnderIce(jc,jb)

          !! Calculate extra salt due to changed layers 
          extra_salt_in_col = 0.0
          DO jk = adj_lev + 1, bt_lev 
            dz_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
              & * stretch_c(jc, jb) 
            dz_new = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
              & * temp_stretch(jc, jb) 
            
            extra_salt_in_col = extra_salt_in_col + &
              & p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) * (dz_old - dz_new)

          END DO
          
          flag = 0
          !! FIXME Test this for salt conservation
          !! Distribute surface fluxes over multiple levels depending on
          !! thickness ratio 
          DO jk = 1, adj_lev
            dz_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
              & * stretch_c(jc, jb) 
            dz_new = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
              & * temp_stretch(jc, jb) 

            dz_ratio = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)/adj_lev_ht

            ! Distribute over 4 levels
!            p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) = p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) &
!              &  + ( p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime ) /dz_new
            !! 
            new_s1  = p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2)*dz_old/dz_new &
              &  +  dz_ratio*p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime/dz_new &
              &  +  dz_ratio*extra_salt_in_col/dz_new
 
            new_s2  = p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2)*dz_old/dz_new &
              &  +  dz_ratio*p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime/dz_new  

            IF (new_s1 .LT. 0) flag = 1 
            IF (new_s2 .LT. 0) flag = 1 

          END DO  

          IF (FLAG .EQ. 0) THEN
          
            DO jk = 1, adj_lev
              dz_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
                & * stretch_c(jc, jb) 
              dz_new = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
                & * temp_stretch(jc, jb) 
  
              dz_ratio = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)/adj_lev_ht
  
              new_s1  = p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2)*dz_old/dz_new &
                &  +  dz_ratio*p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime/dz_new &
                &  +  dz_ratio*extra_salt_in_col/dz_new
   
              p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) = new_s1
              
            END DO  
          ELSE

            DO jk = 1, bt_lev
              dz_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
                & * stretch_c(jc, jb) 
              dz_new = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
                & * temp_stretch(jc, jb) 
  
              new_s1  = p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2)*dz_old/dz_new  
  
              IF (jk .EQ. 1) THEN
                new_s1 = new_s1  +  p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime/dz_new  
              END IF
   
              p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) = new_s1
              
            END DO  

          END IF


          h_new_test =  p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)*temp_stretch(jc, jb)
          p_oce_sfc%top_dilution_coeff(jc,jb) = h_old_test/h_new_test
          
        ENDIF  !  dolic>0
      END DO
    END DO
          
    !! set correct cell thickness under ice
    p_oce_sfc%cellThicknessUnderIce   (:,:) = p_ice%zUnderIce(:,:)

    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)


  END SUBROUTINE apply_surface_fluxes_zstar


  !-------------------------------------------------------------------------
  !>
  !! Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
  !!
  !! Adapted for zstar
  !
  SUBROUTINE apply_surface_fluxes_zstar_v10(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    !
    REAL(wp), INTENT(INOUT) :: eta_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(IN   ) :: stretch_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
 
    ! local variables
    INTEGER               :: jc, jb
    INTEGER               :: i_startidx_c, i_endidx_c
    REAL(wp)              :: sss_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: sst_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zunderice_old
    REAL(wp)              :: zUnderIceIni(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    REAL(wp) :: h_old_test, h_new_test
    
    REAL(wp) :: temp_eta, min_h
    REAL(wp) :: temp_stretch(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    
    INTEGER  :: bt_lev, jk, adj_lev
    REAL(wp) :: d_c, dz_old, dz_new, extra_salt_in_col
    REAL(wp) :: dz_ratio, adj_lev_ht
    
    REAL(wp) :: old_sss, new_s1, new_s2
    INTEGER  :: flag 

    REAL(wp)  :: heatflux_surface_layer ! heatflux into the surface layer
    
    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    
    CHARACTER(LEN=max_char_length), PARAMETER :: str_module = 'ocean_sfc_refactor_zstar'
    
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------
    sst_inter(:,:)    = p_oce_sfc%sst(:,:)
    sss_inter(:,:)    = p_oce_sfc%sss(:,:)
    zUnderIceIni(:,:) = p_oce_sfc%cellThicknessUnderIce (:,:)

    CALL dbg_print('UpdSfcSTART: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcSTART: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)

    !!  Provide total ocean forcing:
    !    - total heat fluxes are aggregated for ice/ocean in ice thermodynamics
    !    - total internal salt flux p_oce_sfc%FrshFlux_TotalIce is calculated in sea ice model
    !    - total freshwater volume forcing
    p_oce_sfc%FrshFlux_VolumeTotal(:,:) = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_VolumeIce (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:) &
      &                                 + p_oce_sfc%FrshFlux_Relax     (:,:)
    ! provide total salinity forcing flux for diagnostics only
    p_oce_sfc%FrshFlux_TotalSalt(:,:)   = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalIce  (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)

    
    ! Apply net surface heat flux to ocean surface (new p_oce_flx%SST)
    IF (no_tracer > 0) THEN

      ! sst-change in surface module after sea-ice thermodynamics using HeatFlux_Total and old freeboard zUnderIceIni
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN

            ! substract the fraction of heatflux used for subsurface heating
            heatflux_surface_layer=p_oce_sfc%HeatFlux_Total(jc,jb)-p_os%p_diag%heatabs(jc,jb)
            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    heatflux_surface_layer*dtime/(clw*rho_ref*zUnderIceIni(jc,jb))

          ENDIF
        ENDDO
      ENDDO

    END IF

    ! apply volume flux to surface elevation
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

      DO jc = i_startidx_c, i_endidx_c
        IF (p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0) THEN

          zunderice_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * &
            &                      stretch_c(jc,jb) 


          h_old_test =  p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * stretch_c(jc,jb)

          eta_c(jc,jb) = eta_c(jc,jb)               &
            &           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
            &           + p_oce_sfc%FrshFlux_TotalIce(jc, jb)*dtime  

          !! Only change the stretching parameter if it is above a certain threshold
          !! This avoids divide by 0 
          temp_stretch(jc, jb) = stretch_c(jc, jb)
          min_h                = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)
          
          !! Update only if height is atleast dz
          if ( d_c  .GT.  min_h ) &
            & temp_stretch(jc, jb) = ( eta_c(jc, jb) + d_c)/( d_c )
 
          !! update zunderice
          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * temp_stretch(jc, jb) 

          new_s1   = ( p_oce_sfc%sss(jc,jb) * zunderice_old + &
           &                         p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime ) / p_ice%zUnderIce(jc,jb)

          p_oce_sfc%sss(jc,jb) = new_s1

          p_oce_sfc%cellThicknessUnderIce(jc,jb) = p_ice%zUnderIce(jc,jb)

        ENDIF  !  dolic>0
      END DO
    END DO
          
    !! set correct cell thickness under ice
    p_oce_sfc%cellThicknessUnderIce   (:,:) = p_ice%zUnderIce(:,:)

    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND v10_zstar: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 1, in_subset=p_patch%cells%owned)


  END SUBROUTINE apply_surface_fluxes_zstar_v10

  !-------------------------------------------------------------------------
  !>
  !! Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
  !!
  !! Adapted for zstar
  !
  SUBROUTINE apply_surface_fluxes_zstar_v11(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    !
    REAL(wp), INTENT(INOUT) :: eta_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(IN   ) :: stretch_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
 
    ! local variables
    INTEGER               :: jc, jb
    INTEGER               :: i_startidx_c, i_endidx_c
    REAL(wp)              :: sss_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: sst_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zunderice_old
    REAL(wp)              :: zUnderIceIni(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    REAL(wp) :: h_old_test, h_new_test
    
    REAL(wp) :: temp_eta, min_h
    REAL(wp) :: temp_stretch(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    
    INTEGER  :: bt_lev, jk, adj_lev
    REAL(wp) :: d_c, dz_old, dz_new, extra_salt_in_col
    REAL(wp) :: dz_ratio, adj_lev_ht
    
    REAL(wp) :: old_sss, new_s1, new_s2, new_sss, thresh_sss
    INTEGER  :: flag 

    REAL(wp)  :: heatflux_surface_layer ! heatflux into the surface layer
    
    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    
    CHARACTER(LEN=max_char_length), PARAMETER :: str_module = 'ocean_sfc_refactor_zstar'
    
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------
    sst_inter(:,:)    = p_oce_sfc%sst(:,:)
    sss_inter(:,:)    = p_oce_sfc%sss(:,:)
    zUnderIceIni(:,:) = p_oce_sfc%cellThicknessUnderIce (:,:)

    CALL dbg_print('UpdSfcSTART: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcSTART: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)

    !!  Provide total ocean forcing:
    !    - total heat fluxes are aggregated for ice/ocean in ice thermodynamics
    !    - total internal salt flux p_oce_sfc%FrshFlux_TotalIce is calculated in sea ice model
    !    - total freshwater volume forcing
    p_oce_sfc%FrshFlux_VolumeTotal(:,:) = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_VolumeIce (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:) &
      &                                 + p_oce_sfc%FrshFlux_Relax     (:,:)
    ! provide total salinity forcing flux for diagnostics only
    p_oce_sfc%FrshFlux_TotalSalt(:,:)   = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalIce  (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)

    
    IF (no_tracer > 0) THEN

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN

            ! substract the fraction of heatflux used for subsurface heating
            heatflux_surface_layer=p_oce_sfc%HeatFlux_Total(jc,jb)-p_os%p_diag%heatabs(jc,jb)
            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    heatflux_surface_layer*dtime/(clw*rho_ref*zUnderIceIni(jc,jb))

          ENDIF
        ENDDO
      ENDDO

    END IF

    CALL dbg_print('UpdSfc: eta-old', eta_c,    str_module, 1, in_subset=p_patch%cells%owned)
    ! apply volume flux to surface elevation
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

      DO jc = i_startidx_c, i_endidx_c
        IF (p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0) THEN

          bt_lev = p_patch_3d%p_patch_1d(1)%dolic_c(jc, jb)      
          d_c    = p_patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, jb)

          !! We are using variable name zUnderIce but here it only refers
          !! to top level thickness
          zunderice_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * &
            &                      stretch_c(jc,jb) 

          temp_eta     = eta_c(jc,jb)              

          eta_c(jc,jb) = eta_c(jc,jb)               &
            &           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
            &           + p_oce_sfc%FrshFlux_TotalIce(jc, jb)*dtime  

          !! Only change the stretching parameter if it is above a certain threshold
          !! This avoids divide by 0 
          temp_stretch(jc, jb) = stretch_c(jc, jb)
          min_h                = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)
          
          !! Update only if height is atleast dz
          if ( d_c  .GT.  min_h ) &
            & temp_stretch(jc, jb) = ( eta_c(jc, jb) + d_c)/( d_c )
 
          !! update zunderice
          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * temp_stretch(jc, jb) 
    
          !! Calculate extra salt due to changed layers 
          extra_salt_in_col = 0.0
          DO jk = 2, bt_lev 
            dz_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
              & * stretch_c(jc, jb) 
            dz_new = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
              & * temp_stretch(jc, jb) 
            
            extra_salt_in_col = extra_salt_in_col + &
              & p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) * (dz_old - dz_new)

          END DO
          
          new_s1   = ( p_oce_sfc%sss(jc,jb) * zunderice_old  &
            &             +  p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime &
            &             +  extra_salt_in_col )/ p_ice%zUnderIce(jc,jb)

          !! Threshold for minimum SSS
          thresh_sss = 0.1_wp 

          flag = 0
          IF (new_s1 .LT. thresh_sss) THEN
            new_sss = thresh_sss
            flag    = 1
            
            extra_salt_in_col = ( new_s1 - new_sss )*p_ice%zUnderIce(jc,jb)

          ELSE
            new_sss = new_s1
          END IF

          !! Restrict SSS to threshold value to prevent negative salinity
          !! This will not be salt conservative but should occur rarely
          p_oce_sfc%sss(jc,jb) = new_sss

!          IF (flag .EQ. 1) THEN
!          
!            dz_new = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, 2, jb)&
!              & * temp_stretch(jc, jb) 
!            
!            new_s1   = p_os%p_prog(nold(1))%tracer(jc, 2, jb, 2) &
!              &             +  extra_salt_in_col /dz_new
!
!            IF (new_s1 .LT. thresh_sss) THEN
!              new_sss = thresh_sss
!              flag    = 2
!              
!              extra_salt_in_col = extra_salt_in_col + ( new_s1 - new_sss)*dz_new
!  
!            ELSE
!              new_sss = new_s1
!            END IF
!
!            p_os%p_prog(nold(1))%tracer(jc,  2, jb, 2) = new_sss
!
!          END IF !! FLAG = 1
!
!          IF (flag .EQ. 2) THEN
!          
!            dz_new = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, 3, jb)&
!              & * temp_stretch(jc, jb) 
!            
!            new_s1   = p_os%p_prog(nold(1))%tracer(jc,  3, jb, 2) &
!              &             +  extra_salt_in_col /dz_new
!
!            !! Limit salinity to thresh_sss
!            !! should not be required so need to check
!            IF (new_s1 .LT. thresh_sss) THEN
!              write(0, *) "=============99999999===================" 
!              write(0, *) bt_lev, extra_salt_in_col, dz_new, p_ice%zUnderIce(jc,jb)
!              write(0, *) "===============HERE=====================" 
!              new_sss = thresh_sss
!            ELSE
!              new_sss = new_s1
!            END IF
!
!            p_os%p_prog(nold(1))%tracer(jc,  3, jb, 2) = new_sss
!
!
!          END IF !! FLAG = 2


        ENDIF  !  dolic>0
      END DO
    END DO
          
    !! set correct cell thickness under ice
    p_oce_sfc%cellThicknessUnderIce   (:,:) = p_ice%zUnderIce(:,:)

    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)

  END SUBROUTINE apply_surface_fluxes_zstar_v11


  !-------------------------------------------------------------------------
  !>
  !! Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
  !!
  !! Adapted for zstar
  !
  SUBROUTINE apply_surface_fluxes_zstar_v12(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    !
    REAL(wp), INTENT(INOUT) :: eta_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(IN   ) :: stretch_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
 
    ! local variables
    INTEGER               :: jc, jb
    INTEGER               :: i_startidx_c, i_endidx_c
    REAL(wp)              :: sss_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: sst_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zunderice_old
    REAL(wp)              :: zUnderIceIni(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    REAL(wp) :: h_old_test, h_new_test
    
    REAL(wp) :: temp_eta, min_h
    REAL(wp) :: temp_stretch(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    
    INTEGER  :: bt_lev, jk, adj_lev
    REAL(wp) :: d_c, dz_old, dz_new, extra_salt_in_col
    REAL(wp) :: dz_ratio, adj_lev_ht
    
    REAL(wp) :: old_sss, new_s1, new_s2, new_sss, thresh_sss
    INTEGER  :: flag 

    REAL(wp)  :: heatflux_surface_layer ! heatflux into the surface layer
    
    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    
    CHARACTER(LEN=max_char_length), PARAMETER :: str_module = 'ocean_sfc_refactor_zstar'
    
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------
    sst_inter(:,:)    = p_oce_sfc%sst(:,:)
    sss_inter(:,:)    = p_oce_sfc%sss(:,:)
    zUnderIceIni(:,:) = p_oce_sfc%cellThicknessUnderIce (:,:)

    CALL dbg_print('UpdSfcSTART: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcSTART: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)

    !!  Provide total ocean forcing:
    !    - total heat fluxes are aggregated for ice/ocean in ice thermodynamics
    !    - total internal salt flux p_oce_sfc%FrshFlux_TotalIce is calculated in sea ice model
    !    - total freshwater volume forcing
    p_oce_sfc%FrshFlux_VolumeTotal(:,:) = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_VolumeIce (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:) &
      &                                 + p_oce_sfc%FrshFlux_Relax     (:,:)
    ! provide total salinity forcing flux for diagnostics only
    p_oce_sfc%FrshFlux_TotalSalt(:,:)   = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalIce  (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)

    
    IF (no_tracer > 0) THEN

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN

            ! substract the fraction of heatflux used for subsurface heating
            heatflux_surface_layer=p_oce_sfc%HeatFlux_Total(jc,jb)-p_os%p_diag%heatabs(jc,jb)
            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    heatflux_surface_layer*dtime/(clw*rho_ref*zUnderIceIni(jc,jb))

          ENDIF
        ENDDO
      ENDDO

    END IF

    CALL dbg_print('UpdSfc: eta-old', eta_c,    str_module, 1, in_subset=p_patch%cells%owned)
    ! apply volume flux to surface elevation
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

      DO jc = i_startidx_c, i_endidx_c
        IF (p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0) THEN

          bt_lev = p_patch_3d%p_patch_1d(1)%dolic_c(jc, jb)      
          d_c    = p_patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, jb)

          !! We are using variable name zUnderIce but here it only refers
          !! to top level thickness
          zunderice_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * &
            &                      stretch_c(jc,jb) 

          temp_eta     = eta_c(jc,jb)              

          eta_c(jc,jb) = eta_c(jc,jb)               &
            &           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
            &           + p_oce_sfc%FrshFlux_TotalIce(jc, jb)*dtime  

          !! Only change the stretching parameter if it is above a certain threshold
          !! This avoids divide by 0 
          temp_stretch(jc, jb) = stretch_c(jc, jb)
          min_h                = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)
          
          !! Update only if height is atleast dz
          if ( d_c  .GT.  min_h ) &
            & temp_stretch(jc, jb) = ( eta_c(jc, jb) + d_c)/( d_c )
 
          !! update zunderice
          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * temp_stretch(jc, jb) 
    
          !! Calculate extra salt due to changed layers 
          extra_salt_in_col = 0.0
          DO jk = 2, bt_lev 
            dz_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
              & * stretch_c(jc, jb) 
            dz_new = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
              & * temp_stretch(jc, jb) 
            
            extra_salt_in_col = extra_salt_in_col + &
              & p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) * (dz_old - dz_new)

          END DO
          
          !! SSS if only the first layer were to be updated
          new_s1   = ( p_oce_sfc%sss(jc,jb) * zunderice_old  &
            &             +  p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime &
            &             +  extra_salt_in_col )/ p_ice%zUnderIce(jc,jb)

          thresh_sss = 0.5_wp

          !! If new_s1 starts approaching 0, distribute salinity over all layers
          IF (new_s1 > thresh_sss) THEN

            p_oce_sfc%sss(jc,jb) = new_s1  

          ELSE

!            write(0, *) "=============99999999===================" 
!            write(0, *) bt_lev, extra_salt_in_col, p_oce_sfc%sss(jc,jb), &
!              & p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime
!!            write(0, *) bt_lev, extra_salt_in_col, p_ice%zUnderIce(jc,jb)
!            write(0, *) "===============HERE=====================" 

            DO jk = 1, bt_lev 
              dz_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
                & * stretch_c(jc, jb) 
              dz_new = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
                & * temp_stretch(jc, jb) 
              
              p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) = & 
                & p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) * (dz_old / dz_new)
  
              if (jk .EQ. 1) THEN
                p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) = p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) &
                  &  +  p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime/dz_new  
              END IF
  
            END DO
          END IF !! choose salinity update
 
        ENDIF  !  dolic>0
      END DO
    END DO
          
    !! set correct cell thickness under ice
    p_oce_sfc%cellThicknessUnderIce   (:,:) = p_ice%zUnderIce(:,:)

    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)

  END SUBROUTINE apply_surface_fluxes_zstar_v12


  
  !>
  !! Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
  !!
  !! Scale salinity and temperature.  
  !
  SUBROUTINE apply_surface_fluxes_zstar_v14(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    !
    REAL(wp), INTENT(INOUT) :: eta_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(IN   ) :: stretch_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
 
    ! local variables
    INTEGER               :: jc, jb
    INTEGER               :: i_startidx_c, i_endidx_c
    REAL(wp)              :: sss_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: sst_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zunderice_old
    REAL(wp)              :: zUnderIceIni(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    REAL(wp) :: h_old_test, h_new_test
    
    REAL(wp) :: temp_eta, min_h
    REAL(wp) :: temp_stretch(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    
    INTEGER  :: bt_lev, jk, adj_lev
    REAL(wp) :: d_c, dz_old, dz_new, extra_salt_in_col
    REAL(wp) :: dz_ratio, adj_lev_ht
    
    REAL(wp) :: old_sss, new_s1, new_s2, new_sss, thresh_sss
    INTEGER  :: flag 

    REAL(wp)  :: heatflux_surface_layer ! heatflux into the surface layer
    
    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    
    CHARACTER(LEN=max_char_length), PARAMETER :: str_module = 'ocean_sfc_refactor_zstar'
    
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------
    sst_inter(:,:)    = p_oce_sfc%sst(:,:)
    sss_inter(:,:)    = p_oce_sfc%sss(:,:)
    zUnderIceIni(:,:) = p_oce_sfc%cellThicknessUnderIce (:,:)

    CALL dbg_print('UpdSfcSTART: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcSTART: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)

    !!  Provide total ocean forcing:
    !    - total heat fluxes are aggregated for ice/ocean in ice thermodynamics
    !    - total internal salt flux p_oce_sfc%FrshFlux_TotalIce is calculated in sea ice model
    !    - total freshwater volume forcing
    p_oce_sfc%FrshFlux_VolumeTotal(:,:) = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_VolumeIce (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:) &
      &                                 + p_oce_sfc%FrshFlux_Relax     (:,:)
    ! provide total salinity forcing flux for diagnostics only
    p_oce_sfc%FrshFlux_TotalSalt(:,:)   = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalIce  (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)

    
    CALL dbg_print('UpdSfc: eta-old', eta_c,    str_module, 1, in_subset=p_patch%cells%owned)
    ! apply volume flux to surface elevation
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

      DO jc = i_startidx_c, i_endidx_c
        IF (p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0) THEN

          bt_lev = p_patch_3d%p_patch_1d(1)%dolic_c(jc, jb)      
          d_c    = p_patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, jb)

          !! We are using variable name zUnderIce but here it only refers
          !! to top level thickness

          temp_eta     = eta_c(jc,jb)              

          eta_c(jc,jb) = eta_c(jc,jb)               &
            &           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
            &           + p_oce_sfc%FrshFlux_TotalIce(jc, jb)*dtime  

          !! Only change the stretching parameter if it is above a certain threshold
          !! This avoids divide by 0 
          temp_stretch(jc, jb) = stretch_c(jc, jb)
          min_h                = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)
          
          !! Update only if height is atleast dz
          if ( d_c  .GT.  min_h ) &
            & temp_stretch(jc, jb) = ( eta_c(jc, jb) + d_c)/( d_c )
 
          !! update zunderice
          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * temp_stretch(jc, jb) 
    
          DO jk = 1, bt_lev 
            dz_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
              & * stretch_c(jc, jb) 
            dz_new = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, jb)&
              & * temp_stretch(jc, jb) 
 
            !! Scale temperature 
            p_os%p_prog(nold(1))%tracer(jc, jk, jb, 1) = & 
              & p_os%p_prog(nold(1))%tracer(jc, jk, jb, 1) * (dz_old / dz_new)
           
            !! Scale salinity
            p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) = & 
              & p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) * (dz_old / dz_new)
  
            !! Add salt flux from ice
            if (jk .EQ. 1) THEN
              p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) = p_os%p_prog(nold(1))%tracer(jc, jk, jb, 2) &
                &  +  p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime/dz_new  
            END IF
  
          END DO
 
        ENDIF  !  dolic>0
      END DO
    END DO
          
    IF (no_tracer > 0) THEN

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN

            ! substract the fraction of heatflux used for subsurface heating
            heatflux_surface_layer=p_oce_sfc%HeatFlux_Total(jc,jb)-p_os%p_diag%heatabs(jc,jb)
            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    heatflux_surface_layer*dtime/(clw*rho_ref*p_ice%zUnderIce(jc,jb))

          ENDIF
        ENDDO
      ENDDO

    END IF

    !! set correct cell thickness under ice
    p_oce_sfc%cellThicknessUnderIce   (:,:) = p_ice%zUnderIce(:,:)

    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)

  END SUBROUTINE apply_surface_fluxes_zstar_v14


  !>
  !! Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
  !!
  !! Update fluxes to apply as boundary condition  
  !
  SUBROUTINE apply_surface_fluxes_zstar_v15(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    !
    REAL(wp), INTENT(INOUT) :: eta_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(IN   ) :: stretch_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
 
    ! local variables
    INTEGER               :: jc, jb
    INTEGER               :: i_startidx_c, i_endidx_c
    REAL(wp)              :: sss_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: sst_inter(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)              :: zunderice_old
    REAL(wp)              :: zUnderIceIni(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    REAL(wp) :: h_old_test, h_new_test
    
    REAL(wp) :: temp_eta, min_h
    REAL(wp) :: temp_stretch(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    
    INTEGER  :: bt_lev, jk, adj_lev
    REAL(wp) :: d_c, dz_old, dz_new, extra_salt_in_col
    REAL(wp) :: dz_ratio, adj_lev_ht
    
    REAL(wp) :: old_sss, new_s1, new_s2, new_sss, thresh_sss
    INTEGER  :: flag 

    REAL(wp)  :: heatflux_surface_layer ! heatflux into the surface layer
    
    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    
    CHARACTER(LEN=max_char_length), PARAMETER :: str_module = 'ocean_sfc_refactor_zstar'
    
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------
    sst_inter(:,:)    = p_oce_sfc%sst(:,:)
    sss_inter(:,:)    = p_oce_sfc%sss(:,:)
    zUnderIceIni(:,:) = p_oce_sfc%cellThicknessUnderIce (:,:)

    CALL dbg_print('UpdSfcSTART: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcSTART: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)

    !!  Provide total ocean forcing:
    !    - total heat fluxes are aggregated for ice/ocean in ice thermodynamics
    !    - total internal salt flux p_oce_sfc%FrshFlux_TotalIce is calculated in sea ice model
    !    - total freshwater volume forcing
    p_oce_sfc%FrshFlux_VolumeTotal(:,:) = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_VolumeIce (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:) &
      &                                 + p_oce_sfc%FrshFlux_Relax     (:,:)
    ! provide total salinity forcing flux for diagnostics only
    p_oce_sfc%FrshFlux_TotalSalt(:,:)   = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalIce  (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:)

    
    CALL dbg_print('UpdSfc: eta-old', eta_c,    str_module, 1, in_subset=p_patch%cells%owned)
    ! apply volume flux to surface elevation
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

      DO jc = i_startidx_c, i_endidx_c
        IF (p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0) THEN

          !! update zunderice
          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) * stretch_c(jc, jb) 
    
        ENDIF  !  dolic>0
      END DO
    END DO
          
    IF (no_tracer > 0) THEN

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN

            ! substract the fraction of heatflux used for subsurface heating
            heatflux_surface_layer=p_oce_sfc%HeatFlux_Total(jc,jb)-p_os%p_diag%heatabs(jc,jb)
            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    heatflux_surface_layer*dtime/(clw*rho_ref*p_ice%zUnderIce(jc,jb))

          ENDIF
        ENDDO
      ENDDO

    END IF

    !! set correct cell thickness under ice
    p_oce_sfc%cellThicknessUnderIce   (:,:) = p_ice%zUnderIce(:,:)

    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)

  END SUBROUTINE apply_surface_fluxes_zstar_v15




  !-------------------------------------------------------------------------
  !
  !>
  !!  Update surface fluxes for ocean forcing. Analytical, OMIP, coupled.
  !!
  !!  OMIP: atmos_fluxes over ice and open water are calculated with the help of
  !!        bulk formulas in calc_omip_budgets_ice and calc_omip_budgets_oce
  !!  Coupled: passes E/P fluxes and heat fluxes over open ocean. Heat fluxes over
  !!           the ice surface are not updated, because they are only used in ice_fast,
  !!           and ice_fast is called within the atmosphere.
  !!
  !! @par Revision History
  !! Initial release by Vladimir Lapin, MPI-M (2016-11)
  !
!<Optimize_Used>
  SUBROUTINE update_atmos_fluxes(p_patch_3D, p_as, atmos_fluxes, p_oce_sfc, p_os, p_ice, this_datetime)

    TYPE(t_patch_3D ),TARGET,  INTENT(IN)       :: p_patch_3D
    TYPE(t_atmos_for_ocean)                     :: p_as
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    TYPE (t_hydro_ocean_state),INTENT(IN)       :: p_os
    TYPE (t_sea_ice),          INTENT(IN)       :: p_ice
    TYPE(datetime), POINTER                     :: this_datetime

    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_surface_refactor:update_atmos_fluxes'
    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells

    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------

    SELECT CASE (iforc_oce)

    CASE (Analytical_Forcing)        !  11      !  Driving the ocean with analytical fluxes

      CALL update_atmos_fluxes_analytical(p_patch_3D, p_os, p_ice, atmos_fluxes,p_oce_sfc)

    CASE (OMIP_FluxFromFile)         !  12      !  Driving the ocean with OMIP fluxes

      !   a) read OMIP data into p_as
      CALL update_flux_fromFile(p_patch_3D, p_as, this_datetime)

      !   b) calculate heat fluxes from p_as
      CALL calc_omip_budgets_oce(p_patch_3d, p_as, p_os, p_ice, atmos_fluxes)

      IF (i_sea_ice >= 1) THEN ! sea ice is on

          CALL calc_omip_budgets_ice(                            &
            &                        p_patch_3d,                                     &  !  input parameter
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

      ELSE   !  no sea ice

      ! apply net surface heat flux in W/m2 for OMIP case, since these fluxes are calculated in calc_omip_budgets_oce
        WHERE (p_patch_3D%lsm_c(:,1,:) <= sea_boundary)
          p_oce_sfc%HeatFlux_ShortWave(:,:) = atmos_fluxes%SWnetw(:,:) ! net SW radiation flux over water
          p_oce_sfc%HeatFlux_LongWave (:,:) = atmos_fluxes%LWnetw(:,:) ! net LW radiation flux over water
          p_oce_sfc%HeatFlux_Sensible (:,:) = atmos_fluxes%sensw (:,:) ! Sensible heat flux over water
          p_oce_sfc%HeatFlux_Latent   (:,:) = atmos_fluxes%latw  (:,:) ! Latent heat flux over water
          ! sum of ocean heat fluxes for ocean boundary condition without ice, generally aggregated in ice thermodynamics
          p_oce_sfc%HeatFlux_Total(:,:) = atmos_fluxes%SWnetw(:,:) + atmos_fluxes%LWnetw(:,:) &
            &                              + atmos_fluxes%sensw(:,:)  + atmos_fluxes%latw(:,:)
        ELSEWHERE
          p_oce_sfc%HeatFlux_ShortWave(:,:) = 0.0_wp
          p_oce_sfc%HeatFlux_LongWave (:,:) = 0.0_wp
          p_oce_sfc%HeatFlux_Sensible (:,:) = 0.0_wp
          p_oce_sfc%HeatFlux_Latent   (:,:) = 0.0_wp
          p_oce_sfc%HeatFlux_Total    (:,:) = 0.0_wp
        ENDWHERE

      ENDIF

      !   c) wind stress is assigned in calc_omip_budgets_oce
      !      over ice: stress_x, stress_y; and over open water: stress_xw, stress_yw

      !   d) freshwater fluxes from p_as

          ! provide evaporation from latent heat flux for OMIP case
          ! under sea ice evaporation is neglected, atmos_fluxes%latw is flux in the absence of sea ice
          p_oce_sfc%FrshFlux_Evaporation(:,:) = atmos_fluxes%latw(:,:) / (alv*rho_ref)

          !  copy variables into atmos_fluxes
          p_oce_sfc%FrshFlux_Runoff(:,:)      = p_as%FrshFlux_Runoff(:,:)

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
          p_oce_sfc%FrshFlux_TotalOcean(:,:) = p_patch_3d%wet_c(:,1,:)*( 1.0_wp-p_ice%concSum(:,:) ) * &
            &  (p_as%FrshFlux_Precipitation(:,:) + p_oce_sfc%FrshFlux_Evaporation(:,:))

    CASE (Coupled_FluxFromAtmo)

      !  Driving the ocean in a coupled mode:
      !  nothing to be done, atmospheric fluxes are provided at the end of time stepping
      !  atmospheric fluxes drive the ocean; fluxes are calculated by atmospheric model
      !  use atmospheric fluxes directly, i.e. no bulk formula as for OMIP is applied

       ! HAMOCC uses p_as to get SW radiation and wind, so we need to copy
       ! the SW radiation onto it in the coupled case
       if(lhamocc) p_as%fswr(:,:) = p_oce_sfc%HeatFlux_ShortWave(:,:)

! heatflux_total(:,:) is provided by coupling interface

      ! these 4 fluxes over open ocean are used in sea ice thermodynamics
      atmos_fluxes%SWnetw (:,:)   = p_oce_sfc%HeatFlux_ShortWave(:,:)
      atmos_fluxes%LWnetw (:,:)   = p_oce_sfc%HeatFlux_LongWave (:,:)
      atmos_fluxes%sensw  (:,:)   = p_oce_sfc%HeatFlux_Sensible (:,:)
      atmos_fluxes%latw   (:,:)   = p_oce_sfc%HeatFlux_Latent   (:,:)

      WHERE ( p_ice%concSum(:,:) > 0._wp) !  corresponding to (1-concSum)*Precip in TotalOcean
   !  WHERE ( ALL( p_ice%hi   (:,:,:) > 0._wp, 2 ) )  !  corresponding to hi>0 in ice_growth_zero
   !  WHERE ( ALL( p_ice%Tsurf(:,:,:) < 0._wp, 2 ) )  !  Tsurf is -1.8 over open water, incorrect specification
        ! SnowFall and liquid rain over ice-covered part of ocean are taken from the atmosphere model
        atmos_fluxes%rpreci(:,:) = p_oce_sfc%FrshFlux_SnowFall(:,:)
        atmos_fluxes%rprecw(:,:) = p_oce_sfc%FrshFlux_Precipitation(:,:) - p_oce_sfc%FrshFlux_SnowFall(:,:)
      ELSEWHERE
        ! not considered in ice_growth_zero
        atmos_fluxes%rpreci(:,:) = 0._wp
        atmos_fluxes%rprecw(:,:) = p_oce_sfc%FrshFlux_Precipitation(:,:)
      ENDWHERE

      ! copy flux for use in TotalOcean, since analytical/omip use p_as:
      !p_as%FrshFlux_Precipitation      = p_oce_sfc%FrshFlux_Precipitation

      ! total water flux over ice-free ocean water: P*(1-C)+E
      !  - whole evaporation over grid-box enters open ocean, this includes evaporation over sea ice covered part
      !  - snowfall is included as (melted) water equivalent
      !  - runoff is added to VolumeTotal below
      p_oce_sfc%FrshFlux_TotalOcean(:,:) = p_patch_3d%wet_c(:,1,:)* &
        &  (( 1.0_wp-p_ice%concSum(:,:) ) * p_oce_sfc%FrshFlux_Precipitation(:,:) + p_oce_sfc%FrshFlux_Evaporation(:,:))

    CASE DEFAULT

      CALL message(TRIM(routine), 'STOP: Ocean Forcing option not implemented' )
      CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION DOES NOT EXIST - TERMINATE')

    END SELECT


    IF (zero_freshwater_flux) THEN
      ! since latw<>0. we must set evap and TotalOcean again to zero:
      p_oce_sfc%FrshFlux_Evaporation  (:,:) = 0.0_wp
      p_oce_sfc%FrshFlux_TotalOcean   (:,:) = 0.0_wp
      p_oce_sfc%FrshFlux_Precipitation(:,:) = 0.0_wp
      p_oce_sfc%FrshFlux_SnowFall     (:,:) = 0.0_wp
      p_oce_sfc%FrshFlux_Evaporation  (:,:) = 0.0_wp
      p_oce_sfc%FrshFlux_Runoff       (:,:) = 0.0_wp
      p_oce_sfc%FrshFlux_TotalOcean   (:,:) = 0.0_wp
      atmos_fluxes%rpreci(:,:)              = 0.0_wp
      atmos_fluxes%rprecw(:,:)              = 0.0_wp
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    !idt_src=5  ! output print level (1-5, fix)
    !  these fluxes are always zero - fluxes over ice-covered area are Qbot, Qtop only
    !CALL dbg_print('aftAtmFluxUpd:atmflx%LWnetIce', atmos_fluxes%LWnet   ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !CALL dbg_print('aftAtmFluxUpd:atmflx%SensIce',  atmos_fluxes%sens    ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !CALL dbg_print('aftAtmFluxUpd:atmflx%LatentIce',atmos_fluxes%lat     ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !CALL dbg_print('aftAtmFluxUpd:atmflx%dsensdT'  ,atmos_fluxes%dsensdT ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !CALL dbg_print('aftAtmFluxUpd:atmflx%dlatdT'   ,atmos_fluxes%dlatdT  ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !CALL dbg_print('aftAtmFluxUpd:atmflx%dLWdT'    ,atmos_fluxes%dLWdt   ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !CALL dbg_print('aftAtmFluxUpd:stress_x'        ,atmos_fluxes%stress_x,str_module,idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------
    CALL dbg_print('aftAtmFluxUpd: Precipitation', p_oce_sfc%FrshFlux_Precipitation,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFluxUpd: Evaporation'  , p_oce_sfc%FrshFlux_Evaporation  ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFluxUpd: SnowFall'     , p_oce_sfc%FrshFlux_SnowFall     ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFluxUpd: Runoff'       , p_oce_sfc%FrshFlux_Runoff       ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFluxUpd: TotalOcean'   , p_oce_sfc%FrshFlux_TotalOcean   ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFluxUpd: rprecw'       , atmos_fluxes%rprecw                ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFluxUpd: rpreci'       , atmos_fluxes%rpreci                ,str_module, 3, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE update_atmos_fluxes

!**********************************************************************
!------------------------------ Analytical ----------------------------
!**********************************************************************

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
  SUBROUTINE update_atmos_fluxes_analytical(p_patch_3D, p_os, p_ice, atmos_fluxes, p_oce_sfc)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)    :: p_patch_3D
    TYPE(t_hydro_ocean_state), INTENT(IN)   :: p_os
    TYPE(t_sea_ice), INTENT(IN)             :: p_ice
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

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_bulk_forcing:update_atmos_fluxes_analytical'
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
      p_oce_sfc%FrshFlux_Precipitation(:,:) = atmos_precip_const

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

    !-----------------------------
    !   varios ad-hoc fixes
    !-----------------------------
    !  needed for old heat flux BC
    p_oce_sfc%HeatFlux_Total=p_oce_sfc%TopBC_Temp_vdiff

    ! provide dLWdt for ice_fast as for OMIP
    atmos_fluxes%dLWdT (:,:,:)  = -4._wp*zemiss_def*stbo*(p_ice%tsurf(:,:,:)+tmelt)**3

    ! provide evaporation from latent heat flux
    p_oce_sfc%FrshFlux_Evaporation(:,:) = atmos_fluxes%latw(:,:) / (alv*rho_ref)

  END SUBROUTINE update_atmos_fluxes_analytical

END MODULE mo_ocean_surface_refactor
