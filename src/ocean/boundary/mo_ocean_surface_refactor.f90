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
    &  limit_elevation, nbgctra, lhamocc

  USE mo_ocean_nml,           ONLY: atmos_flux_analytical_type, relax_analytical_type, &
    &  n_zlev, para_surfRelax_Salt, para_surfRelax_Temp, atmos_precip_const, &  ! atmos_evap_constant
    &  atmos_SWnet_const, atmos_LWnet_const, atmos_lat_const, atmos_sens_const, &
    &  atmos_SWnetw_const, atmos_LWnetw_const, atmos_latw_const, atmos_sensw_const, &
    &  relax_width, forcing_HeatFlux_amplitude, forcing_HeatFlux_base,              &
    &  basin_center_lat, basin_center_lon, basin_width_deg, basin_height_deg

  USE mo_math_constants,      ONLY: pi, deg2rad, rad2deg
  USE mo_physical_constants,  ONLY: rho_ref, alv, tmelt, tf, clw, stbo, zemiss_def
  USE mo_impl_constants,      ONLY: max_char_length, sea_boundary, MIN_DOLIC

  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_ocean_surface_types, ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff

  USE mtime,                  ONLY: datetime, getNoOfSecondsElapsedInDayDateTime
  USE mo_ice_interface,       ONLY: ice_fast_interface, ice_slow_interface
  USE mo_ocean_bulk_forcing,  ONLY: update_surface_relaxation, apply_surface_relaxation, &
                                &   update_flux_fromFile, calc_omip_budgets_ice, calc_omip_budgets_oce, &
                                &   update_ocean_surface_stress, balance_elevation

  USE mo_ocean_diagnostics, ONLY : diag_heat_tendency
  USE mo_name_list_output_init, ONLY: isRegistered

  IMPLICIT NONE
  
  PRIVATE

  ! public interface
  PUBLIC  :: update_ocean_surface_refactor
  ! private routine
  PRIVATE :: apply_surface_fluxes
  PRIVATE :: apply_surface_fluxes_slo
  PRIVATE :: update_atmos_fluxes
  PRIVATE :: update_atmos_fluxes_analytical

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

    IF (no_tracer>=1) p_oce_sfc%sst => p_os%p_prog(nold(1))%tracer(:,1,:,1)
    IF (no_tracer>=2) p_oce_sfc%sss => p_os%p_prog(nold(1))%tracer(:,1,:,2)

    IF (iforc_oce == No_Forcing) RETURN  !  forcing for ocean not defined

    ! save values of ice, snow and temperature for the tendendy and hfbasin diagnostic
    IF ( isRegistered('delta_ice') .OR. isRegistered('delta_snow') .OR. &
           isRegistered('delta_thetao') .OR. &
           isRegistered('global_hfbasin') .OR. isRegistered('atlant_hfbasin') .OR. &
           isRegistered('pacind_hfbasin') ) THEN  

      CALL diag_heat_tendency(p_patch_3d, 1, p_ice,            &
         p_os%p_prog(nold(1))%tracer(:,:,:,1),               &
         p_os%p_diag%delta_ice,                              &
         p_os%p_diag%delta_snow, p_os%p_diag%delta_thetao)

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

        !  (3b) Slow sea ice dynamics and thermodynamics
        CALL ice_slow_interface(p_patch_3D, p_ice, p_oce_sfc, atmos_fluxes, p_os, p_as, p_op_coeff)

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
!    CALL apply_surface_fluxes(p_patch_3D, p_os, p_ice, p_oce_sfc)
    CALL apply_surface_fluxes_slo(p_patch_3D, p_os, p_ice, p_oce_sfc)

    !---------------------------------------------------------------------
    ! (6) Apply volume flux correction
    !---------------------------------------------------------------------
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

    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------

    ! Provide total freshwater volume forcing:
    p_oce_sfc%FrshFlux_VolumeTotal(:,:) = p_oce_sfc%FrshFlux_Runoff    (:,:) &
      &                                 + p_oce_sfc%FrshFlux_VolumeIce (:,:) &
      &                                 + p_oce_sfc%FrshFlux_TotalOcean(:,:) &
      &                                 + p_oce_sfc%FrshFlux_Relax     (:,:)

    ! Note: old freeboard is stored in p_oce_sfc%cellThicknessUnderIce (equiv. to zUnderIceIni in apply_surface_fluxes_slo)

    IF (no_tracer > 0) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, jc, i_bgc_tra) SCHEDULE(dynamic)
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN

            ! (5a) Net heat flux changes sst using old freeboard (Thermodynamic Eq. 1)
            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    p_oce_sfc%HeatFlux_Total(jc,jb)*dtime/(clw*rho_ref*p_oce_sfc%cellThicknessUnderIce(jc,jb))

            ! (5b) Net volume flux (plus snow fall on ice) changes ssh (Thermodynamic Eq. 5)
            p_os%p_prog(nold(1))%h(jc,jb) = p_os%p_prog(nold(1))%h(jc,jb)               &
              &                           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
              &                           + p_ice%totalsnowfall(jc,jb)

            ! (5c) New zUnderIce calculated with NEW h and draftave
            p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) &
              &                    + p_os%p_prog(nold(1))%h(jc,jb) - p_ice%draftave(jc,jb)

            ! (5d) New salinity is calculated from conservation formula:
            !      SSS_new * zUnderIce_new = ( SSS_old * zUnderIceOld + SaltFluxFromIce * dtime )
            p_oce_sfc%sss(jc,jb)   = ( p_oce_sfc%sss(jc,jb) * p_oce_sfc%cellThicknessUnderIce(jc,jb) + &
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

    REAL(wp) :: h_old_test

    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER:: i_bgc_tra
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
            p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
              &                    p_oce_sfc%HeatFlux_Total(jc,jb)*dtime/(clw*rho_ref*zUnderIceIni(jc,jb))
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

          h_old_test=  (p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb)+p_os%p_prog(nold(1))%h(jc,jb))


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
           &        * h_old_test/(p_os%p_prog(nold(1))%h(jc,jb) + p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb))
          ENDDO
          endif

          

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
      CALL calc_omip_budgets_oce(p_patch, p_as, p_os, atmos_fluxes)

      IF (i_sea_ice >= 1) THEN ! sea ice is on

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
