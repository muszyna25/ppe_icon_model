!>
!! Provide an implementation of the sea-ice model.
!!
!! Provide an implementation of the parameters of the surface module (sea ice)
!! used between the atmopshere and the hydrostatic ocean model.
!!
!! @author Peter Korn, MPI
!! @author Dirk Notz, MPI
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
!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_sea_ice_refactor
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2007
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime, ltimer
  USE mo_dynamics_config,     ONLY: nold, nnew
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D, t_patch_vert
  USE mo_exception,           ONLY: finish, message
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref, alf, mu, tf, clw
  USE mo_ocean_nml,           ONLY: no_tracer
  USE mo_sea_ice_nml,         ONLY: i_ice_therm, i_ice_dyn, i_Qio_type, use_constant_tfreez
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_sea_ice,             ONLY: upper_ocean_ts, ice_conc_change, ice_clean_up_thd, ice_zero, energy_content_in_surface
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_util_dbg_prnt,       ONLY: dbg_print
! USE mo_dbg_nml,             ONLY: idbg_mxmn, idbg_val
! USE mo_ice_fem_utils,       ONLY: fem_ice_wrap, ice_advection, ice_ocean_stress

!  USE mo_grid_config,         ONLY: n_dom   ! restrict sea-ice model to the global domain for the time being
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_ice_slow2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ice_slow_slo

  CHARACTER(len=12)           :: str_module    = 'SeaIce_Refac'  ! Output of module for 1 line debug

CONTAINS

  !-------------------------------------------------------------------------------
  !
  !
  !>
  !! !  ice_slow_slo: Ice routines for ocean time step. Calculates average of atmospheric
  ! !           time steps, ice velocity, ice growth rates and updates ice structure accordingly
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !! Rewrite by Stephan Lorenz, MPI-M (2015-01).
  !!
  SUBROUTINE ice_slow_slo(p_patch_3D, p_os, ice, atmos_fluxes, p_op_coeff)
    TYPE(t_patch_3D), TARGET, INTENT(IN)    :: p_patch_3D
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os
    TYPE(t_sea_ice),          INTENT(INOUT) :: ice
    TYPE(t_atmos_fluxes),     INTENT(INOUT) :: atmos_fluxes
    TYPE(t_operator_coeff),   INTENT(IN)    :: p_op_coeff

    TYPE(t_patch),      POINTER :: p_patch
    TYPE(t_patch_vert), POINTER :: p_patch_vert
    TYPE(t_subset_range), POINTER :: all_cells

   !!Local variables
    REAL(wp), DIMENSION (nproma,ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) :: Tfw ! Ocean freezing temperature [C]
    REAL(wp), DIMENSION (nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) :: energyCheck, sst, zuipsnowf
    REAL(wp), POINTER             :: flat(:,:)
    INTEGER                       :: k   !, jb, jc, i_startidx_c, i_endidx_c

    !-------------------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_ice_slow2)

    p_patch      => p_patch_3D%p_patch_2D(1)
    p_patch_vert => p_patch_3D%p_patch_1D(1)
    ! subset range pointer
    all_cells => p_patch%cells%all 
    flat      => p_patch_vert%prism_thick_flat_sfc_c(:,1,:)

    sst(:,:)        = 0.0_wp
    zuipsnowf(:,:)  = 0.0_wp

    CALL ice_zero       (ice)

    ! Save ice thickness at previous time step for calculation of heat and salt
    ! flux into ocean in subroutine upper_ocean_TS
    ice % hiold (:,:,:) = ice%hi(:,:,:)
    ice % hsold (:,:,:) = ice%hs(:,:,:)

    CALL dbg_print('IceSlow: hi before growth' ,ice%hi ,str_module,4, in_subset=p_patch%cells%owned)

    ! needs central variable ice%Tfw and calculation once in ice_slow
    ! freezing temperature of uppermost sea water
    IF ( no_tracer < 2 .OR. use_constant_tfreez ) THEN
      Tfw(:,:,:) = Tf
    ELSE
      DO k=1,ice%kice
        Tfw(:,k,:) = -mu * p_os%p_prog(nold(1))%tracer(:,1,:,2)
      ENDDO
    ENDIF

!---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('GrowZero: heatOceI bef.grow' , ice%heatOceI   , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: Qtop bef. growth'  , ice%Qtop       , str_module, 5, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: Qbot bef. growth'  , ice%Qbot       , str_module, 5, in_subset=p_patch%cells%owned)
!---------------------------------------------------------------------
    
    IF (i_ice_therm /= 3 ) THEN
      ! Heat flux from ocean into ice
      CALL oce_ice_heatflx_slo (p_patch, p_os,ice,Tfw,ice%zHeatOceI)
!!$    ELSE IF ( i_ice_therm == 3) THEN
      ! for i_ice_therm == 3, no ocean-ice heatflx is included!
    END IF

    IF ( i_ice_therm == 1 ) &
      & CALL ice_growth_zero_slo (p_patch, p_os, ice, atmos_fluxes%rpreci)

    sst(:,:) = p_os%p_prog(nold(1))%tracer(:,1,:,1)
    !---  energy  -----------------------------------------------------------
    ! updates should be moved to respective routine ice_growth - done later in upperOceTS
    !  - 2015-01-19: ice growth/melt must not change zunderice yet
 !  ice%draftave (:,:) = (rhos * ice%hs(:,1,:) + rhoi * ice%hi(:,1,:)) * ice%conc(:,1,:) / rho_ref
 !  ice%zUnderIce(:,:) = flat(:,:) + p_os%p_prog(nold(1))%h(:,:) - ice%draftave(:,:)
    energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
      &             ice, sst(:,:), 0, info='AFT GROWTH')

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('IceSlow: hi after growth'   ,ice%hi       ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: Conc. after growth',ice%conc     ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: ice%u bef. dyn'    ,ice%u_prog   ,str_module, 4, in_subset=p_patch%verts%owned)
    CALL dbg_print('IceSlow: ice%v bef. dyn'    ,ice%v_prog   ,str_module, 4, in_subset=p_patch%verts%owned)
    CALL dbg_print('IceSlow: zUnderIce aft.gr',  ice%zUnderIce,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: heatOceI aftgrowth',ice%heatOceI ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: SST    aft. growth',sst          ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: energy aft. growth',energyCheck  ,str_module, 2, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    ! #slo# 2014-12: now use melt/growth energy heatOceI from below to change SST
    !  heatOceI=heatOceI*conc - was done later in upperOceTS - generates inconsistency
    !IF (i_Qio_type < 3) ice%heatOceI(:,1,:) = ice%heatOceI(:,1,:)*ice%conc(:,1,:)
    ! wrong - ice covered part only is used for heat exchange with ocean
    !IF (i_Qio_type == 3) ice%heatOceI(:,1,:) = ice%heatOceI(:,1,:)*ice%conc(:,1,:) 
    ice%heatOceI(:,1,:) = ice%heatOceI(:,1,:)*ice%conc(:,1,:)
    CALL dbg_print('IceSlow: heatOceI aftConCor',ice%heatOceI ,str_module, 4, in_subset=p_patch%cells%owned)

    ! #slo# 2015-01: now move sst-change back to surface module (mo_ocean_bulk)
  ! DO jb = all_cells%start_block, all_cells%end_block
  !   CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
  !   DO jc = i_startidx_c, i_endidx_c
  !     IF ( v_base%lsm_c(jc,1,jb) <= sea_boundary ) THEN
  !       sst(jc,jb) = p_os%p_prog(nold(1))%tracer(jc,1,jb,1)
  !       sst(jc,jb) = sst(jc,jb) + ice%heatOceI(jc,1,jb)*dtime/(clw*rho_ref*ice%zUnderIce(jc,jb))
  !       ! set new sst; heatOceI to zero
  !       p_os%p_prog(nold(1))%tracer(jc,1,jb,1) = sst(jc,jb)
  !       ice%heatOceI(jc,1,jb)                  = 0.0_wp
  !     ENDIF
  !   ENDDO
  ! ENDDO

    energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
      &             ice, sst(:,:), 0, info='AFT UPGROW')
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('IceSlow: heatOceI aftupGrow',ice%heatOceI ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: SST    aft. upGrow',sst          ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: energy aft. upGrow',energyCheck  ,str_module, 4, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    ! Ocean Heat Flux
    CALL upper_ocean_TS(p_patch,p_patch_vert,ice,p_os,atmos_fluxes)

    sst(:,:) = p_os%p_prog(nold(1))%tracer(:,1,:,1)  !  add corresponding flux or calculate zUnderIce
 !  ice%zUnderIce(:,:) = flat(:,:) + p_os%p_prog(nold(1))%h(:,:) &
 !    &                - (rhos * ice%hs(:,1,:) + rhoi * ice%hi(:,1,:)) * ice%conc(:,1,:) / rho_ref
    energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
      &             ice, sst(:,:), 0, info='AFT UOCETS')
    CALL dbg_print('IceSlow: zUnderIce aft.UPTS',ice%zUnderIce,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: energy aft.OceanTS',energyCheck  ,str_module, 2, in_subset=p_patch%cells%owned)

    ! #slo# 2014-12: now use energy atmos_fluxes%HeatFlux_Total from upperOceTS to change SST
    ! #slo# 2015-01: now move sst-change back to surface module (mo_ocean_bulk) using HeatFlux_Total
  ! DO jb = all_cells%start_block, all_cells%end_block
  !   CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
  !   DO jc = i_startidx_c, i_endidx_c
  !     IF ( v_base%lsm_c(jc,1,jb) <= sea_boundary ) THEN
  !       sst(jc,jb) = p_os%p_prog(nold(1))%tracer(jc,1,jb,1)
  !       sst(jc,jb) = sst(jc,jb) + atmos_fluxes%HeatFlux_Total(jc,jb)*dtime/(clw*rho_ref*ice%zUnderIce(jc,jb))
  !       ! set new sst; HeatFlux_Total to zero
  !       p_os%p_prog(nold(1))%tracer(jc,1,jb,1) = sst(jc,jb)
  !       atmos_fluxes%HeatFlux_Total(jc,jb)     = 0.0_wp
  !     ENDIF
  !   ENDDO
  ! ENDDO
    energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
      &             ice, sst(:,:), 0, info='AFT UPUPTS')
    CALL dbg_print('IceSlow: HeatTotal aft.UPTS',atmos_fluxes%HeatFlux_Total,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: SST    aft. upUPTS',sst          ,str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: energy aft. upUPTS',energyCheck  ,str_module, 2, in_subset=p_patch%cells%owned)

    ! Ice Concentration Change
    IF ( i_ice_therm >= 1 ) THEN
        CALL ice_conc_change(p_patch,ice,p_os)
    ENDIF

    ! #slo# 2015-01 - Test: update draft, zunderice now includes totalsnowfall as in ocets, inconsistent with h
    !               - update of draft/zunderice should be done whenever draft is changed
 !  ice%draft       (:,:,:) = (rhos * ice%hs(:,:,:) + rhoi * ice%hi(:,:,:)) / rho_ref
 !  ice%draftave    (:,:)   = sum(ice%draft(:,:,:) * ice%conc(:,:,:),2)
 !  ice%zUnderIce   (:,:)   = p_patch_vert%prism_thick_flat_sfc_c(:,1,:) &
 !    &                            + p_os%p_prog(nold(1))%h(:,:) &
 !    &                            - ice%draftave(:,:)& 
 !    &                            + ice%totalsnowfall(:,:)
    ! - Test: update zunderice only, includes totalsnowfall
    zuipsnowf    (:,:) = flat(:,:) + p_os%p_prog(nold(1))%h(:,:) &
      &                - (rhos * ice%hs(:,1,:) + rhoi * ice%hi(:,1,:)) * ice%conc(:,1,:) / rho_ref &
      &                            + ice%totalsnowfall(:,:)


    sst(:,:) = p_os%p_prog(nold(1))%tracer(:,1,:,1)  !  add corresponding flux or calculate zUnderIce
    ice%zUnderIce(:,:) = flat(:,:) + p_os%p_prog(nold(1))%h(:,:) &
      &                - (rhos * ice%hs(:,1,:) + rhoi * ice%hi(:,1,:)) * ice%conc(:,1,:) / rho_ref
    energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
      &             ice, sst(:,:), 0, info='AFT CONCCH-0')

    CALL dbg_print('IceSlow: energy aft. ConcCh',energyCheck  ,str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: zUnderIce a.ConcCh',ice%zUnderIce,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: zUI+snowf a.ConcCh',zuipsnowf,    str_module, 4, in_subset=p_patch%cells%owned)

    energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
      &             ice, sst(:,:), 1, info='AFT CONCCH-1')
    CALL dbg_print('IceSlow: energy typ1 ConcCh',energyCheck  ,str_module, 4, in_subset=p_patch%cells%owned)

    ! update draftave for comparison to energy using type 0 - difference by snowfall, see above: 
    ice%draftave (:,:) = (rhos * ice%hs(:,1,:) + rhoi * ice%hi(:,1,:)) * ice%conc(:,1,:) / rho_ref
    ice%zUnderIce(:,:) = flat(:,:) + p_os%p_prog(nold(1))%h(:,:) - ice%draftave(:,:)

    energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
      &             ice, sst(:,:), 0, info='AFT CONCCH-corr0')
    CALL dbg_print('IceSlow: energy typ0 draft',energyCheck  ,str_module, 4, in_subset=p_patch%cells%owned)

    ! old call to ice dynamics was here
    ! check after ice dynamics, method 1 since draftave is not yet updated
  ! energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
  !   &             ice, sst(:,:), 1, info='AFT ICEDYN')
  ! CALL dbg_print('IceSlow: energy aftIceAdv',energyCheck  ,str_module, 2, in_subset=p_patch%cells%owned)

    CALL dbg_print('IceSlow: hi    bef.cleanup',ice%hi,       str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: hs    bef.cleanup',ice%hs,       str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: Conc. bef.cleanup',ice%conc     ,str_module, 3, in_subset=p_patch%cells%owned)

    ! the original clean up routine has been split into two: ice_clean_up_dyn, ice_clean_up_thd
    ! ice_clean_up_thd: limit sea ice thickness to seaice_limit of surface layer depth after changes due to the thermodynamic growth/melt
    CALL ice_clean_up_thd( p_patch_3D, ice, atmos_fluxes, p_os )

    !  last check after cleanup with draftave update:
    energyCheck = energy_content_in_surface(p_patch, flat(:,:), p_os%p_prog(nold(1))%h(:,:), &
      &             ice, sst(:,:), 1, info='AFT ICEDYN')
    CALL dbg_print('IceSlow: energy aft CleanUp',energyCheck  ,str_module, 2, in_subset=p_patch%cells%owned)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('IceSlow: hi endOf slow'     ,ice%hi,                 str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: hs endOf slow'     ,ice%hs,                 str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: ConcSumEndOf slow', ice%concSum,            str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: Conc.  EndOf slow', ice%conc,               str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: p_ice%u'           ,ice%u_prog,             str_module, 4, in_subset=p_patch%verts%owned)
    CALL dbg_print('IceSlow: p_ice%v'           ,ice%v_prog,             str_module, 4, in_subset=p_patch%verts%owned)
    CALL dbg_print('IceSlow: p_os%prog(nold)%vn',p_os%p_prog(nold(1))%vn,str_module, 5, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: p_os%prog(nnew)%vn',p_os%p_prog(nnew(1))%vn,str_module, 5, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: p_os%diag%u'       ,p_os%p_diag%u,          str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: p_os%diag%v'       ,p_os%p_diag%v,          str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: atmFlx%windStr-u' ,atmos_fluxes%topBoundCond_windStress_u,str_module,4,in_subset=p_patch%cells%owned)
    CALL dbg_print('IceSlow: atmFlx%windStr-v' ,atmos_fluxes%topBoundCond_windStress_v,str_module,4,in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    IF (ltimer) CALL timer_stop(timer_ice_slow2)

  END SUBROUTINE ice_slow_slo

  !-------------------------------------------------------------------------------
  !  
  !>
  !! ! ice_growth_zero_slo - change ice and snow thickness (Semtner 1976, Appendix)
  !! This function changes:
  !! --- currently not fully --- ice % hs       new snow thickness for each ice category                [m]
  !! ice % hi       new ice  thickness for each ice category                [m]
  !! --- not currently --- ice % evapwi   amount of evaporated water from the mixed layer
  !!                in previously ice covered areas if all ice is gone      [kg/m^3]
  !! ice % heatOceI to contain the energy that is available to the mixed layer
  !!                in previously ice covered areas if all ice is gone      [J]
  !
  ! The counterpart to the  ice_growth subroutine in MPIOM
  !
  !!
  !! @par Revision History
  !! Initial release by Achim Randelhoff
  !! Update and rewrite by Stephan Lorenz, MPI-M (2015-01).
  !!
 
 SUBROUTINE ice_growth_zero_slo(p_patch, p_os, ice, rpreci)
   TYPE(t_patch),             INTENT(IN), TARGET    :: p_patch 
   TYPE(t_hydro_ocean_state), INTENT(IN)            :: p_os
   TYPE (t_sea_ice),          INTENT(INOUT)         :: ice
   REAL(wp),                  INTENT(IN)            :: rpreci(:,:) 
                                   ! water equiv. solid precipitation rate [m/s] DIMENSION (ie,je)

   !!Local variables
    REAL(wp), DIMENSION (nproma,ice%kice, p_patch%alloc_cell_blocks) ::         &
      & Tfw,         & ! Ocean freezing temperature [C]
      & Q_surplus   ! energy surplus during ice growth
    
    REAL(wp) ::      &
      & below_water, & ! Thickness of snow layer below water line           [m]
      & draft          ! depth of ice-ocean interface below sea level       [m]

    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: k, jb, jc, i_startidx_c, i_endidx_c
    
    all_cells            => p_patch%cells%all
    Q_surplus(:,:,:)     =  0.0_wp
    Tfw(:,:,:)           =  0.0_wp

    ! freezing temperature of uppermost sea water
    IF ( no_tracer < 2 .OR. use_constant_tfreez ) THEN
      Tfw(:,:,:) = Tf
    ELSE
      DO k=1,ice%kice
        Tfw(:,k,:) = -mu * p_os%p_prog(nold(1))%tracer(:,1,:,2)
      ENDDO
    ENDIF

!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, k, jc, draft, below_water) SCHEDULE(dynamic)
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO k=1,ice%kice
        DO jc = i_startidx_c,i_endidx_c
          IF (ice%hi(jc,k,jb) > 0._wp) THEN

            ! Add oceanic heat flux to energy available at the bottom of the ice.
            ice%Qbot(jc,k,jb) = ice%Qbot(jc,k,jb) + ice%zHeatOceI(jc,k,jb)

            ! Add snowfall to snow depth
            ! #slo# 2015-01: bugfix: rpreci is rate of snowfall over ice covered area
            ice%hs(jc,k,jb) = ice%hs(jc,k,jb) + rpreci(jc,jb)*dtime*rho_ref/rhos
            ! #slo# 2015-01: bugfix: rpreci is over whole grid-area
            !  this is incorrect, because hs is thickness over ice-covered area only and rpreci is a snowfall rate
            !ice%hs(jc,k,jb) = ice%hs(jc,k,jb) + rpreci(jc,jb)*ice%conc(jc,k,jb)*dtime*rho_ref/rhos
      
            ! for energy flux surplus
            IF ( ice%Qtop(jc,k,jb) > 0.0_wp ) THEN 
              IF  ( ice%hs(jc,k,jb) > 0.0_wp )  THEN ! melt snow where there's snow
                
                ice%hs (jc,k,jb) =  ice%hs(jc,k,jb) - ice%Qtop(jc,k,jb) * dtime / (alf*rhos) 
                ! put remaining heat, if any, into melting ice below
                IF (ice%hs(jc,k,jb) < 0.0_wp) THEN
                  ice%hi(jc,k,jb) = ice%hi(jc,k,jb) + ice%hs(jc,k,jb) * (rhos/rhoi) ! snow thickness loss in ice equivalents
                  ice%hs(jc,k,jb) = 0.0_wp
                ENDIF
                
              ELSE   ! where there's no snow
                ice%hi(jc,k,jb) = ice%hi(jc,k,jb) - ice%Qtop(jc,k,jb) * dtime / (alf*rhoi) 
              ENDIF
            ENDIF
            
            ! bottom melt/freeze
            ice%hi(jc,k,jb) = ice%hi(jc,k,jb) - ice%Qbot(jc,k,jb) * dtime / (alf*rhoi )
            
            ! heat to remove from water
            !  - heatOceI - positive into ocean - positive=downward i.e. same sign convention as HeatFlux_Total into ocean
            !  - zHeatOceI - positive into ice, i.e. positive=upward - melting energy coming from below, from the ocean
            ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) - ice%zHeatOceI(jc,k,jb)
            
            ! hi<0 if melting energy (Qbot+zHeatOceI) is larger than needed to melt all ice and snow, see above
            IF (ice%hi (jc,k,jb) <= 0.0_wp) THEN

              ! remove surplus energy of ice thickness from water
              !  - hi<0, if all ice and snow is already melted
              !  - calculate surplus of heatOceI>0 available for heating of ocean after complete melting
              ! #slo# 2014-11: 3. bugfix: sign error in hi for heatOceI
              !ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) + ice%hi(jc,k,jb)*alf*rhoi/dtime
              ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) - ice%hi(jc,k,jb)*alf*rhoi/dtime

              ! remove latent heat of snow from water
              ! #slo# 2014-11: if there is snow on top of melted ice, hs>0, then heatOceI is reduced by latent heat of snow
              ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) - ice%hs(jc,k,jb)*alf*rhos/dtime
              ! #slo# 2014-11: Attention: if heatOceI is not enough to melt whole snow,
              !                then energy budget is not closed! TODO: check and correct (later, little energy)
              ! IF ( ice%heatOceI(jc,k,jb) >0 ) THEN
              !  - snow is set to rest of ice, since no snow without water is possible
              !   ice%hi(jc,k,jb) = ice%heatOceI(jc,k,jb) * dtime / (alf * rhoi )
              !   ice%heatOceI(jc,k,jb) = 0.0_wp
              ! ELSE  ! melting energy was enough to melt all snow and ice
              
              ice%Tsurf(jc,k,jb) =  Tfw(jc,k,jb)
              !ice%conc (jc,k,jb) = 0.0_wp  !  do not change concentration here, but in ice_conc_change only
              ice%hi   (jc,k,jb) = 0.0_wp
              ice%hs   (jc,k,jb) = 0.0_wp

              ! ENDIF ! ( ice%heatOceI(jc,k,jb) >0 )

            ENDIF

            ! #slo# 2014-11: 2. bugfix: calculation moved down to below recalculaton of hi
            ! #slo# 2015-01: could we update ice%draft here or not?
            draft           = ( rhoi*ice%hi(jc,k,jb) + rhos*ice%hs(jc,k,jb) )/rho_ref
            below_water     = draft - ice%hi(jc,k,jb)  !  thickness to be converted to ice
            
            ! snow -> ice conversion for snow below waterlevel
            ! Currently not quite physical: Snow is pushed together to form new ice, hence snow thickness
            ! decreases more than ice thickness by rhoi/rhos ( analogue to old growth.f90 sea-ice model )
            ! Salt content of snow ice is equal to that of normal ice, salt is removed from the ocean
            ! Temperature of new upper ice is calculated as described in the paragraph below 
            ! Eq. 36
            IF ( below_water > 0.0_wp ) THEN
              ice%snow_to_ice(jc,k,jb) = below_water*rhoi/rhos     ! Thickness of snow that is converted into ice
              ice%hs         (jc,k,jb) = ice%hs(jc,k,jb) - ice%snow_to_ice(jc,k,jb)
              ice%hi         (jc,k,jb) = ice%hi(jc,k,jb) + below_water
            END IF

            IF (ice%hs (jc,k,jb) < 0.0_wp) THEN
               ice % hs(jc,k,jb) = 0.0_wp
               ice % hi(jc,k,jb) = ice%hi(jc,k,jb) + ice%hs(jc,k,jb) * (rhos/rhoi) ! snow thickness loss in ice equivalents
            ENDIF
            
            ! check energy conservation
            ! surplus energy = entering - leaving - latent heat
            !!! what's up with the energy that's put into the ocean?
            ! #slo# 2015-01: snowfall changes energy input - not yet considered
            Q_surplus(jc,k,jb) = &!0.0_wp
              &                   ice%Qbot(jc,k,jb) + ice%Qtop(jc,k,jb) &
              &                   + (ice%hi(jc,k,jb)-ice%hiold(jc,k,jb)) *alf*rhoi/dtime&
              &                   + (ice%hs(jc,k,jb)-ice%hsold(jc,k,jb)) *alf*rhos/dtime

          ELSE  !  hi<=0
            ! #slo# 2014-12: check - heatOceI is set in case of no ice - negative ice possible?
            ice%heatOceI(jc,k,jb) = ice%Qtop(jc,k,jb) + ice%Qbot(jc,k,jb)
            ice%Tsurf(jc,k,jb) =  Tfw(jc,k,jb)
            ice%conc (jc,k,jb) = 0.0_wp
            ice%hi   (jc,k,jb) = 0.0_wp
            ice%hs   (jc,k,jb) = 0.0_wp
          ENDIF

          ! #slo# 2014-12: update zUnderIce better here than in ice_slow?
       !  ice%zUnderIce(:,:) = flat(:,:) + p_os%p_prog(nold(1))%h(:,:) &
       !    &                - (rhos * ice%hs(:,1,:) + rhoi * ice%hi(:,1,:)) * ice%conc(:,1,:) / rho_ref

        END DO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('GrowZero: snow_to_ice', ice%snow_to_ice, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: hi'         , ice%hi         , str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: hs'         , ice%hs         , str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: zHeatOceI'  , ice%zHeatOceI  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: heatOceI '  , ice%heatOceI   , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: Q_surplus'  , Q_surplus      , str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: Qtop'       , ice%Qtop       , str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: Qbot'       , ice%Qbot       , str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: Tsurf'      , ice%Tsurf      , str_module, 4, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------
 
  END SUBROUTINE ice_growth_zero_slo

  !-------------------------------------
  !
  ! oce_ice_heatflx_slo
  !
  ! Calculates the heat flux from the uppermost water layer into the ice.
  !
  ! Currently (as in growth.f90): all energy available in upper ocean grid cell 
  ! is supplied to the ice and the upper ocean temperature is held at the 
  ! freezing point. This is not very physical.
  !
  ! Positive flux upwards.
 
  
  SUBROUTINE oce_ice_heatflx_slo (p_patch, p_os, ice, Tfw, zHeatOceI)
    TYPE(t_patch)            , INTENT(IN), TARGET    :: p_patch
    TYPE(t_hydro_ocean_state), INTENT(IN)  :: p_os
    TYPE(t_sea_ice)          , INTENT(IN)  :: ice
    REAL(wp)                 , INTENT(IN)  :: Tfw(:,:,:)      ! freezing temperature
    REAL(wp)                 , INTENT(OUT) :: zHeatOceI(:,:,:)

    ! Local
    INTEGER :: k, jb, jc, i_startidx_c, i_endidx_c
    TYPE(t_subset_range), POINTER :: all_cells
 !  REAL(wp) :: u_star
    REAL(wp), POINTER  :: sst(:,:)

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice_shared_sr:oce_ice_heatflx'
    
    all_cells => p_patch%cells%all 
    zHeatOceI(:,:,:) = 0.0_wp
    sst => p_os%p_prog(nold(1))%tracer(:,1,:,1)

    ! calculate heat flux from ocean to ice  (zHeatOceI) 
    SELECT CASE ( i_Qio_type )
    CASE (1)
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, k, jc) SCHEDULE(dynamic)
      DO jb = 1,p_patch%nblks_c
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c,i_endidx_c
          DO k=1,ice%kice
            IF (ice%hi(jc,k,jb) > 0._wp) THEN
              ! energy of warm water below ice covered part of grid area only is used for melting
              zHeatOceI(jc,k,jb) = ( sst(jc,jb) - Tfw(jc,k,jb) ) * ice%zUnderIce(jc,jb) * clw*rho_ref/dtime
            ENDIF
          ENDDO
        ENDDO
      END DO
!ICON_OMP_END_PARALLEL_DO

   CASE(2)

! !ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, k, jc) SCHEDULE(dynamic)
  ! DO jb = 1,p_patch%nblks_c
  !   CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
  !   DO jc = i_startidx_c,i_endidx_c
  !     DO k=1,ice%kice
  !       IF (ice%hi(jc,k,jb) > 0._wp) THEN
  !         ! melting energy depends on velocity difference between water and ice, bulk formulation
  !           u_star = SQRT(Cd_io*( (p_os%p_diag%u(jc,1,jb)-ice%u(jc,jb))**2 + &
  !             &         (p_os%p_diag%v(jc,1,jb)-ice%v(jc,jb))**2 ))
  !           zHeatOceI(jc,k,jb) = ( sst(jc,jb) - Tfw(jc,k,jb) ) *rho_ref*clw*Ch_io*u_star
  !       ENDIF
  !     ENDDO
  !   ENDDO
  ! END DO
! !ICON_OMP_END_PARALLEL_DO

    CASE (3)
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, k, jc) SCHEDULE(dynamic)
      DO jb = 1,p_patch%nblks_c
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c,i_endidx_c
          DO k=1,ice%kice
            IF (ice%hi(jc,k,jb) > 0._wp) THEN
              ! ALL energy of warm water over the whole grid area is used for melting ice - divide by concentration
              ! SLO/EO 2014-04-11 - this is wrong, must be accompanied by correction elsewhere,
              !                     since open part of water is still losing heat
              zHeatOceI(jc,k,jb) = ( sst(jc,jb) - Tfw(jc,k,jb) )*ice%zUnderIce(jc,jb)*clw*rho_ref/(dtime*ice%conc(jc,k,jb))
            ENDIF
          ENDDO
        ENDDO
      END DO
!ICON_OMP_END_PARALLEL_DO

    
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Invalid i_Qio_type')
    END SELECT
    
    CALL dbg_print('o-i-heat: SST'       ,sst          ,str_module,4, in_subset=p_patch%cells%owned)
    CALL dbg_print('o-i-heat: Tfw'       ,Tfw          ,str_module,4, in_subset=p_patch%cells%owned)
    CALL dbg_print('o-i-heat: zUnderIce' ,ice%zUnderIce,str_module,4, in_subset=p_patch%cells%owned)
    CALL dbg_print('o-i-heat: zHeatOceI' ,zHeatOceI    ,str_module,3, in_subset=p_patch%cells%owned)
    
  END SUBROUTINE oce_ice_heatflx_slo

END MODULE mo_sea_ice_refactor
