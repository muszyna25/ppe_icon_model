!>
!! Provide interfaces to call sea ice ice model from the ocean and atmosphere.
!!
!! @author Dirk Notz, MPI
!! @author Vladimir Lapin, MPI
!!
!! @par Revision History
!! Original version by Peter Korn, MPI-M (2009)
!! Modified by Vladimir Lapin, MPI-M (2017)
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
MODULE mo_ice_interface
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
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_exception,           ONLY: finish
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_util_dbg_prnt,       ONLY: dbg_print

  USE mo_timer,               ONLY: timer_start, timer_stop, timer_ice_fast,   &
    &                               timers_level, timer_extra40, timer_ice_advection
  USE mtime,                  ONLY: datetime, getDayOfYearFromDateTime

  USE mo_physical_constants,  ONLY: rhoi, ki, Tf, ci
  USE mo_sea_ice_nml,         ONLY: i_ice_therm, i_ice_dyn, i_ice_advec, hci_layer
  USE mo_ocean_nml,           ONLY: atmos_flux_analytical_type, atmos_SWnet_const, atmos_sens_const
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_surface_types, ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff

  USE mo_ice_winton,          ONLY: set_ice_temp_winton
  USE mo_ice_zerolayer,       ONLY: set_ice_temp_zerolayer, set_ice_temp_zerolayer_analytical
  USE mo_ice_slow_thermo,     ONLY: ice_slow_thermo
  USE mo_ice_parameterizations, ONLY: ice_cut_off, set_ice_albedo
  USE mo_ice_fem_interface,   ONLY: ice_fem_interface
  USE mo_ice_advection,       ONLY: ice_advection_upwind

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ice_slow_interface
  PUBLIC :: ice_fast_interface
  PUBLIC :: ice_fast

  CHARACTER(len=12)           :: str_module    = 'IceInterface'  ! Output of module for 1 line debug

CONTAINS

  !-------------------------------------------------------------------------
  !
  !>
  !! ice_slow_interface: calls slow sea ice thermodynamics and dynamics
  !!
  !! This function changes:
  !! p_ice      slow-thermodynamics and dynamics fields of sea ice
  !! p_oce_sfc  heat and fresh-water fluxes, passed to the ocean
  !!
  !! @par Revision History
  !! Initial release by Vladimir Lapin, MPI-M (2016-11)
  !
!<Optimize_Used>
  SUBROUTINE ice_slow_interface(p_patch_3D, p_ice, p_oce_sfc, atmos_fluxes, p_os, p_as, p_op_coeff)

    TYPE(t_patch_3D ),TARGET,   INTENT(IN)      :: p_patch_3D
    TYPE(t_sea_ice),            INTENT(INOUT)   :: p_ice
    TYPE(t_ocean_surface),      INTENT(INOUT)   :: p_oce_sfc
    TYPE(t_atmos_fluxes),       INTENT(IN)      :: atmos_fluxes
    TYPE(t_atmos_for_ocean),    INTENT(INOUT)   :: p_as
    TYPE(t_hydro_ocean_state),  INTENT(IN)      :: p_os
    TYPE(t_operator_coeff),     INTENT(IN)      :: p_op_coeff

    ! Local variables
    TYPE(t_patch),  POINTER :: p_patch

    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('bef.icedyn: hi  ',p_ice%hi    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.icedyn: hs  ',p_ice%hs    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.icedyn: conc',p_ice%conc  ,str_module,2, in_subset=p_patch%cells%owned)

    !---------------------------------------------------------------------
    ! (1) -------------------- Sea ice dynamics --------------------------
    !---------------------------------------------------------------------
    IF ( i_ice_dyn >= 1 ) THEN

      IF (timers_level > 1) CALL timer_start(timer_extra40)

      ! solve for ice velocities (AWI FEM model wrapper)
      CALL ice_fem_interface ( p_patch_3D, p_ice, p_os, p_as, atmos_fluxes, p_op_coeff, p_oce_sfc)

      ! advection
      IF (i_ice_advec == 0) THEN
        IF (ltimer) CALL timer_start(timer_ice_advection)

        CALL ice_advection_upwind( p_patch_3D, p_op_coeff, p_ice )

        IF (ltimer) CALL timer_stop(timer_ice_advection)
      ELSEIF (i_ice_advec == 1) THEN ! nothing to do here
        ! advection on FEM grid is done inside ice_fem_interface
      ENDIF

      ! fix possible overshoots/undershoots after advection (previously, ice_clean_up_dyn)
      CALL ice_cut_off( p_patch, p_ice )

      IF (timers_level > 1) CALL timer_stop(timer_extra40)

    ELSE
      p_ice%u = 0._wp
      p_ice%v = 0._wp
    ENDIF
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('bef.icethm: hi  ',p_ice%hi    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.icethm: hs  ',p_ice%hs    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.icethm: conc',p_ice%conc  ,str_module,2, in_subset=p_patch%cells%owned)

    !---------------------------------------------------------------------
    ! (2) --------------- Slow sea ice thermodynamics --------------------
    !---------------------------------------------------------------------
    CALL ice_slow_thermo(p_patch_3D, p_os, atmos_fluxes, p_ice, p_oce_sfc)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('aft.icethm: hi  ',p_ice%hi    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.icethm: hs  ',p_ice%hs    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.icethm: conc',p_ice%conc  ,str_module,2, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE ice_slow_interface

  !-------------------------------------------------------------------------
  !
  !>
  !! ice_fast_interface: calls fast sea ice thermodynamics in ocean_surface
  !! Note: in coupled runs atmosphere calls ice_fast directly.
  !!
  !! This function changes:
  !! p_ice        - fast-thermodynamics fields (Qtop, Qbot, Tsurf)
  !! atmos_fluxes - ice albedos (albvisdir, albvisdif, albnirdir, albnirdif)
  !!
  !! @par Revision History
  !! Original code (mo_ocean_surface) by Stephan Lorenz, MPI-M (2015-04)
  !! Modified by Vladimir Lapin, MPI-M (2016-11)
  !
!<Optimize_Used>
  SUBROUTINE ice_fast_interface(p_patch, p_ice, atmos_fluxes, this_datetime)

    TYPE(t_patch), TARGET,      INTENT(IN)      :: p_patch
    TYPE(t_sea_ice),            INTENT(INOUT)   :: p_ice
    TYPE(t_atmos_fluxes),       INTENT(INOUT)   :: atmos_fluxes
    TYPE(datetime), POINTER,    INTENT(IN)      :: this_datetime
    !
    ! local variables
    INTEGER               :: jb, i_startidx_c, i_endidx_c

    TYPE(t_subset_range), POINTER :: all_cells

    !-----------------------------------------------------------------------
    all_cells       => p_patch%cells%all
    !---------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(jb, i_startidx_c, i_endidx_c) SCHEDULE(dynamic)
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
            &   p_ice% Tfw  (:,  jb),   &
            &   atmos_fluxes%albvisdir(:,:,jb), &  !  intent(out)
            &   atmos_fluxes%albvisdif(:,:,jb), &  !  intent(out)
            &   atmos_fluxes%albnirdir(:,:,jb), &  !  intent(out)
            &   atmos_fluxes%albnirdif(:,:,jb), &  !  intent(out)
            &   doy=getDayOfYearFromDateTime(this_datetime))
        ENDDO
!ICON_OMP_END_PARALLEL_DO

        ! provide constant heat fluxes for special analytical cases
        ! to-do: move to more appropriate place, should not be done here
        SELECT CASE (atmos_flux_analytical_type)
        CASE (102)
          p_ice%Qtop(:,1,:) = atmos_SWnet_const
          p_ice%Qbot(:,1,:) = 0.0_wp
        CASE (103)
          p_ice%Qtop(:,1,:) = 0.0_wp
          p_ice%Qbot(:,1,:) = atmos_sens_const
        END SELECT

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('ice_fast: hi     ',p_ice%hi       ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('ice_fast: hs     ',p_ice%hs       ,str_module,3, in_subset=p_patch%cells%owned)

    CALL dbg_print('aft.fast: Tsurf  ',p_ice%tsurf    ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.fast: Qtop   ',p_ice%Qtop     ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.fast: Qbot   ',p_ice%Qbot     ,str_module,3, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE ice_fast_interface

  !-------------------------------------------------------------------------
  !>
  !! ! ice_fast: ice thermodynamics at atmospheric time-step.
  !!   Calculates ice/snow surface temp, air-ice fluxes and sets albedos.
  !!
  !! This function changes:
  !! Tsurf and T1, T2 (winton)  - temperature of snow/ice layer(s)
  !! Qtop, Qbot                 - heat flux available for surface/bottom melting
  !! alb{vis/nir}{dir/dif}      - albedos
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !
  SUBROUTINE ice_fast(i_startidx_c, i_endidx_c, nbdim, kice, pdtime, &
            &   Tsurf,          & ! Surface temperature [degC]
            &   T1,             & ! Temperature of upper layer [degC]
            &   T2,             & ! Temperature of lower layer [degC]
            &   hi,             & ! Ice thickness
            &   hs,             & ! Snow thickness
            &   Qtop,           & ! Energy flux available for surface melting [W/m2]
            &   Qbot,           & ! Energy flux available for bottom melting [W/m2]
            &   SWnet,          & ! Net shortwave flux [W/m^2]
            &   nonsolar,       & ! Latent and sensible heat flux and longwave radiation [W/m^2]
            &   dnonsolardT,    & ! Derivative of non-solar fluxes w.r.t. temperature [W/m^2/K]
            &   Tfw,            & ! Freezing temperature of the ocean
            &   albvisdir,      & ! Albedo VIS, direct/parallel
            &   albvisdif,      & ! Albedo VIS, diffuse
            &   albnirdir,      & ! Albedo NIR, direct/parallel
            &   albnirdif,      & ! Albedo NIR, diffuse
            &   doy)              ! Day of the year

    INTEGER, INTENT(IN)    :: i_startidx_c, i_endidx_c, nbdim, kice
    REAL(wp),INTENT(IN)    :: pdtime
    REAL(wp),INTENT(INOUT) :: Tsurf      (nbdim,kice)
    REAL(wp),INTENT(INOUT) :: T1         (nbdim,kice)
    REAL(wp),INTENT(INOUT) :: T2         (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hi         (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hs         (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qtop       (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qbot       (nbdim,kice)
    REAL(wp),INTENT(IN)    :: SWnet      (nbdim,kice)
    REAL(wp),INTENT(IN)    :: nonsolar   (nbdim,kice)
    REAL(wp),INTENT(IN)    :: dnonsolardT(nbdim,kice)
    REAL(wp),INTENT(IN)    :: Tfw        (nbdim)
    REAL(wp),INTENT(OUT)   :: albvisdir  (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: albvisdif  (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: albnirdir  (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: albnirdif  (nbdim,kice)

    INTEGER, OPTIONAL,INTENT(IN)  :: doy

    !-------------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_ice_fast)

    SELECT CASE (i_ice_therm)

    CASE (1)
      CALL set_ice_temp_zerolayer(i_startidx_c, i_endidx_c, nbdim, kice, pdtime, &
                            &   Tsurf, hi, hs, Qtop, Qbot, SWnet, nonsolar, dnonsolardT, Tfw)

    CASE (2)
      CALL set_ice_temp_winton(i_startidx_c, i_endidx_c, nbdim, kice, pdtime, &
                    &   Tsurf, T1, T2, hi, hs, Qtop, Qbot, SWnet, nonsolar, dnonsolardT, Tfw)

    CASE (3)
      IF ( .NOT. PRESENT(doy) ) THEN
        CALL finish(TRIM('mo_ice_interface:ice_fast'),'i_ice_therm = 3 not allowed in this context')
      ENDIF
      CALL set_ice_temp_zerolayer_analytical(i_startidx_c, i_endidx_c, nbdim, kice, &
            &   Tsurf, hi, hs, Qtop, Qbot, Tfw, doy)

    CASE (4)
      WHERE ( hi(:,:) > 0._wp )
      Tsurf=min(0._wp, Tsurf + (SWnet+nonsolar + ki/hi*(Tf-Tsurf)) &
        &               / (ci*rhoi*hci_layer/pdtime-dnonsolardT+ki/hi))
      ELSEWHERE
        Tsurf(:,:) = Tf
      ENDWHERE

    END SELECT

    ! New albedo based on the new surface temperature
    CALL set_ice_albedo(i_startidx_c, i_endidx_c, nbdim, kice, Tsurf, hi, hs, &
      & albvisdir, albvisdif, albnirdir, albnirdif)

    IF (ltimer) CALL timer_stop(timer_ice_fast)

   END SUBROUTINE ice_fast

END MODULE mo_ice_interface
