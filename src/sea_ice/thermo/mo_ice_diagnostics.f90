!>
!! Provide an implementation of the sea-ice model diagnostics.
!!
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
!!
!----------------------------
MODULE mo_ice_diagnostics
  
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_exception,           ONLY: finish
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_dbg_nml,             ONLY: idbg_mxmn, idbg_val
  USE mo_fortran_tools,       ONLY: assign_if_present

  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_ocean_state,         ONLY: v_base

  USE mo_statistics,          ONLY: add_fields

  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref, clw, alf, sice, tf
  USE mo_sea_ice_nml,         ONLY: t_heat_base
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sea_ice_acc

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: salt_in_surface
  PUBLIC :: energy_in_surface
  PUBLIC :: update_ice_statistic, compute_mean_ice_statistics, reset_ice_statistics

  CHARACTER(len=12)           :: str_module = 'IceDiag'  ! Output of module for 1 line debug

CONTAINS

  !-------------------------------------------------------------------------------
  !>
  !! ! salt_in_surface: calculate salt content in the surface layer
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2009)
  FUNCTION salt_in_surface(p_patch, p_ice, sss, computation_type, dbg_lev, info) &
      & RESULT(salt)

    TYPE(t_patch),POINTER, INTENT(IN)                     :: p_patch
    TYPE (t_sea_ice), INTENT(IN)                          :: p_ice
    REAL(wp),DIMENSION(nproma,p_patch%alloc_cell_blocks), INTENT(IN) :: sss
    INTEGER,INTENT(IN), OPTIONAL                          :: computation_type, dbg_lev
    CHARACTER(len=*) , OPTIONAL                           :: info

    ! locals
    REAL(wp), DIMENSION(nproma,p_patch%alloc_cell_blocks) :: salt, saltInSeaice, saltInOcean
    TYPE(t_subset_range), POINTER                         :: all_cells
    REAL(wp), POINTER                                     :: area(:,:)

    INTEGER                                               :: my_computation_type, my_dbg_lev
    CHARACTER(len=20)                                     :: my_info
    INTEGER                                               :: jb, jc, i_startidx_c,i_endidx_c

    !-----------------------------------------------------------------------
    my_dbg_lev          = 4
    CALL assign_if_present(my_dbg_lev, dbg_lev)
    ! exit if actual debug-level < my_dbg_lev
    IF (idbg_mxmn < my_dbg_lev .AND. idbg_val < my_dbg_lev) RETURN

    my_computation_type = 0
    my_info             = 'salt_in_surface'
    CALL assign_if_present(my_computation_type, computation_type)
    CALL assign_if_present(my_info, info)

    all_cells    => p_patch%cells%all
    area         => p_patch%cells%area(:,:)

    !-----------------------------------------------------------------------

    salt(:,:)           = 0.0_wp

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF (all_cells%vertical_levels(jc,jb) < 1) CYCLE
        SELECT CASE (my_computation_type)

        CASE (0)
          saltInSeaice(jc,jb)   = sice * SUM(p_ice%hi(jc,:,jb)*p_ice%conc(jc,:,jb)) * area(jc,jb)
          saltInOcean(jc,jb)    = sss(jc,jb) * p_ice%zUnderIce(jc,jb)               * area(jc,jb)

!        CASE (3) ! use zunderIce for volume in tracer change
!          saltInSeaice(jc,jb)   = sice*rhoi * SUM(p_ice%hi(jc,:,jb)*p_ice%conc(jc,:,jb)) * area(jc,jb)
!
!          sss(jc,jb) = (sss(jc,jb)*zUnderIceOld(jc,jb) - dtime*p_oce_sfc%FrshFlux_TotalSalt(jc,jb)) &
!            &                                           /p_ice%zUnderIce(jc,jb)
!          saltInOcean(jc,jb)    = sss(jc,jb) * p_ice%zUnderIce(jc,jb)*rho_ref * area(jc,jb)

        CASE (5) ! use zunderIce for volume in tracer change, multiply flux with top layer salinity
          saltInSeaice(jc,jb)   = sice*rhoi * SUM(p_ice%hi(jc,:,jb)*p_ice%conc(jc,:,jb)) * area(jc,jb)
          saltInOcean(jc,jb)    = sss(jc,jb) * p_ice%zUnderIce(jc,jb)*rho_ref * area(jc,jb)

        END SELECT

        salt(jc,jb) = saltInSeaice(jc,jb) + saltInOcean(jc,jb)

      END DO
    END DO

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('IceBudget: saltinIce    ', saltInSeaice,    str_module,5, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceBudget: saltInOcean  ', saltInOcean,     str_module,5, in_subset=p_patch%cells%owned)
    CALL dbg_print(TRIM(my_info), salt, str_module, my_dbg_lev, in_subset=all_cells)
    !---------------------------------------------------------------------

  END FUNCTION salt_in_surface

  !-------------------------------------------------------------------------------
  !>
  !! ! energy_in_surface: calculate energy content in the surface layer
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2009)
  FUNCTION energy_in_surface(p_patch, p_ice, ssh, sst, computation_type, dbg_lev, info) &
    & RESULT(energy)

    TYPE(t_patch),POINTER                                 :: p_patch
    TYPE (t_sea_ice), INTENT(IN)                          :: p_ice
    REAL(wp),DIMENSION(nproma,p_patch%alloc_cell_blocks), INTENT(IN) :: ssh, sst
    INTEGER,INTENT(IN), OPTIONAL                          :: computation_type, dbg_lev
    CHARACTER(len=*) , OPTIONAL                           :: info

    ! locals
    TYPE(t_subset_range), POINTER                         :: all_cells
    REAL(wp), DIMENSION(nproma,p_patch%alloc_cell_blocks) :: energy
    INTEGER                                               :: my_computation_type, my_dbg_lev
    CHARACTER(len=20)                                     :: my_info
    INTEGER                                               :: jb, jc, i_startidx_c,i_endidx_c
    REAL(wp)                                              :: zUnderIce, prism_thick_flat

    !-----------------------------------------------------------------------
    my_dbg_lev          = 4
    CALL assign_if_present(my_dbg_lev, dbg_lev)
    ! exit if actual debug-level < my_dbg_lev
    IF (idbg_mxmn < my_dbg_lev .AND. idbg_val < my_dbg_lev) RETURN

    my_computation_type = 0
    my_info             = 'energy_in_surface'
    CALL assign_if_present(my_computation_type, computation_type)
    CALL assign_if_present(my_info, info)

    all_cells    => p_patch%cells%all
    prism_thick_flat = v_base%del_zlev_m(1) ! thickness of the top layer
    !-----------------------------------------------------------------------

    energy(:,:) = 0.0_wp

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF (all_cells%vertical_levels(jc,jb) < 1) CYCLE
        SELECT CASE (my_computation_type)
        CASE (0)
          ! compute energy content of surface layer plus melting energy of ice and snow water equivalent
          !  - relative to arbitrary temperature t_heat_base (e.g. -5C for mostly positive values)
          !  - omit multiplication with area, calculation per unit area, units in Joule/m2
          !  - constant freezing temperature Tf, ice-temperature set to Tf
          !  = (sst-t_heat_base)*zUnderIce*rhow*clw - draftave*rhow*alf + (Tf-t_heat_base)*(hi*rhoi+hs*rhos)*rhow*clw
          energy(jc,jb) = (sst(jc,jb) - t_heat_base) * p_ice%zUnderIce(jc,jb)*rho_ref*clw &
            &                - (p_ice%draftave(jc,jb)*rho_ref*alf) &
            &                + (Tf - t_heat_base)*p_ice%draftave(jc,jb)*rho_ref*clw
        CASE (1)
          !  compute energy content - use zUnderIce and hi, hs, conc
          !  = (sst-t_heat_base)*zUnderIce*rhow*clw - (hi*rhoi+hs*rhos)*alf*conc + (Tf-t_heat_base)*(hi*rhoi+hs*rhos)*conc*clw
          energy(jc,jb) = (sst(jc,jb) - t_heat_base) * p_ice%zUnderIce(jc,jb)*rho_ref*clw &
            &                - ((p_ice%hi(jc,1,jb)*rhoi + p_ice%hs(jc,1,jb)*rhos)*p_ice%conc(jc,1,jb)*alf) &
            &                + (Tf - t_heat_base)*(p_ice%hi(jc,1,jb)*rhoi + p_ice%hs(jc,1,jb)*rhos) &
            &                                *p_ice%conc(jc,1,jb)*clw
        CASE (2)
          !  compute energy content - use hi, hs only, compute local zUnderIce
          !  = (sst-t_heat_base)*zUnderIce*rhow*clw - (hi*rhoi+hs*rhos)*alf*conc + (Tf-t_heat_base)*(hi*rhoi+hs*rhos)*conc*clw
          zUnderIce     = prism_thick_flat+ssh(jc,jb) &
            &                - (rhos * p_ice%hs(jc,1,jb) + rhoi * p_ice%hi(jc,1,jb)) * p_ice%conc(jc,1,jb) / rho_ref
          energy(jc,jb) = (sst(jc,jb) - t_heat_base) *zUnderIce*rho_ref*clw &
            &                - ((p_ice%hi(jc,1,jb)*rhoi + p_ice%hs(jc,1,jb)*rhos)*p_ice%conc(jc,1,jb)*alf) &
            &                + (Tf - t_heat_base)*(p_ice%hi(jc,1,jb)*rhoi + p_ice%hs(jc,1,jb)*rhos) &
            &                               *p_ice%conc(jc,1,jb)*clw
        CASE DEFAULT
          CALL finish ('mo_sea_ice:computation_type','option not supported')
        END SELECT
      END DO
    END DO

    CALL dbg_print(TRIM(my_info), energy, str_module, my_dbg_lev, in_subset=all_cells)

  END FUNCTION energy_in_surface


  !-------------------------------------
  !
  ! Sea ice statistics
  !

    SUBROUTINE update_ice_statistic(p_acc, p_ice, subset)
    TYPE(t_sea_ice_acc),  INTENT(INOUT) :: p_acc
    TYPE(t_sea_ice),      INTENT(IN)    :: p_ice
    TYPE(t_subset_range), INTENT(IN)    :: subset

    CALL add_fields(p_acc%hi  , p_ice%hi  , subset , levels=p_ice%kice)
    CALL add_fields(p_acc%hs  , p_ice%hs  , subset , levels=p_ice%kice)
    CALL add_fields(p_acc%conc, p_ice%conc, subset , levels=p_ice%kice)
    CALL add_fields(p_acc%u   , p_ice%u   , subset)
    CALL add_fields(p_acc%v   , p_ice%v   , subset)
    CALL add_fields(p_acc%HeatOceW , p_ice%HeatOceW , subset)
    CALL add_fields(p_acc%zUnderIce, p_ice%zUnderIce, subset)
    CALL add_fields(p_acc%draftave , p_ice%draftave , subset)
    CALL add_fields(p_acc%Qtop     , p_ice%Qtop     , subset , levels=p_ice%kice)
    CALL add_fields(p_acc%Qbot     , p_ice%Qbot     , subset , levels=p_ice%kice)
    CALL add_fields(p_acc%CondHeat , p_ice%CondHeat , subset , levels=p_ice%kice)
    CALL add_fields(p_acc%HeatOceI , p_ice%HeatOceI , subset , levels=p_ice%kice)
    CALL add_fields(p_acc%zHeatOceI, p_ice%zHeatOceI, subset , levels=p_ice%kice)
  END SUBROUTINE update_ice_statistic

  SUBROUTINE compute_mean_ice_statistics(p_acc,nsteps_since_last_output)
    TYPE(t_sea_ice_acc), INTENT(INOUT) :: p_acc
    INTEGER,INTENT(IN)                 :: nsteps_since_last_output
    p_acc%hi                        = p_acc%hi       /REAL(nsteps_since_last_output,wp)
    p_acc%hs                        = p_acc%hs       /REAL(nsteps_since_last_output,wp)
    p_acc%u                         = p_acc%u        /REAL(nsteps_since_last_output,wp)
    p_acc%v                         = p_acc%v        /REAL(nsteps_since_last_output,wp)
    p_acc%conc                      = p_acc%conc     /REAL(nsteps_since_last_output,wp)
    p_acc%Qtop                      = p_acc%Qtop     /REAL(nsteps_since_last_output,wp)
    p_acc%Qbot                      = p_acc%Qbot     /REAL(nsteps_since_last_output,wp)
    p_acc%CondHeat                  = p_acc%CondHeat /REAL(nsteps_since_last_output,wp)
    p_acc%HeatOceI                  = p_acc%HeatOceI /REAL(nsteps_since_last_output,wp)
    p_acc%HeatOceW                  = p_acc%HeatOceW /REAL(nsteps_since_last_output,wp)
    p_acc%zHeatOceI                 = p_acc%zHeatOceI/REAL(nsteps_since_last_output,wp)
    p_acc%zUnderIce                 = p_acc%zUnderIce/REAL(nsteps_since_last_output,wp)
    p_acc%draftave                  = p_acc%draftave /REAL(nsteps_since_last_output,wp)
  END SUBROUTINE compute_mean_ice_statistics

  SUBROUTINE reset_ice_statistics(p_acc)
    TYPE(t_sea_ice_acc), INTENT(INOUT) :: p_acc
    p_acc%hi                        = 0.0_wp
    p_acc%hs                        = 0.0_wp
    p_acc%u                         = 0.0_wp
    p_acc%v                         = 0.0_wp
    p_acc%conc                      = 0.0_wp
    p_acc%Qtop                      = 0.0_wp
    p_acc%Qbot                      = 0.0_wp
    p_acc%CondHeat                  = 0.0_wp
    p_acc%HeatOceI                  = 0.0_wp
    p_acc%HeatOceW                  = 0.0_wp
    p_acc%zHeatOceI                 = 0.0_wp
    p_acc%zUnderIce                 = 0.0_wp
    p_acc%draftave                  = 0.0_wp
  END SUBROUTINE reset_ice_statistics

  END MODULE mo_ice_diagnostics
