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
#ifndef _OPENMP
#include "consistent_fma.inc"
#endif
!----------------------------
MODULE mo_sea_ice
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
  USE mo_model_domain,        ONLY: t_patch
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: sea_boundary
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref, Tf,        &
    &                               mu, alf, clw, sice

  USE mo_ocean_nml,           ONLY: no_tracer
  USE mo_sea_ice_nml,         ONLY: hnull, hmin,                    &
    &                               leadclose_1, leadclose_2n,      &
    &                               use_constant_tfreez, t_heat_base
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_state,         ONLY: v_base

  USE mo_sea_ice_types,       ONLY: t_sea_ice
  USE mo_ocean_surface_types, ONLY: t_ocean_surface
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_dbg_nml,             ONLY: idbg_mxmn, idbg_val
  USE mo_fortran_tools,       ONLY: assign_if_present

  IMPLICIT NONE

  PRIVATE

  ! public subroutines

  PUBLIC :: ice_conc_change
  PUBLIC :: salt_content_in_surface
  PUBLIC :: energy_content_in_surface

  !to be put into namelist
  !  INTEGER :: i_no_ice_thick_class = 1

  CHARACTER(len=12)           :: str_module    = 'SeaIce'  ! Output of module for 1 line debug

CONTAINS



  !-------------------------------------------------------------------------------
  !
  !
  !>
  !! !! ice_conc_change: Calculates the changes in concentration as well as the grid-cell average
  !                     thickness of new ice forming in open-water areas
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !! Einar Olason, renamed and added support for changing concentration
  !!
  SUBROUTINE ice_conc_change(p_patch,ice, p_os)

    TYPE(t_patch),             INTENT(IN)    :: p_patch
    TYPE (t_sea_ice),          INTENT(INOUT) :: ice
    TYPE(t_hydro_ocean_state), INTENT(IN)    :: p_os

    INTEGER  :: k
 !  REAL(wp) :: sst(nproma,p_patch%alloc_cell_blocks)
 !  REAL(wp) :: sss(nproma,p_patch%alloc_cell_blocks)
    REAL(wp) :: Tfw(nproma,p_patch%alloc_cell_blocks) ! Ocean freezing temperature [C]

 !  REAL(wp) :: leadclose_2n

    CALL dbg_print('IceConcCh: IceConc beg' ,ice%conc, str_module, 4, in_subset=p_patch%cells%owned)

    ! Calculate the sea surface freezing temperature                        [C]
    IF ( no_tracer < 2 .OR. use_constant_tfreez ) THEN
      Tfw(:,:) = Tf
    ELSE
      Tfw(:,:) = -mu * p_os%p_prog(nold(1))%tracer(:,1,:,2)
    ENDIF

    ! This should not be needed
    ! TODO ram - remove all instances of p_patch%cells%area(:,:) and test
    ! See also dynamics_fem/mo_ice_fem_interface.f90
    DO k=1,ice%kice
      ice%vol (:,k,:) = ice%hi(:,k,:)*ice%conc(:,k,:)*p_patch%cells%area(:,:)
      ice%vols(:,k,:) = ice%hs(:,k,:)*ice%conc(:,k,:)*p_patch%cells%area(:,:)
    ENDDO

    CALL dbg_print('IceConcCh: vol  at beg' ,ice%vol , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: vols at beg' ,ice%vols, str_module, 4, in_subset=p_patch%cells%owned)

!ICON_OMP_PARALLEL
!ICON_OMP_WORKSHARE 
    ! Concentration change due to new ice formation
    WHERE ( ice%newice(:,:) > 0._wp .AND. v_base%lsm_c(:,1,:) <= sea_boundary )
      ! New volume - we just preserve volume:
      ! #slo# 2014-12-11: newice is grown over open ocean but already averaged over whole grid area
      !                   newice must not be multiplied by 1-conc
      ice%vol  (:,1,:) = ice%vol(:,1,:) + ice%newice(:,:)*p_patch%cells%area(:,:)

      ! Hibler's way to change the concentration 
      !  - the formulation here uses the default values of leadclose parameters 2 and 3 in MPIOM:
      !    1 and 0 respectively, which recovers the Hibler model: conc=conc+newice/hnull
      ! Fixed 2. April (2014) - we don't need to multiply with 1-A here, like Hibler does, because it's
      ! already included in newice (we use volume, but Hibler growth rate)
      !ice%conc (:,1,:) = min( 1._wp, ice%conc(:,1,:) + ice%newice(:,:)/hnull )

      ! New formulation of leadclose parameter leadclose_2n includes parameters 2 and 3 of MPIOM:
      ! leadclose_2n (=mpiom_leadclose(3)/mpiom_leadclose(2)
      ! standard value of mpiom is: mpiom_leadclose(3)=2. mpiom_leadclose(2)=mpiom_leadclose(3)+1.
      ! i.e. leadclose_2n=2./3. according to mpiom default
      ice%conc(:,1,:) = min( 1._wp, ice%conc(:,1,:) + &
        &                           ice%newice(:,:)/(hnull+leadclose_2n*(ice%hi(:,1,:)-hnull)) )

      ! New ice and snow thickness
      ice%hi   (:,1,:) = ice%vol (:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
      ice%hs   (:,1,:) = ice%vols(:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
      !TODO: Re-calculate temperatures to conserve energy when we change the ice thickness
    ENDWHERE
!ICON_OMP_END_WORKSHARE

#ifndef _OPENMP
    CALL dbg_print('IceConcCh: conc leadcl' ,ice%conc, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hi   leadcl' ,ice%hi  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hs   leadcl' ,ice%hs  , str_module, 4, in_subset=p_patch%cells%owned)
#endif

!ICON_OMP_WORKSHARE
    ! This is where concentration, and thickness change due to ice melt (we must conserve volume)
    ! A.k.a. lateral melt
    WHERE ( ice%hiold(:,1,:) > ice%hi(:,1,:) .AND. ice%hi(:,1,:) > 0._wp )
      ! Hibler's way to change the concentration due to lateral melting (leadclose parameter 1)
      ice%conc(:,1,:) = MAX( 0._wp, ice%conc(:,1,:) &
        &        - ( ice%hiold(:,1,:)-ice%hi(:,1,:) )*ice%conc(:,1,:)*leadclose_1/ice%hiold(:,1,:) )

      ! New ice and snow thickness
      ice%hi  (:,1,:) = ice%vol (:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
      ice%hs  (:,1,:) = ice%vols(:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
      !TODO: Re-calculate temperatures to conserve energy when we change the ice thickness
    ENDWHERE
!ICON_OMP_END_WORKSHARE

#ifndef _OPENMP
    CALL dbg_print('IceConcCh: conc latMlt' ,ice%conc, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hi   latMlt' ,ice%hi  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hs   latMlt' ,ice%hs  , str_module, 4, in_subset=p_patch%cells%owned)
#endif

    ! Ice cannot grow thinner than hmin
    ! Changed 27. March
!ICON_OMP_WORKSHARE
    WHERE ( ice%hi(:,1,:) < hmin .AND. ice%hi(:,1,:) > 0._wp )
      ice%hi  (:,1,:) = hmin
      ice%conc(:,1,:) = ice%vol(:,1,:) / ( ice%hi(:,1,:)*p_patch%cells%area(:,:) )
      ice%hs  (:,1,:) = ice%vols(:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
    ENDWHERE
!ICON_OMP_END_WORKSHARE

!ICON_OMP_WORKSHARE
    WHERE (ice%hi(:,1,:) <= 0._wp)
      ice%Tsurf(:,1,:) = Tfw(:,:)
      ice%T1   (:,1,:) = Tfw(:,:)
      ice%T2   (:,1,:) = Tfw(:,:)
      ice%conc (:,1,:) = 0.0_wp
      ice%hi   (:,1,:) = 0.0_wp
      ice%hs   (:,1,:) = 0.0_wp
      ice%E1   (:,1,:) = 0.0_wp
      ice%E2   (:,1,:) = 0.0_wp
      ice%vol  (:,1,:) = 0.0_wp
    ENDWHERE
!ICON_OMP_END_WORKSHARE
!ICON_OMP_END_PARALLEL

    ice%concSum(:,:)  = SUM(ice%conc(:,:,:),2)

    CALL dbg_print('IceConcCh: IceConc end' ,ice%conc, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hi   at end' ,ice%hi  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hs   at end' ,ice%hs  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: vol  at end' ,ice%vol , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: vols at end' ,ice%vols, str_module, 4, in_subset=p_patch%cells%owned)

  END SUBROUTINE ice_conc_change




  ! compute the salt content in the upper most layer based on the liquid water height from the ice model: zUnderIce
  FUNCTION salt_content_in_surface(p_patch, thickness, p_ice, p_os, surface_fluxes, zUnderIceOld,computation_type, info) &
      & RESULT(salt)
    TYPE(t_patch),POINTER                                 :: p_patch
    REAL(wp),DIMENSION(nproma,p_patch%alloc_cell_blocks), &
      & INTENT(IN)                                        :: thickness,zUnderIceOld
    TYPE (t_sea_ice),       INTENT(INOUT)                 :: p_ice
    TYPE(t_hydro_ocean_state)                             :: p_os
    TYPE(t_ocean_surface)                                 :: surface_fluxes
    INTEGER,INTENT(IN), OPTIONAL                          :: computation_type
    CHARACTER(len=*) , OPTIONAL                           :: info

    ! locals
    REAL(wp), DIMENSION(nproma,p_patch%alloc_cell_blocks) :: salt, salinityDiff
    REAL(wp), DIMENSION(nproma,p_patch%alloc_cell_blocks) :: saltInSeaice, saltInLiquidWater
    INTEGER                                               :: my_computation_type
    CHARACTER(len=20)                                     :: my_info
    TYPE(t_subset_range), POINTER                         :: subset
    INTEGER                                               :: block, cell, cellStart,cellEnd

    my_computation_type = 0
    my_info             = 'BEFORE'

    salinityDiff = 0.0_wp
    salt         = 0.0_wp

    ! exit if actual debug-level is < 2
    IF (idbg_mxmn < 2 .AND. idbg_val < 2) RETURN

    CALL assign_if_present(my_computation_type, computation_type)
    CALL assign_if_present(my_info, info)

    subset => p_patch%cells%owned
    DO block = subset%start_block, subset%end_block
      CALL get_index_range(subset, block, cellStart, cellEnd)
      DO cell = cellStart, cellEnd
        IF (subset%vertical_levels(cell,block) < 1) CYCLE
        SELECT CASE (my_computation_type)
        CASE (0)
          ! compute salt amount in the first layer
          saltInSeaice(cell,block)      = sice &
            &                    * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                    * p_patch%cells%area(cell,block)
          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * p_ice%zUnderIce(cell,block) &
            &                    * p_patch%cells%area(cell,block)
        CASE (1)
          ! compute salt amount in the first layer
          saltInSeaice(cell,block)      = sice &
            &                    * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                    * p_patch%cells%area(cell,block)
    !     saltInSeaice(cell,block)      = 0.0_wp
          salinityDiff(cell,block)      = surface_fluxes%FrshFlux_TotalSalt(cell,block)*(dtime/(thickness(cell,block) &
            &                                                                    + p_os%p_prog(nold(1))%h(cell,block)))
          p_os%p_prog(nold(1))%tracer(cell,1,block,2)  = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                                   + salinityDiff(cell,block)
          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * p_ice%zUnderIce(cell,block) &
            &                    * p_patch%cells%area(cell,block)
    !     saltInLiquidWater(cell,block) = 0.0_wp
        CASE (2)
          ! compute salt amount in the first layer
          saltInSeaice(cell,block)      = sice &
            &                    * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                    * p_patch%cells%area(cell,block)

          salinityDiff(cell,block)      = surface_fluxes%FrshFlux_TotalSalt(cell,block)*(dtime/p_ice%zUnderIce(cell,block))
          p_os%p_prog(nold(1))%tracer(cell,1,block,2)  = p_os%p_prog(nold(1))%tracer(cell,1,block,2) + salinityDiff(cell,block)

          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * p_ice%zUnderIce(cell,block) &
            &                    * p_patch%cells%area(cell,block)

        CASE (3) ! use zunderIce for volume in tracer change
          saltInSeaice(cell,block)      = sice*rhoi &
            &                    * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                    * p_patch%cells%area(cell,block)

          p_os%p_prog(nold(1))%tracer(cell,1,block,2) = (p_os%p_prog(nold(1))%tracer(cell,1,block,2)*zUnderIceOld(cell,block) &
            &                                            - dtime*surface_fluxes%FrshFlux_TotalSalt(cell,block)) &
            &                                           /p_ice%zUnderIce(cell,block)
          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * p_ice%zUnderIce(cell,block)*rho_ref &
            &                    * p_patch%cells%area(cell,block)
        CASE (4) ! use zunderIce for volume in tracer change, multiply flux with top layer salinity
          saltInSeaice(cell,block)      = sice*rhoi &
            &                    * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                    * p_patch%cells%area(cell,block)

          p_os%p_prog(nold(1))%tracer(cell,1,block,2) = (p_os%p_prog(nold(1))%tracer(cell,1,block,2)*zUnderIceOld(cell,block) &
            &                                            -   dtime &
            &                                              * surface_fluxes%FrshFlux_TotalSalt(cell,block) &
            &                                              * p_os%p_prog(nold(1))%tracer(cell,1,block,2)) &
            &                                           /p_ice%zUnderIce(cell,block)
          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * p_ice%zUnderIce(cell,block)*rho_ref &
            &                    * p_patch%cells%area(cell,block)
        CASE (5) ! use zunderIce for volume in tracer change, multiply flux with top layer salinity
          p_ice%zUnderIce(cell,block) = zUnderIceOld(cell,block)
          saltInSeaice(cell,block)      = sice*rhoi &
            &                    * SUM(p_ice%hi(cell,:,block)*p_ice%conc(cell,:,block)) &
            &                    * p_patch%cells%area(cell,block)
        !!DN This is no longer needed since we now update surface salinity
        !directly
        !!DN   p_os%p_prog(nold(1))%tracer(cell,1,block,2) = (p_os%p_prog(nold(1))%tracer(cell,1,block,2)*zUnderIceOld(cell,block) &
        !!DN   &                                            -   dtime &
        !!DN   &                                              * surface_fluxes%FrshFlux_TotalSalt(cell,block) &
        !!DN   &                                              * p_os%p_prog(nold(1))%tracer(cell,1,block,2)) &
        !!DN   &                                           /p_ice%zUnderIce(cell,block)
          saltInLiquidWater(cell,block) = p_os%p_prog(nold(1))%tracer(cell,1,block,2) &
            &                    * p_ice%zUnderIce(cell,block)*rho_ref &
            &                    * p_patch%cells%area(cell,block)
        END SELECT

        salt(cell,block) = saltInSeaice(cell,block) + saltInLiquidWater(cell,block)
      END DO
    END DO

     CALL dbg_print('IceBudget: saltinIce '//TRIM(info)  , &
      &            saltInSeaice , &
      &            str_module, 5, in_subset=p_patch%cells%owned)
     CALL dbg_print('IceBudget: saltinLiquid '//TRIM(info)  , &
      &            saltInLiquidWater , &
      &            str_module, 5, in_subset=p_patch%cells%owned)
     CALL dbg_print('IceBudget: salt '//TRIM(info)  , &
      &            salt , &
      &            str_module, 5, in_subset=p_patch%cells%owned)
     CALL dbg_print('IceBudget: salinityDiff '//TRIM(info)  , &
      &            salinityDiff , &
      &            str_module, 5, in_subset=p_patch%cells%owned)
     CALL dbg_print('IceBudget: zUnderice '//TRIM(info)  , &
      &            p_ice%zUnderIce , &
      &            str_module, 5, in_subset=p_patch%cells%owned)
    !
    ! compute liquid volume in the first layer incl. water prepresentative of sea ice
  END FUNCTION salt_content_in_surface

  ! compute the energy content in the upper most layer based on the liquid water height from the ice model: zUnderIce
  FUNCTION energy_content_in_surface(p_patch, thickness, hold, p_ice, sst, computation_type, info) &
    & RESULT(energy)

    TYPE(t_patch),POINTER                                 :: p_patch
    REAL(wp),DIMENSION(nproma,p_patch%alloc_cell_blocks), INTENT(IN) :: thickness, hold, sst
    TYPE (t_sea_ice), INTENT(IN)                          :: p_ice
    INTEGER,INTENT(IN), OPTIONAL                          :: computation_type
    CHARACTER(len=*) , OPTIONAL                           :: info

    ! locals
    REAL(wp), DIMENSION(nproma,p_patch%alloc_cell_blocks) :: energy
    INTEGER                                               :: my_computation_type
    CHARACTER(len=20)                                     :: my_info
    TYPE(t_subset_range), POINTER                         :: subset
    INTEGER                                               :: block, cell, cellStart,cellEnd
    REAL(wp)                                              :: t_base, zui

    my_computation_type = 0
    my_info             = 'BEFORE'
    energy              = 0.0_wp
  ! t_base              = Tf
  ! t_base              = -5.0_wp
    t_base              = t_heat_base  !  arbitrary temperature basis for calculation of surface heat content

    ! exit if actual debug-level is < 2
    IF (idbg_mxmn < 2 .AND. idbg_val < 2) RETURN

    CALL assign_if_present(my_computation_type, computation_type)
    CALL assign_if_present(my_info, info)

    subset => p_patch%cells%owned
    DO block = subset%start_block, subset%end_block
      CALL get_index_range(subset, block, cellStart, cellEnd)
      DO cell = cellStart, cellEnd
        IF (subset%vertical_levels(cell,block) < 1) CYCLE
        SELECT CASE (my_computation_type)
        CASE (0)
          ! compute energy content of surface layer plus melting energy of ice and snow water equivalent
          !  - relative to arbitrary temperature t_base (e.g. -5C for mostly positive values)
          !  - omit multiplication with area, calculation per unit area, units in Joule/m2
          !  - constant freezing temperature Tf, ice-temperature set to Tf
          !  - use zUnderIce+draftave
          !  = (sst-t_base)*zUnderIce*rhow*clw - draftave*rhow*alf + (Tf-t_base)*(hi*rhoi+hs*rhos)*rhow*clw
          energy(cell,block) = (sst(cell,block) - t_base) * p_ice%zUnderIce(cell,block)*rho_ref*clw &
            &                - (p_ice%draftave(cell,block)*rho_ref*alf) &
            &                + (Tf - t_base)*p_ice%draftave(cell,block)*rho_ref*clw
        CASE (1)
          !  compute energy content - use zUnderIce and hi, hs, conc
          !  = (sst-t_base)*zUnderIce*rhow*clw - (hi*rhoi+hs*rhos)*alf*conc + (Tf-t_base)*(hi*rhoi+hs*rhos)*conc*clw
          energy(cell,block) = (sst(cell,block) - t_base) * p_ice%zUnderIce(cell,block)*rho_ref*clw &
            &                - ((p_ice%hi(cell,1,block)*rhoi + p_ice%hs(cell,1,block)*rhos)*p_ice%conc(cell,1,block)*alf) &
            &                + (Tf - t_base)*(p_ice%hi(cell,1,block)*rhoi + p_ice%hs(cell,1,block)*rhos) &
            &                                *p_ice%conc(cell,1,block)*clw
        CASE (2)
          !  compute energy content - use hi, hs only, compute local zUnderIce
          !  = (sst-t_base)*zUnderIce*rhow*clw - (hi*rhoi+hs*rhos)*alf*conc + (Tf-t_base)*(hi*rhoi+hs*rhos)*conc*clw
          zui                = thickness(cell,block)+hold(cell,block) &
            &                - (rhos * p_ice%hs(cell,1,block) + rhoi * p_ice%hi(cell,1,block)) * p_ice%conc(cell,1,block) / rho_ref
          energy(cell,block) = (sst(cell,block) - t_base) *zui*rho_ref*clw &
            &                - ((p_ice%hi(cell,1,block)*rhoi + p_ice%hs(cell,1,block)*rhos)*p_ice%conc(cell,1,block)*alf) &
            &                + (Tf - t_base)*(p_ice%hi(cell,1,block)*rhoi + p_ice%hs(cell,1,block)*rhos) &
            &                               *p_ice%conc(cell,1,block)*clw
          CONTINUE
        CASE (3)
          write(0,*) " Nothing computed"
        CASE DEFAULT
          CALL finish ('mo_sea_ice:computation_type','option not supported')
        END SELECT
        !salt(cell,block) = saltInSeaice(cell,block) + saltInLiquidWater(cell,block)
      END DO
    END DO

    CALL dbg_print('enContSurf: energy '//TRIM(info),energy,str_module, 5, in_subset=p_patch%cells%owned)

  END FUNCTION energy_content_in_surface

END MODULE mo_sea_ice
