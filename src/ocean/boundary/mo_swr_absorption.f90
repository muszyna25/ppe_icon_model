!>
!! SWR-ABSORPTION SCHEME
!! Shortwave radiation is done in three steps
!! As first step the relative swr absorption factor (swrab) in each level is calculated
!! either for
!! a) a static jerlov type of water (default Type IB)
!!    sbr jerlov_swr_absorption
!!
!! The actual absorption and the update of the temperature is afterwards done
!! again in two steps.
!! a) for the uppermost level together with the rest of the
!!    thermodynamic forcing.
!! b) for the subsurface levels in subsurface_swr_absorption
!!    Heat is distributed over the water column in proportion to the absorption.
!!    No Heat is transpered into sea floor, but is added to the bottom level.
!!
!! @author H.Haak, MPI-M
!! @author U.Mikolajewicz, MPI-M
!!
!!
!! @par Revision History
!! Introduced by H.Haak, MPI-M (2019-06-13)
!!
!!
MODULE mo_swr_absorption
  USE mo_kind,                      ONLY: wp
  USE mo_ocean_nml,                 ONLY: jerlov_atten, jerlov_bluefrac, n_zlev  
  USE mo_physical_constants,        ONLY: rho_ref,  clw
  USE mo_sync,                      ONLY: global_sum_array, enable_sync_checks, disable_sync_checks

  IMPLICIT NONE

  PRIVATE

! Jerlov Water Type values from Kara et. al. (2005)
! jerlov_atten is taken from figure 2 ; jerlov_bluefrac
! is taken from figure 3
!
!   - namelist parameter
!   - default Type IB
!
! type   jerlov_atten [m-1] jerlov_bluefrac [%]
! --------------------------------------------
!  I       0.05                   0.45
!  IA      0.06                   0.41
!  IB      0.08                   0.36
!  II      0.12                   0.28
!  III     0.17                   0.27




  PUBLIC :: jerlov_swr_absorption                      &
       , subsurface_swr_absorption                  &
       , jerlov_atten,jerlov_bluefrac,check_hc2


CONTAINS


  !>
  !! This is the first part of the sw-absorption scheme
  !! ( can to be called before the time loop )
  !! This sbr calculates the relative swr absorption factor (swrab) in each level.
  !!
  !! Numbers are derived for jerlov  type of water
  !!    N.G. Jerlov: Optical Studies of Ocean Waters.
  !!             In: Reports of the Swedish Deep-Sea Expedition, 3:1-59; 1951.


  SUBROUTINE  jerlov_swr_absorption(patch_3d, ocean_state)

    USE mo_model_domain,              ONLY: t_patch, t_patch_3d
    USE mo_ocean_types,               ONLY: t_hydro_ocean_state
    USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
    USE mo_run_config,                ONLY: dtime

    TYPE(t_patch_3d ),TARGET, INTENT(in)              :: patch_3D
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout)  :: ocean_state

    REAL(wp), POINTER :: swsum(:,:)
    REAL(wp), POINTER :: swrab(:,:,:)

    INTEGER  :: blockNo, jc, start_index, end_index, level
  
    TYPE(t_patch), POINTER                   :: patch_2d
    TYPE(t_subset_range), POINTER            :: all_cells


    patch_2d => patch_3D%p_patch_2d(1)
    all_cells => patch_2d%cells%all
 
    swsum => ocean_state%p_diag%swsum
    swrab => ocean_state%p_diag%swrab

    
    !ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index) SCHEDULE(dynamic)
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc =  start_index, end_index

        swsum(jc,blockNo)=jerlov_bluefrac*EXP(-patch_3d%p_patch_1d(1)%zlev_i(2)*jerlov_atten)

        DO level=1,n_zlev

          swrab(jc,level,blockNo)= patch_3d%wet_c(jc,1,blockNo) * &
                (1._wp / swsum(jc,blockNo)) * jerlov_bluefrac * &
               ( EXP(-patch_3d%p_patch_1d(1)%zlev_i(level)*jerlov_atten) &
               - EXP(-patch_3d%p_patch_1d(1)%zlev_i(level+1)*jerlov_atten) )

        END DO ! level

        swsum(jc,blockNo)=swsum(jc,blockNo) * patch_3d%wet_c(jc,1,blockNo)

      END DO
    END DO

  END SUBROUTINE jerlov_swr_absorption

  SUBROUTINE subsurface_swr_absorption(patch_3d, ocean_state)
! !
! !    !>
! !    !! This is part two of the sw-absorption scheme and should be called after
! !    !! the surface thermodynamics -> growth )
! !    !! This sbr calculates the subsurface absorption from the botton to level 2.
! !    !! Absorption in the surface level is done in sbr growth.
! !    !! SWR does not penetrate into Land. Heat from SWR that would penetrate theoreticly
! !    !! below the bottom is added to the bottom layer.
! !

    USE mo_model_domain,              ONLY: t_patch, t_patch_3d
    USE mo_ocean_types,               ONLY: t_hydro_ocean_state
    USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
    USE mo_impl_constants,            ONLY: sea_boundary
    USE mo_dynamics_config,           ONLY: nold
    USE mo_run_config,                ONLY: dtime
    USE mo_parallel_config,           ONLY: nproma


    TYPE(t_patch_3d ),TARGET, INTENT(in)              :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout)  :: ocean_state

    REAL(wp), POINTER :: swsum(:,:)
    REAL(wp), POINTER :: swrab(:,:,:)
    REAL(wp), POINTER :: rsdoabsorb(:,:,:)
    REAL(wp), POINTER :: heatabs(:,:)

    REAL(wp) :: heatabb(nproma, patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: heatabs_t(nproma, patch_3D%p_patch_2d(1)%alloc_cell_blocks)

    REAL(wp) :: dti, cc

    INTEGER  :: blockNo, jc, start_index, end_index, level

    TYPE(t_patch), POINTER                   :: patch_2d
    TYPE(t_subset_range), POINTER            :: all_cells


    patch_2d => patch_3D%p_patch_2d(1)
    all_cells => patch_2d%cells%all

    swsum => ocean_state%p_diag%swsum
    swrab => ocean_state%p_diag%swrab
    heatabs => ocean_state%p_diag%heatabs
    rsdoabsorb => ocean_state%p_diag%rsdoabsorb

    cc = clw * rho_ref
    dti=1.0_wp/dtime

    !ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index,hetabs_t,level,heatabb) SCHEDULE(dynamic)
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc =  start_index, end_index

        heatabs_t(jc,blockno)=patch_3d%wet_c(jc,1,blockNo)*heatabs(jc,blockNo)*dtime/cc
        heatabb(jc,blockno) = 0._wp

       DO level=n_zlev,2,-1

          heatabb(jc,blockno)=heatabb(jc,blockno)+swrab(jc,level,blockNo)*heatabs_t(jc,blockno)

          ! diagnostic for subfurface sw absorbtion
          rsdoabsorb(jc,level,blockNo)=heatabb(jc,blockno)*cc*dti


          ocean_state%p_prog(nold(1))%tracer(jc,level,blockNo,1) = &
               ocean_state%p_prog(nold(1))%tracer(jc,level,blockNo,1) + &
                patch_3d%wet_c(jc,level,blockNo) &    
           * (heatabb(jc,blockno)/patch_3d%p_patch_1d(1)%prism_thick_c(jc,level,blockNo))

          heatabb(jc,blockno) = heatabb(jc,blockno) * ( 1.0_wp - patch_3d%wet_c(jc,level,blockNo))

        END DO

      END DO
    END DO

  END SUBROUTINE subsurface_swr_absorption




  SUBROUTINE check_hc2(patch_3d, ocean_state, p_oce_sfc, hc, hcs)

    USE mo_ocean_surface_types, ONLY: t_ocean_surface
    USE mo_model_domain,              ONLY: t_patch, t_patch_3d
    USE mo_ocean_types,               ONLY: t_hydro_ocean_state
    USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
    USE mo_impl_constants,            ONLY: sea_boundary
    USE mo_dynamics_config,           ONLY: nold
    USE mo_run_config,                ONLY: dtime


    TYPE(t_patch_3d ),TARGET, INTENT(in)              :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout)  :: ocean_state
    TYPE(t_ocean_surface),TARGET, INTENT(in)          :: p_oce_sfc   
    
    REAL(wp) :: cc
    REAL(wp) :: hc(:,:), hcs


    INTEGER  :: blockNo, jc, start_index, end_index, level
 
    TYPE(t_patch), POINTER                   :: patch_2d
    TYPE(t_subset_range), POINTER            :: owned_cells


    patch_2d => patch_3D%p_patch_2d(1)
    owned_cells => patch_2d%cells%owned

    cc = clw * rho_ref

    hc(:,:)=0.0_wp

    DO blockNo = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, blockNo, start_index, end_index)
      DO jc =  start_index, end_index

          hc(jc,blockno) = patch_3D%wet_c(jc,1,blockNo) * cc &
               * ocean_state%p_prog(nold(1))%tracer(jc,1,blockNo,1) &
               * p_oce_sfc%cellThicknessUnderIce(jc,blockno) &
               * patch_2d%cells%area(jc,blockNo)

          DO level = 2, patch_3d%p_patch_1d(1)%dolic_c(jc,blockno)

          hc(jc,blockno) = hc(jc,blockno) &
               + ( patch_3D%wet_c(jc,level,blockNo) * cc &
                 * ocean_state%p_prog(nold(1))%tracer(jc,level,blockNo,1) &
                 * patch_3d%p_patch_1d(1)%prism_thick_c(jc,level,blockNo) &
                 * patch_2d%cells%area(jc,blockNo) )


        END DO

      END DO
    END DO

    CALL disable_sync_checks()
 
    hcs=0.0_wp
    hcs=global_sum_array(hc(:,:))

    CALL enable_sync_checks()

  END SUBROUTINE check_hc2



END MODULE mo_swr_absorption

