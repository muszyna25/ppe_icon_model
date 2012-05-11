MODULE mo_sea_ice_shared_sr
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
  USE mo_exception,           ONLY: finish, message
  USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell, sea_boundary 
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref,ki,ks,Tf,albi,albim,albsm,albs,&
    &                               mu,mus,ci, alf, I_0, alv, albedoW, clw,            &
    &                               cpd, zemiss_def,rd, stbo,tmelt   
  USE mo_math_constants,      ONLY: rad2deg
  USE mo_ocean_nml,           ONLY: no_tracer, init_oce_prog, iforc_oce, &
    &                               FORCING_FROM_FILE_FLUX
  USE mo_oce_state,           ONLY: t_hydro_ocean_state, v_base, ocean_var_list
  USE mo_oce_index,           ONLY: print_mxmn, ipl_src
  USE mo_var_list,            ONLY: add_var
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants
  ! # achim
  USE mo_sea_ice_types, ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: oce_ice_heatflx

CONTAINS

  !-------------------------------------
  !
  ! oce_ice_heatflx
  !
  ! Calculates the heat flux from the uppermost water layer into the ice.
  !
  ! Currently (as in growth.f90): all energy available in upper ocean grid cell 
  ! is supplied to the ice and the upper ocean temperature is held at the 
  ! freezing point. This is not very physical.
  !
  ! Positive flux upwards.
 
  
  SUBROUTINE oce_ice_heatflx (p_os,ice,Tfw,zHeatOceI)
    TYPE(t_hydro_ocean_state), INTENT(IN) :: p_os
    TYPE(t_sea_ice),           INTENT(IN) :: ice
    REAL(wp),                  INTENT(IN):: Tfw(:,:,:) ! freezing temperature
    REAL(wp),                  INTENT(OUT):: zHeatOceI(:,:,:)

    ! Local
    INTEGER :: k ! counter for ice thickness categories

    
    ! calculate heat flux from ocean to ice  (zHeatOceI) 
    DO k=1,ice%kice
      WHERE (ice%isice(:,k,:)) 
        zHeatOceI(:,k,:) = ( p_os%p_prog(nold(1))%tracer(:,1,:,1) - Tfw(:,k,:) ) &
          &                 * ice%zUnderIce(:,:) * clw*rho_ref/dtime
      ENDWHERE
    END DO
  END SUBROUTINE oce_ice_heatflx

END MODULE mo_sea_ice_shared_sr
