!! ===========================================================================================================================
!! Implementation of tides by computation of the Sun's and Moon's full tidal potential
!! This will be used in the pressure gradient calculation
!!
!! Authors: Kai Logemann, Helmholtz-Zentrum Geesthacht, and Leonidas Linardakis, Max-Planck-Institute for Meteorology, Hamburg
!!
!! ===========================================================================================================================

!----------------------------
#include "icon_definitions.inc"
#include "omp_definitions.inc"
#include "iconfor_dsl_definitions.inc"
!----------------------------
MODULE mo_ocean_pressure_bc_conditions
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp, dp
   USE mtime,                    ONLY: datetime
  USE mo_exception,              ONLY: message, finish, warning, message_text
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d
  USE mo_ocean_nml,              ONLY: n_zlev, use_tides,   &
    & use_tides_SAL, atm_pressure_included_in_ocedyn,       &
    & OceanReferenceDensity_inv
  USE mo_physical_constants,     ONLY: grav
  USE mo_grid_subset,            ONLY: t_subset_range, get_index_range
  USE mo_parallel_config,        ONLY: nproma
  USE mo_impl_constants,         ONLY: sea_boundary
  USE mo_ocean_math_operators,   ONLY: grad_fd_norm_oce_2d_3d
  USE mo_util_dbg_prnt,          ONLY: dbg_print, debug_printValue
  USE mo_ocean_tides,            ONLY: calculate_tides_potential  
  USE mo_ocean_surface_types,    ONLY: t_atmos_for_ocean
  USE mo_ocean_types,            ONLY: t_hydro_ocean_state
  USE mo_dynamics_config,        ONLY: nold, nnew


  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: create_pressure_bc_conditions
  
  CHARACTER(LEN=12)  :: str_module = 'mo_ocean_pressure_bc_conditions'  ! Output of module for 1 line debug
  !-------------------------------------------------------------------------
  
CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE create_pressure_bc_conditions(patch_3d,ocean_state, p_as,current_time)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    TYPE(t_atmos_for_ocean)                          :: p_as    
    TYPE(datetime), POINTER                          :: current_time
      
    !------------------------------------------------------------------------
    ! compute tidal potential
    CALL calculate_tides_potential(patch_3d,current_time,ocean_state%p_diag%rho, ocean_state%p_prog(nold(1))%h, &
      ocean_state%p_aux%bc_tides_potential, ocean_state%p_aux%bc_SAL_potential) 
    !------------------------------------------------------------------------
    ! total top potential
    IF (atm_pressure_included_in_ocedyn) THEN
      ocean_state%p_aux%bc_total_top_potential =  &
        & ocean_state%p_aux%bc_tides_potential &
        & + ocean_state%p_aux%bc_SAL_potential &
        & + p_as%pao * OceanReferenceDensity_inv
    
    ELSE IF (use_tides .or. use_tides_SAL) THEN
      ocean_state%p_aux%bc_total_top_potential =  ocean_state%p_aux%bc_tides_potential + ocean_state%p_aux%bc_SAL_potential
    ENDIF
                            
    CALL dbg_print('tides_potential',  ocean_state%p_aux%bc_tides_potential, &
      str_module, 3, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    CALL dbg_print('tides_SAL',      ocean_state%p_aux%bc_SAL_potential, &
      str_module, 3, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
        
  END SUBROUTINE create_pressure_bc_conditions
  !-------------------------------------------------------------------------
  

END MODULE mo_ocean_pressure_bc_conditions
!=============================================================================
