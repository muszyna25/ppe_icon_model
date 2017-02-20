!>
!! Contains the implementation of the semi-implicit Adams-Bashforth timestepping
!! for the ICON ocean model.
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2010/04)
!!  Modified by Stephan Lorenz,     MPI-M (2010-06)
!!   - renaming and adjustment to ocean domain and patch
!!   - implementation of continuity equation for vertical velocities
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ocean_ab_timestepping
  USE mo_ocean_nml,                      ONLY: discretization_scheme
  USE mo_dynamics_config,                ONLY: nold, nnew
  USE mo_sea_ice_types,                  ONLY: t_sfc_flx
  USE mo_model_domain,                   ONLY: t_patch_3D !, t_patch
  USE mo_ext_data_types,                 ONLY: t_external_data
  USE mo_ocean_ab_timestepping_mimetic,  ONLY: solve_free_sfc_ab_mimetic, &
    &  calc_normal_velocity_ab_mimetic, &
    &  calc_vert_velocity_mim_bottomup
  USE mo_ocean_physics_types,            ONLY: t_ho_params
  USE mo_ocean_types,                    ONLY: t_hydro_ocean_state, t_operator_coeff, t_solverCoeff_singlePrecision
  USE mo_exception,                      ONLY: finish!, message_text
  
IMPLICIT NONE

PRIVATE

  INTEGER, PARAMETER :: MIMETIC_TYPE = 1
  INTEGER, PARAMETER :: RBF_TYPE     = 2

  PUBLIC :: solve_free_surface_eq_ab
  PUBLIC :: calc_normal_velocity_ab
  PUBLIC :: calc_vert_velocity
  PUBLIC :: update_time_indices

CONTAINS
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! !  Solves the free surface equation.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE solve_free_surface_eq_ab(patch_3D, ocean_state, external_data, surface_fluxes, &
    & physics_parameters, timestep, operators_coefficients, solverCoeff_sp, return_status)
    TYPE(t_patch_3D ),TARGET, INTENT(INOUT)   :: patch_3D
    TYPE(t_hydro_ocean_state), TARGET         :: ocean_state
    TYPE(t_external_data), TARGET             :: external_data
    TYPE(t_sfc_flx), INTENT(INOUT)            :: surface_fluxes
    TYPE (t_ho_params)                        :: physics_parameters
    INTEGER                                   :: timestep
    TYPE(t_operator_coeff)                    :: operators_coefficients
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp
    INTEGER :: return_status
    
    IF(discretization_scheme==MIMETIC_TYPE)THEN

      CALL solve_free_sfc_ab_mimetic( patch_3D, ocean_state, external_data, surface_fluxes, &
        & physics_parameters, timestep, operators_coefficients, solverCoeff_sp, return_status)

    ELSE
      CALL finish ('calc_vert_velocity: ',' Discretization type not supported !!')
    ENDIF

  END SUBROUTINE solve_free_surface_eq_ab
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Computation of new velocity in Adams-Bashforth timestepping.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE calc_normal_velocity_ab(patch_3D, ocean_state, operators_coefficients, &
    & solverCoeff_sp, external_data, physics_parameters)
    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_operator_coeff)               :: operators_coefficients
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp
    TYPE(t_external_data), TARGET        :: external_data
    TYPE (t_ho_params)                   :: physics_parameters
    !-----------------------------------------------------------------------
    IF(discretization_scheme==MIMETIC_TYPE)THEN

      CALL calc_normal_velocity_ab_mimetic(patch_3D, ocean_state, operators_coefficients)

    ELSE
      CALL finish ('calc_vert_velocity: ',' Discreization type not supported !!')
    ENDIF

  END SUBROUTINE calc_normal_velocity_ab
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Computation of new vertical velocity using continuity equation

  !! Calculate diagnostic vertical velocity from horizontal velocity using the
  !! incommpressibility condition in the continuity equation.
  !! For the case of the semi-implicit-AB scheme the land-sea-mask may be applied
  !! at least after collecting the whole explicit term.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn,   MPI-M (2006).
  !!
!<Optimize:inUse>
  SUBROUTINE calc_vert_velocity(patch_3D, ocean_state, operators_coefficients)
    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_hydro_ocean_state)            :: ocean_state
    TYPE(t_operator_coeff)               :: operators_coefficients
     !-----------------------------------------------------------------------

    !Store current vertical velocity before the new one is calculated
    ocean_state%p_diag%w_old = ocean_state%p_diag%w

    IF(discretization_scheme==MIMETIC_TYPE)THEN

      CALL calc_vert_velocity_mim_bottomup( patch_3D,     &
                                  & ocean_state,                   &
                                  & operators_coefficients )

    ELSE
      CALL finish ('calc_vert_velocity: ',' Discretization type not supported !!')
    ENDIF
  END SUBROUTINE calc_vert_velocity
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE update_time_indices(jg)
    INTEGER, INTENT(IN) :: jg
    INTEGER             :: n_temp
    ! Step 7: Swap time indices before output
    !         half time levels of semi-implicit Adams-Bashforth timestepping are
    !         stored in auxiliary arrays g_n and g_nimd of p_diag%aux
    n_temp    = nold(jg)
    nold(jg)  = nnew(jg)
    nnew(jg)  = n_temp
  END SUBROUTINE update_time_indices
END MODULE mo_ocean_ab_timestepping
