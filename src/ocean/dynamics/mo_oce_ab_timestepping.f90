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
MODULE mo_oce_ab_timestepping
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_ocean_nml,                      ONLY: discretization_scheme
USE mo_dynamics_config,                ONLY: nold, nnew
USE mo_sea_ice_types,                  ONLY: t_sfc_flx
USE mo_model_domain,                   ONLY: t_patch_3D !, t_patch
USE mo_ext_data_types,                 ONLY: t_external_data
USE mo_oce_ab_timestepping_mimetic,    ONLY: solve_free_sfc_ab_mimetic,       &
  &                                          calc_normal_velocity_ab_mimetic, &
  &                                          calc_vert_velocity_mim_bottomup
USE mo_oce_physics,                    ONLY: t_ho_params
USE mo_oce_types,                      ONLY: t_hydro_ocean_state, t_operator_coeff, t_solverCoeff_singlePrecision
USE mo_exception,                      ONLY: finish!, message_text
IMPLICIT NONE

PRIVATE

INTEGER, PARAMETER :: MIMETIC_TYPE = 1
INTEGER, PARAMETER :: RBF_TYPE     = 2
!
! PUBLIC INTERFACE
!
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
!<Optimize_Used>
  SUBROUTINE solve_free_surface_eq_ab(p_patch_3D, p_os, p_ext_data, p_sfc_flx, &
    &                                 p_phys_param, timestep, p_op_coeff, solverCoeff_sp)
    TYPE(t_patch_3D ),TARGET, INTENT(INOUT)   :: p_patch_3D
    TYPE(t_hydro_ocean_state), TARGET             :: p_os
    TYPE(t_external_data), TARGET                 :: p_ext_data
    TYPE(t_sfc_flx), INTENT(INOUT)                :: p_sfc_flx
    TYPE (t_ho_params)                            :: p_phys_param
    INTEGER                                       :: timestep
    TYPE(t_operator_coeff)                        :: p_op_coeff
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp

    IF(discretization_scheme==MIMETIC_TYPE)THEN

      CALL solve_free_sfc_ab_mimetic( p_patch_3D, p_os, p_ext_data, p_sfc_flx, &
        &                            p_phys_param, timestep, p_op_coeff, solverCoeff_sp)

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
!<Optimize_Used>
  SUBROUTINE calc_normal_velocity_ab(p_patch_3D, p_os, p_op_coeff, solverCoeff_sp, p_ext_data, p_phys_param)
    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: p_patch_3D
    TYPE(t_hydro_ocean_state), TARGET    :: p_os
    TYPE(t_operator_coeff)               :: p_op_coeff
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp
    TYPE(t_external_data), TARGET        :: p_ext_data
    TYPE (t_ho_params)                   :: p_phys_param
    !-----------------------------------------------------------------------
    IF(discretization_scheme==MIMETIC_TYPE)THEN

      CALL calc_normal_velocity_ab_mimetic(p_patch_3D, p_os, p_op_coeff, solverCoeff_sp, p_ext_data)

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
!<Optimize_Used>
  SUBROUTINE calc_vert_velocity(p_patch_3D, p_os, p_op_coeff)
    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: p_patch_3D
    TYPE(t_hydro_ocean_state)            :: p_os
    TYPE(t_operator_coeff)               :: p_op_coeff
     !-----------------------------------------------------------------------

    !Store current vertical velocity before the new one is calculated
    p_os%p_diag%w_old = p_os%p_diag%w

    IF(discretization_scheme==MIMETIC_TYPE)THEN

      CALL calc_vert_velocity_mim_bottomup( p_patch_3D,     &
                                  & p_os,                   &
                                  & p_op_coeff )

    ELSE
      CALL finish ('calc_vert_velocity: ',' Discretization type not supported !!')
    ENDIF
  END SUBROUTINE calc_vert_velocity
  !-------------------------------------------------------------------------
!<Optimize_Used>
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
END MODULE mo_oce_ab_timestepping
