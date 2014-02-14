!>
!! @brief Testbed for ocean
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_ocean_testbed

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_datetime,            ONLY: t_datetime
  USE mo_oce_types,           ONLY: t_hydro_ocean_state
  USE mo_oce_physics,         ONLY: t_ho_params
  USE mo_sea_ice_types,       ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, t_sea_ice
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff

  USE mo_run_config,          ONLY: test_mode
  USE mo_grid_config,         ONLY: n_dom

  USE mo_ocean_testbed_dynamics,  ONLY: ocean_test_dynamics

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: ocean_testbed

CONTAINS
  

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_testbed( namelist_filename, shr_namelist_filename, &
    & patch_3d, ocean_state, p_ext_data,          &
    & datetime, p_sfc_flx, p_phys_param,             &
    & p_as, p_atm_f, p_ice, operators_coefficients)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: p_ext_data(n_dom)
    TYPE(t_datetime), INTENT(inout)                  :: datetime
    TYPE(t_sfc_flx)                                  :: p_sfc_flx
    TYPE (t_ho_params)                               :: p_phys_param
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: p_atm_f
    TYPE (t_sea_ice),         INTENT(inout)          :: p_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients

    CHARACTER(LEN=*), PARAMETER ::  method_name = "ocean_testbed"

    SELECT CASE (test_mode)
      CASE (10)
        CALL ocean_test_dynamics( patch_3d, ocean_state, p_ext_data,   &
          & datetime, p_sfc_flx, p_phys_param,             &
          & p_as, p_atm_f, p_ice,operators_coefficients)

      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT


  END SUBROUTINE ocean_testbed
  !-------------------------------------------------------------------------



END MODULE mo_ocean_testbed

