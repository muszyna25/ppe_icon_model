!>
!! @brief Testbed for ocean
!!
!! The test_mode variable of the run_nml defines the mode of running the ocean.
!! When test_mode==0 the normal timesteping (perform_ho_stepping) is called,
!! else the ocean_testbed is called. There are currently 3 testbed modules:
!!
!! test_mode=1-99: ocean_testbed_modules; test ocean modules like advection, diffusion, etc
!!
!! test_mode=99-999: ocean_test_operators; test consistency, properties, round-off errors of operators
!!
!! test_mode=1000-1100: ocean_test_performance; reserved for performance tests of individual routines
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

  USE mo_ocean_testbed_modules,     ONLY: ocean_test_advection
  USE mo_testbed_ocean_performance, ONLY: ocean_test_performance
  USE mo_ocean_testbed_operators,   ONLY: ocean_test_operators

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: ocean_testbed

CONTAINS
  

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_testbed( namelist_filename, shr_namelist_filename, &
    & patch_3d, ocean_state, external_data,          &
    & datetime, surface_fluxes, physics_parameters,             &
    & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice, operators_coefficients)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: external_data(n_dom)
    TYPE(t_datetime), INTENT(inout)                  :: datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: oceans_atmosphere
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: oceans_atmosphere_fluxes
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients

    CHARACTER(LEN=*), PARAMETER ::  method_name = "ocean_testbed"

    SELECT CASE (test_mode)
      CASE (1 : 99)  !  1 - 99 test ocean modules
        CALL ocean_test_advection( patch_3d, ocean_state, external_data,   &
          & datetime, surface_fluxes, physics_parameters,             &
          & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice,operators_coefficients)

      CASE (100 : 999) ! 100 - 999 test ocean operators
        CALL ocean_test_operators( namelist_filename, shr_namelist_filename, &
          & patch_3d, ocean_state, external_data,   &
          & surface_fluxes, physics_parameters,             &
          & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice,operators_coefficients)

      CASE (1000 : 1100) ! 1000 - 1100 performance tests
        CALL ocean_test_performance( namelist_filename, shr_namelist_filename, &
          & patch_3d, ocean_state, external_data,   &
          & surface_fluxes, physics_parameters,             &
          & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice,operators_coefficients)

      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT


  END SUBROUTINE ocean_testbed
  !-------------------------------------------------------------------------



END MODULE mo_ocean_testbed

