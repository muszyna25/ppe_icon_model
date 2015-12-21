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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ocean_testbed

  USE mo_exception,           ONLY: finish
  USE mo_model_domain,        ONLY: t_patch_3d
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state, t_solverCoeff_singlePrecision, t_operator_coeff
  USE mo_ocean_physics,       ONLY: t_ho_params
  USE mo_sea_ice_types,       ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, t_sea_ice
  USE mo_ocean_surface_types, ONLY: t_ocean_surface

  USE mo_run_config,          ONLY: test_mode
  USE mo_grid_config,         ONLY: n_dom

  USE mo_ocean_testbed_modules,     ONLY: ocean_test_modules
  USE mo_testbed_ocean_performance, ONLY: ocean_test_performance
  USE mo_ocean_testbed_operators,   ONLY: ocean_test_operators
  USE mo_ocean_testbed_read,        ONLY: ocean_test_read
  USE mo_ocean_testbed_quads,       ONLY: ocean_test_quads
  USE mo_ocean_math_operators,      ONLY: update_height_depdendent_variables
  USE mtime,                        ONLY: datetime
  
!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: ocean_testbed

CONTAINS
  

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_testbed( namelist_filename, shr_namelist_filename, &
    & patch_3d, ocean_state, external_data,          &
    & this_datetime, surface_fluxes, ocean_surface, physics_parameters,             &
    & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice, operators_coefficients, &
    & solverCoeff_sp)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: external_data(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE (t_ocean_surface)                           :: ocean_surface
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: oceans_atmosphere
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: oceans_atmosphere_fluxes
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp

    CHARACTER(LEN=*), PARAMETER ::  method_name = "ocean_testbed"

    CALL update_height_depdendent_variables( patch_3D, ocean_state(1), external_data(1), operators_coefficients, solverCoeff_sp)
	
    SELECT CASE (test_mode)
      CASE (1 : 99)  !  1 - 99 test ocean modules
        CALL ocean_test_modules( patch_3d, ocean_state,  &
          & this_datetime, surface_fluxes, ocean_surface, physics_parameters,             &
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

      CASE (1101) ! 1000 - 1100 performance tests
        CALL ocean_test_read( namelist_filename, shr_namelist_filename, &
          & patch_3d)

      CASE (1102) ! 1101 -  other tests
        CALL ocean_test_quads( namelist_filename, shr_namelist_filename, &
          & patch_3d)

      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT


  END SUBROUTINE ocean_testbed
  !-------------------------------------------------------------------------



END MODULE mo_ocean_testbed

