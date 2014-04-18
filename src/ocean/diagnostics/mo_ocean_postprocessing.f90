!>
!! @brief postprocessing tools
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
!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_postprocessing

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_io_units,            ONLY: filename_max
  USE mo_run_config,          ONLY: test_mode

  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_grid_config,         ONLY: n_dom
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_datetime,            ONLY: t_datetime
  USE mo_oce_types,           ONLY: t_hydro_ocean_state, t_solverCoeff_singlePrecision, t_operator_coeff
  USE mo_oce_physics,         ONLY: t_ho_params
  USE mo_sea_ice_types,       ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, t_sea_ice
  USE mo_ocean_nml
  USE mo_oce_thermodyn,       ONLY: calc_density, calc_internal_press
  USE mo_dynamics_config,     ONLY: nold, nnew
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

  USE mo_oce_math_operators,  ONLY: calculate_thickness

  USE mo_mpi,                 ONLY: work_mpi_barrier, my_process_is_stdio
  USE mo_timer,               ONLY: init_timer, ltimer, new_timer, timer_start, timer_stop, &
    & print_timer, activate_sync_timers, timers_level, timer_barrier, timer_radiaton_recv
  USE mo_parallel_config,     ONLY: nproma
  USE mo_master_control,      ONLY: get_my_process_name, get_my_model_no

  USE mo_read_interface
  USE mo_statistics
!  USE mo_netcdf_read
!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: ocean_postprocess


CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_postprocess( namelist_filename, shr_namelist_filename, &
    & patch_3d, ocean_state, external_data,          &
    & datetime, surface_fluxes, physics_parameters,             &
    & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice, operators_coefficients, &
    & solverCoeff_sp)

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
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp

    CHARACTER(LEN=*), PARAMETER ::  method_name = "ocean_postprocess"

    SELECT CASE (test_mode)
      CASE (2001)  !  1 - 99 test ocean modules
        CALL  ocean_postprocess_conservation(patch_3d=patch_3d, ocean_state=ocean_state(1), &
          & external_data=external_data(1), operators_coefficients=operators_coefficients, solverCoeff_sp=solverCoeff_sp, &
          & fileName="firstFile.nc", timeIndex=1)
        CALL  ocean_postprocess_conservation(patch_3d=patch_3d, ocean_state=ocean_state(1), &
          & external_data=external_data(1), operators_coefficients=operators_coefficients, solverCoeff_sp=solverCoeff_sp, &
          & fileName="secondFile.nc", timeIndex=4)

      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT

  END SUBROUTINE ocean_postprocess
  !-------------------------------------------------------------------------




  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_postprocess_conservation(patch_3d, ocean_state, &
          & external_data, operators_coefficients, solverCoeff_sp, &
          & fileName, timeIndex)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    TYPE(t_external_data), TARGET, INTENT(in)        :: external_data
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp
    CHARACTER(LEN=*), INTENT(in)         :: fileName
    INTEGER, INTENT(in)                  :: timeIndex

    REAL(wp), POINTER :: T(:,:,:,:), S(:,:,:,:), h(:,:,:)
    REAL(wp), POINTER :: tracers(:,:,:,:)
    REAL(wp) :: meanRho, meanH

    INTEGER :: return_status, stream_id
!    INTEGER :: idx, level, block
!    INTEGER :: start_index, end_index
!    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch),POINTER            :: patch_2d

    CHARACTER(*), PARAMETER :: method_name = "mo_ocean_postprocessing:ocean_postprocess_conservation"

    ! CALL message(method_name,   initialState_InputFileName)
    patch_2d => patch_3d%p_patch_2d(1)
    tracers => ocean_state%p_prog(1)%tracer(:,:,:,:)
    nold(1) = 1
    nnew(1) = 1
    !---------------------------------------------------------------------
    CALL read_onCells_3D_time(                &
      & filename=fileName,                    &
      & variable_name="t_acc",                &
      & start_timestep=timeIndex,             &
      & end_timestep=timeIndex,               &
      & return_pointer=T,                     &
      & patch=patch_2d,                       &
      & return_status=return_status )

    tracers(:,:,:,1) = T(:,:,:,1)

    CALL read_onCells_3D_time(                &
      & filename=fileName,                    &
      & variable_name="s_acc",                &
      & start_timestep=timeIndex,             &
      & end_timestep=timeIndex,               &
      & return_pointer=S,                     &
      & patch=patch_2d,                       &
      & return_status=return_status )

    tracers(:,:,:,2) = S(:,:,:,1)

    CALL read_onCells_2D_time(                &
      & filename=fileName,                    &
      & variable_name="h_acc",                &
      & start_timestep=timeIndex,             &
      & end_timestep=timeIndex,               &
      & return_pointer=h,                     &
      & patch=patch_2d,                       &
      & return_status=return_status )

    ocean_state%p_prog(1)%h(:,:) = h(:,:,1)

    DEALLOCATE(T, S, h)

    !---------------------------------------------------------------------
    ! calclulate volume
    CALL calculate_thickness( patch_3D, ocean_state, external_data, operators_coefficients, solverCoeff_sp)

    CALL calc_density( patch_3d,            &
        & tracers(:,:,:,1:2),&
        & ocean_state%p_diag%rho(:,:,:) )

    meanRho = total_mean(values=ocean_state%p_diag%rho(:,:,:), weights=patch_3d%p_patch_1d(1)%prism_volume(:,:,:), &
      & in_subset=patch_3d%p_patch_2d(1)%cells%owned)

    CALL levels_horizontal_mean(values=ocean_state%p_prog(1)%h(:,:), weights=patch_2D%cells%area(:,:), &
      & in_subset=patch_3d%p_patch_2d(1)%cells%owned, mean=meanH)

    IF (my_process_is_stdio()) &
      write(0,*) TRIM(fileName), " Mean height=", meanH, " Mean density=", meanRho

  END SUBROUTINE ocean_postprocess_conservation
  !-------------------------------------------------------------------------


END MODULE mo_ocean_postprocessing

