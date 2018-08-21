!>
!! Provide an implementation of the ocean tracer nuding functionality.
!!
!!
!! @author Helmuth Haak, MPI
!!
!! @par Revision History
!!  Original version by Helmuth Haak, MPI-M (2018)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!----------------------------
#include "iconfor_dsl_definitions.inc"
#include "omp_definitions.inc"
#include "icon_definitions.inc"
!----------------------------
MODULE mo_ocean_nudging_types
!-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  !USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  !USE mo_impl_constants,            ONLY: sea_boundary, sea, min_dolic
  !USE mo_math_constants,            ONLY: pi
  !USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer,              &
  !  & threshold_min_t, threshold_max_t, threshold_min_s, threshold_max_s, &
  !  & type_3dimrelax_temp, para_3dimrelax_temp,                           &
  !  & type_3dimrelax_salt, para_3dimrelax_salt
  !USE mo_util_dbg_prnt,             ONLY: dbg_print
  !USE mo_parallel_config,           ONLY: nproma
  !USE mo_dynamics_config,           ONLY: nold, nnew
  !USE mo_run_config,                ONLY: dtime, ltimer, debug_check_level
!!  USE mo_ocean_types,               ONLY: t_hydro_ocean_state, t_ocean_tracer
!  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
!  USE mo_exception,                 ONLY: finish !, message_text, message
!  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
!  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range


  IMPLICIT NONE

  TYPE t_ocean_nudge

    ! Variables for 3-dim tracer relaxation:
    onCells :: &
      & data_3dimRelax_Temp, & ! 3-dim temperature relaxation data (T*)
      & forc_3dimRelax_Temp, & ! 3-dim temperature relaxation forcing (1/tau*(T-T*))
      & data_3dimRelax_Salt, & ! 3-dim salinity relaxation data (T*)
      & forc_3dimRelax_Salt, &    ! 3-dim salinity relaxation forcing (1/tau*(T-T*))
      & relax_3dim_coefficient ! 3-dim relaxation coefficient when the relaxation varies

  END TYPE t_ocean_nudge

  PRIVATE

  PUBLIC :: t_ocean_nudge

END MODULE mo_ocean_nudging_types


