!>
!
!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_ocean_testbed_solverMatrix
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp, sp
  USE mo_impl_constants,         ONLY: max_char_length, beta_plane_coriolis,full_coriolis, &
    &                               SEA_BOUNDARY, BOUNDARY, SEA, min_dolic
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d, t_subset_range, t_patch_vert
  USE mo_grid_subset,            ONLY: get_index_range
  USE mo_sync,                   ONLY: sync_patch_array, sync_e, sync_c !, sync_v
  USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_run_config,             ONLY: nsteps, dtime, ltimer, output_mode, test_mode
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ocean_ab_timestepping,  ONLY: solve_free_surface_eq_ab, &
    & calc_normal_velocity_ab,  &
    & calc_vert_velocity,       &
    & update_time_indices
  USE mo_ocean_types,            ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, &
    & t_hydro_ocean_prog
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff! , update_diffusion_matrices
  USE mo_mpi,                    ONLY: my_process_is_stdio, num_work_procs
  USE mo_time_config,            ONLY: time_config
  USE mo_util_dbg_prnt,          ONLY: dbg_print
  USE mo_parallel_config,        ONLY: nproma
  USE mo_timer
  USE mo_ocean_ab_timestepping_mimetic, ONLY : lhs_surface_height

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: createSolverMatrix

  !-------------------------------------------------------------------------
CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE createSolverMatrix( patch_3d, ocean_state, operators_coefficients)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(1)
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients

    REAL(wp) :: x(1:nproma,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: lhs(1:nproma,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    INTEGER :: column,row

    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D

    INTEGER, PARAMETER :: fileNo = 501
    CHARACTER(len=*), PARAMETER :: fileName = 'ocean_matrix.txt'

    CHARACTER(len=*), PARAMETER :: method_name = 'createSolverMatrix'


    IF (nproma /= 1) &
      CALL finish(method_name, "nproma should be = 1")

    IF (num_work_procs /= 1) &
      CALL finish(method_name, "mpi processes should be = 1")

    patch_2D           => patch_3d%p_patch_2d(1)
    cells_in_domain    => patch_2D%cells%in_domain

    open (fileNo, FILE=fileName, STATUS='new')

    DO column = cells_in_domain%start_block, cells_in_domain%end_block
      x   = 0.0_wp
      lhs = 0.0_wp
      x(1,column) = 1.0_wp

      lhs = lhs_surface_height( x, patch_3d, patch_3D%column_thick_e,&
        & operators_coefficients)

      DO row = cells_in_domain%start_block, cells_in_domain%end_block
        IF (lhs(1, row) /= 0.0_wp) &
          write(fileNo, *) "(", row, ",", column, ")", lhs(1, row)
      ENDDO

    ENDDO

    CLOSE(fileNo)

  END SUBROUTINE createSolverMatrix
  !-------------------------------------------------------------------------


END MODULE mo_ocean_testbed_solverMatrix
