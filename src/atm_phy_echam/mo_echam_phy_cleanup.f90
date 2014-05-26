!>
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, 2010-08-16
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_phy_cleanup

  USE mo_echam_phy_memory,   ONLY: destruct_echam_phy_state
  USE mo_echam_phy_config,   ONLY: phy_config => echam_phy_config
  USE mo_echam_conv_config,  ONLY: cleanup_echam_convection
  USE mo_vdiff_solver,       ONLY: cleanup_vdiff_solver

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: cleanup_echam_phy

CONTAINS
  !>
  !! Top-level routine for the cleanup for ECHAM6 physics.
  !! It calls a series of subroutines to deallocate parameter arrays with
  !! "allocatable" attribute.
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2010-08)
  !!
  SUBROUTINE cleanup_echam_phy

    IF (phy_config%lvdiff) CALL cleanup_vdiff_solver  ! Deallocate array "matrix_idx"

    IF (phy_config%lconv) THEN
      CALL cleanup_echam_convection  ! deallocate array "cevapcu"
    END IF

    CALL destruct_echam_phy_state

  END SUBROUTINE cleanup_echam_phy
  !-------------

END MODULE mo_echam_phy_cleanup

