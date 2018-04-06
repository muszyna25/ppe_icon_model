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

  USE mtime,                 ONLY: OPERATOR(>)
  USE mo_grid_config,        ONLY: n_dom
  USE mo_echam_phy_memory,   ONLY: destruct_echam_phy_state
  USE mo_psrad_memory,       ONLY: destruct_psrad_forcing_list
  USE mo_echam_phy_config,   ONLY: echam_phy_tc, dt_zero
  USE mo_echam_cnv_config,   ONLY: dealloc_echam_cnv_config
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

    LOGICAL :: lany
    INTEGER :: jg

    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_vdf > dt_zero)
    END DO
    IF (lany) CALL cleanup_vdiff_solver      ! deallocate array "matrix_idx"

    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_cnv > dt_zero)
    END DO
    IF (lany) CALL dealloc_echam_cnv_config  ! deallocate array "cevapcu"

    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_rad > dt_zero)
    END DO
    IF (lany) CALL destruct_psrad_forcing_list
   
    CALL destruct_echam_phy_state

  END SUBROUTINE cleanup_echam_phy
  !-------------

END MODULE mo_echam_phy_cleanup

