!>
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, 2010-08-16
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_echam_phy_cleanup

  USE mo_echam_phy_config,   ONLY: lconv  => echam_phy_config%lconv, &
                                 & lvdiff => echam_phy_config%lvdiff
  USE mo_echam_conv_nml,     ONLY: cleanup_cuparam
  USE mo_vdiff_solver,       ONLY: cleanup_vdiff_solver

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: cleanup_echam_phy

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

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

    IF (lvdiff) CALL cleanup_vdiff_solver  ! Deallocate array "matrix_idx"

    IF (lconv) THEN
      CALL cleanup_cuparam  ! deallocate array "cevapcu"
    END IF

   !"CALL destruct_echam_phy_state" can be called here instead of in "control_model"

  END SUBROUTINE cleanup_echam_phy
  !-------------

END MODULE mo_echam_phy_cleanup

