!>
!! Contains subroutines for aqua-planet simulation
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI-M (2010-09-21)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ape

  USE mo_kind,                 ONLY: wp
  USE mo_echam_convect_tables, ONLY: tlucua, jptlucu1, jptlucu2, &
    &                                lookuperror, lookupoverflow
  USE mo_physical_constants,   ONLY: vtmpc1

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: aqua_surface

CONTAINS
  !>
  !! Compute saturation specific humidity at the sea surface
  !! from the given SST and surface pressure.
  !!
  SUBROUTINE aqua_surface( kproma, kbdim, ppsfc, ptsfc, pqs )

    INTEGER, INTENT(IN)  :: kbdim, kproma
    REAL(wp),INTENT(IN)  :: ppsfc (kbdim)   !< surface pressure
    REAL(wp),INTENT(IN)  :: ptsfc (kbdim)   !< SST
    REAL(wp),INTENT(INOUT) :: pqs   (kbdim)   !< saturation specific humidity ! out

    INTEGER  :: itemp  !< temperature*1000
    INTEGER  :: jc     !< column index
    REAL(wp) :: zes    !< (saturation vapour pressure)*Rd/Rv/ps

    !-----
    lookupoverflow = .FALSE.

    DO jc = 1,kproma
      itemp   = NINT(ptsfc(jc)*1000._wp)
      IF (itemp<jptlucu1 .OR. itemp>jptlucu2) lookupoverflow = .TRUE.
      itemp   = MAX(MIN(itemp,jptlucu2),jptlucu1)
      zes     = tlucua(itemp)/ppsfc(jc)
      pqs(jc) = zes/(1._wp-vtmpc1*zes)
    ENDDO

    IF (lookupoverflow) CALL lookuperror ('aqua_surface')

  END SUBROUTINE aqua_surface
  !-------------

END MODULE mo_ape

