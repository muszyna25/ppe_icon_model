!>
!! Contains subroutines for aqua-planet simulation
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI-M (2010-09-21)
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
MODULE mo_ape

  USE mo_kind,               ONLY: wp
  USE mo_convect_tables,     ONLY: tlucua, jptlucu1, jptlucu2, &
                                   lookuperror, lookupoverflow
#ifdef __ICON__
  USE mo_physical_constants, ONLY: vtmpc1
#else
  USE mo_constants,          ONLY: vtmpc1
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: aqua_surface

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !! Compute saturation specific humidity at the sea surface
  !! from the given SST and surface pressure.
  !!
  SUBROUTINE aqua_surface( kproma, kbdim, ppsfc, ptsfc, pqs )

    INTEGER, INTENT(IN)  :: kbdim, kproma
    REAL(wp),INTENT(IN)  :: ppsfc (kbdim)   !< surface pressure
    REAL(wp),INTENT(IN)  :: ptsfc (kbdim)   !< SST
    REAL(wp),INTENT(OUT) :: pqs   (kbdim)   !< saturation specific humidity

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

