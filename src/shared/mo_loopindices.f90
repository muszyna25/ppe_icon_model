!>
!! This module contains subroutines needed to determine the start and end.
!!
!! This module contains subroutines needed to determine the start and end
!! indices of do loops for a given patch and block index.
!!
!! @par Revision History
!!  Created by Guenther Zaengl, DWD (2009-03-21)
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
!!
MODULE mo_loopindices
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!

USE mo_kind,            ONLY: wp
USE mo_model_domain,    ONLY: t_patch
USE mo_impl_constants,  ONLY: min_rlcell, min_rledge, min_rlvert
USE mo_parallel_configuration,  ONLY: nproma

IMPLICIT NONE

CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

PUBLIC

CONTAINS


!-------------------------------------------------------------------------
!
!
!
!>
!! Computes the start and end indices of do loops for cell-based variables.
!!
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2009-02-24
!!
SUBROUTINE get_indices_c(p_patch, i_blk, i_startblk, i_endblk, i_startidx, &
                         i_endidx, irl_start, opt_rl_end, opt_chdom)


TYPE(t_patch), INTENT(IN) :: p_patch
INTEGER, INTENT(IN) :: i_blk      ! Current block (variable jb in do loops)
INTEGER, INTENT(IN) :: i_startblk ! Start block of do loop
INTEGER, INTENT(IN) :: i_endblk   ! End block of do loop
INTEGER, INTENT(IN) :: irl_start  ! refin_ctrl level where do loop starts

INTEGER, OPTIONAL, INTENT(IN) :: opt_rl_end ! refin_ctrl level where do loop ends
INTEGER, OPTIONAL, INTENT(IN) :: opt_chdom ! child domain position ID where
                                           ! negative refin_ctrl indices refer to

INTEGER, INTENT(OUT) :: i_startidx, i_endidx ! Start and end indices (jc loop)

! Local variables

INTEGER :: irl_end, i_stdom, i_enddom

IF (PRESENT(opt_rl_end)) THEN
  irl_end = opt_rl_end
ELSE
  irl_end = min_rlcell
ENDIF

IF (PRESENT(opt_chdom)) THEN
  i_stdom = opt_chdom
  i_enddom = opt_chdom
ELSE
  i_stdom = 1
  i_enddom = MAX(1,p_patch%n_childdom)
ENDIF

IF (i_blk == i_startblk) THEN
  i_startidx = MAX(1,p_patch%cells%start_idx(irl_start,i_stdom))
  i_endidx   = nproma
  IF (i_blk == i_endblk) i_endidx = p_patch%cells%end_idx(irl_end,i_enddom)
ELSE IF (i_blk == i_endblk) THEN
  i_startidx = 1
  i_endidx   = p_patch%cells%end_idx(irl_end,i_enddom)
ELSE
  i_startidx = 1
  i_endidx = nproma
ENDIF

END SUBROUTINE get_indices_c

!-------------------------------------------------------------------------
!
!
!
!>
!! Computes the start and end indices of do loops for edge-based variables.
!!
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2009-02-24
!!
SUBROUTINE get_indices_e(p_patch, i_blk, i_startblk, i_endblk, i_startidx, &
                         i_endidx, irl_start, opt_rl_end, opt_chdom)


TYPE(t_patch), INTENT(IN) :: p_patch
INTEGER, INTENT(IN) :: i_blk      ! Current block (variable jb in do loops)
INTEGER, INTENT(IN) :: i_startblk ! Start block of do loop
INTEGER, INTENT(IN) :: i_endblk   ! End block of do loop
INTEGER, INTENT(IN) :: irl_start  ! refin_ctrl level where do loop starts

INTEGER, OPTIONAL, INTENT(IN) :: opt_rl_end ! refin_ctrl level where do loop ends
INTEGER, OPTIONAL, INTENT(IN) :: opt_chdom ! child domain position ID where
                                           ! negative refin_ctrl indices refer to

INTEGER, INTENT(OUT) :: i_startidx, i_endidx ! Start and end indices (je loop)


! Local variables

INTEGER :: irl_end, i_stdom, i_enddom

IF (PRESENT(opt_rl_end)) THEN
  irl_end = opt_rl_end
ELSE
  irl_end = min_rledge
ENDIF

IF (PRESENT(opt_chdom)) THEN
  i_stdom = opt_chdom
  i_enddom = opt_chdom
ELSE
  i_stdom = 1
  i_enddom = MAX(1,p_patch%n_childdom)
ENDIF


IF (i_blk == i_startblk) THEN
  i_startidx = MAX(1,p_patch%edges%start_idx(irl_start,i_stdom))
  i_endidx   = nproma
  IF (i_blk == i_endblk) i_endidx = p_patch%edges%end_idx(irl_end,i_enddom)
ELSE IF (i_blk == i_endblk) THEN
  i_startidx = 1
  i_endidx   = p_patch%edges%end_idx(irl_end,i_enddom)
ELSE
  i_startidx = 1
  i_endidx = nproma
ENDIF

END SUBROUTINE get_indices_e


!-------------------------------------------------------------------------
!
!
!
!>
!! Computes the start and end indices of do loops for cell-based variables.
!!
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2009-02-24
!!
SUBROUTINE get_indices_v(p_patch, i_blk, i_startblk, i_endblk, i_startidx, &
                         i_endidx, irl_start, opt_rl_end, opt_chdom)


TYPE(t_patch), INTENT(IN) :: p_patch
INTEGER, INTENT(IN) :: i_blk      ! Current block (variable jb in do loops)
INTEGER, INTENT(IN) :: i_startblk ! Start block of do loop
INTEGER, INTENT(IN) :: i_endblk   ! End block of do loop
INTEGER, INTENT(IN) :: irl_start  ! refin_ctrl level where do loop starts

INTEGER, OPTIONAL, INTENT(IN) :: opt_rl_end ! refin_ctrl level where do loop ends
INTEGER, OPTIONAL, INTENT(IN) :: opt_chdom ! child domain position ID where
                                           ! negative refin_ctrl indices refer to

INTEGER, INTENT(OUT) :: i_startidx, i_endidx ! Start and end indices (jv loop)

! Local variables

INTEGER :: irl_end, i_stdom, i_enddom

IF (PRESENT(opt_rl_end)) THEN
  irl_end = opt_rl_end
ELSE
  irl_end = min_rlvert
ENDIF

IF (PRESENT(opt_chdom)) THEN
  i_stdom = opt_chdom
  i_enddom = opt_chdom
ELSE
  i_stdom = 1
  i_enddom = MAX(1,p_patch%n_childdom)
ENDIF

IF (i_blk == i_startblk) THEN
  i_startidx = p_patch%verts%start_idx(irl_start,i_stdom)
  i_endidx   = nproma
  IF (i_blk == i_endblk) i_endidx = p_patch%verts%end_idx(irl_end,i_enddom)
ELSE IF (i_blk == i_endblk) THEN
  i_startidx = 1
  i_endidx   = p_patch%verts%end_idx(irl_end,i_enddom)
ELSE
  i_startidx = 1
  i_endidx = nproma
ENDIF

END SUBROUTINE get_indices_v

END MODULE mo_loopindices

