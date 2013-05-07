!> This module defines small arithmetic operations ("post-ops") as
!! post-processing tasks.
!!
!! These post-processing tasks are restricted to point-wise
!! operations (no halo synchronization) of a single field, like
!! value scaling.
!!
!! @note The "post-ops" are performed at output time and DO NOT
!!       MODIFY THE FIELD ITSELF.
!! 
!!
!! @author F. Prill, DWD
!!
!! @par Revision History
!! Initial implementation,            F. Prill, DWD (2013-03-20)
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

MODULE mo_post_op

  USE mo_kind,                ONLY: wp
  USE mo_var_metadata,        ONLY: t_post_op_meta, POST_OP_NONE,    &
    &                               POST_OP_SCALE
  USE mo_exception,           ONLY: finish
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: perform_post_op

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_post_op')

CONTAINS

  !> Performs small arithmetic operations ("post-ops") as
  !  post-processing tasks.
  SUBROUTINE perform_post_op(post_op, field3D)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    REAL(wp),              INTENT(INOUT) :: field3D(:,:,:)
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op"
    INTEGER :: idim(3), l1,l2,l3

    idim = SHAPE(field3D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2,l3), SCHEDULE(runtime)
      DO l3=1,idim(3)
        DO l2=1,idim(2)
          DO l1=1,idim(1)
            field3D(l1,l2,l3) = field3D(l1,l2,l3) * post_op%arg1
          END DO
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
      !
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

  END SUBROUTINE perform_post_op


END MODULE mo_post_op
