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

  INTERFACE perform_post_op
    MODULE PROCEDURE perform_post_op_2D
    MODULE PROCEDURE perform_post_op_3D
  END INTERFACE perform_post_op

CONTAINS

  !> Performs small arithmetic operations ("post-ops") on 3D field as
  !  post-processing tasks.
  SUBROUTINE perform_post_op_2D(post_op, field2D, opt_inverse)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    REAL(wp),              INTENT(INOUT) :: field2D(:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_2D"
    REAL(wp) :: scalfac           ! scale factor
    INTEGER  :: idim(2), l1,l2


    idim = SHAPE(field2D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          scalfac = 1._wp/post_op%arg1
        ELSE
          scalfac = post_op%arg1
        ENDIF
      ELSE
        scalfac = post_op%arg1
      ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2), SCHEDULE(runtime)
      DO l2=1,idim(2)
        DO l1=1,idim(1)
          field2D(l1,l2) = field2D(l1,l2) * scalfac
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
      !
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

  END SUBROUTINE perform_post_op_2D


  !> Performs small arithmetic operations ("post-ops") on 3D field as
  !  post-processing tasks.
  SUBROUTINE perform_post_op_3D(post_op, field3D, opt_inverse)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    REAL(wp),              INTENT(INOUT) :: field3D(:,:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_3D"
    REAL(wp) :: scalfac           ! scale factor
    INTEGER  :: idim(3), l1,l2,l3


    idim = SHAPE(field3D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          scalfac = 1._wp/post_op%arg1
        ELSE
          scalfac = post_op%arg1
        ENDIF
      ELSE
        scalfac = post_op%arg1
      ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2,l3), SCHEDULE(runtime)
      DO l3=1,idim(3)
        DO l2=1,idim(2)
          DO l1=1,idim(1)
            field3D(l1,l2,l3) = field3D(l1,l2,l3) * scalfac
          END DO
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
      !
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

  END SUBROUTINE perform_post_op_3D


END MODULE mo_post_op
