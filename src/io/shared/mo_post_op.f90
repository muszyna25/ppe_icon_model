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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_post_op

  USE mo_kind,                ONLY: dp, sp
  USE mo_var_metadata_types,  ONLY: t_post_op_meta, POST_OP_NONE,    &
    &                               POST_OP_SCALE, POST_OP_LUC
#ifndef __NO_ICON_ATMO__
  USE mo_lnd_nwp_config,      ONLY: convert_luc_ICON2GRIB
#endif
  USE mo_exception,           ONLY: finish
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: perform_post_op

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_post_op')

  INTERFACE perform_post_op
    MODULE PROCEDURE perform_post_op_r2D
    MODULE PROCEDURE perform_post_op_s2D
    MODULE PROCEDURE perform_post_op_i2D
    MODULE PROCEDURE perform_post_op_r3D
    MODULE PROCEDURE perform_post_op_s3D
    MODULE PROCEDURE perform_post_op_i3D
  END INTERFACE perform_post_op

CONTAINS

  !> Performs small arithmetic operations ("post-ops") on 2D REAL field 
  !  as post-processing tasks.
  SUBROUTINE perform_post_op_r2D(post_op, field2D, opt_inverse)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    REAL(dp),              INTENT(INOUT) :: field2D(:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_r2D"
    REAL(dp) :: scalfac           ! scale factor
    INTEGER  :: idim(2), l1,l2


    idim = SHAPE(field2D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          scalfac = 1._dp/post_op%arg1%rval
        ELSE
          scalfac = post_op%arg1%rval
        ENDIF
      ELSE
        scalfac = post_op%arg1%rval
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

  END SUBROUTINE perform_post_op_r2D


  !> Performs small arithmetic operations ("post-ops") on 2D REAL field 
  !  as post-processing tasks.
  SUBROUTINE perform_post_op_s2D(post_op, field2D, opt_inverse)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    REAL(sp),              INTENT(INOUT) :: field2D(:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_s2D"
    REAL(sp) :: scalfac           ! scale factor
    INTEGER  :: idim(2), l1,l2


    idim = SHAPE(field2D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          scalfac = 1._sp/post_op%arg1%rval
        ELSE
          scalfac = post_op%arg1%rval
        ENDIF
      ELSE
        scalfac = post_op%arg1%rval
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

  END SUBROUTINE perform_post_op_s2D


  !> Performs small arithmetic operations ("post-ops") on 2D INTEGER field 
  !  as post-processing tasks.
  SUBROUTINE perform_post_op_i2D(post_op, field2D, opt_inverse)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    INTEGER              , INTENT(INOUT) :: field2D(:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_i2D"
    INTEGER  :: scalfac           ! scale factor
    INTEGER  :: idim(2), l1,l2


    idim = SHAPE(field2D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          CALL finish(routine, "Inverse POST_OP_SCALE not applicable to INTEGER field!")
        ELSE
          scalfac = post_op%arg1%ival
        ENDIF
      ELSE
        scalfac = post_op%arg1%ival
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
#ifndef __NO_ICON_ATMO__
    CASE(POST_OP_LUC)
!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2), SCHEDULE(runtime)
      DO l2=1,idim(2)
        DO l1=1,idim(1)
          field2D(l1,l2) = convert_luc_ICON2GRIB(lc_datbase  = post_op%arg1%ival, &
            &                                    iluc_in     = field2D(l1,l2),    &
            &                                    opt_linverse= opt_inverse)
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
#endif
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

  END SUBROUTINE perform_post_op_i2D


  !> Performs small arithmetic operations ("post-ops") on 3D REAL field 
  !  as post-processing tasks.
  SUBROUTINE perform_post_op_r3D(post_op, field3D, opt_inverse)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    REAL(dp),              INTENT(INOUT) :: field3D(:,:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_r3D"
    REAL(dp) :: scalfac           ! scale factor
    INTEGER  :: idim(3), l1,l2,l3


    idim = SHAPE(field3D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          scalfac = 1._dp/post_op%arg1%rval
        ELSE
          scalfac = post_op%arg1%rval
        ENDIF
      ELSE
        scalfac = post_op%arg1%rval
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

  END SUBROUTINE perform_post_op_r3D


  !> Performs small arithmetic operations ("post-ops") on 3D REAL field 
  !  as post-processing tasks.
  SUBROUTINE perform_post_op_s3D(post_op, field3D, opt_inverse)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    REAL(sp),              INTENT(INOUT) :: field3D(:,:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_s3D"
    REAL(sp) :: scalfac           ! scale factor
    INTEGER  :: idim(3), l1,l2,l3


    idim = SHAPE(field3D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          scalfac = 1._sp/post_op%arg1%rval
        ELSE
          scalfac = post_op%arg1%rval
        ENDIF
      ELSE
        scalfac = post_op%arg1%rval
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

  END SUBROUTINE perform_post_op_s3D


  !> Performs small arithmetic operations ("post-ops") on 3D INTEGER field 
  !  as post-processing tasks.
  SUBROUTINE perform_post_op_i3D(post_op, field3D, opt_inverse)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    INTEGER ,              INTENT(INOUT) :: field3D(:,:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_i3D"
    REAL(dp) :: scalfac           ! scale factor
    INTEGER  :: idim(3), l1,l2,l3


    idim = SHAPE(field3D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          CALL finish(routine, "Inverse POST_OP_SCALE not applicable to INTEGER field!")
        ELSE
          scalfac = post_op%arg1%ival
        ENDIF
      ELSE
        scalfac = post_op%arg1%ival
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
#ifndef __NO_ICON_ATMO__
    CASE(POST_OP_LUC)
!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2,l3), SCHEDULE(runtime)
      DO l3=1,idim(3)
        DO l2=1,idim(2)
          DO l1=1,idim(1)
            field3D(l1,l2,l3) = convert_luc_ICON2GRIB(lc_datbase  = post_op%arg1%ival, &
              &                                       iluc_in     = field3D(l1,l2,l3), &
              &                                       opt_linverse= opt_inverse)
          END DO
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
#endif
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

  END SUBROUTINE perform_post_op_i3D


END MODULE mo_post_op
