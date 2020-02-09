!--------------------------------------------------------------------
!
! Common serialization routines using Serialbox2
!
!--------------------------------------------------------------------

MODULE mo_ser_common

  USE mo_kind,    ONLY: vp, wp
  USE mo_mpi,     ONLY: get_my_mpi_work_id
  IMPLICIT NONE

  LOGICAL :: linitialize = .TRUE.

  PUBLIC :: init

  CONTAINS

  SUBROUTINE init(suffix)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: suffix
    REAL(KIND=8) :: rprecision
    rprecision = 10.0**(-PRECISION(1.0))

    !$ser verbatim    IF (linitialize) THEN

#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser init directory='./ser_data' &
    !$ser&     prefix='reference_'//TRIM(suffix) &
    !$ser&     mpi_rank=get_my_mpi_work_id() &
    !$ser&     rprecision=rprecision &
    !$ser&     rperturb=1.0e-5_8
#else
    !$ser init directory='./ser_data' &
    !$ser&     prefix='current_'//TRIM(suffix) &
    !$ser&     prefix_ref='reference_'//TRIM(suffix) &
    !$ser&     mpi_rank=get_my_mpi_work_id() &
    !$ser&     rprecision=rprecision &
    !$ser&     rperturb=1.0e-5_8
#endif

    !$ser verbatim     linitialize = .FALSE.
    !$ser verbatim     END IF

  END SUBROUTINE init

END MODULE mo_ser_common
