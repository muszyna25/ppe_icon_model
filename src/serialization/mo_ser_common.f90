!--------------------------------------------------------------------
!
! Common serialization routines using Serialbox2
!
!--------------------------------------------------------------------

MODULE mo_ser_common

  USE mo_kind,    ONLY: vp, wp
  USE mo_mpi,     ONLY: get_my_mpi_work_id
  IMPLICIT NONE

  PUBLIC :: init

  CONTAINS

  SUBROUTINE init()
    IMPLICIT NONE
    REAL(KIND=8) :: rprecision
    rprecision = 10.0**(-PRECISION(1.0))

    !$ser verbatim #if defined SERIALIZE_CREATE_REFERENCE 
    !$ser init directory='./ser_data' &
    !$ser&     prefix='reference' &
    !$ser&     mpi_rank=get_my_mpi_work_id() &
    !$ser&     rprecision=rprecision &
    !$ser&     rperturb=1.0e-5_8
    !$ser verbatim #else 
    !$ser init directory='./ser_data' &
    !$ser&     prefix='current' &
    !$ser&     prefix_ref='reference' &
    !$ser&     mpi_rank=get_my_mpi_work_id() &
    !$ser&     rprecision=rprecision &
    !$ser&     rperturb=1.0e-5_8
    !$ser verbatim #endif 
  END SUBROUTINE init

END MODULE mo_ser_common
