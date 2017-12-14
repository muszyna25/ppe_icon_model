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
    CHARACTER(LEN=16)     :: prefix
    REAL(KIND=8) :: rprecision
    rprecision = 10.0**(-PRECISION(1.0))

    !$ser verbatim #if defined SERIALIZE_CREATE_REFERENCE 
    !$ser verbatim prefix = 'reference'
    !$ser verbatim #else 
    !$ser verbatim prefix = 'current'
    !$ser verbatim #endif 
    !$ser init directory='./ser_data' &
    !$ser&     prefix=TRIM(prefix) &
    !$ser&     prefix_ref='reference' &
    !$ser&     mpi_rank=get_my_mpi_work_id() &
    !$ser&     rprecision=rprecision &
    !$ser&     rperturb=1.0e-5_8
  END SUBROUTINE init

END MODULE mo_ser_common
