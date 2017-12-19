!--------------------------------------------------------------------
!
!          Serialization routine for psrad using Serialbox2
!
!--------------------------------------------------------------------

MODULE mo_ser_psrad

  USE mo_kind,        ONLY: vp, wp
  USE mo_ser_common,  ONLY: init

  IMPLICIT NONE
  PUBLIC :: serialize_input
  PUBLIC :: serialize_output

  CONTAINS

  SUBROUTINE serialize_input(a)
    IMPLICIT NONE
    REAL(KIND=wp), DIMENSION(:,:) :: a

    !$ser verbatim call init()
    !$ser savepoint psrad-input
#if defined SERIALIZE_CREATE_REFERENCE 
    !$ser mode write
#elif defined SERIALIZE_PERTURB_REFERENCE
    !$ser mode read-perturb
#elif defined SERIALIZE_READ_REFERENCE
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif 
    !$ser data a=a 

  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(a)
    IMPLICIT NONE
    REAL(KIND=wp), DIMENSION(:,:,:) :: a

    !$ser verbatim call init()
    !$ser savepoint psrad-output
    !$ser mode write
    !$ser data a=a

  END SUBROUTINE serialize_output

END MODULE mo_ser_psrad
