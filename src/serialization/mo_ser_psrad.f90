!--------------------------------------------------------------------
!
!          Serialization routine for psrad using Serialbox2
!
!--------------------------------------------------------------------

MODULE mo_ser_psrad

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE serialize_psrad(a)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:,:) :: a

    !$ser init directory='.' prefix='Serialization_psrad'
    !$ser savepoint sp1
    !$ser mode write
    !$ser data ser_a=a

  END SUBROUTINE serialize_psrad

  SUBROUTINE deserialize_psrad(a)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:,:) :: a

    !$ser init directory='.' prefix='Serialization_psrad-output' prefix_ref='Serialization_psrad'
    !$ser savepoint sp1
    !$ser mode read
    !$ser data ser_a=a
    !$ser mode write
    !$ser data ser_a=a

  END SUBROUTINE deserialize_psrad

  SUBROUTINE deserialize_with_perturb_psrad(a)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:,:) :: a
    REAL(KIND=8) :: rprecision
    rprecision = 10.0**(-PRECISION(1.0))

    !$ser init directory='.' prefix='Serialization_psrad-output' prefix_ref='Serialization_psrad' rprecision=rprecision rperturb=1.0e-5_8
    !$ser savepoint sp1
    !$ser mode read-perturb
    !$ser data ser_a=a

  END SUBROUTINE deserialize_with_perturb_psrad

END MODULE mo_ser_psrad
