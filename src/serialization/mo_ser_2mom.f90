!--------------------------------------------------------------------
!
! Serialization routine for 2 moment microphysics using Serialbox2
!
!--------------------------------------------------------------------

MODULE mo_ser_2mom

  USE mo_kind,        ONLY: vp, wp
  USE mo_ser_common,  ONLY: init

  IMPLICIT NONE
  PUBLIC :: serialize_input
  PUBLIC :: serialize_output

  CONTAINS

  SUBROUTINE serialize_input(dz, rho, pres, w)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: dz, rho, pres, w

    !$ser verbatim call init()
    !$ser savepoint 2mom-input
    !$ser verbatim #if defined SERIALIZE_CREATE_REFERENCE 
    !$ser mode write
    !$ser verbatim #elif defined SERIALIZE_PERTURB_REFERENCE
    !$ser mode read-perturb
    !$ser verbatim #elif defined SERIALIZE_READ_REFERENCE
    !$ser mode read
    !$ser verbatim #else
    !$ser verbatim #error SERIALIZATION MODE IS NOT SET
    !$ser verbatim #endif 
    !$ser data dz=dz rho=rho pres=pres w=w

  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(dz, rho, pres, w)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: dz, rho, pres, w

    !$ser verbatim call init()
    !$ser savepoint 2mom-output
    !$ser mode write
    !$ser data dz=dz rho=rho pres=pres w=w

  END SUBROUTINE serialize_output

END MODULE mo_ser_2mom
