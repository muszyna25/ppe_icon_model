!--------------------------------------------------------------------
!
! Serialization routine for 2 moment microphysics using Serialbox2
!
!--------------------------------------------------------------------

MODULE mo_ser_nh_diffusion

  USE mo_kind,        ONLY: vp, wp
  USE mo_ser_common,  ONLY: init
  IMPLICIT NONE

  LOGICAL :: writeIn = .TRUE.
  LOGICAL :: writeOut = .TRUE.

  PUBLIC :: serialize_input
  PUBLIC :: serialize_output

  CONTAINS

  SUBROUTINE serialize_input(z_temp, z_nabla2_c, z_nabla2_e, z_nabla4_e)
    IMPLICIT NONE
    REAL(KIND=vp), DIMENSION(:,:,:) :: z_temp, z_nabla2_c
    REAL(KIND=wp), DIMENSION(:,:,:) :: z_nabla2_e, z_nabla4_e

    !$ser verbatim IF (writeIn) THEN
    !$ser verbatim call init()
    !$ser savepoint nh_diffusion-input
    !$ser verbatim #if defined SERIALIZE_CREATE_REFERENCE 
    !$ser mode write
    !$ser verbatim #elif defined SERIALIZE_PERTURB_REFERENCE
    !$ser mode read-perturb
    !$ser verbatim #elif defined SERIALIZE_READ_REFERENCE
    !$ser mode read
    !$ser verbatim #else
    !$ser verbatim #error SERIALIZATION MODE IS NOT SET
    !$ser verbatim #endif 
    !$ser data z_temp=z_temp z_nabla2_c=z_nabla2_c &
    !$ser&     z_nabla2_e=z_nabla2_e z_nabla4_e=z_nabla4_e
    !$ser verbatim writeIn = .FALSE.
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(z_temp, z_nabla2_c, z_nabla2_e, z_nabla4_e)
    IMPLICIT NONE
    REAL(KIND=vp), DIMENSION(:,:,:) :: z_temp, z_nabla2_c
    REAL(KIND=wp), DIMENSION(:,:,:) :: z_nabla2_e, z_nabla4_e

    !$ser verbatim if (writeOut) then
    !$ser verbatim call init()
    !$ser savepoint nh_diffusion-output
    !$ser mode write
    !$ser data z_temp=z_temp z_nabla2_c=z_nabla2_c
    !$ser data z_nabla2_e=z_nabla2_e z_nabla4_e=z_nabla4_e
    !$ser verbatim writeOut = .FALSE.
    !$ser verbatim endif

  END SUBROUTINE serialize_output

END MODULE mo_ser_nh_diffusion
