MODULE mo_ocean_solve_lhs_type

! abstract type for lhs-matrix generators

  USE mo_ocean_solve_aux, ONLY: t_destructible
  USE mo_exception, ONLY: finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_lhs_agen, lhs_agen_ptr

  TYPE, ABSTRACT, EXTENDS(t_destructible) :: t_lhs_agen
    LOGICAL :: is_const, use_shortcut
    INTEGER, ALLOCATABLE, DIMENSION(:) :: is_init
  CONTAINS
    PROCEDURE(a_lhs_agen_wp), DEFERRED :: lhs_wp
    PROCEDURE(a_lhs_matrix_shortcut), DEFERRED :: lhs_matrix_shortcut
    GENERIC :: apply => lhs_wp
    GENERIC :: matrix_shortcut => lhs_matrix_shortcut
  END TYPE t_lhs_agen

  ABSTRACT INTERFACE
    SUBROUTINE a_lhs_agen_wp(this, x, ax)
      USE mo_kind, ONLY: wp
      IMPORT t_lhs_agen
      CLASS(t_lhs_agen), INTENT(INOUT) :: this
      REAL(KIND=wp), INTENT(IN) :: x(:,:)
      REAL(KIND=wp), INTENT(OUT) :: ax(:,:)
    END SUBROUTINE a_lhs_agen_wp
    SUBROUTINE a_lhs_matrix_shortcut(this, idx, blk, coeff)
      USE mo_kind, ONLY: wp
      IMPORT t_lhs_agen
      CLASS(t_lhs_agen), INTENT(INOUT) :: this
      INTEGER, INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:,:) :: idx, blk
      REAL(KIND=wp), INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:,:) :: coeff
    END SUBROUTINE a_lhs_matrix_shortcut
  END INTERFACE

CONTAINS

  FUNCTION lhs_agen_ptr(this) RESULT(this_ptr)
    CLASS(t_destructible), TARGET, INTENT(INOUT) :: this
    CLASS(t_lhs_agen), POINTER :: this_ptr

    SELECT TYPE (this)
    CLASS IS (t_lhs_agen)
      this_ptr => this
    CLASS DEFAULT
      NULLIFY(this_ptr)
      CALL finish("trivial_transfer_ptr", "not correct type!")
    END SELECT
  END FUNCTION lhs_agen_ptr

END MODULE mo_ocean_solve_lhs_type
