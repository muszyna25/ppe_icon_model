MODULE mo_primal_flip_flop_lhs

! contains extended lhs-matrix-generator type
! provides the primal flip flop lhs for the "mass matrix inversion" - solve

  USE mo_exception, ONLY: finish
  USE mo_kind, ONLY: wp
  USE mo_ocean_solve_lhs_type, ONLY: t_lhs_agen
  USE mo_ocean_solve_aux, ONLY: t_destructible
  USE mo_model_domain, ONLY: t_patch_3d, t_patch
  USE mo_scalar_product, ONLY: map_edges2edges_viacell_2D_per_level
  USE mo_ocean_types, ONLY: t_operator_coeff

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_primal_flip_flop_lhs, lhs_primal_flip_flop_ptr

  TYPE, EXTENDS(t_lhs_agen) :: t_primal_flip_flop_lhs
    PRIVATE
    TYPE(t_patch_3d), POINTER :: patch_3d => NULL()
    TYPE(t_patch), POINTER :: patch_2d => NULL()
    TYPE(t_operator_coeff), POINTER :: op_coeffs
    INTEGER :: jk
  CONTAINS
    PROCEDURE :: lhs_wp => lhs_primal_flip_flop_wp ! deferred override
    PROCEDURE :: destruct => lhs_primal_flip_flop_destruct ! deferred override
    PROCEDURE :: construct => lhs_primal_flip_flop_construct ! deferred override
    PROCEDURE :: lhs_matrix_shortcut => lhs_primal_flip_flop_matrix_shortcut
  END TYPE t_primal_flip_flop_lhs

CONTAINS

! returns pointer to t_primal_flip_flop_lhs object, if provided a matching object of corresponding abstract type
  FUNCTION lhs_primal_flip_flop_ptr(this) RESULT(this_ptr)
    CLASS(t_destructible), INTENT(IN), TARGET :: this
    CLASS(t_primal_flip_flop_lhs), POINTER :: this_ptr
    
    SELECT TYPE (this)
    CLASS IS (t_primal_flip_flop_lhs)
      this_ptr => this
    CLASS DEFAULT
      NULLIFY(this_ptr)
      CALL finish("primal_flip_flop_ptr()", "not correct type!")
    END SELECT
  END FUNCTION lhs_primal_flip_flop_ptr

! init generator object
  SUBROUTINE lhs_primal_flip_flop_construct(this, patch_3d, op_coeffs, jk)
    CLASS(t_primal_flip_flop_lhs), INTENT(INOUT) :: this
    TYPE(t_patch_3d), TARGET, INTENT(IN) :: patch_3d
    TYPE(t_operator_coeff), TARGET, INTENT(IN) :: op_coeffs
    INTEGER, INTENT(IN) :: jk

    CALL this%destruct()
    this%is_const = .false.
    this%patch_3d => patch_3d
    this%patch_2d => patch_3d%p_patch_2d(1)
    this%op_coeffs => op_coeffs
    this%use_shortcut = .false.
    this%jk = jk
    ALLOCATE(this%is_init(1))
  END SUBROUTINE lhs_primal_flip_flop_construct

! apply operator to vector x
  SUBROUTINE lhs_primal_flip_flop_wp(this, x, ax)
    CLASS(t_primal_flip_flop_lhs), INTENT(INOUT) :: this
    REAL(KIND=wp), INTENT(IN) :: x(:,:)
    REAL(KIND=wp), INTENT(OUT) ::ax(:,:)

    IF(.NOT.ALLOCATED(this%is_init)) &
      CALL finish("lhs_primal_flip_flop_wp()", "not correctly initialized")
    CALL map_edges2edges_viacell_2D_per_level(this%patch_3d, &
      & x(:,:), this%op_coeffs, ax(:,:), this%jk)
  END SUBROUTINE lhs_primal_flip_flop_wp

! clear object internals
  SUBROUTINE lhs_primal_flip_flop_destruct(this)
    CLASS(t_primal_flip_flop_lhs), INTENT(INOUT) :: this

    NULLIFY(this%op_coeffs, this%patch_3d, this%patch_2d)
    IF (ALLOCATED(this%is_init)) DEALLOCATE(this%is_init)
  END SUBROUTINE lhs_primal_flip_flop_destruct

  SUBROUTINE lhs_primal_flip_flop_matrix_shortcut(this, idx, blk, coeff)
    CLASS(t_primal_flip_flop_lhs), INTENT(INOUT) :: this
    INTEGER, INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:,:) :: idx, blk
    REAL(KIND=wp), INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:,:) :: coeff

    CALL finish("lhs_primal_flip_flop_matrix_shortcut", &
      & "not implemented -- go away!")
  END SUBROUTINE lhs_primal_flip_flop_matrix_shortcut

END MODULE mo_primal_flip_flop_lhs
