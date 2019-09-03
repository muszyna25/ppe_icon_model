MODULE mo_ocean_solve_aux

! contains implementation parameters of solver infrastructure
  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  
  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_solve_aux'
  INTEGER, PARAMETER, PUBLIC :: solve_cell = 11, & ! use cell - sync pattern (used in transfer-types)
                                solve_edge = 12, & ! use edge - sync pattern (used in transfer-types)
                                solve_vert = 13    ! use vertex - sync pattern (used in transfer-types)
  INTEGER, PARAMETER, PUBLIC :: solve_gmres = 21, & ! internal ID of GMRES-alike solvers
                                solve_cg = 22, & ! internal ID of CG-alike solvers
                                solve_bcgs = 23, & ! internal ID of BiCG-Stab-alike solvers
                                solve_legacy_gmres = 24, & ! internal ID of legacy GMRES solver
                                solve_mres = 25 ! internal ID of MINRES-alike solvers
  INTEGER, PARAMETER, PUBLIC :: solve_precon_none = 60 ! internal ID of not using any preconditioner
  INTEGER, PARAMETER, PUBLIC :: solve_precon_jac = 61 ! internal ID of Jacobi preconditioner
  INTEGER, PARAMETER, PUBLIC :: solve_trans_scatter = 70, & ! transfer subset mode scatter (every n-th proc is a solve-proc)
                                solve_trans_compact = 71 ! transfer subset mode compact (first p_n_work / n procs are solve procs)
  INTEGER, PARAMETER, PUBLIC :: solve_invalid = -999
  TYPE :: t_ocean_solve_parm
    REAL(KIND=wp) :: tol
    INTEGER :: pt, nr, m, nblk, nblk_a, nidx, nidx_e
    LOGICAL :: use_atol
  CONTAINS
    PROCEDURE :: init => ocean_solve_parm_init
  END TYPE t_ocean_solve_parm

  TYPE, ABSTRACT, PUBLIC :: t_destructible
  CONTAINS
  PROCEDURE(a_destruct), DEFERRED, PUBLIC :: destruct
  END TYPE t_destructible

  ABSTRACT INTERFACE
    SUBROUTINE a_destruct(this)
      IMPORT t_destructible
      CLASS(t_destructible), INTENT(INOUT) :: this
    END SUBROUTINE a_destruct
  END INTERFACE

  PUBLIC :: ocean_solve_clear, t_ocean_solve_parm

  CONTAINS

  SUBROUTINE ocean_solve_clear(this_ptr)
    CLASS(t_destructible), POINTER, INTENT(INOUT) :: this_ptr

    IF (ASSOCIATED(this_ptr)) THEN
      CALL this_ptr%destruct()
      DEALLOCATE(this_ptr)
      NULLIFY(this_ptr)
    END IF
  END SUBROUTINE ocean_solve_clear

  SUBROUTINE ocean_solve_parm_init(this, pt, nr, m, nblk, nblk_a, &
    & nidx, nidx_e, tol, use_atol)
    CLASS(t_ocean_solve_parm), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: pt, nr, m, nblk, nblk_a, nidx, nidx_e
    REAL(KIND=wp), INTENT(IN) :: tol
    LOGICAL :: use_atol

    this%pt = pt
    this%nr = nr
    this%m = m
    this%nblk = nblk
    this%nblk_a = nblk_a
    this%nidx = nidx
    this%nidx_e = nidx_e
    this%tol = tol
    this%use_atol = use_atol
  END SUBROUTINE ocean_solve_parm_init

END MODULE mo_ocean_solve_aux
