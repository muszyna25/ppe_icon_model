#include "icon_definitions.inc"

! contains abstact type for actual solver backends
! this is an abstract interposer layer, 
! in order to use a single interface for all backend solvers

MODULE mo_ocean_solve_backend
  !-------------------------------------------------------------------------
  USE mo_exception, ONLY: finish
  USE mo_kind, ONLY: sp, wp
  USE mo_mpi, ONLY: p_barrier, p_comm_work
  USE mo_ocean_solve_lhs, ONLY: t_lhs
  USE mo_ocean_solve_lhs_type, ONLY: t_lhs_agen
  USE mo_ocean_solve_transfer, ONLY: t_transfer
  USE mo_ocean_solve_aux, ONLY: t_destructible, t_ocean_solve_parm
  USE mo_timer, ONLY: timer_start, timer_stop, new_timer
  USE mo_run_config, ONLY: ltimer
 
  IMPLICIT NONE
  
  PRIVATE
 
  PUBLIC :: t_ocean_solve_backend
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_solve_backend'

! abstract solver backend type
  TYPE, ABSTRACT, EXTENDS(t_destructible) :: t_ocean_solve_backend
! typeIDs, dimensions, max-iters
    INTEGER :: timer_wait
    TYPE(t_ocean_solve_parm) :: par, par_sp
! lhs object
    TYPE(t_lhs) :: lhs
! pointer to communication infrastructure object (pointer, since this may be used by multiple solvers)
    CLASS(t_transfer), POINTER :: trans => NULL()
! tolerances
    REAL(KIND=wp) :: abs_tol_wp, abs_tol_sp
! internal arrays (sp) used in all backends (local portion and solver side portion)
    REAL(KIND=sp), ALLOCATABLE, DIMENSION(:,:) :: x_sp, b_sp
! internal arrays (wp) used in all backends (local portion and solver side portion)
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: x_wp, b_wp
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: res_wp
    REAL(KIND=wp), POINTER, DIMENSION(:,:) :: x_loc_wp, b_loc_wp
    REAL(KIND=wp), POINTER, DIMENSION(:) :: res_loc_wp
! internal arrays (int)
    INTEGER, POINTER, DIMENSION(:) :: niter
    INTEGER, ALLOCATABLE, DIMENSION(:) :: niter_cal
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: x_sp, b_sp
!DIR$ ATTRIBUTES ALIGN : 64 :: x_wp, b_wp, res_wp, x_loc_wp, res_loc_wp, niter, niter_cal
#endif
! interfaces 
  CONTAINS
    PROCEDURE :: dump_matrix => ocean_solve_backend_dump_matrix
    PROCEDURE :: construct => ocean_solve_backend_construct
    PROCEDURE(a_solve_backend_wp), DEFERRED :: doit_wp ! call actual solve (wp)
    PROCEDURE(a_solve_backend_sp), DEFERRED :: doit_sp ! call actual solve (sp)
    PROCEDURE :: solve => ocean_solve_backend_solve
    PROCEDURE :: destruct_commons => ocean_solve_backend_destruct_commons
  END TYPE t_ocean_solve_backend

! abstract interfaces to be declared in extended solver types
  ABSTRACT INTERFACE
    SUBROUTINE a_solve_backend_wp(this)
      IMPORT t_ocean_solve_backend
      CLASS(t_ocean_solve_backend), INTENT(INOUT) :: this
    END SUBROUTINE a_solve_backend_wp
    SUBROUTINE a_solve_backend_sp(this)
      IMPORT t_ocean_solve_backend
      CLASS(t_ocean_solve_backend), INTENT(INOUT) :: this
    END SUBROUTINE a_solve_backend_sp
  END INTERFACE

CONTAINS

! dump lhs or preconditioner matrix
  SUBROUTINE ocean_solve_backend_dump_matrix(this, id, lprecon)
    CLASS(t_ocean_solve_backend), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: id
    LOGICAL, INTENT(IN) :: lprecon
    CHARACTER(LEN=*), PARAMETER :: routine = this_mod_name // &
      & "::ocean_solve_t::ocean_solve_dump_matrix()"

    IF (.NOT.ASSOCIATED(this%trans)) CALL finish(routine, &
      & "ocean_solve_t was not initialized")
    IF (lprecon) THEN
      CALL this%lhs%dump_matrix(id, "ocean_matrix_precon_", .true.)
    ELSE
      CALL this%lhs%dump_matrix(id, "ocean_matrix_lhs_", .false.)
    END IF
  END SUBROUTINE ocean_solve_backend_dump_matrix

! initialize data and arrays
  SUBROUTINE ocean_solve_backend_construct(this, par, par_sp, lhs_agen, trans)
    CLASS(t_ocean_solve_backend), INTENT(INOUT) :: this
    TYPE(t_ocean_solve_parm), INTENT(IN) :: par, par_sp
    CLASS(t_lhs_agen), POINTER, INTENT(IN) :: lhs_agen
    CLASS(t_transfer), POINTER, INTENT(IN) :: trans
    CHARACTER(LEN=*), PARAMETER :: routine = this_mod_name// &
      & "::ocean_solve_t::ocean_solve_backend_construct()"

    IF (ASSOCIATED(this%trans)) CALL finish(routine, &
      & "already initialized!")
    CALL this%destruct()
    this%par = par
    this%par_sp = par_sp
    this%abs_tol_wp = par%tol
    IF (par_sp%nidx .NE. -1) this%abs_tol_sp = REAL(par_sp%tol, sp)
! rhs-pointer
    NULLIFY(this%b_loc_wp)
! iters used arrays
    ALLOCATE(this%niter(2), this%niter_cal(2), this%res_wp(2))
! communication infrastructure
    this%trans => trans
! initialize lhs-object
    CALL this%lhs%construct(par_sp%nidx .EQ. par%nidx, par, lhs_agen, trans)
    this%timer_wait = new_timer("solve_wait")
  END SUBROUTINE ocean_solve_backend_construct

! general solve interface (decides wether to use sp- or wp-variant)
  SUBROUTINE ocean_solve_backend_solve(this, niter, niter_sp, upd)
    CLASS(t_ocean_solve_backend), INTENT(INOUT) :: this
    INTEGER, INTENT(OUT) :: niter, niter_sp
    INTEGER, INTENT(IN) :: upd
    CHARACTER(LEN=*), PARAMETER :: routine = this_mod_name// &
      & '::ocean_solve_t::ocean_solve'
    INTEGER :: sum_it, n_re, n_it

    IF (.NOT.ASSOCIATED(this%trans)) &
      & CALL finish(routine, "solve needs to be initialized")
    IF (upd .EQ. 1) CALL this%lhs%update()
! transfer input fields to internal solver arrays (from worker-PEs to solver-PEs)
    IF (.NOT.ALLOCATED(this%x_wp)) THEN ! must also allocate
      CALL this%trans%into_once(this%x_loc_wp, this%x_wp, 1)
      CALL this%trans%into_once(this%b_loc_wp, this%b_wp, 1)
    ELSE ! transfer/copy only
      CALL this%trans%into(this%x_loc_wp, this%x_wp, this%b_loc_wp, this%b_wp, 1)
    END IF
    this%niter_cal(2) = -2
    IF (this%par_sp%nidx .EQ. this%par%nidx .AND. this%trans%is_solver_pe) THEN
      IF (.NOT.ALLOCATED(this%x_sp)) & ! alloc sp-arrays, if not done, yet
        ALLOCATE(this%x_sp(this%trans%nidx, this%trans%nblk_a), &
          & this%b_sp(this%trans%nidx, this%trans%nblk_a))
! convert wp-input to sp
      this%x_sp(:,:) = REAL(this%x_wp(:,:), sp)
      this%b_sp(:,:) = REAL(this%b_wp(:,:), sp)
      sum_it = 0
      n_re = 0
      n_it = -1
      DO WHILE(n_it .EQ. -1 .AND. n_re .LT. this%par_sp%nr)
        CALL this%doit_sp()
        n_it = this%niter_cal(2)
        sum_it = sum_it + MERGE(n_it, this%par_sp%m, n_it .GT. -1)
        IF (this%abs_tol_sp .LT. this%res_wp(2)) n_it = -1
        n_re = n_re + 1
      END DO ! WHILE(tolerance >= solver_tolerance)
      this%niter_cal(2) = sum_it
      this%x_wp(:,:) = REAL(this%x_sp(:,:), wp)
    END IF
    IF (this%trans%is_solver_pe) THEN
      sum_it = 0
      n_re = 0
      n_it = -1
      DO WHILE(n_it .EQ. -1 .AND. n_re .LT. this%par%nr)
        CALL this%doit_wp()
        n_it = this%niter_cal(1)
        sum_it = sum_it + MERGE(n_it, this%par%m, n_it .GT. -1)
        IF (this%abs_tol_wp .LT. this%res_wp(1)) n_it = -1
        n_re = n_re + 1
      END DO ! WHILE(tolerance >= solver_tolerance)
      this%niter_cal(1) = sum_it
    END IF
! transfer solution and residuals from solver-PEs back onto worker-PEs
    IF (ltimer) CALL timer_start(this%timer_wait)
    IF (ltimer) CALL p_barrier(p_comm_work)
    IF (ltimer) CALL timer_stop(this%timer_wait)
    CALL this%trans%sctr(this%x_wp, this%x_loc_wp)
    CALL this%trans%bcst(this%res_wp, this%res_loc_wp)
    CALL this%trans%bcst(this%niter_cal, this%niter)
    niter = this%niter(1)
    niter_sp = this%niter(2)
  END SUBROUTINE ocean_solve_backend_solve

  SUBROUTINE ocean_solve_backend_destruct_commons(this)
    CLASS(t_ocean_solve_backend), INTENT(INOUT) :: this

    CALL this%lhs%destruct()
    NULLIFY(this%trans)
    IF (ALLOCATED(this%niter_cal)) &
      & DEALLOCATE(this%niter, this%niter_cal, this%res_wp)
    IF (ALLOCATED(this%x_wp)) DEALLOCATE(this%x_wp, this%b_wp)
  END SUBROUTINE ocean_solve_backend_destruct_commons

END MODULE mo_ocean_solve_backend
