#if (defined(_OPENMP) && defined(OCE_SOLVE_OMP))
#include "omp_definitions.inc"
#endif
! contains extension to solver backend type: CG

MODULE mo_ocean_solve_cg

  USE mo_kind, ONLY: wp, sp
  USE mo_exception, ONLY: finish
  USE mo_ocean_solve_backend, ONLY: t_ocean_solve_backend

  IMPLICIT NONE
  
  PRIVATE
 
  PUBLIC :: t_ocean_solve_cg
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_solve_cg'

  TYPE, EXTENDS(t_ocean_solve_backend) :: t_ocean_solve_cg
    PRIVATE
! arrays only used by CG
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: z_wp, d_wp, r_wp, rsq_wp
    REAL(KIND=sp), ALLOCATABLE, DIMENSION(:,:) :: z_sp, d_sp, r_sp, rsq_sp
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: z_wp, d_wp, r_wp, rsq_wp
!DIR$ ATTRIBUTES ALIGN : 64 :: z_sp, d_sp, r_sp, rsq_sp
#endif
! interfaces
  CONTAINS
    PROCEDURE :: doit_wp => ocean_solve_cg_cal_wp ! override deferred
    PROCEDURE :: doit_sp => ocean_solve_cg_cal_sp ! override deferred
    PROCEDURE :: destruct => ocean_solve_cg_destruct ! override deferred
    PROCEDURE, PRIVATE :: recover_arrays_wp => ocean_solve_cg_recover_arrays_wp
    PROCEDURE, PRIVATE :: recover_arrays_sp => ocean_solve_cg_recover_arrays_sp
    GENERIC, PRIVATE :: recover_arrays => recover_arrays_wp, recover_arrays_sp
  END TYPE t_ocean_solve_cg

CONTAINS

! destruct lhs-object and de-allocate arrays
  SUBROUTINE ocean_solve_cg_destruct(this)
    CLASS(t_ocean_solve_cg), INTENT(INOUT) :: this

    CALL this%destruct_commons()
    IF (ALLOCATED(this%z_wp)) DEALLOCATE(this%z_wp, this%d_wp, this%r_wp, &
      & this%rsq_wp)
    IF (ALLOCATED(this%x_sp)) DEALLOCATE(this%x_sp, this%b_sp, this%z_sp, &
      & this%d_sp, this%r_sp, this%rsq_sp)
  END SUBROUTINE ocean_solve_cg_destruct

! get solver arrays (alloc them, if not done so, yet) - wp-variant
  SUBROUTINE ocean_solve_cg_recover_arrays_wp(this, x, b, z, d, r, r2)
    CLASS(t_ocean_solve_cg), INTENT(INOUT), TARGET :: this
    REAL(KIND=wp), INTENT(OUT), POINTER, DIMENSION(:,:) :: &
      & x, b, z, d, r, r2

    IF (.NOT.ALLOCATED(this%z_wp)) THEN
      ALLOCATE(this%z_wp(this%trans%nidx, this%trans%nblk), &
        & this%d_wp(this%trans%nidx, this%trans%nblk_a), &
        & this%r_wp(this%trans%nidx, this%trans%nblk), &
        & this%rsq_wp(this%trans%nidx, this%trans%nblk))
      this%d_wp(:, this%trans%nblk:this%trans%nblk_a) = 0._wp
    END IF
    x => this%x_wp
    b => this%b_wp
    z => this%z_wp
    d => this%d_wp
    r => this%r_wp
    r2 => this%rsq_wp
  END SUBROUTINE ocean_solve_cg_recover_arrays_wp

! actual CG solve (vanilla) - wp-variant
SUBROUTINE ocean_solve_cg_cal_wp(this)
    CLASS(t_ocean_solve_cg), INTENT(INOUT) :: this
    REAL(KIND=wp) :: alpha, beta, dz_glob, tol, tol2, rn, rn_last
    INTEGER :: nidx_e, nblk, iblk, k, m, k_final
    REAL(KIND=wp), POINTER, DIMENSION(:,:), CONTIGUOUS :: &
      & x, b, z, d, r, r2
    LOGICAL :: done

! retrieve extends of vector to solve
    nblk = this%trans%nblk
    nidx_e = this%trans%nidx_e
    m = this%par%m
    k_final = -1
! retrieve arrays
    CALL this%recover_arrays(x, b, z, d, r, r2)
    b(nidx_e+1:, nblk) = 0._wp
! compute initial residual and auxiliary vectors
    CALL this%trans%sync(x)
    CALL this%lhs%apply(x, z)
    z(nidx_e+1:, nblk) = 0._wp
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO iblk = 1, nblk
      r(:, iblk) = b(:, iblk) - z(:, iblk)
      d(:, iblk) = r(:, iblk)
      r2(:, iblk) = r(:, iblk) * r(:, iblk)
    END DO
!ICON_OMP END PARALLEL DO
    CALL this%trans%global_sum(r2, rn)
    ! tolerance
    tol = this%abs_tol_wp
    IF (.NOT.this%par%use_atol) THEN
      tol = this%par%tol * SQRT(rn)
      this%abs_tol_wp = tol
    END IF
    tol2 = tol * tol
    done = .false.
! enter CG iteration
    DO k = 1, m
! check if done
      IF (done) CYCLE
!      IF (this%trans%is_leader_pe) PRINT*,"it, res",k,SQRT(rn)
! check if already reached desired tolerance
      IF (rn .LE. tol2) THEN
        done = .true.
        k_final = k
        CYCLE
      END IF
! correct search direction (in direction of gradient) / update search vector
      IF (k .GT. 1) THEN
        beta = rn / rn_last
        d(nidx_e+1:, nblk) = 0._wp
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
        DO iblk = 1, nblk
          d(:, iblk) = r(:, iblk) + beta * d(:, iblk)
        END DO
!ICON_OMP END PARALLEL DO
      END IF
      CALL this%trans%sync(d)
      CALL this%lhs%apply(d, z)
      d(nidx_e+1:, nblk) = 0._wp
! compute extrapolated location of minimum in direction of d
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        r2(:, iblk) = d(:, iblk) * z(:, iblk)
      END DO
!ICON_OMP END PARALLEL DO
      r2(nidx_e+1:, nblk) = 0._wp
      CALL this%trans%global_sum(r2, dz_glob)
      alpha = rn / dz_glob
! update guess and residuum
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        x(:, iblk) = x(:, iblk) + alpha * d(:, iblk)
        r(:, iblk) = r(:, iblk) - alpha * z(:, iblk)
        r2(:, iblk) = r(:, iblk) * r(:, iblk)
      END DO
!ICON_OMP END PARALLEL DO
      r2(nidx_e+1:, nblk) = 0._wp
! save old and compute new residual norm
      rn_last = rn
      CALL this%trans%global_sum(r2, rn)
    END DO
    this%niter_cal(1) = k_final
    this%res_wp(1) = SQRT(rn)
    CALL this%trans%sync(x)
  END SUBROUTINE ocean_solve_cg_cal_wp

! get solver arrays (alloc them, if not done so, yet) - sp-variant
  SUBROUTINE ocean_solve_cg_recover_arrays_sp(this, x, b, z, d, r, r2)
    CLASS(t_ocean_solve_cg), INTENT(INOUT), TARGET :: this
    REAL(KIND=sp), INTENT(OUT), POINTER, DIMENSION(:,:) :: &
      & x, b, z, d, r, r2

    IF (.NOT.ALLOCATED(this%z_sp)) THEN
      ALLOCATE(this%z_sp(this%trans%nidx, this%trans%nblk), &
        & this%d_sp(this%trans%nidx, this%trans%nblk_a), &
        & this%r_sp(this%trans%nidx, this%trans%nblk), &
        & this%rsq_sp(this%trans%nidx, this%trans%nblk))
      this%d_sp(:, this%trans%nblk:this%trans%nblk_a) = 0._sp
    END IF
    x => this%x_sp
    b => this%b_sp
    z => this%z_sp
    d => this%d_sp
    r => this%r_sp
    r2 => this%rsq_sp
  END SUBROUTINE ocean_solve_cg_recover_arrays_sp

! we should not get here...
  SUBROUTINE ocean_solve_cg_cal_sp(this)
    CLASS(t_ocean_solve_cg), INTENT(INOUT) :: this
!    CALL finish(TRIM(this_mod_name)//":ocean_solve_cg_cal_sp()", &
!      & "single precision cg solver not implemented")
    REAL(KIND=sp) :: alpha, beta, dz_glob, tol, tol2, rn, rn_last
    INTEGER :: nidx_e, nblk, iblk, k, m, k_final
    REAL(KIND=sp), POINTER, DIMENSION(:,:), CONTIGUOUS :: &
      & x, b, z, d, r, r2
    LOGICAL :: done

! retrieve extends of vector to solve
    nblk = this%trans%nblk
    nidx_e = this%trans%nidx_e
    m = this%par%m
    k_final = -1
! retrieve arrays
    CALL this%recover_arrays(x, b, z, d, r, r2)
    b(nidx_e+1:, nblk) = 0._sp
! compute initial residual and auxiliary vectors
    CALL this%trans%sync(x)
    CALL this%lhs%apply(x, z)
    z(nidx_e+1:, nblk) = 0._sp
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO iblk = 1, nblk
      r(:, iblk) = b(:, iblk) - z(:, iblk)
      d(:, iblk) = r(:, iblk)
      r2(:, iblk) = r(:, iblk) * r(:, iblk)
    END DO
!ICON_OMP END PARALLEL DO
    CALL this%trans%global_sum(r2, rn)
    ! tolerance
    tol = REAL(this%abs_tol_sp, sp)
    IF (.NOT.this%par%use_atol) THEN
      tol = this%par%tol * SQRT(rn)
      this%abs_tol_sp = REAL(tol, wp)
    END IF
    tol2 = tol * tol
    done = .false.
! enter CG iteration
    DO k = 1, m
! check if done
      IF (done) CYCLE
!      IF (this%trans%is_leader_pe) PRINT*,"it, res",k,SQRT(rn)
! check if already reached desired tolerance
      IF (rn .LE. tol2) THEN
        done = .true.
        k_final = k
        CYCLE
      END IF
! correct search direction (in direction of gradient) / update search vector
      IF (k .GT. 1) THEN
        beta = rn / rn_last
        d(nidx_e+1:, nblk) = 0._sp
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
        DO iblk = 1, nblk
          d(:, iblk) = r(:, iblk) + beta * d(:, iblk)
        END DO
!ICON_OMP END PARALLEL DO
      END IF
      CALL this%trans%sync(d)
      CALL this%lhs%apply(d, z)
      d(nidx_e+1:, nblk) = 0._sp
! compute extrapolated location of minimum in direction of d
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        r2(:, iblk) = d(:, iblk) * z(:, iblk)
      END DO
!ICON_OMP END PARALLEL DO
      r2(nidx_e+1:, nblk) = 0._sp
      CALL this%trans%global_sum(r2, dz_glob)
      alpha = rn / dz_glob
! update guess and residuum
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        x(:, iblk) = x(:, iblk) + alpha * d(:, iblk)
        r(:, iblk) = r(:, iblk) - alpha * z(:, iblk)
        r2(:, iblk) = r(:, iblk) * r(:, iblk)
      END DO
!ICON_OMP END PARALLEL DO
      r2(nidx_e+1:, nblk) = 0._sp
! save old and compute new residual norm
      rn_last = rn
      CALL this%trans%global_sum(r2, rn)
    END DO
    this%niter_cal(2) = k_final
    this%res_wp(2) = REAL(SQRT(rn), wp)
    CALL this%trans%sync(x)
  END SUBROUTINE ocean_solve_cg_cal_sp

END MODULE mo_ocean_solve_cg
