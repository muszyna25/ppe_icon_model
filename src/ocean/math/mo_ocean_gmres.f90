!>
!!  General comments: provide an implementation of the GMRES solver for.
!!
!!  General comments: provide an implementation of the GMRES solver for
!!  linear systems with the modified Gram-Schmidt orthogonalization.
!!  Implementation as in "Templates for the Solution of Linear Systems:
!!  Building Blocks for Iterative Methods
!!
!! @par Revision History
!!  Original version from http://www.netlib.org/templates/templates.pdf
!!  F90 rewriting by Marco Restelli.
!!  Adapted for use in ICOHDC by Hui Wan, MPI-M (2007-12-15)
!!  Cosmetic changes following the ICON programming guide by Hui Wan,
!!  MPI-M (2007-12-18)
!!  Included blocking by Marco Restelli (2008-09-23)
!!  Dummy argument for preconditioner shifted to last position
!!  and made optional, Marco Giorgetta (2009-02-14)
!!  Inlining of former functions active_dot and norm2 to simplify OpenMP
!!  parallelization, reduction of parallel sections, Guenther Zaengl (2009-06-16)
!!  Adaption of the nonhydrostatic gmres. It works on edges and has no loop over
!!  vertical layers. Almut Gassmann (2010-07-12)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#ifdef __GMRES_OPENMP__
! #include "omp_definitions.inc"
#endif
#include "icon_definitions.inc"
!----------------------------

MODULE mo_ocean_gmres
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp, sp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ltimer
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_timer,               ONLY: timer_start, timer_stop, timers_level, timer_gmres,   &
    & timer_gmres_p_sum, activate_sync_timers
  USE mo_ocean_types,         ONLY: t_operator_coeff, t_solverCoeff_singlePrecision
  USE mo_sync,                ONLY: omp_global_sum_array, global_sum_array
  !  & sync_e, sync_c, sync_v, sync_patch_array
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_mpi,                 ONLY: get_my_global_mpi_id, p_barrier, p_sum, &
    & get_my_mpi_work_communicator
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  
  IMPLICIT NONE
  
  PRIVATE
  
 
  PUBLIC :: ocean_restart_gmres, gmres_oce_old, ocean_restart_gmres_singlePrecesicion
  PUBLIC :: ocean_restart_gmres_e2e, gmres_oce_e2e
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_gmres'
  
  
  !INTERFACE ocean_restart_gmres
  !
  !  MODULE PROCEDURE ocean_restart_gmres
  !  MODULE PROCEDURE gmres_oce_e2e
  !
  !END INTERFACE
  
  !lk  INTERFACE DOT_PRODUCT
  !lk    MODULE PROCEDURE active_dot
  !lk  END INTERFACE
  
  ! module variables used to specify the array bounds in the private
  ! module subroutine and functions (notice that making active_dot an
  ! internal procedure of gmres would not allow for overloading)
  
  
CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !! Based on gmres_oce_old
  !! restart functionality and optimization L.Linardakis, MPIM, 2013
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE ocean_restart_gmres( x,lhs,        &
    & h_e, old_h, patch_3D,      &
    & p_op_coeff, b,                            &
    & absolute_tolerance,  relative_tolerance,  &
    & m, maxiterex,niter,res,                   &
    & preconditioner)
    !
    ! !DESCRIPTION
    !  GMRES solver for linear systems. Notice that, in principle, r and b
    !  can be of any TYPE for which the following operations are defined:
    !   +, -, dot_product
    !
    REAL(wp), INTENT(inout) :: x(:,:)
    REAL(wp), INTENT(in) :: h_e(:,:)
    REAL(wp), INTENT(in) :: old_h(:,:)
    ! patch info needed for calculating lhs
    !TYPE(t_patch), INTENT(IN) :: patch_2D
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    ! index defining the "active" region of the arrays
    ! parameter used in calculating the lhs
    ! REAL(wp), INTENT(in) :: coeff
    TYPE(t_operator_coeff), INTENT(in):: p_op_coeff
    ! right-hand side: same shape as x
    REAL(wp), INTENT(inout) :: b(:,:) ! same size as x
    REAL(wp), INTENT(inout) :: absolute_tolerance
    REAL(wp), INTENT(in) :: relative_tolerance ! (relative or absolute) relative_tolerance
    ! .false. for relative relative_tolerance
    INTEGER,  INTENT(in) :: m         ! maximum number of iterations
    LOGICAL,  INTENT(out) :: maxiterex ! true if reached m iterations
    INTEGER,  INTENT(out) :: niter    ! number of iterations (defined
    ! as the number of evaluation of
    ! the lhs which, due to the
    ! initial evaluation of the
    ! residual, is equal to the
    ! number of Arnoldi iterations +1)
    ! norms of the residual (convergence history); an argument of
    ! dimension at least m is required
    REAL(wp), INTENT(inout) :: res(:) ! (m)
    
    INTERFACE   ! left-hand-side: A*x
      FUNCTION lhs(x, patch_3D, h_e, p_op_coeff) result(ax)
        USE mo_kind, ONLY: wp
        USE mo_model_domain, ONLY: t_patch, t_patch_3d
        USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
        REAL(wp),    INTENT(inout) :: x(:,:)  ! inout for sync
!         REAL(wp), INTENT(in) :: old_h(:,:)
        !TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
        TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
        ! REAL(wp),    INTENT(in) :: coeff
        REAL(wp),    INTENT(in) :: h_e(:,:)
        TYPE(t_operator_coeff),INTENT(in)  :: p_op_coeff
        !     REAL(wp), POINTER :: ax( :, : ) ! same as x
        REAL(wp) :: ax( SIZE(x,1) , SIZE(x,2) ) ! same as x
      END FUNCTION lhs
    END INTERFACE
    
    INTERFACE   ! preconditioner
      SUBROUTINE preconditioner(r, patch_3D, p_op_coeff,h_e)
        USE mo_kind, ONLY: wp
        USE mo_model_domain, ONLY: t_patch, t_patch_3d
        USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
        REAL(wp), INTENT(inout)           :: r(:,:)
        !TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
        TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
        TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
        REAL(wp),    INTENT(in)           :: h_e(:,:)
      END SUBROUTINE preconditioner
    END INTERFACE
    
    OPTIONAL :: preconditioner
    
    ! !LOCAL VARIABLES
    
    LOGICAL :: done
    INTEGER ::     &
      & i,           & ! index for the Arnoldi loop
      & k,           & ! index for the Gram Schmidt orthogonalization
      & nkry           ! dimension of the Krylov space
    REAL(wp) ::    &
      & tol,         & ! effective relative_tolerance (absolute o relative)
      & tol2,        &
      & r(SIZE(x,1),SIZE(x,2)),  & ! residual
      & rn2(m),      &             ! two-norm of the residual
      & v(SIZE(x,1),SIZE(x,2),m),& ! Krylov basis
      & h(m,m),      &             ! Hessemberg matrix
      & hki, hk1i,   &
      & c(m), s(m),  &             ! rotation matrices
      & den, ci
    
    ! REAL(wp), POINTER :: w(:,:)
    REAL(wp) ::  w(SIZE(x,1),SIZE(x,2)) ! new Krylov vector
    
    REAL(wp) :: rrn2, h_aux, rh
    INTEGER :: jb, jk, nlen
    
    INTEGER :: no_of_blocks, pad_nproma, end_nproma
    INTEGER :: my_mpi_work_communicator
    
    REAL(wp) :: sum_aux(patch_3D%p_patch_2d(1)%cells%in_domain%end_block)
    !   REAL(wp) :: sum_x(0:7), sum_w, sum_v
    ! #else
    !   REAL(wp) :: z(SIZE(x,1),SIZE(x,2)) ! needed for global sums in p_test_run
    ! #endif
    
    INTEGER :: mythreadno
    TYPE(t_patch), POINTER :: patch_2d
    CHARACTER(len=*), PARAMETER :: method_name='ocean_restart_gmres'

    !$ INTEGER OMP_GET_THREAD_NUM
    !-------------------------------------------------------------------------
    !  write(0,*) "--------------- gmres --------------------------"
    
    patch_2d =>  patch_3D%p_patch_2d(1)
    my_mpi_work_communicator = get_my_mpi_work_communicator()
    ! 0) set module variables and initialize maxiterex
    no_of_blocks    = patch_2d%cells%in_domain%end_block
    end_nproma      = patch_2d%cells%in_domain%end_index
    pad_nproma      = patch_2d%cells%in_domain%end_index + 1
    
    maxiterex = .FALSE.
    
    v(:,:,:)  = 0.0_wp
    r(:,:)    = 0.0_wp
    b(pad_nproma:nproma, no_of_blocks) = 0.0_wp
 !   w(:, :) = 0.0_wp
 !   b(:, :) = 0.0_wp
 !   x(:, :) = 0.0_wp
    
    ! 1) compute the preconditioned residual
    IF (PRESENT(preconditioner)) CALL preconditioner(x(:,:),patch_3D,p_op_coeff,h_e)
    IF (PRESENT(preconditioner)) CALL preconditioner(b(:,:),patch_3D,p_op_coeff,h_e)
    
    
    w(:, :) = lhs(x(:,:), patch_3D, h_e, p_op_coeff)

    ! w(pad_nproma:nproma, no_of_blocks) = 0.0_wp
    x(pad_nproma:nproma, no_of_blocks) = 0.0_wp
    
    !    write(0,*) "-----------------------------------"
    !    sum_x(1) = SUM(x(:,:))
    !    sum_w    = SUM(w(:,:))
    !    write(0,*) "gmres sum x(0:1), w:", sum_x(0:1), sum_w
    
    start_timer(timer_gmres,2)
    
    mythreadno = 0
! !ICON_OMP_PARALLEL PRIVATE(rrn2, myThreadNo)
!    !$ myThreadNo = OMP_GET_THREAD_NUM()

! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      r(1:nproma,jb) = b(1:nproma,jb) - w(1:nproma,jb)
    ENDDO
! !ICON_OMP_END_DO
    
!     IF (PRESENT(preconditioner)) CALL preconditioner(r(:,:),patch_3D,p_op_coeff,h_e)
    
! !ICON_OMP_DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      sum_aux(jb) = SUM(r(1:nproma,jb) * r(1:nproma,jb))
    ENDDO
    ! sum_aux(no_of_blocks) = SUM(r(1:end_nproma,no_of_blocks) * r(1:end_nproma,no_of_blocks))
! !ICON_OMP_END_DO
    
    IF (mythreadno == 0) THEN
      h_aux  = SUM(sum_aux(1:no_of_blocks))
      start_detail_timer(timer_gmres_p_sum, 5)
      h_aux = p_sum(h_aux, my_mpi_work_communicator)
      stop_detail_timer(timer_gmres_p_sum, 5)
      rn2(1) = SQRT(h_aux)
      !       rn2(1) = SQRT(p_sum(h_aux, my_mpi_work_communicator))
      ! !ICON_OMP FLUSH(rn2(1))
    ENDIF
! !ICON_OMP BARRIER
    
    
    ! 2) compute the first vector of the Krylov space
    IF (rn2(1) /= 0.0_wp) THEN
      rrn2 = 1.0_wp/rn2(1)
! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, no_of_blocks
        v(1:nproma,jb,1) = r(1:nproma,jb) * rrn2
      ENDDO
! !ICON_OMP_END_DO NOWAIT
    ENDIF
! !ICON_OMP_END_PARALLEL

!    CALL dbg_print('1: w', w, method_name, 3, in_subset=patch_2d%cells%owned)
!    CALL dbg_print('1: r', r, method_name, 3, in_subset=patch_2d%cells%owned)
!    write(0,*) "h_aux:", h_aux, " rn2(1):", rn2(1)

    !      write(*,*)  get_my_global_mpi_id(), ':', pad_nproma, nproma, 'r=', r(pad_nproma+1:nproma, no_of_blocks)
    !      write(*,*)  get_my_global_mpi_id(), ':', pad_nproma, nproma, 'v=', v(pad_nproma+1:nproma, no_of_blocks, 1)
    !      CALL p_barrier
    
    IF (rn2(1) == 0.0_wp) THEN ! already done
      !    print*,' gmres: rn2(1)=0.0 '
      ! #slo# 2010-06-15 trying to exit with niter=1 and residual=0.0
      niter = 0
      res(1) = ABS(rn2(1))

      stop_timer(timer_gmres,2)

      RETURN
    ENDIF
    
    ! 3) define the relative_tolerance: can be absolute or relative
    IF (absolute_tolerance > 0.0_wp) THEN
      tol = absolute_tolerance
    ELSE
      tol = relative_tolerance*rn2(1)
      absolute_tolerance = tol
    ENDIF
    tol2 = tol**2
    
    ! 4) Arnoldi loop
    arnoldi: DO i = 1, m-1
      ! 4.1) compute the next (i.e. i+1) Krylov vector
      !      sum_v = SUM(v(:,:,i))
      !       write(0,*) i, " gmres v before lhs:", sum_v, sum_w
      ! write(*,*)  get_my_global_mpi_id(), 'before lhs v(pad):', pad_nproma, nproma, k, 'v=', v(pad_nproma:nproma, no_of_blocks, k)
      w(:,:) = lhs( v(:,:,i), patch_3D, h_e, p_op_coeff )
      ! w(pad_nproma:nproma, no_of_blocks)   = 0.0_wp
      v(pad_nproma:nproma, no_of_blocks,i) = 0.0_wp
      
      !      sum_v = SUM(v(:,:,i))
      !      sum_w = SUM(w(:,:))
      !      write(0,*) i, " gmres sum v, w:", sum_v, sum_w
      
      ! 4.2) Gram-Schmidt orthogonalization
      
      
      IF (PRESENT(preconditioner)) CALL preconditioner(w(:,:),patch_3D,p_op_coeff,h_e)
      
      gs_orth: DO k = 1, i
        
!ICON_OMP_PARALLEL PRIVATE(myThreadNo)
!$   myThreadNo = OMP_GET_THREAD_NUM()
!ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          sum_aux(jb) = SUM(w(1:nproma,jb) * v(1:nproma,jb,k))
        ENDDO
!ICON_OMP_END_DO
        
        ! write(*,*)  get_my_global_mpi_id(), 'v(pad):', pad_nproma, nproma, k, 'v=', v(pad_nproma:nproma, no_of_blocks, k)
        ! write(*,*)  get_my_global_mpi_id(), 'w(pad):', pad_nproma, nproma, 'w=', w(pad_nproma:nproma, no_of_blocks)
        ! CALL p_barrier
        
        IF (mythreadno == 0) THEN
          h_aux = SUM(sum_aux(1:no_of_blocks))
          start_detail_timer(timer_gmres_p_sum,5)
          h_aux = p_sum(h_aux, my_mpi_work_communicator)
          stop_detail_timer(timer_gmres_p_sum,5)
          h(k,i) = h_aux
!ICON_OMP FLUSH(h_aux)
        ENDIF
!ICON_OMP BARRIER
!        CALL dbg_print('anroldi: w', w, method_name, 3, in_subset=patch_2d%cells%owned)
!        CALL dbg_print('anroldi: v', v, method_name, 3, in_subset=patch_2d%cells%owned)
!        write(0,*) "h_aux:", h_aux
        
!ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          w(1:nproma, jb) = w(1:nproma, jb) - h_aux * v(1:nproma, jb, k)
        ENDDO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
      ENDDO gs_orth
      
      !    write(*,*)  get_my_global_mpi_id(), ':', pad_nproma, nproma, 'w=', w(pad_nproma+1:nproma, no_of_blocks)
      !    CALL p_barrier
      !    w(pad_nproma:nproma, no_of_blocks) = 0.0_wp
      !    CALL sync_patch_array(SYNC_C, patch_2D, w(:,:) )
      
! !ICON_OMP_PARALLEL PRIVATE(rh, myThreadNo)
! !$   myThreadNo = OMP_GET_THREAD_NUM()
      ! 4.3) new element for h
! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, no_of_blocks
        sum_aux(jb) = SUM(w(1:nproma,jb) * w(1:nproma,jb))
      ENDDO
! !ICON_OMP_END_DO
      IF (mythreadno == 0) THEN
        h_aux = SUM(sum_aux(1:no_of_blocks))
        start_detail_timer(timer_gmres_p_sum,5)
        h_aux = p_sum(h_aux, my_mpi_work_communicator)
        stop_detail_timer(timer_gmres_p_sum,5)
        h_aux = SQRT(h_aux)
        h(i+1,i) = h_aux
        IF (h_aux < tol2) THEN
          done = .TRUE.
        ELSE
          done = .FALSE.
        ENDIF
      ENDIF
! !ICON_OMP FLUSH(h_aux, done)
! !ICON_OMP BARRIER
      !     write(0,*) i, " gmres  tol2, h_aux", tol2, h_aux
      
      
      ! 4.4) if w is independent from v, add v(:,:,:,i+1)
      IF (.NOT. done) THEN
        rh = 1.0_wp/h_aux
! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          v(1:nproma, jb, i+1) = w(1:nproma,jb) * rh
        ENDDO
! !ICON_OMP_END_DO NOWAIT
      ENDIF
      
! !ICON_OMP_END_PARALLEL
      
      ! 4.5) apply the rotation matrices
      rotation: DO k = 1, i-1
        hki  = h(k  ,i)
        hk1i = h(k+1,i)
        h(k  ,i) =  c(k)*hki + s(k)*hk1i
        h(k+1,i) = -s(k)*hki + c(k)*hk1i
      ENDDO rotation
      
      ! 4.6) compute the new (i.e. i) rotation
      den = SQRT(h(i,i)**2 + h(i+1,i)**2)
      c(i) = h(i  ,i) / den
      s(i) = h(i+1,i) / den
      
      ! 4.7) complete applying the rotation matrices
      h(i,i) =  c(i)*h(i,i) + s(i)*h(i+1,i)
      
      ! 4.8) compute new residual norm
      rn2(i+1) = -s(i)*rn2(i)
      rn2(i  ) =  c(i)*rn2(i)
      
      ! 4.9) check whether we are done
      ! write(0,*) i, " gmres  relative_tolerance bound, residual: ", tol, ABS(rn2(i+1))
      IF ( done .OR. (ABS(rn2(i+1)) < tol) ) EXIT arnoldi
    ENDDO arnoldi
    
    ! 5) check whether we are here because we have exceeded m
    !    Notice that if the arnoldi loop has been terminated with the
    !    exit instruction, we have i=m-1, while if the loop has been
    !    terminated by the loop counter we have i=m (i is incremented
    !    after finishing the do loop)
    IF (i == m) THEN
      maxiterex = .TRUE.
      ! only the first m-1 columns h(i,i) have been set in the arnoldi
      ! loop, but we have i=m. We thus decrement i.
      i = i-1
    ENDIF
    
    nkry = i
    niter = nkry+1
    
    ! 6) compute the coefficient of the Krylov expansion (back sub.)
    krylov: DO i = nkry, 1, -1
      ! coefficients are stored in c for convenience
      ci = rn2(i)
      DO k = nkry, i+1, -1
        ci = ci - h(i,k)*c(k)
      ENDDO
      c(i) = ci/h(i,i)
    ENDDO krylov
    
    ! 7) evaluate the Krylov expansion
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb, i) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      DO i = 1, nkry
        x(:,jb) = x(:,jb) + c(i)*v(:,jb,i)
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL    
    res(1:niter) = ABS(rn2(1:niter))

    stop_timer(timer_gmres,2)
    
  END SUBROUTINE ocean_restart_gmres
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !! Based on gmres_oce_old
  !! restart functionality and optimization L.Linardakis, MPIM, 2013
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE ocean_restart_gmres_singlePrecesicion( x,lhs, old_h, patch_3D, &
    & solverCoefficients, b,                      &
    & tolerance,abstol,m,maxiterex,niter,res)
    !
    ! !DESCRIPTION
    !  GMRES solver for linear systems. Notice that, in principle, r and b
    !  can be of any TYPE for which the following operations are defined:
    !   +, -, dot_product
    !
    REAL(sp), INTENT(inout) :: x(:,:)
    REAL(sp), INTENT(in) :: old_h(:,:)
    ! patch info needed for calculating lhs
    !TYPE(t_patch), INTENT(IN) :: patch_2D
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    ! index defining the "active" region of the arrays
    ! parameter used in calculating the lhs
!    REAL(sp), INTENT(in) :: coeff
    TYPE(t_solverCoeff_singlePrecision), INTENT(in):: solverCoefficients
    ! right-hand side: same shape as x
    REAL(sp), INTENT(inout) :: b(:,:) ! same size as x
    REAL(sp), INTENT(in) :: tolerance ! (relative or absolute) tolerance
    LOGICAL,  INTENT(in) :: abstol    ! .true. for absolute tolerance,
    ! .false. for relative tolerance
    INTEGER,  INTENT(in) :: m         ! maximum number of iterations
    LOGICAL,  INTENT(out) :: maxiterex ! true if reached m iterations
    INTEGER,  INTENT(out) :: niter    ! number of iterations (defined
    ! as the number of evaluation of
    ! the lhs which, due to the
    ! initial evaluation of the
    ! residual, is equal to the
    ! number of Arnoldi iterations +1)
    ! norms of the residual (convergence history); an argument of
    ! dimension at least m is required
    REAL(sp), INTENT(inout) :: res(:) ! (m)

    INTERFACE   ! left-hand-side: A*x
      FUNCTION lhs(x,old_h, patch_3D, sovlerCoeff) result(ax)
        USE mo_kind, ONLY: sp
        USE mo_model_domain, ONLY: t_patch, t_patch_3d
        USE mo_ocean_types, ONLY: t_solverCoeff_singlePrecision
        REAL(sp),    INTENT(inout) :: x(:,:)  ! inout for sync
        REAL(sp), INTENT(in) :: old_h(:,:)
        !TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
        TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
!        REAL(sp),    INTENT(in) :: coeff
        TYPE(t_solverCoeff_singlePrecision),INTENT(in)  :: sovlerCoeff
        REAL(sp) :: ax( SIZE(x,1) , SIZE(x,2) ) ! same as x
      END FUNCTION lhs
    END INTERFACE

    ! !LOCAL VARIABLES

    LOGICAL :: done
    INTEGER ::     &
      & i,           & ! index for the Arnoldi loop
      & k,           & ! index for the Gram Schmidt orthogonalization
      & nkry           ! dimension of the Krylov space
    REAL(sp) ::    &
      & tol,         & ! effective tolerance (absolute o relative)
      & tol2,        &
      & r(SIZE(x,1),SIZE(x,2)),  & ! residual
      & rn2(m),      &             ! two-norm of the residual
      & v(SIZE(x,1),SIZE(x,2),m),& ! Krylov basis
      & h(m,m),      &             ! Hessemberg matrix
      & hki, hk1i,   &
      & c(m), s(m),  &             ! rotation matrices
      & den, ci

    REAL(sp) ::  w(SIZE(x,1),SIZE(x,2)) ! new Krylov vector

    REAL(sp) :: rrn2, h_aux, rh
    INTEGER :: jb, jk, nlen

    INTEGER :: no_of_blocks, pad_nproma, end_nproma
    INTEGER :: my_mpi_work_communicator

    REAL(sp) :: sum_aux(patch_3D%p_patch_2d(1)%cells%in_domain%end_block)

    INTEGER :: mythreadno
    TYPE(t_patch), POINTER :: patch_2d
    CHARACTER(len=*), PARAMETER :: method_name='ocean_restart_gmres'

    !$ INTEGER OMP_GET_THREAD_NUM
    !-------------------------------------------------------------------------
    !  write(0,*) "--------------- gmres --------------------------"

    patch_2d =>  patch_3D%p_patch_2d(1)
    my_mpi_work_communicator = get_my_mpi_work_communicator()
    ! 0) set module variables and initialize maxiterex
    no_of_blocks    = patch_2d%cells%in_domain%end_block
    end_nproma      = patch_2d%cells%in_domain%end_index
    pad_nproma      = patch_2d%cells%in_domain%end_index + 1

    maxiterex = .FALSE.

    v(:,:,:)  = 0.0_sp
    r(:,:)    = 0.0_sp
    b(pad_nproma:nproma, no_of_blocks) = 0.0_sp

    ! 1) compute the preconditioned residual
    ! IF (PRESENT(preconditioner)) CALL preconditioner(x(:,:),patch_3D,p_op_coeff,h_e)
    ! IF (PRESENT(preconditioner)) CALL preconditioner(b(:,:),patch_3D,p_op_coeff,h_e)


    w(:, :) = lhs(x(:,:),old_h, patch_3D, solverCoefficients)

    x(pad_nproma:nproma, no_of_blocks) = 0.0_sp

    !    write(0,*) "-----------------------------------"
    !    sum_x(1) = SUM(x(:,:))
    !    sum_w    = SUM(w(:,:))
    !    write(0,*) "gmres sum x(0:1), w:", sum_x(0:1), sum_w

    start_timer(timer_gmres,2)

    mythreadno = 0
! !ICON_OMP_PARALLEL PRIVATE(rrn2, myThreadNo)
!    !$ myThreadNo = OMP_GET_THREAD_NUM()

! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      r(1:nproma,jb) = b(1:nproma,jb) - w(1:nproma,jb)
    ENDDO
! !ICON_OMP_END_DO

    ! IF (PRESENT(preconditioner)) CALL preconditioner(r(:,:),patch_3D,p_op_coeff,h_e)

! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      sum_aux(jb) = SUM(r(1:nproma,jb) * r(1:nproma,jb))
    ENDDO
    ! sum_aux(no_of_blocks) = SUM(r(1:end_nproma,no_of_blocks) * r(1:end_nproma,no_of_blocks))
! !ICON_OMP_END_DO

    IF (mythreadno == 0) THEN
      h_aux  = SUM(sum_aux(1:no_of_blocks))
      start_detail_timer(timer_gmres_p_sum,5)
      h_aux = p_sum(h_aux, my_mpi_work_communicator)
      stop_detail_timer(timer_gmres_p_sum,5)
      rn2(1) = SQRT(h_aux)
      !       rn2(1) = SQRT(p_sum(h_aux, my_mpi_work_communicator))
      ! !ICON_OMP FLUSH(rn2(1))
    ENDIF
! !ICON_OMP BARRIER


    ! 2) compute the first vector of the Krylov space
    IF (rn2(1) /= 0.0_sp) THEN
      rrn2 = REAL(1.0_wp/rn2(1),sp)
! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, no_of_blocks
        v(1:nproma,jb,1) = r(1:nproma,jb) * rrn2
      ENDDO
! !ICON_OMP_END_DO NOWAIT
    ENDIF
! !ICON_OMP_END_PARALLEL

!    CALL dbg_print('1: w', w, method_name, 3, in_subset=patch_2d%cells%owned)
!    CALL dbg_print('1: r', r, method_name, 3, in_subset=patch_2d%cells%owned)
!    write(0,*) "h_aux:", h_aux, " rn2(1):", rn2(1)

    !      write(*,*)  get_my_global_mpi_id(), ':', pad_nproma, nproma, 'r=', r(pad_nproma+1:nproma, no_of_blocks)
    !      write(*,*)  get_my_global_mpi_id(), ':', pad_nproma, nproma, 'v=', v(pad_nproma+1:nproma, no_of_blocks, 1)
    !      CALL p_barrier

    IF (rn2(1) == 0.0_sp) THEN ! already done
      !    print*,' gmres: rn2(1)=0.0 '
      ! #slo# 2010-06-15 trying to exit with niter=1 and residual=0.0
      niter = 0
      res(1) = ABS(rn2(1))

      IF (ltimer) CALL timer_stop(timer_gmres)

      RETURN
    ENDIF

    ! 3) define the tolerance: can be absolute or relative
    IF (abstol) THEN
      tol = tolerance
    ELSE
      tol = tolerance*rn2(1)
    ENDIF
    tol2 = tol**2

    ! 4) Arnoldi loop
    arnoldi: DO i = 1, m-1
      ! 4.1) compute the next (i.e. i+1) Krylov vector
      !      sum_v = SUM(v(:,:,i))
      !       write(0,*) i, " gmres v before lhs:", sum_v, sum_w
      ! write(*,*)  get_my_global_mpi_id(), 'before lhs v(pad):', pad_nproma, nproma, k, 'v=', v(pad_nproma:nproma, no_of_blocks, k)
      w(:,:) = lhs( v(:,:,i), old_h, patch_3D, solverCoefficients )
      ! w(pad_nproma:nproma, no_of_blocks)   = 0.0_sp
      v(pad_nproma:nproma, no_of_blocks,i) = 0.0_sp

      !      sum_v = SUM(v(:,:,i))
      !      sum_w = SUM(w(:,:))
      !      write(0,*) i, " gmres sum v, w:", sum_v, sum_w

      ! 4.2) Gram-Schmidt orthogonalization


      ! IF (PRESENT(preconditioner)) CALL preconditioner(w(:,:),patch_3D,p_op_coeff,h_e)

      gs_orth: DO k = 1, i

! !ICON_OMP_PARALLEL PRIVATE(myThreadNo)
! !$   myThreadNo = OMP_GET_THREAD_NUM()
! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          sum_aux(jb) = SUM(w(1:nproma,jb) * v(1:nproma,jb,k))
        ENDDO
! !ICON_OMP_END_DO

        ! write(*,*)  get_my_global_mpi_id(), 'v(pad):', pad_nproma, nproma, k, 'v=', v(pad_nproma:nproma, no_of_blocks, k)
        ! write(*,*)  get_my_global_mpi_id(), 'w(pad):', pad_nproma, nproma, 'w=', w(pad_nproma:nproma, no_of_blocks)
        ! CALL p_barrier

        IF (mythreadno == 0) THEN
          h_aux = SUM(sum_aux(1:no_of_blocks))
          start_detail_timer(timer_gmres_p_sum,5)
          h_aux = p_sum(h_aux, my_mpi_work_communicator)
          stop_detail_timer(timer_gmres_p_sum,5)
          h(k,i) = h_aux
! !ICON_OMP FLUSH(h_aux)
        ENDIF
! !ICON_OMP BARRIER
!        CALL dbg_print('anroldi: w', w, method_name, 3, in_subset=patch_2d%cells%owned)
!        CALL dbg_print('anroldi: v', v, method_name, 3, in_subset=patch_2d%cells%owned)
!        write(0,*) "h_aux:", h_aux

! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          w(1:nproma, jb) = w(1:nproma, jb) - h_aux * v(1:nproma, jb, k)
        ENDDO
! !ICON_OMP_END_DO NOWAIT
! !ICON_OMP_END_PARALLEL
      ENDDO gs_orth

      !    write(*,*)  get_my_global_mpi_id(), ':', pad_nproma, nproma, 'w=', w(pad_nproma+1:nproma, no_of_blocks)
      !    CALL p_barrier
      !    w(pad_nproma:nproma, no_of_blocks) = 0.0_wp
      !    CALL sync_patch_array(SYNC_C, patch_2D, w(:,:) )

!ICON_OMP_PARALLEL PRIVATE(myThreadNo, rh)
!$   myThreadNo = OMP_GET_THREAD_NUM()
      ! 4.3) new element for h
!ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, no_of_blocks
        sum_aux(jb) = SUM(w(1:nproma,jb) * w(1:nproma,jb))
      ENDDO
!ICON_OMP_END_DO
      IF (mythreadno == 0) THEN
        h_aux = SUM(sum_aux(1:no_of_blocks))
        start_detail_timer(timer_gmres_p_sum,5)
        h_aux = p_sum(h_aux, my_mpi_work_communicator)
        stop_detail_timer(timer_gmres_p_sum,5)
        h_aux = SQRT(h_aux)
        h(i+1,i) = h_aux
        IF (h_aux < tol2) THEN
          done = .TRUE.
        ELSE
          done = .FALSE.
        ENDIF
      ENDIF
!ICON_OMP FLUSH(h_aux, done)
!ICON_OMP BARRIER
      !     write(0,*) i, " gmres  tol2, h_aux", tol2, h_aux

      ! 4.4) if w is independent from v, add v(:,:,:,i+1)
      IF (.NOT. done) THEN
        rh = REAL(1.0_wp/h_aux,sp)
!ICON_OMP_DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          v(1:nproma, jb, i+1) = w(1:nproma,jb) * rh
        ENDDO
!ICON_OMP_END_DO NOWAIT
      ENDIF

!ICON_OMP_END_PARALLEL

      ! 4.5) apply the rotation matrices
      rotation: DO k = 1, i-1
        hki  = h(k  ,i)
        hk1i = h(k+1,i)
        h(k  ,i) =  c(k)*hki + s(k)*hk1i
        h(k+1,i) = -s(k)*hki + c(k)*hk1i
      ENDDO rotation

      ! 4.6) compute the new (i.e. i) rotation
      den = SQRT(h(i,i)**2 + h(i+1,i)**2)
      c(i) = h(i  ,i) / den
      s(i) = h(i+1,i) / den

      ! 4.7) complete applying the rotation matrices
      h(i,i) =  c(i)*h(i,i) + s(i)*h(i+1,i)

      ! 4.8) compute new residual norm
      rn2(i+1) = -s(i)*rn2(i)
      rn2(i  ) =  c(i)*rn2(i)

      ! 4.9) check whether we are done
      ! write(0,*) i, " gmres  tolerance bound, residual: ", tol, ABS(rn2(i+1))
      IF ( done .OR. (ABS(rn2(i+1)) < tol) ) EXIT arnoldi
    ENDDO arnoldi

    ! 5) check whether we are here because we have exceeded m
    !    Notice that if the arnoldi loop has been terminated with the
    !    exit instruction, we have i=m-1, while if the loop has been
    !    terminated by the loop counter we have i=m (i is incremented
    !    after finishing the do loop)
    IF (i == m) THEN
      maxiterex = .TRUE.
      ! only the first m-1 columns h(i,i) have been set in the arnoldi
      ! loop, but we have i=m. We thus decrement i.
      i = i-1
    ENDIF

    nkry = i
    niter = nkry+1

    ! 6) compute the coefficient of the Krylov expansion (back sub.)
    krylov: DO i = nkry, 1, -1
      ! coefficients are stored in c for convenience
      ci = rn2(i)
      DO k = nkry, i+1, -1
        ci = ci - h(i,k)*c(k)
      ENDDO
      c(i) = ci/h(i,i)
    ENDDO krylov

    ! 7) evaluate the Krylov expansion
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb, i) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      DO i = 1, nkry
        x(:,jb) = x(:,jb) + c(i)*v(:,jb,i)
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    res(1:niter) = ABS(rn2(1:niter))

    stop_timer(timer_gmres,2)

  END SUBROUTINE ocean_restart_gmres_singlePrecesicion
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !! Based on gmres_oce_old
  !! restart functionality and optimization L.Linardakis, MPIM, 2013
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE ocean_restart_gmres_e2e( x,lhs,        &
    & patch_3D,  &
    & p_op_coeff, level,b,                            &
    & absolute_tolerance,  relative_tolerance,  &
    & m, maxiterex,niter,res,                   &
    & preconditioner)
    !
    ! !DESCRIPTION
    !  GMRES solver for linear systems. Notice that, in principle, r and b
    !  can be of any TYPE for which the following operations are defined:
    !   +, -, dot_product
    !
    REAL(wp), INTENT(inout) :: x(:,:)
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    TYPE(t_operator_coeff), INTENT(in)     :: p_op_coeff
    INTEGER                                :: level
    REAL(wp), INTENT(inout)                :: b(:,:) ! same size as x
    REAL(wp), INTENT(inout)                :: absolute_tolerance
    REAL(wp), INTENT(in)                   :: relative_tolerance ! (relative or absolute) relative_tolerance
    INTEGER,  INTENT(in)                   :: m         ! maximum number of iterations
    LOGICAL,  INTENT(out)                  :: maxiterex ! true if reached m iterations
    INTEGER,  INTENT(out)                  :: niter    ! number of iterations (defined
    ! as the number of evaluation of the lhs which, due to the
    ! initial evaluation of the residual, is equal to the
    ! number of Arnoldi iterations +1) norms of the residual (convergence history); an argument of
    ! dimension at least m is required
    REAL(wp), INTENT(inout)                   :: res(:) ! (m)
    
    INTERFACE   ! left-hand-side: A*x
      FUNCTION lhs(x,patch_3D, p_op_coeff, level) result(ax)
        USE mo_kind, ONLY: wp
        USE mo_model_domain, ONLY: t_patch, t_patch_3d
        USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
        REAL(wp),    INTENT(inout)             :: x(:,:)  ! inout for sync
        TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
        TYPE(t_operator_coeff),INTENT(in)      :: p_op_coeff
        INTEGER                                :: level
        REAL(wp) :: ax( SIZE(x,1) , SIZE(x,2) ) ! same as x
      END FUNCTION lhs
    END INTERFACE
    
    INTERFACE   ! preconditioner
      SUBROUTINE preconditioner(r, patch_3D, p_op_coeff,h_e)
        USE mo_kind, ONLY: wp
        USE mo_model_domain, ONLY: t_patch, t_patch_3d
        USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
        REAL(wp), INTENT(inout)           :: r(:,:)
        !TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
        TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
        TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
        REAL(wp),    INTENT(in)           :: h_e(:,:)
      END SUBROUTINE preconditioner
    END INTERFACE
    
    OPTIONAL :: preconditioner
    
    ! !LOCAL VARIABLES
    
    LOGICAL :: done
    INTEGER ::     &
      & i,           & ! index for the Arnoldi loop
      & k,           & ! index for the Gram Schmidt orthogonalization
      & nkry           ! dimension of the Krylov space
    REAL(wp) ::    &
      & tol,         & ! effective relative_tolerance (absolute o relative)
      & tol2,        &
      & r(SIZE(x,1),SIZE(x,2)),  & ! residual
      & rn2(m),      &             ! two-norm of the residual
      & v(SIZE(x,1),SIZE(x,2),m),& ! Krylov basis
      & h(m,m),      &             ! Hessemberg matrix
      & hki, hk1i,   &
      & c(m), s(m),  &             ! rotation matrices
      & den, ci
    
    ! REAL(wp), POINTER :: w(:,:)
    REAL(wp) ::  w(SIZE(x,1),SIZE(x,2)) ! new Krylov vector
    
    REAL(wp) :: rrn2, h_aux, rh
    INTEGER :: jb, jk, nlen
    
    INTEGER :: no_of_blocks, pad_nproma, end_nproma
    INTEGER :: my_mpi_work_communicator
    
    REAL(wp) :: sum_aux(patch_3D%p_patch_2d(1)%edges%in_domain%end_block)
    INTEGER :: mythreadno
    TYPE(t_patch), POINTER :: patch_2d
    CHARACTER(len=*), PARAMETER :: method_name='ocean_restart_gmres_e2e'

    !$ INTEGER OMP_GET_THREAD_NUM
    !-------------------------------------------------------------------------
    !  write(0,*) "--------------- gmres --------------------------"
    
    patch_2d =>  patch_3D%p_patch_2d(1)
    my_mpi_work_communicator = get_my_mpi_work_communicator()
    ! 0) set module variables and initialize maxiterex
    no_of_blocks    = patch_2d%edges%in_domain%end_block
    end_nproma      = patch_2d%edges%in_domain%end_index
    pad_nproma      = patch_2d%edges%in_domain%end_index + 1
    
    maxiterex = .FALSE.
    
    v(:,:,:)  = 0.0_wp
    r(:,:)    = 0.0_wp
    b(pad_nproma:nproma, no_of_blocks) = 0.0_wp
 !   w(:, :) = 0.0_wp
 !   b(:, :) = 0.0_wp
 !   x(:, :) = 0.0_wp
    
    ! 1) compute the preconditioned residual
    !IF (PRESENT(preconditioner)) CALL preconditioner(x(:,:),patch_3D,p_op_coeff,h_e)
    !IF (PRESENT(preconditioner)) CALL preconditioner(b(:,:),patch_3D,p_op_coeff,h_e)
    
    
    !w(:, :) = lhs(x(:,:),old_h, patch_3D, h_e, thickness_c, p_op_coeff)
    w(:, :) =lhs(x,patch_3D, p_op_coeff, level)

    ! w(pad_nproma:nproma, no_of_blocks) = 0.0_wp
    x(pad_nproma:nproma, no_of_blocks) = 0.0_wp
    
    !    write(0,*) "-----------------------------------"
    !    sum_x(1) = SUM(x(:,:))
    !    sum_w    = SUM(w(:,:))
    !    write(0,*) "gmres sum x(0:1), w:", sum_x(0:1), sum_w
    
    IF (ltimer) CALL timer_start(timer_gmres)
    
    mythreadno = 0
! !ICON_OMP_PARALLEL PRIVATE(rrn2, myThreadNo)
!    !$ myThreadNo = OMP_GET_THREAD_NUM()

! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      r(1:nproma,jb) = b(1:nproma,jb) - w(1:nproma,jb)
    ENDDO
! !ICON_OMP_END_DO
    
!     IF (PRESENT(preconditioner)) CALL preconditioner(r(:,:),patch_3D,p_op_coeff,h_e)
    
! !ICON_OMP_DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      sum_aux(jb) = SUM(r(1:nproma,jb) * r(1:nproma,jb))
    ENDDO
    ! sum_aux(no_of_blocks) = SUM(r(1:end_nproma,no_of_blocks) * r(1:end_nproma,no_of_blocks))
! !ICON_OMP_END_DO
    
    IF (mythreadno == 0) THEN
      h_aux  = SUM(sum_aux(1:no_of_blocks))
      IF (activate_sync_timers) CALL timer_start(timer_gmres_p_sum)
      h_aux = p_sum(h_aux, my_mpi_work_communicator)
      IF (activate_sync_timers) CALL timer_stop(timer_gmres_p_sum)
      rn2(1) = SQRT(h_aux)
      !       rn2(1) = SQRT(p_sum(h_aux, my_mpi_work_communicator))
      ! !ICON_OMP FLUSH(rn2(1))
    ENDIF
! !ICON_OMP BARRIER
    
    
    ! 2) compute the first vector of the Krylov space
    IF (rn2(1) /= 0.0_wp) THEN
      rrn2 = 1.0_wp/rn2(1)
! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, no_of_blocks
        v(1:nproma,jb,1) = r(1:nproma,jb) * rrn2
      ENDDO
! !ICON_OMP_END_DO NOWAIT
    ENDIF
! !ICON_OMP_END_PARALLEL

!    CALL dbg_print('1: w', w, method_name, 3, in_subset=patch_2d%cells%owned)
!    CALL dbg_print('1: r', r, method_name, 3, in_subset=patch_2d%cells%owned)
!    write(0,*) "h_aux:", h_aux, " rn2(1):", rn2(1)

    !      write(*,*)  get_my_global_mpi_id(), ':', pad_nproma, nproma, 'r=', r(pad_nproma+1:nproma, no_of_blocks)
    !      write(*,*)  get_my_global_mpi_id(), ':', pad_nproma, nproma, 'v=', v(pad_nproma+1:nproma, no_of_blocks, 1)
    !      CALL p_barrier
    
    IF (rn2(1) == 0.0_wp) THEN ! already done
      !    print*,' gmres: rn2(1)=0.0 '
      ! #slo# 2010-06-15 trying to exit with niter=1 and residual=0.0
      niter = 0
      res(1) = ABS(rn2(1))

      IF (ltimer) CALL timer_stop(timer_gmres)

      RETURN
    ENDIF
    
    ! 3) define the relative_tolerance: can be absolute or relative
    IF (absolute_tolerance > 0.0_wp) THEN
      tol = absolute_tolerance
    ELSE
      tol = relative_tolerance*rn2(1)
      absolute_tolerance = tol
    ENDIF
    tol2 = tol**2
    
    ! 4) Arnoldi loop
    arnoldi: DO i = 1, m-1
      ! 4.1) compute the next (i.e. i+1) Krylov vector
      !      sum_v = SUM(v(:,:,i))
      !       write(0,*) i, " gmres v before lhs:", sum_v, sum_w
      ! write(*,*)  get_my_global_mpi_id(), 'before lhs v(pad):', pad_nproma, nproma, k, 'v=', v(pad_nproma:nproma, no_of_blocks, k)
      !w(:,:) = lhs( v(:,:,i), old_h, patch_3D, h_e, thickness_c, p_op_coeff )
      w(:, :) =lhs(v(:,:,i),patch_3D, p_op_coeff, level)

      ! w(pad_nproma:nproma, no_of_blocks)   = 0.0_wp
      v(pad_nproma:nproma, no_of_blocks,i) = 0.0_wp
      
      !      sum_v = SUM(v(:,:,i))
      !      sum_w = SUM(w(:,:))
      !      write(0,*) i, " gmres sum v, w:", sum_v, sum_w
      
      ! 4.2) Gram-Schmidt orthogonalization
      
      
      !IF (PRESENT(preconditioner)) CALL preconditioner(w(:,:),patch_3D,p_op_coeff,h_e)
      
      gs_orth: DO k = 1, i
        
!ICON_OMP_PARALLEL PRIVATE(myThreadNo)
       !$   myThreadNo = OMP_GET_THREAD_NUM()
!ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          sum_aux(jb) = SUM(w(1:nproma,jb) * v(1:nproma,jb,k))
        ENDDO
!ICON_OMP_END_DO
        
        ! write(*,*)  get_my_global_mpi_id(), 'v(pad):', pad_nproma, nproma, k, 'v=', v(pad_nproma:nproma, no_of_blocks, k)
        ! write(*,*)  get_my_global_mpi_id(), 'w(pad):', pad_nproma, nproma, 'w=', w(pad_nproma:nproma, no_of_blocks)
        ! CALL p_barrier
        
        IF (mythreadno == 0) THEN
          h_aux = SUM(sum_aux(1:no_of_blocks))
          IF (activate_sync_timers) CALL timer_start(timer_gmres_p_sum)
          h_aux = p_sum(h_aux, my_mpi_work_communicator)
          IF (activate_sync_timers) CALL timer_stop(timer_gmres_p_sum)
          h(k,i) = h_aux
!ICON_OMP FLUSH(h_aux)
        ENDIF
!ICON_OMP BARRIER
!        CALL dbg_print('anroldi: w', w, method_name, 3, in_subset=patch_2d%cells%owned)
!        CALL dbg_print('anroldi: v', v, method_name, 3, in_subset=patch_2d%cells%owned)
!        write(0,*) "h_aux:", h_aux
        
!ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          w(1:nproma, jb) = w(1:nproma, jb) - h_aux * v(1:nproma, jb, k)
        ENDDO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
      ENDDO gs_orth
      
      !    write(*,*)  get_my_global_mpi_id(), ':', pad_nproma, nproma, 'w=', w(pad_nproma+1:nproma, no_of_blocks)
      !    CALL p_barrier
      !    w(pad_nproma:nproma, no_of_blocks) = 0.0_wp
      !    CALL sync_patch_array(SYNC_C, patch_2D, w(:,:) )
      
! !ICON_OMP_PARALLEL PRIVATE(rh, myThreadNo)
! !$   myThreadNo = OMP_GET_THREAD_NUM()
      ! 4.3) new element for h
! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, no_of_blocks
        sum_aux(jb) = SUM(w(1:nproma,jb) * w(1:nproma,jb))
      ENDDO
! !ICON_OMP_END_DO
      IF (mythreadno == 0) THEN
        h_aux = SUM(sum_aux(1:no_of_blocks))
        IF (activate_sync_timers) CALL timer_start(timer_gmres_p_sum)
        h_aux = p_sum(h_aux, my_mpi_work_communicator)
        IF (activate_sync_timers) CALL timer_stop(timer_gmres_p_sum)
        h_aux = SQRT(h_aux)
        h(i+1,i) = h_aux
        IF (h_aux < tol2) THEN
          done = .TRUE.
        ELSE
          done = .FALSE.
        ENDIF
      ENDIF
! !ICON_OMP FLUSH(h_aux, done)
! !ICON_OMP BARRIER
      !     write(0,*) i, " gmres  tol2, h_aux", tol2, h_aux
      
      
      ! 4.4) if w is independent from v, add v(:,:,:,i+1)
      IF (.NOT. done) THEN
        rh = 1.0_wp/h_aux
! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          v(1:nproma, jb, i+1) = w(1:nproma,jb) * rh
        ENDDO
! !ICON_OMP_END_DO NOWAIT
      ENDIF
      
! !ICON_OMP_END_PARALLEL
      
      ! 4.5) apply the rotation matrices
      rotation: DO k = 1, i-1
        hki  = h(k  ,i)
        hk1i = h(k+1,i)
        h(k  ,i) =  c(k)*hki + s(k)*hk1i
        h(k+1,i) = -s(k)*hki + c(k)*hk1i
      ENDDO rotation
      
      ! 4.6) compute the new (i.e. i) rotation
      den = SQRT(h(i,i)**2 + h(i+1,i)**2)
      c(i) = h(i  ,i) / den
      s(i) = h(i+1,i) / den
      
      ! 4.7) complete applying the rotation matrices
      h(i,i) =  c(i)*h(i,i) + s(i)*h(i+1,i)
      
      ! 4.8) compute new residual norm
      rn2(i+1) = -s(i)*rn2(i)
      rn2(i  ) =  c(i)*rn2(i)
      
      ! 4.9) check whether we are done
      ! write(0,*) i, " gmres  relative_tolerance bound, residual: ", tol, ABS(rn2(i+1))
      IF ( done .OR. (ABS(rn2(i+1)) < tol) ) EXIT arnoldi
    ENDDO arnoldi
    
    ! 5) check whether we are here because we have exceeded m
    !    Notice that if the arnoldi loop has been terminated with the
    !    exit instruction, we have i=m-1, while if the loop has been
    !    terminated by the loop counter we have i=m (i is incremented
    !    after finishing the do loop)
    IF (i == m) THEN
      maxiterex = .TRUE.
      ! only the first m-1 columns h(i,i) have been set in the arnoldi
      ! loop, but we have i=m. We thus decrement i.
      i = i-1
    ENDIF
    
    nkry = i
    niter = nkry+1
    
    ! 6) compute the coefficient of the Krylov expansion (back sub.)
    krylov: DO i = nkry, 1, -1
      ! coefficients are stored in c for convenience
      ci = rn2(i)
      DO k = nkry, i+1, -1
        ci = ci - h(i,k)*c(k)
      ENDDO
      c(i) = ci/h(i,i)
    ENDDO krylov
    
    ! 7) evaluate the Krylov expansion
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb, i) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      DO i = 1, nkry
        x(:,jb) = x(:,jb) + c(i)*v(:,jb,i)
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    IF (ltimer) CALL timer_stop(timer_gmres)
    
    res(1:niter) = ABS(rn2(1:niter))
    
  END SUBROUTINE ocean_restart_gmres_e2e
  !-------------------------------------------------------------------------
  
  
  
  
  
  
  
  !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  Original version from http://www.netlib.org/templates/templates.pdf
  !!  F90 rewriting by Marco Restelli.
  !!  Adapted for use in ICOHDC by Hui Wan, MPI-M (2007-12-15)
  !!  Included blocking by Marco Restelli (2008-09-22)
  !!  MPI Parallelization by Rainer Johanni (2009-11-11)
  !! interpolation state removed from parameter list and adapted to 2D arrays, Peter Korn (2010-04)
  !! @par
  !! inital guess, overwritten with the solution
  !!
  !-----------------------------------------------------------------
  SUBROUTINE gmres_oce_old( x,lhs,h_e, old_h, patch_3D, &
    & p_op_coeff, b,                      &
    & tolerance,abstol,m,maxiterex,niter,res,    &
    & preconditioner)
    !
    ! !DESCRIPTION
    !  GMRES solver for linear systems. Notice that, in principle, r and b
    !  can be of any TYPE for which the following operations are defined:
    !   +, -, dot_product
    !
    REAL(wp), INTENT(inout) :: x(:,:)
    REAL(wp), INTENT(in) :: h_e(:,:)
    REAL(wp), INTENT(in) :: old_h(:,:)
    ! patch info needed for calculating lhs
    !TYPE(t_patch), INTENT(IN) :: patch_2D
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    ! index defining the "active" region of the arrays
    ! parameter used in calculating the lhs
!    REAL(wp), INTENT(in) :: coeff
    TYPE(t_operator_coeff), INTENT(in):: p_op_coeff
    ! right-hand side: same shape as x
    REAL(wp), INTENT(inout) :: b(:,:) ! same size as x
    REAL(wp), INTENT(in) :: tolerance ! (relative or absolute) tolerance
    LOGICAL,  INTENT(in) :: abstol    ! .true. for absolute tolerance,
    ! .false. for relative tolerance
    INTEGER,  INTENT(in) :: m         ! maximum number of iterations
    LOGICAL,  INTENT(out) :: maxiterex ! true if reached m iterations
    INTEGER,  INTENT(out) :: niter    ! number of iterations (defined
    ! as the number of evaluation of
    ! the lhs which, due to the
    ! initial evaluation of the
    ! residual, is equal to the
    ! number of Arnoldi iterations +1)
    ! norms of the residual (convergence history); an argument of
    ! dimension at least m is required
    REAL(wp), INTENT(inout) :: res(:) ! (m)
    
    INTERFACE   ! left-hand-side: A*x
      FUNCTION lhs(x, patch_3D, h_e, p_op_coeff) result(ax)
        USE mo_kind, ONLY: wp
        USE mo_model_domain, ONLY: t_patch, t_patch_3d
        USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
        REAL(wp),    INTENT(inout) :: x(:,:)  ! inout for sync
!         REAL(wp), INTENT(in) :: old_h(:,:)
        !TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
        TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
!        REAL(wp),    INTENT(in) :: coeff
        REAL(wp),    INTENT(in) :: h_e(:,:)
        TYPE(t_operator_coeff),INTENT(in)  :: p_op_coeff
        !     REAL(wp), POINTER :: ax( :, : ) ! same as x
        REAL(wp) :: ax( SIZE(x,1) , SIZE(x,2) ) ! same as x
      ENDFUNCTION lhs
    END INTERFACE
    
    INTERFACE   ! preconditioner
      SUBROUTINE preconditioner(r, patch_3D, p_op_coeff,h_e)
        USE mo_kind, ONLY: wp
        USE mo_model_domain, ONLY: t_patch, t_patch_3d
        USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
        REAL(wp), INTENT(inout)           :: r(:,:)
        !TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
        TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
        TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
        REAL(wp),    INTENT(in)           :: h_e(:,:)
      END SUBROUTINE preconditioner
    END INTERFACE
    
    OPTIONAL :: preconditioner
    
    
    LOGICAL :: done
    INTEGER ::     &
      & i,           & ! index for the Arnoldi loop
      & k,           & ! index for the Gram Schmidt orthogonalization
      & nkry           ! dimension of the Krylov space
    REAL(wp) ::    &
      & tol,         & ! effective tolerance (absolute o relative)
      & tol2,        &
      & r(SIZE(x,1),SIZE(x,2)),  & ! residual
      & rn2(m),      &             ! two-norm of the residual
      & v(SIZE(x,1),SIZE(x,2),m),& ! Krylov basis
      & h(m,m),      &             ! Hessemberg matrix
      & hki, hk1i,   &
      & c(m), s(m),  &             ! rotation matrices
      & den, ci
    
    ! REAL(wp), POINTER :: w(:,:)
    REAL(wp) ::  w(SIZE(x,1),SIZE(x,2)) ! new Krylov vector
    
    REAL(wp) :: rrn2, rn2_aux, h_aux, rh
    INTEGER :: jb, jk, nlen
    
    INTEGER :: no_of_blocks, end_nproma
    
    REAL(wp) :: z(SIZE(x,1),SIZE(x,2)) ! needed for global sums in p_test_run
    
    INTEGER :: mythreadno
    TYPE(t_patch), POINTER :: patch_2d
    CHARACTER(len=*), PARAMETER :: method_name='gmres_oce_old'
    !$ INTEGER OMP_GET_THREAD_NUM
    !-------------------------------------------------------------------------
    patch_2d =>  patch_3D%p_patch_2d(1)
    
    ! 0) set module variables and initialize maxiterex
    !>
    !!
    no_of_blocks = patch_2d%cells%in_domain%end_block
    end_nproma   = patch_2d%cells%in_domain%end_index
    
    maxiterex = .FALSE.
    z(:,:)    = 0._wp 
    v(:,:,:)  = 0.0_wp
    r(:,:)    = 0.0_wp
    
    ! 1) compute the preconditioned residual
    IF (PRESENT(preconditioner)) CALL preconditioner(x(:,:),patch_3D,p_op_coeff,h_e)
    IF (PRESENT(preconditioner)) CALL preconditioner(b(:,:),patch_3D,p_op_coeff,h_e)
    w(:,:) = lhs(x(:,:), patch_3D, h_e, p_op_coeff)
    
    IF (ltimer) CALL timer_start(timer_gmres)
    
    mythreadno = 0
! !ICON_OMP_PARALLEL PRIVATE(rn2_aux, rrn2, myThreadNo)
!     !$ myThreadNo = OMP_GET_THREAD_NUM()
!     
!     
! !ICON_OMP_DO PRIVATE(jb,nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      IF (jb /= no_of_blocks) THEN
        nlen = nproma
      ELSE
        nlen = end_nproma
        r(nlen+1:nproma,jb) = 0.0_wp
      ENDIF
      r(1:nlen,jb) = b(1:nlen,jb)-w(1:nlen,jb)
    ENDDO
! !ICON_OMP_END_DO
    
    IF (PRESENT(preconditioner)) CALL preconditioner(r(:,:),patch_3D,p_op_coeff,h_e)
    
! !ICON_OMP_DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      IF (jb /= no_of_blocks) THEN
        z(:,jb) = r(:,jb)*r(:,jb)
      ELSE
        z(1:end_nproma,jb) = r(1:end_nproma,jb)*r(1:end_nproma,jb)
      ENDIF
      WHERE(.NOT.patch_2d%cells%decomp_info%owner_mask(:,jb)) z(:,jb) = 0.0_wp
    ENDDO
! !ICON_OMP_END_DO
    
    rn2_aux = SQRT(omp_global_sum_array(z(:,1:no_of_blocks)))
    
    IF (mythreadno == 0) rn2(1) = rn2_aux
    
    ! 2) compute the first vector of the Krylov space
    IF (rn2_aux /= 0.0_wp) THEN
      rrn2 = 1.0_wp/rn2_aux
! !ICON_OMP_DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, no_of_blocks
        v(:,jb,1) = r(:,jb)*rrn2
      ENDDO
! !ICON_OMP_END_DO NOWAIT
    ENDIF
    
! !ICON_OMP_END_PARALLEL
!    CALL dbg_print('1: w', w, method_name, 3, in_subset=patch_2d%cells%owned)
!    CALL dbg_print('1: r', r, method_name, 3, in_subset=patch_2d%cells%owned)
!    write(0,*) " rn2(1):", rn2(1)
    
    IF (rn2(1) == 0.0_wp) THEN ! already done
      !    print*,' gmres: rn2(1)=0.0 '
      ! #slo# 2010-06-15 trying to exit with niter=1 and residual=0.0
      niter = 0
      !     niter = 1
      res(1) = ABS(rn2(1))
      IF (ltimer) CALL timer_stop(timer_gmres)
      RETURN
    ENDIF
    
    ! 3) define the tolerance: can be absolute or relative
    IF (abstol) THEN
      tol = tolerance
    ELSE
      tol = tolerance*rn2(1)
    ENDIF
    tol2 = tol**2
    
    ! 4) Arnoldi loop
    arnoldi: DO i = 1, m-1
      ! 4.1) compute the next (i.e. i+1) Krylov vector
      w(:, :) = lhs( v(:,:,i),patch_3D,h_e, p_op_coeff )
      
      ! 4.2) Gram-Schmidt orthogonalization
      
! !ICON_OMP_PARALLEL PRIVATE(rh, h_aux, myThreadNo, k)
      !$   myThreadNo = OMP_GET_THREAD_NUM()
      
      IF (PRESENT(preconditioner)) CALL preconditioner(w(:,:),patch_3D,p_op_coeff,h_e)
      
      gs_orth: DO k = 1, i
        
! !ICON_OMP_DO PRIVATE(jb, jk, nlen) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          IF (jb /= no_of_blocks) THEN
            z(:,jb) = w(:,jb)*v(:,jb,k)
          ELSE
            z(1:end_nproma,jb) = w(1:end_nproma,jb) * v(1:end_nproma,jb,k)
          ENDIF
          WHERE(.NOT.patch_2d%cells%decomp_info%owner_mask(:,jb)) z(:,jb) = 0.0_wp
        ENDDO
! !ICON_OMP_END_DO
        h_aux = omp_global_sum_array(z)
!        CALL dbg_print('anroldi: w', w, method_name, 3, in_subset=patch_2d%cells%owned)
!        CALL dbg_print('anroldi: v', v, method_name, 3, in_subset=patch_2d%cells%owned)
!        write(0,*) "h_aux:", h_aux

        IF (mythreadno == 0) h(k,i) = h_aux
        
! !ICON_OMP_DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          IF (jb /= no_of_blocks) THEN
            nlen = nproma
          ELSE
            nlen = end_nproma
            w(nlen+1:nproma,jb) = 0._wp
          ENDIF
          w(1:nlen,jb) = w(1:nlen,jb) - h_aux*v(1:nlen,jb,k)
        ENDDO
! !ICON_OMP_END_DO
        
      ENDDO gs_orth
      
      ! 4.3) new element for h
      
! !ICON_OMP_DO PRIVATE(jb, jk, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, no_of_blocks
        IF (jb /= no_of_blocks) THEN
          z(:,jb) = w(:,jb) * w(:,jb)
        ELSE
          z(1:end_nproma,jb) = w(1:end_nproma,jb)*w(1:end_nproma,jb)
        ENDIF
        WHERE(.NOT.patch_2d%cells%decomp_info%owner_mask(:,jb)) z(:,jb) = 0.0_wp
      ENDDO
! !ICON_OMP_END_DO
      
      h_aux = SQRT(omp_global_sum_array(z))
      
      IF (mythreadno == 0) h(i+1,i) = h_aux
      
      IF (h_aux < tol2) THEN
        done = .TRUE.
      ELSE
        done = .FALSE.
      ENDIF
      
      ! 4.4) if w is independent from v, add v(:,:,:,i+1)
      IF (.NOT. done) THEN
        rh = 1.0_wp/h_aux
! !ICON_OMP_DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          v(:,jb,i+1) = w(:,jb)*rh
        ENDDO
! !ICON_OMP_END_DO NOWAIT
      ENDIF
      
! !ICON_OMP_END_PARALLEL
      
      ! 4.5) apply the rotation matrices
      rotation: DO k = 1, i-1
        hki  = h(k  ,i)
        hk1i = h(k+1,i)
        h(k  ,i) =  c(k)*hki + s(k)*hk1i
        h(k+1,i) = -s(k)*hki + c(k)*hk1i
      ENDDO rotation
      
      ! 4.6) compute the new (i.e. i) rotation
      den = SQRT(h(i,i)**2 + h(i+1,i)**2)
      c(i) = h(i  ,i) / den
      s(i) = h(i+1,i) / den
      
      ! 4.7) complete applying the rotation matrices
      h(i,i) =  c(i)*h(i,i) + s(i)*h(i+1,i)
      
      ! 4.8) compute new residual norm
      rn2(i+1) = -s(i)*rn2(i)
      rn2(i  ) =  c(i)*rn2(i)
      
      ! 4.9) check whether we are done
      IF ( done .OR. (ABS(rn2(i+1)) < tol) ) EXIT arnoldi
    ENDDO arnoldi
    
    ! 5) check whether we are here because we have exceeded m
    !    Notice that if the arnoldi loop has been terminated with the
    !    exit instruction, we have i=m-1, while if the loop has been
    !    terminated by the loop counter we have i=m (i is incremented
    !    after finishing the do loop)
    IF (i == m) THEN
      maxiterex = .TRUE.
      ! only the first m-1 columns h(i,i) have been set in the arnoldi
      ! loop, but we have i=m. We thus decrement i.
      i = i-1
    ENDIF
    
    nkry = i
    niter = nkry+1
    
    ! 6) compute the coefficient of the Krylov expansion (back sub.)
    krylov: DO i = nkry, 1, -1
      ! coefficients are stored in c for convenience
      ci = rn2(i)
      DO k = nkry, i+1, -1
        ci = ci - h(i,k)*c(k)
      ENDDO
      c(i) = ci/h(i,i)
    ENDDO krylov
    
    ! 7) evaluate the Krylov expansion
! !ICON_OMP_PARALLEL
! !ICON_OMP_DO PRIVATE(jb, i) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      DO i = 1, nkry
        x(:,jb) = x(:,jb) + c(i)*v(:,jb,i)
      ENDDO
    ENDDO
! !ICON_OMP_END_DO NOWAIT
! !ICON_OMP_END_PARALLEL
    IF (ltimer) CALL timer_stop(timer_gmres)
    
    res(1:niter) = ABS(rn2(1:niter))
    
  END SUBROUTINE gmres_oce_old
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  Original version from http://www.netlib.org/templates/templates.pdf
  !!  F90 rewriting by Marco Restelli.
  !!  Adapted for use in ICOHDC by Hui Wan, MPI-M (2007-12-15)
  !!  Included blocking by Marco Restelli (2008-09-22)
  !!  MPI Parallelization by Rainer Johanni (2009-11-11)
  !! interpolation state removed from parameter list and adapted to 2D arrays, Peter Korn (2010-04)
  !! @par
  !! inital guess, overwritten with the solution
  !!
  SUBROUTINE gmres_oce_e2e( x,lhs, &
    & patch_3D,                  &
    & level,                         &
    & p_op_coeff, b,               &
    & tolerance,abstol,m,maxiterex,niter,res,    &
    & preconditioner)
    !
    ! !DESCRIPTION
    !  GMRES solver for linear systems. Notice that, in principle, r and b
    !  can be of any TYPE for which the following operations are defined:
    !   +, -, dot_product
    !
    REAL(wp), INTENT(inout) :: x(:,:)
    INTEGER,     INTENT(in) :: level
    ! REAL(wp), INTENT(IN) :: thickness_c(:,:)
    ! REAL(wp), INTENT(IN) :: old_h(:,:)
    ! patch info needed for calculating lhs
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    TYPE(t_operator_coeff), INTENT(in):: p_op_coeff
    ! right-hand side: same shape as x
    REAL(wp), INTENT(inout) :: b(:,:) ! same size as x
    REAL(wp), INTENT(in) :: tolerance ! (relative or absolute) tolerance
    LOGICAL,  INTENT(in) :: abstol    ! .true. for absolute tolerance,
    ! .false. for relative tolerance
    INTEGER,  INTENT(in) :: m         ! maximum number of iterations
    LOGICAL,  INTENT(out) :: maxiterex ! true if reached m iterations
    INTEGER,  INTENT(out) :: niter    ! number of iterations (defined
    ! as the number of evaluation of
    ! the lhs which, due to the
    ! initial evaluation of the
    ! residual, is equal to the
    ! number of Arnoldi iterations +1)
    ! norms of the residual (convergence history); an argument of
    ! dimension at least m is required
    REAL(wp), INTENT(inout) :: res(:) ! (m)
    
      ! !   FUNCTION lhs(x,old_h, patch_2D,p_patch_3D,coeff, h_e, thickness_c, p_op_coeff) RESULT(ax)
      ! !     USE mo_kind, ONLY: wp
      ! !     USE mo_model_domain, ONLY: t_patch, t_patch_3D
      ! !     USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
      ! !     REAL(wp),    INTENT(inout) :: x(:,:)  ! inout for sync
      ! !     REAL(wp), INTENT(IN) :: old_h(:,:)
      ! !     TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
      ! !     TYPE(t_patch_3D ),TARGET, INTENT(INOUT)   :: p_patch_3D
      ! !     REAL(wp),    INTENT(in) :: coeff
      ! !     REAL(wp),    INTENT(in) :: h_e(:,:)
      ! !     REAL(wp),    INTENT(in) :: thickness_c(:,:)
      ! !     TYPE(t_operator_coeff),INTENT(IN)  :: p_op_coeff
      ! !     REAL(wp) :: ax( SIZE(x,1) , SIZE(x,2) ) ! same as x
      ! !   ENDFUNCTION lhs
    INTERFACE   ! left-hand-side: A*x
      FUNCTION lhs(x,patch_3D, p_op_coeff, level) result(ax)
        USE mo_kind, ONLY: wp
        USE mo_model_domain, ONLY: t_patch, t_patch_3d
        USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
        REAL(wp),    INTENT(inout)             :: x(:,:)  ! inout for sync
        TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
        TYPE(t_operator_coeff),INTENT(in)      :: p_op_coeff
        INTEGER                                :: level
        REAL(wp) :: ax( SIZE(x,1) , SIZE(x,2) ) ! same as x
      END FUNCTION lhs
    END INTERFACE
    
    INTERFACE   ! preconditioner
      SUBROUTINE preconditioner(r,patch_2d, patch_3D, p_op_coeff,h_e)
        USE mo_kind, ONLY: wp
        USE mo_model_domain, ONLY: t_patch, t_patch_3d
        USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
        REAL(wp), INTENT(inout)           :: r(:,:)
        TYPE(t_patch), TARGET, INTENT(in) :: patch_2d
        TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
        TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
        REAL(wp),    INTENT(in)           :: h_e(:,:)
      END SUBROUTINE preconditioner
    END INTERFACE
    
    OPTIONAL :: preconditioner
    
    ! !LOCAL VARIABLES
    
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: done
    INTEGER ::     &
      & i,           & ! index for the Arnoldi loop
      & k,           & ! index for the Gram Schmidt orthogonalization
      & nkry           ! dimension of the Krylov space
    REAL(wp) ::    &
      & tol,         & ! effective tolerance (absolute o relative)
      & tol2,        &
      & r(SIZE(x,1),SIZE(x,2)),  & ! residual
      & rn2(m),      &             ! two-norm of the residual
      & v(SIZE(x,1),SIZE(x,2),m),& ! Krylov basis
      & w(SIZE(x,1),SIZE(x,2)),  & ! new Krylov vector
      & h(m,m),      &             ! Hessemberg matrix
      & hki, hk1i,   &
      & c(m), s(m),  &             ! rotation matrices
      & den, ci
    
    REAL(wp) :: rrn2, rn2_aux, h_aux, rh
    INTEGER :: jb, jk, nlen
    
    INTEGER :: no_of_blocks, end_nproma
    
    REAL(wp) :: z(SIZE(x,1),SIZE(x,2)) ! needed for global sums in p_test_run
    
    INTEGER :: mythreadno
    !$ INTEGER OMP_GET_THREAD_NUM
    !-------------------------------------------------------------------------
    patch_2d =>  patch_3D%p_patch_2d(1)
    
    
    ! 0) set module variables and initialize maxiterex
    no_of_blocks = patch_2d%edges%in_domain%end_block
    end_nproma   = patch_2d%edges%in_domain%end_index
    
    maxiterex = .FALSE.
    z(:,:) = 0._wp
    v(:,:,:)  = 0.0_wp
    r(:,:)    = 0.0_wp
    
    ! 1) compute the preconditioned residual
    
    !w(:,:) = lhs(x(:,:),old_h, patch_2D, p_patch_3D,coeff, h_e, thickness_c, p_op_coeff)
    w(:, :) =lhs(x,patch_3D, p_op_coeff, level)
    
    IF (ltimer) CALL timer_start(timer_gmres)
    
    mythreadno = 0
! !ICON_OMP_PARALLEL PRIVATE(rn2_aux, rrn2, myThreadNo)
    !$ myThreadNo = OMP_GET_THREAD_NUM()
    
    
! !ICON_OMP_DO PRIVATE(jb,nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      IF (jb /= no_of_blocks) THEN
        nlen = nproma
      ELSE
        nlen = end_nproma
        r(nlen+1:nproma,jb) = 0.0_wp
      ENDIF
      r(1:nlen,jb) = b(1:nlen,jb)-w(1:nlen,jb)
    ENDDO
! !ICON_OMP_END_DO
    
    
! !ICON_OMP_DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      IF (jb /= no_of_blocks) THEN
        z(:,jb) = r(:,jb)*r(:,jb)
      ELSE
        z(1:end_nproma,jb) = r(1:end_nproma,jb)*r(1:end_nproma,jb)
      ENDIF
      WHERE(.NOT.patch_2d%edges%decomp_info%owner_mask(:,jb)) z(:,jb) = 0.0_wp
    ENDDO
! !ICON_OMP_END_DO
    
    rn2_aux = SQRT(omp_global_sum_array(z(:,1:no_of_blocks)))
    
    
    IF (mythreadno == 0) rn2(1) = rn2_aux
    
    ! 2) compute the first vector of the Krylov space
    IF (rn2_aux /= 0.0_wp) THEN
      rrn2 = 1.0_wp/rn2_aux
! !ICON_OMP_DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, no_of_blocks
        v(:,jb,1) = r(:,jb)*rrn2
      ENDDO
! !ICON_OMP_END_DO NOWAIT
    ENDIF
    
! !ICON_OMP_END_PARALLEL
    
    IF (rn2(1) == 0.0_wp) THEN ! already done
      !    print*,' gmres: rn2(1)=0.0 '
      ! #slo# 2010-06-15 trying to exit with niter=1 and residual=0.0
      niter = 0
      !     niter = 1
      res(1) = ABS(rn2(1))
      IF (ltimer) CALL timer_stop(timer_gmres)
      RETURN
    ENDIF
    
    ! 3) define the tolerance: can be absolute or relative
    IF (abstol) THEN
      tol = tolerance
    ELSE
      tol = tolerance*rn2(1)
    ENDIF
    tol2 = tol**2
    
    ! 4) Arnoldi loop
    arnoldi: DO i = 1, m-1
      ! 4.1) compute the next (i.e. i+1) Krylov vector
      !w(:,:) = lhs( v(:,:,i),old_h, patch_2D,p_patch_3D,coeff,h_e, thickness_c, p_op_coeff )
       w(:, :) =lhs(v(:,:,i),patch_3D, p_op_coeff, level)
    
      ! 4.2) Gram-Schmidt orthogonalization
      
! !ICON_OMP_PARALLEL PRIVATE(rh, h_aux, myThreadNo, k)
      !$   myThreadNo = OMP_GET_THREAD_NUM()
      
      
      gs_orth: DO k = 1, i
        
! !ICON_OMP_DO PRIVATE(jb, jk, nlen) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          IF (jb /= no_of_blocks) THEN
            z(:,jb) = w(:,jb)*v(:,jb,k)
          ELSE
            z(1:end_nproma,jb) = w(1:end_nproma,jb) * v(1:end_nproma,jb,k)
          ENDIF
          WHERE(.NOT.patch_2d%edges%decomp_info%owner_mask(:,jb)) z(:,jb) = 0.0_wp
        ENDDO
! !ICON_OMP_END_DO
        h_aux = omp_global_sum_array(z)
        
        IF (mythreadno == 0) h(k,i) = h_aux
        
! !ICON_OMP_DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          IF (jb /= no_of_blocks) THEN
            nlen = nproma
          ELSE
            nlen = end_nproma
            w(nlen+1:nproma,jb) = 0._wp
          ENDIF
          w(1:nlen,jb) = w(1:nlen,jb) - h_aux*v(1:nlen,jb,k)
        ENDDO
! !ICON_OMP_END_DO
        
      ENDDO gs_orth
      
      ! 4.3) new element for h
! !ICON_OMP_DO PRIVATE(jb, jk, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, no_of_blocks
        IF (jb /= no_of_blocks) THEN
          z(:,jb) = w(:,jb) * w(:,jb)
        ELSE
          z(1:end_nproma,jb) = w(1:end_nproma,jb)*w(1:end_nproma,jb)
        ENDIF
        WHERE(.NOT.patch_2d%edges%decomp_info%owner_mask(:,jb)) z(:,jb) = 0.0_wp
      ENDDO
! !ICON_OMP_END_DO
      
      h_aux = SQRT(omp_global_sum_array(z))
      
      IF (mythreadno == 0) h(i+1,i) = h_aux
      
      IF (h_aux < tol2) THEN
        done = .TRUE.
      ELSE
        done = .FALSE.
      ENDIF
      
      ! 4.4) if w is independent from v, add v(:,:,:,i+1)
      IF (.NOT. done) THEN
        rh = 1.0_wp/h_aux
! !ICON_OMP_DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          v(:,jb,i+1) = w(:,jb)*rh
        ENDDO
! !ICON_OMP_END_DO NOWAIT
      ENDIF
      
! !ICON_OMP_END_PARALLEL
      
      ! 4.5) apply the rotation matrices
      rotation: DO k = 1, i-1
        hki  = h(k  ,i)
        hk1i = h(k+1,i)
        h(k  ,i) =  c(k)*hki + s(k)*hk1i
        h(k+1,i) = -s(k)*hki + c(k)*hk1i
      ENDDO rotation
      
      ! 4.6) compute the new (i.e. i) rotation
      den = SQRT(h(i,i)**2 + h(i+1,i)**2)
      c(i) = h(i  ,i) / den
      s(i) = h(i+1,i) / den
      
      ! 4.7) complete applying the rotation matrices
      h(i,i) =  c(i)*h(i,i) + s(i)*h(i+1,i)
      
      ! 4.8) compute new residual norm
      rn2(i+1) = -s(i)*rn2(i)
      rn2(i  ) =  c(i)*rn2(i)
      
      ! 4.9) check whether we are done
      IF ( done .OR. (ABS(rn2(i+1)) < tol) ) EXIT arnoldi
    ENDDO arnoldi
    
    ! 5) check whether we are here because we have exceeded m
    !    Notice that if the arnoldi loop has been terminated with the
    !    exit instruction, we have i=m-1, while if the loop has been
    !    terminated by the loop counter we have i=m (i is incremented
    !    after finishing the do loop)
    IF (i == m) THEN
      maxiterex = .TRUE.
      ! only the first m-1 columns h(i,i) have been set in the arnoldi
      ! loop, but we have i=m. We thus decrement i.
      i = i-1
    ENDIF
    
    nkry = i
    niter = nkry+1
    
    ! 6) compute the coefficient of the Krylov expansion (back sub.)
    krylov: DO i = nkry, 1, -1
      ! coefficients are stored in c for convenience
      ci = rn2(i)
      DO k = nkry, i+1, -1
        ci = ci - h(i,k)*c(k)
      ENDDO
      c(i) = ci/h(i,i)
    ENDDO krylov
    
    ! 7) evaluate the Krylov expansion
! !ICON_OMP_PARALLEL
! !ICON_OMP_DO PRIVATE(jb, i) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      DO i = 1, nkry
        x(:,jb) = x(:,jb) + c(i)*v(:,jb,i)
      ENDDO
    ENDDO
! !ICON_OMP_END_DO NOWAIT
! !ICON_OMP_END_PARALLEL

    IF (ltimer) CALL timer_stop(timer_gmres)
    
    res(1:niter) = ABS(rn2(1:niter))
    
  END SUBROUTINE gmres_oce_e2e  
  !-------------------------------------------------------------------------
  
END MODULE mo_ocean_gmres
