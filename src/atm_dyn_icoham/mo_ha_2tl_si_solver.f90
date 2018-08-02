!>
!! General comments: provide an implementation of the GMRES solver for
!! linear systems with the modified Gram-Schmidt orthogonalization.
!! Implementation as in "Templates for the Solution of Linear Systems:
!! Building Blocks for Iterative Methods".
!! This module is used by the 2-time-level semi-implicit time stepping
!! scheme of the ICON hydrostatic atmospheric dynamical core.
!!
!! @par Revision History
!! Original version from http://www.netlib.org/templates/templates.pdf
!! F90 rewriting by Marco Restelli.
!! Adapted for use in ICOHDC by Hui Wan, MPI-M (2007-12-15)
!! Cosmetic changes following the ICON programming guide by Hui Wan,
!! MPI-M (2007-12-18)
!! Included blocking by Marco Restelli (2008-09-23)
!! Dummy argument for preconditioner shifted to last position
!! and made optional, Marco Giorgetta (2009-02-14)
!! Inlining of former functions active_dot and norm2 to simplify OpenMP
!! parallelization, reduction of parallel sections, Guenther Zaengl (2009-06-16)
!! Modification for the 2-time-level semi-implicit time stepping scheme
!! by Hui Wan (MPI-M, 2009-12).
!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ha_2tl_si_solver

  USE mo_kind,         ONLY: wp
  USE mo_model_domain, ONLY: t_patch
  USE mo_timer,        ONLY: timer_start, timer_stop, timer_gmres
  USE mo_intp_data_strc, ONLY: t_int_state
  USE mo_sync,         ONLY: omp_global_sum_array
  USE mo_parallel_config,  ONLY: nproma
  USE mo_run_config,   ONLY: ltimer

  IMPLICIT NONE

  PRIVATE


! !PUBLIC MEMBER FUNCTIONS/SUBROUTINES
  PUBLIC :: solver    ! solver
! !DEFINED PARAMETERS:
  CHARACTER(len=*), PARAMETER :: this_mod_name = 'm_ha_2tlsi_solver'

! module variables used to specify the array bounds in the private
! module subroutine and functions (notice that making active_dot an
! internal procedure of gmres would not allow for overloading)


CONTAINS

 !>
 !!
 !! @par Revision History
 !!  Original version from http://www.netlib.org/templates/templates.pdf
 !!  F90 rewriting by Marco Restelli.
 !!  Adapted for use in ICOHDC by Hui Wan, MPI-M (2007-12-15)
 !!  Included blocking by Marco Restelli (2008-09-22)
 !!  MPI Parallelization by Rainer Johanni (2009-11-11)
 !! @par
 !!
 SUBROUTINE solver(x, &
                   lhs,temp,rdelp_c,delp_e,rdlnpr,rdalpha,cgradps,   &
                   curr_patch,p_int,nblks,npromz,coeff,b,  &
                   tolerance,abstol,m,maxiterex,niter,res, &
                   preconditioner)

  ! inital guess, overwritten with the solution
  REAL(wp), INTENT(INOUT) :: x(:,:,:)
  ! temperature- and pressure-related quantities needed for the lhs
  REAL(wp), INTENT(IN) :: temp(:,:,:), rdelp_c(:,:,:), delp_e(:,:,:), &
                          rdlnpr(:,:,:), rdalpha(:,:,:), cgradps(:,:,:)
  ! patch info needed for calculating lhs
  TYPE(t_patch), INTENT(INOUT) :: curr_patch
  ! interpolation state (needed for divergence averaging)
  TYPE(t_int_state), INTENT(IN) :: p_int
  ! index defining the "active" region of the arrays
  INTEGER, INTENT(IN) :: nblks, npromz
  ! parameter used in calculating the lhs
  REAL(wp), INTENT(IN) :: coeff
  ! right-hand side: same shape as x
  REAL(wp), INTENT(IN) :: b(:,:,:) ! same size as x
  REAL(wp), INTENT(IN) :: tolerance ! (relative or absolute) tolerance
  LOGICAL,  INTENT(IN) :: abstol    ! .true. for absolute tolerance,
                                    ! .false. for relative tolerance
  INTEGER,  INTENT(IN) :: m         ! maximum number of iterations
  LOGICAL,  INTENT(OUT) :: maxiterex ! true if reached m iterations
  INTEGER,  INTENT(OUT) :: niter    ! number of iterations (defined
                                    ! as the number of evaluation of
                                    ! the lhs which, due to the
                                    ! initial evaluation of the
                                    ! residual, is equal to the
                                    ! number of Arnoldi iterations +1)
  ! norms of the residual (convergence history); an argument of
  ! dimension at least m is required
  REAL(wp), INTENT(INOUT) :: res(:) ! (m)

  INTERFACE   ! left-hand-side: A*x
    FUNCTION lhs(x,temp,rdelp_c,delp_e,rdlnpr,rdalpha,cgradps, &
                 curr_patch,p_int,coeff)  RESULT(ax)
      USE mo_kind, ONLY: wp
      USE mo_model_domain, ONLY: t_patch
      USE mo_intp_data_strc,ONLY: t_int_state
      REAL(wp),    INTENT(in) :: x(:,:,:)
      REAL(wp),    INTENT(in) :: temp(:,:,:), rdelp_c(:,:,:), delp_e(:,:,:)
      REAL(wp),    INTENT(in) :: rdlnpr(:,:,:), rdalpha(:,:,:), cgradps(:,:,:)
      TYPE(t_patch), INTENT(inout) :: curr_patch
      TYPE(t_int_state), INTENT(in) :: p_int
      REAL(wp),    INTENT(in) :: coeff
      REAL(wp) :: ax( SIZE(x,1) , SIZE(x,2) , SIZE(x,3) ) ! same as x
    ENDFUNCTION lhs
  END INTERFACE

  INTERFACE   ! preconditioner
   SUBROUTINE preconditioner(r)
     USE mo_kind, ONLY: wp
     REAL(wp), INTENT(inout) :: r(:,:,:)
   END SUBROUTINE preconditioner
  END INTERFACE

  OPTIONAL :: preconditioner

! !LOCAL VARIABLES

  LOGICAL :: done
  INTEGER :: &
   i,           & ! index for the Arnoldi loop
   k,           & ! index for the Gram Schmidt orthogonalization
   nkry           ! dimension of the Krylov space
  REAL(wp) :: &
   tol,         & ! effective tolerance (absolute o relative)
   tol2,        &
   r(SIZE(x,1),SIZE(x,2),SIZE(x,3)),  & ! residual
   rn2(m),      & ! two-norm of the residual
   v(SIZE(x,1),SIZE(x,2),SIZE(x,3),m),& ! Krylov basis
   w(SIZE(x,1),SIZE(x,2),SIZE(x,3)),  & ! new Krylov vector
   h(m,m),      & ! Hessemberg matrix
   hki, hk1i,   &
   c(m), s(m),  & ! rotation matrices
   den, ci

  REAL(wp):: rrn2, rn2_aux, h_aux, rh
  
! NOMPI is disabled for checking the bit-reproducability with mpi versions
! this is not efficient though, and probably will be re-intronduced
#ifdef NOMPI_DISABLED
  REAL(wp):: sum_aux(nblks)
#endif
  INTEGER :: jb, nlen, ndim2

  INTEGER :: mnblks, mnpromz

#ifndef NOMPI_DISABLED
  REAL(wp) :: z(SIZE(x,1),SIZE(x,3)) ! needed for global sums
  INTEGER :: jk
#endif

  INTEGER :: myThreadNo
!$ INTEGER OMP_GET_THREAD_NUM
!-------------------------------------------------------------------------


   ! 0) set module variables and initialize maxiterex

   mnblks = nblks
   mnpromz = npromz
#ifndef NOMPI_DISABLED
   z(:,:) = 0._wp
#endif

   maxiterex = .FALSE.
   ndim2 = SIZE(x,2)

   ! 1) compute the preconditioned residual

   w(:,:,:) = lhs(x(:,:,:),temp,rdelp_c,delp_e, &
                  rdlnpr,rdalpha,cgradps,      &
                  curr_patch,p_int,coeff)

   IF (ltimer) CALL timer_start(timer_gmres)

   myThreadNo = 0
!$OMP PARALLEL PRIVATE(rn2_aux, rrn2, myThreadNo)
!$ myThreadNo = OMP_GET_THREAD_NUM()

!$OMP DO PRIVATE(jb,nlen) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = 1, mnblks
     IF (jb /= mnblks) THEN
       nlen = nproma
     ELSE
       nlen = mnpromz
       r(nlen+1:nproma,:,jb) = 0._wp
     ENDIF
      r(1:nlen,:,jb) = b(1:nlen,:,jb)-w(1:nlen,:,jb)
   ENDDO
!$OMP END DO

   IF (PRESENT(preconditioner)) CALL preconditioner(r(:,:,:))

#ifdef NOMPI_DISABLED
!$OMP DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
     DO jb = 1, mnblks
       IF (jb /= mnblks) THEN
         sum_aux(jb) = SUM(r(:,:,jb)*r(:,:,jb))
       ELSE
         sum_aux(jb) = SUM( r(1:mnpromz,:,jb)*r(1:mnpromz,:,jb) )
       ENDIF
     ENDDO
!$OMP END DO

   rn2_aux = SQRT(SUM(sum_aux))

#else
!$OMP DO PRIVATE(jb, jk, nlen) ICON_OMP_DEFAULT_SCHEDULE
     DO jb = 1, mnblks
       IF (jb /= mnblks) THEN
         nlen = nproma
       ELSE
         nlen = mnpromz
       ENDIF
       z(1:nlen,jb) = r(1:nlen,1,jb)*r(1:nlen,1,jb)
       DO jk = 2,ndim2
         z(1:nlen,jb) = z(1:nlen,jb) + r(1:nlen,jk,jb)*r(1:nlen,jk,jb)
       ENDDO
       WHERE(.NOT.curr_patch%edges%decomp_info%owner_mask(:,jb)) z(:,jb) = 0._wp
     ENDDO
!$OMP END DO

   rn2_aux = SQRT(omp_global_sum_array(z))
#endif

   IF (myThreadNo == 0) rn2(1) = rn2_aux


   ! 2) compute the first vector of the Krylov space
   IF (rn2_aux /= 0.0_wp) THEN
      rrn2 = 1.0_wp/rn2_aux
!$OMP DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, mnblks
        v(:,:,jb,1) = r(:,:,jb)*rrn2
      ENDDO
!$OMP END DO NOWAIT
   ENDIF

!$OMP END PARALLEL

   IF (rn2(1) == 0.0_wp) THEN ! already done
     niter = 0
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
     w(:,:,:) = lhs( v(:,:,:,i), temp, rdelp_c, delp_e, &
                     rdlnpr, rdalpha, cgradps,          &
                     curr_patch, p_int, coeff )

     ! 4.2) Gram-Schmidt orthogonalization

!$OMP PARALLEL PRIVATE(rh, h_aux, myThreadNo, k)
!$   myThreadNo = OMP_GET_THREAD_NUM()

     IF (PRESENT(preconditioner)) CALL preconditioner(r(:,:,:))

     gs_orth: DO k = 1, i

#ifdef NOMPI_DISABLED
!$OMP DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
     DO jb = 1, mnblks
       IF (jb /= mnblks) THEN
         sum_aux(jb) = SUM(w(:,:,jb)*v(:,:,jb,k))
       ELSE
         sum_aux(jb) = SUM( w(1:mnpromz,:,jb)*v(1:mnpromz,:,jb,k) )
       ENDIF
     ENDDO
!$OMP END DO

     h_aux = SUM(sum_aux)
#else

!$OMP DO PRIVATE(jb, jk, nlen) ICON_OMP_DEFAULT_SCHEDULE
     DO jb = 1, mnblks
       IF (jb /= mnblks) THEN
         nlen = nproma
       ELSE
         nlen = mnpromz
       ENDIF
       z(1:nlen,jb) = w(1:nlen,1,jb)*v(1:nlen,1,jb,k)
       DO jk = 2,ndim2
         z(1:nlen,jb) = z(1:nlen,jb) + w(1:nlen,jk,jb)*v(1:nlen,jk,jb,k)
       ENDDO
       WHERE(.NOT.curr_patch%edges%decomp_info%owner_mask(:,jb)) z(:,jb) = 0._wp
     ENDDO
!$OMP END DO

     h_aux = omp_global_sum_array(z)
#endif

     IF (myThreadNo == 0) h(k,i) = h_aux

!$OMP DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
     DO jb = 1, mnblks
       IF (jb /= mnblks) THEN
         nlen = nproma
       ELSE
         nlen = mnpromz
         w(nlen+1:nproma,:,jb) = 0._wp
       ENDIF
       w(1:nlen,:,jb) = w(1:nlen,:,jb) - h_aux*v(1:nlen,:,jb,k)
     ENDDO
!$OMP END DO

   ENDDO gs_orth

     ! 4.3) new element for h

#ifdef NOMPI_DISABLED
!$OMP DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
     DO jb = 1, mnblks
       IF (jb /= mnblks) THEN
         sum_aux(jb) = SUM(w(:,:,jb)*w(:,:,jb))
       ELSE
         sum_aux(jb) = SUM( w(1:mnpromz,:,jb)*w(1:mnpromz,:,jb) )
       ENDIF
     ENDDO
!$OMP END DO

     h_aux = SQRT(SUM(sum_aux))

#else
!$OMP DO PRIVATE(jb, jk, nlen) ICON_OMP_DEFAULT_SCHEDULE
     DO jb = 1, mnblks
       IF (jb /= mnblks) THEN
         nlen = nproma
       ELSE
         nlen = mnpromz
       ENDIF
       z(1:nlen,jb) = w(1:nlen,1,jb)*w(1:nlen,1,jb)
       DO jk = 2,ndim2
         z(1:nlen,jb) = z(1:nlen,jb) + w(1:nlen,jk,jb)*w(1:nlen,jk,jb)
       ENDDO
       WHERE(.NOT.curr_patch%edges%decomp_info%owner_mask(:,jb)) z(:,jb) = 0._wp
     ENDDO
!$OMP END DO

     h_aux = SQRT(omp_global_sum_array(z))
#endif

     IF (myThreadNo == 0) h(i+1,i) = h_aux

     IF (h_aux < tol2) THEN
       done = .TRUE.
     ELSE
       done = .FALSE.
     ENDIF

     ! 4.4) if w is independent from v, add v(:,:,:,i+1)
     IF (.NOT. done) THEN
       rh = 1.0_wp/h_aux
!$OMP DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
       DO jb = 1, mnblks
         v(:,:,jb,i+1) = w(:,:,jb)*rh
      ENDDO
!$OMP END DO NOWAIT
     ENDIF

!$OMP END PARALLEL

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
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = 1, mnblks
     DO i = 1, nkry
       x(:,:,jb) = x(:,:,jb) + c(i)*v(:,:,jb,i)
    ENDDO
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
   IF (ltimer) CALL timer_stop(timer_gmres)

   res(1:niter) = ABS(rn2(1:niter))

 END SUBROUTINE solver


!-------------------------------------------------------------------------
END MODULE mo_ha_2tl_si_solver
!EOP
