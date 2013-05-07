
! LL: xlc has trouble optimizing routines with implicit shaped parameters
#ifdef __xlC__
@process nohot
! @PROCESS NOOPTIMIZE
#endif
!>
!!   Contains the implementation of various mathematical algorithms.
!!
!!   Contains the implementation of various mathematical algorithms
!!   used by the shallow water model.
!!
!! @par Revision History
!!  Developed  by Luca Bonaventura (2002-4).
!!  Modified to ProTeX-style by  Luca Bonaventura and Thomas Heinze (2004).
!!  Adapted to new data structure by Thomas Heinze,
!!  Peter Korn and Luca Bonaventura (2005).
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!!  Modification by Peter Korn and Luca Bonaventura (2006-07-21):
!!  - moved here linear algebra subroutines and other auxiliary functions
!!  - added some Protex documentation and great cleanup
!!  Modification by Thomas Heinze (2006-09-07):
!!  - added functions barycenter and bary_center
!!  Modification by Thomas Heinze (2006-10-19):
!!  - added functions vector_product and integral_over_triangle
!!  Modification by Tobias Ruppert and Thomas Heinze (2006-11-14):
!!  - added functions solve_chol and choldec
!!  - renamed function solve to solve_lu
!!  Modification by Thomas Heinze (2006-11-16):
!!  - added functions cvec2gvec and gvec2cvec
!!  Modification by Peter Korn, MPI-M, (2006-11-23):
!!  - replaced vertex_index by vertex_idx
!!  Modification by Hui Wan (2007-02-22):
!!  - functions barycenter and bary_center removed.
!!    (These two functions used TYPE grid, thus should not be put
!!     in to shr_general.)
!!  - function ll2xyz removed because it had been replaced by
!!    gc2cc and was not used anymore.
!!  - type cartesian_coordinates and type geographical_coordinates
!!    moved to this module from mo_model_domain
!!  - functions func_f and dxg moved from mo_functions to here
!!  Problem related to the range of longitude in the spherical coordnate
!!  identified by Th. Heinze and P. Ripodas (2007-02-26). Correction in 'cc2gc'
!!  made by H. Wan (2007-02-27).
!!  Modification by Thomas Heinze, DWD, (2007-07-26):
!!  - including all the improvements of Tobias Ruppert's diploma thesis
!!  - several changes according to the programming guide
!!  - functions func_f and dxg moved from here to mo_interpolation
!!  Modification by Daniel Reinert, DWD, (2009-07-20)
!!  - added subroutine qrdec
!!  Modification by Daniel Reinert, DWD, (2012-04-04)
!!  - added function which can be used to check whether to line segments
!!    intersect (in 2D cartesian system)
!!  Modification by Daniel Reinert, DWD, (2012-04-05)
!!  - added function which computes the intersection point between 2 lines
!!    (2D cartesian)
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!
MODULE mo_math_utility_solvers
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_parallel_config,     ONLY: nproma
  
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
      
  PUBLIC :: solve_lu
  PUBLIC :: solve_chol
  PUBLIC :: solve_chol_v
  PUBLIC :: ludec
  PUBLIC :: qrdec
  PUBLIC :: choldec
  PUBLIC :: choldec_v
  PUBLIC :: inv_mat
    
  
CONTAINS
  
  
  !-------------------------------------------------------------------------
  !>
  !!  Home made routine for LU decomposition.
  !!
  !!  Home made routine for LU decomposition
  !!  to be substituted by call from linear algebra library;
  !!  provides input matrix for routine <i>solve</i>;<br>
  !!  @f$k\_dim@f$  : matrix dimension <br>
  !!  @f$p\_a@f$    : input matrix, overwritten by lu decomposition <br>
  !!  @f$k\_pivot@f$: integer index array recording row pivot line exchanges
  !!
  !! @par Revision History
  !! Developed  by P.Korn (2006).
  !! Revised by Th.Heinze, DWD, (2006-11-01):
  !! - changed stop criterion from 1.e-08 to 1.e-12
  !!
  SUBROUTINE ludec (k_dim, p_a, k_pivot)
    !
    INTEGER, INTENT(in)     :: k_dim              ! matrix dimension
    
    REAL(wp), INTENT(inout) :: p_a(k_dim,k_dim)   ! input matrix
    
    INTEGER, INTENT(out)    :: k_pivot(k_dim)     ! pivoting index
    
    REAL(wp)                :: z_rmax(k_dim)      ! reciprocal of max value for
    ! each row
    REAL(wp)                :: z_temp             ! temporary variable
    REAL(wp)                :: z_max              ! max value
    REAL(wp)                :: z_sum, z_dum       ! variables generated by matrix
    ! entries
    INTEGER :: ji, jj, jk         ! integer over dimension
    INTEGER :: imax
    
    !-----------------------------------------------------------------------
    
    imax = 0
    
    DO ji = 1, k_dim
      z_max = 0._wp
      DO jj = 1, k_dim
        z_temp = ABS(p_a(ji,jj))
        IF(z_temp >= z_max) THEN
          z_max = z_temp
        ENDIF
      ENDDO
#if defined(SINGULARITY_CHECKS)
      IF (z_max <= 1.e-12_wp) THEN
        CALL message ('mo_math_utilities:ludec', &
          & 'error in matrix inversion, nearly singular matrix')
        CALL finish  ('mo_math_utilities:ludec', &
          & 'error in matrix inversion, nearly singular matrix')
      ENDIF
#endif
      
      z_rmax(ji) = 1._wp / z_max
    ENDDO
    
    DO jj = 1, k_dim
      DO ji = 1, jj-1
        z_sum = p_a(ji,jj)
        DO jk = 1, ji-1
          z_sum = z_sum - p_a(ji,jk) * p_a(jk,jj)
        ENDDO
        p_a(ji,jj) = z_sum
      ENDDO
      
      z_max = 0._wp
      
      DO ji = jj, k_dim
        z_sum = p_a(ji,jj)
        DO jk = 1, jj-1
          z_sum = z_sum - p_a(ji,jk) * p_a(jk,jj)
        ENDDO
        p_a(ji,jj) = z_sum
        z_dum = z_rmax(ji) * ABS(z_sum)
        IF (z_dum >= z_max) THEN
          z_max = z_dum
          imax  = ji
        ENDIF
      ENDDO
      
      IF (jj /= imax) THEN
        
        DO jk=1,k_dim
          z_dum        = p_a(imax,jk)
          p_a(imax,jk) = p_a(jj,jk)
          p_a(jj,jk)   = z_dum
        ENDDO
        z_rmax(imax) = z_rmax(jj)
        
      ENDIF
      
      k_pivot(jj) = imax
#if defined(SINGULARITY_CHECKS)
      IF (ABS(p_a(jj,jj)) <= 1.e-12_wp) THEN
        CALL message ('mo_math_utilities:ludec', &
          & 'error  in matrix inversion, nearly singular matrix 2')
        CALL finish  ('mo_math_utilities:ludec', &
          & 'error  in matrix inversion, nearly singular matrix 2')
      ENDIF
#endif
      IF (jj /= k_dim) THEN
        z_dum = 1._wp / p_a(jj,jj)
        DO ji = jj + 1, k_dim
          p_a(ji,jj) = p_a(ji,jj) * z_dum
        ENDDO
      ENDIF
    ENDDO
    
  END SUBROUTINE ludec
  !-------------------------------------------------------------------------  
  
  !-------------------------------------------------------------------------
  !>
  !!  QR decomposition of (m x n) matrix  Applies modified Gram-Schmidt algorithm.
  !!
  !!  QR decomposition of (m x n) matrix  Applies modified Gram-Schmidt algorithm
  !!  @f$p\_inmat@f$ : input matrix (m x n)<br>
  !!  @f$p\_qmat@f$  : Q matrix (orthogonal) (m x n) <br>
  !!  @f$p\_rmat@f$  : R matrix (upper triangular) (n x n)
  !!
  !! @par Revision History
  !! Developed and tested by Daniel Reinert, DWD (2009-11-26)
  !!
  SUBROUTINE qrdec (nrow, ncol, i_startidx, i_endidx, p_inmat, p_qmat, p_rmat)
    !
    ! !LITERATURE
    ! Ruenger et al. (2005): Comparison of different parallel modified
    !                        Gram-Schmidt Algorithms. Lecture Notes in
    !                        Computer Science, 3648, 826-836
    !
    INTEGER, INTENT(in)  :: nrow          ! number of rows (m) of matrix
    ! p\_inmat
    
    INTEGER, INTENT(in)  :: ncol          ! number of columns (n) of matrix
    ! p\_inmat
    
    INTEGER, INTENT(in)  :: i_startidx,  & ! start and and indices for
      & i_endidx       ! loop over cells
    
    REAL(wp), INTENT(in) :: p_inmat(:,:,:)   ! matrix for which QR-
    ! factorization is calculated
    
    REAL(wp), INTENT(out)  :: p_qmat(nproma,nrow,ncol)  ! Q matrix
    REAL(wp), INTENT(out)  :: p_rmat(nproma,ncol,ncol)  ! R matrix
    
    INTEGER :: i, k, jc                   ! loop indices
    REAL(wp) :: z_q(nproma,nrow)          ! column of input matrix
    REAL(wp) :: d_product(nproma)
    
    !-----------------------------------------------------------------------    
    ! Modified Gram-Schmidt algorithm
    !
    p_qmat(i_startidx:i_endidx,:,:) = p_inmat(i_startidx:i_endidx,:,:)
    
    DO i=1, ncol
      
      DO jc = i_startidx, i_endidx
        z_q(jc,1:nrow)      = p_qmat(jc,1:nrow,i)                          ! ith column
        p_rmat(jc,i,i) = SQRT(DOT_PRODUCT(z_q(jc,1:nrow),z_q(jc,1:nrow)))  ! vector length
        z_q(jc,1:nrow)      = z_q(jc,1:nrow)/p_rmat(jc,i,i)                ! normalization
        p_qmat(jc,1:nrow,i) = z_q(jc,1:nrow)
        
        ! loop over all old columns of the input array
        DO k=i+1,ncol
          d_product(jc) = DOT_PRODUCT(p_qmat(jc,1:nrow,k),z_q(jc,1:nrow))
          p_qmat(jc,1:nrow,k) = p_qmat(jc,1:nrow,k) - d_product(jc) * z_q(jc,1:nrow)
          p_rmat(jc,i,k) = d_product(jc)
        ENDDO
      ENDDO
      
    ENDDO
    
  END SUBROUTINE qrdec
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  !! Given a positive-definite symmetric matrix @f$p\_a(1:k\_dim,1:k\_dim)@f$, with.
  !!
  !! Given a positive-definite symmetric matrix @f$p\_a(1:k\_dim,1:k\_dim)@f$, with
  !! physical dimension @f$k\_dim@f$, this routine constructs its Cholesky
  !! decomposition, @f$p\_a = LL^T@f$. On input, only the upper triangle of @f$p\_a@f$
  !! need be given; it is not modified. The Cholesky factor @f$L@f$ is returned
  !! in the lower triangle of @f$p\_a@f$, except for its diagonal elements which are
  !! returned in @f$p\_diag(1:k\_dim)@f$.
  !!
  !! @par Revision History
  !!  Original version by Tobias Ruppert and Thomas Heinze, DWD (2006-11-14)
  !!
  SUBROUTINE choldec (k_dim, p_a, p_diag)
    !
    
    INTEGER, INTENT(in)       :: k_dim                ! matrix dimension
    
    REAL(wp), INTENT(inout)   :: p_a(:,:)             ! input matrix
    
    REAL(wp), INTENT(out)     :: p_diag(:)            ! diagonal of output matrix
    
    REAL(wp)                  :: z_sum                ! variable generated by
    ! matrix entries
    INTEGER :: ji, jj, jk           ! integer over dimension
    
    !-----------------------------------------------------------------------   
    DO ji = 1, k_dim
      DO jj = ji, k_dim
        z_sum = p_a(ji,jj)
        
        DO jk = ji-1, 1, -1
          z_sum = z_sum - p_a(ji,jk) * p_a(jj,jk)
        ENDDO
        
        IF (ji == jj) THEN
#if defined(SINGULARITY_CHECKS)
          IF (z_sum <= 1.e-12_wp) THEN
            CALL message ('mo_math_utilities:choldec',                           &
              & 'error in matrix inversion, nearly singular matrix')
            CALL finish  ('mo_math_utilities:choldec',                           &
              & 'error in matrix inversion, nearly singular matrix')
          ELSE
#endif
            p_diag(ji) = SQRT(z_sum)
#if defined(SINGULARITY_CHECKS)
          ENDIF
#endif
        ELSE
          p_a(jj,ji) = z_sum / p_diag(ji)
        ENDIF
      ENDDO
    ENDDO
    
  END SUBROUTINE choldec
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Vectorized version of choldec.
  !!
  !! Vectorized version of choldec
  !! Given a positive-definite symmetric matrix @f$p\_a(1:k\_dim,1:k\_dim)@f$, with
  !! physical dimension @f$k\_dim@f$, this routine constructs its Cholesky
  !! decomposition, @f$p\_a = LL^T@f$. On input, only the upper triangle of @f$p\_a@f$
  !! need be given; it is not modified. The Cholesky factor @f$L@f$ is returned
  !! in the lower triangle of @f$p\_a@f$, except for its diagonal elements which are
  !! returned in @f$p\_diag(1:k\_dim)@f$.
  !!
  !! @par Revision History
  !!  Original version by Tobias Ruppert and Thomas Heinze, DWD (2006-11-14)
  !!  Vectorized version by Guenther Zaengl, DWD (2009-04-20)
  !!
#ifdef __SX__
  SUBROUTINE choldec_v (istart, iend, k_dim, maxdim, p_a, p_diag)
    
    INTEGER, INTENT(in)       :: istart, iend         ! start and end index    
    INTEGER, INTENT(in)       :: k_dim(:)                ! matrix dimension    
    INTEGER, INTENT(in)       :: maxdim                  ! maximum matrix dimension
    REAL(wp), INTENT(inout)   :: p_a(:,:,:)             ! input matrix    
    REAL(wp), INTENT(out)     :: p_diag(:,:)            ! diagonal of output matrix
    
    REAL(wp)                  :: z_sum(nproma)        ! variable generated by
    ! matrix entries
    INTEGER :: jc, ji, jj, jk           ! integer over dimension
    
    !-----------------------------------------------------------------------    
    DO ji = 1, maxdim
      DO jj = ji, maxdim
        DO jc = istart, iend
          IF ((ji > k_dim(jc)) .OR. (jj > k_dim(jc))) CYCLE
          z_sum(jc) = p_a(jc,ji,jj)
        ENDDO
        DO jk = ji-1, 1, -1
          DO jc = istart, iend
            IF ((ji > k_dim(jc)) .OR. (jj > k_dim(jc))) CYCLE
            z_sum(jc) = z_sum(jc) - p_a(jc,ji,jk) * p_a(jc,jj,jk)
          ENDDO
        ENDDO
        
        IF (ji == jj) THEN
          DO jc = istart, iend
            IF (ji > k_dim(jc)) CYCLE
            p_diag(jc,ji) = SQRT(z_sum(jc))
          ENDDO
        ELSE
          DO jc = istart, iend
            IF ((ji > k_dim(jc)) .OR. (jj > k_dim(jc))) CYCLE
            p_a(jc,jj,ji) = z_sum(jc) / p_diag(jc,ji)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    
  END SUBROUTINE choldec_v
#else
  ! non-vecorized version
  SUBROUTINE choldec_v (istart, iend, k_dim, p_a, p_diag)
    !    
    INTEGER, INTENT(in)       :: istart, iend         ! start and end index    
    INTEGER, INTENT(in)       :: k_dim(:)                ! matrix dimension    
    REAL(wp), INTENT(inout)   :: p_a(:,:,:)             ! input matrix   
    REAL(wp), INTENT(out)     :: p_diag(:,:)            ! diagonal of output matrix
    
    REAL(wp)                  :: z_sum                ! variable generated by
    ! matrix entries
    INTEGER :: jc, ji, jj, jk           ! integer over dimension
    
    !-----------------------------------------------------------------------    
    DO jc = istart, iend
      DO ji = 1, k_dim(jc)
        DO jj = ji, k_dim(jc)
          z_sum = p_a(jc,ji,jj)
          
          DO jk = ji-1, 1, -1
            z_sum = z_sum - p_a(jc,ji,jk) * p_a(jc,jj,jk)
          ENDDO
          
          IF (ji == jj) THEN
#if defined(SINGULARITY_CHECKS)
            IF (z_sum <= 1.e-12_wp) THEN
              WRITE (*,*) p_a
              CALL message ('mo_math_utilities:choldec',                           &
                & 'error in matrix inversion, nearly singular matrix')
              CALL finish  ('mo_math_utilities:choldec',                           &
                & 'error in matrix inversion, nearly singular matrix')
            ELSE
#endif
              p_diag(jc,ji) = SQRT(z_sum)
#if defined(SINGULARITY_CHECKS)
            ENDIF
#endif
          ELSE
            p_a(jc,jj,ji) = z_sum / p_diag(jc,ji)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    
  END SUBROUTINE choldec_v
#endif
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Compute the inverse of matrix a (use only for VERY small matrices!).
  !!
  !!
  !! @par Revision History
  !! Original version by Marco Restelli (2007-11-22)
  !!
  SUBROUTINE inv_mat (a)
    
    REAL(wp), INTENT(inout)   :: a(:,:)     ! input matrix
    
    INTEGER :: k_pivot(SIZE(a,1)), j
    REAL(wp) :: b(SIZE(a,1),SIZE(a,1)), e(SIZE(a,1))
    !-----------------------------------------------------------------------
    
    b = a
    CALL ludec(SIZE(a,1),b,k_pivot)
    
    DO j=1,SIZE(a,1)
      e = 0.0_wp
      e(j) = 1.0_wp
      CALL solve_lu(SIZE(a,1),b,e,a(:,j),k_pivot)
    ENDDO
    
  END SUBROUTINE inv_mat
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Solves the set of @f$k\_dim@f$ linear equations @f$p\_a \cdot p\_x = p\_b@f$, where.
  !!
  !! Solves the set of @f$k\_dim@f$ linear equations @f$p\_a \cdot p\_x = p\_b@f$, where
  !! @f$p\_a@f$ is a positive-definite symmetric matrix with physical dimension @f$k\_dim@f$.
  !! @f$p\_a@f$ and @f$p\_diag@f$ are input as the output of the routine {\\tt choldec}.
  !! Only the lower triangle of @f$p\_a@f$ is accessed. @f$p\_b(1:k\_dim)@f$ is input as
  !! the right-hand side vector. The solution vector is returned in @f$p\_x(1:k\_dim)@f$.
  !! @f$p\_a@f$ and @f$p\_diag@f$ are not modified and can be left in place for successive
  !! calls with different right-hand sides @f$p\_b@f$. @f$p\_b@f$ is not modified unless
  !! you identify @f$p\_b@f$ and @f$p\_x@f$ in the calling sequence, which is allowed.
  !!
  !! @par Revision History
  !!  Original version by Tobias Ruppert and Thomas Heinze, DWD (2006-11-14)
  !!
  SUBROUTINE solve_chol (k_dim, p_a, p_diag, p_b, p_x)
    !
    INTEGER,  INTENT(in)   :: k_dim         ! matrix dimension
    REAL(wp), INTENT(in)   :: p_a(:,:)      ! (k_dim,k_dim) matrix
    REAL(wp), INTENT(in)   :: p_diag(:)     ! (k_dim) diagonal of matrix
    
    REAL(wp), INTENT(inout):: p_b(:)        ! (k_dim) input vector    
    REAL(wp), INTENT(inout):: p_x(:)        ! (k_dim) output vector
    
    REAL(wp)               :: z_sum         ! variable generated by
    ! matrix entries
    INTEGER :: ji, jj        ! integer over dimensions
    
    !-----------------------------------------------------------------------    
    ! Solve L*y = b, storing y in x.    
    DO ji = 1, k_dim
      z_sum = p_b(ji)
      DO jj = ji - 1, 1, -1
        z_sum = z_sum - p_a(ji,jj) * p_x(jj)
      ENDDO
      p_x(ji) = z_sum / p_diag(ji)
    ENDDO
    
    ! Solve L^T * x = y
    
    DO ji = k_dim, 1, -1
      z_sum = p_x(ji)
      DO jj = ji+1, k_dim
        z_sum = z_sum - p_a(jj,ji) * p_x(jj)
      ENDDO
      p_x(ji) = z_sum / p_diag(ji)
    ENDDO
    
  END SUBROUTINE solve_chol
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Solves the set of @f$k\_dim@f$ linear equations @f$p\_a \cdot p\_x = p\_b@f$, where.
  !!
  !! Solves the set of @f$k\_dim@f$ linear equations @f$p\_a \cdot p\_x = p\_b@f$, where
  !! @f$p\_a@f$ is a positive-definite symmetric matrix with physical dimension @f$k\_dim@f$.
  !! @f$p\_a@f$ and @f$p\_diag@f$ are input as the output of the routine {\\tt choldec}.
  !! Only the lower triangle of @f$p\_a@f$ is accessed. @f$p\_b(1:k\_dim)@f$ is input as
  !! the right-hand side vector. The solution vector is returned in @f$p\_x(1:k\_dim)@f$.
  !! @f$p\_a@f$ and @f$p\_diag@f$ are not modified and can be left in place for successive
  !! calls with different right-hand sides @f$p\_b@f$. @f$p\_b@f$ is not modified unless
  !! you identify @f$p\_b@f$ and @f$p\_x@f$ in the calling sequence, which is allowed.
  !!
  !! @par Revision History
  !!  Original version by Tobias Ruppert and Thomas Heinze, DWD (2006-11-14)
  !!
#ifdef __SX__
  SUBROUTINE solve_chol_v (istart, iend, k_dim, maxdim, p_a, p_diag, p_b, p_x)
    !
    INTEGER,  INTENT(in)   :: istart, iend    ! start and end indices
    INTEGER,  INTENT(in)   :: k_dim(:)        ! (nproma) matrix dimension
    INTEGER,  INTENT(in)   :: maxdim          ! maximum matrix dimension
    REAL(wp), INTENT(in)   :: p_a(:,:,:)      ! (nproma,k_dim,k_dim) matrix
    REAL(wp), INTENT(in)   :: p_diag(:,:)     ! (nproma,k_dim) diagonal of matrix
    
    REAL(wp), INTENT(inout):: p_b(:,:)        ! (nproma,k_dim) input vector
    
    REAL(wp), INTENT(inout):: p_x(:,:)        ! (k_dim,nproma) output vector
    
    REAL(wp)               :: z_sum(nproma) ! variable generated by
    ! matrix entries
    INTEGER :: jc, ji, jj    ! integer over dimensions
    
    !-----------------------------------------------------------------------
    ! Solve L*y = b, storing y in x.    
    DO ji = 1,maxdim
      DO jc = istart, iend
        IF (ji > k_dim(jc)) CYCLE
        z_sum(jc) = p_b(jc,ji)
      ENDDO
      DO jj = ji - 1, 1, -1
        DO jc = istart, iend
          IF (ji > k_dim(jc)) CYCLE
          z_sum(jc) = z_sum(jc) - p_a(jc,ji,jj) * p_x(jj,jc)
        ENDDO
      ENDDO
      DO jc = istart, iend
        IF (ji > k_dim(jc)) CYCLE
        p_x(ji,jc) = z_sum(jc) / p_diag(jc,ji)
      ENDDO
    ENDDO
    
    ! Solve L^T * x = y    
    DO ji = maxdim, 1, -1
      DO jc = istart, iend
        IF (ji > k_dim(jc)) CYCLE
        z_sum(jc) = p_x(ji,jc)
      ENDDO
      DO jj = ji+1, maxdim
        DO jc = istart, iend
          IF ((ji > k_dim(jc)) .OR. (jj > k_dim(jc))) CYCLE
          z_sum(jc) = z_sum(jc) - p_a(jc,jj,ji) * p_x(jj,jc)
        ENDDO
      ENDDO
      DO jc = istart, iend
        IF (ji > k_dim(jc)) CYCLE
        p_x(ji,jc) = z_sum(jc) / p_diag(jc,ji)
      ENDDO
    ENDDO
  END SUBROUTINE solve_chol_v
  
#else
  ! non-vecotrized version
  SUBROUTINE solve_chol_v (istart, iend, k_dim, p_a, p_diag, p_b, p_x)
    !
    INTEGER,  INTENT(in)   :: istart, iend    ! start and end indices
    INTEGER,  INTENT(in)   :: k_dim(:)        ! (nproma) matrix dimension
    REAL(wp), INTENT(in)   :: p_a(:,:,:)      ! (nproma,k_dim,k_dim) matrix
    REAL(wp), INTENT(in)   :: p_diag(:,:)     ! (nproma,k_dim) diagonal of matrix
    
    REAL(wp), INTENT(inout):: p_b(:,:)        ! (nproma,k_dim) input vector
    
    REAL(wp), INTENT(inout):: p_x(:,:)        ! (k_dim,nproma) output vector
    
    REAL(wp)               :: z_sum         ! variable generated by
    ! matrix entries
    INTEGER :: jc, ji, jj    ! integer over dimensions
    
    !-----------------------------------------------------------------------
    ! Solve L*y = b, storing y in x.    
    DO jc = istart, iend
      DO ji = 1, k_dim(jc)
        z_sum = p_b(jc,ji)
        DO jj = ji - 1, 1, -1
          z_sum = z_sum - p_a(jc,ji,jj) * p_x(jj,jc)
        ENDDO
        p_x(ji,jc) = z_sum / p_diag(jc,ji)
      ENDDO
      
      ! Solve L^T * x = y      
      DO ji = k_dim(jc), 1, -1
        z_sum = p_x(ji,jc)
        DO jj = ji+1, k_dim(jc)
          z_sum = z_sum - p_a(jc,jj,ji) * p_x(jj,jc)
        ENDDO
        p_x(ji,jc) = z_sum / p_diag(jc,ji)
      ENDDO
    ENDDO
    
  END SUBROUTINE solve_chol_v
  !-------------------------------------------------------------------------
#endif
  
  !-------------------------------------------------------------------------
  !>
  !!  home made routine for backsubstitution.
  !!
  !!  home made routine for backsubstitution
  !!  to be substituted by call from linear algebra library
  !!  dim: matrix dimension
  !!  a: input matrix  lu decomposed by ludec (otherwise it does not work)
  !!  b: right hand side, row pivot  exchanges to be performed using ind
  !!  ind: integer index array recording row pivot  exchanges
  !!
  SUBROUTINE solve_lu (k_dim, p_a, p_b, p_x, k_pivot)
    !
    INTEGER,  INTENT(in)    :: k_dim                   ! matrix dimension
    INTEGER,  INTENT(in)    :: k_pivot(k_dim)          ! pivoting index
    REAL(wp), INTENT(in)    :: p_a(k_dim,k_dim)        ! matrix
    
    REAL(wp), INTENT(inout) :: p_b(k_dim)              ! input vector
    
    REAL(wp), INTENT(out)   :: p_x(k_dim)              ! output vector
    
    INTEGER :: ji, jj, ip              ! integer over dimensions
    REAL(wp)                :: z_sum                   ! sum
    
    !-----------------------------------------------------------------------   
    ! forward    
    DO ji = 1, k_dim
      ip      = k_pivot(ji)
      z_sum   = p_b(ip)
      p_b(ip) = p_b(ji)
      DO jj = 1, ji-1
        z_sum = z_sum - p_a(ji,jj) * p_x(jj)
      ENDDO
      p_x(ji) = z_sum
    ENDDO
    
    ! backward    
    DO ji = k_dim, 1, -1
      z_sum = p_x(ji)
      DO jj = ji + 1, k_dim
        z_sum = z_sum - p_a(ji,jj) * p_x(jj)
      ENDDO
      p_x(ji) = z_sum / p_a(ji,ji)
    ENDDO
    
  END SUBROUTINE solve_lu
  !-------------------------------------------------------------------------
    
  
END MODULE mo_math_utility_solvers

