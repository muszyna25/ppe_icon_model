!! This algorithm computes the interpolation weights for remapping
!! with radial basis functions (RBFs).
!!
!! @note Restrictions:
!!       - Algorithm is restricted to Gaussian kernels.
!!       - The RBF stencil is formed by four lon-lat boxes. The
!!         centers of these boxes form the lon-lat area where the edge
!!         midpoint lies.
!!
!! @author F. Prill, DWD
!!
MODULE mo_remap_weights_rbf

!$  USE OMP_LIB

  USE mo_kind,                 ONLY: wp
  USE mo_parallel_config,      ONLY: nproma
  USE mo_exception,            ONLY: finish
  USE mo_communication,        ONLY: blk_no, idx_no
  USE mo_math_constants,       ONLY: pi_180
  USE mo_math_utilities,       ONLY: t_geographical_coordinates, t_cartesian_coordinates, &
    &                                gvec2cvec, gc2cc
  USE mo_math_utility_solvers, ONLY: choldec_v, solve_chol_v
  USE mo_remap_shared,         ONLY: t_grid, GRID_TYPE_REGULAR, GRID_TYPE_ICON,  &
    &                                dist_deg
  USE mo_remap_config,         ONLY: MAX_NSTENCIL_RBF, dbg_level, rbf_vec_dim
  USE mo_remap_intp,           ONLY: t_intp_data
  USE mo_remap_grid_regular,   ONLY: get_containing_dual_cell

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: prepare_interpolation_rbf_vec

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_weights_rbf')


CONTAINS

  !-------------------------------------------------------------------------   
  !>
  !! Gaussian kernel for RBF interpolation.
  !! 
  !! @f$\phi(r)=e^{-r^2}@f$
  !! 
  !! @par Revision History
  !! Developed by L. Bonaventura  (2004).
  !! 
  FUNCTION gaussi (p_x, p_scale)  RESULT (p_rbf_val)
    REAL(wp) , INTENT(in) :: p_x             ! radial distance
    REAL(wp) , INTENT(in) :: p_scale         ! scale parameter
    
    REAL(wp)              :: p_rbf_val       ! RBF value
    
    !-----------------------------------------------------------------------  
    
    p_rbf_val = p_x / p_scale
    p_rbf_val = -1._wp * p_rbf_val * p_rbf_val
    p_rbf_val = EXP(p_rbf_val)
  END FUNCTION gaussi


  !> Set RBF interpolation stencil for a single destination point.
  !
  !  The RBF stencil is formed by four lon-lat boxes. The centers of
  !  these boxes form the lon-lat area where the edge midpoint lies
  !  (in the vicinity of the poles we do a slightly different
  !  treatment).
  !
  SUBROUTINE set_rbf_stencil(p_in, grid, je, jb, intp_data_u, intp_data_v)
    TYPE (t_geographical_coordinates), INTENT(IN) :: p_in
    TYPE (t_grid),     INTENT(INOUT) :: grid         !< source grid partition
    INTEGER,           INTENT(IN)    :: je,jb        !< index of destination point
    TYPE(t_intp_data), INTENT(INOUT) :: intp_data_u  !< interpolation coefficients, u comp. (result)
    TYPE(t_intp_data), INTENT(INOUT) :: intp_data_v  !< interpolation coefficients, v comp. (result)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::set_rbf_stencil')
    INTEGER :: global_idx(4), i, g_idx, g_blk

    IF (MAX_NSTENCIL_RBF < rbf_vec_dim) &
      &   CALL finish(routine, "Allocated stencil size too small!")

    ! search for surrounding lon-lat boxes, exploit regular grid
    ! structure:
    global_idx(:) = get_containing_dual_cell(grid%regular_grid, p_in)

    ! data sites for u,v component are identical:
    intp_data_u%nidx(je,jb) = rbf_vec_dim
    intp_data_v%nidx(je,jb) = rbf_vec_dim
    DO i=1,rbf_vec_dim
      g_idx = idx_no(global_idx(i))
      g_blk = blk_no(global_idx(i))
      intp_data_u%iidx(i,je,jb) = g_idx
      intp_data_u%iblk(i,je,jb) = g_blk
      intp_data_v%iidx(i,je,jb) = g_idx
      intp_data_v%iblk(i,je,jb) = g_blk
    END DO
  END SUBROUTINE set_rbf_stencil


  !> Computes RBF interpolation weights; implementation for vector
  !  field tangential to the sphere.
  !
  !  @todo insert normalization of weights
  !
  SUBROUTINE prepare_interpolation_rbf_vec(src_grid, dst_grid, intp_data_u, intp_data_v, &
    &                                      lcompute_vt, rbf_vec_scale)
    TYPE (t_grid),     INTENT(INOUT) :: src_grid      !< source grid partition
    TYPE (t_grid),     INTENT(INOUT) :: dst_grid      !< destination grid partition
    TYPE(t_intp_data), INTENT(INOUT) :: intp_data_u   !< interpolation coefficients, u component (result)
    TYPE(t_intp_data), INTENT(INOUT) :: intp_data_v   !< interpolation coefficients, v component (result)
    LOGICAL,           INTENT(IN)    :: lcompute_vt   !< Flag: .TRUE., if the tangential wind shall be computed
    REAL(wp),          INTENT(IN)    :: rbf_vec_scale !< RBF scaling factor
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::prepare_interpolation_rbf')   
    INTEGER, PARAMETER :: N = 2*rbf_vec_dim

    INTEGER  :: jb,je,start_idx,end_idx,                     &
      &         j1_idx,j1_blk, is1, is2, j1, j2,             &
      &         kdim(nproma)
    REAL(wp) :: mat_A(nproma, N, N), rhs(nproma,N),          &
      &         e1(3,nproma,N), e2(3,nproma,N),              &
      &         z_diag(nproma, N),                           &
      &         rbf_vec_coeff(N,nproma),                     &
      &         phi, z_dist, z_norm, n_e1, n_e2, mat_R(3,3)
    TYPE (t_cartesian_coordinates)    :: cc(nproma,N),       &
      &         vec_x(nproma), cc_x(nproma), dfr, n_x
    TYPE (t_geographical_coordinates) :: gc(nproma,N), gc_x

    ! ------------------
    ! consistency checks
    ! ------------------

    ! test, if the direction of interpolation is currently supported:
    IF (.NOT. (src_grid%structure == GRID_TYPE_REGULAR) .AND.  &
      &       (dst_grid%structure == GRID_TYPE_ICON)) THEN
      CALL finish(routine, "Not yet implemented!")
    END IF

    IF (dbg_level >= 2) WRITE (0,*) "# computing RBF interpolation weights."

    ! -------------------
    ! stencil computation
    ! -------------------

    ! get nearest neighbor data point for each edge midpoint of the
    ! triangular mesh, set RBF stencil indices
!$OMP PARALLEL PRIVATE(jb,je, start_idx,end_idx)
!$OMP DO
    DO jb = 1,dst_grid%p_patch%nblks_e
      start_idx = 1
      end_idx   = nproma
      IF (jb == dst_grid%p_patch%nblks_e) end_idx = dst_grid%p_patch%npromz_e
      DO je = start_idx, end_idx
        CALL set_rbf_stencil(dst_grid%p_patch%edges%center(je,jb), src_grid, je, jb, &
          &                  intp_data_u, intp_data_v)
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    ! --------------------------
    ! computation of RBF weights
    ! --------------------------

    kdim(:) = N

    ! blocks loop: set of destination points (= edge midpoints):
!$OMP PARALLEL PRIVATE(jb,je, j1,j2, is1,is2, start_idx,end_idx,mat_A,rhs,gc_x,   &
!$OMP                  cc_x,mat_R,n_x,vec_x,j1_idx,j1_blk,gc,cc,z_norm,e1,e2,dfr, &
!$OMP                  z_dist,phi,n_e1,n_e2,z_diag,rbf_vec_coeff)
!$OMP DO
    DO jb = 1,dst_grid%p_patch%nblks_e
      start_idx = 1
      end_idx   = nproma
      IF (jb == dst_grid%p_patch%nblks_e) end_idx = dst_grid%p_patch%npromz_e

      mat_A(:,:,:) = 0._wp
      rhs(:,:)     = 0._wp

      ! -----------------------------------------------------
      ! precomputations for destination point (edge midpoint)
      ! -----------------------------------------------------

      DO je = start_idx, end_idx
        ! evaluation point x (edge midpoint)
        gc_x%lon = dst_grid%p_patch%edges%center(je,jb)%lon * pi_180
        gc_x%lat = dst_grid%p_patch%edges%center(je,jb)%lat * pi_180
        cc_x(je) = gc2cc(gc_x)
        
        ! precompute tangential vector (if required)
        IF (lcompute_vt) THEN
          ! indices loop: set of destination points (= edge midpoints):
          
          ! rotation matrix: "vn" -> "vt" (for optional computation of "vt")
          mat_R(1,1:3) = (/ cc_x(je)%x(1)*cc_x(je)%x(1),                 &
            &               cc_x(je)%x(1)*cc_x(je)%x(2) - cc_x(je)%x(3), &
            &               cc_x(je)%x(1)*cc_x(je)%x(3) + cc_x(je)%x(2) /) 
          mat_R(2,1:3) = (/ cc_x(je)%x(1)*cc_x(je)%x(2) + cc_x(je)%x(3), &
            &               cc_x(je)%x(2)*cc_x(je)%x(2),                 &
            &               cc_x(je)%x(2)*cc_x(je)%x(3) - cc_x(je)%x(1) /) 
          mat_R(3,1:3) = (/ cc_x(je)%x(1)*cc_x(je)%x(3) - cc_x(je)%x(2), &
            &               cc_x(je)%x(2)*cc_x(je)%x(3) + cc_x(je)%x(1), &
            &               cc_x(je)%x(3)*cc_x(je)%x(3)  /) 
          
          n_x = dst_grid%p_patch%edges%primal_cart_normal(je,jb) 

!CDIR BEGIN EXPAND=3
          vec_x(je)%x(1) = DOT_PRODUCT(mat_R(1,1:3), n_x%x(1:3))
          vec_x(je)%x(2) = DOT_PRODUCT(mat_R(2,1:3), n_x%x(1:3))
          vec_x(je)%x(3) = DOT_PRODUCT(mat_R(3,1:3), n_x%x(1:3))
!CDIR END
        ELSE
          vec_x(je) = dst_grid%p_patch%edges%primal_cart_normal(je,jb)
        END IF
      END DO

      ! --------------------------------------------
      ! compute right hand side of the linear system
      ! --------------------------------------------

      DO j1 = 1,rbf_vec_dim
        ! indices loop: set of destination points (= edge midpoints):
        DO je = start_idx, end_idx

          ! data site j1
          ! ------------
          
          ! get index/block of data site j1
          ! (data sites for u,v component are identical)
          j1_idx = intp_data_u%iidx(j1,je,jb)
          j1_blk = intp_data_u%iblk(j1,je,jb)
          ! get lon/lat coordinates on unit sphere
          gc(je,j1) = src_grid%p_patch%cells%center(j1_idx,j1_blk)
          gc(je,j1)%lon = gc(je,j1)%lon * pi_180
          gc(je,j1)%lat = gc(je,j1)%lat * pi_180
          ! get Cartesian coordinates
          cc(je,j1) = gc2cc(gc(je,j1))

          ! compute zonal unit vector:
          CALL gvec2cvec(1._wp,0._wp, gc(je,j1)%lon, gc(je,j1)%lat, &
            &            e1(1,je,j1), e1(2,je,j1), e1(3,je,j1))
          z_norm   = SQRT( DOT_PRODUCT( e1(:,je,j1), e1(:,je,j1)) )
          e1(:,je,j1) = 1._wp/z_norm * e1(:,je,j1)
          ! compute meridional unit vector:
          CALL gvec2cvec(0._wp,1._wp, gc(je,j1)%lon, gc(je,j1)%lat, &
            &            e2(1,je,j1), e2(2,je,j1), e2(3,je,j1))
          z_norm   = SQRT( DOT_PRODUCT(e2(:,je,j1),e2(:,je,j1)) )
          e2(:,je,j1) = 1._wp/z_norm * e2(:,je,j1)

          ! evaluate scalar RBF
          ! -------------------

          dfr%x(:) = cc(je,j1)%x(:) - cc_x(je)%x(:)
          z_dist = SQRT( DOT_PRODUCT(dfr%x, dfr%x) )
          phi = gaussi(z_dist, rbf_vec_scale)

          ! compute RHS
          ! -----------
          !
          ! ( vector (n^T * B), where B is the output matrix in R^{3 x 2N} )
          n_e1 = DOT_PRODUCT(vec_x(je)%x(:), e1(:,je,j1))
          n_e2 = DOT_PRODUCT(vec_x(je)%x(:), e2(:,je,j1))

          is1 = 2*j1-1
          rhs(je,is1:(is1+1)) = (/ phi*n_e1, phi*n_e2 /)

        END DO ! je
      END DO ! j1

      ! ----------------------------------
      ! compute RBF interpolation matrix A 
      ! ----------------------------------
      !
      ! ( symmetric matrix A in R^{2N x 2N}, where we have N data
      !   sites and components u,v at each site )
      !
      ! note: it is sufficient to compute only half of the entries,
      ! since only the **upper** triangle of A is accessed.
      DO j1 = 1,rbf_vec_dim
        DO j2 = j1,rbf_vec_dim
          is1 = 2*j1 - 1
          is2 = 2*j2 - 1
          ! indices loop: set of destination points (= edge midpoints):
          DO je = start_idx, end_idx
            
            ! evaluate scalar RBF:
            dfr%x(:) = cc(je,j1)%x(:) - cc(je,j2)%x(:)
            z_dist = SQRT( DOT_PRODUCT(dfr%x, dfr%x) )
            phi = gaussi(z_dist, rbf_vec_scale)
            ! fill 2x2 subblock in the A matrix:
            mat_A(je, is1,     is2:(is2+1) ) = &
              & (/ phi*DOT_PRODUCT(e1(:,je,j1),e1(:,je,j2)), phi*DOT_PRODUCT(e1(:,je,j1),e2(:,je,j2)) /)
            mat_A(je, (is1+1), is2:(is2+1) ) = &
              & (/ phi*DOT_PRODUCT(e2(:,je,j1),e1(:,je,j2)), phi*DOT_PRODUCT(e2(:,je,j1),e2(:,je,j2)) /)
          END DO ! je
        END DO ! j2
      END DO ! j1

      ! ---------------------------------------------
      ! solve linear system by Cholesky decomposition
      ! ---------------------------------------------

#ifdef __SX__
!CDIR NOIEXPAND
      CALL choldec_v(start_idx,end_idx, kdim,N, mat_A, z_diag)
!CDIR NOIEXPAND
      CALL solve_chol_v(start_idx,end_idx, kdim,N, mat_A, z_diag, rhs, rbf_vec_coeff)
#else
      CALL choldec_v(start_idx,end_idx, kdim, mat_A, z_diag)
      CALL solve_chol_v(start_idx,end_idx, kdim, mat_A, z_diag, rhs, rbf_vec_coeff)
#endif

      ! ---------------------------------------------
      ! distribute coeffients (interpolation weights)
      ! ---------------------------------------------

      DO j1 = 1,rbf_vec_dim
        is1 = 2*j1 - 1
        ! indices loop: set of destination points (= edge midpoints):
        DO je = start_idx, end_idx
          intp_data_u%wgt(j1,je,jb) = rbf_vec_coeff(is1,  je)
          intp_data_v%wgt(j1,je,jb) = rbf_vec_coeff(is1+1,je)
        END DO ! je
      END DO ! j1
    END DO ! jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE prepare_interpolation_rbf_vec

END MODULE mo_remap_weights_rbf
