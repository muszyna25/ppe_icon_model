#if (defined(_OPENMP) && defined(OCE_SOLVE_OMP))
#include "omp_definitions.inc"
#endif

MODULE mo_surface_height_lhs

! contains extended lhs-matrix-generator type
! provides surface height lhs for free ocean surface solve

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: finish
  USE mo_ocean_solve_lhs_type, ONLY: t_lhs_agen
  USE mo_ocean_solve_aux, ONLY: t_destructible
  USE mo_ocean_types, ONLY: t_solverCoeff_singlePrecision, t_operator_coeff
  USE mo_grid_subset, ONLY: t_subset_range, get_index_range
  USE mo_ocean_nml, ONLY: select_lhs, select_lhs_matrix, ab_gam, ab_beta, &
    & l_edge_based, l_lhs_direct
  USE mo_communication, ONLY: exchange_data
  USE mo_parallel_config, ONLY: nproma
  USE mo_impl_constants, ONLY: sea_boundary
  USE mo_run_config, ONLY: dtime, debug_check_level
  USE mo_physical_constants, ONLY: grav
  USE mo_ocean_math_operators, ONLY: &
    & grad_fd_norm_oce_2d_3d, div_oce_2D_onTriangles_onBlock, &
    & div_oce_2D_general_onBlock
  USE mo_scalar_product, ONLY: map_edges2edges_viacell_3d_const_z
  USE mo_model_domain, ONLY: t_patch_3d, t_patch
  USE mo_mpi, ONLY: my_process_is_mpi_parallel

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_surface_height_lhs, lhs_surface_height_ptr

  TYPE, EXTENDS(t_lhs_agen) :: t_surface_height_lhs
    PRIVATE
    TYPE(t_patch_3d), POINTER :: patch_3d => NULL()
    TYPE(t_patch), POINTER :: patch_2d => NULL()
    REAL(wp), POINTER :: thickness_e_wp(:,:)
    TYPE(t_operator_coeff), POINTER :: op_coeffs_wp => NULL()
    TYPE(t_solverCoeff_singlePrecision), POINTER :: op_coeffs_sp => NULL()
    REAL(wp), ALLOCATABLE, DIMENSION(:,:), PRIVATE :: z_grad_h_wp, z_e_wp
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: z_grad_h_wp, z_e_wp
#endif
  CONTAINS
    PROCEDURE :: lhs_wp => lhs_surface_height_wp
    PROCEDURE :: construct => lhs_surface_height_construct
    PROCEDURE :: destruct => lhs_surface_height_destruct
    PROCEDURE, PRIVATE :: internal_matrix_wp => lhs_surface_height_ab_mim_matrix_wp
    PROCEDURE, PRIVATE :: internal_wp => lhs_surface_height_ab_mim_wp
    PROCEDURE :: lhs_matrix_shortcut => lhs_surface_height_ab_mim_matrix_shortcut
  END TYPE t_surface_height_lhs

  INTEGER, PARAMETER :: topLevel = 1

CONTAINS

! returns pointer to t_primal_flip_flop_lhs object, if provided a matching
! object of corresponding abstract type
  FUNCTION lhs_surface_height_ptr(this) RESULT(this_ptr)
    CLASS(t_destructible), INTENT(IN), TARGET :: this
    CLASS(t_surface_height_lhs), POINTER :: this_ptr

    SELECT TYPE (this)
    CLASS IS (t_surface_height_lhs)
      this_ptr => this
    CLASS DEFAULT
      NULLIFY(this_ptr)
      CALL finish("surface_height_lhs_ptr", "not correct type!")
    END SELECT
  END FUNCTION lhs_surface_height_ptr

!init generator object
  SUBROUTINE lhs_surface_height_construct(this, patch_3d, thick_e, &
      & op_coeffs_wp, op_coeffs_sp)
    CLASS(t_surface_height_lhs), INTENT(INOUT) :: this
    TYPE(t_patch_3d), POINTER, INTENT(IN) :: patch_3d
    REAL(wp), POINTER, INTENT(IN) :: thick_e(:,:)
    TYPE(t_operator_coeff), TARGET, INTENT(IN) :: op_coeffs_wp
    TYPE(t_solverCoeff_singlePrecision), TARGET, INTENT(IN) :: op_coeffs_sp

    CALL this%destruct()
    this%patch_3d => patch_3d
    this%patch_2d => patch_3d%p_patch_2d(1)
    this%thickness_e_wp => thick_e
    this%op_coeffs_wp => op_coeffs_wp
    this%op_coeffs_sp => op_coeffs_sp
    this%is_const = .false.
    this%use_shortcut = (select_lhs .GT. select_lhs_matrix .AND. select_lhs .LE. select_lhs_matrix + 1)
    IF (this%patch_2d%cells%max_connectivity .NE. 3 .AND. .NOT.l_lhs_direct) &
      & CALL finish("t_surface_height_lhs::lhs_surface_height_construct", &
      &  "internal matrix implementation only works with triangular grids!")
    ALLOCATE(this%is_init(1))
  END SUBROUTINE lhs_surface_height_construct

! interface routine clear object internals
  SUBROUTINE lhs_surface_height_destruct(this)
    CLASS(t_surface_height_lhs), INTENT(INOUT) :: this

    NULLIFY(this%patch_3d, this%patch_2d, this%thickness_e_wp)
    NULLIFY(this%op_coeffs_wp, this%op_coeffs_sp)
    IF (ALLOCATED(this%z_e_wp)) DEALLOCATE(this%z_e_wp, this%z_grad_h_wp)
    IF (ALLOCATED(this%is_init)) DEALLOCATE(this%is_init)
  END SUBROUTINE lhs_surface_height_destruct

! interface routine for the left hand side computation
  SUBROUTINE lhs_surface_height_wp(this, x, ax)
    CLASS(t_surface_height_lhs), INTENT(INOUT) :: this
    REAL(wp), INTENT(IN) :: x(:,:)
    REAL(wp), INTENT(OUT) :: ax(:,:)

    IF (this%use_shortcut) &
      & CALL finish("t_surface_height_lhs::lhs_surface_height_wp", &
        & "should not be here because of shortcut!")
! alloc internal arrays, if needed to
! if using operator implementation of lhs
    IF (select_lhs .NE. select_lhs_matrix) THEN
      IF (ALLOCATED(this%z_e_wp)) THEN
        IF (nproma .NE. SIZE(this%z_e_wp, 1) .OR. &
            & this%patch_2d%nblks_e .NE. SIZE(this%z_e_wp, 2)) &
            & DEALLOCATE(this%z_e_wp, this%z_grad_h_wp)
      END IF
      IF (.NOT.ALLOCATED(this%z_e_wp)) &
        & ALLOCATE(this%z_e_wp(nproma, this%patch_2d%nblks_e), &
            & this%z_grad_h_wp(nproma, this%patch_2d%nblks_e))
    END IF
! call approriate backend depending on select_lhs choice
    IF (select_lhs == select_lhs_matrix) THEN
      CALL this%internal_matrix_wp(x, ax)
    ELSE
      CALL this%internal_wp(x, ax)
    ENDIF
  END SUBROUTINE lhs_surface_height_wp

! internal backend routine to compute surface height lhs -- "matrix" implementation
  SUBROUTINE lhs_surface_height_ab_mim_matrix_wp(this, x, lhs)
    CLASS(t_surface_height_lhs), INTENT(INOUT) :: this
    REAL(wp), INTENT(IN), CONTIGUOUS :: x(:,:)
    REAL(wp), INTENT(OUT), CONTIGUOUS :: lhs(:,:)
    INTEGER :: start_index, end_index, jc, blkNo, ico
    TYPE(t_subset_range), POINTER :: cells_in_domain
    REAL(wp), POINTER, DIMENSION(:,:,:), CONTIGUOUS :: lhs_coeffs
    REAL(wp) :: xco(9)
    INTEGER, POINTER, DIMENSION(:,:,:), CONTIGUOUS :: idx, blk

    cells_in_domain => this%patch_2D%cells%in_domain
    lhs_coeffs => this%op_coeffs_wp%lhs_all
    idx => this%op_coeffs_wp%lhs_CellToCell_index
    blk => this%op_coeffs_wp%lhs_CellToCell_block
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, xco, ico) ICON_OMP_DEFAULT_SCHEDULE
    DO blkNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blkNo, start_index, end_index)
      lhs(:,blkNo) = 0.0_wp
      DO jc = start_index, end_index
        IF(.NOT.(this%patch_3d%lsm_c(jc,1,blkNo) > sea_boundary)) THEN
          FORALL(ico = 1:9) xco(ico) = x(idx(ico, jc, blkNo), blk(ico, jc, blkNo))
          lhs(jc,blkNo) = x(jc,blkNo) * lhs_coeffs(0, jc, blkNo) + &
            & SUM(xco(:) * lhs_coeffs(1:9, jc, blkNo))
        END IF
      END DO
    END DO ! blkNo
!ICON_OMP_END_PARALLEL_DO
    IF (debug_check_level > 20) THEN
      DO blkNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blkNo, start_index, end_index)
        DO jc = start_index, end_index
          IF(this%patch_3d%lsm_c(jc,1,blkNo) > sea_boundary) THEN
            IF (lhs(jc,blkNo) /= 0.0_wp) &
              & CALL finish("lhs_surface_height_ab_mim", "lhs(jc,blkNo) /= 0 on land")
          ENDIF
        END DO
      END DO
    ENDIF
  END SUBROUTINE lhs_surface_height_ab_mim_matrix_wp

! internal backend routine to compute surface height lhs -- "operator" implementation
  SUBROUTINE lhs_surface_height_ab_mim_wp(this, x, lhs)
    CLASS(t_surface_height_lhs), INTENT(INOUT) :: this
    REAL(wp), INTENT(IN), CONTIGUOUS :: x(:,:)
    REAL(wp), INTENT(OUT), CONTIGUOUS :: lhs(:,:)
    REAL(wp) :: gdt2_inv, gam_times_beta
    INTEGER :: start_index, end_index, jc, blkNo, je
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain

    cells_in_domain => this%patch_2D%cells%in_domain
    edges_in_domain => this%patch_2D%edges%in_domain
    gdt2_inv = REAL(1.0_wp / (grav*(dtime)**2),wp)
    gam_times_beta = REAL(ab_gam * ab_beta, wp)
    lhs(1:nproma,cells_in_domain%end_block:) = 0.0_wp
    IF (l_edge_based) THEN
      !Step 1) Calculate gradient of iterated height.
      CALL grad_fd_norm_oce_2d_3d( x, this%patch_2D, this%op_coeffs_wp%grad_coeff(:,1,:), &
        & this%z_grad_h_wp(:,:), subset_range=this%patch_2D%edges%gradIsCalculable)
      DO blkNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blkNo, start_index, end_index)
        DO je = start_index, end_index
          this%z_e_wp(je,blkNo) = this%z_grad_h_wp(je,blkNo) * this%thickness_e_wp(je,blkNo)
        END DO
      END DO
    ELSE  !IF(.NOT.l_edge_based)THEN
      CALL grad_fd_norm_oce_2d_3d( x, this%patch_2D, this%op_coeffs_wp%grad_coeff(:,1,:), &
        & this%z_grad_h_wp(:,:), subset_range=this%patch_2D%edges%gradIsCalculable)
      IF (this%patch_2d%cells%max_connectivity /= 3 .AND. my_process_is_mpi_parallel()) &
        & CALL exchange_data(this%patch_2D%comm_pat_e, this%z_grad_h_wp)
      CALL map_edges2edges_viacell_3d_const_z( this%patch_3d, &
        & this%z_grad_h_wp(:,:), this%op_coeffs_wp, this%z_e_wp(:,:))
    END IF ! l_edge_based
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO blkNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blkNo, start_index, end_index)
      IF (this%patch_2d%cells%max_connectivity .EQ. 3) THEN
        CALL div_oce_2D_onTriangles_onBlock(this%z_e_wp, this%patch_2D, &
          & this%op_coeffs_wp%div_coeff, lhs(:,blkNo), level=topLevel, &
          & blockNo=blkNo, start_index=start_index, end_index=end_index)
      ELSE
        CALL div_oce_2D_general_onBlock(this%z_e_wp, this%patch_2D, &
          & this%op_coeffs_wp%div_coeff, lhs(:,blkNo), level=topLevel, &
          & blockNo=blkNo, start_index=start_index, end_index=end_index)
      END IF
      !Step 4) Finalize LHS calculations
      DO jc = start_index, end_index
        lhs(jc,blkNo) = x(jc,blkNo) * gdt2_inv - gam_times_beta * lhs(jc,blkNo)
      END DO
    END DO ! blkNo
!ICON_OMP_END_PARALLEL_DO
    IF (debug_check_level > 20) THEN
      DO blkNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blkNo, start_index, end_index)
        DO jc = start_index, end_index
          IF(this%patch_3d%lsm_c(jc,1,blkNo) > sea_boundary) THEN
            IF (lhs(jc,blkNo) /= 0.0_wp) &
              & CALL finish("lhs_surface_height_ab_mim",&
                  & "lhs(jc,blkNo) /= 0 on land")
          END IF
        END DO
      END DO
    ENDIF
  END SUBROUTINE lhs_surface_height_ab_mim_wp

  SUBROUTINE lhs_surface_height_ab_mim_matrix_shortcut(this, idx, blk, coeff)
    CLASS(t_surface_height_lhs), INTENT(INOUT) :: this
    INTEGER, INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:,:) :: idx, blk
    REAL(KIND=wp), INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:,:) :: coeff
    INTEGER :: nidx, nblk, nnz, iidx, iblk, inz
    INTEGER, POINTER, DIMENSION(:,:,:), CONTIGUOUS :: opc_idx, opc_blk
    REAL(wp), POINTER, DIMENSION(:,:,:), CONTIGUOUS :: opc_coeff

    IF (.NOT.this%use_shortcut) &
      & CALL finish( &
        & "t_surface_height_lhs::lhs_surface_height_ab_mim_matrix_shortcut", &
        & "wrong turn!")    
    nidx = SIZE(coeff, 1)
    nblk = SIZE(coeff, 2)
    nnz = 10
    opc_coeff => this%op_coeffs_wp%lhs_all
    IF (.NOT.ALLOCATED(idx)) THEN
      DEALLOCATE(coeff)
      ALLOCATE(idx(nidx, nblk, nnz), blk(nidx, nblk, nnz), &
        & coeff(nidx, nblk, nnz))
      opc_idx => this%op_coeffs_wp%lhs_CellToCell_index
      opc_blk => this%op_coeffs_wp%lhs_CellToCell_block
!ICON_OMP PARALLEL PRIVATE(iblk, iidx, inz)
!ICON_OMP DO SCHEDULE(GUIDED)
      DO iblk = 1, nblk
        idx(:, iblk, 1) = (/(iidx, iidx = 1, nidx)/)
        blk(:, iblk, 1) = iblk
        coeff(:, iblk, 1) = opc_coeff(0, :, iblk)
      END DO
!ICON_OMP END DO NOWAIT
!ICON_OMP DO SCHEDULE(GUIDED)
      DO inz = 2, nnz
        DO iblk = 1, nblk
          DO iidx = 1, nidx
            IF (ABS(opc_coeff(inz - 1, iidx, iblk)) .EQ. 0._wp) THEN
              idx(iidx, iblk, inz) = idx(iidx, iblk, 1)
              blk(iidx, iblk, inz) = blk(iidx, iblk, 1)
              coeff(iidx, iblk, inz) = 0._wp
            ELSE
              idx(iidx, iblk, inz) = opc_idx(inz - 1, iidx, iblk)
              blk(iidx, iblk, inz) = opc_blk(inz - 1, iidx, iblk)
              coeff(iidx, iblk, inz) = opc_coeff(inz - 1, iidx, iblk)
            END IF
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
!ICON_OMP END PARALLEL
    ELSE
!ICON_OMP PARALLEL DO PRIVATE(iblk)
      DO inz = 1, nnz
        DO iblk = 1, nblk
          coeff(:, iblk, inz) = opc_coeff(inz - 1, :, iblk)
        END DO
      END DO
!ICON_OMP END PARALLEL DO
    END IF
  END SUBROUTINE lhs_surface_height_ab_mim_matrix_shortcut

END MODULE mo_surface_height_lhs
