!! Subroutines for MPI communication.
!!
!! @author F. Prill, DWD
!!
MODULE mo_remap_sync

#ifdef __ICON__
  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_communication,      ONLY: idx_no, blk_no
#else
  USE mo_utilities,    ONLY: wp, SUCCESS, nproma, idx_no,             &
    &                        blk_no, finish
#endif

  USE mo_mpi,          ONLY: p_real_dp, p_comm_work, p_n_work, p_int, &
    &                        get_my_mpi_work_id, p_gather, p_gatherv, &
    &                        p_scatterv, p_wait, p_isend, p_recv,     &
    &                        p_clear_request
  USE mo_remap_config, ONLY: dbg_level
  USE mo_remap_shared, ONLY: t_grid
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: allocate_gather_c, allocate_gather_e
  PUBLIC :: finalize_gather
  PUBLIC :: scatter_field2D 
  PUBLIC :: gather_field2D, gather_field3D
  PUBLIC :: t_gather

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_sync')

  !> Data structure for local-to-global communication of distributed
  !> field data
  !
  TYPE t_gather
    INTEGER,  ALLOCATABLE :: global_idx(:), &   !< mapping: (global, PE-ordered array) -> (global array)
      &                      local_size(:), &   !< local field size
      &                      displs(:)          !< offset in global, PE-ordered array
    INTEGER :: local_dim = -1
    INTEGER :: total_dim = -1
    INTEGER :: rank0     = -1
    INTEGER :: nblks     = -1
    INTEGER :: npromz    = -1
  END TYPE t_gather

  INTERFACE gather_field2D
    MODULE PROCEDURE gather_field2D_real
    MODULE PROCEDURE gather_field2D_int
  END INTERFACE

CONTAINS

  !> Prepare communication data structures for cell-based grid.
  !
  SUBROUTINE allocate_gather_c(grid, rank0, gather_c)
    TYPE (t_grid),  INTENT(IN)      :: grid             !< local grid partition
    INTEGER,        INTENT(IN)      :: rank0            !< root PE, where data is collected.
    TYPE(t_gather), INTENT(INOUT)   :: gather_c         !< communication pattern
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::allocate_gather_c'
    INTEGER :: ierrstat, i, this_pe

    IF (dbg_level >= 5) WRITE (0,*) "# allocate communication pattern"
    this_pe = get_my_mpi_work_id()
    gather_c%rank0     = rank0
    gather_c%local_dim = grid%p_patch%n_patch_cells

#if !defined(NOMPI)
    ALLOCATE(gather_c%displs(p_n_work), gather_c%local_size(p_n_work), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! communicate local field sizes:
    gather_c%local_size(:) = 0
    CALL p_gather(gather_c%local_dim, gather_c%local_size, gather_c%rank0, p_comm_work)

    IF (this_pe == rank0) THEN
      gather_c%total_dim = SUM(gather_c%local_size)
      gather_c%nblks     = blk_no(grid%p_patch%n_patch_cells_g)
      gather_c%npromz    = grid%p_patch%n_patch_cells_g - nproma*(gather_c%nblks-1)
    ELSE
      gather_c%total_dim = 0
      gather_c%nblks     = 0
      gather_c%npromz    = 0
    END IF
    ALLOCATE(gather_c%global_idx(gather_c%total_dim), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! compute array displacements (offsets):
    gather_c%displs(:) = gather_c%local_size(:)
    DO i=2,p_n_work
      gather_c%displs(i) = gather_c%displs(i) + gather_c%displs(i-1)
    END DO
    DO i=p_n_work,2,-1
      gather_c%displs(i) = gather_c%displs(i-1)
    END DO
    gather_c%displs(1) = 0

    ! call variable MPI gather operations:
    ! communicate local-to-global index mappings:
    CALL p_gatherv(grid%p_patch%cells%glb_index, gather_c%local_dim,          &
      &            gather_c%global_idx, gather_c%local_size, gather_c%displs, &
      &            gather_c%rank0, p_comm_work)
#else
    ! proper initialization for non-MPI compilation:
    ALLOCATE(gather_c%global_idx(grid%p_patch%n_patch_cells),  &
      &      gather_c%local_size(1), gather_c%displs(1), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    gather_c%total_dim     = gather_c%local_dim
    gather_c%local_size(1) = gather_c%local_dim
    gather_c%nblks         = blk_no(grid%p_patch%n_patch_cells_g)
    gather_c%npromz        = grid%p_patch%n_patch_cells_g - nproma*(gather_c%nblks-1)
    gather_c%displs(1)     = 0
    gather_c%global_idx(:) = grid%p_patch%cells%glb_index(:)
#endif
  END SUBROUTINE allocate_gather_c


  !> Prepare communication data structures for edge-based grid.
  !
  SUBROUTINE allocate_gather_e(grid, rank0, gather_e)
    TYPE (t_grid),  INTENT(IN)      :: grid             !< local grid partition
    INTEGER,        INTENT(IN)      :: rank0            !< root PE, where data is collected.
    TYPE(t_gather), INTENT(INOUT)   :: gather_e         !< communication pattern
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::allocate_gather_e'
    INTEGER :: ierrstat, i, this_pe

    IF (dbg_level >= 5) WRITE (0,*) "# allocate communication pattern"
    this_pe = get_my_mpi_work_id()
    gather_e%rank0     = rank0
    gather_e%local_dim = grid%p_patch%n_patch_edges

#if !defined(NOMPI)
    ALLOCATE(gather_e%displs(p_n_work), gather_e%local_size(p_n_work), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! communicate local field sizes:
    gather_e%local_size(:) = 0
    CALL p_gather(gather_e%local_dim, gather_e%local_size, gather_e%rank0, p_comm_work)

    IF (this_pe == rank0) THEN
      gather_e%total_dim = SUM(gather_e%local_size)
      gather_e%nblks     = blk_no(grid%p_patch%n_patch_edges_g)
      gather_e%npromz    = grid%p_patch%n_patch_edges_g - nproma*(gather_e%nblks-1)
    ELSE
      gather_e%total_dim = 0
      gather_e%nblks     = 0
      gather_e%npromz    = 0
    END IF
    ALLOCATE(gather_e%global_idx(gather_e%total_dim), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! compute array displacements (offsets):
    gather_e%displs(:) = gather_e%local_size(:)
    DO i=2,p_n_work
      gather_e%displs(i) = gather_e%displs(i) + gather_e%displs(i-1)
    END DO
    DO i=p_n_work,2,-1
      gather_e%displs(i) = gather_e%displs(i-1)
    END DO
    gather_e%displs(1) = 0

    ! call variable MPI gather operations:
    ! communicate local-to-global index mappings:
    CALL p_gatherv(grid%p_patch%edges%glb_index, gather_e%local_dim,          &
      &            gather_e%global_idx, gather_e%local_size, gather_e%displs, &
      &            gather_e%rank0, p_comm_work)
#else
    ! proper initialization for non-MPI compilation:
    ALLOCATE(gather_e%global_idx(grid%p_patch%n_patch_edges),  &
      &      gather_e%local_size(1), gather_e%displs(1), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    gather_e%total_dim     = gather_e%local_dim
    gather_e%local_size(1) = gather_e%local_dim
    gather_e%nblks         = blk_no(grid%p_patch%n_patch_edges_g)
    gather_e%npromz        = grid%p_patch%n_patch_edges_g - nproma*(gather_e%nblks-1)
    gather_e%displs(1)     = 0
    gather_e%global_idx(:) = grid%p_patch%edges%glb_index(:)
#endif
  END SUBROUTINE allocate_gather_e


  !> Free communication data structures.
  !
  RECURSIVE SUBROUTINE finalize_gather(gather, &
    &                                  opt_gather2, opt_gather3, opt_gather4)
    TYPE(t_gather), INTENT(INOUT) :: gather
    TYPE(t_gather), INTENT(INOUT), OPTIONAL :: opt_gather2, opt_gather3, opt_gather4
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::finalize_gather')
    INTEGER :: ierrstat

    ! the following recursion makes this more readable on the caller
    ! side (we can call this SR with more than one object as an
    ! argument):
    IF (PRESENT(opt_gather4)) THEN
      CALL finalize_gather(opt_gather4)
      CALL finalize_gather(gather, opt_gather2, opt_gather3)
    ELSE IF (PRESENT(opt_gather3)) THEN
      CALL finalize_gather(opt_gather3, opt_gather4)
      CALL finalize_gather(gather, opt_gather2)
    ELSE IF (PRESENT(opt_gather2)) THEN
      CALL finalize_gather(opt_gather2)
      CALL finalize_gather(gather)
    ELSE
      ! clean up
      IF (dbg_level >= 10) WRITE (0,*) "# finalize gather/scatter pattern."
      IF (ALLOCATED(gather%displs))  DEALLOCATE(gather%displs, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
      IF (ALLOCATED(gather%global_idx))  DEALLOCATE(gather%global_idx, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
      IF (ALLOCATED(gather%local_size))  DEALLOCATE(gather%local_size, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
  END SUBROUTINE finalize_gather


  !> Scatter a global 2D field data to working PEs.
  !
  SUBROUTINE scatter_field2D(gather, field_in, field_out)
    REAL(wp),         INTENT(IN)    :: field_in(:,:)    !< input field (only required for PE rank0)
    REAL(wp),         INTENT(INOUT) :: field_out(:,:)   !< local output field (size may differ on PEs)
    TYPE(t_gather),   INTENT(IN)    :: gather           !< communication pattern
#if !defined(NOMPI)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::scatter_field2D')
    INTEGER  :: ierrstat, i, jc,jb
    REAL(wp), ALLOCATABLE :: send_buf(:)

    IF (dbg_level >= 11) &
      & WRITE (0,*) "# scatter data from PE ", gather%rank0
    ! allocate temporary data structures:
    ALLOCATE(send_buf(gather%total_dim), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! map values to the right positions in the send buffer:
    IF (get_my_mpi_work_id() == gather%rank0) THEN
!$OMP PARALLEL DO PRIVATE(jc,jb)
      DO i=1,gather%total_dim
        jc = idx_no(gather%global_idx(i))
        jb = blk_no(gather%global_idx(i))
        send_buf(i) = field_in(jc,jb)
      END DO
!$OMP END PARALLEL DO
    END IF ! if rank0

    ! communicate field data:
    CALL p_scatterv(send_buf, gather%local_size, gather%displs,   &
      &             field_out, gather%local_dim,                  &
      &             gather%rank0, p_comm_work)

    ! clean up
    DEALLOCATE(send_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#else
    field_out(:,:) = field_in(:,:)
#endif
  END SUBROUTINE scatter_field2D


  !> Gather distributed 2D field data from working PEs in a global
  !  field on rank0.
  !
  SUBROUTINE gather_field2D_real(gather, field_in, field_out)
    REAL(wp),          INTENT(IN)    :: field_in(:,:)    !< input field (size may differ on PEs)
    REAL(wp),          INTENT(INOUT) :: field_out(:,:)   !< global output field
    TYPE(t_gather),    INTENT(IN)    :: gather           !< communication pattern
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::gather_field2D')
#if !defined(NOMPI)
    INTEGER  :: ierrstat, i, jc,jb
    REAL(wp), ALLOCATABLE :: recv_buf(:)

    IF (dbg_level >= 11) WRITE (0,*) "# gather data on PE ", gather%rank0
    ! allocate temporary data structures:
    ALLOCATE(recv_buf(gather%total_dim), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! communicate field data:
    CALL p_gatherv(field_in, gather%local_dim,                 &
      &            recv_buf, gather%local_size, gather%displs, &
      &            gather%rank0, p_comm_work)

    ! map received values to the right (global) positions:
!$OMP PARALLEL DO PRIVATE(jc,jb)
    DO i=1,gather%total_dim
      jc = idx_no(gather%global_idx(i))
      jb = blk_no(gather%global_idx(i))
      field_out(jc,jb) = recv_buf(i)
    END DO
!$OMP END PARALLEL DO

    ! clean up
    DEALLOCATE(recv_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#else
    IF (SIZE(field_out,1) /= SIZE(field_in,1)) CALL finish(routine, "Mismatch dimension 1")
    IF (SIZE(field_out,2) /= SIZE(field_in,2)) CALL finish(routine, "Mismatch dimension 2")
    field_out(:,:) = field_in(:,:)
#endif
  END SUBROUTINE gather_field2D_real


  !> Gather distributed 2D field data from working PEs in a global
  !  field on rank0: INTEGER variant
  !
  SUBROUTINE gather_field2D_int(gather, field_in, field_out)
    INTEGER,           INTENT(IN)    :: field_in(:,:)    !< input field (size may differ on PEs)
    INTEGER,           INTENT(INOUT) :: field_out(:,:)   !< global output field
    TYPE(t_gather),    INTENT(IN)    :: gather           !< communication pattern
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::gather_field2D')
#if !defined(NOMPI)
    INTEGER  :: ierrstat, i, jc,jb
    INTEGER, ALLOCATABLE :: recv_buf(:)

    IF (dbg_level >= 11) WRITE (0,*) "# gather data on PE ", gather%rank0
    ! allocate temporary data structures:
    ALLOCATE(recv_buf(gather%total_dim), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! communicate field data:
    CALL p_gatherv(field_in, gather%local_dim,                 &
      &            recv_buf, gather%local_size, gather%displs, &
      &            gather%rank0, p_comm_work)

    ! map received values to the right (global) positions:
!$OMP PARALLEL DO PRIVATE(jc,jb)
    DO i=1,gather%total_dim
      jc = idx_no(gather%global_idx(i))
      jb = blk_no(gather%global_idx(i))
      field_out(jc,jb) = recv_buf(i)
    END DO
!$OMP END PARALLEL DO

    ! clean up
    DEALLOCATE(recv_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#else
    IF (SIZE(field_out,1) /= SIZE(field_in,1)) CALL finish(routine, "Mismatch dimension 1")
    IF (SIZE(field_out,2) /= SIZE(field_in,2)) CALL finish(routine, "Mismatch dimension 2")
    field_out(:,:) = field_in(:,:)
#endif
  END SUBROUTINE gather_field2D_int


  !> Gather distributed 3D field data from working PEs in a global
  !  field on rank0.
  !
  SUBROUTINE gather_field3D(gather, nlev, field_in, field_out)
    TYPE(t_gather),    INTENT(IN)    :: gather             !< communication pattern
    INTEGER,           INTENT(IN)    :: nlev               !< no. of levels (dim 2)
    REAL(wp),          INTENT(IN)    :: field_in(:,:,:)    !< input field (size may differ on PEs)
    REAL(wp),          INTENT(INOUT) :: field_out(:,:,:)   !< global output field
#if !defined(NOMPI)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::gather_field3D')
    INTEGER  :: ierrstat, jc,jk,jb, iglobal_idx, i
    REAL(wp), ALLOCATABLE :: recv_buf(:)
    INTEGER :: ilocal_size(p_n_work), idispls(p_n_work), shape_snd(3)
    REAL(wp) :: send_buf(UBOUND(field_in,2),UBOUND(field_in,1),UBOUND(field_in,3))

    IF (dbg_level >= 11) WRITE (0,*) "# gather data on PE ", gather%rank0
    ! allocate temporary data structures:
    ALLOCATE(recv_buf(gather%total_dim*nlev), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! communicate field data:
    ilocal_size(:) = gather%local_size(:) * nlev
    idispls(:)     = gather%displs(:) * nlev
    shape_snd = (/ UBOUND(field_in,2),UBOUND(field_in,1),UBOUND(field_in,3) /)
    send_buf(:,:,:) = RESHAPE(field_in, shape_snd, order=(/2,1,3/) )

    CALL p_gatherv(send_buf, gather%local_dim*nlev, &
      &            recv_buf, ilocal_size, idispls,  &
      &            gather%rank0, p_comm_work)

    ! map received values to the right (global) positions:
!$OMP PARALLEL DO PRIVATE(jc,jk,jb,iglobal_idx)
    DO i=1,gather%total_dim
      jc = idx_no(gather%global_idx(i))
      jb = blk_no(gather%global_idx(i))
      DO jk=1,nlev
        iglobal_idx = jk + (i-1)*nlev
        field_out(jc,jk,jb) = recv_buf(iglobal_idx)
      END DO
    END DO
!$OMP END PARALLEL DO

    ! clean up
    DEALLOCATE(recv_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#else
    field_out(:,:,:) = field_in(:,:,:)
#endif
  END SUBROUTINE gather_field3D

END MODULE mo_remap_sync
