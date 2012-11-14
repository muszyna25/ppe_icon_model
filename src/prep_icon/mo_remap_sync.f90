!! Subroutines for MPI communication.
!!
MODULE mo_remap_sync

  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_communication,      ONLY: idx_no, blk_no
  USE mo_remap_config,       ONLY: dbg_level
  USE mo_mpi,                ONLY: p_real_dp, p_comm_work, p_n_work, p_int, &
    &                              get_my_mpi_work_id, p_gather, p_gatherv, &
    &                              p_scatterv, p_wait, p_isend, p_recv,     &
    &                              p_clear_request
  USE mo_remap_shared,       ONLY: t_grid
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: allocate_gather_c, finalize_gather_c
  PUBLIC :: scatter_field2D_c, gather_field2D_c
  PUBLIC :: gather_field3D_c
  PUBLIC :: igather_field2D_c
  PUBLIC :: t_gather_c

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_sync')  

  !> Data structure for local-to-global communication of distributed
  !> field data
  !
  TYPE t_gather_c
    INTEGER,  ALLOCATABLE :: global_idx(:), &   !< mapping: (global, PE-ordered array) -> (global array)
      &                      local_size(:), &   !< local field size
      &                      displs(:)          !< offset in global, PE-ordered array
    INTEGER :: local_dim, total_dim, rank0, nblks_c, npromz_c
    
    INTEGER :: mpi_request

    !> mpi_isend buffer
    REAL(wp), ALLOCATABLE :: isend_buf(:)
  END TYPE t_gather_c

CONTAINS

  !> Prepare communication data structures.
  !
  SUBROUTINE allocate_gather_c(grid, rank0, gather_c)
    TYPE (t_grid), INTENT(IN)       :: grid             !< local grid partition
    INTEGER,       INTENT(IN)       :: rank0            !< root PE, where data is collected.
    TYPE(t_gather_c), INTENT(INOUT) :: gather_c         !< communication pattern
#if !defined(NOMPI)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::allocate_gather_c')
    INTEGER :: ierrstat, i, this_pe
    
    IF (dbg_level >= 5) &
      &   WRITE (0,*) "# allocate communication pattern"
    this_pe = get_my_mpi_work_id()
    gather_c%rank0     = rank0
    gather_c%local_dim = grid%p_patch%n_patch_cells
    ALLOCATE(gather_c%displs(p_n_work), gather_c%local_size(p_n_work), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! communicate local field sizes:
    gather_c%local_size(:) = 0
    CALL p_gather(gather_c%local_dim, gather_c%local_size, gather_c%rank0, p_comm_work)

    IF (this_pe == rank0) THEN
      gather_c%total_dim = SUM(gather_c%local_size)
      gather_c%nblks_c   = blk_no(grid%p_patch%n_patch_cells_g)
      gather_c%npromz_c  = grid%p_patch%n_patch_cells_g - nproma*(gather_c%nblks_c-1)
    ELSE
      gather_c%total_dim = 0
      gather_c%nblks_c   = 0
      gather_c%npromz_c  = 0
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


    ! prepare mpi_isend:
    CALL p_clear_request(gather_c%mpi_request)
    ALLOCATE(gather_c%isend_buf(gather_c%local_dim), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE for isend failed!")
#endif
  END SUBROUTINE allocate_gather_c


  !> Free communication data structures.
  !
  RECURSIVE SUBROUTINE finalize_gather_c(gather_c, opt_gather_c2, opt_gather_c3, opt_gather_c4)
    TYPE(t_gather_c), INTENT(INOUT) :: gather_c
    TYPE(t_gather_c), INTENT(INOUT), OPTIONAL :: opt_gather_c2, opt_gather_c3, opt_gather_c4
#if !defined(NOMPI)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::finalize_gather_c')
    INTEGER :: ierrstat

    ! the following recursion makes this more readable on the caller
    ! side:
    IF (PRESENT(opt_gather_c4)) THEN
      CALL finalize_gather_c(opt_gather_c4)
      CALL finalize_gather_c(gather_c, opt_gather_c2, opt_gather_c3)
    ELSE IF (PRESENT(opt_gather_c3)) THEN
      CALL finalize_gather_c(opt_gather_c3, opt_gather_c4)
      CALL finalize_gather_c(gather_c, opt_gather_c2)
    ELSE IF (PRESENT(opt_gather_c2)) THEN
      CALL finalize_gather_c(opt_gather_c2)
      CALL finalize_gather_c(gather_c)
    ELSE

      ! clean up
      IF (dbg_level >= 10) WRITE (0,*) "# finalize gather/scatter pattern."
      DEALLOCATE(gather_c%displs, gather_c%global_idx, gather_c%local_size, &
        &        gather_c%isend_buf, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
#endif
  END SUBROUTINE finalize_gather_c


  !> Scatter a global 2D field data to working PEs.
  !
  SUBROUTINE scatter_field2D_c(gather_c, field_in, field_out)
    REAL(wp),         INTENT(IN)    :: field_in(:,:)    !< input field (only required for PE rank0)
    REAL(wp),         INTENT(INOUT) :: field_out(:,:)   !< local output field (size may differ on PEs)
    TYPE(t_gather_c), INTENT(IN)    :: gather_c         !< communication pattern
#if !defined(NOMPI)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::gather_field2D_c')
    INTEGER  :: ierrstat, i, jc,jb
    REAL(wp), ALLOCATABLE :: send_buf(:)

    IF (dbg_level >= 11) &
      & WRITE (0,*) "# scatter data from PE ", gather_c%rank0
    ! allocate temporary data structures:
    ALLOCATE(send_buf(gather_c%total_dim), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! map values to the right positions in the send buffer:
    IF (get_my_mpi_work_id() == gather_c%rank0) THEN
!$OMP PARALLEL DO PRIVATE(jc,jb)
      DO i=1,gather_c%total_dim
        jc = idx_no(gather_c%global_idx(i))
        jb = blk_no(gather_c%global_idx(i))
        send_buf(i) = field_in(jc,jb)
      END DO
!$OMP END PARALLEL DO
    END IF ! if rank0

    ! communicate field data:
    CALL p_scatterv(send_buf, gather_c%local_size, gather_c%displs,   &
      &             field_out, gather_c%local_dim,                    &
      &             gather_c%rank0, p_comm_work)

    ! clean up
    DEALLOCATE(send_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#else
    field_out(:,:) = field_in(:,:)
#endif
  END SUBROUTINE scatter_field2D_c


  !> Gather distributed 2D field data from working PEs in a global
  !  field on rank0.
  !
  SUBROUTINE gather_field2D_c(gather_c, field_in, field_out)
    REAL(wp),          INTENT(IN)    :: field_in(:,:)    !< input field (size may differ on PEs)
    REAL(wp),          INTENT(INOUT) :: field_out(:,:)   !< global output field
    TYPE(t_gather_c),  INTENT(IN)    :: gather_c         !< communication pattern
#if !defined(NOMPI)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::gather_field2D_c')
    INTEGER  :: ierrstat, i, jc,jb
    REAL(wp), ALLOCATABLE :: recv_buf(:)

    IF (dbg_level >= 11) WRITE (0,*) "# gather data on PE ", gather_c%rank0
    ! allocate temporary data structures:
    ALLOCATE(recv_buf(gather_c%total_dim), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! communicate field data:
    CALL p_gatherv(field_in, gather_c%local_dim,                   &
      &            recv_buf, gather_c%local_size, gather_c%displs, &
      &            gather_c%rank0, p_comm_work)

    ! map received values to the right (global) positions:
!$OMP PARALLEL DO PRIVATE(jc,jb)
    DO i=1,gather_c%total_dim
      jc = idx_no(gather_c%global_idx(i))
      jb = blk_no(gather_c%global_idx(i))
      field_out(jc,jb) = recv_buf(i)
    END DO
!$OMP END PARALLEL DO

    ! clean up
    DEALLOCATE(recv_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#else
    field_out(:,:) = field_in(:,:)
#endif
  END SUBROUTINE gather_field2D_c


  !> Gather distributed 3D field data from working PEs in a global
  !  field on rank0.
  !
  SUBROUTINE gather_field3D_c(gather_c, nlev, field_in, field_out)
    TYPE(t_gather_c),  INTENT(IN)    :: gather_c           !< communication pattern
    INTEGER,           INTENT(IN)    :: nlev               !< no. of levels (dim 2)
    REAL(wp),          INTENT(IN)    :: field_in(:,:,:)    !< input field (size may differ on PEs)
    REAL(wp),          INTENT(INOUT) :: field_out(:,:,:)   !< global output field
#if !defined(NOMPI)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::gather_field2D_c')
    INTEGER  :: ierrstat, jc,jk,jb, iendidx, iglobal_idx
    REAL(wp), ALLOCATABLE :: recv_buf(:)
    INTEGER :: ilocal_size(p_n_work), idispls(p_n_work)

    IF (dbg_level >= 11) WRITE (0,*) "# gather data on PE ", gather_c%rank0
    ! allocate temporary data structures:
    ALLOCATE(recv_buf(gather_c%total_dim*nlev), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! communicate field data:
    ilocal_size(:) = gather_c%local_size(:) * nlev
    idispls(:)     = gather_c%displs(:) * nlev
    CALL p_gatherv(field_in, gather_c%local_dim*nlev, &
      &            recv_buf, ilocal_size, idispls,    &
      &            gather_c%rank0, p_comm_work)

    ! map received values to the right (global) positions:
!$OMP PARALLEL DO PRIVATE(jc,jk,jb,iglobal_idx)
    DO jb=1,gather_c%nblks_c
      DO jk=1,nlev
        iendidx = nproma
        IF (jb == gather_c%nblks_c) iendidx = gather_c%npromz_c
        DO jc=1,iendidx
          iglobal_idx = (jb-1)*nproma + (jk-1)*nproma*nlev + jc
          field_out(jc,jk,jb) = recv_buf(iglobal_idx)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! clean up
    DEALLOCATE(recv_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#else
    field_out(:,:,:) = field_in(:,:,:)
#endif
  END SUBROUTINE gather_field3D_c


  !> Gather distributed 2D field data from working PEs in a global
  !  field on rank0 with a number of ISEND/RECV calls.
  !
  SUBROUTINE igather_field2D_c(gather_c, field_in, field_out)
    REAL(wp),          INTENT(IN)    :: field_in(:,:)    !< input field (size may differ on PEs)
    REAL(wp),          INTENT(INOUT) :: field_out(:,:)   !< global output field
    TYPE(t_gather_c),  INTENT(INOUT) :: gather_c         !< communication pattern
#if !defined(NOMPI)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::igather_field2D_c')
    INTEGER  :: ierrstat, i, jc,jb
    REAL(wp), ALLOCATABLE :: recv_buf(:)

    IF (dbg_level >= 11) WRITE (0,*) "# gather data on PE ", gather_c%rank0
    ! allocate temporary data structures:
    ALLOCATE(recv_buf(gather_c%total_dim), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! communicate field data:
    gather_c%isend_buf(:) = RESHAPE(field_in, (/ gather_c%local_dim /))
    CALL p_wait(gather_c%mpi_request)
    CALL p_isend(gather_c%isend_buf, gather_c%rank0, 0, &
      &          gather_c%local_dim, p_comm_work, gather_c%mpi_request)

    IF (get_my_mpi_work_id() == gather_c%rank0) THEN
      DO i=1,SIZE(gather_c%local_size)
        CALL p_recv(recv_buf, (i-1), 0, gather_c%local_size(i), &
          &         p_comm_work, gather_c%displs(i)+1)
      END DO
    END IF


    ! map received values to the right (global) positions:
!$OMP PARALLEL DO PRIVATE(jc,jb)
    DO i=1,gather_c%total_dim
      jc = idx_no(gather_c%global_idx(i))
      jb = blk_no(gather_c%global_idx(i))
      field_out(jc,jb) = recv_buf(i)
    END DO
!$OMP END PARALLEL DO

    ! clean up
    DEALLOCATE(recv_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#else
    field_out(:,:) = field_in(:,:)
#endif
  END SUBROUTINE igather_field2D_c

END MODULE mo_remap_sync
