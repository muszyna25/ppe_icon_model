! -----------------------------------------------------------------------
! MPI global reshuffle, 05/2016 F. Prill
!
! This program works on a distributed-memory array, where each PE's
! part is defined by global indices "owner_idx" (local size!).  The
! "reshuffle" operation means that each PE writes a number of "nsend"
! values (where "nsend" may be different for each PE) to the array at
! global indices "glb_idx". This involves, of course, some
! send/receive operations, since the destination indices may be stored
! on a different PE. This implementation of the "reshuffle" operation
! involves no global-size arrays, but it costs 2x MPI_ALLTOALL
! operations and 4x MPI_ALLTOALLV operations.
!
! IN:  glb_idx(1,...,nsend)    : global indices
!      values(1,...,nsend)     : values to send
!      owner_idx(1,...,nlocal) : global indices owned by this PE
! OUT: recv_buf(1,...,nlocal)  : buffer for received values
!
! TODO[FP]
! - Non-MPI version missing for this subroutine.

! -----------------------------------------------------------------------

MODULE mo_reshuffle
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE t_regular_partition
    INTEGER :: n_indices, n_pes, irank
  CONTAINS
    PROCEDURE :: local_size => regular_partition_local_size
    PROCEDURE :: glb2local  => regular_partition_glb2local
    PROCEDURE :: glbidx2pe  => regular_partition_glbidx2pe
  END TYPE t_regular_partition

CONTAINS

  INTEGER FUNCTION regular_partition_local_size(this)
    CLASS(t_regular_partition), INTENT(IN) :: this
    regular_partition_local_size = ((this%n_indices+this%n_pes-1) / this%n_pes)
  END FUNCTION regular_partition_local_size

  INTEGER FUNCTION regular_partition_glb2local(this,iidx)
    CLASS(t_regular_partition), INTENT(IN) :: this
    INTEGER,                    INTENT(IN) :: iidx
    regular_partition_glb2local = iidx - this%local_size()*this%irank 
  END FUNCTION regular_partition_glb2local

  INTEGER FUNCTION regular_partition_glbidx2pe(this,iidx)
    CLASS(t_regular_partition), INTENT(IN) :: this
    INTEGER,                    INTENT(IN) :: iidx

    regular_partition_glbidx2pe = -1
    IF ((iidx > 0) .AND. (iidx <= this%n_indices) .AND. (this%n_indices >= this%n_pes)) THEN
      regular_partition_glbidx2pe = (iidx-1)/this%local_size()
    END IF
  END FUNCTION regular_partition_glbidx2pe

  SUBROUTINE swap_int(a, i,j, permutation)
    INTEGER,  INTENT(INOUT)           :: a(:)           !< array for in-situ sorting
    INTEGER,  INTENT(IN)              :: i,j            !< indices to be exchanged
    INTEGER,  INTENT(INOUT), OPTIONAL :: permutation(:) !< (optional) permutation of indices
    ! local variables
    INTEGER :: t, t_p

    t    = a(i)
    a(i) = a(j)
    a(j) = t
    IF (PRESENT(permutation)) THEN
      t_p            = permutation(i)
      permutation(i) = permutation(j)
      permutation(j) = t_p
    END IF
  END SUBROUTINE swap_int

  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of INTEGER values.
  ! 
  !  Ordering after the sorting process: smallest...largest.
  !
  RECURSIVE SUBROUTINE quicksort_int(a, permutation, l_in, r_in)
    INTEGER,  INTENT(INOUT)           :: a(:)           !< array for in-situ sorting
    INTEGER,  INTENT(INOUT), OPTIONAL :: permutation(:) !< (optional) permutation of indices
    INTEGER,  INTENT(IN),    OPTIONAL :: l_in,r_in      !< left, right partition indices
    ! local variables
    INTEGER :: i,j,l,r,t_p,t,v,m

    IF (PRESENT(l_in)) THEN
      l = l_in
    ELSE
      l = 1
    END IF
    IF (PRESENT(r_in)) THEN
      r = r_in
    ELSE
      r = SIZE(a,1)
    END IF
    IF (r>l) THEN
      i = l-1
      j = r
      
      ! median-of-three selection of partitioning element
      IF ((r-l) > 3) THEN 
        m = (l+r)/2
        IF (a(l)>a(m))  CALL swap_int(a, l,m, permutation)
        IF (a(l)>a(r)) THEN
          CALL swap_int(a, l,r, permutation)
        ELSE IF (a(r)>a(m)) THEN
          CALL swap_int(a, r,m, permutation)
        END IF
      END IF

      v = a(r)
      LOOP : DO
        CNTLOOP1 : DO
          i = i+1
          IF (a(i) >= v) EXIT CNTLOOP1
        END DO CNTLOOP1
        CNTLOOP2 : DO
          j = j-1
          IF ((a(j) <= v) .OR. (j==1)) EXIT CNTLOOP2
        END DO CNTLOOP2
        t    = a(i)
        a(i) = a(j)
        a(j) = t
        IF (PRESENT(permutation)) THEN
          t_p            = permutation(i)
          permutation(i) = permutation(j)
          permutation(j) = t_p
        END IF
        IF (j <= i) EXIT LOOP
      END DO LOOP
      a(j) = a(i)
      a(i) = a(r)
      a(r) = t
      IF (PRESENT(permutation)) THEN
        permutation(j) = permutation(i)
        permutation(i) = permutation(r)
        permutation(r) = t_p
      END IF
      CALL quicksort_int(a,permutation,l,i-1)
      CALL quicksort_int(a,permutation,i+1,r)
    END IF
  END SUBROUTINE quicksort_int

  SUBROUTINE calc_displs(i_pe, displs)
    INTEGER, INTENT(IN)    :: i_pe(:)    !< i_pe(i) contains the PE where to send/receive index "i" to/from
    INTEGER, INTENT(INOUT) :: displs(0:) !< displs(p) contains the start offset for PE "p"
    INTEGER :: i, counts(0:(SIZE(displs)-1))

    counts(:) = 0
    DO i=1,SIZE(i_pe)
      counts(i_pe(i)) = counts(i_pe(i)) + 1
    END DO
    displs(:) = 0
    DO i=1,SIZE(displs)
      displs(i) = displs(i-1) + counts(i-1)
    END DO
  END SUBROUTINE calc_displs

  ! -----------------------------------------------------------------------
  ! This subroutine works on a distributed-memory array, where each
  ! PE's part is defined by global indices "owner_idx" (local size!).
  ! The "reshuffle" operation means that each PE writes a number of
  ! values (which may be different for each PE) to the array at global
  ! indices "glb_idx". This involves, of course, some send/receive
  ! operations, since the destination indices may be stored on a
  ! different PE. This implementation of the "reshuffle" operation
  ! involves no global-size arrays.
  !
  ! MPI_ALLTOALL operations needed:
  ! 1x MPI_ALLTOALL : send/receive counts
  ! 1x MPI_ALLTOALLV: send/receive indices
  ! 2x MPI_ALLTOALLV: send/receive values
  !
  ! IN:  glb_idx(1,...,nsend)    : global indices
  !      values(1,...,nsend)     : values to send
  !      owner_idx(1,...,nlocal) : global indices owned by this PE
  ! OUT: recv_buf(1,...,nlocal)  : buffer for received values
  ! -----------------------------------------------------------------------
  SUBROUTINE reshuffle(in_glb_idx, in_values, owner_idx, nglb_indices, communicator, out_values, ierr)
    INTEGER, INTENT(IN)    :: in_glb_idx(:)    !< global indices to which "values" correspond
    INTEGER, INTENT(IN)    :: in_values(:)     !< values to send
    INTEGER, INTENT(IN)    :: owner_idx(:)     !< array indices "owned" by local PE
    INTEGER, INTENT(IN)    :: nglb_indices     !< total size of distributed array
    INTEGER, INTENT(IN)    :: communicator     !< MPI comm.
    INTEGER, INTENT(INOUT) :: out_values(:)    !< resulting local part of distributed array
    INTEGER, INTENT(OUT)   :: ierr             !< error return code 
    ! local variables
    INTEGER                   :: i, j, nsend, nlocal, rank, isize, ierr, communicator,     &
      &                          npairs_recv, local_idx, npairs_recv_owner                 
    INTEGER, ALLOCATABLE      :: icounts(:), irecv(:), irecv_idx(:), recv_vals(:),         &
      &                          isend_idx(:), send_vals(:), i_pe(:), glb_idx(:),          &
      &                          values(:), permutation(:), send_displs(:),                &
      &                          recv_displs(:), reg_partition_buf(:),                     &
      &                          reordered_owner_idx(:), isendbuf(:), i_pe_owner(:),       &
      &                          icounts_owner(:), irecv_owner(:), icounts_buf(:,:),       &
      &                          irecv_buf(:,:), permutation_owner(:), irecv_idx_owner(:), &
      &                          send_displs_owner(:), recv_displs_owner(:),               &
      &                          send_counts2(:), recv_counts2(:), send_displs2(:),        &
      &                          recv_displs2(:), tmp_idx(:), irecv_tmp(:)
    TYPE(t_regular_partition) :: reg_partition
        
    ierr   = 0
    nlocal = SIZE(owner_idx)
    nsend  = SIZE(in_glb_idx)

    ! consistency checks:
    IF (SIZE(out_values) /= nlocal) THEN
      WRITE (0,*) "SIZE(out_values) /= nlocal"
      ierr = 1
    END IF
    IF (SIZE(in_values) /= nsend) THEN
      WRITE (0,*) "SIZE(in_values) /= nsend"
      ierr = 1
    END IF
    IF ((MAXVAL(in_glb_idx) > nglb_indices) .OR. (MINVAL(in_glb_idx) < 1)) THEN
      WRITE (0,*) "MAXVAL(in_glb_idx) > nglb_indices"
      ierr = 1
    END IF
    IF ((MAXVAL(owner_idx) > nglb_indices)  .OR. (MINVAL(owner_idx) < 1)) THEN
      WRITE (0,*) "MAXVAL(in_glb_idx) > nglb_indices"
      ierr = 1
    END IF

    CALL MPI_COMM_RANK(communicator, rank, ierr)      ; IF (ierr /= 0) WRITE (*,*) "MPI Error!"
    CALL MPI_COMM_SIZE(communicator, isize, ierr)     ; IF (ierr /= 0) WRITE (*,*) "MPI Error!"

    ! 1)   intermediate, regular partitioning of global index space:
    reg_partition%n_indices = nglb_indices
    reg_partition%n_pes     = isize
    reg_partition%irank     = rank

    ! 1.0) for every "glb_idx" entry compute the destination PE,
    !      for every "owner_idx" entry compute the source PE
    ALLOCATE(i_pe(nsend), i_pe_owner(nlocal))
    DO i=1,nsend
      i_pe(i) = reg_partition%glbidx2pe(in_glb_idx(i))
    END DO
    DO i=1,nlocal
      i_pe_owner(i) = reg_partition%glbidx2pe(owner_idx(i))
    END DO

    ! 1.2) each PE "a" sends to every other PE "b" the no. of
    !      "glb_idx" entries for "b",
    !      each PE "a" sends to every other PE "b" the no. of
    !      "owner_idx" entries he wants to receive
    ALLOCATE(icounts(0:(isize-1)), irecv(0:(isize-1)), &
      &      icounts_owner(0:(isize-1)), irecv_owner(0:(isize-1)))
    icounts(:)        = 0
    irecv(:)          = 0
    DO i=1,nsend
      icounts(i_pe(i)) = icounts(i_pe(i)) + 1
    END DO
    icounts_owner(:)  = 0
    irecv_owner(:)    = 0
    DO i=1,nlocal
      icounts_owner(i_pe_owner(i)) = icounts_owner(i_pe_owner(i)) + 1
    END DO

    ! we use a temp buffer to exchange both index arrays in one step
    ALLOCATE(icounts_buf(2,0:(isize-1)), irecv_buf(2,0:(isize-1)))
    icounts_buf(1,:) = icounts(:)
    icounts_buf(2,:) = icounts_owner(:)
    irecv_buf(1,:)   = irecv(:)
    irecv_buf(2,:)   = irecv_owner(:)
    CALL MPI_ALLTOALL(icounts_buf, 2, MPI_INTEGER, irecv_buf, 2, MPI_INTEGER, communicator, ierr)
    IF (ierr /= 0) WRITE (*,*) "MPI Error!"
    icounts(:)       = icounts_buf(1,:)
    icounts_owner(:) = icounts_buf(2,:)
    irecv(:)         = irecv_buf(1,:)
    irecv_owner(:)   = irecv_buf(2,:)

    ! 1.1) reorder "glb_idx" and "values" according to the rank of the
    !      destination PE ;
    !      reorder "owner_idx" according to the rank of the source PE
    ALLOCATE(permutation(nsend), glb_idx(nsend), values(nsend), &
      &      permutation_owner(nlocal), reordered_owner_idx(nlocal))
    permutation(:)       = (/ ( i, i=1,nsend) /)
    permutation_owner(:) = (/ ( i, i=1,nlocal) /)
    CALL quicksort_int(i_pe, permutation)
    CALL quicksort_int(i_pe_owner, permutation_owner)
    glb_idx(:)             = in_glb_idx(permutation(:))
    values(:)              = in_values(permutation(:))
    reordered_owner_idx(:) = owner_idx(permutation_owner(:))
 
    ! 1.3a) PE "a" sends "glb_idx", "owner_idx" to "b"
    npairs_recv       = SUM(irecv)
    npairs_recv_owner = SUM(irecv_owner)

    ALLOCATE(irecv_idx(npairs_recv),         recv_vals(npairs_recv),         &
      &      send_displs(0:(isize-1)),       recv_displs(0:(isize-1)),       &
      &      send_displs_owner(0:(isize-1)), recv_displs_owner(0:(isize-1)), &
      &      send_displs2(0:(isize-1)),      recv_displs2(0:(isize-1)),      &
      &      irecv_idx_owner(npairs_recv_owner))
    CALL calc_displs(i_pe, send_displs)
    recv_displs(0) = 0
    DO i=1,SIZE(recv_displs)
      recv_displs(i) = recv_displs(i-1) + irecv(i-1)
    END DO
    CALL calc_displs(i_pe_owner, send_displs_owner)
    recv_displs_owner(0) = 0
    DO i=1,SIZE(recv_displs_owner)
      recv_displs_owner(i) = recv_displs_owner(i-1) + irecv_owner(i-1)
    END DO

    ! we use a temp buffer to send of "glb_idx" and "owner_idx" in one
    ! step:
    ALLOCATE(send_counts2(0:(isize-1)), recv_counts2(0:(isize-1)), &
      &      tmp_idx(nsend+nlocal), irecv_tmp(npairs_recv+npairs_recv_owner))
    send_counts2(:) = icounts(:)     + icounts_owner(:)
    recv_counts2(:) = irecv(:)       + irecv_owner(:)
    send_displs2(:) = send_displs(:) + send_displs_owner(:)
    recv_displs2(:) = recv_displs(:) + recv_displs_owner(:)
    j = 1
    DO i=0,(isize-1)
      tmp_idx(j:(j+icounts(i)-1))       = glb_idx( (1+send_displs(i)):(send_displs(i)+icounts(i)) )
      j = j+icounts(i)
      tmp_idx(j:(j+icounts_owner(i)-1)) = reordered_owner_idx( (1+send_displs_owner(i)):(send_displs_owner(i)+icounts_owner(i)) )
      j = j+icounts_owner(i)
    END DO
    CALL MPI_ALLTOALLV(tmp_idx, send_counts2, send_displs2, MPI_INTEGER, irecv_tmp, &
      &                recv_counts2, recv_displs2, MPI_INTEGER, communicator, ierr)
    IF (ierr /= 0) WRITE (*,*) "MPI Error!"
    j = 1
    DO i=0,(isize-1)
      irecv_idx( (1+recv_displs(i)):(recv_displs(i)+irecv(i)) ) = irecv_tmp(j:(j+irecv(i)-1))
      j = j+irecv(i)
      irecv_idx_owner( (1+recv_displs_owner(i)):(recv_displs_owner(i)+irecv_owner(i)) ) = irecv_tmp(j:(j+irecv_owner(i)-1)) 
      j = j+irecv_owner(i)
    END DO

    ! 1.3b) PE "a" sends values to "b"
    CALL MPI_ALLTOALLV(values, icounts, send_displs, MPI_INTEGER, recv_vals, &
      &                irecv, recv_displs, MPI_INTEGER, communicator, ierr)
    IF (ierr /= 0) WRITE (*,*) "MPI Error!"

    ! 1.4) each PE inserts the received index/value pairs into "its"
    !      part of the global index space (wrt. to the regular
    !      partition).
    ALLOCATE(reg_partition_buf(reg_partition%local_size()))
    reg_partition_buf(:) = 0
    DO i=1,npairs_recv
      local_idx = reg_partition%glb2local(irecv_idx(i))
      reg_partition_buf(local_idx) = recv_vals(i)
    END DO

    ! 2.4) collect send buffer based on received "owner_idx"
    ALLOCATE(isendbuf(npairs_recv_owner))
    DO i=1,npairs_recv_owner
      isendbuf(i) = reg_partition_buf( reg_partition%glb2local(irecv_idx_owner(i)) )
    END DO

    ! 2.5) PE "b" sends values back to "a"
    DEALLOCATE(irecv_idx_owner)
    ALLOCATE(irecv_idx_owner(SUM(icounts_owner)))
    CALL MPI_ALLTOALLV(isendbuf, irecv_owner, recv_displs_owner, MPI_INTEGER, irecv_idx_owner, &
      &                icounts_owner, send_displs_owner, MPI_INTEGER, communicator, ierr)
    IF (ierr /= 0) WRITE (*,*) "MPI Error!"
    
    ! 2.5) each PE inserts the received index/values pairs into its
    !      part of the global index space
    out_values(permutation_owner(1:nlocal)) = irecv_idx_owner(1:nlocal)
    
    ! 3)   clean-up
    DEALLOCATE(icounts, irecv, i_pe, permutation, glb_idx, values, send_displs, recv_displs, &
      &        irecv_idx, recv_vals, reg_partition_buf, reordered_owner_idx, isendbuf,       &
      &        i_pe_owner, icounts_owner, irecv_owner, icounts_buf, irecv_buf,               &
      &        permutation_owner, irecv_idx_owner, send_displs_owner, recv_displs_owner,     &
      &        send_counts2, recv_counts2, tmp_idx, irecv_tmp)
  END SUBROUTINE reshuffle

END MODULE mo_reshuffle
