!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_util_uuid

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_SIGNED_CHAR, C_NULL_CHAR, &
    &                                    C_DOUBLE, C_INT, C_PTR, C_F_POINTER, C_LOC
  USE MPI

  IMPLICIT NONE 

  PRIVATE

  PUBLIC :: t_uuid

  PUBLIC :: uuid_generate
  PUBLIC :: uuid_generate_parallel

  PUBLIC :: compare_uuid

  PUBLIC :: uuid_parse
  PUBLIC :: uuid_unparse

  PUBLIC :: uuid2char
  PUBLIC :: char2uuid

  PUBLIC :: clear_uuid

  PUBLIC :: OPERATOR(==)

  PUBLIC :: UUID_STRING_LENGTH
  PUBLIC :: UUID_DATA_LENGTH
  PUBLIC :: UUID_EQUAL, UUID_EQUAL_LIMITED_ACCURACY, UUID_UNEQUAL

  INTEGER, PARAMETER :: UUID_STRING_LENGTH = 36
  INTEGER, PARAMETER :: UUID_DATA_LENGTH   = 16

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_uuid'

  TYPE, BIND(C) :: t_uuid
    INTEGER(C_SIGNED_CHAR) :: data(16)
  END TYPE t_uuid

  TYPE, BIND(C) :: t_context
    INTEGER(C_INT) :: f64(2)
    INTEGER(C_INT) :: f48(2)
    INTEGER(C_INT) :: tl64(2)
    INTEGER(C_INT) :: tl48(2)
    INTEGER(C_INT) :: max_zero
  END TYPE t_context

  ENUM, BIND(c)
    ENUMERATOR :: UUID_EQUAL                  = 0
    ENUMERATOR :: UUID_EQUAL_LIMITED_ACCURACY = 1
    ENUMERATOR :: UUID_UNEQUAL                = 2
  END ENUM

  INTERFACE OPERATOR (==)
    MODULE PROCEDURE uuid_compare
  END INTERFACE OPERATOR (==)

  INTERFACE
    SUBROUTINE my_uuid_generate(val, nval, uuid) BIND(C,NAME='uuid_generate')
      IMPORT :: C_DOUBLE, C_INT, t_uuid
      REAL(c_double),   INTENT(IN)         :: val(*)
      INTEGER(c_int),   INTENT(IN), VALUE  :: nval
      TYPE(t_uuid),     INTENT(OUT)        :: uuid
    END SUBROUTINE my_uuid_generate
  END INTERFACE

  INTERFACE
    INTEGER(C_INT) FUNCTION my_compare_uuid(uuid_A, uuid_B, min_difference) BIND(C,NAME='compare_UUID')
      IMPORT :: t_uuid, C_INT, C_DOUBLE
      TYPE(t_uuid),     INTENT(IN), VALUE  :: uuid_A, uuid_B
      REAL(c_double),   INTENT(OUT)        :: min_difference
    END FUNCTION my_compare_uuid
  END INTERFACE

  INTERFACE
    SUBROUTINE my_uuid_unparse(uuid_string, uuid) BIND(C,NAME='uuid_unparse')
      IMPORT :: C_CHAR, t_uuid
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: uuid_string
      TYPE(t_uuid),                    INTENT(in)  :: uuid
    END SUBROUTINE my_uuid_unparse
  END INTERFACE

  INTERFACE
    SUBROUTINE my_uuid_parse(uuid, uuid_string) BIND(C,NAME='uuid_parse')
      IMPORT :: C_CHAR, t_uuid
      TYPE(t_uuid),                    INTENT(out) :: uuid
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in)  :: uuid_string
    END SUBROUTINE my_uuid_parse
  END INTERFACE

  INTERFACE
    SUBROUTINE my_encode_uuid(fingerprint, uuid) BIND(C,NAME='encode_uuid')
      IMPORT :: t_uuid, C_PTR
      TYPE(C_PTR), VALUE        :: fingerprint
      TYPE(t_uuid), INTENT(out) :: uuid
    END SUBROUTINE my_encode_uuid
  END INTERFACE

  INTERFACE
    FUNCTION my_uuid_scan_data(val, nval) BIND(C,NAME='uuid_scan_data') RESULT(fingerprint)
      IMPORT :: C_PTR, C_DOUBLE, C_INT
      REAL(C_DOUBLE),     INTENT(IN) :: val(*)
      INTEGER(C_INT), INTENT(IN), VALUE :: nval
      TYPE(C_PTR) :: fingerprint
    END FUNCTION my_uuid_scan_data
  END INTERFACE

  INTERFACE
    FUNCTION my_concat_fingerprints(context0, context1) BIND(C,NAME='concat_fingerprints') RESULT(fingerprint)
      IMPORT :: C_PTR
      TYPE(C_PTR), VALUE :: context0
      TYPE(C_PTR), VALUE :: context1
      TYPE(C_PTR) :: fingerprint
    END FUNCTION my_concat_fingerprints
  END INTERFACE

  INTERFACE
    SUBROUTINE deallocate_c(ptr) BIND(C,NAME='deallocate_fingerprint') 
      IMPORT :: C_PTR
      TYPE(C_PTR), VALUE :: ptr
    END SUBROUTINE deallocate_c
  END INTERFACE

CONTAINS
 
  SUBROUTINE uuid_parse(uuid_string, uuid)
    CHARACTER(len=*), INTENT(in)  :: uuid_string
    TYPE(t_uuid),     INTENT(out) :: uuid
    CALL my_uuid_parse(uuid, TRIM(uuid_string)//C_NULL_CHAR)
  END SUBROUTINE uuid_parse

  SUBROUTINE uuid_unparse(uuid, uuid_string)
    TYPE(t_uuid),     INTENT(in)  :: uuid
    CHARACTER(len=*), INTENT(out) :: uuid_string
    CALL my_uuid_unparse(uuid_string, uuid)
  END SUBROUTINE uuid_unparse

  FUNCTION uuid_compare(uuid1, uuid2)
    LOGICAL :: uuid_compare
    TYPE(t_uuid),     INTENT(in)  :: uuid1
    TYPE(t_uuid),     INTENT(in)  :: uuid2
    uuid_compare = .TRUE.
!CDIR NOVECTOR
    IF (ANY(uuid1%data /= uuid2%data)) THEN
      uuid_compare = .FALSE.
    ENDIF
  END FUNCTION uuid_compare
  
  SUBROUTINE uuid2char(uuid, string)
    TYPE(t_uuid), INTENT(in) :: uuid
    CHARACTER(len=1), INTENT(out) :: string(16)
    string = TRANSFER(uuid%data, string)
  END SUBROUTINE uuid2char

  SUBROUTINE char2uuid(string, uuid)
    CHARACTER(len=1), INTENT(in) :: string(16)
    TYPE(t_uuid), INTENT(out) :: uuid
    uuid%data  = TRANSFER(string, uuid%data)
  END SUBROUTINE char2uuid

  SUBROUTINE uuid_generate(val, nval, uuid)
    REAL(C_DOUBLE), INTENT(IN)  :: val(*)
    INTEGER(C_INT), INTENT(IN)  :: nval
    TYPE(t_uuid),   INTENT(out) :: uuid
    CALL my_uuid_generate(val, nval, uuid)
  END SUBROUTINE uuid_generate

  INTEGER(C_INT) FUNCTION compare_uuid(uuid_A, uuid_B) 
    TYPE(t_uuid),   INTENT(IN)  :: uuid_A, uuid_B    
    REAL(C_DOUBLE) :: min_difference
    compare_uuid = my_compare_uuid(uuid_A, uuid_B, min_difference) 
  END FUNCTION compare_uuid


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
    INTEGER,  INTENT(INOUT)           :: permutation(:) !< (optional) permutation of indices
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

        t_p            = permutation(i)
        permutation(i) = permutation(j)
        permutation(j) = t_p

        IF (j <= i) EXIT LOOP
      END DO LOOP
      a(j) = a(i)
      a(i) = a(r)
      a(r) = t

      permutation(j) = permutation(i)
      permutation(i) = permutation(r)
      permutation(r) = t_p
      CALL quicksort_int(a,permutation,l,i-1)
      CALL quicksort_int(a,permutation,i+1,r)
    END IF
  END SUBROUTINE quicksort_int


  ! Parallel computation of UUID.
  !
  ! In contrast to the sequential version "uuid_generate", here the
  ! data may be spread over different PEs in an unordered way. Without
  ! using a global array, the UUID is then generated by concatenating
  ! independently computed fingerprints.
  !
  SUBROUTINE uuid_generate_parallel(comm, in_val, in_glbidx, glb_nval, uuid)
    INTEGER,        INTENT(IN)  :: comm
    REAL(C_DOUBLE), INTENT(IN)  :: in_val(:)
    INTEGER,        INTENT(IN)  :: in_glbidx(:)
    INTEGER,        INTENT(IN)  :: glb_nval
    TYPE(t_uuid),   INTENT(OUT) :: uuid
    ! local variables
    INTEGER, PARAMETER :: ROOTPE = 0
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':uuid_generate_parallel'
    TYPE(C_PTR)                          :: fingerprint
    TYPE(C_PTR)                          :: fingerprint_merge, new_fingerprint
    TYPE(t_context), TARGET, ALLOCATABLE :: fingerprint_i(:)
    INTEGER                              :: i, p_error, ipe, npes, target_pe, nval_local, chunksize
    TYPE(t_context), POINTER             :: ptr
    INTEGER, ALLOCATABLE                 :: sendbuf(:), recvbuf(:,:), permutation(:), glbidx_sorted(:), &
      &                                     glbidx_local(:), glbidx(:)
    INTEGER, ALLOCATABLE                 :: sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
    REAL(C_DOUBLE), ALLOCATABLE          :: val_sorted(:), val_local(:)
    INTEGER(C_INT)                       :: nval
    LOGICAL                              :: lcontiguous

    IF (SIZE(in_glbidx) /= SIZE(in_val)) THEN
      WRITE (0,*) routine, ": Size of input fields does not match!"
      RETURN
    END IF

    CALL MPI_COMM_SIZE(comm, npes, p_error)
    CALL MPI_COMM_RANK(comm, ipe,  p_error)

    nval      = SIZE(in_val)
    chunksize = (glb_nval+npes/2)/npes

    ! first, we need to reorder the scattered global indices without
    ! the use of a global-size array. We need to build an
    ! MPI_ALLTOALL, where the global index "glbidx" ends up at PE
    ! glbidx/chunksize:

    ! --- build local copies of "val", "glbidx", sorted by "glbidx":
    ALLOCATE(val_sorted(nval), permutation(nval), glbidx_sorted(nval))
    glbidx = in_glbidx
    permutation = (/ ( i, i=1,nval ) /)
    CALL quicksort_int(glbidx, permutation)
    val_sorted    = in_val(permutation)
    glbidx_sorted = glbidx
    DEALLOCATE(permutation)

    ! --- build a "displacement array", i.e. the first index which is
    ! --- going to a specific PE:
    ALLOCATE(sdispls(npes), sendcounts(npes))
    sdispls       = -1
    sendcounts(:) =  0
    DO i=1,SIZE(glbidx_sorted)
      target_pe = (glbidx_sorted(i)-1)/chunksize +1
      sendcounts(target_pe) = sendcounts(target_pe) + 1
      IF (sdispls(target_pe) == -1)  sdispls(target_pe) = i - 1
    END DO
    WHERE (sdispls(:) < 0)  sdispls(:) = 0

    ! --- each PE tells each other PE how many items it indends to
    ! --- send:
    ALLOCATE(recvcounts(npes), rdispls(npes))
    CALL MPI_ALLTOALL(sendcounts, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, comm, p_error)

    ! --- build a "displacement array" for the receiver side, i.e. the
    ! --- first index where to receive a block from another PE:
    rdispls(:) = 0
    DO i=1,npes-1
      rdispls(i+1) = rdispls(i) + recvcounts(i)
    END DO

    ! --- consistency check: test if our PE is going to receive at
    ! --- most "chunksize" items:
    nval_local = SUM(recvcounts)
    IF (nval_local > chunksize) THEN
      WRITE (0,*) routine, ": Internal error! #local values=", nval_local, "; chunksize=", chunksize
      RETURN
    END IF

    ! --- now perform the ALLTOALLV operation for the data and the
    ! --- global indices:

    ALLOCATE(glbidx_local(nval_local), val_local(nval_local))
    CALL MPI_ALLTOALLV(glbidx_sorted, sendcounts, sdispls, MPI_INTEGER, &
      &                glbidx_local, recvcounts, rdispls, MPI_INTEGER, comm, p_error)
    CALL MPI_ALLTOALLV(val_sorted, sendcounts, sdispls, MPI_DOUBLE, &
      &                val_local, recvcounts, rdispls, MPI_DOUBLE, comm, p_error)

    ALLOCATE(permutation(nval_local))
    permutation = (/ ( i, i=1,nval_local ) /)
    CALL quicksort_int(glbidx_local, permutation)
    val_local    = val_local(permutation)
    DEALLOCATE(permutation)

    ! --- consistency check: test if this PE has received a contiguous
    ! --- chunk of indices:
    lcontiguous = .TRUE.
    DO i=2,nval_local
      IF (glbidx_local(i) /= glbidx_local(i-1)+1) THEN
        lcontiguous = .FALSE.
        EXIT
      END IF
    END DO
    IF (.NOT. lcontiguous) THEN
      WRITE (0,*) routine, ": Internal error, PE has not received a contiguous chunk of indices!"
      RETURN
    END IF

    ! each PE creates "its" part as an independent fingerprint:
    fingerprint = my_uuid_scan_data(val_local, nval_local)

    ! gather all fingerprints on PE#0:
    CALL C_F_POINTER(fingerprint, ptr)
    ALLOCATE(recvbuf(9,npes), sendbuf(9))
    sendbuf(:) = (/ ptr%f64(1),  ptr%f64(2),  &
      &             ptr%f48(1),  ptr%f48(2),  &
      &             ptr%tl64(1), ptr%tl64(2), &
      &             ptr%tl48(1), ptr%tl48(2), &
      &             ptr%max_zero /)
    CALL MPI_GATHER(sendbuf, SIZE(sendbuf), MPI_INTEGER, &
      &             recvbuf, SIZE(sendbuf), MPI_INTEGER, ROOTPE, comm, p_error)

    IF (ipe == ROOTPE) THEN
      ALLOCATE(fingerprint_i(npes))
      DO i=1,npes
        fingerprint_i(i)%f64(1)   = recvbuf(1,i)
        fingerprint_i(i)%f64(2)   = recvbuf(2,i)
        fingerprint_i(i)%f48(1)   = recvbuf(3,i)
        fingerprint_i(i)%f48(2)   = recvbuf(4,i)
        fingerprint_i(i)%tl64(1)  = recvbuf(5,i)
        fingerprint_i(i)%tl64(2)  = recvbuf(6,i)
        fingerprint_i(i)%tl48(1)  = recvbuf(7,i)
        fingerprint_i(i)%tl48(2)  = recvbuf(8,i)
        fingerprint_i(i)%max_zero = recvbuf(9,i)
      END DO
      DEALLOCATE(sendbuf, recvbuf)
      
      ! concatenate them:
      fingerprint_merge = C_LOC(fingerprint_i(1))
      DO i=2,npes
        new_fingerprint   = my_concat_fingerprints(fingerprint_merge, C_LOC(fingerprint_i(i)));
        IF (i>2)  CALL deallocate_c(fingerprint_merge)
        fingerprint_merge = new_fingerprint
      END DO
      
      CALL my_encode_uuid(fingerprint_merge, uuid)
      CALL deallocate_c(fingerprint_merge)
      DEALLOCATE(fingerprint_i)
    END IF
    DEALLOCATE(val_sorted, glbidx_sorted, sdispls, sendcounts, &
      &        glbidx_local, val_local, recvcounts, rdispls, glbidx)
  END SUBROUTINE uuid_generate_parallel


  SUBROUTINE clear_uuid(uuid)
    TYPE(t_uuid), INTENT(inout) :: uuid
    uuid%data(:) = INT(0, C_SIGNED_CHAR)
  END SUBROUTINE clear_uuid

END MODULE mo_util_uuid
