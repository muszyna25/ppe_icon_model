!! Module containing basic sorting algorithms (required until
!! a Fortran STL has been invented).
!!
!! Initial revision: F. Prill, DWD (10/2012)
!!
MODULE mo_util_sort

#ifdef __ICON__
  USE mo_kind,   ONLY: wp
#else
  USE mo_utilities, ONLY: wp
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: quicksort
  PUBLIC :: insertion_sort

  ! generic interface for in-situ QuickSort sorting routine
  INTERFACE quicksort
    MODULE PROCEDURE quicksort_real
    MODULE PROCEDURE quicksort_int
  END INTERFACE

  INTERFACE insertion_sort
    MODULE PROCEDURE insertion_sort_int
  END INTERFACE insertion_sort

CONTAINS

  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of REAL values.
  !
  RECURSIVE SUBROUTINE quicksort_real(a, permutation, l_in, r_in)
    REAL(wp), INTENT(INOUT)           :: a(:)           !< array for in-situ sorting
    INTEGER,  INTENT(INOUT), OPTIONAL :: permutation(:) !< (optional) permutation of indices
    INTEGER , INTENT(IN),    OPTIONAL :: l_in,r_in      !< left, right partition indices
    ! local variables
    INTEGER  :: i,j,l,r,t_p
    REAL(wp) :: t,v

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
      v   = a(r)
      i   = l-1
      j   = r
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
      CALL quicksort(a,permutation,l,i-1)
      CALL quicksort(a,permutation,i+1,r)
    END IF
  END SUBROUTINE quicksort_real


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
      CALL quicksort(a,permutation,l,i-1)
      CALL quicksort(a,permutation,i+1,r)
    END IF
  END SUBROUTINE quicksort_int

  SUBROUTINE insertion_sort_int(a)
    INTEGER, INTENT(inout) :: a(:)

    INTEGER :: t, h
    INTEGER :: i, n

    n = SIZE(a)

    DO i = 2, n
      t = a(i)
      DO h = i - 1, 1, -1
        IF (t >= a(h)) EXIT
        a(h + 1) = a(h)
      END DO
      a(h + 1) = t
    END DO
  END SUBROUTINE insertion_sort_int

END MODULE mo_util_sort
