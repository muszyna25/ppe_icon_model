!! Module containing basic sorting algorithms (required until
!! a Fortran STL has been invented).
!!
!! Initial revision: F. Prill, DWD (10/2012)
!!
MODULE mo_util_sort

  USE mo_kind,   ONLY: wp
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: quicksort

  ! generic interface for in-situ QuickSort sorting routine
  INTERFACE quicksort
    MODULE PROCEDURE quicksort_real
    MODULE PROCEDURE quicksort_int
  END INTERFACE
  
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


  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of INTEGER values.
  !
  RECURSIVE SUBROUTINE quicksort_int(a, permutation, l_in, r_in)
    INTEGER,  INTENT(INOUT)           :: a(:)           !< array for in-situ sorting
    INTEGER,  INTENT(INOUT), OPTIONAL :: permutation(:) !< (optional) permutation of indices
    INTEGER , INTENT(IN),    OPTIONAL :: l_in,r_in      !< left, right partition indices
    ! local variables
    INTEGER :: i,j,l,r,t_p,t,v

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
  END SUBROUTINE quicksort_int

END MODULE mo_util_sort
