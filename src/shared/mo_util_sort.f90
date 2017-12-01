!! Module containing basic sorting algorithms (required until
!! a Fortran STL has been invented).
!!
!! Initial revision: F. Prill, DWD (10/2012)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_util_sort
  USE mo_exception, ONLY: finish
  USE mo_impl_constants, ONLY: SUCCESS

#ifdef __ICON__
  USE mo_kind,   ONLY: wp
#else
  USE mo_utilities, ONLY: wp
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: quicksort
  PUBLIC :: insertion_sort
  PUBLIC :: t_Permutation

  ! generic interface for in-situ QuickSort sorting routine
  INTERFACE quicksort
    MODULE PROCEDURE quicksort_real
    MODULE PROCEDURE quicksort_int
  END INTERFACE

  INTERFACE insertion_sort
    MODULE PROCEDURE insertion_sort_int
  END INTERFACE insertion_sort

  ! A simple tool to reorder DATA.
  TYPE t_Permutation
    INTEGER, ALLOCATABLE :: indexTranslation(:), inverseTranslation(:)
  CONTAINS
    PROCEDURE :: construct => permutation_construct
    PROCEDURE :: permute => permutation_permute
    PROCEDURE :: reverse => permutation_reverse
    PROCEDURE :: destruct => permutation_destruct
  END TYPE t_Permutation

  CHARACTER(*), PARAMETER :: modname = "mo_util_sort"

CONTAINS

  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of REAL values.
  ! 
  !  Ordering after the sorting process: smallest...largest.
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

  ! Constructs a permutation on the length of the supplied INTEGER
  ! array such that permutating that same array results IN a sorted
  ! array.  There IS no requirement on the supplied array itself: It
  ! may USE non-consecutive numbers AND/OR contain repetitions, AND
  ! the number range IS NOT restricted IN ANY way.  IF two OR more
  ! entries of the supplied array are equal, the resulting FUNCTION IS
  ! NOT injective: permute() will ignore all but one of each input
  ! corresponding to the same order VALUE, AND its output array IS
  ! expected to be smaller than the input array.  reverse() will
  ! output more points than it has input by multiplexing each entry to
  ! all positions that correspond to its VALUE IN the order array.
  !     Thus, the sequence
  !         CALL permutation%reverse(a, b);
  !         CALL permutation%permute(b, c);
  !     IS guaranteed to RESULT IN `a == c`. However, the sequence
  !         CALL permutation%permute(a, b);
  !         CALL permutation%reverse(b, c);
  !     does NOT provide the same guarantee as information IS lost within the smaller array `b`.
  SUBROUTINE permutation_construct(me, order)
    CLASS(t_Permutation), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: order(:)

    INTEGER :: inputSize, outputSize, i, error
    INTEGER, ALLOCATABLE :: sorted(:), inversePermutation(:), outputReductionTranslation(:)
    CHARACTER(*), PARAMETER :: routine = modname//":permutation_construct"

    ! ALLOCATE memory AND handle the trivial CASE that the permutation SIZE IS zero
    inputSize = SIZE(order)
    ALLOCATE(me%indexTranslation(inputSize), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    IF(inputSize == 0) THEN
        ! normaly, this IS ALLOCATED later, but we should ALLOCATE it before returning
        ALLOCATE(me%inverseTranslation(inputSize), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        RETURN
    END IF

    ALLOCATE(inversePermutation(inputSize), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    ALLOCATE(sorted(inputSize), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    ALLOCATE(outputReductionTranslation(inputSize), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

    ! compute an inverse permutation first as this IS easily achieved using the quicksort() FUNCTION
    DO i = 1, inputSize
        inversePermutation(i) = i
    END DO
    sorted(:) = order(:)
    CALL quicksort(sorted, inversePermutation)

    ! derive a forward permutation from the inverse translation
    ! this IS NOT the FINAL forward translation as that has to remove the duplicates
    me%indexTranslation(:) = -1
    DO i = 1, inputSize
        me%indexTranslation(inversePermutation(i)) = i
    END DO

    ! assert that every element of indexTranslation was initialized correctly
    IF(ANY(me%indexTranslation < 1)) CALL finish(routine, "assertion failed")

    ! determine how we have to reduce the output
    outputReductionTranslation(1) = 1
    DO i = 2, inputSize
        IF(sorted(i-1) == sorted(i)) THEN
            outputReductionTranslation(i) = outputReductionTranslation(i-1)
        ELSE
            outputReductionTranslation(i) = outputReductionTranslation(i-1) + 1
        END IF
    END DO
    outputSize = outputReductionTranslation(inputSize)

    ! update the forward translation accordingly
    DO i = 1, inputSize
        me%indexTranslation(i) = outputReductionTranslation(me%indexTranslation(i))
    END DO

    ! derive the backward translation from the forward translation
    ALLOCATE(me%inverseTranslation(outputSize), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    me%inverseTranslation(:) = -1
    DO i = 1, inputSize
        me%inverseTranslation(me%indexTranslation(i)) = i
    END DO

    ! assert that every element of inverseTranslation was initialized correctly
    IF(ANY(me%inverseTranslation < 1)) CALL finish(routine, "assertion failed")

    ! cleanup
    DEALLOCATE(sorted, inversePermutation, outputReductionTranslation)
  END SUBROUTINE permutation_construct

  SUBROUTINE permutation_permute(me, input, output)
    CLASS(t_Permutation), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: input(:)
    INTEGER, INTENT(OUT) :: output(:)

    INTEGER :: i
    CHARACTER(*), PARAMETER :: routine = modname//":permutation_permute"
    CHARACTER(len=1000) :: message_text = ''

    ! check preconditions
    IF(SIZE(input) /= SIZE(me%indexTranslation)) THEN
      WRITE (message_text, *)  "illegal argument: input size (= ", SIZE(input), ") ", &
        & "must match input size of the permutation (= ", SIZE(me%indexTranslation), ")"
      CALL finish(routine, message_text)
    END IF
    IF(SIZE(output) /= SIZE(me%inverseTranslation)) THEN
      WRITE (message_text, *)  "illegal argument: output size (= ", SIZE(output), ") ", &
        & "must match output size of the permutation (= ", SIZE(me%inverseTranslation), &
        & ") (which may be smaller than the input size of the permutation)"
      CALL finish(routine, message_text)
    END IF

    DO i = 1, SIZE(output)
        output(i) = input(me%inverseTranslation(i))
    END DO
  END SUBROUTINE permutation_permute

  SUBROUTINE permutation_reverse(me, input, output)
    CLASS(t_Permutation), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: input(:)
    INTEGER, INTENT(OUT) :: output(:)

    INTEGER :: i
    CHARACTER(*), PARAMETER :: routine = modname//":permutation_reverse"
    CHARACTER(len=1000) :: message_text = ''

    ! check preconditions
    IF(SIZE(output) /= SIZE(me%indexTranslation)) THEN
      WRITE (message_text, *) "illegal argument: output size (= ", SIZE(output), ") ", &
        & "must match input size of the permutation (= ", SIZE(me%indexTranslation), ")"
      CALL finish(routine, message_text)
    END IF
    IF(SIZE(input) /= SIZE(me%inverseTranslation)) THEN
      WRITE (message_text, *) "illegal argument: input size (= ", SIZE(input), ") ", &
        & "must match output size of the permutation (= ", SIZE(me%inverseTranslation), &
        & ") (which may be smaller than the input size of the permutation)"
      CALL finish(routine, message_text)
    END IF

    DO i = 1, SIZE(output)
        output(i) = input(me%indexTranslation(i))
    END DO
  END SUBROUTINE permutation_reverse

  SUBROUTINE permutation_destruct(me)
    CLASS(t_Permutation), INTENT(INOUT) :: me

    DEALLOCATE(me%indexTranslation, me%inverseTranslation)
  END SUBROUTINE permutation_destruct

END MODULE mo_util_sort
