!>
!! This module contains often-used Fortran language constructs.
!!
!! The small functions and subroutines in this module should depend
!! only on most basic types and should not call other model-specific
!! subroutines.
!!
!! @par Revision History
!!    Initial implementation : F. Prill, DWD (2012-07-04)
!!    moved routines from "mo_var_list.f90"
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_fortran_tools

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: finish
  USE mo_var_metadata_types,      ONLY: VARNAME_LEN

  IMPLICIT NONE
  PRIVATE


  TYPE t_ptr_2d3d
    REAL(wp),POINTER :: p_3d(:,:,:)  ! REAL pointer to 3D (spatial) array
    REAL(wp),POINTER :: p_2d(:,:)    ! REAL pointer to 2D (spatial) array
  END TYPE t_ptr_2d3d

  TYPE t_ptr_i2d3d
    INTEGER,POINTER :: p_3d(:,:,:)  ! INTEGER pointer to 3D (spatial) array
    INTEGER,POINTER :: p_2d(:,:)    ! INTEGER pointer to 2D (spatial) array
  END TYPE t_ptr_i2d3d

  ! Type to pass pointer arrays to convection and turbulent diffusion subroutines
  TYPE t_ptr_tracer
    REAL(wp), POINTER :: ptr(:,:)
    INTEGER           :: idx_tracer
  END TYPE t_ptr_tracer

  PUBLIC :: assign_if_present
  PUBLIC :: t_ptr_2d3d
  PUBLIC :: t_ptr_i2d3d
  PUBLIC :: t_ptr_tracer
  PUBLIC :: swap
  PUBLIC :: resize_arr_c1d

  INTERFACE assign_if_present
    MODULE PROCEDURE assign_if_present_character
    MODULE PROCEDURE assign_if_present_logical
    MODULE PROCEDURE assign_if_present_logicals
    MODULE PROCEDURE assign_if_present_integer
    MODULE PROCEDURE assign_if_present_integers
    MODULE PROCEDURE assign_if_present_real
  END INTERFACE assign_if_present

  INTERFACE swap
    MODULE PROCEDURE swap_int
  END INTERFACE swap

CONTAINS

  ! private routines to assign values if actual parameters are present
  !
  SUBROUTINE assign_if_present_character (y,x)
    CHARACTER(len=*), INTENT(inout)        :: y
    CHARACTER(len=*), INTENT(in) ,OPTIONAL :: x
    IF (.NOT. PRESENT(x)) RETURN
    IF ( x == ' ' )       RETURN
    y = x
  END SUBROUTINE assign_if_present_character


  SUBROUTINE assign_if_present_logical (y,x)
    LOGICAL, INTENT(inout)        :: y
    LOGICAL, INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_logical


  SUBROUTINE assign_if_present_logicals (y,x)
    LOGICAL, INTENT(inout)        :: y(:)
    LOGICAL, INTENT(in) ,OPTIONAL :: x(:)
    INTEGER :: n
    IF (PRESENT(x)) THEN
      n = MIN(SIZE(x), SIZE(y))
      y(1:n) = x(1:n)
    ENDIF
  END SUBROUTINE assign_if_present_logicals


  SUBROUTINE assign_if_present_integer (y,x)
    INTEGER, INTENT(inout)        :: y
    INTEGER, INTENT(in) ,OPTIONAL :: x
    IF (.NOT. PRESENT(x)) RETURN
    IF ( x == -HUGE(x)  ) RETURN
    y = x
  END SUBROUTINE assign_if_present_integer


  SUBROUTINE assign_if_present_integers (y,x)
    INTEGER, INTENT(inout)        :: y (:)
    INTEGER, INTENT(in) ,OPTIONAL :: x (:)
    INTEGER :: n
    IF (PRESENT(x)) THEN
      n = MIN(SIZE(x), SIZE(y))
      y(1:n) = x(1:n)
    ENDIF
  END SUBROUTINE assign_if_present_integers


  SUBROUTINE assign_if_present_real (y,x)
    REAL(wp), INTENT(inout)        :: y
    REAL(wp), INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    IF ( x == -HUGE(x) ) RETURN
    y = x
  END SUBROUTINE assign_if_present_real


  !>
  !! Swap content of two Integers
  !!
  !! Swap content of two Integers
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-10-28)
  !!
  SUBROUTINE swap_int(a,b)
    INTEGER, INTENT(INOUT) :: a
    INTEGER, INTENT(INOUT) :: b

    ! local variables
    INTEGER :: temp
  !-----------------------------
    temp = a
    a    = b
    b    = temp
  END SUBROUTINE swap_int



  !>
  !! Expand array by given size
  !!
  !! Expand a 1D character array by given size.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-10-28)
  !!
  SUBROUTINE resize_arr_c1d(arr,nelem)
    ! GCC 4.9.0 complained about CHARACTER(:); Cray did not!
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE, INTENT(INOUT) :: arr(:)   ! array to be resized
    INTEGER                  , INTENT(IN)    :: nelem    ! number of elements to expand
    !
    ! local variables
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: tmp_arr(:)
    INTEGER :: istat                   ! status
    !-----------------------------

    ! If arr has not yet been allocated, do it.
    IF (.NOT. ALLOCATED(arr)) THEN
      ALLOCATE(arr(1), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_fortran_tools:resize_arr_c1d', &
             &      ' initial allocation of array arr failed')
      ENDIF
    ELSE
      ! check for appropriate nelem
      IF ( nelem < 0) CALL finish('mo_fortran_tools:resize_arr_c1d', &
        &                         ' nelem must be > 0')

      ! allocate temporary array of size SIZE(arr)+nelem
      ALLOCATE(tmp_arr(1:SIZE(arr)+nelem), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_fortran_tools:resize_char_arr_1d', &
             &      'allocation of array tmp_arr failed')
      ENDIF

      ! copy 
      tmp_arr(1:SIZE(arr)) = arr(1:SIZE(arr))

      CALL move_alloc(tmp_arr, arr)
      ! now arr has been resized to the size of tmp_arr, 
      ! and tmp_arr is deallocated.
    ENDIF

  END SUBROUTINE resize_arr_c1d

END MODULE mo_fortran_tools
