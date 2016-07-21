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

  USE mo_kind,                    ONLY: wp, sp, dp, ik4 => i4
  USE mo_exception,               ONLY: finish
  USE mo_impl_constants,          ONLY: VARNAME_LEN

  IMPLICIT NONE

  PUBLIC :: t_Destructible
  PUBLIC :: assign_if_present
  PUBLIC :: t_ptr_2d3d
  PUBLIC :: t_ptr_i2d3d
  PUBLIC :: t_ptr_tracer
  PUBLIC :: copy, init, swap, var_scale, negative2zero
  PUBLIC :: init_zero_contiguous_dp, init_zero_contiguous_sp
  PUBLIC :: resize_arr_c1d

  PRIVATE


  !> Just a small base CLASS for anything that needs a destructor.
  !> XXX: This will become unnecessary once all relevant compilers support the FINAL keyword.
  TYPE, ABSTRACT :: t_Destructible
  CONTAINS
    PROCEDURE(interface_destructor), DEFERRED :: destruct
  END TYPE t_Destructible

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

  INTERFACE assign_if_present
    MODULE PROCEDURE assign_if_present_character
    MODULE PROCEDURE assign_if_present_logical
    MODULE PROCEDURE assign_if_present_logicals
    MODULE PROCEDURE assign_if_present_integer
    MODULE PROCEDURE assign_if_present_integers
    MODULE PROCEDURE assign_if_present_real
  END INTERFACE assign_if_present

  !> this is meant to make it easier for compilers to circumvent
  !! temporaries as are too often created in a(:, :, :) = b(:, :, :)
  !! uses omp orphaning
  INTERFACE copy
    MODULE PROCEDURE copy_2d_dp
    MODULE PROCEDURE copy_3d_dp
    MODULE PROCEDURE copy_4d_dp
    MODULE PROCEDURE copy_5d_dp
    MODULE PROCEDURE copy_5d_sp
    MODULE PROCEDURE copy_5d_spdp
    MODULE PROCEDURE copy_2d_i4
    MODULE PROCEDURE copy_3d_i4
    MODULE PROCEDURE copy_5d_i4
  END INTERFACE copy

  INTERFACE init
    MODULE PROCEDURE init_zero_1d_dp
    MODULE PROCEDURE init_zero_2d_dp
    MODULE PROCEDURE init_zero_2d_i4
    MODULE PROCEDURE init_zero_3d_dp
    MODULE PROCEDURE init_zero_3d_sp
    MODULE PROCEDURE init_zero_3d_i4
    MODULE PROCEDURE init_zero_4d_i4
    MODULE PROCEDURE init_zero_4d_dp
    MODULE PROCEDURE init_zero_4d_sp
    MODULE PROCEDURE init_3d_dp
    MODULE PROCEDURE init_3d_spdp
    MODULE PROCEDURE init_5d_dp
    MODULE PROCEDURE init_5d_i4
  END INTERFACE init

  INTERFACE negative2zero
    MODULE PROCEDURE negative2zero_4d_dp
  END INTERFACE negative2zero

  INTERFACE var_scale
    MODULE PROCEDURE var_scale_3d_dp
  END INTERFACE var_scale

  INTERFACE swap
    MODULE PROCEDURE swap_int
  END INTERFACE swap

  ABSTRACT INTERFACE
    !> destructor interface
    SUBROUTINE interface_destructor(me)
        IMPORT t_Destructible
        CLASS(t_Destructible), INTENT(INOUT) :: me
    END SUBROUTINE interface_destructor
  END INTERFACE

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

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_2d_dp(src, dest)
    REAL(dp), INTENT(in) :: src(:, :)
    REAL(dp), INTENT(out) :: dest(:, :)
    INTEGER :: i1, i2, m1, m2
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
!$omp do collapse(2)
    DO i2 = 1, m2
      DO i1 = 1, m1
        dest(i1, i2) = src(i1, i2)
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE copy_2d_dp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_3d_dp(src, dest)
    REAL(dp), INTENT(in) :: src(:, :, :)
    REAL(dp), INTENT(out) :: dest(:, :, :)
    INTEGER :: i1, i2, i3, m1, m2, m3
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
#ifdef _CRAYFTN
!$omp do
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          dest(i1, i2, i3) = src(i1, i2, i3)
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE copy_3d_dp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_4d_dp(src, dest)
    REAL(dp), INTENT(in) :: src(:, :, :, :)
    REAL(dp), INTENT(out) :: dest(:, :, :, :)
    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
    m4 = SIZE(dest, 4)
#ifdef _CRAYFTN
!$omp do
#else
!$omp do collapse(4)
#endif
    DO i4 = 1, m4
      DO i3 = 1, m3
        DO i2 = 1, m2
          DO i1 = 1, m1
            dest(i1, i2, i3, i4) = src(i1, i2, i3, i4)
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE copy_4d_dp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_5d_dp(src, dest)
    REAL(dp), INTENT(in) :: src(:, :, :, :, :)
    REAL(dp), INTENT(out) :: dest(:, :, :, :, :)
    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
    m4 = SIZE(dest, 4)
    m5 = SIZE(dest, 5)
!$omp do collapse(5)
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              dest(i1, i2, i3, i4, i5) = src(i1, i2, i3, i4, i5)
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE copy_5d_dp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_5d_sp(src, dest)
    REAL(sp), INTENT(in) :: src(:, :, :, :, :)
    REAL(sp), INTENT(out) :: dest(:, :, :, :, :)
    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
    m4 = SIZE(dest, 4)
    m5 = SIZE(dest, 5)
!$omp do collapse(5)
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              dest(i1, i2, i3, i4, i5) = src(i1, i2, i3, i4, i5)
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE copy_5d_sp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_5d_spdp(src, dest)
    REAL(sp), INTENT(in) :: src(:, :, :, :, :)
    REAL(dp), INTENT(out) :: dest(:, :, :, :, :)
    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
    m4 = SIZE(dest, 4)
    m5 = SIZE(dest, 5)
!$omp do collapse(5)
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              dest(i1, i2, i3, i4, i5) = REAL(src(i1, i2, i3, i4, i5),KIND=dp)
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE copy_5d_spdp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_2d_i4(src, dest)
    INTEGER(ik4), INTENT(in) :: src(:, :)
    INTEGER(ik4), INTENT(out) :: dest(:, :)
    INTEGER :: i1, i2, m1, m2
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
!$omp do collapse(2)
    DO i2 = 1, m2
      DO i1 = 1, m1
        dest(i1, i2) = src(i1, i2)
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE copy_2d_i4

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_5d_i4(src, dest)
    INTEGER(ik4), INTENT(in) :: src(:, :, :, :, :)
    INTEGER(ik4), INTENT(out) :: dest(:, :, :, :, :)
    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
    m4 = SIZE(dest, 4)
    m5 = SIZE(dest, 5)
!$omp do collapse(5)
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              dest(i1, i2, i3, i4, i5) = src(i1, i2, i3, i4, i5)
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE copy_5d_i4

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_3d_i4(src, dest)
    INTEGER(ik4), INTENT(in) :: src(:, :, :)
    INTEGER(ik4), INTENT(out) :: dest(:, :, :)
    INTEGER :: i1, i2, i3, m1, m2, m3
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
!$omp do collapse(3)
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          dest(i1, i2, i3) = src(i1, i2, i3)
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE copy_3d_i4


  SUBROUTINE init_zero_3d_dp(init_var)
    REAL(dp), INTENT(out) :: init_var(:, :, :)
    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
#ifdef _CRAYFTN
!$omp do
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = 0.0_dp
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE init_zero_3d_dp

  SUBROUTINE init_zero_4d_dp(init_var)
    REAL(dp), INTENT(out) :: init_var(:, :, :, :)
    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)
#ifdef _CRAYFTN
!$omp do
#else
!$omp do collapse(4)
#endif
    DO i4 = 1, m4
      DO i3 = 1, m3
        DO i2 = 1, m2
          DO i1 = 1, m1
            init_var(i1, i2, i3, i4) = 0.0_dp
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE init_zero_4d_dp

  SUBROUTINE init_zero_4d_sp(init_var)
    REAL(sp), INTENT(out) :: init_var(:, :, :, :)
    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)
#ifdef _CRAYFTN
!$omp do
#else
!$omp do collapse(4)
#endif
    DO i4 = 1, m4
      DO i3 = 1, m3
        DO i2 = 1, m2
          DO i1 = 1, m1
            init_var(i1, i2, i3, i4) = 0.0_sp
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE init_zero_4d_sp

  SUBROUTINE init_zero_3d_sp(init_var)
    REAL(sp), INTENT(out) :: init_var(:, :, :)
    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
!$omp do collapse(3)
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = 0.0_sp
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE init_zero_3d_sp

  SUBROUTINE init_zero_2d_i4(init_var)
    INTEGER(ik4), INTENT(out) :: init_var(:, :)
    INTEGER :: i1, i2, m1, m2

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
!$omp do collapse(2)
    DO i2 = 1, m2
      DO i1 = 1, m1
        init_var(i1, i2) = 0_ik4
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE init_zero_2d_i4

  SUBROUTINE init_zero_3d_i4(init_var)
    INTEGER(ik4), INTENT(out) :: init_var(:, :, :)
    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
!$omp do collapse(3)
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = 0_ik4
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE init_zero_3d_i4

  SUBROUTINE init_zero_4d_i4(init_var)
    INTEGER(ik4), INTENT(out) :: init_var(:, :, :, :)
    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)
!$omp do collapse(4)
    DO i4 = 1, m4
      DO i3 = 1, m3
        DO i2 = 1, m2
          DO i1 = 1, m1
            init_var(i1, i2, i3, i4) = 0_ik4
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE init_zero_4d_i4

  SUBROUTINE init_zero_1d_dp(init_var)
    REAL(dp), INTENT(out) :: init_var(:)
    INTEGER :: i1, m1

    m1 = SIZE(init_var, 1)
!$omp do
    DO i1 = 1, m1
      init_var(i1) = 0.0_dp
    END DO
!$omp end do nowait
  END SUBROUTINE init_zero_1d_dp

  SUBROUTINE init_zero_2d_dp(init_var)
    REAL(dp), INTENT(out) :: init_var(:, :)
    INTEGER :: i1, i2, m1, m2

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
!$omp do collapse(2)
    DO i2 = 1, m2
      DO i1 = 1, m1
        init_var(i1, i2) = 0.0_dp
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE init_zero_2d_dp

  SUBROUTINE init_3d_dp(init_var, init_val)
    REAL(dp), INTENT(out) :: init_var(:, :, :)
    REAL(dp), INTENT(in) :: init_val

    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
!$omp do collapse(3)
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = init_val
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE init_3d_dp

  SUBROUTINE init_3d_spdp(init_var, init_val)
    REAL(sp), INTENT(out) :: init_var(:, :, :)
    REAL(dp), INTENT(in) :: init_val

    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
!$omp do collapse(3)
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = REAL(init_val,KIND=sp)
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE init_3d_spdp

  SUBROUTINE init_5d_dp(init_var, init_val)
    REAL(dp), INTENT(out) :: init_var(:, :, :, :, :)
    REAL(dp), INTENT(in) :: init_val

    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)
    m5 = SIZE(init_var, 5)
!$omp do collapse(5)
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              init_var(i1, i2, i3, i4, i5) = init_val
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE init_5d_dp

  SUBROUTINE init_5d_i4(init_var, init_val)
    INTEGER(ik4), INTENT(out) :: init_var(:, :, :, :, :)
    INTEGER(ik4), INTENT(in) :: init_val

    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)
    m5 = SIZE(init_var, 5)
!$omp do collapse(5)
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              init_var(i1, i2, i3, i4, i5) = init_val
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE init_5d_i4


  SUBROUTINE var_scale_3d_dp(var, scale_val)
    REAL(dp), INTENT(inout) :: var(:, :, :)
    REAL(dp), INTENT(in) :: scale_val

    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(var, 1)
    m2 = SIZE(var, 2)
    m3 = SIZE(var, 3)
!$omp do collapse(3)
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          var(i1, i2, i3) = var(i1, i2, i3) * scale_val
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE var_scale_3d_dp

  SUBROUTINE negative2zero_4d_dp(var)
    REAL(dp), INTENT(inout) :: var(:, :, :, :)

    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4
    REAL(dp) :: v

    m1 = SIZE(var, 1)
    m2 = SIZE(var, 2)
    m3 = SIZE(var, 3)
    m4 = SIZE(var, 4)
!$omp do collapse(4)
    DO i4 = 1, m4
      DO i3 = 1, m3
        DO i2 = 1, m2
          DO i1 = 1, m1
            v = var(i1, i2, i3, i4)
            var(i1, i2, i3, i4) = (ABS(v) + v) * 0.5_dp
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE negative2zero_4d_dp

  SUBROUTINE init_zero_contiguous_dp(var, n)
    INTEGER, INTENT(in) :: n
    REAL(dp), INTENT(out) :: var(n)

    INTEGER :: i
!$omp do
    DO i = 1, n
      var(i) = 0.0_dp
    END DO
!$omp end do nowait
  END SUBROUTINE init_zero_contiguous_dp

  SUBROUTINE init_zero_contiguous_sp(var, n)
    INTEGER, INTENT(in) :: n
    REAL(sp), INTENT(out) :: var(n)

    INTEGER :: i
!$omp do
    DO i = 1, n
      var(i) = 0.0_sp
    END DO
!$omp end do nowait
  END SUBROUTINE init_zero_contiguous_sp


END MODULE mo_fortran_tools
