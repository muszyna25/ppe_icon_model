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

  USE mo_kind,                    ONLY: wp, sp, vp, dp, ik4 => i4
  USE mo_exception,               ONLY: finish
  USE mo_impl_constants,          ONLY: SUCCESS
  USE mo_impl_constants,          ONLY: VARNAME_LEN
#ifdef _OPENACC
  USE mo_mpi,                     ONLY: i_am_accel_node
#endif
  USE iso_c_binding, ONLY: c_ptr, c_f_pointer, c_loc

  IMPLICIT NONE

  PUBLIC :: t_Destructible
  PUBLIC :: assign_if_present
  PUBLIC :: t_ptr_2d3d, t_ptr_2d3d_vp
  PUBLIC :: assign_if_present_allocatable
  PUBLIC :: alloc
  PUBLIC :: ensureSize
  PUBLIC :: t_alloc_character
  PUBLIC :: t_ptr_1d
  PUBLIC :: t_ptr_1d_ptr_1d
  PUBLIC :: t_ptr_2d, t_ptr_2d_sp, t_ptr_2d_int
  PUBLIC :: t_ptr_3d, t_ptr_3d_sp
  PUBLIC :: t_ptr_i2d3d
  PUBLIC :: t_ptr_tracer
  PUBLIC :: copy, init, swap, negative2zero
  PUBLIC :: var_scale, var_add
  PUBLIC :: init_zero_contiguous_dp, init_zero_contiguous_sp
  PUBLIC :: resize_arr_c1d
  PUBLIC :: DO_DEALLOCATE
  PUBLIC :: DO_PTR_DEALLOCATE
  PUBLIC :: insert_dimension

  PRIVATE

  !> Just a small base CLASS for anything that needs a destructor.
  !> XXX: This will become unnecessary once all relevant compilers support the FINAL keyword.
  TYPE, ABSTRACT :: t_Destructible
  CONTAINS
    PROCEDURE(interface_destructor), DEFERRED :: destruct
  END TYPE t_Destructible

  TYPE t_alloc_character
    CHARACTER(:), ALLOCATABLE :: a
  END TYPE t_alloc_character

  TYPE t_ptr_1d
    REAL(wp),POINTER :: p(:)  ! pointer to 1D (spatial) array
  END TYPE t_ptr_1d

  TYPE t_ptr_1d_ptr_1d
    TYPE(t_ptr_1d), POINTER :: p(:)  ! pointer to a 1D array of pointers to 1D (spatial) arrays
  END TYPE t_ptr_1d_ptr_1d

  TYPE t_ptr_2d
    REAL(dp),POINTER :: p(:,:)  ! pointer to 2D (spatial) array
  END TYPE t_ptr_2d

  TYPE t_ptr_2d_sp
    REAL(sp),POINTER :: p(:,:)  ! pointer to 2D (spatial) array
  END TYPE t_ptr_2d_sp

  TYPE t_ptr_2d_int
    INTEGER,POINTER :: p(:,:)  ! pointer to 2D (spatial) array
  END TYPE t_ptr_2d_int

  TYPE t_ptr_3d
    REAL(dp),POINTER :: p(:,:,:)  ! pointer to 3D (spatial) array
  END TYPE t_ptr_3d

  TYPE t_ptr_3d_sp
    REAL(sp),POINTER :: p(:,:,:)  ! pointer to 3D (spatial) array
  END TYPE t_ptr_3d_sp

  TYPE t_ptr_2d3d
    REAL(wp),POINTER :: p_3d(:,:,:)  ! REAL pointer to 3D (spatial) array
    REAL(wp),POINTER :: p_2d(:,:)    ! REAL pointer to 2D (spatial) array
  END TYPE t_ptr_2d3d

  TYPE t_ptr_2d3d_vp
    REAL(vp),POINTER :: p_3d(:,:,:)  ! REAL pointer to 3D (spatial) array
    REAL(vp),POINTER :: p_2d(:,:)    ! REAL pointer to 2D (spatial) array
  END TYPE t_ptr_2d3d_vp


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
    MODULE PROCEDURE assign_if_present_real_sp
  END INTERFACE assign_if_present

  INTERFACE assign_if_present_allocatable
    MODULE PROCEDURE assign_if_present_logical_allocatable_1d
    MODULE PROCEDURE assign_if_present_integer_allocatable
    MODULE PROCEDURE assign_if_present_integer_allocatable_1d
    MODULE PROCEDURE assign_if_present_real_allocatable
    MODULE PROCEDURE assign_if_present_real_allocatable_1d
  END INTERFACE assign_if_present_allocatable

  ! This allocates an array, adjusting the allocation SIZE to 1 IF the given SIZE IS zero OR less, AND checking for allocation failure.
  INTERFACE alloc
    MODULE PROCEDURE alloc_int_1d
    MODULE PROCEDURE alloc_double_1d
    MODULE PROCEDURE alloc_single_1d
  END INTERFACE alloc

  ! This handles the recuring CASE of growing a buffer to match possibly increasing needs.
  ! We USE a POINTER to pass the buffer because that allows us to avoid an extra copy when reallocating the buffer.
  ! The association status of the POINTER that IS passed IN must be defined.
  INTERFACE ensureSize
    MODULE PROCEDURE ensureSize_dp_1d
  END INTERFACE ensureSize

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

  INTERFACE var_add
    MODULE PROCEDURE var_addc_3d_dp
  END INTERFACE var_add

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


  ! auxiliary routines
  INTERFACE DO_DEALLOCATE
    MODULE PROCEDURE DO_DEALLOCATE_r4D
    MODULE PROCEDURE DO_DEALLOCATE_r3D
    MODULE PROCEDURE DO_DEALLOCATE_r2D
    MODULE PROCEDURE DO_DEALLOCATE_r1D
    MODULE PROCEDURE DO_DEALLOCATE_i3D
    MODULE PROCEDURE DO_DEALLOCATE_i2D
  END INTERFACE DO_DEALLOCATE

  INTERFACE DO_PTR_DEALLOCATE
    MODULE PROCEDURE DO_PTR_DEALLOCATE_r3D
    MODULE PROCEDURE DO_PTR_DEALLOCATE_r2D
  END INTERFACE DO_PTR_DEALLOCATE


  CHARACTER(LEN = *), PARAMETER :: modname = "mo_fortran_tools"

#if defined( _OPENACC )
#if defined(__FORTRAN_TOOLS_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
  LOGICAL, PARAMETER ::  acc_validate = .FALSE.   ! Only activate for validation phase
#endif


  INTERFACE insert_dimension
    MODULE PROCEDURE insert_dimension_r_wp_6_5, insert_dimension_r_wp_6_5_s
    MODULE PROCEDURE insert_dimension_r_sp_6_5, insert_dimension_r_sp_6_5_s
    MODULE PROCEDURE insert_dimension_i4_6_5, insert_dimension_i4_6_5_s
  END INTERFACE insert_dimension
CONTAINS

  ! routines to assign values if actual parameters are present
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

  SUBROUTINE assign_if_present_real_sp (y,x)
    REAL(sp), INTENT(inout)        :: y
    REAL(sp), INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    IF ( x == -HUGE(x) ) RETURN
    y = x
  END SUBROUTINE assign_if_present_real_sp

  SUBROUTINE assign_if_present_logical_allocatable_1d(y, x)
    LOGICAL, ALLOCATABLE, INTENT(INOUT) :: y(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: x(:)

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_logical_allocatable_1d"

    IF(.NOT.PRESENT(x)) RETURN
    IF(ALLOCATED(y)) THEN
        IF(SIZE(y) /= SIZE(x)) DEALLOCATE(y)
    END IF
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y(SIZE(x)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y(:) = x(:)
  END SUBROUTINE assign_if_present_logical_allocatable_1d

  SUBROUTINE assign_if_present_integer_allocatable(y, x)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: y
    INTEGER, OPTIONAL, INTENT(IN) :: x

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_integer_allocatable"

    IF(.NOT.PRESENT(x)) RETURN
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y = x
  END SUBROUTINE assign_if_present_integer_allocatable

  SUBROUTINE assign_if_present_integer_allocatable_1d(y, x)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: y(:)
    INTEGER, OPTIONAL, INTENT(IN) :: x(:)

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_integer_allocatable_1d"

    IF(.NOT.PRESENT(x)) RETURN
    IF(ALLOCATED(y)) THEN
        IF(SIZE(y) /= SIZE(x)) DEALLOCATE(y)
    END IF
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y(SIZE(x)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y(:) = x(:)
  END SUBROUTINE assign_if_present_integer_allocatable_1d

  SUBROUTINE assign_if_present_real_allocatable(y, x)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: y
    REAL(wp), OPTIONAL, INTENT(IN) :: x

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_real_allocatable"

    IF(.NOT.PRESENT(x)) RETURN
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y = x
  END SUBROUTINE assign_if_present_real_allocatable

  SUBROUTINE assign_if_present_real_allocatable_1d(y, x)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: y(:)
    REAL(wp), OPTIONAL, INTENT(IN) :: x(:)

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_real_allocatable_1d"

    IF(.NOT.PRESENT(x)) RETURN
    IF(ALLOCATED(y)) THEN
        IF(SIZE(y) /= SIZE(x)) DEALLOCATE(y)
    END IF
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y(SIZE(x)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y(:) = x(:)
  END SUBROUTINE assign_if_present_real_allocatable_1d

  SUBROUTINE alloc_int_1d(array, allocSize)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: array(:)
    INTEGER, VALUE :: allocSize

    INTEGER :: error
    CHARACTER(*), PARAMETER :: routine = modname//":alloc_int_1d"

    IF(allocSize < 1) allocSize = 1
    IF(ALLOCATED(array)) THEN
        IF(SIZE(array) == allocSize) RETURN
        DEALLOCATE(array)
    END IF
    ALLOCATE(array(allocSize), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
  END SUBROUTINE alloc_int_1d

  SUBROUTINE alloc_double_1d(array, allocSize)
    REAL(dp), ALLOCATABLE, INTENT(INOUT) :: array(:)
    INTEGER, VALUE :: allocSize

    INTEGER :: error
    CHARACTER(*), PARAMETER :: routine = modname//":alloc_double_1d"

    IF(allocSize < 1) allocSize = 1
    IF(ALLOCATED(array)) THEN
        IF(SIZE(array) == allocSize) RETURN
        DEALLOCATE(array)
    END IF
    ALLOCATE(array(allocSize), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
  END SUBROUTINE alloc_double_1d

  SUBROUTINE alloc_single_1d(array, allocSize)
    REAL(sp), ALLOCATABLE, INTENT(INOUT) :: array(:)
    INTEGER, VALUE :: allocSize

    INTEGER :: error
    CHARACTER(*), PARAMETER :: routine = modname//":alloc_single_1d"

    IF(allocSize < 1) allocSize = 1
    IF(ALLOCATED(array)) THEN
        IF(SIZE(array) == allocSize) RETURN
        DEALLOCATE(array)
    END IF
    ALLOCATE(array(allocSize), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
  END SUBROUTINE alloc_single_1d

  SUBROUTINE ensureSize_dp_1d(buffer, requiredSize)
    REAL(wp), POINTER, INTENT(INOUT) :: buffer(:)
    INTEGER, VALUE ::requiredSize

    REAL(wp), POINTER :: newBuffer(:)
    INTEGER :: oldSize, error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":ensureSize_dp_1d"

    IF(ASSOCIATED(buffer)) THEN
        oldSize = SIZE(buffer, 1)
        IF(oldSize >= requiredSize) RETURN  ! nothing to DO IF it's already big enough
        requiredSize = MAX(requiredSize, 2*oldSize) ! avoid quadratic complexity

        ALLOCATE(newBuffer(requiredSize), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")

        newBuffer(1:oldSize) = buffer(1:oldSize)
        newBuffer(oldSize + 1:requiredSize) = 0.0

        DEALLOCATE(buffer)
        buffer => newBuffer
        newBuffer => NULL()
    ELSE
        ALLOCATE(buffer(requiredSize), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")

        buffer(1:requiredSize) = 0.0
    END IF
  END SUBROUTINE ensureSize_dp_1d


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
#ifdef _OPENACC
!$ACC DATA PCOPYIN( src ), PCOPYOUT( dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( src ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC PARALLEL PRESENT( src, dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(2)
#else
!$omp do collapse(2)
#endif
    DO i2 = 1, m2
      DO i1 = 1, m1
        dest(i1, i2) = src(i1, i2)
      END DO
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( dest ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif

  END SUBROUTINE copy_2d_dp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_3d_dp(src, dest)
    REAL(dp), INTENT(in) :: src(:, :, :)
    REAL(dp), INTENT(out) :: dest(:, :, :)
    INTEGER :: i1, i2, i3, m1, m2, m3
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
#ifdef _OPENACC
!$ACC DATA PCOPYIN( src ), PCOPYOUT( dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( src ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC PARALLEL PRESENT( src, dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(3)
#else
#ifdef _CRAYFTN
!$omp do
#else
!$omp do collapse(3)
#endif
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          dest(i1, i2, i3) = src(i1, i2, i3)
        END DO
      END DO
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( dest ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
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
#ifdef _OPENACC
!$ACC DATA PCOPYIN( src ), PCOPYOUT( dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( src ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC PARALLEL PRESENT( src, dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(4)
#else
#ifdef _CRAYFTN
!$omp do
#else
!$omp do collapse(4)
#endif
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
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( dest ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
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
#ifdef _OPENACC
!$ACC DATA PCOPYIN( src ), PCOPYOUT( dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( src ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC PARALLEL PRESENT( src, dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(5)
#else
!$omp do collapse(5)
#endif
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
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( dest ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
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
#ifdef _OPENACC
!$ACC DATA PCOPYIN( src ), PCOPYOUT( dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( src ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC PARALLEL PRESENT( src, dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(5)
#else
!$omp do collapse(5)
#endif
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
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( dest ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
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
#ifdef _OPENACC
!$ACC DATA PCOPYIN( src ), PCOPYOUT( dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( src ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC PARALLEL PRESENT( src, dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(5)
#else
!$omp do collapse(5)
#endif
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
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( dest ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE copy_5d_spdp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_2d_i4(src, dest)
    INTEGER(ik4), INTENT(in) :: src(:, :)
    INTEGER(ik4), INTENT(out) :: dest(:, :)
    INTEGER :: i1, i2, m1, m2
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
#ifdef _OPENACC
!$ACC DATA PCOPYIN( src ), PCOPYOUT( dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( src ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC PARALLEL PRESENT( src, dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(2)
#else
!$omp do collapse(2)
#endif
    DO i2 = 1, m2
      DO i1 = 1, m1
        dest(i1, i2) = src(i1, i2)
      END DO
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( dest ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
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

#ifdef _OPENACC
!$ACC DATA PCOPYIN( src ), PCOPYOUT( dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( src ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC PARALLEL PRESENT( src, dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(5)
#else
!$omp do collapse(5)
#endif
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
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( dest ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE copy_5d_i4

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_3d_i4(src, dest)
    INTEGER(ik4), INTENT(in) :: src(:, :, :)
    INTEGER(ik4), INTENT(out) :: dest(:, :, :)
    INTEGER :: i1, i2, i3, m1, m2, m3
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
#ifdef _OPENACC
!$ACC DATA PCOPYIN( src ), PCOPYOUT( dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( src ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC PARALLEL PRESENT( src, dest ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(3)
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
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( dest ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE copy_3d_i4


  SUBROUTINE init_zero_3d_dp(init_var)
    REAL(dp), INTENT(out) :: init_var(:, :, :)
    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
#ifdef _OPENACC
!$ACC DATA PCOPYOUT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(3)
#else
#ifdef _CRAYFTN
!$omp do
#else
!$omp do collapse(3)
#endif
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = 0.0_dp
        END DO
      END DO
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( init_var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE init_zero_3d_dp

  SUBROUTINE init_zero_4d_dp(init_var)
    REAL(dp), INTENT(out) :: init_var(:, :, :, :)
    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)
#ifdef _OPENACC
!$ACC DATA PCOPYOUT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(4)
#else
#ifdef _CRAYFTN
!$omp do
#else
!$omp do collapse(4)
#endif
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
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( init_var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE init_zero_4d_dp

  SUBROUTINE init_zero_4d_sp(init_var)
    REAL(sp), INTENT(out) :: init_var(:, :, :, :)
    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)
#ifdef _OPENACC
!$ACC DATA PCOPYOUT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(4)
#else
#ifdef _CRAYFTN
!$omp do
#else
!$omp do collapse(4)
#endif
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
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( init_var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE init_zero_4d_sp

  SUBROUTINE init_zero_3d_sp(init_var)
    REAL(sp), INTENT(out) :: init_var(:, :, :)
    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
#ifdef _OPENACC
!$ACC DATA PCOPYOUT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = 0.0_sp
        END DO
      END DO
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( init_var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif

  END SUBROUTINE init_zero_3d_sp

  SUBROUTINE init_zero_2d_i4(init_var)
    INTEGER(ik4), INTENT(out) :: init_var(:, :)
    INTEGER :: i1, i2, m1, m2

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)

#ifdef _OPENACC
!$ACC DATA PCOPYOUT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(2)
#else
!$omp do collapse(2)
#endif
    DO i2 = 1, m2
      DO i1 = 1, m1
        init_var(i1, i2) = 0_ik4
      END DO
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( init_var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE init_zero_2d_i4

  SUBROUTINE init_zero_3d_i4(init_var)
    INTEGER(ik4), INTENT(out) :: init_var(:, :, :)
    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
#ifdef _OPENACC
!$ACC DATA PCOPYOUT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = 0_ik4
        END DO
      END DO
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( init_var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE init_zero_3d_i4

  SUBROUTINE init_zero_4d_i4(init_var)
    INTEGER(ik4), INTENT(out) :: init_var(:, :, :, :)
    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)
#ifdef _OPENACC
!$ACC DATA PCOPY( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(4)
#else
!$omp do collapse(4)
#endif
    DO i4 = 1, m4
      DO i3 = 1, m3
        DO i2 = 1, m2
          DO i1 = 1, m1
            init_var(i1, i2, i3, i4) = 0_ik4
          END DO
        END DO
      END DO
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( init_var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE init_zero_4d_i4

  SUBROUTINE init_zero_1d_dp(init_var)
    REAL(dp), INTENT(out) :: init_var(:)
    INTEGER :: i1, m1

    m1 = SIZE(init_var, 1)
#ifdef _OPENACC
!$ACC DATA PCOPYOUT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP 
#else
!$omp do
#endif
    DO i1 = 1, m1
      init_var(i1) = 0.0_dp
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( init_var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE init_zero_1d_dp

  SUBROUTINE init_zero_2d_dp(init_var)
    REAL(dp), INTENT(out) :: init_var(:, :)
    INTEGER :: i1, i2, m1, m2

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
#ifdef _OPENACC
!$ACC DATA PCOPYOUT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(2)
#else
!$omp do collapse(2)
#endif
    DO i2 = 1, m2
      DO i1 = 1, m1
        init_var(i1, i2) = 0.0_dp
      END DO
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( init_var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE init_zero_2d_dp

  SUBROUTINE init_3d_dp(init_var, init_val)
    REAL(dp), INTENT(out) :: init_var(:, :, :)
    REAL(dp), INTENT(in) :: init_val

    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
#ifdef _OPENACC
!$ACC DATA PCOPYOUT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = init_val
        END DO
      END DO
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( init_var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE init_3d_dp

  SUBROUTINE init_3d_spdp(init_var, init_val)
    REAL(sp), INTENT(out) :: init_var(:, :, :)
    REAL(dp), INTENT(in) :: init_val

    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
#ifdef _OPENACC
!$ACC DATA PCOPYOUT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = REAL(init_val,KIND=sp)
        END DO
      END DO
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( init_var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
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
#ifdef _OPENACC
!$ACC DATA PCOPYOUT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(5)
#else
!$omp do collapse(5)
#endif
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
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( init_var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
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
#ifdef _OPENACC
!$ACC DATA PCOPYOUT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( init_var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(5)
#else
!$omp do collapse(5)
#endif
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
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( init_var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE init_5d_i4


  SUBROUTINE var_scale_3d_dp(var, scale_val)
    REAL(dp), INTENT(inout) :: var(:, :, :)
    REAL(dp), INTENT(in) :: scale_val

    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(var, 1)
    m2 = SIZE(var, 2)
    m3 = SIZE(var, 3)

#ifdef _OPENACC
!$ACC DATA PCOPYOUT( var ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC PARALLEL PRESENT( var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          var(i1, i2, i3) = var(i1, i2, i3) * scale_val
        END DO
      END DO
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE var_scale_3d_dp


  ! add a constant value to a 3D field
  SUBROUTINE var_addc_3d_dp(var, add_val)
    REAL(dp), INTENT(inout) :: var(:, :, :)
    REAL(dp), INTENT(in) :: add_val

    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(var, 1)
    m2 = SIZE(var, 2)
    m3 = SIZE(var, 3)
!$omp do collapse(3)
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          var(i1, i2, i3) = var(i1, i2, i3) + add_val
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE var_addc_3d_dp


  SUBROUTINE negative2zero_4d_dp(var)
    REAL(dp), INTENT(inout) :: var(:, :, :, :)

    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4
    REAL(dp) :: v

    m1 = SIZE(var, 1)
    m2 = SIZE(var, 2)
    m3 = SIZE(var, 3)
    m4 = SIZE(var, 4)

#ifdef _OPENACC
!$ACC DATA PCOPY( var ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC PARALLEL PRESENT( var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP COLLAPSE(4)
#else
!$omp do collapse(4)
#endif
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
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE negative2zero_4d_dp

  SUBROUTINE init_zero_contiguous_dp(var, n)
    INTEGER, INTENT(in) :: n
    REAL(dp), INTENT(out) :: var(n)

    INTEGER :: i
#ifdef _OPENACC
!$ACC DATA PCOPYOUT( var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP
#else
!$omp do
#endif
    DO i = 1, n
      var(i) = 0.0_dp
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE init_zero_contiguous_dp

  SUBROUTINE init_zero_contiguous_sp(var, n)
    INTEGER, INTENT(in) :: n
    REAL(sp), INTENT(out) :: var(n)

    INTEGER :: i
#ifdef _OPENACC
!$ACC DATA PCOPYOUT( var ), IF( i_am_accel_node .AND. acc_on )
!$ACC PARALLEL PRESENT( var ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP
#else
!$omp do
#endif
    DO i = 1, n
      var(i) = 0.0_sp
    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
!$ACC UPDATE HOST( var ), IF( i_am_accel_node .AND. acc_on .AND. acc_validate )
!$ACC END DATA
#else
!$omp end do nowait
#endif
  END SUBROUTINE init_zero_contiguous_sp

  SUBROUTINE insert_dimension_r_wp_6_5_s(ptr_out, ptr_in, in_shape, new_dim_rank)
    INTEGER, INTENT(in) :: in_shape(5), new_dim_rank
    REAL(wp), POINTER, INTENT(out) :: ptr_out(:,:,:,:,:,:)
    REAL(wp), TARGET, INTENT(in) :: ptr_in(in_shape(1),in_shape(2),&
         in_shape(3),in_shape(4),in_shape(5))
    INTEGER :: out_shape(6), i
    TYPE(c_ptr) :: cptr
    out_shape(1:5) = SHAPE(ptr_in)
    cptr = C_LOC(ptr_in)
    DO i = 6, new_dim_rank, -1
      out_shape(i) = out_shape(i-1)
    END DO
    out_shape(new_dim_rank) = 1
    CALL C_F_POINTER(cptr, ptr_out, out_shape)
  END SUBROUTINE insert_dimension_r_wp_6_5_s

  ! insert dimension of size 1 (so that total array size remains the
  ! same but an extra dimension is inserted into the shape)
  SUBROUTINE insert_dimension_r_wp_6_5(ptr_out, ptr_in, new_dim_rank)
    REAL(wp), POINTER, INTENT(out) :: ptr_out(:,:,:,:,:,:)
    REAL(wp), TARGET, INTENT(in) :: ptr_in(:,:,:,:,:)
    INTEGER, INTENT(in) :: new_dim_rank
    INTEGER :: in_shape(5)
    in_shape = SHAPE(ptr_in)
    CALL insert_dimension(ptr_out, ptr_in, in_shape, new_dim_rank)
  END SUBROUTINE insert_dimension_r_wp_6_5

  SUBROUTINE insert_dimension_r_sp_6_5_s(ptr_out, ptr_in, in_shape, new_dim_rank)
    INTEGER, INTENT(in) :: in_shape(5), new_dim_rank
    REAL(sp), POINTER, INTENT(out) :: ptr_out(:,:,:,:,:,:)
    REAL(sp), TARGET, INTENT(in) :: ptr_in(in_shape(1),in_shape(2),&
         in_shape(3),in_shape(4),in_shape(5))
    INTEGER :: out_shape(6), i
    TYPE(c_ptr) :: cptr
    out_shape(1:5) = SHAPE(ptr_in)
    cptr = C_LOC(ptr_in)
    DO i = 6, new_dim_rank, -1
      out_shape(i) = out_shape(i-1)
    END DO
    out_shape(new_dim_rank) = 1
    CALL C_F_POINTER(cptr, ptr_out, out_shape)
  END SUBROUTINE insert_dimension_r_sp_6_5_s

  ! insert dimension of size 1 (so that total array size remains the
  ! same but an extra dimension is inserted into the shape)
  SUBROUTINE insert_dimension_r_sp_6_5(ptr_out, ptr_in, new_dim_rank)
    REAL(sp), POINTER, INTENT(out) :: ptr_out(:,:,:,:,:,:)
    REAL(sp), TARGET, INTENT(in) :: ptr_in(:,:,:,:,:)
    INTEGER, INTENT(in) :: new_dim_rank
    INTEGER :: in_shape(5)
    in_shape = SHAPE(ptr_in)
    CALL insert_dimension(ptr_out, ptr_in, in_shape, new_dim_rank)
  END SUBROUTINE insert_dimension_r_sp_6_5

  SUBROUTINE insert_dimension_i4_6_5_s(ptr_out, ptr_in, in_shape, new_dim_rank)
    INTEGER, INTENT(in) :: in_shape(5), new_dim_rank
    INTEGER(ik4), POINTER, INTENT(out) :: ptr_out(:,:,:,:,:,:)
    INTEGER(ik4), TARGET, INTENT(in) :: ptr_in(in_shape(1),in_shape(2),&
         in_shape(3),in_shape(4),in_shape(5))
    INTEGER :: out_shape(6), i
    TYPE(c_ptr) :: cptr
    out_shape(1:5) = SHAPE(ptr_in)
    cptr = C_LOC(ptr_in)
    DO i = 6, new_dim_rank, -1
      out_shape(i) = out_shape(i-1)
    END DO
    out_shape(new_dim_rank) = 1
    CALL C_F_POINTER(cptr, ptr_out, out_shape)
  END SUBROUTINE insert_dimension_i4_6_5_s

  ! insert dimension of size 1 (so that total array size remains the
  ! same but an extra dimension is inserted into the shape)
  SUBROUTINE insert_dimension_i4_6_5(ptr_out, ptr_in, new_dim_rank)
    INTEGER(ik4), POINTER, INTENT(out) :: ptr_out(:,:,:,:,:,:)
    INTEGER(ik4), TARGET, INTENT(in) :: ptr_in(:,:,:,:,:)
    INTEGER, INTENT(in) :: new_dim_rank
    INTEGER :: in_shape(5)
    in_shape = SHAPE(ptr_in)
    CALL insert_dimension(ptr_out, ptr_in, in_shape, new_dim_rank)
  END SUBROUTINE insert_dimension_i4_6_5


  ! AUXILIARY ROUTINES FOR DEALLOCATION

  SUBROUTINE DO_DEALLOCATE_r4D(object)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: object(:,:,:,:)
    INTEGER :: ierrstat
    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_DEALLOCATE_r4D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_DEALLOCATE_R4D

  SUBROUTINE DO_DEALLOCATE_r3D(object)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: object(:,:,:)
    INTEGER :: ierrstat
    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_DEALLOCATE_r3D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_DEALLOCATE_R3D

  SUBROUTINE DO_DEALLOCATE_r2D(object)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: object(:,:)
    INTEGER :: ierrstat
    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_DEALLOCATE_r2D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_DEALLOCATE_R2D

  SUBROUTINE DO_DEALLOCATE_r1D(object)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: object(:)
    INTEGER :: ierrstat
    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_DEALLOCATE_r1D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_DEALLOCATE_R1D

  SUBROUTINE DO_DEALLOCATE_i3D(object)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: object(:,:,:)
    INTEGER :: ierrstat
    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_DEALLOCATE_i3D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_DEALLOCATE_i3D

  SUBROUTINE DO_DEALLOCATE_i2D(object)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: object(:,:)
    INTEGER :: ierrstat
    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_DEALLOCATE_i2D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_DEALLOCATE_i2D

  SUBROUTINE DO_PTR_DEALLOCATE_r3D(object)
    REAL(wp), POINTER, INTENT(INOUT) :: object(:,:,:)
    INTEGER :: ierrstat
    IF (ASSOCIATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_PTR_DEALLOCATE_r3D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_PTR_DEALLOCATE_R3D

  SUBROUTINE DO_PTR_DEALLOCATE_r2D(object)
    REAL(wp), POINTER, INTENT(INOUT) :: object(:,:)
    INTEGER :: ierrstat
    IF (ASSOCIATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_PTR_DEALLOCATE_r2D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_PTR_DEALLOCATE_R2D



END MODULE mo_fortran_tools
