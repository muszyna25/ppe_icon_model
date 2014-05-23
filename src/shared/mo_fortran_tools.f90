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
  PUBLIC :: t_ptr_tracer!,pcen,ptenc


  INTERFACE assign_if_present
    MODULE PROCEDURE assign_if_present_character
    MODULE PROCEDURE assign_if_present_logical
    MODULE PROCEDURE assign_if_present_logicals
    MODULE PROCEDURE assign_if_present_integer
    MODULE PROCEDURE assign_if_present_integers
    MODULE PROCEDURE assign_if_present_real
  END INTERFACE assign_if_present

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

END MODULE mo_fortran_tools
