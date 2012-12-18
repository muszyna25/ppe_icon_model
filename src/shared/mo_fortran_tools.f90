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
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!
MODULE mo_fortran_tools

  USE mo_kind,                    ONLY: wp
  IMPLICIT NONE
  PRIVATE


  TYPE t_ptr_2d3d
    REAL(wp),POINTER :: p_3d(:,:,:)  ! pointer to 3D (spatial) array
    REAL(wp),POINTER :: p_2d(:,:)    ! pointer to 2D (spatial) array
  END TYPE t_ptr_2d3d
  
  TYPE t_ptr_tracer                   ! Type to pass pointer arrays to convection subroutines
    REAL(wp), POINTER :: ptr(:,:)
    INTEGER           :: idx_tracer    
  END TYPE t_ptr_tracer

  PUBLIC :: assign_if_present
  PUBLIC :: t_ptr_2d3d
  PUBLIC :: t_ptr_tracer!,pcen,ptenc


  INTERFACE assign_if_present  ! purely internal
    MODULE PROCEDURE assign_if_present_character
    MODULE PROCEDURE assign_if_present_logical
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
