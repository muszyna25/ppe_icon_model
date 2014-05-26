!>
!!   Contains basic math types
!!
!! @par Revision History
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
MODULE mo_math_types
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: t_cartesian_coordinates
  PUBLIC :: t_geographical_coordinates
  PUBLIC :: t_line
  PUBLIC :: t_tangent_vectors
    
  ! cartesian coordinate class
  TYPE t_cartesian_coordinates
    REAL(wp) :: x(3)
  END TYPE t_cartesian_coordinates
  
  ! geographical coordinate class
  TYPE t_geographical_coordinates
    REAL(wp) :: lon
    REAL(wp) :: lat
  END TYPE t_geographical_coordinates

  ! the two coordinates on the tangent plane
  TYPE t_tangent_vectors
    REAL(wp) :: v1
    REAL(wp) :: v2
  END TYPE t_tangent_vectors
  
  ! line class
  TYPE t_line
    TYPE(t_geographical_coordinates) :: p1
    TYPE(t_geographical_coordinates) :: p2
  END TYPE t_line
  
    
END MODULE mo_math_types

