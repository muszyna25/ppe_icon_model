!>
!!  Contains the definition of basic used defined data types.
!!
!!  Contains the definition of basic used defined data types
!!  to store the gridpoint coordinates either in cartesian or
!!  spherical coordinates. It also contains a number of (mostly
!!  ELEMENTAL) functions used to compute vectors and geometric quantities
!!  needed by the grid generator.
!!
!! @par Revision History
!!  Initial version  by Luis Kornblueh (2004)
!!  Modified to include tangent vectors and Protex headers
!!  by Luca Bonaventura (2005)
!! @par
!!  Guenther Zaengl, DWD, 2008-10-10:
!!  Add subroutine gvec2cvec (copied from mo_math_utilities)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_base_geometry

  USE mo_kind, ONLY: wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: x_rot_angle,y_rot_angle,z_rot_angle

  REAL(wp) :: x_rot_angle,y_rot_angle,z_rot_angle

END MODULE mo_base_geometry
!----------------------------------------------------------------------------









