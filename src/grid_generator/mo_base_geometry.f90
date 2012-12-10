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
MODULE mo_base_geometry

  USE mo_kind, ONLY: wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: x_rot_angle,y_rot_angle,z_rot_angle

  REAL(wp) :: x_rot_angle,y_rot_angle,z_rot_angle

END MODULE mo_base_geometry
!----------------------------------------------------------------------------









