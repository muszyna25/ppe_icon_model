!>
!!  Contains the definition of basic structures and geometry parameters 
!!  These are included in the grid/patch info
!!
!! @par Revision History
!!  Initial version  by Leonidas Linardakis, MPIM (2012-12)
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
MODULE mo_grid_geometry_info

  USE mo_kind, ONLY: wp
  USE mo_math_utilities,  ONLY: t_geographical_coordinates, t_cartesian_coordinates, &
    & t_tangent_vectors

  IMPLICIT NONE

  PRIVATE

  ! public parameters
  PUBLIC :: sphere_geometry, planar_torus_geometry

  ! public structures
  PUBLIC :: t_planar_torus_geometry_info

  ! public methods
  PUBLIC :: copy_planar_torus_info

  ! -----------------------------
  ! types of grid geometries
  INTEGER, PARAMETER ::  sphere_geometry       = 1
  INTEGER, PARAMETER ::  planar_torus_geometry = 2
  
  ! -----------------------------
  ! types of grids
  INTEGER, PARAMETER ::  triangular_grid = 3
  INTEGER, PARAMETER ::  hexagonal_grid  = 6

  !--------------------------------------------------------------
  !> Holds the planar torus gemoetry parameters
  TYPE t_planar_torus_geometry_info
    REAL(wp) :: cell_edge_length
    TYPE(t_cartesian_coordinates) :: center
    REAL(wp) :: length
    REAL(wp) :: height
  END TYPE t_planar_torus_geometry_info
  
CONTAINS

  SUBROUTINE copy_planar_torus_info(from_planar_torus_info, to_planar_torus_info)
    TYPE(t_planar_torus_geometry_info) :: from_planar_torus_info, to_planar_torus_info

    to_planar_torus_info%cell_edge_length = from_planar_torus_info%cell_edge_length
    to_planar_torus_info%center%x         = from_planar_torus_info%center%x
    to_planar_torus_info%length           = from_planar_torus_info%length
    to_planar_torus_info%height           = from_planar_torus_info%height
  
  END SUBROUTINE copy_planar_torus_info
    

END MODULE mo_grid_geometry_info
!----------------------------------------------------------------------------









