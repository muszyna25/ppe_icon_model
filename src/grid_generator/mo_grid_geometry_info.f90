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
  USE mo_physical_constants, ONLY: earth_radius
  USE mo_math_constants,     ONLY: pi
  USE mo_exception,          ONLY: message_text, message, finish, warning

#ifndef NOMPI
  ! The USE statement below lets this module use the routines from
  ! mo_read_netcdf_parallel where only 1 processor is reading
  ! and broadcasting the results  
  USE mo_read_netcdf_parallel, ONLY:                &
    & nf_nowrite, nf_global, nf_noerr, nf_strerror,  &
    & nf_open            => p_nf_open,               &
    & nf_close           => p_nf_close,              &
    & nf_inq_dimid       => p_nf_inq_dimid,          &
    & nf_inq_dimlen      => p_nf_inq_dimlen,         &
    & nf_inq_varid       => p_nf_inq_varid,          &
    & nf_get_att_text    => p_nf_get_att_text,       &
    & nf_get_att_int     => p_nf_get_att_int,        &
    & nf_get_att_double  => p_nf_get_att_double,     &
    & nf_get_var_int     => p_nf_get_var_int,        &
    & nf_get_var_double  => p_nf_get_var_double
#endif
  

  IMPLICIT NONE

  PRIVATE
#ifdef NOMPI
  INCLUDE 'netcdf.inc'
#endif

  ! public parameters
  PUBLIC :: sphere_geometry, planar_torus_geometry
  PUBLIC :: triangular_cell, hexagonal_cell
  PUBLIC :: cut_off_grid, refined_bisection_grid, dualy_refined_grid

  ! public structures  
  PUBLIC :: t_grid_geometry_info

  ! public methods
  PUBLIC :: set_default_geometry_info, copy_grid_geometry_info, &
    & set_grid_geometry_derived_info, read_geometry_info

  ! -----------------------------
  ! types of grid geometries
  INTEGER, PARAMETER ::  sphere_geometry       = 1
  INTEGER, PARAMETER ::  planar_torus_geometry = 2
  
  ! -----------------------------
  ! types of grids
  INTEGER, PARAMETER ::  triangular_cell = 3
  INTEGER, PARAMETER ::  hexagonal_cell  = 6

  ! -----------------------------
  ! types of grid creation
  INTEGER, PARAMETER ::  cut_off_grid = 1
  INTEGER, PARAMETER ::  refined_bisection_grid = 2
  INTEGER, PARAMETER ::  dualy_refined_grid = 3
  
  !--------------------------------------------------------------
  INTEGER, PARAMETER ::  undefined    = -1
  !--------------------------------------------------------------
  !> Holds the grid geometry parameters
  TYPE t_grid_geometry_info
    INTEGER :: cell_type          ! triangular_cell, etc
    INTEGER :: geometry_type      ! sphere_geometry, etc
    !> The creation process of the grid (cut_off, refined, etc, see parameters).
    INTEGER :: grid_creation_process 
    !> The grid optimization procces
    INTEGER :: grid_optimization_process
        
    TYPE(t_cartesian_coordinates) :: center
    REAL(wp) :: mean_edge_length  ! (meters)
    REAL(wp) :: mean_cell_area    ! (meters^2)
    REAL(wp) :: domain_length     ! (meters)
    REAL(wp) :: domain_height     ! (meters)

    !> Sphere parameters 
    REAL(wp) :: sphere_radius
    
    !> derived info, the following can be calclulated from the previous
    REAL(wp) :: mean_characteristic_length ! the sqrt(mean_cell_area)
    
  END TYPE t_grid_geometry_info
  
CONTAINS

  !------------------------------------------------------------------------
  !>
  ! The default is sphere triangular geometry
  ! Note the cell characteristice (mean area and lenght) are unknown
  ! and set to 0
  SUBROUTINE set_default_geometry_info(to_geometry_info)
    TYPE(t_grid_geometry_info) :: to_geometry_info

    to_geometry_info%cell_type                  = triangular_cell
    to_geometry_info%geometry_type              = sphere_geometry
    to_geometry_info%grid_creation_process      = undefined
    to_geometry_info%grid_optimization_process  = undefined
    to_geometry_info%center %x(:)               = 0.0_wp
    to_geometry_info%mean_edge_length           = 0.0_wp
    to_geometry_info%mean_cell_area             = 0.0_wp
    to_geometry_info%domain_length              = 2.0_wp * pi * earth_radius
    to_geometry_info%domain_height              = 2.0_wp * pi * earth_radius
    to_geometry_info%sphere_radius              = earth_radius
    to_geometry_info%mean_characteristic_length = 0.0_wp

  END SUBROUTINE set_default_geometry_info
  !------------------------------------------------------------------------
    
  !------------------------------------------------------------------------
  !>
  SUBROUTINE copy_grid_geometry_info(from_geometry_info, to_geometry_info)
    TYPE(t_grid_geometry_info) :: from_geometry_info, to_geometry_info

    to_geometry_info%cell_type                  = from_geometry_info%cell_type
    to_geometry_info%geometry_type              = from_geometry_info%geometry_type
    to_geometry_info%grid_creation_process      = from_geometry_info%grid_creation_process
    to_geometry_info%grid_optimization_process  = from_geometry_info%grid_optimization_process
    to_geometry_info%center                     = from_geometry_info%center
    to_geometry_info%mean_edge_length           = from_geometry_info%mean_edge_length
    to_geometry_info%mean_cell_area             = from_geometry_info%mean_cell_area
    to_geometry_info%domain_length              = from_geometry_info%domain_length
    to_geometry_info%domain_height              = from_geometry_info%domain_height
    to_geometry_info%sphere_radius              = from_geometry_info%sphere_radius
    to_geometry_info%mean_characteristic_length = from_geometry_info%mean_characteristic_length

  END SUBROUTINE copy_grid_geometry_info
  !------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !> 
  SUBROUTINE set_grid_geometry_derived_info( to_geometry_info )
    TYPE(t_grid_geometry_info) :: to_geometry_info

    ! derived geometry parameters
    to_geometry_info%mean_characteristic_length = SQRT(to_geometry_info%mean_cell_area)
!     write(0,*) "------------------------------------------------"
!     write(0,*) "mean_cell_area:", to_geometry_info%mean_cell_area
!     write(0,*) "mean_characteristic_length:", to_geometry_info%mean_characteristic_length
!     write(0,*) "------------------------------------------------"
    
  END SUBROUTINE set_grid_geometry_derived_info
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  INTEGER FUNCTION read_geometry_info(ncid, geometry_info)
    INTEGER, INTENT(in) :: ncid
    TYPE(t_grid_geometry_info) :: geometry_info
    
    INTEGER :: netcd_status
    CHARACTER(*), PARAMETER :: method_name = "read_geometry_info"
            
    CALL set_default_geometry_info(geometry_info)
    read_geometry_info = -1
    
    netcd_status = nf_get_att_int(ncid, nf_global,'grid_geometry', geometry_info%geometry_type)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","grid_geometry")
      RETURN
    ENDIF
        
    netcd_status = nf_get_att_double(ncid, nf_global,'mean_edge_length', &
      & geometry_info%mean_edge_length)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","mean_edge_length")
      RETURN
    ENDIF
    
    netcd_status = nf_get_att_double(ncid, nf_global,'mean_cell_area', &
      & geometry_info%mean_cell_area)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","mean_cell_area")
      RETURN
    ENDIF
    
    netcd_status = nf_get_att_double(ncid, nf_global,'domain_length', &
      & geometry_info%domain_length)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","domain_length")
      RETURN
    ENDIF
    
    netcd_status = nf_get_att_double(ncid, nf_global,'domain_height', &
      & geometry_info%domain_height)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","domain_height")
      RETURN
    ENDIF

    netcd_status = nf_get_att_double(ncid, nf_global,'sphere_radius', geometry_info%sphere_radius)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","sphere_radius")
      RETURN
    ENDIF
    
    netcd_status = nf_get_att_double(ncid, nf_global,'domain_cartesian_center', &
      & geometry_info%center%x )
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","domain_cartesian_center")
      RETURN
    ENDIF
    
    ! return status ok
    read_geometry_info = 0
    
  END FUNCTION read_geometry_info
  !-------------------------------------------------------------------------

END MODULE mo_grid_geometry_info
!----------------------------------------------------------------------------









