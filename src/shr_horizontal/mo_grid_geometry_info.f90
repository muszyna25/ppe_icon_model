!>
!!  Contains the definition of basic structures and geometry parameters 
!!  These are included in the grid/patch info
!!
!! @par Revision History
!!  Initial version  by Leonidas Linardakis, MPIM (2012-12)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_grid_geometry_info

  USE mo_kind, ONLY: wp
  USE mo_math_types,  ONLY: t_cartesian_coordinates
  USE mo_physical_constants, ONLY: earth_radius
  USE mo_math_constants,     ONLY: pi
  USE mo_exception,          ONLY: finish

#if !( defined (NOMPI) || defined (__ICON_GRID_GENERATOR__))
  ! The USE statement below lets this module use the routines from
  ! mo_netcdf_parallel where only 1 processor is reading and
  ! broadcasting the results  
  USE mo_netcdf_parallel, ONLY:                        &
!     & nf_nowrite, nf_global, nf_noerr, nf_strerror,  &
!     & nf_open            => p_nf_open,               &
!     & nf_close           => p_nf_close,              &
!     & nf_inq_dimid       => p_nf_inq_dimid,          &
!     & nf_inq_dimlen      => p_nf_inq_dimlen,         &
!     & nf_inq_varid       => p_nf_inq_varid,          &
!     & nf_get_att_text    => p_nf_get_att_text,       &
!     & p_nf_get_att_int     => p_nf_get_att_int,        &
!     & p_nf_get_att_double  => p_nf_get_att_double,     &
!     & nf_get_var_int     => p_nf_get_var_int,        &
!     & nf_get_var_double  => p_nf_get_var_double
    & p_nf_get_att_int, p_nf_get_att_double
#endif
  

  IMPLICIT NONE

  PRIVATE
  INCLUDE 'netcdf.inc'


  ! public parameters
  PUBLIC :: sphere_geometry, planar_torus_geometry
  PUBLIC :: triangular_cell, hexagonal_cell
  PUBLIC :: cut_off_grid, refined_bisection_grid, dualy_refined_grid

  ! public structures  
  PUBLIC :: t_grid_geometry_info

  ! public methods
  PUBLIC :: set_default_geometry_info, copy_grid_geometry_info, &
    & set_grid_geometry_derived_info, read_geometry_info,       &
    & parallel_read_geometry_info, write_geometry_info,         &
    & get_resolution_string

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
    REAL(wp) :: mean_edge_length       ! (meters)
    REAL(wp) :: mean_dual_edge_length  ! (meters)
    REAL(wp) :: mean_cell_area         ! (meters^2)
    REAL(wp) :: mean_dual_cell_area    ! (meters^2)
    REAL(wp) :: domain_length          ! (meters)
    REAL(wp) :: domain_height          ! (meters)

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
    to_geometry_info%mean_dual_edge_length      = 0.0_wp
    to_geometry_info%mean_cell_area             = 0.0_wp
    to_geometry_info%mean_dual_cell_area        = 0.0_wp
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
    to_geometry_info%mean_dual_edge_length      = from_geometry_info%mean_dual_edge_length
    to_geometry_info%mean_cell_area             = from_geometry_info%mean_cell_area
    to_geometry_info%mean_dual_cell_area        = from_geometry_info%mean_dual_cell_area
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
  FUNCTION get_resolution_string(geometry_info) result(resolution_string)
    TYPE(t_grid_geometry_info) :: geometry_info

    CHARACTER(len=16) :: resolution_string

    CALL set_grid_geometry_derived_info(geometry_info)
    IF (geometry_info%mean_characteristic_length > 5000.0_wp) THEN
      WRITE(resolution_string,'(i4.4,a)') NINT(geometry_info%mean_characteristic_length / 1000.0_wp), "km"
    ELSE
      WRITE(resolution_string,'(i4.4,a)') NINT(geometry_info%mean_characteristic_length), "m"
    ENDIF

!    write(0,*) SQRT(geometry_info%mean_cell_area) / 1000.0_wp, resolution_string
!    stop

  END FUNCTION get_resolution_string
  !-------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
  INTEGER FUNCTION parallel_read_geometry_info(ncid, geometry_info)
    INTEGER, INTENT(in) :: ncid
    TYPE(t_grid_geometry_info) :: geometry_info
    
#if ( defined (NOMPI) || defined (__ICON_GRID_GENERATOR__))
    parallel_read_geometry_info = read_geometry_info(ncid, geometry_info)
#else
    
    INTEGER :: netcd_status
    CHARACTER(*), PARAMETER :: method_name = "read_geometry_info"
            
    CALL set_default_geometry_info(geometry_info)
    parallel_read_geometry_info = -1
    
    netcd_status = p_nf_get_att_int(ncid, nf_global,'grid_geometry', geometry_info%geometry_type)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","grid_geometry")
      RETURN
    ENDIF
    
    netcd_status = p_nf_get_att_int(ncid, nf_global,'grid_cell_type', geometry_info%cell_type)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","grid_geometry")
      RETURN
    ENDIF
        
    netcd_status = p_nf_get_att_double(ncid, nf_global,'mean_edge_length', &
      & geometry_info%mean_edge_length)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","mean_edge_length")
      RETURN
    ENDIF
    
    netcd_status = p_nf_get_att_double(ncid, nf_global,'mean_dual_edge_length', &
      & geometry_info%mean_dual_edge_length)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","mean_dual_edge_length")
      RETURN
    ENDIF
    
    netcd_status = p_nf_get_att_double(ncid, nf_global,'mean_cell_area', &
      & geometry_info%mean_cell_area)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","mean_cell_area")
      RETURN
    ENDIF
    
    netcd_status = p_nf_get_att_double(ncid, nf_global,'mean_dual_cell_area', &
      & geometry_info%mean_dual_cell_area)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","mean_dual_cell_area")
      RETURN
    ENDIF
    
    netcd_status = p_nf_get_att_double(ncid, nf_global,'domain_length', &
      & geometry_info%domain_length)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","domain_length")
      RETURN
    ENDIF
    
    netcd_status = p_nf_get_att_double(ncid, nf_global,'domain_height', &
      & geometry_info%domain_height)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","domain_height")
      RETURN
    ENDIF

    netcd_status = p_nf_get_att_double(ncid, nf_global,'sphere_radius', geometry_info%sphere_radius)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","sphere_radius")
      RETURN
    ENDIF
    
    netcd_status = p_nf_get_att_double(ncid, nf_global,'domain_cartesian_center', &
      & geometry_info%center%x )
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","domain_cartesian_center")
      RETURN
    ENDIF
    
    ! return status ok
    parallel_read_geometry_info = 0
#endif

  END FUNCTION parallel_read_geometry_info
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
    
    netcd_status = nf_get_att_int(ncid, nf_global,'grid_cell_type', geometry_info%cell_type)
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
    
    netcd_status = nf_get_att_double(ncid, nf_global,'mean_dual_edge_length', &
      & geometry_info%mean_dual_edge_length)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","mean_dual_edge_length")
      RETURN
    ENDIF
    
    netcd_status = nf_get_att_double(ncid, nf_global,'mean_cell_area', &
      & geometry_info%mean_cell_area)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","mean_cell_area")
      RETURN
    ENDIF
    
    netcd_status = nf_get_att_double(ncid, nf_global,'mean_dual_cell_area', &
      & geometry_info%mean_dual_cell_area)
    IF (netcd_status /= nf_noerr) THEN
!       CALL finish("Cannot read","mean_dual_cell_area")
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
  
  !-------------------------------------------------------------------------
  INTEGER FUNCTION write_geometry_info(ncid, geometry_info)
    INTEGER, INTENT(in) :: ncid
    TYPE(t_grid_geometry_info) :: geometry_info

    write_geometry_info = -1

    CALL nf(nf_put_att_int      (ncid, nf_global, 'grid_geometry', nf_int, 1,     &
      & geometry_info%geometry_type))
    
    CALL nf(nf_put_att_int      (ncid, nf_global, 'grid_cell_type', nf_int, 1,     &
      & geometry_info%cell_type))
    
    CALL nf(nf_put_att_double(ncid, nf_global, 'mean_edge_length' , nf_double, 1, &
      & geometry_info%mean_edge_length))
      
    CALL nf(nf_put_att_double(ncid, nf_global, 'mean_dual_edge_length' , nf_double, 1, &
      & geometry_info%mean_dual_edge_length))
      
    CALL nf(nf_put_att_double  (ncid, nf_global, 'mean_cell_area' , nf_double, 1, &
      & geometry_info%mean_cell_area))
      
    CALL nf(nf_put_att_double  (ncid, nf_global, 'mean_dual_cell_area' , nf_double, 1, &
      & geometry_info%mean_dual_cell_area))
      
    CALL nf(nf_put_att_double   (ncid, nf_global, 'domain_length' , nf_double, 1, &
      & geometry_info%domain_length))
    
    CALL nf(nf_put_att_double   (ncid, nf_global, 'domain_height' , nf_double, 1, &
      & geometry_info%domain_height))
    
    CALL nf(nf_put_att_double   (ncid, nf_global, 'sphere_radius' , nf_double, 1, &
      & geometry_info%sphere_radius))
    
    CALL nf(nf_put_att_double  (ncid, nf_global, 'domain_cartesian_center', nf_double, 3, &
      & geometry_info%center%x))
      
    write_geometry_info = 0

  END FUNCTION write_geometry_info
  !-------------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE nf(return_status)
    INTEGER, INTENT(in) :: return_status

    IF (return_status /= nf_noerr) THEN
      CALL finish('mo_io_grid netCDF error', nf_strerror(return_status))
    ENDIF

  END SUBROUTINE nf
  !-------------------------------------------------------------------------

END MODULE mo_grid_geometry_info
!----------------------------------------------------------------------------









