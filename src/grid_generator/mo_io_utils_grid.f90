!>
!!   Set of IO utils for reading/writing grid info 
!!
!!
!! @par Revision History
!!  Developed  by Leonidas Linardakis, MPIM, 20012-12
!!  Modifications by Thomas Heinze (2005-08-29):
!!
!! @par Copyright
!! 2002-2008 by DWD and MPI-M
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
!! $Id: n/a$
!!
MODULE mo_io_utils_grid

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: find_next_free_unit, filename_max
  USE mo_exception,          ONLY: message_text, message, finish, warning
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates, t_geographical_coordinates
  USE mo_grid_geometry_info, ONLY: sphere_geometry, planar_torus_geometry, &
    & t_planar_torus_geometry_info

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: read_planar_torus_info, write_planar_torus_info, nf
  !--------------------------------------------------------------------

CONTAINS


  !--------------------------------------------------------------------
  SUBROUTINE nf(return_status)
    INTEGER, INTENT(in) :: return_status

    IF (return_status /= nf_noerr) THEN
      CALL finish('mo_io_grid netCDF error', nf_strerror(return_status))
    ENDIF

  END SUBROUTINE nf
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE read_planar_torus_info(ncid, planar_torus_info)
    INTEGER, INTENT(in) :: ncid
    TYPE(t_planar_torus_geometry_info) :: planar_torus_info
    
    INTEGER :: netcd_status
    CHARACTER(*), PARAMETER :: method_name = "read_planar_torus_info"
        
    netcd_status = nf_get_att_double(ncid, nf_global,'planar_torus_cell_edge_length', &
      & planar_torus_info%cell_edge_length)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "planar_torus_cell_edge_length")
    
    netcd_status = nf_get_att_double(ncid, nf_global,'planar_torus_length', &
      & planar_torus_info%length)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "planar_torus_length")
    
    netcd_status = nf_get_att_double(ncid, nf_global,'planar_torus_height', &
      & planar_torus_info%height)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "planar_torus_height")

    netcd_status = nf_get_att_double(ncid, nf_global,'planar_torus_cartesian_center', &
      & planar_torus_info%center%x )
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "planar_torus_cartesian_center")
    
  END SUBROUTINE read_planar_torus_info
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE write_planar_torus_info(ncid, planar_torus_info)
    INTEGER, INTENT(in) :: ncid
    TYPE(t_planar_torus_geometry_info) :: planar_torus_info

    CALL nf(nf_put_att_double  (ncid, nf_global, 'planar_torus_cell_edge_length' , nf_double, 1, &
      & planar_torus_info%cell_edge_length))
      
    CALL nf(nf_put_att_double  (ncid, nf_global, 'planar_torus_length' , nf_double, 1, &
      & planar_torus_info%length))
    
    CALL nf(nf_put_att_double  (ncid, nf_global, 'planar_torus_height' , nf_double, 1, &
      & planar_torus_info%height))
    
    CALL nf(nf_put_att_double  (ncid, nf_global, 'planar_torus_cartesian_center', nf_double, 3, &
      & planar_torus_info%center%x))

  END SUBROUTINE write_planar_torus_info
  !-------------------------------------------------------------------------

END MODULE mo_io_utils_grid



