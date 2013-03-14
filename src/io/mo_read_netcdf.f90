!>
!! This module provides basic methods for reading 
!! a NetCDF file in a parallel or sequential in a transparent way.
!!
!! @par Revision History
!! Initial version by Leonidas Linardakis, March 2013
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
MODULE mo_read_netcdf

  USE mo_kind
  USE mo_mpi
  USE mo_gather_scatter,     ONLY: scatter_cells
  USE mo_model_domain,       ONLY: t_patch
  USE mo_exception,          ONLY: message_text, message, warning, finish, em_warn
  USE mo_impl_constants,     ONLY: success, max_char_length

IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: read_netcdf_cells_2D

  INTERFACE read_netcdf_cells_2D
    MODULE PROCEDURE read_netcdf_REAL_CELLS_2D
  END INTERFACE

!-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION read_netcdf_REAL_CELLS_2D(filename, variable_name, fill_array, patch)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp), POINTER :: fill_array(:,:)
    TYPE(t_patch) :: patch

    INTEGER :: total_number_of_cells
    INTEGER :: ncid, varid, var_type, var_dims
    INTEGER :: var_size(NF_MAX_VAR_DIMS)
    INTEGER :: return_status
    REAL(wp), POINTER :: tmp_array(:)
    
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_netcdf:read_netcdf_REAL_CELLS_2D'

    total_number_of_cells = patch%n_patch_cells_g
    
    IF( my_process_is_mpi_workroot()  ) THEN
      ncid = netcdf_open_input(filename)
      CALL nf(netcdf_inq_var(ncid, variable_name, varid, var_type, var_dims, var_size))
      
      IF (var_dims /= 1 .OR. var_size(1) /= total_number_of_cells) &
        & CALL finish(method_name, "Dimensions mismatch")
    ENDIF
    
    ALLOCATE( tmp_array(total_number_of_cells), stat=return_status )
    IF (return_status /= success) THEN
      CALL finish (method_name, 'ALLOCATE( tmp_array )')
    ENDIF
    
    IF( my_process_is_mpi_workroot()) THEN
      CALL nf(nf_get_var_double(ncid, varid, tmp_array(:)))
      CALL nf(netcdf_close(ncid))
    ENDIF
    
    CALL scatter_cells(tmp_array, fill_array, patch)
    
    DEALLOCATE(tmp_array)    
                              
  END FUNCTION read_netcdf_REAL_CELLS_2D
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_open_input(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: ncid

    IF( my_process_is_mpi_workroot()  ) THEN
        CALL nf(nf_open(TRIM(filename), nf_nowrite, ncid))
    ELSE
        ncid = -1 ! set it to an invalid value
    ENDIF

    netcdf_open_input = ncid
    
  END FUNCTION netcdf_open_input
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_close(ncid)
    INTEGER, INTENT(IN) :: ncid

    netcdf_close = -1
    IF( my_process_is_mpi_workroot()  ) THEN
        netcdf_close = nf_close(ncid)
    ENDIF

  END FUNCTION netcdf_close
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_inq_var(ncid, name, varid, var_type, var_dims, var_size)
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: name
    
    INTEGER, INTENT(OUT) :: varid, var_type, var_dims
    INTEGER, INTENT(OUT) :: var_size(NF_MAX_VAR_DIMS)

    INTEGER :: number_of_attributes  
    CHARACTER, POINTER :: NULL_CHAR_POINTER
    INTEGER :: res

    netcdf_inq_var = -1
    IF ( .NOT. my_process_is_mpi_workroot() ) RETURN

    NULLIFY(NULL_CHAR_POINTER)

    netcdf_inq_var = nf_inq_varid(ncid, name, varid)
    CALL nf(netcdf_inq_var)
    netcdf_inq_var = nf_inq_var (ncid, varid, NULL_CHAR_POINTER, var_type, var_dims, var_size, number_of_attributes)
    
    CALL nf(netcdf_inq_var)

  END FUNCTION netcdf_inq_var
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE nf(STATUS, warnonly, silent)
    
    INTEGER, INTENT(in)           :: STATUS
    LOGICAL, INTENT(in), OPTIONAL :: warnonly
    LOGICAL, INTENT(in), OPTIONAL :: silent
    
    LOGICAL :: lwarnonly, lsilent
    
    lwarnonly = .FALSE.
    lsilent   = .FALSE.
    IF(PRESENT(warnonly)) lwarnonly = .TRUE.
    IF(PRESENT(silent))   lsilent   = silent
    
    IF (lsilent) RETURN
    IF (STATUS /= nf_noerr) THEN
      IF (lwarnonly) THEN
        CALL message('mo_read_netcdf error', nf_strerror(STATUS), &
          & level=em_warn)
      ELSE
        CALL finish('mo_read_netcdf', nf_strerror(STATUS))
      ENDIF
    ENDIF
    
  END SUBROUTINE nf
  !-------------------------------------------------------------------------

END MODULE mo_read_netcdf
