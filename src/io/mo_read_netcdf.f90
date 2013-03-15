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
  USE mo_parallel_config,    ONLY: nproma
  USE mo_io_units,           ONLY: filename_max

IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: read_netcdf_cells_2D

  INTERFACE read_netcdf_cells_2D
    MODULE PROCEDURE read_netcdf_REAL_CELLS_2D_filename
    MODULE PROCEDURE read_netcdf_REAL_CELLS_2D_fileid
  END INTERFACE

  INTEGER, PARAMETER :: MAX_VAR_DIMS = NF_MAX_VAR_DIMS
!-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION read_netcdf_REAL_CELLS_2D_filename(filename, variable_name, fill_array, patch)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp), POINTER            :: fill_array(:,:)
    TYPE(t_patch)                :: patch

    INTEGER :: ncid
    INTEGER :: return_status    
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_netcdf:read_netcdf_REAL_CELLS_2D_filename'

    ncid = netcdf_open_input(filename)
    read_netcdf_REAL_CELLS_2D_filename = &
      & read_netcdf_REAL_CELLS_2D_fileid(ncid, variable_name, fill_array, patch)
    return_status = netcdf_close(ncid)
                              
  END FUNCTION read_netcdf_REAL_CELLS_2D_filename
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION read_netcdf_REAL_CELLS_2D_fileid(ncid, variable_name, fill_array, patch)
    INTEGER, INTENT(IN)          :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp), POINTER            :: fill_array(:,:)
    TYPE(t_patch)                :: patch

    INTEGER :: total_number_of_cells
    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    INTEGER :: return_status
    REAL(wp), POINTER :: tmp_array(:)
    
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_read_netcdf:read_netcdf_REAL_CELLS_2D_fileid'

    total_number_of_cells = patch%n_patch_cells_g
    
    IF( my_process_is_mpi_workroot()  ) THEN
      CALL nf(netcdf_inq_var(ncid, variable_name, varid, var_type, var_dims, var_size))
      
      IF (var_dims /= 1 .OR. var_size(1) /= total_number_of_cells) THEN
        write(0,*) "var_dims = ", var_dims, " var_size=", var_size, &
          & " total_number_of_cells=", total_number_of_cells
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF
      
    ENDIF
    
    ALLOCATE( tmp_array(total_number_of_cells), stat=return_status )
    IF (return_status /= success) THEN
      CALL finish (method_name, 'ALLOCATE( tmp_array )')
    ENDIF
    
    IF( my_process_is_mpi_workroot()) THEN
      CALL nf(nf_get_var_double(ncid, varid, tmp_array(:)))
    ENDIF

    IF (.NOT. ASSOCIATED(fill_array)) THEN
      ALLOCATE( fill_array(nproma, patch%nblks_c), stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( fill_array )')
      ENDIF
    ENDIF
    
    CALL scatter_cells(tmp_array, fill_array, patch)
    
    DEALLOCATE(tmp_array)    
                              
  END FUNCTION read_netcdf_REAL_CELLS_2D_fileid
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
    INTEGER, INTENT(OUT) :: var_size(MAX_VAR_DIMS)

    INTEGER  :: number_of_attributes  
    INTEGER :: var_dims_reference(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: check_var_name
    INTEGER :: i, return_status

    netcdf_inq_var = -1
    IF ( .NOT. my_process_is_mpi_workroot() ) RETURN


    netcdf_inq_var = nf_inq_varid(ncid, name, varid)
    CALL nf(netcdf_inq_var)
    netcdf_inq_var = nf_inq_var (ncid, varid, check_var_name, var_type, var_dims, &
      & var_dims_reference, number_of_attributes)
    DO i=1, var_dims
!       return_status = nf_inq_dimlen(ncid, var_dims_reference(i), var_size(i))
      CALL nf(nf_inq_dimlen(ncid, var_dims_reference(i), var_size(i)))
    ENDDO
!     write(0,*) " Read var_dims, var_size:",  var_dims, var_size
!     write(0,*) " check_var_name:",  check_var_name
!     write(0,*) " name:", name    
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
