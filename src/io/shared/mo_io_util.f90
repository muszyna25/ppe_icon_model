!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Guenther Zaengl, DWD
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial Revision by Daniel Reinert, DWD, 2012-03-22
!! - some IO-routines, which might be of future use, moved here from 
!!   the outdated output module mo_io_vlist. 
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_io_util

  USE mo_exception,             ONLY: finish
  USE mo_cdi,                   ONLY: FILETYPE_NC, FILETYPE_NC2, FILETYPE_NC4,         &
    &                                 FILETYPE_GRB, FILETYPE_GRB2
  USE mo_util_string,           ONLY: tolower
  USE mo_read_interface,        ONLY: nf

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'

  PUBLIC :: get_filetype
  PUBLIC :: get_file_extension
  PUBLIC :: read_netcdf_int_1d

  ! module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_io_util'

CONTAINS

  !-------------------------------------------------------------------------
  !> @return One of CDI's FILETYPE\_XXX constants. Possible values: 2
  !          (=FILETYPE\_GRB2), 4 (=FILETYPE\_NC2)
  !
  !  The file type is determined by the setting of the "filetype"
  !  namelist parameter in "initicon_nml". If this parameter has not
  !  been set, we try to determine the file type by its extension
  !  "*.grb*" or ".nc".
  !
  FUNCTION get_filetype(filename)
    INTEGER :: get_filetype
    CHARACTER(LEN=*), INTENT(IN) :: filename
    ! local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//'::get_filetype'
    INTEGER :: idx
    
    idx = INDEX(tolower(filename),'.nc')
    IF (idx==0) THEN
      idx = INDEX(tolower(filename),'.grb')
      IF (idx==0) THEN
        CALL finish(routine, "File type could not be determined!  File: " // trim (filename))
      ELSE
        get_filetype = FILETYPE_GRB2
      END IF
    ELSE
      get_filetype = FILETYPE_NC2
    END IF
  END FUNCTION get_filetype


  !> @return file extension string corresponding to file type
  !
  FUNCTION get_file_extension(filetype) RESULT(extn)
    CHARACTER(LEN=16)   :: extn
    INTEGER, INTENT(IN) :: filetype
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::get_file_extension"

    SELECT CASE (filetype)
    CASE (FILETYPE_NC)
      CALL finish(routine,'netCDF classic not supported')
    CASE (FILETYPE_NC2, FILETYPE_NC4)
      ! this is ok, both formats can write more than 2GB files
      extn = '.nc'
    CASE (FILETYPE_GRB)
      CALL finish(routine,'GRIB1 not supported')
    CASE (FILETYPE_GRB2)
      extn = '.grb'
    CASE default
      CALL finish(routine,'unknown output_type')
    END SELECT
  END FUNCTION get_file_extension


  !----------------------------------------------------------------
  !> Reads integer variables from NetCDF file
  !
  !  Note: opens and closes file and allocates output variables.
  !----------------------------------------------------------------
  SUBROUTINE read_netcdf_int_1d(filename, varname1, var1, opt_varname2, opt_var2, &
    &                           opt_attname, opt_attvar1, opt_attvar2) 
    CHARACTER(len=*),                INTENT(IN)    :: filename      ! NetCDF file name
    CHARACTER(len=*),                INTENT(IN)    :: varname1      ! variable name string
    INTEGER, ALLOCATABLE,            INTENT(INOUT) :: var1(:)       ! output data
    CHARACTER(len=*),     OPTIONAL,  INTENT(IN)    :: opt_varname2  ! variable name string
    INTEGER, ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: opt_var2(:)   ! output data
    CHARACTER(len=*),     OPTIONAL,  INTENT(IN)    :: opt_attname   ! optional variable attribute name
    INTEGER,              OPTIONAL,  INTENT(INOUT) :: opt_attvar1   ! attribute value for variable 1
    INTEGER,              OPTIONAL,  INTENT(INOUT) :: opt_attvar2   ! attribute value for variable 2
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine    = "read_netcdf_int_1d"
    INTEGER :: dimID, varid, ndims, dimids(1), dimlen, ncfileID
    LOGICAL :: l_exist

    ! consistency checks
    IF ((PRESENT(opt_var2) .AND. .NOT. PRESENT(opt_varname2)) .OR.  &
      & (PRESENT(opt_varname2) .AND. .NOT. PRESENT(opt_var2))) THEN
      CALL finish(routine, "Internal error!")
    END IF
    IF ((PRESENT(opt_attvar1) .AND. .NOT. PRESENT(opt_attname)) .OR. &
      & (PRESENT(opt_attvar2) .AND. .NOT. PRESENT(opt_attname))) THEN
      CALL finish(routine, "Internal error!")
    END IF

    ! --- open NetCDF file
    INQUIRE (FILE=filename, EXIST=l_exist)
    IF (.NOT. l_exist) THEN
      CALL finish(routine, 'file "'//TRIM(filename)//'" not found!')
    END IF
    CALL nf(nf_open(TRIM(FILENAME), NF_NOWRITE, ncfileID), routine)

    ! ----------------------
    ! --- variable "var1"
    ! ----------------------

    ! --- find out about variable dimensions:
    CALL nf(nf_inq_varid(ncfileID, TRIM(varname1), varid), routine)
    CALL nf(nf_inq_varndims(ncfileID, varID, ndims), routine)
    IF (ndims /= 1)  CALL finish(routine, "Variable '"//TRIM(varname1)//"' has more than one dimension!")
    CALL nf(nf_inq_vardimid(ncfileID, varID, dimids), routine)
    CALL nf(nf_inq_dimlen(ncfileID, dimids, dimlen), routine)

    ! --- allocate output variable, read data
    ALLOCATE(var1(dimlen))
    CALL nf(nf_get_var_int(ncfileID, varID, var1), routine)

    ! --- optional: read integer attribute
    IF (PRESENT(opt_attvar1)) THEN
      CALL nf(nf_get_att_int(ncfileID, varID, opt_attname, opt_attvar1), routine)
    END IF

    ! -------------------------
    ! --- variable "opt_var2"
    ! -------------------------

    IF (PRESENT(opt_var2)) THEN

      ! --- find out about variable dimensions:
      CALL nf(nf_inq_varid(ncfileID, TRIM(opt_varname2), varid), routine)
      CALL nf(nf_inq_varndims(ncfileID, varID, ndims), routine)
      IF (ndims /= 1)  CALL finish(routine, "Variable '"//TRIM(opt_varname2)//"' has more than one dimension!")
      CALL nf(nf_inq_vardimid(ncfileID, varID, dimids), routine)
      CALL nf(nf_inq_dimlen(ncfileID, dimids, dimlen), routine)
      
      ! --- allocate output variable, read data
      ALLOCATE(opt_var2(dimlen))
      CALL nf(nf_get_var_int(ncfileID, varID, opt_var2), routine)

      ! --- optional: read integer attribute
      IF (PRESENT(opt_attvar2)) THEN
        CALL nf(nf_get_att_int(ncfileID, varID, opt_attname, opt_attvar2), routine)
      END IF

    END IF

    ! --- close NetCDF file
    CALL nf(nf_close(ncfileID), routine)

  END SUBROUTINE read_netcdf_int_1d


END MODULE mo_io_util

