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

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: get_filetype
  PUBLIC :: get_file_extension

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

END MODULE mo_io_util

