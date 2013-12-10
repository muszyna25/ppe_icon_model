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
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_io_util

  USE mo_exception,             ONLY: finish
  USE mo_cdi_constants,         ONLY: FILETYPE_NC, FILETYPE_NC2, FILETYPE_NC4,         &
    &                                 FILETYPE_GRB, FILETYPE_GRB2
  USE mo_util_string,           ONLY: tolower

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: get_filetype
  PUBLIC :: get_file_extension

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

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
    CHARACTER(len=*), PARAMETER :: routine = 'mo_nh_initicon:get_filetype'
    INTEGER :: idx
    
    idx = INDEX(tolower(filename),'.nc')
    IF (idx==0) THEN
      idx = INDEX(tolower(filename),'.grb')
      IF (idx==0) THEN
        CALL finish(routine, "File type could not be determined!")
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

