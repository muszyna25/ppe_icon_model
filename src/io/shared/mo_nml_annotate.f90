!>
!! Utility module: print namelist to file, annotating all changed values.
!!
!! The implementation is more or less a "poor man's approach". We write
!! and re-read the namelists in text form, since Fortran is missing the
!! language features for looping over namelist elements.
!!
!! @par Revision History
!!  by F. Prill, DWD (2013-06-13)
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_nml_annotate

#ifdef __ICON__  
  USE mo_io_units,    ONLY: find_next_free_unit, filename_max
  USE mo_exception,   ONLY: finish
#else
  USE mo_utilities,   ONLY: finish, find_next_free_unit,        &
    &                       filename_max
#endif
  USE mo_util_file,   ONLY: util_filesize, util_tmpnam, util_unlink
  USE mo_util_nml,    ONLY: util_annotate_nml
  
  IMPLICIT NONE
  PRIVATE
  

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_nml_annotate')

  INTEGER :: tmpnml = -1   !< file handle for temporary text file with defaults and settings
  
  PUBLIC :: temp_defaults
  PUBLIC :: temp_settings
  PUBLIC :: log_nml_settings 
  
CONTAINS
  
  !> Opens a new temporary text file.
  !
  FUNCTION temp_defaults()
    INTEGER :: temp_defaults
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::temp_defaults'
    LOGICAL                     :: lopen
    CHARACTER(len=filename_max) :: filename
    INTEGER                     :: flen

    IF (tmpnml == -1) THEN
      lopen = .FALSE.
    ELSE
      INQUIRE(tmpnml, OPENED=lopen)
    END IF
    IF (.NOT. lopen) THEN
      flen   = util_tmpnam(filename, filename_max)
      tmpnml = find_next_free_unit(10,100)
      IF (tmpnml < 0) THEN
        CALL finish(routine, "Failed call to find_next_free_unit.")
      END IF
      OPEN(UNIT=tmpnml, FILE=filename(1:flen), &
        &  ACCESS='sequential',   action='write', delim='apostrophe')
    END IF
    temp_defaults = tmpnml
  END FUNCTION temp_defaults
  
  
  !> Opens a new temporary text file.
  !
  FUNCTION temp_settings()
    INTEGER :: temp_settings
    ! Local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::temp_settings'
    LOGICAL :: lopen

    INQUIRE(tmpnml, OPENED=lopen)
    IF (.NOT. lopen) CALL finish(routine, "Internal error!")
    temp_settings = tmpnml
  END FUNCTION temp_settings
  
  
  !> Read defaults and settings into string buffer, compare them and do the
  !  print-out.
  !
  SUBROUTINE log_nml_settings(dst_filename)
    CHARACTER(LEN=*), INTENT(IN) ::dst_filename !< log file name (optional)
    ! Local variables
    CHARACTER(len=filename_max) :: tmp_filename
    INTEGER                     :: iret
    LOGICAL                     :: lopen

    INQUIRE(tmpnml, NAME=tmp_filename)
    INQUIRE(tmpnml, OPENED=lopen)
    IF (lopen) THEN
      CLOSE (tmpnml)
      tmpnml = -1
      ! call a C routine which reads and parses the temporary
      ! namelists with a finite state machine:
      iret = util_annotate_nml(TRIM(tmp_filename), TRIM(dst_filename))
      iret = util_unlink(TRIM(tmp_filename))
    END IF
  END SUBROUTINE log_nml_settings
  
END MODULE mo_nml_annotate
