!>
!!   * open_nml         - open namelist file
!!   * close_nml        - close namelist file
!!   * position_nml     - position namelist file for reading
!!   * open_nml_output  - open a new ASCII file to record the
!!                        actual values of the namelist variables.
!!   * close_nml_output - close that ASCII file.
!!
!! @par Revision History
!!   L. Kornblueh, MPI, March 2001, original source for ECHAM
!!   L. Kornblueh, MPI, January 2003, added error message for OPEN
!!  Modification by Hui Wan, MPI-M (2007-01-17):
!!   - added the ProTeX-style preamble
!!  Modification by L. Kornblueh, MPI, April 2005, added close_nml
!!  Modification by Thomas Heinze, DWD, (2005-05-03):
!!  - replaced FUNCTION position_nml by one provided by L.Kornblueh
!!  Modification by Hui Wan, MPI-M  (2009-03-17):
!!  - added subroutines open_nml_output and close_nml_output
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
MODULE mo_namelist
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!

  USE mo_util_string,   ONLY: tolower
  USE mo_exception,     ONLY: finish
  USE mo_io_units,      ONLY: nnml, nnml_output, find_next_free_unit

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

! !PUBLIC ENTITIES
  PUBLIC :: position_nml           ! function  : position namelist file
  PUBLIC :: open_nml               ! subroutine: open default namelist file
  PUBLIC :: close_nml              ! subroutine: close namelist file
  PUBLIC :: POSITIONED, MISSING,  &! return values from position_nml
            LENGTH_ERROR, READ_ERROR
  PUBLIC :: open_nml_output
  PUBLIC :: close_nml_output

!  return values of function 'position_nml':

  INTEGER, PARAMETER :: POSITIONED   =  0 ! file pointer set to namelist group
  INTEGER, PARAMETER :: MISSING      =  1 ! namelist group is missing
  INTEGER, PARAMETER :: LENGTH_ERROR =  2
  INTEGER, PARAMETER :: READ_ERROR   =  3

CONTAINS

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!

  !>
  !!  opens the namelist file.
  !!
  !!
  SUBROUTINE open_nml (file)


    CHARACTER(len=*) ,INTENT(in) :: file

    INTEGER :: istat

!-------------------------------------------------------------------------

    OPEN (nnml, file=file, iostat=istat, status='old', action='read', &
         delim='apostrophe')

    IF (istat /= 0) THEN
      CALL finish ('open_nml','Could not open '//TRIM(file))
    END IF

  END SUBROUTINE open_nml

!-------------------------------------------------------------------------
!
!
  !>
  !!               close namelist file.
  !!
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, March 2001
  !!
  SUBROUTINE close_nml
!

    INTEGER :: istat

!-------------------------------------------------------------------------

    CLOSE (nnml, IOSTAT=istat)

    IF (istat /= 0) THEN
      CALL finish ('close_nml','Could not close namelist file.')
    END IF

  END SUBROUTINE close_nml

!-------------------------------------------------------------------------
!

  !>
  !!               set file pointer to begin of namelist for reading
  !!               namelist /name/ (case independent).
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, March 2001
  !!
  SUBROUTINE position_nml (name, lrewind, status)
!
    CHARACTER(len=*), INTENT(in)            :: name    ! namelist group name
    LOGICAL,          INTENT(in)  ,OPTIONAL :: lrewind ! default: true

    INTEGER,          INTENT(out) ,OPTIONAL :: status  ! error return value

    CHARACTER(len=256) :: yline    ! line read
    CHARACTER(len=256) :: test     ! uppercase namelist group name
    INTEGER            :: stat     ! local copy of status variable
    INTEGER            :: ios      ! status variable from read operation
    LOGICAL            :: lrew     ! local copy of rewind flag
    INTEGER            :: len_name ! length of requested namelist group name
    CHARACTER          :: ytest    ! character to test for delimiter
    CHARACTER(len=12)  :: code     ! error code printed
    INTEGER            :: ind      ! index from index routine
    INTEGER            :: indc     ! index of comment character (!)

!-------------------------------------------------------------------------

    lrew  = .TRUE.   ; IF (PRESENT(lrewind)) lrew  = lrewind
    stat  =  MISSING
    code  = 'MISSING'

    len_name = LEN_TRIM (name)

    IF (len_name > LEN(test)) THEN
      stat =  LENGTH_ERROR
      code = 'LENGTH_ERROR'
    ENDIF

    test = tolower (name)

    ! Reposition file at beginning:

    IF (lrew) REWIND (nnml)

    ! Search start of namelist

    DO
      IF (stat /= MISSING) EXIT

      yline = ' '

      READ (nnml,'(a)',IOSTAT=ios) yline
      IF (ios < 0) THEN
        EXIT  ! MISSING
      ELSE IF (ios > 0) THEN
        stat =  READ_ERROR
        code = 'READ_ERROR'
        EXIT
      END IF

      yline = tolower (yline)

      ind = INDEX(yline,'&'//TRIM(test))

      IF (ind == 0) CYCLE

      indc = INDEX(yline,'!')

      IF (indc > 0 .AND. indc < ind) CYCLE

      ! test for delimiter

      ytest = yline(ind+len_name+1:ind+len_name+1)

      IF ( (LGE(ytest,'0') .AND. LLE(ytest,'9')) .OR. &
           (LGE(ytest,'a') .AND. LLE(ytest,'z')) .OR. &
            ytest == '_'                         .OR. &
           (LGE(ytest,'A') .AND. LLE(ytest,'Z'))) THEN
        CYCLE
      ELSE
        stat = POSITIONED
        BACKSPACE (nnml)
        EXIT
      END IF
    ENDDO

    IF (PRESENT (status)) status = stat
      SELECT CASE (stat)
      CASE (POSITIONED)
        RETURN
      CASE (MISSING)
        IF (PRESENT (status)) RETURN
      END SELECT

    CALL finish ('position_nml','namelist /'//TRIM(test)//'/ '//code)

  END SUBROUTINE position_nml

!-------------------------------------------------------------------------
!

  !>
  !!  opens an ASCII file which will later contain all the namelist variables
  !!  and their actual values used in the simulation.
  !!
  SUBROUTINE open_nml_output (file)


    CHARACTER(len=*) ,INTENT(in) :: file

    INTEGER :: istat

!-------------------------------------------------------------------------

    nnml_output = find_next_free_unit(10,20)

    OPEN (nnml_output, FILE=TRIM(file), IOSTAT=istat, DELIM='apostrophe')

    IF (istat /= 0) THEN
      CALL finish ('open_nml_output','Could not open '//TRIM(file))
    END IF

  END SUBROUTINE open_nml_output

!-------------------------------------------------------------------------
!
!
  !>
  !!  close the ASCII output that contains all the namelist
  !!  variables and their actual values used in the simulation.
  !!
  !! @par Revision History
  !!  Hui Wan, MPI-M, 2009-03-17
  !!
  SUBROUTINE close_nml_output
!

    INTEGER :: istat

!-------------------------------------------------------------------------

    CLOSE (nnml_output, IOSTAT=istat)

    IF (istat /= 0) THEN
      CALL finish ('close_nml_output','Could not close the namelist output')
    END IF

  END SUBROUTINE close_nml_output

END MODULE mo_namelist
