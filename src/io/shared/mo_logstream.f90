!>
!! Utility module: open log (ASCII) file where messages or other
!! files' contents may be dumped, annotated by a time stamp and MPI
!! rank info.
!!
!! @author F. Prill, DWD
!!
!!
!! @par Revision History
!! Initial revision: 2018-10-09 : F. Prill, DWD
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_logstream

  USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH
  USE mo_io_units,       ONLY: find_next_free_unit
  USE mo_mpi,            ONLY: get_my_global_mpi_id

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_logstream
  PUBLIC :: logstream
  PUBLIC :: SUCCESS, ERROR_NOT_OPEN, ERROR_FILE_NOT_FOUND, ERROR_FILE_NOT_READABLE

  ! module name string (for debugging purposes):
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_logstream'

  ! log stream status:
  ENUM, BIND(C)
    ENUMERATOR ::                &
    STATUS_UNDEFINED        = 0, &
    STATUS_OPEN             = 1, &
    STATUS_CLOSED           = 2
  END ENUM

  ! log stream errors:
  ENUM, BIND(C)
    ENUMERATOR ::                 &
    SUCCESS                 =  0, &
    ERROR_NOT_OPEN          = -1, &
    ERROR_FILE_NOT_FOUND    = -2, &
    ERROR_FILE_NOT_READABLE = -3 
  END ENUM


  !> log stream object: log (ASCII) file where messages or other
  !  files' contents may be dumped, annotated by a time stamp and MPI
  !  rank info.
  !
  TYPE :: t_logstream
    INTEGER                        :: status = STATUS_UNDEFINED !< Flag: .TRUE. = opened
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: filename                  !< name of the log stream file.
    INTEGER                        :: funit                     !< output file unit.
  CONTAINS
    PROCEDURE :: init        => logstream_init                  !< initialize, open logstream file.
    PROCEDURE :: finalize    => logstream_finalize              !< destructor.

    PROCEDURE :: append      => logstream_append                !< append message string.
    PROCEDURE :: append_file => logstream_append_file           !< append some file's contents to log file.
  END TYPE t_logstream

  ! local instance of logstream object
  TYPE(t_logstream) :: logstream

  ! temporary variable for time stamp construction:
  INTEGER, PARAMETER               :: MAX_TIMESTAMP_LEN =  40
  CHARACTER(LEN=MAX_TIMESTAMP_LEN) :: tmp_timestamp
  
  ! line buffer (file read-in):
  INTEGER, PARAMETER               :: BUFFER_LEN        = 512
  CHARACTER(len=BUFFER_LEN)        :: buffer

CONTAINS

  ! Constructor.
  !
  SUBROUTINE logstream_init(logstream, filename_prefix)
    CLASS(t_logstream), INTENT(INOUT) :: logstream         !< log stream object
    CHARACTER(LEN=*),   INTENT(IN)    :: filename_prefix   !< output file name (prefix).
    ! local variables
    INTEGER :: pe_id
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: filename

    ! generate log stream file from prefix "filename" and PE number:
    pe_id = get_my_global_mpi_id()
    WRITE (filename, '(a,a,i0)') TRIM(filename_prefix), ".PE_", pe_id

    logstream%funit  = find_next_free_unit(10,100)
    OPEN (unit=logstream%funit, file=TRIM(filename), action="write",status="replace")
    logstream%status = STATUS_OPEN

    CALL logstream%append("Log stream opened.")
  END SUBROUTINE logstream_init


  ! Destructor.
  !
  SUBROUTINE logstream_finalize(logstream)
    CLASS(t_logstream), INTENT(INOUT) :: logstream  !< log stream object

    CALL logstream%append("Closing log stream")

    CLOSE(logstream%funit)
    logstream%status = STATUS_CLOSED
  END SUBROUTINE logstream_finalize


  ! Append message string.
  !
  SUBROUTINE logstream_append(logstream, msg, opt_ierr)
    CLASS(t_logstream), INTENT(INOUT) :: logstream  !< log stream object.
    CHARACTER(LEN=*),   INTENT(IN)    :: msg        !< message string.
    INTEGER, INTENT(OUT), OPTIONAL    :: opt_ierr   !< error code (0=SUCCESS).
    ! local variables
    INTEGER :: ierr

    ierr = SUCCESS
    IF (logstream%status == STATUS_OPEN) THEN
      CALL get_timestamp(logstream, tmp_timestamp)
      WRITE (logstream%funit, '(2a)')  tmp_timestamp, TRIM(msg)
    ELSE
      ierr = ERROR_NOT_OPEN
    END IF

    IF (PRESENT(opt_ierr))  opt_ierr = ierr
  END SUBROUTINE logstream_append


  ! Append some file's contents to log file.
  !
  SUBROUTINE logstream_append_file(logstream, filename, opt_ierr)
    CLASS(t_logstream), INTENT(INOUT) :: logstream  !< log stream object.
    CHARACTER(LEN=*),   INTENT(IN)    :: filename   !< source file name.
    INTEGER, INTENT(OUT), OPTIONAL    :: opt_ierr   !< error code (0=SUCCESS).
    ! local variables
    LOGICAL :: l_exist
    INTEGER :: io_error, iunit, ierr

    ierr = SUCCESS
    IF (logstream%status /= STATUS_OPEN) THEN
      ierr = ERROR_NOT_OPEN
      IF (PRESENT(opt_ierr))  opt_ierr = ierr
      RETURN
    END IF
    INQUIRE (FILE=filename, EXIST=l_exist)
    IF (.NOT.l_exist) THEN
      ierr = ERROR_FILE_NOT_FOUND
      IF (PRESENT(opt_ierr))  opt_ierr = ierr
      RETURN
    END IF

    ! print time stamp:
    CALL get_timestamp(logstream, tmp_timestamp)
    WRITE (logstream%funit, '(2a)') tmp_timestamp, &
      &                             "----------------------------------------------------------------"

    ! print input file's contents:
    iunit = find_next_free_unit(10,100)
    OPEN(unit=iunit, file=filename, status='old',action='read', iostat=io_error) 
    IF ( io_error /= 0) THEN
      ierr = ERROR_FILE_NOT_READABLE
      IF (PRESENT(opt_ierr))  opt_ierr = ierr
      RETURN
    END IF

    WRITE (logstream%funit, '(4a)') tmp_timestamp, "Contents of file '", TRIM(filename), "':"
    DO
      READ (iunit, '(A)', iostat=io_error) buffer
      IF (io_error /= 0)  EXIT

      WRITE (logstream%funit, '(2a)') tmp_timestamp, TRIM(buffer)
    END DO

    CLOSE(iunit) 
    IF (PRESENT(opt_ierr))  opt_ierr = ierr
  END SUBROUTINE logstream_append_file


  ! Auxiliary routine: generate time stamp preceding each log message.
  !
  SUBROUTINE get_timestamp(logstream, tstamp)
    TYPE(t_logstream), INTENT(INOUT) :: logstream  !< log stream object.
    CHARACTER(LEN=*),  INTENT(INOUT) :: tstamp     !< time stamp (result).
    ! local variables
    INTEGER           :: date_time(8), pe_id
    CHARACTER(LEN=10) :: b

    ! retrieve PE number:
    pe_id = get_my_global_mpi_id()
    ! retrieve date and time:
    CALL DATE_AND_TIME(b, values=date_time)
    ! "date_time" holds year, month of the year, day of the month,
    ! time diff wrt. UTC, hour of day, minutes, seconds, milliseconds.

    ! construct time stamp:
    tstamp = ""
    WRITE (tstamp,'(4(a,i0),4(a,i2),a)') &
      &       "PE ", pe_id, " :: ", date_time(1), "-", date_time(2), "-", date_time(3), &
      &       " ", date_time(5), ":", date_time(6), ":", date_time(7), " :: "
  END SUBROUTINE get_timestamp

END MODULE mo_logstream
