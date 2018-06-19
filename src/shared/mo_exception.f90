!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! - Initial version taken from ECHAM6;
!! - Subroutines print_status and print_value taken from mo_submodel of
!!   ECHAM6 by Hui Wan (2010-07-14).
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_exception

  USE mo_io_units, ONLY: nerr, nlog, filename_max
  USE mo_mpi,      ONLY: run_is_global_mpi_parallel, abort_mpi, my_process_is_stdio, &
    & get_my_global_mpi_id, get_my_mpi_work_id, get_glob_proc0, proc_split, comm_lev
  USE mo_kind,     ONLY: wp
  USE mo_impl_constants,  ONLY: MAX_CHAR_LENGTH
  

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: message_text
  PUBLIC :: message, warning, finish, print_status, print_value
  PUBLIC :: em_none, em_info, em_warn, em_error, em_param, em_debug
  PUBLIC :: open_log, close_log
  PUBLIC :: debug_messages_on, debug_messages_off
  PUBLIC :: number_of_warnings, number_of_errors
  PUBLIC :: get_filename_noext
  PUBLIC :: msg_timestamp

  INTEGER, PARAMETER :: em_none  = 0   !< normal message
  INTEGER, PARAMETER :: em_info  = 1   !< informational message
  INTEGER, PARAMETER :: em_warn  = 2   !< warning message: number of warnings counted
  INTEGER, PARAMETER :: em_error = 3   !< error message: number of errors counted
  INTEGER, PARAMETER :: em_param = 4   !< report parameter value
  INTEGER, PARAMETER :: em_debug = 5   !< debugging message

  CHARACTER(len=MAX_CHAR_LENGTH) :: message_text = ''

  LOGICAL :: l_debug     = .FALSE.
  LOGICAL :: l_log       = .FALSE.

  !> Flag. If .TRUE., precede output messages by time stamp.
  !  This parameters is set by mo_run_nml
  LOGICAL :: msg_timestamp


  INTEGER :: number_of_warnings  = 0
  INTEGER :: number_of_errors    = 0

  INTERFACE print_value            !< report on a parameter value
    MODULE PROCEDURE print_lvalue  !< logical
    MODULE PROCEDURE print_ivalue  !< integer
    MODULE PROCEDURE print_rvalue  !< real
  END INTERFACE

CONTAINS
  !>
  SUBROUTINE debug_messages_on
    l_debug = .TRUE.
  END SUBROUTINE debug_messages_on
  !-------------
  !>
  SUBROUTINE debug_messages_off
    l_debug = .FALSE.
  END SUBROUTINE debug_messages_off
  !-------------
  
  !-------------
  !>
  !!
  FUNCTION get_filename_noext(name) result(filename)
    CHARACTER (len=*), INTENT(in) :: name    
    CHARACTER(LEN=filename_max) :: filename

    INTEGER :: end_name

    end_name = INDEX(name, '.', back=.true.)
    IF (end_name > 0) THEN
      filename = name(1:end_name-1)
    ELSE
      filename = TRIM(name)
    ENDIF

  END FUNCTION get_filename_noext
  !-------------
  
  !-------------
  !>
  !! @brief Finish model simulation and report the reason
  !!
  SUBROUTINE finish (name, text, exit_no)

#ifdef __INTEL_COMPILER
    USE ifcore
#endif

    CHARACTER(len=*), INTENT(in)           :: name
    CHARACTER(len=*), INTENT(in), OPTIONAL :: text
    INTEGER,          INTENT(in), OPTIONAL :: exit_no

    INTEGER           :: iexit

#ifndef __STANDALONE
    EXTERNAL :: util_exit
#if ! (defined (__SX__) || defined (__INTEL_COMPILER) || defined (__xlC__))
    EXTERNAL :: util_backtrace
#endif
#endif

    WRITE (nerr,'(/,80("="),/)')
    IF (l_log) WRITE (nlog,'(/,80("="),/)')

    IF (PRESENT(exit_no)) THEN
       iexit = exit_no
    ELSE
       iexit = 1     ! POSIX defines this as EXIT_FAILURE
    END IF

    IF (PRESENT(text)) THEN
      IF (iexit == 1) THEN
        WRITE (nerr,'(a,a,a,a)') 'FATAL ERROR in ', TRIM(name), ': ', TRIM(text)
      ELSE
        WRITE (nerr,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
      ENDIF
      IF (l_log) WRITE (nlog,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
    ELSE
      IF (iexit == 1) THEN
        WRITE (nerr,'(a,a)') 'FATAL ERROR in ', TRIM(name)
      ELSE
        WRITE (nerr,'(1x,a)') TRIM(name)
      ENDIF
      IF (l_log) WRITE (nlog,'(a,a)') TRIM(name)
    ENDIF

!     IF (my_process_is_mpi_parallel()) THEN
      WRITE (nerr,'(1x,a,i0)') 'FINISH called from PE: ', get_my_global_mpi_id()
      IF (l_log) WRITE (nlog,'(1x,a,i0)') 'FINISH called from PE: ', get_my_global_mpi_id()
!     ENDIF

    WRITE (nerr,'(/,80("-"),/,/)')
    IF (l_log) WRITE (nlog,'(/,80("-"),/,/)')

#ifdef __INTEL_COMPILER
    CALL tracebackqq
#else
#ifdef __xlC__
    CALL xl__trbk
#else
#ifdef __SX__
    CALL mesput('Traceback: ', 11, 1)
#else
#ifndef __STANDALONE
    CALL util_backtrace
#endif
#endif
#endif
#endif

    WRITE (nerr,'(/,80("="),/)')
    IF (l_log) WRITE (nlog,'(/,80("="),/)')

    IF (run_is_global_mpi_parallel()) THEN
      CALL abort_mpi
    ELSE
#ifndef __STANDALONE
       CALL util_exit(iexit)
#else
       STOP 'mo_exception: finish ..'
#endif
    END IF

  END SUBROUTINE finish
  !-------------
  

  !-------------
  !>
  !!
  SUBROUTINE warning (name, text)

    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text

    CALL message (name, text, level=em_warn)

  END SUBROUTINE warning
  !-------------
  

  !-------------
  !>
  !!
  SUBROUTINE message (name, text, out, level, all_print, adjust_right)

    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text
    INTEGER,           INTENT(in), OPTIONAL :: out
    INTEGER,           INTENT(in), OPTIONAL :: level
    LOGICAL,           INTENT(in), OPTIONAL :: all_print
    LOGICAL,           INTENT(in), OPTIONAL :: adjust_right

    INTEGER :: iout
    INTEGER :: ilevel
    LOGICAL :: lprint
    LOGICAL :: ladjust
    LOGICAL :: lactive

    CHARACTER(len=32) :: prefix

    CHARACTER(len=LEN(message_text)) :: write_text

    CHARACTER(len=10) :: ctime
    CHARACTER(len=8)  :: cdate
 
    IF (PRESENT(all_print)) THEN
      lprint = all_print
    ELSE
      lprint = .FALSE.
    ENDIF

    IF (PRESENT(adjust_right)) THEN
      ladjust = adjust_right
    ELSE
      ladjust = .FALSE.
    ENDIF

    IF (PRESENT(out)) THEN
      iout = out
    ELSE
      iout = nerr
    END IF

    IF (PRESENT(level)) THEN
      ilevel = level
    ELSE
      ilevel = em_none
    END IF

    SELECT CASE (ilevel)
    CASE (em_none)  ; prefix = '        '
    CASE (em_info)  ; prefix = 'INFO   :'
    CASE (em_warn)  ; prefix = 'WARNING:' ; number_of_warnings  = number_of_warnings+1
    CASE (em_error) ; prefix = 'ERROR  :' ; number_of_errors    = number_of_errors+1
    CASE (em_param) ; prefix = '---     '
    CASE (em_debug) ; prefix = 'DEBUG  :'
    END SELECT

    IF (.NOT. ladjust) THEN
      message_text = TRIM(ADJUSTL(text))
    ELSE
      message_text = TRIM(text)
    ENDIF
    IF (name /= '')  THEN
      message_text = TRIM(name) // ': ' // TRIM(message_text)
    ENDIF
    IF (ilevel > em_none) THEN
      message_text = TRIM(prefix) // ' ' // TRIM(message_text)
    ENDIF

    IF (run_is_global_mpi_parallel() .AND. &
      & (l_debug .OR. ilevel == em_warn .OR. ilevel == em_error)) THEN
      WRITE(write_text,'(1x,a,i6,a,a)') 'PE ', get_my_global_mpi_id(), ' ', &
        & TRIM(message_text)
      lprint = .TRUE.
    ELSE

      IF ( proc_split ) THEN
        IF ( (comm_lev>0)  .AND.   &
        & (get_my_mpi_work_id() == get_glob_proc0()) ) THEN
        WRITE(write_text,*) 'PROC_SPLIT: PE ', get_glob_proc0(), ' ', &
          & TRIM(message_text)
        ELSE
          write_text = message_text
        ENDIF
      ELSE
        write_text = message_text
      END IF

    END IF

    ! conditions under which this process actually prints out a message:
    ! - printout explicity requested
    lactive = lprint
    ! - or: this process is stdio process
    lactive = lactive .OR. my_process_is_stdio()
    ! - or: this process is root process in processor splitting
    IF (proc_split) &
     &lactive = lactive .OR. &
       &       (proc_split .AND. (get_my_mpi_work_id() == get_glob_proc0()))

    IF (lactive) THEN
      IF (msg_timestamp) THEN
        CALL DATE_AND_TIME(date=cdate,time=ctime)
        WRITE(iout,'(a)') '['//cdate//' '//ctime//'] '//TRIM(write_text)
        IF (l_log) WRITE(nlog,'(a)') '['//cdate//' '//ctime//']: '//TRIM(write_text)
      ELSE
        WRITE(iout,'(1x,a)') TRIM(write_text)
        IF (l_log) WRITE(nlog,'(1x,a)') TRIM(write_text)
      END IF
    END IF

  END SUBROUTINE message
  !-------------
  

  !-------------
  !>
  !!
  SUBROUTINE print_status (mstring, flag)

    CHARACTER(len=*), intent(in)   :: mstring
    LOGICAL, intent(in)            :: flag

    IF ( flag ) THEN
       WRITE(message_text,'(a60,1x,":",a)') mstring,'active'
    ELSE
       WRITE(message_text,'(a60,1x,":",a)') mstring,'*not* active'
    END IF
    CALL message('', message_text, level=em_param)

  END SUBROUTINE print_status
  !-------------
  !>
  !! Report the value of a logical, integer or real variable
  !! Convenience routine interfaced by print_value(mstring, value)
  !!
  SUBROUTINE print_lvalue (mstring, lvalue)

    CHARACTER(len=*), intent(in)   :: mstring
    LOGICAL, intent(in)            :: lvalue

    IF (lvalue) THEN
      WRITE(message_text,'(a60,1x,": ",a)') mstring,'TRUE'
    ELSE
      WRITE(message_text,'(a60,1x,": ",a)') mstring,'FALSE'
    END IF
    CALL message('', message_text, level=em_param)

  END SUBROUTINE print_lvalue
  !-------------
  !>
  !!
  SUBROUTINE print_ivalue (mstring, ivalue)

    CHARACTER(len=*), intent(in)   :: mstring
    INTEGER, intent(in)            :: ivalue

    WRITE(message_text,'(a60,1x,":",i10)') mstring, ivalue
    CALL message('', message_text, level=em_param)

  END SUBROUTINE print_ivalue
  !-------------
  !>
  !!
  SUBROUTINE print_rvalue (mstring, rvalue)

    CHARACTER(len=*), intent(in)   :: mstring
    REAL(wp), intent(in)           :: rvalue

    WRITE(message_text,'(a60,1x,":",g12.5)') mstring, rvalue
    CALL message('', message_text, level=em_param)

  END SUBROUTINE print_rvalue
  !-------------
  !>
  !!
  SUBROUTINE open_log (logfile_name)

    CHARACTER(len=*), INTENT(in) :: logfile_name
    LOGICAL                      :: l_opened

    INQUIRE (UNIT=nlog,OPENED=l_opened)

    IF (l_opened) THEN
      WRITE (message_text,'(a)') 'log file unit has been used already.'
      CALL message ('open_log', message_text)
      WRITE (message_text,'(a)') 'Close unit and reopen for log file.'
      CALL message ('open_log', message_text, level=em_warn)
      CLOSE (nlog)
    ENDIF

    OPEN (nlog,file=TRIM(logfile_name))

    l_log = .TRUE.

  END SUBROUTINE open_log
  !-------------
  !>
  !!
  SUBROUTINE close_log
    LOGICAL :: l_opened

    INQUIRE (UNIT=nlog,OPENED=l_opened)
    IF (l_opened) THEN
      CLOSE (nlog)
    ENDIF

    l_log = .FALSE.

  END SUBROUTINE close_log
  !-------------

END MODULE mo_exception
