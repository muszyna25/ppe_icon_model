!> Routines for handling of regular output steps and ready file events
!> on multiple I/O PEs.
!!
!! Note that this module contains the abstract implementation of the
!! output event handling, while the actual control of the trigger
!! mechanism is located in the module "mo_output_event_control".
!!
!! Detailed description:
!!
!! 1. The writing of regular output is triggered at so-called "output
!!    event steps". Additional tasks for such an output event step may
!!    be the opening of a new output file or the closing of an output
!!    file. There may be several output events, each with a different
!!    number of output event steps.
!!
!! 2. Output event steps happen at regular intervals. These are given
!!    by an interval size and the time stamps for begin and end.  One
!!    may think of each output event as being handled by its own
!!    PE. The completion of an output event step is communicated via
!!    non-blocking MPI messages to the root PE, which keeps track of
!!    the overall event status.
!!
!! 3. Different output events/PEs may be joined together to form a
!!    single new event. Then not every PE necessarily participates in
!!    every event step. The purpose of joining/grouping events is that
!!    certain actions can be triggered for the group, e.g. writing
!!    "ready files".
!!
!!    A "ready file" is a technique for handling dependencies between
!!    the NWP processes at DWD: When a program - parallel or
!!    sequential, shell script or binary - produces some output which
!!    is necessary for other running applications, then the completion
!!    of the write process signals this by creating a small file
!!    (size: a few bytes). Only when this file exists, the second
!!    program starts reading its input data. Implicity, this assumes
!!    that a file system creates (and closes) files in the same order
!!    as they are written by the program.
!!
!!
!! @note This event handling is *static*, i.e. all event occurrences
!!       are pre-defined during the initialization. This is necessary
!!       for generating a complete list of "ready" files and output
!!       files at the beginning of the simulation.
!!
!!
!! MPI ROLES
!! ---------
!!
!! The output event handling requires communication between the I/O
!! PEs and a root I/O MPI rank "ROOT_OUTEVENT" - usually chosen as
!! rank 0 from the I/O PEs.
!!
!! Root PE: Parallel communication is necessary
!!
!! - During the setup phase:
!!
!!   The root I/O MPI rank asks all participating I/O PEs for their
!!   output event info. This information about these events is
!!   forwarded (broadcasted) to the worker PEs. Afterwards, each of
!!   the I/O and worker PEs generates a unified output event,
!!   indicating which PE performs a write process at which step. Thus,
!!   every PE with the exception of the asynchronous restart PEs has
!!   the same view on the sequence of output events.
!!
!!   see FUNCTION union_of_all_events
!!
!! - During the loop of output steps:
!!
!!   The root I/O rank places MPI_IRECV calls, thereby asking all
!!   participating PEs for acknowledging the completion of the current
!!   output step.
!!
!!   see  SUBROUTINE trigger_output_step_irecv
!!        FUNCTION is_output_step_complete
!!
!! All I/O PEs: Parallel communication is performed
!!
!! - During setup:
!!
!!   Output events are generated and the necessary meta-data for
!!   replication on the root I/O PE is send via MPI to the root
!!   rank. An event name which is different from the string
!!   DEFAULT_EVENT_NAME ( = "default") leads to the creation of ready
!!   files during the output loop.
!!
!!   see FUNCTION new_parallel_output_event
!!
!! - During the loop of output steps:
!!
!!   The I/O PEs acknowledge the completion of their respective output
!!   step by sending non-blocking MPI messages.  The non-blocking
!!   messages have unique MPI tags that are computed as i_tag =
!!   SENDRECV_TAG_OUTEVENT + this_pe + (local_event_no-1)*nranks
!!
!!   see SUBROUTINE pass_output_step
!!
!!
!! @author F. Prill, DWD
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2013-09-17)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! -----------------------------------------------------------------------------------
MODULE mo_output_event_handler

  ! actual method (MPI-2)
#ifndef NOMPI
#if !defined (__SUNPRO_F95)
  USE mpi
#endif
#endif

  USE mo_kind,                   ONLY: i8
  USE mo_impl_constants,         ONLY: SUCCESS, MAX_TIME_INTERVALS
  USE mo_exception,              ONLY: finish
  USE mo_io_units,               ONLY: FILENAME_MAX, find_next_free_unit
  USE mo_util_string,            ONLY: int2string, remove_whitespace
  USE mo_mpi,                    ONLY: p_int,                                               &
    &                                  p_pack_int, p_pack_string, p_pack_bool, p_pack_real, &
    &                                  p_unpack_int, p_unpack_string, p_unpack_bool,        &
    &                                  p_unpack_real, p_send_packed, p_irecv_packed,        &
    &                                  p_wait, p_bcast, get_my_global_mpi_id,               &
    &                                  my_process_is_mpi_test, p_pe,                        &
    &                                  my_process_is_mpi_workroot,                          &
    &                                  p_comm_rank, p_comm_size
  USE mtime,                     ONLY: MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN,         &
    &                                  datetime, timedelta,  newTimedelta,                  &
    &                                  deallocateDatetime, datetimeToString,                &
    &                                  newDatetime, OPERATOR(>=), OPERATOR(==),&
    &                                  OPERATOR(<), OPERATOR(<=), &
    &                                  OPERATOR(>), OPERATOR(+), OPERATOR(/=),              &
    &                                  deallocateTimedelta, newJulianDay, JulianDay,        &
    &                                  getJulianDayFromDatetime, getDatetimeFromJulianDay
  USE mo_output_event_types,     ONLY: t_sim_step_info, t_event_step_data,                  &
    &                                  t_event_step, t_output_event, t_par_output_event,    &
    &                                  MAX_FILENAME_STR_LEN, MAX_EVENT_NAME_STR_LEN,        &
    &                                  DEFAULT_EVENT_NAME
  USE mo_name_list_output_types, ONLY: t_fname_metadata
  USE mo_util_table,             ONLY: initialize_table, finalize_table, add_table_column,  &
    &                                  set_table_entry, print_table, t_table
  USE mo_name_list_output_config,   ONLY: use_async_name_list_io
  IMPLICIT NONE

  ! public subroutines + functions:
  PRIVATE

#ifndef NOMPI
#if defined (__SUNPRO_F95)
  INCLUDE "mpif.h"
#endif
#endif

  ! initialization + destruction
  PUBLIC :: new_parallel_output_event
  PUBLIC :: complete_event_setup
  PUBLIC :: union_of_all_events
  PUBLIC :: deallocate_output_event
  PUBLIC :: blocking_wait_for_irecvs
  ! inquiry functions
  PUBLIC :: get_current_date
  PUBLIC :: get_current_step
  PUBLIC :: get_current_filename
  PUBLIC :: get_current_jfile
  PUBLIC :: check_write_readyfile
  PUBLIC :: check_open_file
  PUBLIC :: check_close_file
  PUBLIC :: is_event_root_pe
  PUBLIC :: is_output_step
  PUBLIC :: is_output_step_complete
  PUBLIC :: is_output_event_finished
  ! handshake routines (called after steps)
  PUBLIC :: pass_output_step
#ifndef NOMPI
  PUBLIC :: trigger_output_step_irecv
#endif
  ! auxiliary functions
  PUBLIC :: print_output_event
  PUBLIC :: set_event_to_simstep
  PUBLIC :: strip_from_modifiers

  INTERFACE deallocate_output_event
    MODULE PROCEDURE deallocate_output_event
    MODULE PROCEDURE deallocate_par_output_event
  END INTERFACE

  INTERFACE is_output_event_finished
    MODULE PROCEDURE is_output_event_finished
    MODULE PROCEDURE is_par_output_event_finished
  END INTERFACE

  INTERFACE is_output_step
    MODULE PROCEDURE is_output_step
    MODULE PROCEDURE is_par_output_step
  END INTERFACE

  INTERFACE get_current_date
    MODULE PROCEDURE get_current_date
    MODULE PROCEDURE get_current_date_par
  END INTERFACE

  INTERFACE get_current_step
    MODULE PROCEDURE get_current_step
    MODULE PROCEDURE get_current_step_par
  END INTERFACE

  INTERFACE print_output_event
    MODULE PROCEDURE print_output_event
    MODULE PROCEDURE print_par_output_event
  END INTERFACE

  INTERFACE set_event_to_simstep
    MODULE PROCEDURE set_event_to_simstep
    MODULE PROCEDURE set_event_to_simstep_par
  END INTERFACE


  !---------------------------------------------------------------
  ! constants

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_output_event_handler'

  !> maximum buffer size for sending event meta-data (MPI_PACK)
  INTEGER, PARAMETER :: MAX_BUF_SIZE =    (    MAX_EVENT_NAME_STR_LEN    &
    &                                      + 9*(MAX_DATETIME_STR_LEN+1)  &
    &                                      + MAX_TIMEDELTA_STR_LEN       &
    &                                      + 3*FILENAME_MAX              &
    &                                      + 1024 )

  !> MPI message tag for setup of output events
  INTEGER, PARAMETER :: SENDRECV_TAG_SETUP    = 1002

  !> MPI message tag for output event handshake
  INTEGER, PARAMETER :: SENDRECV_TAG_OUTEVENT = 1001

  !> MPI rank of root PE handling ready files
  INTEGER, PARAMETER :: ROOT_OUTEVENT         = 0

  !> Internal switch for debugging output
  LOGICAL, PARAMETER :: ldebug                = .FALSE.

  !> Max. no. of steps printed out to stderr
  INTEGER, PARAMETER :: MAX_PRINTOUT          = 2000


  !---------------------------------------------------------------
  ! local list event with event meta-data
  !---------------------------------------------------------------

  !> event meta-data: this is only required for I/O PE #0, where we have to
  !  keep a local list of output events.
  TYPE t_event_data_local
    CHARACTER(LEN=MAX_EVENT_NAME_STR_LEN)        :: name                 !< output event name
    CHARACTER(LEN=MAX_DATETIME_STR_LEN+1)        :: begin_str(MAX_TIME_INTERVALS) !< date-time stamp + modifier
    CHARACTER(LEN=MAX_DATETIME_STR_LEN+1)        :: end_str(MAX_TIME_INTERVALS)   !< date-time stamp + modifier
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)         :: intvl_str(MAX_TIME_INTERVALS)
    !> Flag. If .TRUE. the last step is always written
    LOGICAL                                      :: l_output_last
    !> definitions for conversion "time stamp -> simulation step"
    TYPE(t_sim_step_info)                        :: sim_step_info
    !> additional meta-data for generating output filename
    TYPE(t_fname_metadata)                       :: fname_metadata
    !> this event's MPI tag
    INTEGER                                      :: i_tag
  END TYPE t_event_data_local

  !> Maximum length of local event meta-data list
  INTEGER, PARAMETER :: LOCAL_NMAX_EVENT_LIST = 100

  !> local list of output events
  TYPE(t_event_data_local) :: event_list_local(LOCAL_NMAX_EVENT_LIST)

  !> length of local list of output events
  INTEGER :: ievent_list_local = 0


CONTAINS

  !---------------------------------------------------------------
  ! ROUTINES FOR DEALLOCATION
  !---------------------------------------------------------------

  !> Deallocate t_event_step data structure.
  !  @author F. Prill, DWD
  !
  SUBROUTINE deallocate_event_step(event_step)
    TYPE(t_event_step), INTENT(INOUT) :: event_step
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::deallocate_event_step"
    INTEGER :: ierrstat
    DEALLOCATE(event_step%event_step_data, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE deallocate_event_step


  !> Deallocate t_output_event data structure.
  ! @author F. Prill, DWD
  !
  SUBROUTINE deallocate_output_event(event)
    TYPE(t_output_event), POINTER :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::deallocate_output_event"
    INTEGER :: ierrstat,i
    IF (.NOT. ASSOCIATED(event)) RETURN
    DO i=1,event%n_event_steps
      CALL deallocate_event_step(event%event_step(i))
    END DO
    event%n_event_steps = 0
    DEALLOCATE(event%event_step, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    DEALLOCATE(event, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE deallocate_output_event


  !> Deallocate t_par_output_event data structure.
  !
  !  @note This subroutine recursively deallocates the complete
  !        singly-linked list rooted at @p event.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE deallocate_par_output_event(event)
    TYPE(t_par_output_event), POINTER :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::deallocate_par_output_event"
    INTEGER :: ierrstat

    TYPE(t_par_output_event), POINTER :: ev, next

    ierrstat = success
    ev => event
    DO WHILE (ASSOCIATED(ev))
      next => ev%next
      CALL deallocate_output_event(ev%output_event)
      IF (ALLOCATED(ev%irecv_req)) THEN
        DEALLOCATE(ev%output_event, ev%irecv_req, ev%irecv_buf, STAT=ierrstat)
      END IF
      IF (ierrstat == success) DEALLOCATE(ev, STAT=ierrstat)
      IF (ierrstat /= success) CALL finish (routine, 'DEALLOCATE failed.')
      ev => next
    END DO
    NULLIFY(event)
  END SUBROUTINE deallocate_par_output_event


  !---------------------------------------------------------------
  ! ROUTINES FOR PRINTING EVENT DESCRIPTIONS
  !---------------------------------------------------------------

  !> Screen print-out of an output event.
  !  @author F. Prill, DWD
  !
  !  @note If an (optional) destination file is provided, we
  !        implicitly assume that this file has already been opened by
  !        the called and will be closed by the caller.
  !
  SUBROUTINE print_output_event(event, opt_dstfile)
    TYPE(t_output_event), POINTER             :: event           !< output event data structure
    INTEGER, OPTIONAL,             INTENT(IN) :: opt_dstfile     !< optional destination ASCII file unit
    ! local variables
    INTEGER                     :: i, j, irow, dst
    TYPE(t_table)               :: table
    TYPE(t_event_step), POINTER :: event_step

    dst = 0
    IF (PRESENT(opt_dstfile)) dst = opt_dstfile

    IF (ldebug) THEN
      WRITE (0,*) "PE ", get_my_global_mpi_id(), ": event%n_event_steps = ", event%n_event_steps
    END IF

    IF (check_write_readyfile(event)) THEN
      WRITE (dst,'(3a)') 'output "', TRIM(event%event_data%name), '", writes ready files:'
    ELSE
      WRITE (dst,'(3a)') 'output "', TRIM(event%event_data%name), '", does not write ready files:'
    END IF
#ifdef __SX__
    IF (dst == 0) THEN
      WRITE (dst,'(a)') "output on SX9 has been shortened, cf. 'output_schedule.txt'."
    END IF
#endif

    ! do not print, if the table is excessively long
    IF (event%n_event_steps > MAX_PRINTOUT) THEN
      WRITE (dst,*) "detailed print-out of output steps has been omitted, cf. 'output_schedule.txt'."
      RETURN
    END IF

    ! classic output
    !    DO i=1,event%n_event_steps
    !      CALL print_event_step(event%event_step(i))
    !    END DO

    ! table-based output
    CALL initialize_table(table)
    CALL add_table_column(table, "model step")
    CALL add_table_column(table, "model date")
    CALL add_table_column(table, "filename")
    CALL add_table_column(table, "I/O PE")
    CALL add_table_column(table, "output date")
#ifdef __SX__
    IF (dst /= 0) THEN
#endif
      ! do not add the file-part column on the NEC SX9, because we have a
      ! line limit of 132 characters there
      CALL add_table_column(table, "#")
#ifdef __SX__
    END IF
#endif
    CALL add_table_column(table, "open")
    CALL add_table_column(table, "close")
    irow = 0
    DO i=1,event%n_event_steps
      event_step => event%event_step(i)
      DO j=1,event_step%n_pes
        irow = irow + 1
        IF (j==1) THEN
          CALL set_table_entry(table,irow,"model step", int2string(event_step%i_sim_step))
          CALL set_table_entry(table,irow,"model date", TRIM(event_step%exact_date_string))
        ELSE
          CALL set_table_entry(table,irow,"model step", " ")
          CALL set_table_entry(table,irow,"model date", " ")
        END IF
#ifdef __SX__
        ! save some characters on SX:
        IF ((LEN_TRIM(event_step%event_step_data(j)%filename_string) > 20) .AND. (dst == 0)) THEN
          CALL set_table_entry(table,irow,"filename",    TRIM(event_step%event_step_data(j)%filename_string(1:20)//"..."))
        ELSE
          CALL set_table_entry(table,irow,"filename",    TRIM(event_step%event_step_data(j)%filename_string))
        END IF
#else
        CALL set_table_entry(table,irow,"filename",    TRIM(event_step%event_step_data(j)%filename_string))
#endif
        CALL set_table_entry(table,irow,"I/O PE",      int2string(event_step%event_step_data(j)%i_pe))

        CALL set_table_entry(table,irow,"output date", TRIM(event_step%event_step_data(j)%datetime_string))
#ifdef __SX__
        IF (dst /= 0) THEN
#endif
          ! do not add the file-part column on the NEC SX9, because we have a
          ! line limit of 132 characters there
          CALL set_table_entry(table,irow,"#",           &
            & TRIM(int2string(event_step%event_step_data(j)%jfile))//"."//TRIM(int2string(event_step%event_step_data(j)%jpart)))
#ifdef __SX__
        END IF
#endif
        ! append "+ open" or "+ close" according to event step data:
        IF (event_step%event_step_data(j)%l_open_file) THEN
          CALL set_table_entry(table,irow,"open", "x")
        ELSE
          CALL set_table_entry(table,irow,"open", " ")
        END IF
        IF (event_step%event_step_data(j)%l_close_file) THEN
          CALL set_table_entry(table,irow,"close","x")
        ELSE
          CALL set_table_entry(table,irow,"close"," ")
        END IF
      END DO
    END DO
    CALL print_table(table, opt_delimiter='   ', opt_dstfile=dst)
    CALL finalize_table(table)
  END SUBROUTINE print_output_event


  !> Screen print-out of a parallel output event.
  !
  !  @note This subroutine recursively prints the complete
  !        singly-linked list rooted at @p event.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE print_par_output_event(event, opt_filename, opt_dstfile)
    TYPE(t_par_output_event), TARGET              :: event
    CHARACTER(LEN=*), OPTIONAL,        INTENT(IN) :: opt_filename    !< name of ASCII file (optional)
    INTEGER,          OPTIONAL,        INTENT(IN) :: opt_dstfile     !< optional destination ASCII file unit
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::print_par_output_event"
    INTEGER :: dst, ierrstat
    TYPE(t_par_output_event), POINTER :: ev

    ! consistency check:
    IF (PRESENT(opt_dstfile) .AND. PRESENT(opt_filename)) THEN
      CALL finish(routine, "Routine was called with both a file unit and a filename to open!")
    END IF

    ! do not print, if the table is excessively long
    dst = 0
    IF (event%output_event%n_event_steps > MAX_PRINTOUT) THEN
      WRITE (dst,*) "detailed print-out of output steps has been omitted, cf. 'output_schedule.txt'."
      RETURN
    END IF

    ! open ASCII output file (if necessary):
    IF (PRESENT(opt_dstfile)) THEN
      dst = opt_dstfile
    ELSE IF (PRESENT(opt_filename)) THEN
      dst = find_next_free_unit(10,100)
      OPEN (dst, file=TRIM(opt_filename), status='REPLACE', iostat=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'OPEN failed.')
    END IF

    ev => event
    DO WHILE (ASSOCIATED(ev))
      WRITE (dst,'(a)') "" ! newline
      CALL print_output_event(ev%output_event, opt_dstfile=dst)
      ev => ev%next
    END DO
    WRITE (dst,'(a)') "" ! newline

    ! close ASCII output file (if necessary):
    IF (PRESENT(opt_filename)) THEN
      CLOSE (dst, iostat=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'CLOSE failed.')
    END IF
  END SUBROUTINE print_par_output_event


  !> Screen print-out of a single output event step.
  !  @author F. Prill, DWD
  !
  SUBROUTINE print_event_step(event_step)
    TYPE(t_event_step), INTENT(IN) :: event_step
    ! local variables
    CHARACTER(LEN=7), PARAMETER :: SUFFIX_STR(3) = (/ "+ open ", "+ close", "       "/)
    INTEGER :: i, suffix1, suffix2

    WRITE (0,'(a,i0,a,a)') "   model step ", event_step%i_sim_step,      &
      &                    ", exact date: ", TRIM(event_step%exact_date_string)
    DO i=1,event_step%n_pes
      suffix1 = 3
      suffix2 = 3
      ! append "+ open" or "+ close" according to event step data:
      IF (event_step%event_step_data(i)%l_open_file)  suffix1 = 1
      IF (event_step%event_step_data(i)%l_close_file) suffix2 = 2
      WRITE (0,'(a,a,a,i0,a,a,a,a,a)') &
        & "      output to '",                                       &
        & TRIM(event_step%event_step_data(i)%filename_string),       &
        & "' (I/O PE ", event_step%event_step_data(i)%i_pe,          &
        & ", ", TRIM(event_step%event_step_data(i)%datetime_string), &
        & ") ", SUFFIX_STR(suffix1), SUFFIX_STR(suffix2)
    END DO
  END SUBROUTINE print_event_step


  !---------------------------------------------------------------
  ! ROUTINES FOR CREATING NEW OUTPUT EVENTS
  !---------------------------------------------------------------

  !> Create a simple output event, happening at regular intervals.
  !
  !  These are given by an interval size @p intvl_str and the time
  !  stamps for begin and end, @p begin_str and @p end_str.
  !
  !  The begin and end time stamps may contain "modifier symbols",
  !  e.g. ">2014-06-01T00:00:00.000", where ">" means that output
  !  should not include the start date.
  !
  !  Note that this subroutine generates the output event steps but
  !  does not map these time stamps onto the corresponding simulation
  !  steps. This simulation-specific task is performed by a different
  !  subroutine, which is given as a function parameter
  !  "fct_time2simstep".
  !
  !  @author F. Prill, DWD
  !
  FUNCTION new_output_event(evd, i_pe, fct_time2simstep,    &
    &                       fct_generate_filenames) RESULT(p_event)
    TYPE(t_output_event),                POINTER :: p_event
    TYPE(t_event_data_local),            INTENT(IN)  :: evd
    INTEGER,                             INTENT(IN)  :: i_pe                 !< rank of participating PE

    INTERFACE
      !> As an argument of this function, the user must provide a
      !!  conversion "time stamp -> simulation step"
      SUBROUTINE fct_time2simstep(nstrings, date_string, sim_step_info, &
        &                         result_steps, result_exactdate)
        USE mo_output_event_types, ONLY: t_sim_step_info
        INTEGER,               INTENT(IN)    :: nstrings             !< no. of string to convert
        CHARACTER(len=*),      INTENT(IN)    :: date_string(:)       !< array of ISO 8601 time stamp strings
        TYPE(t_sim_step_info), INTENT(IN)    :: sim_step_info        !< definitions: time step size, etc.
        INTEGER,               INTENT(INOUT) :: result_steps(:)      !< resulting step indices
        CHARACTER(LEN=*),      INTENT(INOUT) :: result_exactdate(:)  !< resulting (exact) time step strings
      END SUBROUTINE fct_time2simstep

      !> As an argument of this function, the user must provide a
      !! function for generating output file names
      SUBROUTINE fct_generate_filenames(nstrings, date_string, sim_steps, &
        &                             sim_step_info, fname_metadata, &
        &                             skipped_dates, result_fnames)
        USE mo_output_event_types,     ONLY: t_sim_step_info, t_event_step_data
        USE mo_name_list_output_types, ONLY: t_fname_metadata

        INTEGER,                 INTENT(IN)    :: nstrings           !< no. of string to convert
        CHARACTER(len=*),        INTENT(IN)    :: date_string(:)     !< array of ISO 8601 time stamp strings
        INTEGER,                 INTENT(IN)    :: sim_steps(:)       !< array of corresponding simulation steps
        TYPE(t_sim_step_info),   INTENT(IN)    :: sim_step_info      !< definitions: time step size, etc.
        TYPE(t_fname_metadata),  INTENT(IN)    :: fname_metadata     !< additional meta-data for generating output filename
        INTEGER,                 INTENT(IN)    :: skipped_dates
        TYPE(t_event_step_data), INTENT(out)   :: result_fnames(SIZE(date_string))
      END SUBROUTINE fct_generate_filenames
    END INTERFACE

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::new_output_event"
    !> Max. no. of event steps (used for local array sizes)
    INTEGER, PARAMETER :: INITIAL_NEVENT_STEPS = 256 ! tested to be a good compromise

    TYPE(datetime),  POINTER :: mtime_date, mtime_begin, mtime_end, mtime_restart, &
      &                         sim_end, mtime_dom_start, mtime_dom_end, run_start
    TYPE(timedelta), POINTER :: delta
    INTEGER                  :: ierrstat, i, n_event_steps, &
      &                         nintvls, iintvl, skipped_dates, &
      &                         old_size, new_size
    LOGICAL                  :: l_active, l_append_step
    CHARACTER(len=MAX_DATETIME_STR_LEN), ALLOCATABLE :: mtime_date_string(:), tmp_string(:)
    INTEGER,                             ALLOCATABLE :: mtime_sim_steps(:)
    CHARACTER(len=MAX_DATETIME_STR_LEN), ALLOCATABLE :: mtime_exactdate(:)
    TYPE(t_event_step_data),             ALLOCATABLE :: filename_metadata(:)
    TYPE(t_event_step_data),             POINTER     :: step_data
    CHARACTER(len=MAX_DATETIME_STR_LEN+1)            :: dt_string
    CHARACTER(len=MAX_DATETIME_STR_LEN)              :: begin_str2(MAX_TIME_INTERVALS)
    CHARACTER(len=MAX_DATETIME_STR_LEN)              :: end_str2(MAX_TIME_INTERVALS)
    LOGICAL                                          :: incl_begin(MAX_TIME_INTERVALS)
    LOGICAL                                          :: incl_end(MAX_TIME_INTERVALS)
    CHARACTER                                        :: char
    TYPE(julianday), POINTER :: tmp_jd => NULL()
    TYPE tmp_container
      INTEGER(i8) :: day
      INTEGER(i8) :: ms
    END TYPE tmp_container

    TYPE(tmp_container), ALLOCATABLE, TARGET :: mtime_date_container_a(:)
    TYPE(tmp_container), ALLOCATABLE, TARGET :: mtime_date_container_b(:)  
    TYPE(tmp_container), ALLOCATABLE :: tmp(:)
    TYPE(tmp_container), ALLOCATABLE :: mtime_date_uniq(:)
    
    TYPE(tmp_container), POINTER :: mtime_date_container(:) => NULL()
    
    INTEGER, ALLOCATABLE :: indices_to_use(:)
    INTEGER :: remaining_intvls, iselected_intvl
    
    INTEGER :: n_event_steps_a, n_event_steps_b, remaining_event_steps
    
    ! allocate event data structure
    ALLOCATE(p_event, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    n_event_steps                    = 0
    p_event%i_event_step             = 1
    p_event%event_data%name          = evd%name
    p_event%event_data%sim_start     = evd%sim_step_info%sim_start

    ! count the number of different time intervals for this event (usually 1)
    DO nintvls = 0, max_time_intervals-1
      IF (LEN_TRIM(evd%begin_str(nintvls+1)) == 0) EXIT
    END DO

    ! status output
    IF (ldebug) THEN
      WRITE (0,*) "PE ",get_my_global_mpi_id(), ":"
      WRITE (0,'(7a)') 'Defining output event "', &
           TRIM(evd%begin_str(1)), '", "', &
           TRIM(evd%end_str(1)), '", "', &
           TRIM(evd%intvl_str(1)), '"'
      DO i=2,nintvls
        WRITE (0,'(7a)') ' +  "', &
             TRIM(evd%begin_str(i)), '", "', &
             TRIM(evd%end_str(i)), '", "',   &
             TRIM(evd%intvl_str(i)), '"'
      END DO
      WRITE (0,*) 'Simulation bounds:    "'//TRIM(evd%sim_step_info%sim_start)//'", "'//TRIM(evd%sim_step_info%sim_end)//'"'
      WRITE (0,*) 'restart bound: "'//TRIM(evd%sim_step_info%restart_time)
    END IF

    ! set some dates used later:
    sim_end     => newDatetime(TRIM(evd%sim_step_info%sim_end))
    run_start   => newDatetime(TRIM(evd%sim_step_info%run_start))

    ! Domains (and their output) can be activated and deactivated
    ! during the simulation. This is determined by the parameters
    ! "dom_start_time" and "dom_end_time". Therefore, we must create
    ! a corresponding event.
    mtime_dom_start => newDatetime(TRIM(evd%sim_step_info%dom_start_time))
    mtime_dom_end   => newDatetime(TRIM(evd%sim_step_info%dom_end_time)) ! this carries the domain-specific namelist value "end_time"

    ! To avoid further case discriminations, sim_end is set to the minimum of the simulation end time
    ! and the time at which a nested domain is turned off
    IF (sim_end > mtime_dom_end) THEN
      sim_end => newDatetime(TRIM(evd%sim_step_info%dom_end_time))
    ENDIF

    ! Compute the end time wrt. "dt_restart": It might be that the
    ! simulation end is limited by this parameter
    mtime_restart   => newDatetime(TRIM(evd%sim_step_info%restart_time))

    ! loop over the event occurrences
    
    ALLOCATE(mtime_date_container_a(256), stat=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ALLOCATE(mtime_date_container_b(256), stat=ierrstat)   
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! The begin and end time stamps may contain "modifier symbols",
    ! e.g. ">2014-06-01T00:00:00.000", where ">" means that output
    ! should not include the start date.
    !
    ! At this point, we strip these modifiers from the begin and end
    ! time stamps and store this information in logical flags.
    incl_begin(:) = .TRUE.
    incl_end(:)   = .TRUE.
    DO iintvl=1,nintvls
      ! begin time stamp
      dt_string = evd%begin_str(iintvl)
      char      = dt_string(1:1)
      SELECT CASE(char)
      CASE ('>')
        incl_begin(iintvl) = .FALSE.
        dt_string = dt_string(2:)
      END SELECT
      begin_str2(iintvl) = dt_string

      ! end time stamp
      dt_string = evd%end_str(iintvl)
      char      = dt_string(1:1)
      SELECT CASE(char)
      CASE ('<')
        incl_end(iintvl) = .FALSE.
        dt_string = dt_string(2:)
      END SELECT
      end_str2(iintvl) = dt_string
    END DO

    allocate(indices_to_use(nintvls))
    
    CALL remove_duplicate_intervals(begin_str2, end_str2, evd%intvl_str, nintvls, &
         &                          indices_to_use, remaining_intvls)
    
    ! there may be multiple starts/ends/intervals (usually only one):

    mtime_date_container => mtime_date_container_a
    
    tmp_jd => newJulianday(0_i8, 0_i8);
    
    n_event_steps = 0
    skipped_dates = 0

    INTERVAL_LOOP: DO iselected_intvl = 1, remaining_intvls    
      iintvl = indices_to_use(iselected_intvl)
      
      mtime_begin => newDatetime(TRIM(begin_str2(iintvl)))
      IF (.NOT. ASSOCIATED(mtime_begin))  THEN
        CALL finish(routine, "date-time conversion error: "//TRIM(begin_str2(iintvl)))
      END IF
      mtime_end   => newDatetime(TRIM(end_str2(iintvl)))
      IF (.NOT. ASSOCIATED(mtime_end))  THEN
        CALL finish(routine, "date-time conversion error: "//TRIM(end_str2(iintvl)))
      END IF

      mtime_date  => newDatetime(mtime_begin)
      delta       => newTimedelta(TRIM(evd%intvl_str(iintvl))) ! create a time delta
      IF (mtime_end >= mtime_begin) THEN
        EVENT_LOOP: DO
          IF ((mtime_date    >= run_start)                         .AND. &
            & (sim_end       >= mtime_date)                        .AND. &
            & (mtime_restart >= mtime_date)                        .AND. &
            & (incl_begin(iintvl) .OR. (mtime_date > mtime_begin)) .AND. &
            & (incl_end(iintvl)   .OR. (mtime_end  > mtime_date)) )  THEN

            IF  (mtime_date >= mtime_dom_start) THEN

              n_event_steps = n_event_steps + 1

              old_size = SIZE(mtime_date_container)
              IF (n_event_steps > old_size) THEN
                ALLOCATE(tmp(2*old_size), stat=ierrstat)
                IF (ierrstat /= 0) STOP 'allocate failed'
                tmp(1:old_size) = mtime_date_container(:)
                IF (ASSOCIATED(mtime_date_container, mtime_date_container_a)) THEN
                  CALL MOVE_ALLOC(tmp, mtime_date_container_a)
                  mtime_date_container => mtime_date_container_a
                ELSE
                  CALL MOVE_ALLOC(tmp, mtime_date_container_b)
                  mtime_date_container => mtime_date_container_b
                ENDIF
              ENDIF

              CALL getJulianDayFromDatetime(mtime_date, tmp_jd)
              mtime_date_container(n_event_steps)%day = tmp_jd%day
              mtime_date_container(n_event_steps)%ms = tmp_jd%ms

            ELSE
              ! we skip an output date when the domain is not yet
              ! "active" - this leads to a file count offset.
              skipped_dates = skipped_dates+1
            END IF
          END IF

          IF (ldebug)  WRITE (0,*) get_my_global_mpi_id(), ": adding time delta."
          mtime_date = mtime_date + delta

          l_active = mtime_date <= mtime_end .AND.   &
            &        mtime_date <= sim_end   .AND.   &
            &        mtime_date <= mtime_restart
          IF (.NOT. l_active) EXIT EVENT_LOOP
        END DO EVENT_LOOP
      END IF

      IF (iintvl == 1) THEN
        n_event_steps_a = n_event_steps
        mtime_date_container => mtime_date_container_b
        n_event_steps = 0
        CYCLE interval_loop
      ELSE
        n_event_steps_b = n_event_steps
        n_event_steps = 0      
        
        IF (n_event_steps_b > 0) THEN
          CALL merge2SortedAndRemoveDuplicates(mtime_date_container_a, n_event_steps_a, &
            &                               mtime_date_container_b, n_event_steps_b, &
            &                               mtime_date_uniq, remaining_event_steps)

          IF (remaining_event_steps > SIZE(mtime_date_container_a)) THEN
            ALLOCATE(tmp(SIZE(mtime_date_container_a)), stat=ierrstat)
            IF (ierrstat /= 0) STOP 'allocate failed'
            tmp(1:remaining_event_steps) = mtime_date_uniq(1:remaining_event_steps)
            ! FIXME: this is probably a bug and the deallocation of
            ! mtime_date_container_a needs to be framed with
            ! 1.
            ! l_assoc = ASSOCIATED(mtime_date_container, mtime_date_container_a)
            ! ...
            CALL MOVE_ALLOC(tmp, mtime_date_container_a)
            ! ... and
            ! 2.
            ! mtime_date_container => mtime_date_container_a
          ELSE
            mtime_date_container_a(1:remaining_event_steps) = mtime_date_uniq(1:remaining_event_steps)
          ENDIF
          
          n_event_steps_a = remaining_event_steps
          
          DEALLOCATE(mtime_date_uniq)
        END IF
      ENDIF
      
    END DO INTERVAL_LOOP

    ! copy back results into original data structures
    
    n_event_steps = n_event_steps_a
    ! to prevent a potential reallocation in next step add 1 element
    ALLOCATE(mtime_date_string(n_event_steps+1), STAT=ierrstat)    
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! copy mtime_date_container_a to mtime_date_string 
    DO i = 1, n_event_steps
      tmp_jd => newJulianday(mtime_date_container_a(i)%day, mtime_date_container_a(i)%ms);      
      CALL getDatetimeFromJulianDay(tmp_jd, mtime_date)
      CALL datetimeToString(mtime_date, mtime_date_string(i))
    END DO

    ! Optional: Append the last event time step
    IF (evd%l_output_last .AND. mtime_date < sim_end .AND. mtime_end >= sim_end) THEN
      ! check, that we do not duplicate the last time step:
      l_append_step = .TRUE.
      IF (n_event_steps > 0) THEN
        mtime_date => newDatetime(TRIM(mtime_date_string(n_event_steps)))
        IF (mtime_date == sim_end)  l_append_step = .FALSE.
      END IF
      IF (l_append_step) THEN 
        n_event_steps = n_event_steps + 1
        old_size = SIZE(mtime_date_string)
        IF (n_event_steps > old_size) THEN
          ! resize buffer
          new_size = old_size + INITIAL_NEVENT_STEPS
          ALLOCATE(tmp_string(new_size), STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
          tmp_string(1:old_size) = mtime_date_string
          CALL MOVE_ALLOC(tmp_string, mtime_date_string)
        END IF
        CALL datetimeToString(sim_end, mtime_date_string(n_event_steps))
        IF (ldebug) THEN
          WRITE (0,*) get_my_global_mpi_id(), ": ", &
               &      n_event_steps, ": output event '", mtime_date_string(n_event_steps), "'"
        END IF
      END IF
    END IF
    
    ALLOCATE(mtime_sim_steps(SIZE(mtime_date_string)),   &
         &   mtime_exactdate(SIZE(mtime_date_string)),   &
         &   filename_metadata(SIZE(mtime_date_string)), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    
    CALL fct_time2simstep(n_event_steps, mtime_date_string, &
      &                   evd%sim_step_info, mtime_sim_steps, mtime_exactdate)

    ! remove all those event steps which have no corresponding simulation
    ! step (mtime_sim_steps(i) < 0):
    DO i=1,n_event_steps
      IF (mtime_sim_steps(i) < 0)  EXIT
    END DO
    n_event_steps = (i-1)
    
    IF (n_event_steps > 0) THEN
      CALL fct_generate_filenames(n_event_steps, mtime_date_string,       &
        &                   mtime_sim_steps, evd%sim_step_info, evd%fname_metadata, &
        &                   skipped_dates, filename_metadata)
    END IF
    
    ! from this list of time stamp strings: generate the event steps
    ! for this event
    p_event%n_event_steps = n_event_steps
    ALLOCATE(p_event%event_step(n_event_steps), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    DO i=1,n_event_steps
      p_event%event_step(i)%exact_date_string  = mtime_exactdate(i)
      p_event%event_step(i)%i_sim_step         = mtime_sim_steps(i)
      p_event%event_step(i)%n_pes              = 1
      ALLOCATE(p_event%event_step(i)%event_step_data(p_event%event_step(i)%n_pes), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      step_data => p_event%event_step(i)%event_step_data(1)
      step_data%i_pe            = i_pe
      step_data%datetime_string = mtime_date_string(i)
      step_data%i_tag           = evd%i_tag
      step_data%filename_string = filename_metadata(i)%filename_string
      step_data%jfile           = filename_metadata(i)%jfile
      step_data%jpart           = filename_metadata(i)%jpart
      step_data%l_open_file     = filename_metadata(i)%l_open_file
      step_data%l_close_file    = filename_metadata(i)%l_close_file
    END DO
    IF (ldebug) THEN
      WRITE (0,*) routine, ": defined event ",                            &
        &         TRIM(p_event%event_step(1)%event_step_data(1)%filename_string), "; tag = ", &
        &         p_event%event_step(1)%event_step_data(1)%i_tag
    END IF

    ! clean up
    CALL deallocateDatetime(mtime_begin)
    CALL deallocateDatetime(mtime_end)
    CALL deallocateDatetime(mtime_date)
    CALL deallocateDatetime(mtime_dom_start)
    CALL deallocateDatetime(mtime_dom_end)
    CALL deallocateDatetime(mtime_restart)
    CALL deallocateDatetime(sim_end)
    CALL deallocateDatetime(run_start)
    CALL deallocateTimedelta(delta)

    DEALLOCATE(mtime_date_string, mtime_sim_steps, &
      &        mtime_exactdate, filename_metadata, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  CONTAINS

    SUBROUTINE merge2SortedAndRemoveDuplicates(InputArray1, nsize_IA1, &
         &                                     InputArray2, nsize_IA2, &
         &                                     OutputArray, nsize_OA)
      TYPE(tmp_container), INTENT(in) :: InputArray1(:)
      TYPE(tmp_container), INTENT(in) :: InputArray2(:)
      TYPE(tmp_container), ALLOCATABLE, INTENT(out) :: OutputArray(:)
      INTEGER, INTENT(in) :: nsize_IA1
      INTEGER, INTENT(in) :: nsize_IA2
      INTEGER, INTENT(out) :: nsize_OA
      
      INTEGER(i8) :: diff
      
      INTEGER :: n, na, nb
      INTEGER :: i, j, k


      IF (nsize_IA1 == 0 .AND. nsize_IA2 == 0) THEN
        CALL finish('merge2SortedAndRemoveDuplicates', 'Invalid size of both input arrays.')
      ENDIF

      IF (nsize_IA1 == 0) THEN
        nsize_OA = nsize_IA2
        ALLOCATE(OutputArray(nsize_IA2))
        OutputArray(:) = InputArray2(:nsize_IA2)
        RETURN
      ENDIF
      IF (nsize_IA2 == 0) THEN
        nsize_OA = nsize_IA1
        ALLOCATE(OutputArray(nsize_IA1))
        OutputArray(:) = InputArray1(:nsize_IA1)
        RETURN
      ENDIF
      
      na = nsize_IA1
      nb = nsize_IA2 
      n = na + nb
      
      ALLOCATE(OutputArray(n))
      
      i = 1
      j = 1
      k = 1

      ! handle special case for k == 1 to be Fortran conforming in start-up step

      diff = 86400000_i8 * (InputArray1(i)%day - InputArray2(j)%day) + InputArray1(i)%ms - InputArray2(j)%ms
      IF (diff < 0_i8) THEN
        OutputArray(1) = InputArray1(1)
        i = 2
        k = 2
      ELSE IF (diff > 0_i8) THEN
        OutputArray(1) = InputArray2(1)
        j = 2
        k = 2
      ELSE
        OutputArray(1) = InputArray1(1)
        i = 2
        j = 2
        k = 2
      ENDIF
      
      DO WHILE (i <= na .AND. j <= nb)
        diff = 86400000_i8 * (InputArray1(i)%day - InputArray2(j)%day) + InputArray1(i)%ms - InputArray2(j)%ms
        IF (diff < 0_i8) THEN
          IF ((InputArray1(i)%day /= OutputArray(k-1)%day) .OR. (InputArray1(i)%ms /= OutputArray(k-1)%ms)) THEN
            OutputArray(k) = InputArray1(i)
            k = k+1
          ENDIF
          i = i+1
        ELSE IF (diff > 0_i8) THEN
          IF ((InputArray2(j)%day /= OutputArray(k-1)%day) .OR. (InputArray2(j)%ms /= OutputArray(k-1)%ms)) THEN
            OutputArray(k) = InputArray2(j)
            k = k+1
          ENDIF
          j = j+1
        ELSE
          IF ((InputArray1(i)%day /= OutputArray(k-1)%day) .OR. (InputArray1(i)%ms /= OutputArray(k-1)%ms)) THEN
            OutputArray(k) = InputArray1(i)
            k = k+1
          ENDIF
          i = i+1
          j = j+1
        ENDIF
      ENDDO
      
      DO WHILE (i <= na)
        IF ((InputArray1(i)%day /= OutputArray(k-1)%day) .OR. (InputArray1(i)%ms /= OutputArray(k-1)%ms)) THEN
          OutputArray(k) = InputArray1(i)
          k = k+1
          i = i+1
        ELSE
          i = i+1
        ENDIF
      ENDDO
      
      DO WHILE (j <= nb)
        IF ((InputArray2(j)%day /= OutputArray(k-1)%day) .OR. (InputArray2(j)%ms /= OutputArray(k-1)%ms)) THEN
          OutputArray(k) = InputArray2(j)
          k = k+1
          j = j+1
        ELSE
          j = j+1
        ENDIF
      ENDDO
      
      nsize_OA = k-1
      
    END SUBROUTINE merge2SortedAndRemoveDuplicates

    SUBROUTINE remove_duplicate_intervals(starts, ends, intvls, n, indices_to_use, remaining)
      CHARACTER(len=*), INTENT(in) :: starts(:)
      CHARACTER(len=*), INTENT(in) :: ends(:)
      CHARACTER(len=*), INTENT(in) :: intvls(:)
      INTEGER, INTENT(in) :: n         
      INTEGER, INTENT(out) :: indices_to_use(n)
      INTEGER, INTENT(out) :: remaining
      
      INTEGER :: i, j
      
      remaining = 1
      indices_to_use(1) = 1
      
      OUTER_LOOP: DO i = 2, n
        DO j = 1, remaining
          IF (TRIM(starts(j)) == TRIM(starts(i))) THEN
            IF (TRIM(ends(j)) == TRIM(ends(i))) THEN
              IF (TRIM(intvls(j)) == TRIM(intvls(i))) THEN
                ! found a match so start looking again
                CYCLE OUTER_LOOP
              END IF
            END IF
          END IF
        END DO
        ! no match found so add it to the output
        remaining = remaining + 1
        indices_to_use(remaining) = i
      END DO OUTER_LOOP
      
    END SUBROUTINE remove_duplicate_intervals
    
  END FUNCTION new_output_event


  !> Create a simple *parallel* output event, happening at regular intervals.
  !
  !  This subroutine calls the local version "new_output_event", adds
  !  data structures for parallel communication and launches a
  !  non-blocking MPI send to the root I/O PE.
  !
  !  The begin and end time stamps may contain "modifier symbols",
  !  e.g. ">2014-06-01T00:00:00.000", where ">" means that output
  !  should not include the start date.
  !
  !  @author F. Prill, DWD
  !
  FUNCTION new_parallel_output_event(name, begin_str, end_str, intvl_str,                               &
    &                                l_output_last, sim_step_info, fname_metadata, fct_time2simstep,    &
    &                                fct_generate_filenames, local_event_no, icomm) RESULT(p_event)
    TYPE(t_par_output_event), POINTER :: p_event
    CHARACTER(LEN=*),                      INTENT(IN)  :: name                 !< output event name
    !> start time stamp + modifier
    CHARACTER(len=MAX_DATETIME_STR_LEN+1), INTENT(IN)  :: begin_str(MAX_TIME_INTERVALS)
    !> start time stamp + modifier
    CHARACTER(len=MAX_DATETIME_STR_LEN+1), INTENT(IN)  :: end_str(MAX_TIME_INTERVALS)
    CHARACTER(len=MAX_DATETIME_STR_LEN),   INTENT(IN)  :: intvl_str(MAX_TIME_INTERVALS)
    !> Flag. If .TRUE. the last step is always written
    LOGICAL,                               INTENT(IN)  :: l_output_last
    !> definitions for conversion "time stamp -> simulation step"
    TYPE(t_sim_step_info),                 INTENT(IN)  :: sim_step_info
    !> additional meta-data for generating output filename
    TYPE(t_fname_metadata),                INTENT(IN)  :: fname_metadata
    INTEGER,                               INTENT(IN)  :: local_event_no       !< local index of this event on local PE
    INTEGER,                               INTENT(IN)  :: icomm                !< MPI communicator


    INTERFACE
      !> As an argument of this function, the user must provide a
      !! conversion "time stamp -> simulation step"
      SUBROUTINE fct_time2simstep(nstrings, date_string, sim_step_info, &
        &                         result_steps, result_exactdate)
        USE mo_output_event_types, ONLY: t_sim_step_info
        INTEGER,                   INTENT(IN)    :: nstrings             !< no. of string to convert
        CHARACTER(len=*),          INTENT(IN)    :: date_string(:)       !< array of ISO 8601 time stamp strings
        TYPE(t_sim_step_info),     INTENT(IN)    :: sim_step_info        !< definitions: time step size, etc.
        INTEGER,                   INTENT(INOUT) :: result_steps(:)      !< resulting step indices
        CHARACTER(LEN=*),          INTENT(INOUT) :: result_exactdate(:)  !< resulting (exact) time step strings
      END SUBROUTINE fct_time2simstep

      !> As an argument of this function, the user must provide a
      !! function for generating output file names
      SUBROUTINE fct_generate_filenames(nstrings, date_string, sim_steps, &
        &                               sim_step_info, fname_metadata, &
        &                               skipped_dates, result_fnames)
        USE mo_output_event_types,     ONLY: t_sim_step_info, t_event_step_data
        USE mo_name_list_output_types, ONLY: t_fname_metadata

        INTEGER,                   INTENT(IN)    :: nstrings           !< no. of string to convert
        CHARACTER(len=*),          INTENT(IN)    :: date_string(:)     !< array of ISO 8601 time stamp strings
        INTEGER,                   INTENT(IN)    :: sim_steps(:)       !< array of corresponding simulation steps
        TYPE(t_sim_step_info),     INTENT(IN)    :: sim_step_info      !< definitions: time step size, etc.
        TYPE(t_fname_metadata),    INTENT(IN)    :: fname_metadata     !< additional meta-data for generating output filename
        INTEGER,                   INTENT(IN)    :: skipped_dates
        TYPE(t_event_step_data),   INTENT(OUT)   :: result_fnames(SIZE(date_string))
      END SUBROUTINE fct_generate_filenames
    END INTERFACE

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::new_parallel_output_event"
    INTEGER :: ierrstat, this_pe, i_tag, nranks, i, nintvls
    TYPE(t_event_data_local) :: evd

    ! determine this PE's MPI rank wrt. the given MPI communicator:
    this_pe = 0
    nranks  = 1
#ifndef NOMPI
    IF (icomm /= MPI_COMM_NULL) THEN
      this_pe = p_comm_rank(icomm)
      nranks = p_comm_size(icomm)
      IF (ldebug) THEN
        WRITE (0,*) "PE ",get_my_global_mpi_id(), ": local rank is ", this_pe, "; icomm has size ", nranks
      END IF
    END IF
#endif

    ! compute i_tag ID st. it stays unique even for multiple events
    ! running on the same I/O PE:
    i_tag = SENDRECV_TAG_OUTEVENT + this_pe + (local_event_no-1)*nranks

    ! allocate parallel event data structure
    ALLOCATE(p_event, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    NULLIFY(p_event%next)

    ! count the number of different time intervals for this event (usually 1)
    DO nintvls = 0, MAX_TIME_INTERVALS-1
      IF (LEN_TRIM(begin_str(nintvls+1)) == 0) EXIT
    END DO

    ! create the local (non-parallel) output event data structure:
    evd%name = name
    evd%begin_str = begin_str
    evd%end_str = end_str
    evd%intvl_str = intvl_str
    evd%l_output_last = l_output_last
    evd%sim_step_info = sim_step_info
    evd%fname_metadata = fname_metadata
    evd%i_tag = i_tag
    p_event%output_event => new_output_event(evd, this_pe,         &
      &                                      fct_time2simstep,     &
      &                                      fct_generate_filenames)

    IF (ldebug) THEN
      WRITE (0,*) "PE ", get_my_global_mpi_id(), ": created event with ", &
        &      p_event%output_event%n_event_steps, " steps."
    END IF

    ! set the other MPI-related data fields:
    p_event%icomm         = icomm
    p_event%iroot         = ROOT_OUTEVENT
    p_event%irecv_nreq    = 0
    p_event%isend_req     = 0
#ifndef NOMPI
    p_event%isend_req     = MPI_REQUEST_NULL
    IF (this_pe /= ROOT_OUTEVENT) THEN
      ! now send all the meta-data to the root PE
      IF (icomm /= MPI_COMM_NULL) THEN
        IF (ldebug) THEN
          WRITE (0,*) "new_parallel_output_event: PE ", get_my_global_mpi_id(), ": send event data: ", &
            &      TRIM(begin_str(1)), TRIM(end_str(1)), TRIM(intvl_str(1))
          DO i=2,nintvls
            WRITE (0,*) " + PE ", get_my_global_mpi_id(), ": send event data: ", &
              &      TRIM(begin_str(i)), TRIM(end_str(i)), TRIM(intvl_str(i))
          END DO
        END IF
        CALL send_event_data(evd, icomm, ROOT_OUTEVENT)
      END IF
    ELSE
#endif
      ! I/O PE #0: we keep a local list of event meta-data
      ievent_list_local = ievent_list_local + 1
      event_list_local(ievent_list_local)%name               = name
      event_list_local(ievent_list_local)%begin_str(:)       = begin_str(:)
      event_list_local(ievent_list_local)%end_str(:)         = end_str(:)
      event_list_local(ievent_list_local)%intvl_str(:)       = intvl_str(:)
      event_list_local(ievent_list_local)%l_output_last      = l_output_last
      event_list_local(ievent_list_local)%sim_step_info      = sim_step_info
      event_list_local(ievent_list_local)%fname_metadata     = fname_metadata
      event_list_local(ievent_list_local)%i_tag              = i_tag
#ifndef NOMPI
    END IF
#endif
  END FUNCTION new_parallel_output_event


  !---------------------------------------------------------------
  ! ROUTINES FOR JOINING OUTPUT EVENTS
  !---------------------------------------------------------------

  !> Receives event meta-data from all PEs within the given MPI
  !> communicator; creates the union of these events.
  !
  !  Optional: Broadcast events via an inter-communicator, eg. to
  !            worker PEs
  !
  !  @author F. Prill, DWD
  !
  FUNCTION union_of_all_events(fct_time2simstep, fct_generate_filenames, icomm, &
    &                          opt_broadcast_comm, opt_broadcast_root)
    TYPE(t_par_output_event), POINTER :: union_of_all_events
    INTEGER, OPTIONAL,      INTENT(IN)  :: icomm                       !< MPI communicator for intra-I/O communication
    INTEGER, OPTIONAL,      INTENT(IN)  :: opt_broadcast_comm          !< MPI communicator for broadcast IO->workers
    INTEGER, OPTIONAL,      INTENT(IN)  :: opt_broadcast_root          !< MPI rank (broadcast source)

    INTERFACE
      !> As an argument of this function, the user must provide a
      !! conversion "time stamp -> simulation step"
      SUBROUTINE fct_time2simstep(nstrings, date_string, sim_step_info, &
        &                         result_steps, result_exactdate)
        USE mo_output_event_types, ONLY: t_sim_step_info
        INTEGER,                  INTENT(IN)    :: nstrings             !< no. of string to convert
        CHARACTER(len=*),         INTENT(IN)    :: date_string(:)       !< array of ISO 8601 time stamp strings
        TYPE(t_sim_step_info),    INTENT(IN)    :: sim_step_info        !< definitions: time step size, etc.
        INTEGER,                  INTENT(INOUT) :: result_steps(:)      !< resulting step indices
        CHARACTER(LEN=*),         INTENT(INOUT) :: result_exactdate(:)  !< resulting (exact) time step strings
      END SUBROUTINE fct_time2simstep

      !> As an argument of this function, the user must provide a
      !! function for generating output file names
      SUBROUTINE fct_generate_filenames(nstrings, date_string, sim_steps, &
        &                               sim_step_info, fname_metadata, &
        &                               skipped_dates, result_fnames)
        USE mo_output_event_types,     ONLY: t_sim_step_info, t_event_step_data
        USE mo_name_list_output_types, ONLY: t_fname_metadata

        INTEGER,                   INTENT(IN)    :: nstrings           !< no. of string to convert
        CHARACTER(len=*),          INTENT(IN)    :: date_string(:)     !< array of ISO 8601 time stamp strings
        INTEGER,                   INTENT(IN)    :: sim_steps(:)       !< array of corresponding simulation steps
        TYPE(t_sim_step_info),     INTENT(IN)    :: sim_step_info      !< definitions: time step size, etc.
        TYPE(t_fname_metadata),    INTENT(IN)    :: fname_metadata     !< additional meta-data for generating output filename
        INTEGER,                   INTENT(IN)    :: skipped_dates
        TYPE(t_event_step_data),  INTENT(OUT)    :: result_fnames(SIZE(date_string))
      END SUBROUTINE fct_generate_filenames
    END INTERFACE

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::union_of_all_events"
    TYPE(t_par_output_event), POINTER     :: par_event, last_node
    TYPE(t_output_event),     POINTER     :: ev2
    TYPE(t_event_data_local)              :: evd
    LOGICAL                               :: lrecv
    INTEGER                               :: i_pe, nranks, ierrstat,            &
      &                                      this_pe, nbcast_ranks, i
    LOGICAL                               :: lbroadcast

    IF (ldebug)  WRITE (0,*) routine, " enter."
#ifndef NOMPI
    lbroadcast = PRESENT(opt_broadcast_root) .AND. &
         &       PRESENT(opt_broadcast_comm) .AND. &
         &       (.NOT. my_process_is_mpi_test())
#else
    lbroadcast = .FALSE.
#endif

    NULLIFY(union_of_all_events)
    ! get the number of ranks in this MPI communicator
    nranks  =  1
    this_pe = -1
#ifndef NOMPI
    IF (icomm /= MPI_COMM_NULL) THEN
      this_pe = p_comm_rank(icomm)
      nranks = p_comm_size(icomm)
      IF (ldebug) THEN
         WRITE (0,*) "PE ",get_my_global_mpi_id(), ": local rank is ", this_pe, "; icomm has size ", nranks
         IF (lbroadcast) THEN
            nbcast_ranks = p_comm_size(opt_broadcast_comm)
            WRITE (0,*) "PE ",get_my_global_mpi_id(), ": local rank is ", this_pe, "; bcast comm has size ", nbcast_ranks, &
                 &      ", root is ", opt_broadcast_root
            lbroadcast = (nbcast_ranks > 1)
         END IF
      END IF
    END IF
#endif
    IF (lbroadcast)  CALL p_bcast(nranks, opt_broadcast_root, opt_broadcast_comm)

    ! loop over all PEs of the I/O communicator:
    DO i_pe = 0,(nranks-1)
      RECEIVE_LOOP : DO
        lrecv = .FALSE.
#ifndef NOMPI
        ! receive event meta-data from participating I/O PEs:
        IF (this_pe == i_pe) THEN
          lrecv = .FALSE.
        ELSE
          IF (this_pe == ROOT_OUTEVENT) THEN
            lrecv = receive_event_data(evd, icomm, i_pe, SENDRECV_TAG_SETUP)
          END IF
        END IF
#endif
        ! I/O PE #0: we keep a local list of event meta-data
        IF (.NOT. lrecv .AND. (ievent_list_local > 0)) THEN
          lrecv = .TRUE.
          evd = event_list_local(ievent_list_local)
          ievent_list_local = ievent_list_local - 1
          IF (ldebug) THEN
            DO i=1,MAX_TIME_INTERVALS
              WRITE (0,*) "Taking event ", TRIM(evd%begin_str(i)), " / ", TRIM(evd%end_str(i)), " / ", &
                &         TRIM(evd%intvl_str(i)), " from local list."
            END DO
            WRITE (0,*) ievent_list_local, " entries are left."
          END IF
        END IF
        ! forward the event meta-data to the worker PEs
        IF (lbroadcast) THEN
          lrecv = broadcast_event_data(evd,                                    &
            &                          opt_broadcast_comm, opt_broadcast_root, &
            &                          lrecv)
        END IF
        
        IF (.NOT. lrecv) EXIT RECEIVE_LOOP

        ! create the event steps from the received meta-data:
        ev2 => new_output_event(evd, i_pe, &
          &                     fct_time2simstep, fct_generate_filenames)

        ! find the parallel output event in the linked list that
        ! matches the event name
        NULLIFY(last_node)
        par_event => union_of_all_events
        DO WHILE (ASSOCIATED(par_event))
          IF (par_event%output_event%event_data%name == ev2%event_data%name) EXIT
          last_node => par_event
          par_event => par_event%next
        END DO
        ! if there is no such parallel output event, then create one
        ! at the end of the linked list:
        IF (.NOT. ASSOCIATED(par_event)) THEN
          IF (.NOT. ASSOCIATED(last_node)) THEN
            ALLOCATE(union_of_all_events, STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
            par_event => union_of_all_events
          ELSE
            ALLOCATE(last_node%next, STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
            par_event => last_node%next
          END IF
          ! set the MPI-related data fields:
          par_event%icomm         = icomm
          par_event%iroot         = ROOT_OUTEVENT
          par_event%irecv_nreq    = 0
          NULLIFY(par_event%next)
          NULLIFY(par_event%output_event)
        END IF
        ! increase MPI message tag ID st. it stays unique even for
        ! multiple events running on the same I/O PE:
        IF (ASSOCIATED(par_event%output_event)) THEN
          ! create union of the two events:
          CALL merge_events(par_event%output_event, ev2)
          CALL deallocate_output_event(ev2)
        ELSE
          par_event%output_event => ev2
        END IF
        IF (ldebug) THEN
          WRITE (0,*) "PE ", get_my_global_mpi_id(), ": n_event_steps = ", &
            & union_of_all_events%output_event%n_event_steps
        END IF
      END DO RECEIVE_LOOP
    END DO
    IF (ldebug)  WRITE (0,*) routine, " done."
  END FUNCTION union_of_all_events


  !> Merge compatible event into another.
  !
  !  Different output events/PEs may be joined together to form a
  !  single new event. Then not every PE necessarily participates in
  !  every event step. The purpose of joining/grouping events is that
  !  certain actions can be triggered for the group, e.g. writing
  !  "ready files".
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE merge_events(event1, event2)
    TYPE(t_output_event), INTENT(INOUT) :: event1
    TYPE(t_output_event), INTENT(IN) :: event2
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::merge_event"
    INTEGER                             :: i1, i2, ierrstat, i_sim_step1, i_sim_step2, &
      &                                    max_sim_step, j1, j2, nsteps
    LOGICAL :: copy_i1
    CHARACTER(LEN=MAX_FILENAME_STR_LEN) :: filename_string1
    TYPE(t_event_step), ALLOCATABLE :: new_steps(:)
    ! allocate event data structure
    event1%i_event_step             = 1
    IF (event1%event_data%name /= event2%event_data%name) THEN
      i1 = LEN_TRIM(event1%event_data%name)
      event1%event_data%name(i1:) = ","//event2%event_data%name
    END IF
    IF (event1%event_data%sim_start /= event2%event_data%sim_start) THEN
      CALL finish(routine, "Simulation start dates do not match!")
    END IF

    ! consistency check: test, if any of the filenames in event1
    ! occurs in event2 (avoid duplicate names)
    DO i1=1,event1%n_event_steps
      DO j1=1,event1%event_step(i1)%n_pes
        IF (.NOT. event1%event_step(i1)%event_step_data(j1)%l_open_file) CYCLE
        filename_string1 = event1%event_step(i1)%event_step_data(j1)%filename_string
        DO i2=1,event2%n_event_steps
          DO j2=1,event2%event_step(i2)%n_pes
            IF (.NOT. event2%event_step(i2)%event_step_data(j2)%l_open_file) CYCLE
            IF (event2%event_step(i2)%event_step_data(j2)%filename_string == filename_string1) THEN
              ! found a duplicate filename:
              CALL finish(routine, "Error! Ambiguous output file name: '"//TRIM(filename_string1)//"'")
            END IF
          END DO
        END DO
      END DO
    END DO

    ! choose a maximum simulation step number as abort criterion:
    max_sim_step = 0
    IF (event1%n_event_steps > 0)  max_sim_step = MAX(max_sim_step, event1%event_step(event1%n_event_steps)%i_sim_step)
    IF (event2%n_event_steps > 0)  max_sim_step = MAX(max_sim_step, event2%event_step(event2%n_event_steps)%i_sim_step)
    max_sim_step = max_sim_step + 1

    nsteps = 0
    i1 = 1
    i2 = 1
    DO WHILE (i1 <= event1%n_event_steps .OR. i2 <= event2%n_event_steps)
      IF (i1 > event1%n_event_steps) THEN
        i_sim_step1 = max_sim_step
      ELSE
        i_sim_step1 = event1%event_step(i1)%i_sim_step
      END IF
      IF (i2 > event2%n_event_steps) THEN
        i_sim_step2 = max_sim_step
      ELSE
        i_sim_step2 = event2%event_step(i2)%i_sim_step
      END IF
      nsteps = nsteps+1
      ! copy event step i1 from event1:
      i1 = i1 + MERGE(1, 0, i_sim_step2 >= i_sim_step1)
      ! copy event step i2 from event2:
      i2 = i2 + MERGE(1, 0, i_sim_step2 <= i_sim_step1)
    END DO

    ! now create the new event:
    ALLOCATE(new_steps(nsteps), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    nsteps = 0
    i1 = 1
    i2 = 1
    DO WHILE (i1 <= event1%n_event_steps .OR. i2 <= event2%n_event_steps)
      IF (i1 > event1%n_event_steps) THEN
        i_sim_step1 = max_sim_step
      ELSE
        i_sim_step1 = event1%event_step(i1)%i_sim_step
      END IF
      IF (i2 > event2%n_event_steps) THEN
        i_sim_step2 = max_sim_step
      ELSE
        i_sim_step2 = event2%event_step(i2)%i_sim_step
      END IF
      nsteps = nsteps+1
      copy_i1 = i_sim_step1 <= i_sim_step2
      IF (copy_i1) THEN
        ! copy event step i1 from event1:
        new_steps(nsteps)%i_sim_step = event1%event_step(i1)%i_sim_step
        new_steps(nsteps)%exact_date_string &
          &    = event1%event_step(i1)%exact_date_string
        new_steps(nsteps)%n_pes = event1%event_step(i1)%n_pes
        CALL MOVE_ALLOC(from=event1%event_step(i1)%event_step_data, &
          &             to=new_steps(nsteps)%event_step_data)
        i1 = i1 + 1
      END IF
      IF (i_sim_step2 <= i_sim_step1) THEN
        ! copy event step i2 from event2:
        CALL append_event_step(new_steps(nsteps), event2%event_step(i2), l_create=(.NOT. copy_i1))
        i2 = i2 + 1
      END IF
    END DO
    CALL MOVE_ALLOC(from=new_steps, to=event1%event_step)
    event1%n_event_steps = nsteps
  END SUBROUTINE merge_events


  !> Appends the data of an event step to the data of another event step.
  !  @author F. Prill, DWD
  !
  SUBROUTINE append_event_step(dst_event_step, src_event_step, l_create)
    TYPE(t_event_step), INTENT(INOUT)  :: dst_event_step  !< source event step
    TYPE(t_event_step), INTENT(IN)     :: src_event_step  !< destination event step
    LOGICAL,            INTENT(IN)     :: l_create        !< Flag. If .TRUE., the the destination event step is created
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::append_event_step"
    TYPE(t_event_step_data), ALLOCATABLE  :: tmp_event_step_data(:)
    INTEGER :: iold, ierrstat, i1, i2, inew, iadd

    IF (l_create) THEN
      dst_event_step%i_sim_step = src_event_step%i_sim_step
      dst_event_step%exact_date_string = src_event_step%exact_date_string
    ELSE
      IF (src_event_step%i_sim_step /= dst_event_step%i_sim_step) &
        &   CALL finish(routine, "sim_steps do not match!")
      IF (src_event_step%exact_date_string /= dst_event_step%exact_date_string) &
        &   CALL finish(routine, "exact_date_strings do not match!")
    END IF
    ! event_step%event_step_data
    IF (l_create) THEN
      dst_event_step%n_pes = src_event_step%n_pes
      ALLOCATE(dst_event_step%event_step_data(dst_event_step%n_pes), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      dst_event_step%event_step_data(1:dst_event_step%n_pes) = &
        &   src_event_step%event_step_data(1:src_event_step%n_pes)
    ELSE
      iold = dst_event_step%n_pes
      iadd = src_event_step%n_pes
      inew = iold + iadd
      ALLOCATE(tmp_event_step_data(inew), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      tmp_event_step_data(1:iold) = dst_event_step%event_step_data(1:iold)
      tmp_event_step_data(iold+1:inew) = &
           & src_event_step%event_step_data(1:iadd)
      CALL MOVE_ALLOC(from=tmp_event_step_data, &
        &             to=dst_event_step%event_step_data)
      dst_event_step%n_pes = inew

      ! consistency check: test, if any of the filenames in group1
      ! occurs in group2 (avoid duplicate names)
      DO i1=1,iold
        DO i2=(iold+1),(iold+src_event_step%n_pes)
          IF (dst_event_step%event_step_data(i1)%filename_string &
            == dst_event_step%event_step_data(i2)%filename_string) THEN
            ! found a duplicate filename:
            CALL finish(routine, "Ambiguous output file name: '"//TRIM(dst_event_step%event_step_data(i1)%filename_string)//"'")
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE append_event_step


  !---------------------------------------------------------------
  ! ROUTINES FOR TESTING / INQURING OUTPUT EVENTS
  !---------------------------------------------------------------

  !> @return .TRUE. if the output event has exceeded its final step.
  !  @author F. Prill, DWD
  !
  FUNCTION is_output_event_finished(event)
    LOGICAL :: is_output_event_finished
    TYPE(t_output_event), INTENT(IN) :: event
    is_output_event_finished = (event%i_event_step > event%n_event_steps)
  END FUNCTION is_output_event_finished


  !> @return .TRUE. if the output event has exceeded its final step.
  !  @author F. Prill, DWD
  !
  FUNCTION is_par_output_event_finished(event)
    LOGICAL :: is_par_output_event_finished
    TYPE(t_par_output_event), POINTER             :: event

    IF (.NOT. ASSOCIATED(event)) THEN
      is_par_output_event_finished = .TRUE.
    ELSE
      is_par_output_event_finished = is_output_event_finished(event%output_event)
    END IF
  END FUNCTION is_par_output_event_finished


  !> @return .TRUE. if the given step index matches the step index of
  !          the current event step, ie. if the current event step is
  !          active.
  !
  !  @author F. Prill, DWD
  !
  FUNCTION is_output_step(event, jstep)
    LOGICAL :: is_output_step
    TYPE(t_output_event), INTENT(IN) :: event  !< output event
    INTEGER,              INTENT(IN) :: jstep  !< given step index
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::is_output_step"
    INTEGER :: istep

    IF (is_output_event_finished(event)) THEN
      is_output_step = .FALSE.
    ELSE
      istep = event%i_event_step
      is_output_step = (event%event_step(istep)%i_sim_step == jstep)
    END IF
  END FUNCTION is_output_step


  !> @return .TRUE. if the given step index matches the step index of
  !          the current event step, ie. if the current event step is
  !          active.
  !
  !  @note This function recursively checks the complete
  !        singly-linked list rooted at @p event.
  !
  !  @author F. Prill, DWD
  !
  RECURSIVE FUNCTION is_par_output_step(event, jstep)
    LOGICAL :: is_par_output_step
    TYPE(t_par_output_event), POINTER                :: event  !< output event
    INTEGER,                           INTENT(IN)    :: jstep  !< given step index
    ! local variables
    LOGICAL :: ret, ret_local

    ret = .FALSE.
    IF (ASSOCIATED(event)) THEN
      ! first, check other output events in linked list:
      IF (ASSOCIATED(event%next)) THEN
        ret = ret .OR. is_output_step(event%next, jstep)
      END IF
      ret_local = is_output_step(event%output_event, jstep)
      ret = ret .OR. ret_local
    END IF
    is_par_output_step = ret
  END FUNCTION is_par_output_step


  !> @return .TRUE. if all participants of the parallel output event
  !>         have acknowledged the completion of the current step.
  !
  !  @note We do not use the wrapper routines from "mo_mpi" here,
  !        since we need direct control over the non-blocking P2P
  !        requests.
  !  @author F. Prill, DWD
  !
  FUNCTION is_output_step_complete(event) RESULT(ret)
    LOGICAL :: ret
    TYPE(t_par_output_event), INTENT(inout) :: event
    ! local variables
#ifndef NOMPI
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::is_output_step_complete"
    INTEGER              :: ierrstat
    ! INTEGER, ALLOCATABLE :: irecv_status(:,:)

    ret = .TRUE.
    IF (event%irecv_nreq /= 0) THEN
      ! ALLOCATE(irecv_status(MPI_STATUS_SIZE,event%irecv_nreq), STAT=ierrstat)
      ! IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      CALL MPI_TESTALL(event%irecv_nreq, event%irecv_req, ret, &
        &              mpi_statuses_ignore, ierrstat)
      IF (ierrstat /= mpi_success) CALL finish (routine, 'Error in MPI_TESTALL.')
      ! DEALLOCATE(irecv_status, STAT=ierrstat)
      ! IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    END IF
#else
    ret = .TRUE.
#endif
  END FUNCTION is_output_step_complete


  !> @return current date-time stamp string.
  !  @author F. Prill, DWD
  !
  FUNCTION get_current_date(event)
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: get_current_date
    TYPE(t_output_event), INTENT(in) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::get_current_date"
    INTEGER :: istep

    istep = MIN(event%i_event_step, event%n_event_steps)
    ! For events that are the union of multiple other events we return
    ! the *model date-time string*:
    IF (event%event_step(istep)%n_pes > 1) THEN
      get_current_date = event%event_step(istep)%exact_date_string
    ELSE
      get_current_date = event%event_step(istep)%event_step_data(1)%datetime_string
    END IF
  END FUNCTION get_current_date


  !> @return current date-time stamp string.
  !  @author F. Prill, DWD
  !
  FUNCTION get_current_date_par(event)
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: get_current_date_par
    TYPE(t_par_output_event), POINTER :: event
    IF (.NOT. ASSOCIATED(event)) THEN
      get_current_date_par = ""
    ELSE
      get_current_date_par =  get_current_date(event%output_event)
    END IF
  END FUNCTION get_current_date_par


  !> @return current output step.
  !  @author F. Prill, DWD
  !
  FUNCTION get_current_step(event)
    INTEGER :: get_current_step
    TYPE(t_output_event), POINTER :: event
    get_current_step = MIN(event%i_event_step, event%n_event_steps)
  END FUNCTION get_current_step


  !> @return current output step.
  !  @author F. Prill, DWD
  !
  FUNCTION get_current_step_par(event)
    INTEGER :: get_current_step_par
    TYPE(t_par_output_event), POINTER :: event
    IF (.NOT. ASSOCIATED(event)) THEN
      get_current_step_par = -1
    ELSE
      get_current_step_par =  get_current_step(event%output_event)
    END IF
  END FUNCTION get_current_step_par


  !> @return current filename string.
  !  @author F. Prill, DWD
  !
  FUNCTION get_current_filename(event)
    CHARACTER(LEN=MAX_FILENAME_STR_LEN) :: get_current_filename
    TYPE(t_par_output_event), INTENT(IN) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::get_current_filename"
    INTEGER :: istep

    istep = event%output_event%i_event_step
    ! We cannot check events that are the union of multiple other
    ! events:
    IF (event%output_event%event_step(istep)%n_pes > 1) THEN
      WRITE (0,*) "Event step ", istep,  "(", TRIM(event%output_event%event_step(istep)%exact_date_string), ")", &
        &         ", shared by ", event%output_event%event_step(istep)%n_pes, " PEs."
      CALL finish(routine, "Error! Multi-part event step!")
    END IF
    get_current_filename = TRIM(event%output_event%event_step(istep)%event_step_data(1)%filename_string)
  END FUNCTION get_current_filename


  !> @return current file number.
  !  @author F. Prill, DWD
  !
  FUNCTION get_current_jfile(event)
    INTEGER :: get_current_jfile
    TYPE(t_par_output_event), INTENT(IN) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::get_current_jfile"
    INTEGER :: istep

    istep = MIN(event%output_event%i_event_step, event%output_event%n_event_steps)
    IF (istep == 0) THEN
      get_current_jfile = 0
      RETURN
    END IF
    ! We cannot check events that are the union of multiple other
    ! events:
    IF (event%output_event%event_step(istep)%n_pes > 1) THEN
      WRITE (0,*) "Event step ", istep,  "(", TRIM(event%output_event%event_step(istep)%exact_date_string), ")", &
        &         ", shared by ", event%output_event%event_step(istep)%n_pes, " PEs."
      CALL finish(routine, "Error! Multi-part event step!")
    END IF
    get_current_jfile = event%output_event%event_step(istep)%event_step_data(1)%jfile
  END FUNCTION get_current_jfile


  !> @return .TRUE. if this PE should write a ready file for the given
  !          event.
  !  @author F. Prill, DWD
  !
  FUNCTION check_write_readyfile(event)
    LOGICAL :: check_write_readyfile
    TYPE(t_output_event), POINTER :: event

    IF (ASSOCIATED(event)) THEN
      check_write_readyfile = (TRIM(event%event_data%name) /= DEFAULT_EVENT_NAME) .AND.  &
        &                     (event%n_event_steps > 0)
    ELSE
      check_write_readyfile = .FALSE.
    END IF
  END FUNCTION check_write_readyfile


  !> @return .TRUE. if current event step has the "open file" flag
  !          enabled.
  !
  !  @author F. Prill, DWD
  !
  FUNCTION check_open_file(event)
    LOGICAL :: check_open_file
    TYPE(t_par_output_event), INTENT(IN) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::check_open_file"
    INTEGER :: istep

    istep = event%output_event%i_event_step
    ! We cannot check events that are the union of multiple other
    ! events:
    IF (event%output_event%event_step(istep)%n_pes > 1) THEN
      WRITE (0,*) "Event step ", istep,  "(", TRIM(event%output_event%event_step(istep)%exact_date_string), ")", &
        &         ", shared by ", event%output_event%event_step(istep)%n_pes, " PEs."
      CALL finish(routine, "Error! Multi-part event step!")
    END IF
    check_open_file = event%output_event%event_step(istep)%event_step_data(1)%l_open_file
  END FUNCTION check_open_file


  !> @return .TRUE. if current event step has the "close file" flag
  !          enabled.
  !
  !  @author F. Prill, DWD
  !
  FUNCTION check_close_file(event)
    LOGICAL :: check_close_file
    TYPE(t_par_output_event), INTENT(IN) :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::check_close_file"
    INTEGER :: istep

    istep = event%output_event%i_event_step
    ! We cannot check events that are the union of multiple other
    ! events:
    IF (event%output_event%event_step(istep)%n_pes > 1) THEN
      WRITE (0,*) "Event step ", istep,  "(", TRIM(event%output_event%event_step(istep)%exact_date_string), ")", &
        &         ", shared by ", event%output_event%event_step(istep)%n_pes, " PEs."
      CALL finish(routine, "Error! Multi-part event step!")
    END IF
    check_close_file = event%output_event%event_step(istep)%event_step_data(1)%l_close_file
  END FUNCTION check_close_file


  !> @return .TRUE. if this PE is the event's root PE.  If no event
  !          data structure has been provided as an argument to this
  !          function, we test for ROOT_OUTEVENT.
  !
  !  @author F. Prill, DWD
  !
  FUNCTION is_event_root_pe(opt_event, opt_icomm)
    LOGICAL :: is_event_root_pe
    TYPE(t_par_output_event), POINTER, OPTIONAL :: opt_event
    INTEGER, INTENT(IN),               OPTIONAL :: opt_icomm
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::is_event_root_pe"
#ifndef NOMPI
    INTEGER :: this_pe, p_comm
#endif

    ! at least one of the two arguments must be provided:
    IF (.NOT. PRESENT(opt_event) .AND. .NOT. PRESENT(opt_icomm)) THEN
      CALL finish(routine, "Internal error!")
    END IF

#ifndef NOMPI
    IF (PRESENT(opt_event)) THEN
      p_comm = opt_event%icomm
    ELSE
      p_comm = opt_icomm
    END IF
    ! determine this PE's MPI rank wrt. the given MPI communicator and
    ! return if this PE is not root PE for the given event:
    this_pe = p_comm_rank(p_comm)
    IF (PRESENT(opt_event)) THEN
      is_event_root_pe = (this_pe == opt_event%iroot)
    ELSE
      is_event_root_pe = (this_pe == ROOT_OUTEVENT)
    END IF
#else
    is_event_root_pe = .TRUE.
#endif
  END FUNCTION is_event_root_pe


  !---------------------------------------------------------------
  ! ROUTINES PERFORMING DATA TRANSFER TO ROOT PE DURING SETUP
  !---------------------------------------------------------------

  !> MPI-send event data to root PE.
  !
  !  All data necessary to create the event dates must be transmitted.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE send_event_data(evd, icomm, dst_rank)
    TYPE(t_event_data_local),            INTENT(IN)  :: evd
    INTEGER,                             INTENT(IN)  :: icomm                !< MPI communicator
    INTEGER,                             INTENT(IN)  :: dst_rank             !< MPI destination rank
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::send_event_data"
    INTEGER :: nitems, ierrstat, position
    CHARACTER, ALLOCATABLE :: buffer(:)   !< MPI buffer for packed

    ! allocate message buffer
    ALLOCATE(buffer(MAX_BUF_SIZE), stat=ierrstat)
    IF (ierrstat /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed')
    position = 0

    ! prepare an MPI message:
    !
    nitems = 1
    CALL p_pack_int(nitems, buffer, position, icomm)
    CALL pack_metadata(buffer, position, evd, icomm)

    ! send packed message:
    CALL p_send_packed(buffer, dst_rank, SENDRECV_TAG_SETUP, position, icomm)

    ! clean up
    DEALLOCATE(buffer, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE send_event_data


  !> Send an empty record of event data to root PE st. the receiver
  !> knows that no more events will be transmitted.
  SUBROUTINE complete_event_setup(icomm)
    INTEGER,                INTENT(IN)  :: icomm                 !< MPI communicator
    ! local variables
#ifndef NOMPI
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::complete_event_setup"
    INTEGER :: nitems, ierrstat, position, this_pe
    CHARACTER, ALLOCATABLE :: buffer(:)   !< MPI buffer for packed

    IF (icomm /= MPI_COMM_NULL) THEN
      this_pe = p_comm_rank(icomm)

      IF (this_pe /= ROOT_OUTEVENT) THEN
        ! allocate message buffer
        ALLOCATE(buffer(MAX_BUF_SIZE), stat=ierrstat)
        IF (ierrstat /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed')
        position = 0
        ! prepare an empty MPI message:
        nitems = 0
        CALL p_pack_int(nitems, buffer, position, icomm)
        
        ! send packed message:
        CALL p_send_packed(buffer, ROOT_OUTEVENT, SENDRECV_TAG_SETUP, position, icomm)
        
        ! clean up
        DEALLOCATE(buffer, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
      END IF
    END IF
#endif
  END SUBROUTINE complete_event_setup


  !> MPI-receive event data.
  !
  !  @return .FALSE. if no more event data is to be received from the
  !          given PE.
  !
  !  @author F. Prill, DWD
  !
  FUNCTION receive_event_data(evd, icomm, isrc, isendrecv_tag)
    LOGICAL :: receive_event_data
    TYPE(t_event_data_local),              INTENT(INOUT)  :: evd
    !> MPI communicator
    INTEGER,                               INTENT(IN)     :: icomm
    !> MPI rank of the sending PE
    INTEGER,                               INTENT(IN)     :: isrc
    !> MPI tag for this messages isend/irecv communication
    INTEGER,                               INTENT(IN)     :: isendrecv_tag
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::receive_event_data"
    INTEGER :: nitems, ierrstat, position, i
    CHARACTER, ALLOCATABLE :: buffer(:)   !< MPI buffer for packed

    ! allocate message buffer
    ALLOCATE(buffer(MAX_BUF_SIZE), stat=ierrstat)
    IF (ierrstat /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed')
    position = 0

    ! receive packed message:
    CALL p_irecv_packed(buffer, isrc, isendrecv_tag, MAX_BUF_SIZE, icomm)
    ! wait for message to arrive:
    CALL p_wait()

    ! unpack MPI message:
    !
    nitems = 1
    CALL p_unpack_int(buffer, position, nitems, icomm)
    IF (nitems == 0) THEN
      receive_event_data = .FALSE.
    ELSE
      receive_event_data = .TRUE.
      CALL unpack_metadata(buffer, position, evd, icomm)
    END IF

    IF (ldebug) THEN
      DO i=1,MAX_TIME_INTERVALS
        WRITE (0,*) "received event data: ", TRIM(evd%begin_str(i)), TRIM(evd%end_str(i)), TRIM(evd%intvl_str(i))
      END DO
    END IF

    ! clean up
    DEALLOCATE(buffer, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END FUNCTION receive_event_data


  !> MPI-broadcast event data.
  !
  !  @return .FALSE. if no more event data is to be received from the
  !          given PE.
  !
  !  @author F. Prill, DWD
  !
  FUNCTION broadcast_event_data(evd, icomm, iroot, l_no_end_message)
    LOGICAL :: broadcast_event_data
    TYPE(t_event_data_local),              INTENT(INOUT)  :: evd
    !> MPI communicator
    INTEGER,                               INTENT(IN)     :: icomm
    !> MPI broadcast root rank
    INTEGER,                               INTENT(IN)     :: iroot
    !> Flag. .FALSE. if "end message" shall be broadcasted
    LOGICAL,                               INTENT(IN)     :: l_no_end_message
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::broadcast_event_data"
    INTEGER :: nitems, ierrstat, position, this_pe, i
    CHARACTER, ALLOCATABLE :: buffer(:)   !< MPI buffer for packed

    ! allocate message buffer
    ALLOCATE(buffer(MAX_BUF_SIZE), stat=ierrstat)
    IF (ierrstat /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed')

    ! determine this PE's MPI rank wrt. the given MPI communicator:
    this_pe = -1
#ifndef NOMPI
    IF (icomm /= MPI_COMM_NULL) this_pe = p_comm_rank(icomm)
#endif

    ! prepare an MPI message:
    IF (iroot == this_pe) THEN
      position = 0
      IF (l_no_end_message) THEN
        ! normal message:
        nitems = 1
        CALL p_pack_int(nitems, buffer, position, icomm)
        CALL pack_metadata(buffer, position, evd, icomm)
        IF (ldebug) THEN
          DO i=1,MAX_TIME_INTERVALS
            WRITE (0,*) "PE ", get_my_global_mpi_id(), ": send event data: ", &
              &         TRIM(evd%begin_str(i)), TRIM(evd%end_str(i)), TRIM(evd%intvl_str(i))
          END DO
        END IF
      ELSE
        ! empty, "end message":
        nitems = 0
        CALL p_pack_int(nitems, buffer, position, icomm)
        IF (ldebug) THEN
          WRITE (0,*) "PE ", get_my_global_mpi_id(), ": send end message after "
        END IF
      END IF
    END IF

#ifndef NOMPI
    ! broadcast packed message:
    CALL MPI_BCAST(buffer, MAX_BUF_SIZE, MPI_PACKED, iroot, icomm, ierrstat)
    IF (ierrstat /= 0) CALL finish (routine, 'Error in MPI_BCAST.')
#endif

    ! unpack message:
    IF (iroot /= this_pe) THEN
      position = 0
      nitems   = 1
      CALL p_unpack_int(buffer, position, nitems, icomm)
      IF (nitems == 0) THEN
        broadcast_event_data = .FALSE.
      ELSE
        broadcast_event_data = .TRUE.
        CALL unpack_metadata(buffer, position, evd, icomm)
      END IF
    ELSE
      broadcast_event_data = l_no_end_message
    END IF

    ! clean up
    DEALLOCATE(buffer, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END FUNCTION broadcast_event_data


  !> Utility routine: Unpack MPI message buffer with event meta-data.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE pack_metadata(buffer, position, evd, icomm)
    CHARACTER,                             INTENT(INOUT)  :: buffer(:)             !< MPI buffer for packed
    INTEGER,                               INTENT(INOUT)  :: position              !< MPI buffer position
    TYPE(t_event_data_local),              INTENT(IN)     :: evd
    INTEGER,                               INTENT(IN)     :: icomm                !< MPI communicator
    ! local variables
    INTEGER :: i

    ! encode event name string
    CALL p_pack_string(TRIM(evd%name),                           buffer, position, icomm)
    DO i=1,MAX_TIME_INTERVALS
      ! encode event begin and end string
      CALL p_pack_string(TRIM(evd%begin_str(i)),                 buffer, position, icomm)
      CALL p_pack_string(TRIM(evd%end_str(i)),                   buffer, position, icomm)
      ! encode event interval string
      CALL p_pack_string(TRIM(evd%intvl_str(i)),                 buffer, position, icomm)
    END DO
    ! encode flag "evd%l_output_last":
    CALL p_pack_bool(evd%l_output_last,                          buffer, position, icomm)
    ! encode t_sim_step_info data
    CALL p_pack_string(TRIM(evd%sim_step_info%sim_start),        buffer, position, icomm)
    CALL p_pack_string(TRIM(evd%sim_step_info%sim_end),          buffer, position, icomm)
    CALL p_pack_real(evd%sim_step_info%dtime,                    buffer, position, icomm)
    CALL p_pack_string(TRIM(evd%sim_step_info%run_start),        buffer, position, icomm)
    CALL p_pack_string(TRIM(evd%sim_step_info%restart_time),     buffer, position, icomm)
    CALL p_pack_int(evd%sim_step_info%jstep0,                    buffer, position, icomm)
    CALL p_pack_string(evd%sim_step_info%dom_start_time,         buffer, position, icomm)
    CALL p_pack_string(evd%sim_step_info%dom_end_time,           buffer, position, icomm)
    ! encode fname_metadata data
    CALL p_pack_int(evd%fname_metadata%steps_per_file,           buffer, position, icomm)
    CALL p_pack_bool(evd%fname_metadata%steps_per_file_inclfirst,buffer, position, icomm)
    CALL p_pack_string(TRIM(evd%fname_metadata%file_interval),   buffer, position, icomm)
    CALL p_pack_int(evd%fname_metadata%phys_patch_id,            buffer, position, icomm)
    CALL p_pack_int(evd%fname_metadata%ilev_type,                buffer, position, icomm)
    CALL p_pack_string(TRIM(evd%fname_metadata%filename_format), buffer, position, icomm)
    CALL p_pack_string(TRIM(evd%fname_metadata%filename_pref),   buffer, position, icomm)
    CALL p_pack_string(TRIM(evd%fname_metadata%extn),            buffer, position, icomm)
    CALL p_pack_int(evd%fname_metadata%jfile_offset,             buffer, position, icomm)
    CALL p_pack_int(evd%fname_metadata%npartitions,              buffer, position, icomm)
    CALL p_pack_int(evd%fname_metadata%ifile_partition,          buffer, position, icomm)
    ! encode this event's MPI tag
    CALL p_pack_int(evd%i_tag,                                   buffer, position, icomm)
  END SUBROUTINE pack_metadata


  !> Utility routine: Unpack MPI message buffer with event meta-data.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE unpack_metadata(buffer, position, evd, icomm)
    CHARACTER,                             INTENT(INOUT)  :: buffer(:)             !< MPI buffer for packed
    INTEGER,                               INTENT(INOUT)  :: position              !< MPI buffer position
    TYPE(t_event_data_local),              INTENT(INOUT) :: evd
    INTEGER,                               INTENT(IN)     :: icomm                !< MPI communicator
    ! local variables
    INTEGER :: i

    ! decode event name string
    CALL p_unpack_string(buffer, position, evd%name,                                     icomm)
    DO i=1,(MAX_TIME_INTERVALS)
      ! decode event begin and end string                                                            
      CALL p_unpack_string(buffer, position, evd%begin_str(i),                           icomm)
      CALL p_unpack_string(buffer, position, evd%end_str(i),                             icomm)
      ! decode event interval string                                                           
      CALL p_unpack_string(buffer, position, evd%intvl_str(i),                           icomm)
    END DO
    ! decode flag "l_output_last":                                                                 
    CALL p_unpack_bool(  buffer, position, evd%l_output_last,                            icomm)
    ! decode t_sim_step_info data                                                                  
    CALL p_unpack_string(buffer, position, evd%sim_step_info%sim_start,                  icomm)
    CALL p_unpack_string(buffer, position, evd%sim_step_info%sim_end,                    icomm)
    CALL p_unpack_real(  buffer, position, evd%sim_step_info%dtime,                      icomm)
    CALL p_unpack_string(buffer, position, evd%sim_step_info%run_start,                  icomm)
    CALL p_unpack_string(buffer, position, evd%sim_step_info%restart_time,               icomm)
    CALL p_unpack_int(   buffer, position, evd%sim_step_info%jstep0,                     icomm)
    CALL p_unpack_string(buffer, position, evd%sim_step_info%dom_start_time,             icomm)
    CALL p_unpack_string(buffer, position, evd%sim_step_info%dom_end_time,               icomm)
    ! decode fname_metadata data
    CALL p_unpack_int(   buffer, position, evd%fname_metadata%steps_per_file,            icomm)
    CALL p_unpack_bool(  buffer, position, evd%fname_metadata%steps_per_file_inclfirst,  icomm)
    CALL p_unpack_string(buffer, position, evd%fname_metadata%file_interval,             icomm)
    CALL p_unpack_int(   buffer, position, evd%fname_metadata%phys_patch_id,             icomm)
    CALL p_unpack_int(   buffer, position, evd%fname_metadata%ilev_type,                 icomm)
    CALL p_unpack_string(buffer, position, evd%fname_metadata%filename_format,           icomm)
    CALL p_unpack_string(buffer, position, evd%fname_metadata%filename_pref,             icomm)
    CALL p_unpack_string(buffer, position, evd%fname_metadata%extn,                      icomm)
    CALL p_unpack_int(   buffer, position, evd%fname_metadata%jfile_offset,              icomm)
    CALL p_unpack_int(   buffer, position, evd%fname_metadata%npartitions,               icomm)
    CALL p_unpack_int(   buffer, position, evd%fname_metadata%ifile_partition,           icomm)
    ! decode this event's MPI tag                                                                 
    CALL p_unpack_int(   buffer, position, evd%i_tag,                                    icomm)
  END SUBROUTINE unpack_metadata


  !---------------------------------------------------------------
  ! ROUTINES PERFORMING HANDSHAKE DURING EVENT STEPS
  !---------------------------------------------------------------

  !> Declare output step as being completed by this participating PE.
  !
  !  @note We do not use the wrapper routines from "mo_mpi" here,
  !        since we need direct control over the non-blocking P2P
  !        requests.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE pass_output_step(event)
    TYPE(t_par_output_event), POINTER :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::pass_output_step"
    LOGICAL :: step_is_not_finished
#ifndef NOMPI
    INTEGER :: ierrstat, impi_status(MPI_STATUS_SIZE), istep, i_tag
#endif

    step_is_not_finished = .NOT. is_output_event_finished(event)

#ifndef NOMPI
    IF (use_async_name_list_io .AND. step_is_not_finished &
         .AND. (event%icomm /= MPI_COMM_NULL)) THEN
      ! wait for the last ISEND to be processed:
      IF (ldebug) WRITE (0,*) p_pe, ": waiting for request handle."
      CALL MPI_WAIT(event%isend_req, impi_status, ierrstat)
      IF (ierrstat /= mpi_success) CALL finish (routine, 'Error in MPI_WAIT.')
      IF (ldebug) WRITE (0,*) p_pe, ": waiting for request handle done."
      ! launch a new non-blocking send:
      istep = event%output_event%i_event_step
      i_tag = event%output_event%event_step(istep)%event_step_data(1)%i_tag
      event%isend_buf = istep
      IF (ldebug) THEN
        WRITE (0,*) routine, ": sending message ", i_tag, " from ", get_my_global_mpi_id(), &
          &         " to ", event%iroot
      END IF
      CALL MPI_IBSEND(event%isend_buf, 1, p_int, event%iroot, i_tag, &
        &            event%icomm, event%isend_req, ierrstat)
      IF (ierrstat /= mpi_success) CALL finish (routine, 'Error in MPI_ISEND.')
      IF (ldebug) THEN
        WRITE (0,*) "pass ", event%output_event%i_event_step, &
          &         " (",  event%output_event%event_step(istep)%event_step_data(1)%jfile, &
          &         " / ", event%output_event%event_step(istep)%event_step_data(1)%jpart, &
          &         " )"
      END IF
    END IF
#endif

    IF (step_is_not_finished) THEN
      ! increment step counter
      event%output_event%i_event_step = event%output_event%i_event_step + 1
    END IF
  END SUBROUTINE pass_output_step


#ifndef NOMPI
  !> Place MPI_IRECV calls, thereby asking all participating PEs for
  !> acknowledging the completion of the output step.
  !
  !  @note We do not use the wrapper routines from "mo_mpi" here,
  !        since we need direct control over the non-blocking P2P
  !        requests.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE trigger_output_step_irecv(event)
    TYPE(t_par_output_event), POINTER :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::trigger_output_step_irecv"
    INTEGER                     :: ierrstat, cur_step, i, i_pe, i_tag
    TYPE(t_event_step), POINTER :: event_step

    IF (.NOT. ASSOCIATED(event))  RETURN
    ! determine this PE's MPI rank wrt. the given MPI communicator and
    ! return if this PE is not root PE for the given event:
    IF (.NOT. is_event_root_pe(event))  RETURN

    ! wait for the last IRECVs to be processed:
    CALL wait_for_pending_irecvs(event)

    IF (ALLOCATED(event%irecv_req)) THEN
      DEALLOCATE(event%irecv_req, event%irecv_buf, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
      event%irecv_nreq = 0
    END IF
    IF (.NOT. is_output_event_finished(event)) THEN
      ! launch a couple of new non-blocking receives:
      cur_step         =  event%output_event%i_event_step
      event_step       => event%output_event%event_step(cur_step)
      event%irecv_nreq =  event_step%n_pes

      ALLOCATE(event%irecv_req(event%irecv_nreq), event%irecv_buf(event%irecv_nreq), STAT=ierrstat)
      event%irecv_req(:) = MPI_REQUEST_NULL
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      DO i=1,event%irecv_nreq
        i_pe  = event_step%event_step_data(i)%i_pe
        i_tag = event_step%event_step_data(i)%i_tag
        IF (ldebug) THEN
          WRITE (0,*) routine, ": launching IRECV ", i_tag, " on ", get_my_global_mpi_id(), &
            &         " to ", i_pe
        END IF
        CALL MPI_IRECV(event%irecv_buf(i), 1, p_int, i_pe, &
          &            i_tag, event%icomm, event%irecv_req(i), ierrstat)
        IF (ierrstat /= mpi_success) CALL finish (routine, 'Error in MPI_IRECV.')
      END DO
    END IF
  END SUBROUTINE trigger_output_step_irecv
#endif


#ifndef NOMPI
  !> Utility routine: blocking wait for pending "output completed"
  !> messages.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE wait_for_pending_irecvs(event)
    TYPE(t_par_output_event), POINTER :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::wait_for_pending_irecs"
    INTEGER                     :: ierrstat
    ! INTEGER, ALLOCATABLE        :: irecv_status(:,:)

    ! wait for the last IRECVs to be processed:
    IF (ALLOCATED(event%irecv_req)) THEN
      ! ALLOCATE(irecv_status(MPI_STATUS_SIZE,event%irecv_nreq), STAT=ierrstat)
      ! IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      CALL MPI_WAITALL(event%irecv_nreq, event%irecv_req, mpi_statuses_ignore, ierrstat)
      IF (ierrstat /= mpi_success) CALL finish(routine, 'Error in MPI_WAITALL.')
      ! DEALLOCATE(irecv_status, STAT=ierrstat)
      ! IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
      ! increment event step counter:
      event%output_event%i_event_step = event%output_event%i_event_step + 1
    END IF
  END SUBROUTINE wait_for_pending_irecvs
#endif


  !> Utility routine: blocking wait for pending "output completed"
  !> messages.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE blocking_wait_for_irecvs(event)
    TYPE(t_par_output_event), POINTER :: event
    ! local variables
#ifndef NOMPI
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::blocking_wait_for_irecvs"
    INTEGER                           :: ierrstat, nreq, ireq
    INTEGER, ALLOCATABLE              :: irecv_req(:)
    TYPE(t_par_output_event), POINTER :: ev

    ! count the number of request handles
    ev   =>  event
    nreq = 0
    DO WHILE (ASSOCIATED(ev))
      IF (ALLOCATED(ev%irecv_req)) nreq = nreq + ev%irecv_nreq
      ev => ev%next
    END DO
    IF (ldebug) THEN
      WRITE (0,*) "Total ", nreq, " IRECV request handles."
    END IF
    IF (nreq == 0) RETURN

    ! collect the request handles
    ALLOCATE(irecv_req(nreq), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ev   =>  event
    ireq = 1
    DO WHILE (ASSOCIATED(ev))
      IF (ALLOCATED(ev%irecv_req)) &
        irecv_req(ireq:(ireq+ev%irecv_nreq-1)) = ev%irecv_req(1:ev%irecv_nreq)
      ireq = ireq + ev%irecv_nreq
      ev => ev%next
    END DO

    ! wait for the last IRECVs to be processed:
    CALL MPI_WAITALL(nreq, irecv_req, mpi_statuses_ignore, ierrstat)
    IF (ierrstat /= mpi_success) CALL finish (routine, 'Error in MPI_WAITALL.')
    DEALLOCATE(irecv_req, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    ! clear the request handles
    ev => event
    DO WHILE (ASSOCIATED(ev))
      IF (ALLOCATED(ev%irecv_req)) &
        ev%irecv_req(:) = MPI_REQUEST_NULL
      ev => ev%next
    END DO
#endif
  END SUBROUTINE blocking_wait_for_irecvs


  !> Spool event state fast-forward to a given event step.
  !
  !  This functionality is, e.g., required for resume after restart.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE set_event_to_simstep(event, jstep, l_isrestart, lrecover_open_file)
    TYPE(t_output_event),   INTENT(INOUT), TARGET :: event              !< output event data structure
    INTEGER,                INTENT(IN)            :: jstep              !< simulation step
    LOGICAL,                INTENT(IN)            :: l_isrestart        !< .TRUE. if this is a restart run
    !> Flag. If true, we test for an existing file from previous runs
    LOGICAL,                INTENT(IN)            :: lrecover_open_file
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_event_to_simstep"
    CHARACTER(LEN=10) :: jpart_str
    INTEGER :: ev_step, istep, n_pes, i_pe
    TYPE(t_event_step_data), POINTER :: event_step_data

    ev_step = 0
    DO istep=1,event%n_event_steps
      IF (event%event_step(istep)%i_sim_step < jstep)  ev_step = istep
    END DO
    event%i_event_step = ev_step + 1

    IF (event%i_event_step <= event%n_event_steps) THEN
      istep = event%i_event_step
      n_pes = event%event_step(istep)%n_pes

      DO i_pe=1,n_pes
        event_step_data => event%event_step(istep)%event_step_data(i_pe)

        ! Spooling forward the event status to a given step means that
        ! we have to set the "open" flag.
        event_step_data%l_open_file = .TRUE.

        ! Test if the opening of the file at step "jstep" would
        ! destroy a file of the same name that has been created in a
        ! previous run (i.e. we are in a restart run).
        !
        ! This happens in the following two situations:
        !
        ! - Each file contains more than one output step and "jstep"
        !   corresponds to an output step "in the middle of the file".
        !
        ! - We are in a restart run and we are dealing with the first
        !   output file, which also contains the initial state.

        IF ((event_step_data%jpart > 1) .OR.         &
          & (l_isrestart .AND. (event_step_data%jfile == 1))) THEN
          ! Resuming after a restart means that we have to open the
          ! file for output though this has not been planned
          ! initially. We must find a unique suffix then for this new
          ! file (otherwise we would overwrite the file from the last
          ! step) and we must inform the user about this incident.
          IF (.NOT. lrecover_open_file) THEN
            ! simply throw an error message
            CALL finish(routine, "Attempt to overwrite existing file after restart!")
          ELSE
            ! otherwise: modify file name s.t. the new, resumed file
            ! is clearly distinguishable: We append "_<part>+"
            jpart_str = int2string(event_step_data%jpart)
            IF (my_process_is_mpi_workroot()) THEN
              WRITE (0,*) "Modify filename ", TRIM(event_step_data%filename_string), " to ", &
                &      TRIM(event_step_data%filename_string)//"_part_"//TRIM(jpart_str)//"+",  &
                &      " after restart."
            END IF
            CALL modify_filename(event, TRIM(event_step_data%filename_string), &
              &       TRIM(event_step_data%filename_string)//"_part_"//TRIM(jpart_str)//"+", &
              &       start_step=istep)
          END IF
        END IF

      END DO
    END IF
  END SUBROUTINE set_event_to_simstep


  !> Spool event state fast-forward to a given event step.
  !  Implementation for parallel events.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE set_event_to_simstep_par(event, jstep, l_isrestart, lrecover_open_file)
    TYPE(t_par_output_event), TARGET    :: event              !< output event data structure
    INTEGER,                  INTENT(IN) :: jstep              !< simulation step
    LOGICAL,                  INTENT(IN) :: l_isrestart        !< .TRUE. if this is a restart run
    LOGICAL,                  INTENT(IN) :: lrecover_open_file !< Flag. If true, we test for an existing file from previous runs
    TYPE(t_par_output_event), POINTER :: ev
    ev => event
    DO WHILE (ASSOCIATED(ev))
      CALL set_event_to_simstep(ev%output_event, jstep, l_isrestart, lrecover_open_file)
      ev => ev%next
    END DO
  END SUBROUTINE set_event_to_simstep_par


  !> Utility routine: Loop over all event steps of a given event and
  !> modify all matching filenames.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE modify_filename(event, old_name, new_name, start_step)
    TYPE(t_output_event), INTENT(INOUT), TARGET :: event              !< output event data structure    
    CHARACTER(LEN=*),     INTENT(IN)            :: old_name, new_name !< old file name string and replacement
    INTEGER,              INTENT(IN)            :: start_step         !< event step where to start searching
    ! local variables
    INTEGER :: istep, n_pes, i_pe, ipart
    TYPE(t_event_step_data), POINTER :: event_step_data

    ipart = 0
    DO istep=start_step,event%n_event_steps
      n_pes = event%event_step(istep)%n_pes
      DO i_pe=1,n_pes
        event_step_data => event%event_step(istep)%event_step_data(i_pe)        
        IF (TRIM(event_step_data%filename_string) == TRIM(old_name)) THEN
          event_step_data%filename_string = TRIM(new_name)
          ipart = ipart + 1
          event_step_data%jpart = ipart
        END IF
      END DO
    END DO
  END SUBROUTINE modify_filename


  !> Utility routine: Strip date-time stamp (string) from modifiers,
  !  e.g. ">", "<".
  !
  !  @author F. Prill, DWD
  !
  FUNCTION strip_from_modifiers(dt_string)
    CHARACTER(LEN=*), INTENT(IN) :: dt_string
    CHARACTER(LEN=LEN_TRIM(dt_string)) :: strip_from_modifiers
    ! local variables
    CHARACTER :: char
    CHARACTER(LEN=LEN(dt_string)) :: str

    str  = TRIM(remove_whitespace(TRIM(dt_string)))
    char = str(1:1)
    SELECT CASE(char)
    CASE ('>','<')
      str = str(2:)
    END SELECT
    strip_from_modifiers = TRIM(str)
  END FUNCTION strip_from_modifiers

END MODULE mo_output_event_handler
