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
!! -----------------------------------------------------------------------------------
MODULE mo_output_event_handler

  ! actual method (MPI-2)
#ifndef NOMPI
#if !defined (__SUNPRO_F95)
  USE mpi
#endif
#endif

  USE mo_impl_constants,         ONLY: SUCCESS
  USE mo_exception,              ONLY: finish, message
  USE mo_kind,                   ONLY: wp
  USE mo_io_units,               ONLY: FILENAME_MAX, find_next_free_unit
  USE mo_util_string,            ONLY: int2string
  USE mo_mpi,                    ONLY: p_int, p_pe_work, p_comm_work,                       &
    &                                  p_pack_int, p_pack_string, p_pack_bool, p_pack_real, &
    &                                  p_unpack_int, p_unpack_string, p_unpack_bool,        &
    &                                  p_unpack_real, p_send_packed, p_irecv_packed, p_wait,&
    &                                  p_bcast, get_my_global_mpi_id
  USE mo_fortran_tools,          ONLY: assign_if_present
  USE mtime,                     ONLY: MAX_DATETIME_STR_LEN,                                &
    &                                  MAX_TIMEDELTA_STR_LEN, PROLEPTIC_GREGORIAN,          &
    &                                  event, datetime, newEvent, timedelta,                &
    &                                  setCalendar, resetCalendar, newTimedelta,            &
    &                                  deallocateDatetime, datetimeToString,                &
    &                                  deallocateEvent, newDatetime, OPERATOR(>=),          &
    &                                  OPERATOR(>), OPERATOR(+), deallocateTimedelta
  USE mo_mtime_extensions,       ONLY: isCurrentEventActive, getPTStringFromMS
  USE mo_output_event_types,     ONLY: t_sim_step_info, t_event_data, t_event_step_data,    &
    &                                  t_event_step, t_output_event, t_par_output_event,    &
    &                                  MAX_FILENAME_STR_LEN, MAX_EVENT_NAME_STR_LEN,        &
    &                                  DEFAULT_EVENT_NAME
  USE mo_output_event_control,   ONLY: compute_matching_sim_steps, generate_output_filenames
  USE mo_name_list_output_types, ONLY: t_fname_metadata
  USE mo_util_table,             ONLY: initialize_table, finalize_table, add_table_column,  &
    &                                  set_table_entry, print_table, t_table

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
  PUBLIC :: wait_for_pending_irecvs
  ! inquiry functions
  PUBLIC :: get_current_date
  PUBLIC :: get_current_step
  PUBLIC :: get_current_filename
  PUBLIC :: check_write_readyfile
  PUBLIC :: check_open_file
  PUBLIC :: check_close_file
  PUBLIC :: is_event_root_pe
  PUBLIC :: is_output_step
  PUBLIC :: is_output_step_complete
  PUBLIC :: is_output_event_finished
  ! handshake routines (called after steps)
  PUBLIC :: pass_output_step
  PUBLIC :: trigger_output_step_irecv
  ! auxiliary functions
  PUBLIC :: print_output_event
  PUBLIC :: set_event_to_simstep


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
  INTEGER, PARAMETER :: MAX_BUF_SIZE =    (    MAX_EVENT_NAME_STR_LEN &
    &                                      + 7*MAX_DATETIME_STR_LEN   &
    &                                      + MAX_TIMEDELTA_STR_LEN    &
    &                                      + 3*FILENAME_MAX           &
    &                                      + 1024 )

  !> MPI message tag for setup of output events
  INTEGER, PARAMETER :: SENDRECV_TAG_SETUP    = 1002

  !> MPI message tag for output event handshake
  INTEGER, PARAMETER :: SENDRECV_TAG_OUTEVENT = 1001

  !> MPI rank of root PE handling ready files
  INTEGER, PARAMETER :: ROOT_OUTEVENT         = 0

  !> Internal switch for debugging output
  LOGICAL, PARAMETER :: ldebug  = .FALSE.


#ifdef NOMPI
  !---------------------------------------------------------------
  ! non-MPI runs: local list event with event meta-data
  !---------------------------------------------------------------

  !> event meta-data: this is only required for non-MPI runs, where we have to
  !  keep a local list of output events.
  TYPE t_event_data_nompi
    CHARACTER(LEN=MAX_EVENT_NAME_STR_LEN)        :: name                 !< output event name
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)          :: begin_str, end_str
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)         :: intvl_str
    LOGICAL                                      :: l_output_last        !< Flag. If .TRUE. the last step is always written
    TYPE(t_sim_step_info)                        :: sim_step_info        !< definitions for conversion "time stamp -> simulation step"
    TYPE(t_fname_metadata)                       :: fname_metadata       !< additional meta-data for generating output filename
    INTEGER                                      :: icomm                !< MPI communicator
    INTEGER                                      :: dst_rank             !< MPI destination rank
  END TYPE t_event_data_nompi

  !> Maximum length of local event meta-data list: only required for non-MPI
  !  runs
  INTEGER, PARAMETER :: NOMPI_NMAX_EVENT_LIST = 100

  !> local list of output events: only required for non-MPI runs
  TYPE(t_event_data_nompi) :: event_list_nompi(NOMPI_NMAX_EVENT_LIST)

  !> length of local list of output events: only required for non-MPI runs
  INTEGER :: ievent_list_nompi = 0
#endif

 
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
  RECURSIVE SUBROUTINE deallocate_par_output_event(event)
    TYPE(t_par_output_event), POINTER :: event
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::deallocate_par_output_event"
    INTEGER :: ierrstat

    IF (.NOT. ASSOCIATED(event)) RETURN
    IF (ASSOCIATED(event%next)) THEN
      CALL deallocate_output_event(event%next)
      NULLIFY(event%next)
    END IF

    CALL deallocate_output_event(event%output_event)
    IF (ALLOCATED(event%irecv_req)) THEN
      DEALLOCATE(event%output_event, event%irecv_req, event%irecv_buf, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    END IF
    DEALLOCATE(event, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
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

    IF (check_write_readyfile(event)) THEN
      WRITE (dst,'(a,a,a)') 'output "', TRIM(event%event_data%name), '", writes ready files:'
    ELSE
      WRITE (dst,'(a,a,a)') 'output "', TRIM(event%event_data%name), '", does not write ready files:'
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
        CALL set_table_entry(table,irow,"filename",    TRIM(event_step%event_step_data(j)%filename_string))
        CALL set_table_entry(table,irow,"I/O PE",      int2string(event_step%event_step_data(j)%i_pe))
        CALL set_table_entry(table,irow,"output date", TRIM(event_step%event_step_data(j)%datetime_string))
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
  RECURSIVE SUBROUTINE print_par_output_event(event, opt_filename, opt_dstfile)
    TYPE(t_par_output_event), POINTER             :: event
    CHARACTER(LEN=*), OPTIONAL,        INTENT(IN) :: opt_filename    !< name of ASCII file (optional)
    INTEGER,          OPTIONAL,        INTENT(IN) :: opt_dstfile     !< optional destination ASCII file unit
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::print_par_output_event"
    INTEGER :: dst, ierrstat

    ! consistency check:
    IF (PRESENT(opt_dstfile) .AND. PRESENT(opt_filename)) THEN
      CALL finish(routine, "Routine was called with both a file unit and a filename to open!")
    END IF

    ! open ASCII output file (if necessary):
    dst = 0
    IF (PRESENT(opt_dstfile)) THEN
      dst = opt_dstfile
    ELSE
      IF (PRESENT(opt_filename)) THEN
        dst = find_next_free_unit(10,100)
        OPEN (dst, file=TRIM(opt_filename), status='REPLACE', iostat=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish (routine, 'OPEN failed.')
      END IF
    END IF

    WRITE (dst,*) " " ! newline
    CALL print_output_event(event%output_event, opt_dstfile=dst)
    IF (ASSOCIATED(event%next)) THEN
      CALL print_par_output_event(event%next, opt_dstfile=dst)
    END IF
    WRITE (dst,*) " " ! newline

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
  !  Note that this subroutine generates the output event steps but
  !  does not map these time stamps onto the corresponding simulation
  !  steps. This simulation-specific task is performed by a different
  !  subroutine, which is given as a function parameter
  !  "fct_time2simstep".
  !
  !  @author F. Prill, DWD
  !
  FUNCTION new_output_event(name, i_pe, i_tag, begin_str, end_str, intvl_str, l_output_last, sim_step_info,  &
    &                       fname_metadata, fct_time2simstep, fct_generate_filenames) RESULT(p_event)
    TYPE(t_output_event),   POINTER :: p_event
    CHARACTER(LEN=*),       INTENT(IN)  :: name                 !< output event name
    INTEGER,                INTENT(IN)  :: i_pe                 !< rank of participating PE
    INTEGER,                INTENT(IN)  :: i_tag                !< tag, e.g. for MPI isend/irecv messages
    CHARACTER(len=*),       INTENT(IN)  :: begin_str, end_str
    CHARACTER(len=*),       INTENT(IN)  :: intvl_str
    LOGICAL,                INTENT(IN)  :: l_output_last        !< Flag. If .TRUE. the last step is always written
    TYPE(t_sim_step_info),  INTENT(IN)  :: sim_step_info        !< definitions for conversion "time stamp -> simulation step"
    TYPE(t_fname_metadata), INTENT(IN)  :: fname_metadata       !< additional meta-data for generating output filename

    !> As an argument of this function, the user must provide a
    !  conversion "time stamp -> simulation step"
    INTERFACE
      SUBROUTINE fct_time2simstep(nstrings, date_string, sim_step_info, &
        &                         result_steps, result_exactdate)
        USE mo_output_event_types, ONLY: t_sim_step_info
        INTEGER,              INTENT(IN)    :: nstrings             !< no. of string to convert
        CHARACTER(len=*),     INTENT(IN)    :: date_string(:)       !< array of ISO 8601 time stamp strings
        TYPE(t_sim_step_info),INTENT(IN)    :: sim_step_info        !< definitions: time step size, etc.
        INTEGER,              INTENT(INOUT) :: result_steps(:)      !< resulting step indices
        CHARACTER(LEN=*),     INTENT(INOUT) :: result_exactdate(:)  !< resulting (exact) time step strings 
      END SUBROUTINE fct_time2simstep
    END INTERFACE

    !> As an argument of this function, the user must provide a
    !  function for generating output file names
    INTERFACE
      FUNCTION fct_generate_filenames(nstrings, date_string, sim_steps, &
        &                             sim_step_info, fname_metadata)  RESULT(result_fnames)
        USE mo_output_event_types,     ONLY: t_sim_step_info, t_event_step_data
        USE mo_name_list_output_types, ONLY: t_fname_metadata

        INTEGER,                INTENT(IN)    :: nstrings           !< no. of string to convert
        CHARACTER(len=*),       INTENT(IN)    :: date_string(:)     !< array of ISO 8601 time stamp strings
        INTEGER,                INTENT(IN)    :: sim_steps(:)       !< array of corresponding simulation steps
        TYPE(t_sim_step_info),  INTENT(IN)    :: sim_step_info      !< definitions: time step size, etc.
        TYPE(t_fname_metadata), INTENT(IN)    :: fname_metadata     !< additional meta-data for generating output filename
        TYPE(t_event_step_data) :: result_fnames(SIZE(date_string))
      END FUNCTION fct_generate_filenames
    END INTERFACE

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::new_output_event"
    !> Max. no. of event steps (used for local array sizes)
    INTEGER, PARAMETER :: MAX_NEVENT_STEPS = 1000

    TYPE(datetime),  POINTER :: mtime_date, mtime_begin, mtime_end, sim_end, &
      &                         mtime_dom_start
    TYPE(event),     POINTER :: mtime_event
    TYPE(timedelta), POINTER :: delta
    INTEGER                  :: ierrstat, i, n_event_steps  
    LOGICAL                  :: l_active
    CHARACTER(len=MAX_DATETIME_STR_LEN)  :: mtime_date_string(MAX_NEVENT_STEPS)
    INTEGER                              :: mtime_sim_steps(MAX_NEVENT_STEPS)
    CHARACTER(len=MAX_DATETIME_STR_LEN)  :: mtime_exactdate(MAX_NEVENT_STEPS)
    TYPE(t_event_step_data)              :: filename_metadata(MAX_NEVENT_STEPS)

    ! allocate event data structure
    ALLOCATE(p_event, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    n_event_steps                    = 0
    p_event%i_event_step             = 1
    p_event%event_data%name          = TRIM(name)

    ! initialize event with the "mtime" library:
    CALL setCalendar(PROLEPTIC_GREGORIAN)
    mtime_event => newEvent(TRIM(name),                  &  ! name
      &                     TRIM(begin_str),             &  ! ref.date
      &                     TRIM(begin_str),             &  ! first
      &                     TRIM(end_str),               &  ! last
      &                     TRIM(intvl_str))                ! interval

    ! status output
    IF (ldebug) THEN
      WRITE (0,*) 'Defining output event "'//TRIM(begin_str)//'", "'//TRIM(end_str)//'", "'//TRIM(intvl_str)//'"'
      WRITE (0,*) 'Simulation bounds:    "'//TRIM(sim_step_info%sim_start)//'", "'//TRIM(sim_step_info%sim_end)//'"'
    END IF

    ! set some dates used later:
    sim_end     => newDatetime(TRIM(sim_step_info%sim_end))
    mtime_begin => newDatetime(TRIM(begin_str))
    mtime_end   => newDatetime(TRIM(end_str))
    ! Domains (and their output) can be activated and deactivated
    ! during the simulation. This is determined by the parameters
    ! "dom_start_time" and "dom_end_time". Therefore, we must create
    ! a corresponding event.
    mtime_dom_start => newDatetime(TRIM(sim_step_info%dom_start_time))

    ! loop over the event occurrences    
    mtime_date => mtime_begin
    delta      => newTimedelta(TRIM(intvl_str)) ! create a time delta 
    n_event_steps = 0
    EVENT_LOOP: DO
      IF ((mtime_date >= mtime_begin) .AND. &
        & (mtime_date >= mtime_dom_start)) THEN
        n_event_steps = n_event_steps + 1
        IF (n_event_steps > SIZE(mtime_date_string)) THEN
          CALL finish(routine, "Internal error: step buffer size exceeded!")
        END IF
        CALL datetimeToString(mtime_date, mtime_date_string(n_event_steps))
        IF (ldebug) THEN
          WRITE (0,*) "Adding step ", n_event_steps, ": ", &
            &         mtime_date_string(n_event_steps)
        END IF
      END IF
      mtime_date = mtime_date + delta
      l_active = isCurrentEventActive(mtime_event, mtime_date)
      IF (.NOT. l_active) EXIT EVENT_LOOP
      ! Optional: Append the last event time step
      IF (l_output_last .AND. .NOT. (mtime_date >= mtime_end) .AND. (mtime_date >= sim_end)) THEN
        n_event_steps = n_event_steps + 1
        IF (n_event_steps > SIZE(mtime_date_string)) THEN
          CALL finish(routine, "Internal error: step buffer size exceeded!")
        END IF
        CALL datetimeToString(sim_end, mtime_date_string(n_event_steps))
        IF (ldebug) THEN
          WRITE (0,*) n_event_steps, ": output event '", mtime_date_string(n_event_steps), "'"
        END IF
        EXIT EVENT_LOOP
      END IF
    END DO EVENT_LOOP

    CALL fct_time2simstep(n_event_steps, mtime_date_string,                            &
      &                   sim_step_info, mtime_sim_steps, mtime_exactdate)

    ! remove all those event steps which have no corresponding simulation
    ! step (mtime_sim_steps(i) < 0):
    do i=1,n_event_steps
       if (mtime_sim_steps(i) < 0)  EXIT
    end do
    n_event_steps = (i-1)

    filename_metadata = fct_generate_filenames(n_event_steps, mtime_date_string,       &
      &                   mtime_sim_steps, sim_step_info, fname_metadata)

    ! from this list of time stamp strings: generate the event steps
    ! for this event
    p_event%n_event_steps = n_event_steps
    ALLOCATE(p_event%event_step(p_event%n_event_steps), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    DO i=1,p_event%n_event_steps
      p_event%event_step(i)%exact_date_string  = TRIM(mtime_exactdate(i))
      p_event%event_step(i)%i_sim_step         = mtime_sim_steps(i)
      p_event%event_step(i)%n_pes              = 1
      ALLOCATE(p_event%event_step(i)%event_step_data(p_event%event_step(i)%n_pes), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      p_event%event_step(i)%event_step_data(1)%i_pe            = i_pe
      p_event%event_step(i)%event_step_data(1)%datetime_string = TRIM(mtime_date_string(i))
      p_event%event_step(i)%event_step_data(1)%i_tag           = i_tag
      p_event%event_step(i)%event_step_data(1)%filename_string = TRIM(filename_metadata(i)%filename_string)
      p_event%event_step(i)%event_step_data(1)%l_open_file     = filename_metadata(i)%l_open_file 
      p_event%event_step(i)%event_step_data(1)%l_close_file    = filename_metadata(i)%l_close_file
    END DO

    ! clean up
    CALL deallocateEvent(mtime_event)
    CALL deallocateDatetime(mtime_begin)
    CALL deallocateDatetime(mtime_end)
    CALL deallocateDatetime(mtime_dom_start)
    CALL deallocateDatetime(sim_end)
    CALL deallocateTimedelta(delta)
    CALL resetCalendar()
  END FUNCTION new_output_event


  !> Create a simple *parallel* output event, happening at regular intervals.
  ! 
  !  This subroutine calls the local version "new_output_event", adds
  !  data structures for parallel communication and launches a
  !  non-blocking MPI send to the root I/O PE.
  !
  !  @author F. Prill, DWD
  !
  FUNCTION new_parallel_output_event(name, begin_str, end_str, intvl_str, l_output_last, sim_step_info, &
    &                                fname_metadata, fct_time2simstep, fct_generate_filenames,          &
    &                                local_event_no, icomm) RESULT(p_event)
    TYPE(t_par_output_event),   POINTER :: p_event
    CHARACTER(LEN=*),       INTENT(IN)  :: name                 !< output event name
    CHARACTER(len=*),       INTENT(IN)  :: begin_str, end_str
    CHARACTER(len=*),       INTENT(IN)  :: intvl_str
    LOGICAL,                INTENT(IN)  :: l_output_last        !< Flag. If .TRUE. the last step is always written
    TYPE(t_sim_step_info),  INTENT(IN)  :: sim_step_info        !< definitions for conversion "time stamp -> simulation step"
    TYPE(t_fname_metadata), INTENT(IN)  :: fname_metadata       !< additional meta-data for generating output filename
    INTEGER,                INTENT(IN)  :: local_event_no       !< local index of this event on local PE
    INTEGER,                INTENT(IN)  :: icomm                !< MPI communicator

 
    !> As an argument of this function, the user must provide a
    !  conversion "time stamp -> simulation step"
    INTERFACE
      SUBROUTINE fct_time2simstep(nstrings, date_string, sim_step_info, &
        &                         result_steps, result_exactdate)
        USE mo_output_event_types, ONLY: t_sim_step_info
        INTEGER,              INTENT(IN)    :: nstrings             !< no. of string to convert
        CHARACTER(len=*),     INTENT(IN)    :: date_string(:)       !< array of ISO 8601 time stamp strings
        TYPE(t_sim_step_info),INTENT(IN)    :: sim_step_info        !< definitions: time step size, etc.
        INTEGER,              INTENT(INOUT) :: result_steps(:)      !< resulting step indices
        CHARACTER(LEN=*),     INTENT(INOUT) :: result_exactdate(:)  !< resulting (exact) time step strings 
      END SUBROUTINE fct_time2simstep
    END INTERFACE

    !> As an argument of this function, the user must provide a
    !  function for generating output file names
    INTERFACE
      FUNCTION fct_generate_filenames(nstrings, date_string, sim_steps, &
        &                             sim_step_info, fname_metadata)  RESULT(result_fnames)
        USE mo_output_event_types,     ONLY: t_sim_step_info, t_event_step_data
        USE mo_name_list_output_types, ONLY: t_fname_metadata

        INTEGER,                INTENT(IN)    :: nstrings           !< no. of string to convert
        CHARACTER(len=*),       INTENT(IN)    :: date_string(:)     !< array of ISO 8601 time stamp strings
        INTEGER,                INTENT(IN)    :: sim_steps(:)       !< array of corresponding simulation steps
        TYPE(t_sim_step_info),  INTENT(IN)    :: sim_step_info      !< definitions: time step size, etc.
        TYPE(t_fname_metadata), INTENT(IN)    :: fname_metadata     !< additional meta-data for generating output filename
        TYPE(t_event_step_data) :: result_fnames(SIZE(date_string))
      END FUNCTION fct_generate_filenames
    END INTERFACE

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::new_parallel_output_event"
    INTEGER :: ierrstat, this_pe, i_tag, nranks

    ! determine this PE's MPI rank wrt. the given MPI communicator:
#ifndef NOMPI
    IF (icomm /= MPI_COMM_NULL) THEN
      CALL MPI_COMM_RANK(icomm, this_pe, ierrstat)
      IF (ierrstat /= 0) CALL finish (routine, 'Error in MPI_COMM_RANK.')
      CALL MPI_COMM_SIZE (icomm, nranks, ierrstat)
      IF (ierrstat /= 0) CALL finish (routine, 'Error in MPI_COMM_SIZE.')
      IF (ldebug) THEN
         write (0,*) "PE ",get_my_global_mpi_id(), ": local rank is ", this_pe, "; icomm has size ", nranks
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

    ! create the local (non-parallel) output event data structure:
    p_event%output_event => new_output_event(name, this_pe, i_tag , begin_str, end_str, intvl_str, l_output_last, &
      &                                      sim_step_info, fname_metadata,                                       &
      &                                      fct_time2simstep, fct_generate_filenames)

    ! set the other MPI-related data fields:
    p_event%icomm         = icomm
    p_event%iroot         = ROOT_OUTEVENT
    p_event%irecv_nreq    = 0
    p_event%isend_req     = 0
#ifndef NOMPI
    p_event%isend_req     = MPI_REQUEST_NULL
    ! now send all the meta-data to the root PE
    IF (icomm /= MPI_COMM_NULL) THEN
      if (ldebug) then
         write (0,*) "new_parallel_output_event: PE ", get_my_global_mpi_id(), ": send event data: ", &
              &      trim(begin_str), trim(end_str), trim(intvl_str)
      end if
      CALL send_event_data(name, begin_str, end_str, intvl_str, l_output_last, sim_step_info, &
        &                  fname_metadata, icomm, ROOT_OUTEVENT)
    END IF
#else
    ! non-MPI runs: we keep a local list of event meta-data
    ievent_list_nompi = ievent_list_nompi + 1
    event_list_nompi(ievent_list_nompi)%name            = name
    event_list_nompi(ievent_list_nompi)%begin_str       = begin_str
    event_list_nompi(ievent_list_nompi)%end_str         = end_str
    event_list_nompi(ievent_list_nompi)%intvl_str       = intvl_str
    event_list_nompi(ievent_list_nompi)%l_output_last   = l_output_last
    event_list_nompi(ievent_list_nompi)%sim_step_info   = sim_step_info
    event_list_nompi(ievent_list_nompi)%fname_metadata  = fname_metadata
    event_list_nompi(ievent_list_nompi)%icomm           = icomm
    event_list_nompi(ievent_list_nompi)%dst_rank        = ROOT_OUTEVENT
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

    !> As an argument of this function, the user must provide a
    !  conversion "time stamp -> simulation step"
    INTERFACE
      SUBROUTINE fct_time2simstep(nstrings, date_string, sim_step_info, &
        &                         result_steps, result_exactdate)
        USE mo_output_event_types, ONLY: t_sim_step_info
        INTEGER,              INTENT(IN)    :: nstrings             !< no. of string to convert
        CHARACTER(len=*),     INTENT(IN)    :: date_string(:)       !< array of ISO 8601 time stamp strings
        TYPE(t_sim_step_info),INTENT(IN)    :: sim_step_info        !< definitions: time step size, etc.
        INTEGER,              INTENT(INOUT) :: result_steps(:)      !< resulting step indices
        CHARACTER(LEN=*),     INTENT(INOUT) :: result_exactdate(:)  !< resulting (exact) time step strings 
      END SUBROUTINE fct_time2simstep
    END INTERFACE

    !> As an argument of this function, the user must provide a
    !  function for generating output file names
    INTERFACE
      FUNCTION fct_generate_filenames(nstrings, date_string, sim_steps, &
        &                             sim_step_info, fname_metadata)  RESULT(result_fnames)
        USE mo_output_event_types,     ONLY: t_sim_step_info, t_event_step_data
        USE mo_name_list_output_types, ONLY: t_fname_metadata

        INTEGER,                INTENT(IN)    :: nstrings           !< no. of string to convert
        CHARACTER(len=*),       INTENT(IN)    :: date_string(:)     !< array of ISO 8601 time stamp strings
        INTEGER,                INTENT(IN)    :: sim_steps(:)       !< array of corresponding simulation steps
        TYPE(t_sim_step_info),  INTENT(IN)    :: sim_step_info      !< definitions: time step size, etc.
        TYPE(t_fname_metadata), INTENT(IN)    :: fname_metadata     !< additional meta-data for generating output filename
        TYPE(t_event_step_data) :: result_fnames(SIZE(date_string))
      END FUNCTION fct_generate_filenames
    END INTERFACE

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::union_of_all_events"
    TYPE(t_par_output_event), POINTER     :: par_event, last_node
    TYPE(t_output_event),     POINTER     :: ev1, ev2
    CHARACTER(LEN=MAX_EVENT_NAME_STR_LEN) :: recv_name                 !< output event name
    CHARACTER(len=MAX_DATETIME_STR_LEN)   :: recv_begin_str, recv_end_str
    CHARACTER(len=MAX_DATETIME_STR_LEN)   :: recv_intvl_str
    LOGICAL                               :: lrecv
    LOGICAL                               :: recv_l_output_last        !< Flag. If .TRUE. the last step is always written
    TYPE(t_sim_step_info)                 :: recv_sim_step_info        !< definitions for conversion "time stamp -> simulation step"
    TYPE(t_fname_metadata)                :: recv_fname_metadata       !< additional meta-data for generating output filename
    INTEGER                               :: i_pe, nranks, ierrstat, &
      &                                      i_tag, this_pe, nbcast_ranks
    LOGICAL                               :: lbroadcast

    lbroadcast = PRESENT(opt_broadcast_root) .AND. &
         &       PRESENT(opt_broadcast_comm)

    NULLIFY(union_of_all_events)
    ! get the number of ranks in this MPI communicator
    nranks  =  1
    this_pe = -1
#ifndef NOMPI
    IF (icomm /= MPI_COMM_NULL) THEN
      CALL MPI_COMM_RANK(icomm, this_pe, ierrstat)
      IF (ierrstat /= 0) CALL finish (routine, 'Error in MPI_COMM_RANK.')
      CALL MPI_COMM_SIZE (icomm, nranks, ierrstat)
      IF (ierrstat /= 0) CALL finish (routine, 'Error in MPI_COMM_SIZE.')
      IF (ldebug) THEN
         write (0,*) "PE ",get_my_global_mpi_id(), ": local rank is ", this_pe, "; icomm has size ", nranks
         if (lbroadcast) then
            CALL MPI_COMM_SIZE (opt_broadcast_comm, nbcast_ranks, ierrstat)
            IF (ierrstat /= 0) CALL finish (routine, 'Error in MPI_COMM_SIZE.')
            write (0,*) "PE ",get_my_global_mpi_id(), ": local rank is ", this_pe, "; bcast comm has size ", nbcast_ranks, &
                 &      ", root is ", opt_broadcast_root
         end if
      END IF
    END IF
#endif
    IF (lbroadcast)  CALL p_bcast(nranks, opt_broadcast_root, opt_broadcast_comm)

    ! loop over all PEs of the I/O communicator:
    DO i_pe = 0,(nranks-1)
      i_tag = SENDRECV_TAG_OUTEVENT + i_pe
      lrecv = .TRUE.
      RECEIVE_LOOP : DO
#ifndef NOMPI
        ! receive event meta-data from participating I/O PEs:
        IF (this_pe == ROOT_OUTEVENT) THEN
          lrecv = receive_event_data(recv_name, recv_begin_str, recv_end_str, recv_intvl_str, recv_l_output_last, &
            &                        recv_sim_step_info, recv_fname_metadata, icomm, i_pe, SENDRECV_TAG_SETUP)
        END IF
        ! forward the event meta-data to the worker PEs
        IF (lbroadcast) THEN
          lrecv =  broadcast_event_data(recv_name, recv_begin_str, recv_end_str,                    &
            &                           recv_intvl_str, recv_l_output_last, recv_sim_step_info,     &
            &                           recv_fname_metadata, opt_broadcast_comm,                    &
            &                           opt_broadcast_root, lrecv)          
        END IF
#else
        ! non-MPI runs: we keep a local list of event meta-data
        lrecv = (ievent_list_nompi > 0)
        IF (ievent_list_nompi > 0) THEN
          recv_name           = event_list_nompi(ievent_list_nompi)%name
          recv_begin_str      = event_list_nompi(ievent_list_nompi)%begin_str
          recv_end_str        = event_list_nompi(ievent_list_nompi)%end_str
          recv_intvl_str      = event_list_nompi(ievent_list_nompi)%intvl_str
          recv_l_output_last  = event_list_nompi(ievent_list_nompi)%l_output_last
          recv_sim_step_info  = event_list_nompi(ievent_list_nompi)%sim_step_info
          recv_fname_metadata = event_list_nompi(ievent_list_nompi)%fname_metadata
          ievent_list_nompi = ievent_list_nompi - 1
          IF (ldebug) THEN
            WRITE (0,*) "Taking event ", TRIM(recv_begin_str), " / ", TRIM(recv_end_str), " / ", TRIM(recv_intvl_str), &
              &         " from local list."
          END IF
        END IF
#endif

        IF (.NOT. lrecv) EXIT RECEIVE_LOOP

        ! create the event steps from the received meta-data:
        ev2 => new_output_event(TRIM(recv_name), i_pe, i_tag, TRIM(recv_begin_str),              &
          &                     TRIM(recv_end_str), TRIM(recv_intvl_str),                        &
          &                     recv_l_output_last, recv_sim_step_info, recv_fname_metadata,     &
          &                     fct_time2simstep, fct_generate_filenames)

        ! find the parallel output event in the linked list that
        ! matches the event name
        NULLIFY(last_node)
        par_event => union_of_all_events
        DO
          IF (.NOT. ASSOCIATED(par_event)) EXIT
          IF (TRIM(par_event%output_event%event_data%name) == TRIM(ev2%event_data%name)) EXIT
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
        i_tag = i_tag + nranks 
        IF (ASSOCIATED(par_event%output_event)) THEN
          ev1 => par_event%output_event
          ! create union of the two events:
          par_event%output_event => event_union(ev1,ev2)
          CALL deallocate_output_event(ev1)
          CALL deallocate_output_event(ev2)
        ELSE
          par_event%output_event => ev2
        END IF
      END DO RECEIVE_LOOP
    END DO
  END FUNCTION union_of_all_events


  !> Create a new output event from two joint events.
  !
  !  Different output events/PEs may be joined together to form a
  !  single new event. Then not every PE necessarily participates in
  !  every event step. The purpose of joining/grouping events is that
  !  certain actions can be triggered for the group, e.g. writing
  !  "ready files".
  !
  !  @author F. Prill, DWD
  !
  FUNCTION event_union(event1, event2) RESULT(p_event)
    TYPE(t_output_event), POINTER :: p_event
    TYPE(t_output_event), INTENT(IN) :: event1, event2
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::event_union"
    INTEGER              :: i1, i2, ierrstat, i, i_sim_step1, i_sim_step2, &
      &                     max_sim_step

    ! allocate event data structure
    ALLOCATE(p_event, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    p_event%n_event_steps            = 0
    p_event%i_event_step             = 1
    IF (TRIM(event1%event_data%name) == TRIM(event2%event_data%name)) THEN
      p_event%event_data%name        = TRIM(event1%event_data%name)
    ELSE
      p_event%event_data%name        = TRIM(event1%event_data%name)//","//TRIM(event2%event_data%name)
    END IF

    ! loop over the two events, determine the size of the union:
    i2 = 1
    DO i1=1,event1%n_event_steps
      p_event%n_event_steps = p_event%n_event_steps + 1
      ! find step in second event matching i1:
      STEP_LOOP1 : DO
        IF (i2 > event2%n_event_steps)  EXIT STEP_LOOP1
        IF (event2%event_step(i2)%i_sim_step > event1%event_step(i1)%i_sim_step)  EXIT STEP_LOOP1
        ! avoid double-counting joint event steps:
        IF (.NOT. (event1%event_step(i1)%i_sim_step == event2%event_step(i2)%i_sim_step)) THEN
          p_event%n_event_steps = p_event%n_event_steps + 1
        END IF
        i2 = i2 + 1
      END DO STEP_LOOP1
    END DO 

    ! choose a maximum simulation step number as abort criterion:
    max_sim_step = 0
    IF (event1%n_event_steps > 0)  max_sim_step = MAX(max_sim_step, event1%event_step(event1%n_event_steps)%i_sim_step)
    IF (event2%n_event_steps > 0)  max_sim_step = MAX(max_sim_step, event2%event_step(event2%n_event_steps)%i_sim_step)
    max_sim_step = max_sim_step + 1

    i  = 0
    i1 = 1
    i2 = 1
    DO
      IF ((i1 > event1%n_event_steps) .AND. (i2 > event2%n_event_steps))  EXIT
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
      IF (i_sim_step2 > i_sim_step1) THEN
        ! copy event step i1 from event1:
        i = i+1
        i1 = i1 + 1
      ELSE IF (i_sim_step2 < i_sim_step1) THEN
        ! copy event step i2 from event2:
        i = i+1
        i2 = i2 + 1
      ELSE
        i = i+1
        ! join event steps:
        i1 = i1 + 1
        i2 = i2 + 1
      END IF
    END DO
    p_event%n_event_steps = i

    ! now create the new event:
    ALLOCATE(p_event%event_step(p_event%n_event_steps), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    i  = 0
    i1 = 1
    i2 = 1
    DO
      IF ((i1 > event1%n_event_steps) .AND. (i2 > event2%n_event_steps))  EXIT
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
      IF (i_sim_step2 > i_sim_step1) THEN
        ! copy event step i1 from event1:
        i = i+1
        CALL append_event_step(p_event%event_step(i), event1%event_step(i1), l_create=.TRUE.)
        i1 = i1 + 1
      ELSE IF (i_sim_step2 < i_sim_step1) THEN
        ! copy event step i2 from event2:
        i = i+1
        CALL append_event_step(p_event%event_step(i), event2%event_step(i2), l_create=.TRUE.)
        i2 = i2 + 1
      ELSE
        i = i+1
        ! join event steps:
        CALL append_event_step(p_event%event_step(i), event1%event_step(i1), l_create=.TRUE.)
        CALL append_event_step(p_event%event_step(i), event2%event_step(i2), l_create=.FALSE.)
        i1 = i1 + 1
        i2 = i2 + 1
      END IF
    END DO
  END FUNCTION event_union


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
    INTEGER :: iold, ierrstat, i1, i2

    ! event_step%i_sim_step
    IF (l_create) THEN
      dst_event_step%i_sim_step = src_event_step%i_sim_step
    ELSE
      IF (src_event_step%i_sim_step /= dst_event_step%i_sim_step) &
        &   CALL finish(routine, "sim_steps do not match!")
    END IF
    ! event_step%exact_date_string
    IF (l_create) THEN
      dst_event_step%exact_date_string = src_event_step%exact_date_string
    ELSE
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
      ALLOCATE(tmp_event_step_data(dst_event_step%n_pes), STAT=ierrstat)
      iold = dst_event_step%n_pes
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      tmp_event_step_data(1:dst_event_step%n_pes) = &
        &   dst_event_step%event_step_data(1:dst_event_step%n_pes)
      DEALLOCATE(dst_event_step%event_step_data, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

      dst_event_step%n_pes = dst_event_step%n_pes + src_event_step%n_pes
      ALLOCATE(dst_event_step%event_step_data(dst_event_step%n_pes), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      dst_event_step%event_step_data(1:iold) = tmp_event_step_data(1:iold)
      dst_event_step%event_step_data((iold+1):(iold+src_event_step%n_pes)) = &
        &   src_event_step%event_step_data(1:src_event_step%n_pes)

      ! consistency check: test, if any of the filenames in group1
      ! occurs in group2 (avoid duplicate names)
      DO i1=1,iold
         DO i2=(iold+1),(iold+src_event_step%n_pes)
            IF (dst_event_step%event_step_data(i1)%filename_string == dst_event_step%event_step_data(i2)%filename_string) THEN
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

    if (is_output_event_finished(event)) then
       is_output_step = .false.
    else
       istep = event%i_event_step
       is_output_step = (event%event_step(istep)%i_sim_step == jstep)
    end if
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
    TYPE(t_par_output_event), POINTER                :: event
    ! local variables
#ifndef NOMPI
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::is_output_step_complete"
    INTEGER              :: ierrstat
    INTEGER, ALLOCATABLE :: irecv_status(:,:)

    ret = .TRUE.
    IF (event%irecv_nreq == 0)  RETURN
    ALLOCATE(irecv_status(MPI_STATUS_SIZE,event%irecv_nreq), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    CALL MPI_TESTALL(event%irecv_nreq, event%irecv_req, ret, &
      &              irecv_status, ierrstat)
    IF (ierrstat /= 0) CALL finish (routine, 'Error in MPI_TESTALL.')
    DEALLOCATE(irecv_status, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
#else
    ret = .TRUE.
#endif
  END FUNCTION is_output_step_complete


  !> @return current date-time stamp string.
  !  @author F. Prill, DWD
  !
  FUNCTION get_current_date(event)
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: get_current_date
    TYPE(t_output_event), POINTER :: event
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


  !> @return .TRUE. if this PE should write a ready file for the given
  !          event.
  !  @author F. Prill, DWD
  !
  FUNCTION check_write_readyfile(event)
    LOGICAL :: check_write_readyfile
    TYPE(t_output_event), POINTER :: event

    IF (ASSOCIATED(event)) THEN
      check_write_readyfile = (TRIM(event%event_data%name) /= DEFAULT_EVENT_NAME)
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
    INTEGER :: ierrstat, this_pe, p_comm
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
    CALL MPI_COMM_RANK(p_comm, this_pe, ierrstat)
    IF (ierrstat /= 0) CALL finish (routine, 'Error in MPI_COMM_RANK.')
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
  SUBROUTINE send_event_data(name, begin_str, end_str, intvl_str, l_output_last, sim_step_info, &
    &                        fname_metadata, icomm, dst_rank)
    CHARACTER(LEN=*),       INTENT(IN)  :: name                 !< output event name
    CHARACTER(len=*),       INTENT(IN)  :: begin_str, end_str
    CHARACTER(len=*),       INTENT(IN)  :: intvl_str
    LOGICAL,                INTENT(IN)  :: l_output_last        !< Flag. If .TRUE. the last step is always written
    TYPE(t_sim_step_info),  INTENT(IN)  :: sim_step_info        !< definitions for conversion "time stamp -> simulation step"
    TYPE(t_fname_metadata), INTENT(IN)  :: fname_metadata       !< additional meta-data for generating output filename
    INTEGER,                INTENT(IN)  :: icomm                !< MPI communicator
    INTEGER,                INTENT(IN)  :: dst_rank             !< MPI destination rank
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
    CALL p_pack_int(nitems,         buffer, MAX_BUF_SIZE, position, icomm)
    CALL pack_metadata(buffer, position, name, begin_str, end_str, intvl_str, l_output_last, &
      &                sim_step_info, fname_metadata, icomm)

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
    INTEGER :: nitems, ierrstat, position
    CHARACTER, ALLOCATABLE :: buffer(:)   !< MPI buffer for packed

    IF (icomm /= MPI_COMM_NULL) THEN
      ! allocate message buffer
      ALLOCATE(buffer(MAX_BUF_SIZE), stat=ierrstat)  
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed')
      position = 0
      ! prepare an empty MPI message:
      nitems = 0
      CALL p_pack_int(nitems, buffer, MAX_BUF_SIZE, position, icomm)

      ! send packed message:
      CALL p_send_packed(buffer, ROOT_OUTEVENT, SENDRECV_TAG_SETUP, position, icomm)
      ! clean up
      DEALLOCATE(buffer, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')    
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
  FUNCTION receive_event_data(name, begin_str, end_str, intvl_str, l_output_last, sim_step_info, &
    &                         fname_metadata, icomm, isrc, isendrecv_tag)
    LOGICAL :: receive_event_data
    CHARACTER(LEN=MAX_EVENT_NAME_STR_LEN), INTENT(INOUT) :: name                 !< output event name
    CHARACTER(len=MAX_DATETIME_STR_LEN)  , INTENT(INOUT) :: begin_str, end_str
    CHARACTER(len=MAX_DATETIME_STR_LEN)  , INTENT(INOUT) :: intvl_str
    LOGICAL,                INTENT(INOUT)  :: l_output_last        !< Flag. If .TRUE. the last step is always written
    TYPE(t_sim_step_info),  INTENT(INOUT)  :: sim_step_info        !< definitions for conversion "time stamp -> simulation step"
    TYPE(t_fname_metadata), INTENT(INOUT)  :: fname_metadata       !< additional meta-data for generating output filename
    INTEGER,                INTENT(IN)     :: icomm                !< MPI communicator
    INTEGER,                INTENT(IN)     :: isrc                 !< MPI rank of the sending PE
    INTEGER,                INTENT(IN)     :: isendrecv_tag        !< MPI tag for this messages isend/irecv communication
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::receive_event_data"
    INTEGER :: nitems, ierrstat, position
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
    CALL p_unpack_int(buffer, MAX_BUF_SIZE, position, nitems, icomm)
    IF (nitems == 0) THEN
      receive_event_data = .FALSE.
    ELSE
      receive_event_data = .TRUE.
      CALL unpack_metadata(buffer, position, name, begin_str, end_str, intvl_str, l_output_last, &
        &                  sim_step_info, fname_metadata, icomm)
    END IF

    if (ldebug) then
       write (0,*) "received event data: ", trim(begin_str), trim(end_str), trim(intvl_str)
    end if

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
  FUNCTION broadcast_event_data(name, begin_str, end_str, intvl_str, l_output_last, sim_step_info, &
    &                           fname_metadata, icomm, iroot, l_no_end_message)
    LOGICAL :: broadcast_event_data
    CHARACTER(LEN=MAX_EVENT_NAME_STR_LEN), INTENT(INOUT) :: name                 !< output event name
    CHARACTER(len=MAX_DATETIME_STR_LEN)  , INTENT(INOUT) :: begin_str, end_str
    CHARACTER(len=MAX_DATETIME_STR_LEN)  , INTENT(INOUT) :: intvl_str
    LOGICAL,                INTENT(INOUT)  :: l_output_last        !< Flag. If .TRUE. the last step is always written
    TYPE(t_sim_step_info),  INTENT(INOUT)  :: sim_step_info        !< definitions for conversion "time stamp -> simulation step"
    TYPE(t_fname_metadata), INTENT(INOUT)  :: fname_metadata       !< additional meta-data for generating output filename
    INTEGER,                INTENT(IN)     :: icomm                !< MPI communicator
    INTEGER,                INTENT(IN)     :: iroot                !< MPI broadcast root rank
    LOGICAL,                INTENT(IN)     :: l_no_end_message     !< Flag. .FALSE. if "end message" shall be broadcasted
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::broadcast_event_data"
    INTEGER :: nitems, ierrstat, position, this_pe
    CHARACTER, ALLOCATABLE :: buffer(:)   !< MPI buffer for packed

    ! allocate message buffer
    ALLOCATE(buffer(MAX_BUF_SIZE), stat=ierrstat)  
    IF (ierrstat /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed')

    ! determine this PE's MPI rank wrt. the given MPI communicator:
    this_pe = -1
#ifndef NOMPI
    IF (icomm /= MPI_COMM_NULL) THEN
      CALL MPI_COMM_RANK(icomm, this_pe, ierrstat)
      IF (ierrstat /= 0) CALL finish (routine, 'Error in MPI_COMM_RANK.')
    END IF
#endif

    ! prepare an MPI message:
    IF (iroot == this_pe) THEN
      position = 0
      IF (l_no_end_message) THEN
        ! normal message:
        nitems = 1
        CALL p_pack_int(nitems,         buffer, MAX_BUF_SIZE, position, icomm)
        CALL pack_metadata(buffer, position, name, begin_str, end_str, intvl_str, l_output_last, &
          &                sim_step_info, fname_metadata, icomm)
      ELSE
        ! empty, "end message":
        nitems = 0
        CALL p_pack_int(nitems, buffer, MAX_BUF_SIZE, position, icomm)
      END IF
      if (ldebug) then
         write (0,*) "PE ", get_my_global_mpi_id(), ": send event data: ", trim(begin_str), trim(end_str), trim(intvl_str)
      end if
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
      CALL p_unpack_int(buffer, MAX_BUF_SIZE, position, nitems, icomm)
      IF (nitems == 0) THEN
        broadcast_event_data = .FALSE.
      ELSE
        broadcast_event_data = .TRUE.
        CALL unpack_metadata(buffer, position, name, begin_str, end_str, intvl_str, l_output_last, &
          &                  sim_step_info, fname_metadata, icomm)
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
  SUBROUTINE pack_metadata(buffer, position, name, begin_str, end_str, intvl_str, l_output_last, &
    &                      sim_step_info, fname_metadata, icomm)
    CHARACTER,                             INTENT(INOUT)  :: buffer(:)             !< MPI buffer for packed
    INTEGER,                               INTENT(INOUT)  :: position              !< MPI buffer position
    CHARACTER(LEN=MAX_EVENT_NAME_STR_LEN), INTENT(IN)     :: name                  !< output event name
    CHARACTER(len=MAX_DATETIME_STR_LEN)  , INTENT(IN)     :: begin_str, end_str
    CHARACTER(len=MAX_DATETIME_STR_LEN)  , INTENT(IN)     :: intvl_str
    LOGICAL,                               INTENT(IN)     :: l_output_last        !< Flag. If .TRUE. the last step is always written
    TYPE(t_sim_step_info),                 INTENT(IN)     :: sim_step_info        !< definitions for conversion "time stamp -> simulation step"
    TYPE(t_fname_metadata),                INTENT(IN)     :: fname_metadata       !< additional meta-data for generating output filename
    INTEGER,                               INTENT(IN)     :: icomm                !< MPI communicator

    ! encode event name string
    CALL p_pack_string(TRIM(name),                           buffer, MAX_BUF_SIZE, position, icomm)
    ! encode event begin and end string
    CALL p_pack_string(TRIM(begin_str),                      buffer, MAX_BUF_SIZE, position, icomm)
    CALL p_pack_string(TRIM(end_str),                        buffer, MAX_BUF_SIZE, position, icomm)
    ! encode event interval string
    CALL p_pack_string(TRIM(intvl_str),                      buffer, MAX_BUF_SIZE, position, icomm)
    ! encode flag "l_output_last":
    CALL p_pack_bool(l_output_last,                          buffer, MAX_BUF_SIZE, position, icomm)
    ! encode t_sim_step_info data
    CALL p_pack_string(TRIM(sim_step_info%sim_start),        buffer, MAX_BUF_SIZE, position, icomm)
    CALL p_pack_string(TRIM(sim_step_info%sim_end),          buffer, MAX_BUF_SIZE, position, icomm)
    CALL p_pack_real(sim_step_info%dtime,                    buffer, MAX_BUF_SIZE, position, icomm)
    CALL p_pack_int(sim_step_info%iadv_rcf,                  buffer, MAX_BUF_SIZE, position, icomm)
    CALL p_pack_int(sim_step_info%jstep0,                    buffer, MAX_BUF_SIZE, position, icomm)
    CALL p_pack_string(sim_step_info%dom_start_time,         buffer, MAX_BUF_SIZE, position, icomm)
    CALL p_pack_string(sim_step_info%dom_end_time,           buffer, MAX_BUF_SIZE, position, icomm)
    ! encode fname_metadata data
    CALL p_pack_int(fname_metadata%steps_per_file,           buffer, MAX_BUF_SIZE, position, icomm)
    CALL p_pack_string(TRIM(fname_metadata%file_interval),   buffer, MAX_BUF_SIZE, position, icomm)
    CALL p_pack_int(fname_metadata%phys_patch_id,            buffer, MAX_BUF_SIZE, position, icomm)
    CALL p_pack_int(fname_metadata%ilev_type,                buffer, MAX_BUF_SIZE, position, icomm)
    CALL p_pack_string(TRIM(fname_metadata%filename_format), buffer, MAX_BUF_SIZE, position, icomm)
    CALL p_pack_string(TRIM(fname_metadata%filename_pref),   buffer, MAX_BUF_SIZE, position, icomm)
    CALL p_pack_string(TRIM(fname_metadata%extn),            buffer, MAX_BUF_SIZE, position, icomm)
  END SUBROUTINE pack_metadata


  !> Utility routine: Unpack MPI message buffer with event meta-data.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE unpack_metadata(buffer, position, name, begin_str, end_str, intvl_str, l_output_last, &
    &                        sim_step_info, fname_metadata, icomm)
    CHARACTER,                             INTENT(INOUT)  :: buffer(:)             !< MPI buffer for packed
    INTEGER,                               INTENT(INOUT)  :: position              !< MPI buffer position
    CHARACTER(LEN=MAX_EVENT_NAME_STR_LEN), INTENT(INOUT)  :: name                  !< output event name
    CHARACTER(len=MAX_DATETIME_STR_LEN)  , INTENT(INOUT)  :: begin_str, end_str
    CHARACTER(len=MAX_DATETIME_STR_LEN)  , INTENT(INOUT)  :: intvl_str
    LOGICAL,                               INTENT(INOUT)  :: l_output_last        !< Flag. If .TRUE. the last step is always written
    TYPE(t_sim_step_info),                 INTENT(INOUT)  :: sim_step_info        !< definitions for conversion "time stamp -> simulation step"
    TYPE(t_fname_metadata),                INTENT(INOUT)  :: fname_metadata       !< additional meta-data for generating output filename
    INTEGER,                               INTENT(IN)     :: icomm                !< MPI communicator

    ! decode event name string
    CALL p_unpack_string(buffer, MAX_BUF_SIZE, position, name,                           icomm)
    ! decode event begin and end string
    CALL p_unpack_string(buffer, MAX_BUF_SIZE, position, begin_str,                      icomm)
    CALL p_unpack_string(buffer, MAX_BUF_SIZE, position, end_str,                        icomm)
    ! decode event interval string
    CALL p_unpack_string(buffer, MAX_BUF_SIZE, position, intvl_str,                      icomm)
    ! decode flag "l_output_last":
    CALL p_unpack_bool(  buffer, MAX_BUF_SIZE, position, l_output_last,                  icomm)
    ! decode t_sim_step_info data
    CALL p_unpack_string(buffer, MAX_BUF_SIZE, position, sim_step_info%sim_start,        icomm)
    CALL p_unpack_string(buffer, MAX_BUF_SIZE, position, sim_step_info%sim_end,          icomm)
    CALL p_unpack_real(  buffer, MAX_BUF_SIZE, position, sim_step_info%dtime,            icomm)
    CALL p_unpack_int(   buffer, MAX_BUF_SIZE, position, sim_step_info%iadv_rcf,         icomm)
    CALL p_unpack_int(   buffer, MAX_BUF_SIZE, position, sim_step_info%jstep0,           icomm)
    CALL p_unpack_string(buffer, MAX_BUF_SIZE, position, sim_step_info%dom_start_time,   icomm)
    CALL p_unpack_string(buffer, MAX_BUF_SIZE, position, sim_step_info%dom_end_time,     icomm)
    ! decode fname_metadata data
    CALL p_unpack_int(   buffer, MAX_BUF_SIZE, position, fname_metadata%steps_per_file,  icomm)
    CALL p_unpack_string(buffer, MAX_BUF_SIZE, position, fname_metadata%file_interval,   icomm)
    CALL p_unpack_int(   buffer, MAX_BUF_SIZE, position, fname_metadata%phys_patch_id,   icomm)
    CALL p_unpack_int(   buffer, MAX_BUF_SIZE, position, fname_metadata%ilev_type,       icomm)
    CALL p_unpack_string(buffer, MAX_BUF_SIZE, position, fname_metadata%filename_format, icomm)
    CALL p_unpack_string(buffer, MAX_BUF_SIZE, position, fname_metadata%filename_pref,   icomm)
    CALL p_unpack_string(buffer, MAX_BUF_SIZE, position, fname_metadata%extn,            icomm)
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
#ifndef NOMPI
    INTEGER :: ierrstat, impi_status(MPI_STATUS_SIZE), istep, i_tag

    IF (.NOT. is_output_event_finished(event) .AND. &
      & (event%icomm /= MPI_COMM_NULL)) THEN
      ! wait for the last ISEND to be processed:
      CALL MPI_WAIT(event%isend_req, impi_status, ierrstat)
      IF (ierrstat /= 0) CALL finish (routine, 'Error in MPI_WAIT.')
      ! launch a new non-blocking send:
      istep = event%output_event%i_event_step
      i_tag = event%output_event%event_step(istep)%event_step_data(1)%i_tag
      event%isend_buf = istep
      CALL MPI_ISEND(event%isend_buf, 1, p_int, event%iroot, i_tag, &
        &            event%icomm, event%isend_req, ierrstat)
      IF (ierrstat /= 0) CALL finish (routine, 'Error in MPI_ISEND.')
    END IF
#endif
    IF (.NOT. is_output_event_finished(event)) THEN
      ! increment step counter
      event%output_event%i_event_step = event%output_event%i_event_step + 1
    END IF
  END SUBROUTINE pass_output_step


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
#ifndef NOMPI
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
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      DO i=1,event%irecv_nreq
        i_pe  = event_step%event_step_data(i)%i_pe
        i_tag = event_step%event_step_data(i)%i_tag
        CALL MPI_IRECV(event%irecv_buf(i), 1, p_int, i_pe, &
          &            i_tag, event%icomm, event%irecv_req(i), ierrstat)
        IF (ierrstat /= 0) CALL finish (routine, 'Error in MPI_IRECV.')
      END DO
    END IF
#else
    ! increment event step counter:
    event%output_event%i_event_step = event%output_event%i_event_step + 1
#endif
  END SUBROUTINE trigger_output_step_irecv


  !> Utility routine: blocking wait for pending "output completed"
  !> messages.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE wait_for_pending_irecvs(event)
    TYPE(t_par_output_event), POINTER :: event
    ! local variables
#ifndef NOMPI
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::wait_for_pending_irecs"    
    INTEGER                     :: ierrstat
    INTEGER, ALLOCATABLE        :: irecv_status(:,:)

    ! wait for the last IRECVs to be processed:
    IF (ALLOCATED(event%irecv_req)) THEN
      ALLOCATE(irecv_status(MPI_STATUS_SIZE,event%irecv_nreq), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      CALL MPI_WAITALL(event%irecv_nreq, event%irecv_req, irecv_status, ierrstat)
      IF (ierrstat /= 0) CALL finish (routine, 'Error in MPI_WAITALL.')
      DEALLOCATE(irecv_status, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
      ! increment event step counter:
      event%output_event%i_event_step = event%output_event%i_event_step + 1
    END IF
#endif
  END SUBROUTINE wait_for_pending_irecvs


  !> Spool event state fast-forward to a given event step.
  ! 
  !  This functionality is, e.g., required for resume after restart.
  !
  !  @author F. Prill, DWD
  !
  SUBROUTINE set_event_to_simstep(event, jstep, lforce_open_file)
    TYPE(t_output_event), INTENT(INOUT) :: event            !< output event data structure
    INTEGER,              INTENT(IN)    :: jstep            !< simulation step
    LOGICAL,              INTENT(IN)    :: lforce_open_file !< Flag. If true, the open-Flag is set for the reached event step.
    ! local variables
    INTEGER :: ev_step, istep, n_pes
   
    ev_step = 0
    DO istep=1,event%n_event_steps
      IF (event%event_step(istep)%i_sim_step <= jstep)  ev_step = istep
    END DO
    event%i_event_step = ev_step + 1
    IF ((event%i_event_step <= event%n_event_steps) .AND. &
      & lforce_open_file) THEN
      n_pes = event%event_step(event%i_event_step)%n_pes
      istep = event%i_event_step
      event%event_step(istep)%event_step_data(1:n_pes)%l_open_file = .TRUE.
    END IF
  END SUBROUTINE set_event_to_simstep


  !> Spool event state fast-forward to a given event step.
  !  Implementation for parallel events.
  !
  !  @author F. Prill, DWD
  !
  RECURSIVE SUBROUTINE set_event_to_simstep_par(event, jstep, lforce_open_file)
    TYPE(t_par_output_event), POINTER    :: event            !< output event data structure
    INTEGER,                  INTENT(IN) :: jstep            !< simulation step
    LOGICAL,                  INTENT(IN) :: lforce_open_file !< Flag. If true, the open-Flag is set for the reached event step.
   
    IF (.NOT. ASSOCIATED(event)) RETURN
    IF (ASSOCIATED(event%next)) THEN
      CALL set_event_to_simstep_par(event%next, jstep, lforce_open_file)
    END IF
    CALL set_event_to_simstep(event%output_event, jstep, lforce_open_file)
  END SUBROUTINE set_event_to_simstep_par

END MODULE mo_output_event_handler
