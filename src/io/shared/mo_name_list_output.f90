!! -------------------------------------------------------------------------
!>
!! Module handling synchronous and asynchronous output; supporting
!! multiple I/O PEs and horizontal interpolation.
!!
!! @author R. Johanni
!!
!! @par Revision History
!! Initial implementation  by  R. Johanni  (2011)
!! Major changes: F. Prill, DWD (2012-2013)
!!
!! @todo In asynchronous I/O mode, windows are created but not freed
!!
!! @todo Several fields are allocated but not freed at the end of the
!!       simulation. A pseudo-destructor should be implemented!
!!
!! @note: The spelling "name_list" (with underscore) is intended to make
!!        clear that this does not pertain to a FORTRAN namelist but rather
!!        to a list of names of output variables
!!
!! Define USE_CRAY_POINTER for platforms having problems with ISO_C_BINDING
!! BUT understand CRAY pointers
!!   #define USE_CRAY_POINTER
!!
!! -------------------------------------------------------------------------
!!
!! MPI roles in asynchronous communication:
!! 
!! - Compute PEs create local memory windows, buffering all variables
!!   for all output files (for the local horizontal grid partition).
!!
!! - Asynchronous I/O servers create trivial local memory windows of
!!   size 1.
!!
!! - Additionally, when writing, the asynchronous I/O servers allocate
!!   a 3D buffer for a single variable. This temporary field serves as
!!   a target buffer for the one-sided MPI_GET operation.
!!
!! Transfer of meta-info:
!!
!!  Since parts of the variable's "info-field" TYPE(t_var_metadata) may change
!!  during simulation, the following mechanism updates the metadata on the
!!  asynchronous output PEs:
!!  For each output file, a separate MPI window is created on work PE#0, where 
!!  the work root stores the current variable meta-info. This is then retrieved
!!  via an additional MPI_GET by the I/O PEs.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! -------------------------------------------------------------------------
MODULE mo_name_list_output

  ! constants
  USE mo_kind,                      ONLY: wp, i8, dp, sp
  USE mo_impl_constants,            ONLY: max_dom, SUCCESS, MAX_TIME_LEVELS
  USE mo_cdi_constants              ! We need all
  ! utility functions
  USE mo_io_units,                  ONLY: FILENAME_MAX, find_next_free_unit
  USE mo_exception,                 ONLY: finish, message, message_text
  USE mo_util_string,               ONLY: t_keyword_list, associate_keyword, with_keywords
  USE mo_dictionary,                ONLY: dict_finalize
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_write_output, ltimer,      &
    &                                     print_timer
  USE mo_name_list_output_gridinfo, ONLY: write_grid_info_grb2, GRID_INFO_NONE
  ! config
  USE mo_master_nml,                ONLY: model_base_dir
  USE mo_grid_config,               ONLY: n_dom
  USE mo_run_config,                ONLY: msg_level
  USE mo_io_config,                 ONLY: lkeep_in_sync
  USE mo_parallel_config,           ONLY: nproma, p_test_run, use_dp_mpi2io, num_io_procs
  USE mo_name_list_output_config,   ONLY: use_async_name_list_io
  ! data types
  USE mo_var_metadata_types,        ONLY: t_var_metadata, POST_OP_SCALE
  USE mo_name_list_output_types,    ONLY: t_output_name_list, t_output_file, t_reorder_info,        &
  &                                       msg_io_start, msg_io_done, msg_io_shutdown, all_events
  USE mo_output_event_types,        ONLY: t_sim_step_info, t_par_output_event
  ! parallelization
  USE mo_communication,             ONLY: exchange_data, t_comm_gather_pattern, idx_no, blk_no
  USE mo_mpi,                       ONLY: p_send, p_recv, p_barrier, stop_mpi,                      &
    &                                     p_mpi_wtime, p_irecv, p_wait, p_test, p_isend,            &
    &                                     p_comm_work, p_real_dp, p_real_sp, p_int,                 &
    &                                     my_process_is_stdio, my_process_is_mpi_test,              &
    &                                     my_process_is_mpi_workroot,                               &
    &                                     my_process_is_io, my_process_is_mpi_ioroot,               &
    &                                     process_mpi_all_test_id, process_mpi_all_workroot_id,     &
    &                                     num_work_procs, p_pe, p_pe_work, p_work_pe0, p_io_pe0
  ! calendar operations
  USE mtime,                        ONLY: datetime, newDatetime, deallocateDatetime,                &
    &                                     PROLEPTIC_GREGORIAN, setCalendar,                         &
    &                                     timedelta, newTimedelta, deallocateTimedelta
  USE mo_mtime_extensions,          ONLY: getTimeDeltaFromDateTime
  ! output scheduling
  USE mo_output_event_handler,      ONLY: is_output_step, check_open_file, check_close_file,        &
    &                                     pass_output_step, get_current_filename,                   &
    &                                     get_current_date, trigger_output_step_irecv,              &
    &                                     is_output_step_complete, is_output_event_finished,        &
    &                                     check_write_readyfile, blocking_wait_for_irecvs
  ! output initialization
  USE mo_name_list_output_init,     ONLY: init_name_list_output, setup_output_vlist,                &
    &                                     varnames_dict, out_varnames_dict,                         &
    &                                     output_file, patch_info, lonlat_info
  USE mo_name_list_output_metadata, ONLY: metainfo_write_to_memwin, metainfo_get_from_memwin,       &
    &                                     metainfo_get_size
  USE mo_util_cdi,                  ONLY: set_timedependent_GRIB2_keys
  ! post-ops
  USE mo_post_op,                   ONLY: perform_post_op

#ifndef __NO_ICON_ATMO__
  USE mo_dynamics_config,           ONLY: nnow, nnow_rcf, nnew, nnew_rcf
  USE mo_meteogram_output,          ONLY: meteogram_init, meteogram_finalize
  USE mo_meteogram_config,          ONLY: meteogram_output_config
  USE mo_intp_data_strc,            ONLY: lonlat_grid_list
#endif

  IMPLICIT NONE

  PRIVATE

  !-----------------------------------------------------------------
  ! include NetCDF headers (direct NetCDF library calls are required
  ! for output of grid information).
  !-----------------------------------------------------------------
  INCLUDE 'netcdf.inc'

  PUBLIC :: write_name_list_output
  PUBLIC :: close_name_list_output
  PUBLIC :: istime4name_list_output
  PUBLIC :: name_list_io_main_proc

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output'

  !> constant for better readability
  INTEGER, PARAMETER :: WAIT_UNTIL_FINISHED = -1

  !> Internal switch for debugging output
  LOGICAL, PARAMETER :: ldebug  = .FALSE.


CONTAINS


  !------------------------------------------------------------------------------------------------
  !> open_output_file:
  !  Opens a output file and sets its vlist
  !
  !  Please note that this routine is only executed on one processor
  !  (for a specific file) and thus all calls to message get the
  !  all_print=.TRUE. argument so that the messages really appear in
  !  the log.
  !
  SUBROUTINE open_output_file(of)

    TYPE(t_output_file), INTENT(INOUT) :: of
    ! local variables:
    CHARACTER(LEN=*), PARAMETER       :: routine = modname//"::open_output_file"
    CHARACTER(LEN=filename_max)       :: filename

    ! open file:
    filename = TRIM(get_current_filename(of%out_event))
    of%cdiFileID = streamOpenWrite(TRIM(filename), of%output_type)

    IF (of%cdiFileID < 0) THEN
      WRITE(message_text,'(a)') cdiStringError(of%cdiFileID)
      CALL message('',message_text,all_print=.TRUE.)
      CALL finish (routine, 'open failed on '//TRIM(filename))
    ELSE
      CALL message (routine, 'opened '//TRIM(filename),all_print=.TRUE.)
    END IF

    ! assign the vlist (which must have ben set before)
    CALL streamDefVlist(of%cdiFileID, of%cdiVlistID)

    ! set cdi internal time index to 0 for writing time slices in netCDF
    of%cdiTimeIndex = 0

  END SUBROUTINE open_output_file


  !------------------------------------------------------------------------------------------------
  !> Close all name_list files
  !
  SUBROUTINE close_name_list_output()
    ! local variables
    INTEGER :: i

#ifndef NOMPI
#ifndef __NO_ICON_ATMO__
    IF (use_async_name_list_io    .AND.  &
      & .NOT. my_process_is_io()  .AND.  &
      & .NOT. my_process_is_mpi_test()) THEN
      !-- compute PEs (senders):
      
      CALL compute_wait_for_async_io(jstep=WAIT_UNTIL_FINISHED)
      CALL compute_shutdown_async_io()

    ELSE
#endif
#endif
      !-- asynchronous I/O PEs (receiver):
      DO i = 1, SIZE(output_file)
        IF (output_file(i)%cdiFileID >= 0) THEN
          CALL close_output_file(output_file(i))
          CALL destroy_output_vlist(output_file(i))
        END IF
      ENDDO
#ifndef NOMPI
#ifndef __NO_ICON_ATMO__
    ENDIF
#endif
#endif

    DEALLOCATE(output_file)

    ! destroy variable name dictionaries:
    CALL dict_finalize(varnames_dict)
    CALL dict_finalize(out_varnames_dict)

  END SUBROUTINE close_name_list_output


  !------------------------------------------------------------------------------------------------
  !> Close output stream and the associated file,
  !  destroy all vlist related data for this file
  !
  SUBROUTINE close_output_file(of)
    TYPE (t_output_file), INTENT(INOUT) :: of
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::close_output_file"
    LOGICAL :: is_output_process

    ! GRB2 format: define geographical longitude, latitude as special
    ! variables "RLON", "RLAT"
    is_output_process = my_process_is_io() .OR. &
      &                 ((.NOT. use_async_name_list_io) .AND. my_process_is_mpi_workroot())

    IF(of%cdiFileID /= CDI_UNDEFID) THEN
#ifndef __NO_ICON_ATMO__
      IF ((of%name_list%output_grid)                                      .AND. &
        & (patch_info(of%phys_patch_id)%grid_info_mode /= GRID_INFO_NONE) .AND. &
        & is_output_process                                               .AND. &
        & (of%name_list%filetype == FILETYPE_GRB2)) THEN
        CALL write_grid_info_grb2(of, patch_info)
      END IF
#endif

      CALL streamClose(of%cdiFileID)
    END IF

    of%cdiFileID = CDI_UNDEFID
  END SUBROUTINE close_output_file


  !------------------------------------------------------------------------------------------------
  !> Close output stream and the associated file,
  !  destroy all vlist related data for this file
  !
  SUBROUTINE destroy_output_vlist(of)
    TYPE (t_output_file), INTENT(INOUT) :: of
    ! local variables
    INTEGER :: j, vlistID

    vlistID = of%cdiVlistID
    IF(vlistID /= CDI_UNDEFID) THEN
      IF(of%cdiCellGridID   /= CDI_UNDEFID) CALL gridDestroy(of%cdiCellGridID)
      IF(of%cdiEdgeGridID   /= CDI_UNDEFID) CALL gridDestroy(of%cdiEdgeGridID)
      IF(of%cdiVertGridID   /= CDI_UNDEFID) CALL gridDestroy(of%cdiVertGridID)
      IF(of%cdiLonLatGridID /= CDI_UNDEFID) CALL gridDestroy(of%cdiLonLatGridID)
      IF(of%cdiTaxisID      /= CDI_UNDEFID) CALL taxisDestroy(of%cdiTaxisID)
      DO j = 1, SIZE(of%cdiZaxisID)
        IF(of%cdiZaxisID(j) /= CDI_UNDEFID) CALL zaxisDestroy(of%cdiZaxisID(j))
      ENDDO
      CALL vlistDestroy(vlistID)
    ENDIF

    of%cdiVlistID      = CDI_UNDEFID
    of%cdiCellGridID   = CDI_UNDEFID
    of%cdiEdgeGridID   = CDI_UNDEFID
    of%cdiVertGridID   = CDI_UNDEFID
    of%cdiLonLatGridID = CDI_UNDEFID
    of%cdiTaxisID      = CDI_UNDEFID
    of%cdiZaxisID(:)   = CDI_UNDEFID

  END SUBROUTINE destroy_output_vlist


  !------------------------------------------------------------------------------------------------
  !> Loop over all output_name_list's, write the ones for which output is due
  !  This routine also cares about opening the output files the first time
  !  and reopening the files after a certain number of steps.
  !
  SUBROUTINE write_name_list_output(jstep, opt_lhas_output)
    INTEGER,           INTENT(IN)   :: jstep             !< model step
    LOGICAL, OPTIONAL, INTENT(OUT)  :: opt_lhas_output   !< (Optional) Flag: .TRUE. if this async I/O PE has written during this step.
    ! local variables
    CHARACTER(LEN=*), PARAMETER  :: routine = modname//"::write_name_list_output"
    INTEGER                           :: i, j, idate, itime, iret
    TYPE(t_output_name_list), POINTER :: p_onl
    TYPE(datetime),           POINTER :: mtime_datetime
    CHARACTER(LEN=filename_max+100)   :: text
    TYPE(t_par_output_event), POINTER :: ev
    INTEGER                           :: noutput_pe_list, list_idx
    INTEGER                           :: output_pe_list(MAX(1,num_io_procs))

    IF (ltimer) CALL timer_start(timer_write_output)
#ifndef NOMPI
#ifndef __NO_ICON_ATMO__
    IF  (use_async_name_list_io      .AND.  &
      &  .NOT. my_process_is_io()    .AND.  &
      &  .NOT. my_process_is_mpi_test()) THEN
      
      ! If asynchronous I/O is enabled, the compute PEs have to make
      ! sure that the I/O PEs are ready with the last output step
      ! before writing data into the I/O memory window.  This routine
      ! (write_name_list_output) is also called from the I/O PEs, but
      ! in this case the calling routine cares about the flow control.
      CALL compute_wait_for_async_io(jstep)

    ENDIF
#endif
#endif

    IF (PRESENT(opt_lhas_output))  opt_lhas_output = .FALSE.

    ! during the following loop, we collect a list of all I/O PEs for
    ! which output is performed:
    output_pe_list(:) = -1
    noutput_pe_list   =  0

    ! Go over all output files
    OUTFILE_LOOP : DO i=1,SIZE(output_file)

      ! Skip this output file if it is not due for output!
      IF (.NOT. is_output_step(output_file(i)%out_event, jstep))  CYCLE OUTFILE_LOOP

      ! -------------------------------------------------
      ! Check if files have to be (re)opened
      ! -------------------------------------------------

      IF (check_open_file(output_file(i)%out_event) .AND.  &
        & (output_file(i)%io_proc_id == p_pe)) THEN 
        IF (output_file(i)%cdiVlistId == CDI_UNDEFID)  &
          &  CALL setup_output_vlist(output_file(i))
        CALL open_output_file(output_file(i))
      END IF

      ! -------------------------------------------------
      ! Do the output
      ! -------------------------------------------------

      ! Notify user
#ifndef __SX__
      IF (output_file(i)%io_proc_id == p_pe) THEN
        WRITE(text,'(a,a,a,a,a,i0)') &
          & 'Output to ',TRIM(get_current_filename(output_file(i)%out_event)),        &
          & ' at simulation time ', TRIM(get_current_date(output_file(i)%out_event)), &
          & ' by PE ', p_pe
        CALL message(routine, text,all_print=.TRUE.)
      END IF
#endif
      IF (output_file(i)%io_proc_id == p_pe) THEN
        ! convert time stamp string into
        ! year/month/day/hour/minute/second values using the mtime
        ! library:
        CALL setCalendar(PROLEPTIC_GREGORIAN)
        mtime_datetime => newDatetime(TRIM(get_current_date(output_file(i)%out_event)))
        
        idate = cdiEncodeDate(INT(mtime_datetime%date%year),   &
          &                   INT(mtime_datetime%date%month),  &
          &                   INT(mtime_datetime%date%day))
        itime = cdiEncodeTime(INT(mtime_datetime%time%hour),   &
          &                   INT(mtime_datetime%time%minute), &
          &                   INT(mtime_datetime%time%second))
        CALL deallocateDatetime(mtime_datetime)
        CALL taxisDefVdate(output_file(i)%cdiTaxisID, idate)
        CALL taxisDefVtime(output_file(i)%cdiTaxisID, itime)
        iret = streamDefTimestep(output_file(i)%cdiFileId, output_file(i)%cdiTimeIndex)
        output_file(i)%cdiTimeIndex = output_file(i)%cdiTimeIndex + 1
      END IF

      p_onl => output_file(i)%name_list
      IF(my_process_is_io()) THEN
#ifndef NOMPI
        IF(output_file(i)%io_proc_id == p_pe) THEN
          CALL io_proc_write_name_list(output_file(i), check_open_file(output_file(i)%out_event))
          IF (PRESENT(opt_lhas_output))  opt_lhas_output = .TRUE.
        ENDIF
#endif 
      ELSE
        CALL write_name_list(output_file(i), check_open_file(output_file(i)%out_event))
      ENDIF

      ! -------------------------------------------------
      ! Check if files have to be closed
      ! -------------------------------------------------

      IF (check_close_file(output_file(i)%out_event) .AND.  &
        & (output_file(i)%io_proc_id == p_pe)) THEN 
        CALL close_output_file(output_file(i))
        CALL message (routine, 'closed '//TRIM(get_current_filename(output_file(i)%out_event)),all_print=.TRUE.)
      END IF

      ! -------------------------------------------------
      ! add I/O PE of output file to the "output_list"
      ! -------------------------------------------------

      list_idx = -1
      DO j=1,noutput_pe_list
        IF (output_pe_list(j) == output_file(i)%io_proc_id) THEN 
          list_idx = j
          EXIT
        END IF
      END DO
      IF (list_idx == -1) THEN
        noutput_pe_list = noutput_pe_list + 1
        list_idx   = noutput_pe_list
      END IF
      output_pe_list(list_idx) = output_file(i)%io_proc_id

      ! -------------------------------------------------
      ! hand-shake protocol: step finished!
      ! -------------------------------------------------
      CALL pass_output_step(output_file(i)%out_event)

    ENDDO OUTFILE_LOOP

    ! If asynchronous I/O is enabled, the compute PEs can now start
    ! the I/O PEs
#ifndef NOMPI
    IF (use_async_name_list_io  .AND. &
      & .NOT.my_process_is_io() .AND. &
      & .NOT.my_process_is_mpi_test()) THEN
      CALL compute_start_async_io(jstep, output_pe_list, noutput_pe_list)
    END IF
#endif
    ! Handle incoming "output step completed" messages: After all
    ! participating I/O PE's have acknowledged the completion of their
    ! write processes, we trigger a "ready file" on the first I/O PE.
    IF (.NOT.my_process_is_mpi_test()) THEN
       IF ((      use_async_name_list_io .AND. my_process_is_mpi_ioroot()) .OR.  &
            & (.NOT. use_async_name_list_io .AND. my_process_is_mpi_workroot())) THEN
          ev => all_events
          HANDLE_COMPLETE_STEPS : DO
             IF (.NOT. ASSOCIATED(ev)) EXIT HANDLE_COMPLETE_STEPS
             IF (.NOT. is_output_step_complete(ev) .OR.  &
                  & is_output_event_finished(ev)) THEN 
                ev => ev%next
                CYCLE HANDLE_COMPLETE_STEPS
             END IF
             !--- write ready file
             IF (check_write_readyfile(ev%output_event))  CALL write_ready_file(ev)
             ! launch a non-blocking request to all participating PEs to
             ! acknowledge the completion of the next output event
             CALL trigger_output_step_irecv(ev)
          END DO HANDLE_COMPLETE_STEPS
       END IF
    END IF
    IF (ltimer) CALL timer_stop(timer_write_output)
    IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": write_name_list_output done."
  END SUBROUTINE write_name_list_output


  !------------------------------------------------------------------------------------------------
  !> Create a "ready file"
  !
  !  A "ready file" is a technique for handling dependencies between
  !  the NWP processes at DWD: When a program - parallel or
  !  sequential, shell script or binary - produces some output which
  !  is necessary for other running applications, then the completion
  !  of the write process signals this by creating a small file (size:
  !  a few bytes). Only when this file exists, the second program
  !  starts reading its input data. Implicity, this assumes that a
  !  file system creates (and closes) files in the same order as they
  !  are written by the program.
  !
  SUBROUTINE write_ready_file(ev)
    TYPE(t_par_output_event), POINTER :: ev
    ! local variables
    CHARACTER(LEN=*), PARAMETER         :: routine = modname//"::write_ready_file"
    CHARACTER(LEN=FILENAME_MAX)         :: rdy_filename
    CHARACTER(LEN=FILENAME_MAX)         :: forecast_delta_str
    TYPE(datetime),  POINTER            :: mtime_begin, mtime_date
    TYPE(timedelta), POINTER            :: forecast_delta
    INTEGER                             :: iunit
    TYPE (t_keyword_list), POINTER      :: keywords     => NULL()

    CALL setCalendar(PROLEPTIC_GREGORIAN)

    ! compute current forecast time (delta):
    mtime_date     => newDatetime(TRIM(get_current_date(ev)))
    mtime_begin    => newDatetime(TRIM(ev%output_event%event_data%sim_start))
    forecast_delta => newTimedelta("P01D")
    CALL getTimeDeltaFromDateTime(mtime_date, mtime_begin, forecast_delta)
    WRITE (forecast_delta_str,'(4(i2.2))') forecast_delta%day, forecast_delta%hour, &
      &                                    forecast_delta%minute, forecast_delta%second 
    CALL deallocateDatetime(mtime_date)
    CALL deallocateDatetime(mtime_begin)
    CALL deallocateTimedelta(forecast_delta)

    ! substitute tokens in ready file name
    CALL associate_keyword("<path>",            TRIM(model_base_dir),            keywords)
    CALL associate_keyword("<datetime>",        TRIM(get_current_date(ev)),      keywords)
    CALL associate_keyword("<ddhhmmss>",        TRIM(forecast_delta_str),        keywords)
    rdy_filename = TRIM(with_keywords(keywords, ev%output_event%event_data%name))

    IF ((      use_async_name_list_io .AND. my_process_is_mpi_ioroot()) .OR.  &
      & (.NOT. use_async_name_list_io .AND. my_process_is_stdio())) THEN
      WRITE (0,*) 'Write ready file "', TRIM(rdy_filename), '"'
    END IF

    ! actually create ready file:
    iunit = find_next_free_unit(10,20)
    OPEN (iunit, file=TRIM(rdy_filename), form='formatted')
    WRITE(iunit, '(A)') 'ready'
    CLOSE(iunit)
  END SUBROUTINE write_ready_file


  !------------------------------------------------------------------------------------------------
  !> Write an output name list. Called by non-IO PEs.
  !
  SUBROUTINE write_name_list(of, l_first_write)

#ifndef NOMPI
#ifdef  __SUNPRO_F95
  INCLUDE "mpif.h"
#else
    USE mpi, ONLY: MPI_LOCK_EXCLUSIVE, MPI_MODE_NOCHECK
#endif
#endif

    TYPE (t_output_file), INTENT(INOUT), TARGET :: of
    LOGICAL,              INTENT(IN)            :: l_first_write
    ! local variables:
    CHARACTER(LEN=*), PARAMETER                 :: routine = modname//"::write_name_list"
    REAL(wp),         PARAMETER                 :: SYNC_ERROR_PRINT_TOL = 1e-13_wp
    INTEGER,          PARAMETER                 :: iUNKNOWN = 0
    INTEGER,          PARAMETER                 :: iINTEGER = 1
    INTEGER,          PARAMETER                 :: iREAL    = 2
                                               
    INTEGER                                     :: tl, i_dom, i_log_dom, i, iv, jk, n_points, &
      &                                            nlevs, nblks, nindex, mpierr, lonlat_id,   &
      &                                            idata_type
    INTEGER(i8)                                 :: ioff
    TYPE (t_var_metadata), POINTER              :: info
    TYPE(t_reorder_info),  POINTER              :: p_ri
    REAL(wp),          ALLOCATABLE              :: r_ptr(:,:,:)
    INTEGER,           ALLOCATABLE              :: i_ptr(:,:,:)
    REAL(wp),          ALLOCATABLE              :: r_out_recv(:,:)
    REAL(wp),              POINTER              :: r_out_wp(:,:)
    INTEGER,           ALLOCATABLE              :: r_out_int(:,:)
    REAL(sp),          ALLOCATABLE, TARGET      :: r_out_sp(:,:)
    REAL(dp),          ALLOCATABLE, TARGET      :: r_out_dp(:,:)
    TYPE(t_comm_gather_pattern), POINTER        :: p_pat
    LOGICAL                                     :: l_error
    LOGICAL                                     :: have_GRIB

    ! Offset in memory window for async I/O
    ioff = 0_i8

    i_dom = of%phys_patch_id
    i_log_dom = of%log_patch_id

    tl = 0 ! to prevent warning
    
    ! ---------------------------------------------------------
    ! PE#0 : store variable meta-info to be accessed by I/O PEs
    ! ---------------------------------------------------------
#ifndef NOMPI
    ! In case of async IO: Lock own window before writing to it
    IF (      use_async_name_list_io   .AND.  &
      & .NOT. my_process_is_mpi_test() .AND.  &
      &       my_process_is_mpi_workroot()) THEN
      CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, p_pe_work, MPI_MODE_NOCHECK, of%mem_win%mpi_win_metainfo, mpierr)
    END IF
#endif

    DO iv = 1, of%num_vars
      info => of%var_desc(iv)%info
      CALL metainfo_write_to_memwin(of%mem_win, iv, info)
    END DO

#ifndef NOMPI
    ! In case of async IO: Done writing to memory window, unlock it
    IF (      use_async_name_list_io   .AND.  &
      & .NOT. my_process_is_mpi_test() .AND.  &
      &       my_process_is_mpi_workroot()) THEN
      CALL MPI_Win_unlock(p_pe_work, of%mem_win%mpi_win_metainfo, mpierr)
    END IF
#endif

#ifndef NOMPI
    ! In case of async IO: Lock own window before writing to it
    IF(use_async_name_list_io .AND. .NOT.my_process_is_mpi_test()) THEN
      CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, p_pe_work, MPI_MODE_NOCHECK, of%mem_win%mpi_win, mpierr)
    END IF
#endif

    ! ----------------------------------------------------
    ! Go over all name list variables for this output file
    ! ----------------------------------------------------
    DO iv = 1, of%num_vars

      info => of%var_desc(iv)%info
!dbg      write(0,*)'>>>>>>>>>>>>>>>>>>>>>>>>> VARNAME: ',info%name

      ! inspect time-constant variables only if we are writing the
      ! first step in this file:
      IF ((info%isteptype == TSTEP_CONSTANT) .AND. .NOT. l_first_write) CYCLE

      ! Check if first dimension of array is nproma.
      ! Otherwise we got an array which is not suitable for this output scheme.
      IF(info%used_dimensions(1) /= nproma) &
        CALL finish(routine,'1st dim is not nproma: '//TRIM(info%name))

      idata_type = iUNKNOWN

      ! For time level dependent elements: set time level and check if
      ! time level is present:
        ! set a default time level (which is not used anyway, but must
        ! be a valid array subscript):
      tl = 1
#ifndef __NO_ICON_ATMO__
      IF (.NOT. ASSOCIATED(of%var_desc(iv)%r_ptr)  .AND.    &
        & .NOT. ASSOCIATED(of%var_desc(iv)%i_ptr)) THEN
        SELECT CASE (info%tlev_source)
          CASE(0); tl = nnow(i_log_dom)
          CASE(1); tl = nnow_rcf(i_log_dom)
          CASE(2); tl = nnew(i_log_dom)
          CASE(3); tl = nnew_rcf(i_log_dom)
          CASE DEFAULT
            CALL finish(routine,'Unsupported tlev_source')
        END SELECT
        IF(tl<=0 .OR. tl>max_time_levels) &
          CALL finish(routine, 'Illegal time level in nnow()/nnow_rcf()')
        ! Check if present
        IF (.NOT. ASSOCIATED(of%var_desc(iv)%tlev_rptr(tl)%p)   .AND.   &
          & .NOT. ASSOCIATED(of%var_desc(iv)%tlev_iptr(tl)%p)) THEN
          CALL finish(routine,'Actual timelevel not in '//TRIM(info%name))
        END IF
      ENDIF
#endif

      IF (info%lcontained) THEN
        nindex = info%ncontained
      ELSE
        nindex = 1
      ENDIF

      ! determine, if this is a REAL or an INTEGER variable:
      IF (ASSOCIATED(of%var_desc(iv)%r_ptr) .OR.  &
        & ASSOCIATED(of%var_desc(iv)%tlev_rptr(tl)%p)) THEN
        idata_type = iREAL
      ELSE IF (ASSOCIATED(of%var_desc(iv)%i_ptr) .OR.  &
        & ASSOCIATED(of%var_desc(iv)%tlev_iptr(tl)%p)) THEN
        idata_type = iINTEGER
      END IF

      SELECT CASE (info%ndims)
      CASE (1)
        CALL message(routine, info%name)
        CALL finish(routine,'1d arrays not handled yet.')
      CASE (2)
        ! 2D fields: Make a 3D copy of the array
        IF (idata_type == iREAL)    ALLOCATE(r_ptr(info%used_dimensions(1),1,info%used_dimensions(2)))
        IF (idata_type == iINTEGER) ALLOCATE(i_ptr(info%used_dimensions(1),1,info%used_dimensions(2)))

        IF (ASSOCIATED(of%var_desc(iv)%r_ptr)) THEN
          r_ptr(:,1,:) = of%var_desc(iv)%r_ptr(:,:,nindex,1,1)
        ELSE IF (ASSOCIATED(of%var_desc(iv)%i_ptr)) THEN
          i_ptr(:,1,:) = of%var_desc(iv)%i_ptr(:,:,nindex,1,1)
        ELSE IF (ASSOCIATED(of%var_desc(iv)%tlev_rptr(tl)%p)) THEN
          r_ptr(:,1,:) = of%var_desc(iv)%tlev_rptr(tl)%p(:,:,nindex,1,1)
        ELSE IF (ASSOCIATED(of%var_desc(iv)%tlev_iptr(tl)%p)) THEN
          i_ptr(:,1,:) = of%var_desc(iv)%tlev_iptr(tl)%p(:,:,nindex,1,1)
        ELSE
          CALL finish(routine, "Internal error!")
        ENDIF
      CASE (3)
        ! 3D fields: Here we could just set a pointer to the
        ! array... if there were no post-ops
        IF (idata_type == iREAL)    ALLOCATE(r_ptr(info%used_dimensions(1), &
          &                                        info%used_dimensions(2), &
          &                                        info%used_dimensions(3)))
        IF (idata_type == iINTEGER) ALLOCATE(i_ptr(info%used_dimensions(1), &
          &                                        info%used_dimensions(2), &
          &                                        info%used_dimensions(3)))

        IF(ASSOCIATED(of%var_desc(iv)%r_ptr)) THEN
          r_ptr = of%var_desc(iv)%r_ptr(:,:,:,nindex,1)
        ELSE IF (ASSOCIATED(of%var_desc(iv)%i_ptr)) THEN
          i_ptr = of%var_desc(iv)%i_ptr(:,:,:,nindex,1)
        ELSE IF (ASSOCIATED(of%var_desc(iv)%tlev_rptr(tl)%p)) THEN
          r_ptr = of%var_desc(iv)%tlev_rptr(tl)%p(:,:,:,nindex,1)
        ELSE IF (ASSOCIATED(of%var_desc(iv)%tlev_iptr(tl)%p)) THEN
          i_ptr = of%var_desc(iv)%tlev_iptr(tl)%p(:,:,:,nindex,1)
        ELSE
          CALL finish(routine, "Internal error!")
        ENDIF
      CASE (4)
        CALL message(routine, info%name)
        CALL finish(routine,'4d arrays not handled yet.')
      CASE (5)
        CALL message(routine, info%name)
        CALL finish(routine,'5d arrays not handled yet.')
      CASE DEFAULT
        CALL message(routine, info%name)
        CALL finish(routine,'dimension not set.')
      END SELECT

      ! --------------------------------------------------------
      ! Perform post-ops (small arithmetic operations on fields)
      ! --------------------------------------------------------

      IF (of%var_desc(iv)%info%post_op%ipost_op_type == POST_OP_SCALE) THEN
        IF (idata_type == iINTEGER) CALL finish(routine, "Not yet implemented!")
        CALL perform_post_op(of%var_desc(iv)%info%post_op, r_ptr)
      END IF

      IF(info%ndims == 2) THEN
        nlevs = 1
      ELSE
        nlevs = info%used_dimensions(2)
      ENDIF

      ! Get pointer to appropriate reorder_info

      SELECT CASE (info%hgrid)
      CASE (GRID_UNSTRUCTURED_CELL)
        p_ri => patch_info(i_dom)%cells
        p_pat => patch_info(i_dom)%p_pat_c
      CASE (GRID_UNSTRUCTURED_EDGE)
        p_ri => patch_info(i_dom)%edges
        p_pat => patch_info(i_dom)%p_pat_e
      CASE (GRID_UNSTRUCTURED_VERT)
        p_ri => patch_info(i_dom)%verts
        p_pat => patch_info(i_dom)%p_pat_v
#ifndef __NO_ICON_ATMO__
      CASE (GRID_REGULAR_LONLAT)
        lonlat_id = info%hor_interp%lonlat_id
        p_ri  => lonlat_info(lonlat_id, i_log_dom)%ri
        p_pat => lonlat_grid_list(lonlat_id)%p_pat(i_log_dom)
#endif
      CASE default
        CALL finish(routine,'unknown grid type')
      END SELECT

      IF (.NOT.use_async_name_list_io .OR. my_process_is_mpi_test()) THEN

        ! -------------------
        ! No asynchronous I/O
        ! -------------------
        !
        ! gather the array on stdio PE and write it out there

        n_points = p_ri%n_glb
        nblks = (n_points-1)/nproma + 1

        IF (idata_type == iREAL) THEN

          ALLOCATE(r_out_dp(MERGE(n_points, 0, &
            &                      my_process_is_mpi_workroot()), nlevs))
          r_out_wp => r_out_dp
          r_out_wp(:,:) = 0._wp
          CALL exchange_data(r_ptr(:,:,:), r_out_wp(:,:), p_pat)

        ELSE IF (idata_type == iINTEGER) THEN

          ALLOCATE(r_out_int(MERGE(n_points, 0, &
            &                      my_process_is_mpi_workroot()), nlevs))
          r_out_int(:,:) = 0
          CALL exchange_data(i_ptr(:,:,:), r_out_int(:,:), p_pat)

        END IF

        have_GRIB = of%output_type == FILETYPE_GRB .OR. of%output_type == FILETYPE_GRB2
        IF(my_process_is_mpi_workroot()) THEN

          IF (use_dp_mpi2io .or. have_GRIB) THEN
            IF (.NOT. ALLOCATED(r_out_dp)) &
              ALLOCATE(r_out_dp(n_points, nlevs))
            IF (idata_type == iINTEGER) &
              r_out_dp(:,:) = REAL(r_out_int(:,:), dp)
          ELSE
            ALLOCATE(r_out_sp(n_points, nlevs))
            IF (idata_type == iREAL) THEN
              r_out_sp(:,:) = REAL(r_out_wp(:,:), sp)
            ELSE IF (idata_type == iINTEGER) THEN
              r_out_sp(:,:) = REAL(r_out_int(:,:), sp)
            END IF
          END IF

          ! ------------------
          ! case of a test run
          ! ------------------
          !
          ! compare results on worker PEs and test PE
          IF (p_test_run  .AND.  use_dp_mpi2io) THEN
            ! Currently we don't do the check for REAL*4, we would need
            ! p_send/p_recv for this type
            IF(.NOT. my_process_is_mpi_test()) THEN
              ! Send to test PE
              CALL p_send(r_out_dp, process_mpi_all_test_id, 1)
            ELSE
              ! Receive result from parallel worker PEs
              ALLOCATE(r_out_recv(n_points,nlevs))
              CALL p_recv(r_out_recv, process_mpi_all_workroot_id, 1)
              ! check for correctness
              l_error = .FALSE.
              DO jk = 1, nlevs
                DO i = 1, n_points
                  IF (r_out_recv(i,jk) /= r_out_dp(i,jk)) THEN
                    ! do detailed print-out only for "large" errors:
                    IF (ABS(r_out_recv(i,jk) - r_out_dp(i,jk)) > SYNC_ERROR_PRINT_TOL) THEN
                      WRITE (0,*) 'Sync error test PE/worker PEs for ', TRIM(info%name)
                      WRITE (0,*) "global pos (", idx_no(i), ",", blk_no(i),")"
                      WRITE (0,*) "level", jk, "/", nlevs
                      WRITE (0,*) "vals: ", r_out_recv(i,jk), r_out_dp(i,jk)
                      l_error = .TRUE.
                    END IF
                  END IF
                ENDDO
              END DO
              if (l_error)   CALL finish(routine,"Sync error!")
              DEALLOCATE(r_out_recv)
            ENDIF
          ENDIF

        ENDIF

        ! set some GRIB2 keys that may have changed during simulation
        IF  (( of%output_type == FILETYPE_GRB2 ) .AND.  &
          &  ( my_process_is_mpi_workroot() )    .AND. &
          &  ( .NOT. my_process_is_mpi_test() )) THEN
          CALL set_timedependent_GRIB2_keys(of%cdiVlistID, info%cdiVarID, info,      &
            &                               TRIM(of%out_event%output_event%event_data%sim_start), &
            &                               TRIM(get_current_date(of%out_event)))
        END IF

        ! ----------
        ! write data
        ! ----------

        IF (my_process_is_mpi_workroot() .AND. .NOT. my_process_is_mpi_test()) THEN
          IF (use_dp_mpi2io .or. have_GRIB) THEN
            CALL streamWriteVar(of%cdiFileID, info%cdiVarID, r_out_dp, 0)
          ELSE
            CALL streamWriteVarF(of%cdiFileID, info%cdiVarID, r_out_sp, 0)
          ENDIF
          IF (lkeep_in_sync) CALL streamSync(of%cdiFileID)
        ENDIF

        ! clean up
        IF (ALLOCATED(r_out_int)) DEALLOCATE(r_out_int)
        IF (ALLOCATED(r_out_sp)) DEALLOCATE(r_out_sp)
        IF (ALLOCATED(r_out_dp)) DEALLOCATE(r_out_dp)

      ELSE

        ! ------------------------
        ! Asynchronous I/O is used
        ! ------------------------

        ! just copy the OWN DATA points to the memory window
        DO jk = 1, nlevs
          IF (use_dp_mpi2io) THEN
            IF (idata_type == iREAL) THEN
              DO i = 1, p_ri%n_own
                of%mem_win%mem_ptr_dp(ioff+INT(i,i8)) = REAL(r_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i)),dp)
              ENDDO
            END IF
            IF (idata_type == iINTEGER) THEN
              DO i = 1, p_ri%n_own
                of%mem_win%mem_ptr_dp(ioff+INT(i,i8)) = REAL(i_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i)),dp)
              ENDDO
            END IF
          ELSE
            IF (idata_type == iREAL) THEN
              DO i = 1, p_ri%n_own
                of%mem_win%mem_ptr_sp(ioff+INT(i,i8)) = REAL(r_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i)),sp)
              ENDDO
            END IF
            IF (idata_type == iINTEGER) THEN
              DO i = 1, p_ri%n_own
                of%mem_win%mem_ptr_sp(ioff+INT(i,i8)) = REAL(i_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i)),sp)
              ENDDO
            END IF
          END IF
          ioff = ioff + INT(p_ri%n_own,i8)
        END DO ! nlevs

      END IF

      ! clean up
      IF (ALLOCATED(r_ptr)) DEALLOCATE(r_ptr)
      IF (ALLOCATED(i_ptr)) DEALLOCATE(i_ptr)

    ENDDO

#ifndef NOMPI
    ! In case of async IO: Done writing to memory window, unlock it
    IF(use_async_name_list_io .AND. .NOT.my_process_is_mpi_test()) THEN
      CALL MPI_Win_unlock(p_pe_work, of%mem_win%mpi_win, mpierr)
    END IF
#endif

  END SUBROUTINE write_name_list


  !------------------------------------------------------------------------------------------------
  !> Returns if it is time for the next output step
  !  Please note:
  !  This function returns .TRUE. whenever the next output time of any name list
  !  is reached at the simulation step @p jstep.
  !
  FUNCTION istime4name_list_output(jstep) 
    LOGICAL :: istime4name_list_output
    INTEGER, INTENT(IN)   :: jstep            ! simulation time step
    ! local variables
    INTEGER :: i
    LOGICAL :: ret, ret_local

    ret = .FALSE.
    IF (ALLOCATED(output_file)) THEN
       ! note: there may be cases where no output namelist has been
       ! defined. thus we must check if "output_file" has been
       ! allocated.
       DO i = 1, SIZE(output_file)
          ret_local = is_output_step(output_file(i)%out_event, jstep)
          ret = ret .OR. ret_local
       END DO
    END IF
    istime4name_list_output = ret
  END FUNCTION istime4name_list_output


  !------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------
  ! The following routines are only needed for asynchronous IO
  !------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------

#ifdef NOMPI
  ! Just define the entry point of name_list_io_main_proc, it will never be called

  SUBROUTINE name_list_io_main_proc(sim_step_info, isample)
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info
    INTEGER,                INTENT(in) :: isample
  END SUBROUTINE name_list_io_main_proc

#else


  !-------------------------------------------------------------------------------------------------
  !> Main routine for I/O PEs.
  !  Please note that this routine never returns.
  !
  SUBROUTINE name_list_io_main_proc(sim_step_info, isample)
    !> Data structure containing all necessary data for mapping an
    !  output time stamp onto a corresponding simulation step index.
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info
    INTEGER,                INTENT(in) :: isample
    ! local variables:

#ifndef __NO_ICON_ATMO__
    LOGICAL             :: done, l_complete, lhas_output, &
      &                    lset_timers_for_idle_pe
    INTEGER             :: jg, jstep, i
    TYPE(t_par_output_event), POINTER :: ev

    ! setup of meteogram output
    DO jg =1,n_dom
      IF (meteogram_output_config(jg)%lenabled) THEN
        CALL meteogram_init(meteogram_output_config(jg), jg)
      END IF
    END DO

    ! Initialize name list output, this is a collective call for all PEs
    CALL init_name_list_output(sim_step_info)

    ! Tell the compute PEs that we are ready to work
    IF (ANY(output_file(:)%io_proc_id == p_pe)) THEN 
      CALL async_io_send_handshake(0)
    END IF

    ! Enter I/O loop
    DO
      ! skip loop, if this output PE is idle:
      IF (ALL(output_file(:)%io_proc_id /= p_pe)) EXIT

      ! Wait for a message from the compute PEs to start
      CALL async_io_wait_for_start(done, jstep)
      IF(done) EXIT ! leave loop, we are done

      ! perform I/O
      CALL write_name_list_output(jstep, opt_lhas_output=lhas_output)

      ! Inform compute PEs that we are done, if this I/O PE has
      ! written output:
      IF (lhas_output)  CALL async_io_send_handshake(jstep)

      ! Handle final pending "output step completed" messages: After
      ! all participating I/O PE's have acknowledged the completion of
      ! their write processes, we trigger a "ready file" on the first
      ! I/O PE.
      IF (.NOT.my_process_is_mpi_test()  .AND.  &
        & use_async_name_list_io .AND. my_process_is_mpi_ioroot()) THEN

        ! Go over all output files
        l_complete = .TRUE.
        OUTFILE_LOOP : DO i=1,SIZE(output_file)
          l_complete = l_complete .AND. is_output_event_finished(output_file(i)%out_event)
        END DO OUTFILE_LOOP

        IF (l_complete) THEN
          IF (ldebug)   WRITE (0,*) p_pe, ": wait for fellow I/O PEs..."
          WAIT_FINAL : DO           
            CALL blocking_wait_for_irecvs(all_events)
            ev => all_events
            l_complete = .TRUE.
            HANDLE_COMPLETE_STEPS : DO
              IF (.NOT. ASSOCIATED(ev)) EXIT HANDLE_COMPLETE_STEPS
              
              !--- write ready file
              IF (check_write_readyfile(ev%output_event))  CALL write_ready_file(ev)
              IF (.NOT. is_output_event_finished(ev)) THEN 
                l_complete = .FALSE.
                IF (is_output_step_complete(ev)) THEN 
                  CALL trigger_output_step_irecv(ev)
                END IF
              END IF
              ev => ev%next
            END DO HANDLE_COMPLETE_STEPS
            IF (l_complete) EXIT WAIT_FINAL
          END DO WAIT_FINAL
          IF (ldebug)  WRITE (0,*) p_pe, ": Finalization sequence"
        END IF
      END IF

    ENDDO

    ! Finalization sequence:
    CALL close_name_list_output

    ! finalize meteogram output
    DO jg = 1, n_dom
      IF (meteogram_output_config(jg)%lenabled)  CALL meteogram_finalize(jg)
    END DO

    DO jg = 1, max_dom
      DEALLOCATE(meteogram_output_config(jg)%station_list)
    END DO


    ! Purely idle output PEs: Empty calls of timer start/stop. For
    ! this pathological case it is important to call the same timers
    ! as the "normal" output PEs. Otherwise we will get a deadlock
    ! situation when computing the global sums for these timers.
    lset_timers_for_idle_pe = .FALSE.
    IF (ltimer .AND. ALLOCATED(output_file)) THEN
      IF (ALL(output_file(:)%io_proc_id /= p_pe))  lset_timers_for_idle_pe = .TRUE.
    END IF
    IF (ltimer .AND. .NOT. ALLOCATED(output_file)) lset_timers_for_idle_pe = .TRUE.
    IF (lset_timers_for_idle_pe) THEN
        CALL timer_start(timer_write_output)
        CALL timer_stop(timer_write_output)
    END IF

    IF (ltimer) CALL print_timer
    ! Shut down MPI
    CALL stop_mpi
    STOP
#endif

  END SUBROUTINE name_list_io_main_proc
  !------------------------------------------------------------------------------------------------


  !------------------------------------------------------------------------------------------------
  !> Output routine on the IO PEs
  !
  !  @note This subroutine is called by asynchronous I/O PEs only.
  !
  SUBROUTINE io_proc_write_name_list(of, l_first_write)

#ifdef __SUNPRO_F95
    INCLUDE "mpif.h"
#else
    USE mpi, ONLY: MPI_ADDRESS_KIND, MPI_LOCK_SHARED, MPI_MODE_NOCHECK
#endif

    TYPE (t_output_file), TARGET, INTENT(IN) :: of
    LOGICAL                     , INTENT(IN) :: l_first_write

    CHARACTER(LEN=*), PARAMETER    :: routine = modname//"::io_proc_write_name_list"

    INTEGER                        :: nval, nlev_max, iv, jk, nlevs, mpierr, nv_off, np, i_dom, &
      &                               lonlat_id, i_log_dom, ierrstat,                           &
      &                               dst_start, dst_end, src_start, src_end
    INTEGER(KIND=MPI_ADDRESS_KIND) :: ioff(0:num_work_procs-1)
    INTEGER                        :: voff(0:num_work_procs-1), nv_off_np(0:num_work_procs)
    REAL(sp), ALLOCATABLE          :: var1_sp(:), var3_sp(:)
    REAL(dp), ALLOCATABLE          :: var1_dp(:), var3_dp(:)
    TYPE (t_var_metadata), POINTER :: info
    TYPE (t_var_metadata)          :: updated_info
    TYPE(t_reorder_info) , POINTER :: p_ri
    LOGICAL                        :: have_GRIB
    INTEGER, ALLOCATABLE           :: bufr_metainfo(:)

    !-- for timing
    CHARACTER(len=10)              :: ctime
    REAL(dp)                       :: t_get, t_write, t_copy, t_0, mb_get, mb_wr

  !------------------------------------------------------------------------------------------------
#if defined (__SX__) && !defined (NOMPI)
! It may be necessary that var1 is in global memory on NEC
! (Note: this is only allowed when we compile with MPI.)
!CDIR GM_ARRAY(var1)
#endif

    CALL date_and_time(TIME=ctime)
    WRITE (0, '(a,i0,a)') '#################### I/O PE ',p_pe,' starting I/O at '//ctime

    t_get   = 0.d0
    t_write = 0.d0
    t_copy  = 0.d0
    mb_get  = 0.d0
    mb_wr   = 0.d0


    ! Get maximum number of data points in a slice and allocate tmp variables

    i_dom = of%phys_patch_id
    nval = MAX(patch_info(i_dom)%cells%n_glb, &
               patch_info(i_dom)%edges%n_glb, &
               patch_info(i_dom)%verts%n_glb)
#ifndef __NO_ICON_ATMO__
    ! take also the lon-lat grids into account
    DO iv = 1, of%num_vars
      info => of%var_desc(iv)%info
      IF (info%hgrid == GRID_REGULAR_LONLAT) THEN
        lonlat_id = info%hor_interp%lonlat_id
        i_log_dom = of%log_patch_id
        p_ri  => lonlat_info(lonlat_id, i_log_dom)%ri
        nval = MAX(nval, p_ri%n_glb)
      END IF
    END DO
#endif

    nlev_max = 1
    DO iv = 1, of%num_vars
      info => of%var_desc(iv)%info
      IF(info%ndims == 3) nlev_max = MAX(nlev_max, info%used_dimensions(2))
    ENDDO

    IF (use_dp_mpi2io) THEN
      ALLOCATE(var1_dp(nval*nlev_max), STAT=ierrstat)
    ELSE
      ALLOCATE(var1_sp(nval*nlev_max), STAT=ierrstat)
    ENDIF
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! retrieve info object from PE#0 (via a separate MPI memory
    ! window)
    ALLOCATE(bufr_metainfo(of%num_vars*metainfo_get_size()), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    CALL MPI_Win_lock(MPI_LOCK_SHARED, 0, MPI_MODE_NOCHECK, of%mem_win%mpi_win_metainfo, mpierr)
    CALL MPI_Get(bufr_metainfo, SIZE(bufr_metainfo), p_int, 0, &
      &          0_MPI_ADDRESS_KIND, SIZE(bufr_metainfo), p_int, of%mem_win%mpi_win_metainfo, mpierr)
    CALL MPI_Win_unlock(0, of%mem_win%mpi_win_metainfo, mpierr)

    ! Go over all name list variables for this output file

    ioff(:) = 0_MPI_ADDRESS_KIND
    DO iv = 1, of%num_vars
      ! POINTER to this variable's meta-info
      info => of%var_desc(iv)%info
      ! get also an update for this variable's meta-info (separate object)
      CALL metainfo_get_from_memwin(bufr_metainfo, iv, updated_info)

      IF ( of%output_type == FILETYPE_GRB2 ) THEN
        CALL set_timedependent_GRIB2_keys(of%cdiVlistID, info%cdiVarID, updated_info, &
          &                               TRIM(of%out_event%output_event%event_data%sim_start),    &
          &                               TRIM(get_current_date(of%out_event)))
      END IF

      ! inspect time-constant variables only if we are writing the
      ! first step in this file:
      IF ((info%isteptype == TSTEP_CONSTANT) .AND. .NOT. l_first_write) CYCLE

      IF(info%ndims == 2) THEN
        nlevs = 1
      ELSE
        nlevs = info%used_dimensions(2)
      ENDIF
      ! Get pointer to appropriate reorder_info

      SELECT CASE (info%hgrid)
        CASE (GRID_UNSTRUCTURED_CELL)
          p_ri => patch_info(of%phys_patch_id)%cells
        CASE (GRID_UNSTRUCTURED_EDGE)
          p_ri => patch_info(of%phys_patch_id)%edges
        CASE (GRID_UNSTRUCTURED_VERT)
          p_ri => patch_info(of%phys_patch_id)%verts

#ifndef __NO_ICON_ATMO__
        CASE (GRID_REGULAR_LONLAT)
          lonlat_id = info%hor_interp%lonlat_id
          i_log_dom = of%log_patch_id
          p_ri  => lonlat_info(lonlat_id, i_log_dom)%ri
#endif

        CASE DEFAULT
          CALL finish(routine,'unknown grid type')
      END SELECT

      ! Retrieve part of variable from every worker PE using MPI_Get

      nv_off  = 0
      DO np = 0, num_work_procs-1

        IF(p_ri%pe_own(np) == 0) CYCLE

        nval = p_ri%pe_own(np)*nlevs ! Number of words to transfer

        t_0 = p_mpi_wtime()
        CALL MPI_Win_lock(MPI_LOCK_SHARED, np, MPI_MODE_NOCHECK, of%mem_win%mpi_win, mpierr)

        IF (use_dp_mpi2io) THEN
          CALL MPI_Get(var1_dp(nv_off+1), nval, p_real_dp, np, ioff(np), &
            &          nval, p_real_dp, of%mem_win%mpi_win, mpierr)
        ELSE
          CALL MPI_Get(var1_sp(nv_off+1), nval, p_real_sp, np, ioff(np), &
            &          nval, p_real_sp, of%mem_win%mpi_win, mpierr)
        ENDIF

        CALL MPI_Win_unlock(np, of%mem_win%mpi_win, mpierr)
        t_get  = t_get  + p_mpi_wtime() - t_0
        mb_get = mb_get + nval

        ! Update the offset in var1
        nv_off = nv_off + nval

        ! Update the offset in the memory window on compute PEs
        ioff(np) = ioff(np) + INT(nval,i8)

      ENDDO

      ! compute the total offset for each PE
      nv_off       = 0
      nv_off_np(0) = 0
      DO np = 0, num_work_procs-1
        voff(np)        = nv_off
        nval            = p_ri%pe_own(np)*nlevs
        nv_off          = nv_off + nval
        nv_off_np(np+1) = nv_off_np(np) + p_ri%pe_own(np)
      END DO

      ! var1 is stored in the order in which the variable was stored on compute PEs,
      ! get it back into the global storage order

      t_0 = p_mpi_wtime() ! performance measurement
      have_GRIB = of%output_type == FILETYPE_GRB .OR. of%output_type == FILETYPE_GRB2
      IF (use_dp_mpi2io .OR. have_GRIB) THEN
        ALLOCATE(var3_dp(p_ri%n_glb), STAT=ierrstat) ! Must be allocated to exact size
      ELSE
        ALLOCATE(var3_sp(p_ri%n_glb), STAT=ierrstat) ! Must be allocated to exact size
      ENDIF
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      t_copy = t_copy + p_mpi_wtime() - t_0 ! performance measurement

      LevelLoop: DO jk = 1, nlevs

        t_0 = p_mpi_wtime() ! performance measurement

        IF (use_dp_mpi2io) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(dst_start, dst_end, src_start, src_end)
          DO np = 0, num_work_procs-1
            dst_start = nv_off_np(np)+1
            dst_end   = nv_off_np(np+1)
            src_start = voff(np)+1
            src_end   = voff(np)+p_ri%pe_own(np)
            voff(np)  = src_end
            var3_dp(p_ri%reorder_index(dst_start:dst_end)) = var1_dp(src_start:src_end)
          ENDDO
!$OMP END DO
!$OMP END PARALLEL
        ELSE

          IF (have_GRIB) THEN
            ! ECMWF GRIB-API/CDI has only a double precision interface at the date of coding this
!$OMP PARALLEL
!$OMP DO PRIVATE(dst_start, dst_end, src_start, src_end)
            DO np = 0, num_work_procs-1
              dst_start = nv_off_np(np)+1
              dst_end   = nv_off_np(np+1)
              src_start = voff(np)+1
              src_end   = voff(np)+p_ri%pe_own(np)
              voff(np)  = src_end
              var3_dp(p_ri%reorder_index(dst_start:dst_end)) = REAL(var1_sp(src_start:src_end), dp)
            ENDDO
!$OMP END DO
!$OMP END PARALLEL
          ELSE
!$OMP PARALLEL
!$OMP DO PRIVATE(dst_start, dst_end, src_start, src_end)
            DO np = 0, num_work_procs-1
              dst_start = nv_off_np(np)+1
              dst_end   = nv_off_np(np+1)
              src_start = voff(np)+1
              src_end   = voff(np)+p_ri%pe_own(np)
              voff(np)  = src_end
              var3_sp(p_ri%reorder_index(dst_start:dst_end)) = var1_sp(src_start:src_end)
            ENDDO
!$OMP END DO
!$OMP END PARALLEL
          END IF
        ENDIF
        t_copy = t_copy + p_mpi_wtime() - t_0 ! performance measurement
        ! Write calls (via CDIs) of the asynchronous I/O PEs:
        t_0 = p_mpi_wtime() ! performance measurement
        IF (use_dp_mpi2io .OR. have_GRIB) THEN
          CALL streamWriteVarSlice(of%cdiFileID, info%cdiVarID, jk-1, var3_dp, 0)
          mb_wr = mb_wr + REAL(SIZE(var3_dp), wp)
        ELSE
          CALL streamWriteVarSliceF(of%cdiFileID, info%cdiVarID, jk-1, var3_sp, 0)
          mb_wr = mb_wr + REAL(SIZE(var3_sp),wp)
        ENDIF
        t_write = t_write + p_mpi_wtime() - t_0 ! performance measurement
      ENDDO LevelLoop

      IF (use_dp_mpi2io .OR. have_GRIB) THEN
        DEALLOCATE(var3_dp, STAT=ierrstat)
      ELSE
        DEALLOCATE(var3_sp, STAT=ierrstat)
      ENDIF
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    ENDDO ! Loop over output variables

    IF (use_dp_mpi2io) THEN
      DEALLOCATE(var1_dp, STAT=ierrstat)
    ELSE
      DEALLOCATE(var1_sp, STAT=ierrstat)
    ENDIF
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    !
    !-- timing report
    !
    CALL date_and_time(TIME=ctime)
    WRITE (0, '(a,i0,a)') '#################### I/O PE ',p_pe,' done at '//ctime
    ! Convert mb_get/mb_wr to MB
    IF (use_dp_mpi2io) THEN
      mb_get = mb_get*8*1.d-6
    ELSE
      mb_get = mb_get*4*1.d-6
    ENDIF
    mb_wr  = mb_wr*4*1.d-6 ! 4 byte since dp output is implicitly converted to sp
    ! writing this message causes a runtime error on the NEC because formatted output to stdio/stderr is limited to 132 chars
#ifndef __SX__
    IF (msg_level >= 12) THEN
      WRITE (0,'(10(a,f10.3))') &  ! remark: CALL message does not work here because it writes only on PE0
           & ' Got ',mb_get,' MB, time get: ',t_get,' s [',mb_get/MAX(1.e-6_wp,t_get), &
           & ' MB/s], time write: ',t_write,' s [',mb_wr/MAX(1.e-6_wp,t_write),        &
           & ' MB/s], times copy: ',t_copy,' s'
   !   CALL message('',message_text)
    ENDIF
#endif

    DEALLOCATE(bufr_metainfo, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  END SUBROUTINE io_proc_write_name_list


  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  ! Flow control routines between compute and IO procs ...
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------
  ! ... called on IO procs:

  !-------------------------------------------------------------------------------------------------
  !> Send a message to the compute PEs that the I/O is ready. The
  !  counterpart on the compute side is compute_wait_for_async_io
  !
  SUBROUTINE async_io_send_handshake(jstep)
    INTEGER, INTENT(IN) :: jstep
    ! local variables
    REAL(wp) :: msg
    TYPE(t_par_output_event), POINTER :: ev

    IF (ldebug) &
         & WRITE (0,*) "pe ", p_pe, ": async_io_send_handshake, jstep=", jstep

    ! --- Send a message from this I/O PE to the compute PE #0
    !
    ! Note: We have to do this in a non-blocking fashion in order to
    !       receive "ready file" messages.
    CALL p_wait() 
    msg = REAL(msg_io_done, wp)
    CALL p_isend(msg, p_work_pe0, 0)

    ! --- I/O PE #0  :  take care of ready files    
    IF(p_pe_work == 0) THEN
      DO 
        IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": trigger, async_io_send_handshake"
        ev => all_events
        HANDLE_COMPLETE_STEPS : DO
          IF (.NOT. ASSOCIATED(ev)) EXIT HANDLE_COMPLETE_STEPS
          IF (.NOT. is_output_step_complete(ev) .OR.  &
            & is_output_event_finished(ev)) THEN 
            ev => ev%next
            CYCLE HANDLE_COMPLETE_STEPS
          END IF
          !--- write ready file
          IF (check_write_readyfile(ev%output_event))  CALL write_ready_file(ev)
          ! launch a non-blocking request to all participating PEs to
          ! acknowledge the completion of the next output event
          CALL trigger_output_step_irecv(ev)
        END DO HANDLE_COMPLETE_STEPS
        IF (p_test()) EXIT
      END DO
    END IF
    CALL p_wait() 

  END SUBROUTINE async_io_send_handshake


  !-------------------------------------------------------------------------------------------------
  !> async_io_wait_for_start: Wait for a message from work PEs that we
  !  should start I/O or finish.  The counterpart on the compute side is
  !  compute_start_async_io/compute_shutdown_async_io
  !
  SUBROUTINE async_io_wait_for_start(done, jstep)
    LOGICAL, INTENT(OUT)          :: done ! flag if we should shut down
    INTEGER, INTENT(OUT)          :: jstep
    ! local variables
    REAL(wp) :: msg(2)
    TYPE(t_par_output_event), POINTER :: ev

    ! Set output parameters to default values
    done  = .FALSE.
    jstep = -1

    ! Receive message that we may start I/O (or should finish)
    ! 
    ! If this I/O PE will write output in this step, or if it has
    ! finished all its tasks and waits for the shutdown message,
    ! launch a non-blocking receive request to compute PE #0:
    !
    ! Note: We have to do this in a non-blocking fashion in order to
    !       receive "ready file" messages.
    !
    CALL p_wait()
    CALL p_irecv(msg, p_work_pe0, 0)
    
    IF(p_pe_work == 0) THEN
      DO 
        ev => all_events
        HANDLE_COMPLETE_STEPS : DO
          IF (.NOT. ASSOCIATED(ev)) EXIT HANDLE_COMPLETE_STEPS
          IF (.NOT. is_output_step_complete(ev) .OR.  &
            & is_output_event_finished(ev)) THEN 
            ev => ev%next
            CYCLE HANDLE_COMPLETE_STEPS
          END IF
          !--- write ready file
          IF (check_write_readyfile(ev%output_event))  CALL write_ready_file(ev)
          ! launch a non-blocking request to all participating PEs to
          ! acknowledge the completion of the next output event
          CALL trigger_output_step_irecv(ev)
        END DO HANDLE_COMPLETE_STEPS

        IF (p_test()) EXIT
      END DO
    END IF

    CALL p_wait()
    
    SELECT CASE(INT(msg(1)))
    CASE(msg_io_start)
      jstep = INT(msg(2))
    CASE(msg_io_shutdown)
      done = .TRUE.
    CASE DEFAULT
      ! Anything else is an error
      CALL finish(modname, 'I/O PE: Got illegal I/O tag')
    END SELECT
  END SUBROUTINE async_io_wait_for_start


  !-------------------------------------------------------------------------------------------------
  ! ... called on compute procs:

  !-------------------------------------------------------------------------------------------------
  !> compute_wait_for_async_io: Wait for a message that the I/O is ready
  !  The counterpart on the I/O side is io_send_handshake
  !
  SUBROUTINE compute_wait_for_async_io(jstep)
    INTEGER, INTENT(IN) :: jstep         !< model step
    ! local variables
    REAL(wp) :: msg
    INTEGER  :: i,j, nwait_list, wait_idx
    INTEGER  :: wait_list(num_io_procs)

    ! Compute PE #0 receives message from I/O PEs
    !
    ! Note: We only need to wait for those I/O PEs which are involved
    !       in the current step.
    IF (p_pe_work==0) THEN
      IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": compute_wait_for_async_io, jstep=",jstep
      IF (jstep == WAIT_UNTIL_FINISHED) THEN
        CALL p_wait()
      ELSE
        wait_list(:) = -1
        nwait_list   =  0

        ! Go over all output files, collect IO PEs
        OUTFILE_LOOP : DO i=1,SIZE(output_file)
          ! Skip this output file if it is not due for output!
          IF (.NOT. is_output_step(output_file(i)%out_event, jstep))  CYCLE OUTFILE_LOOP
          wait_idx = -1
          DO j=1,nwait_list
            IF (wait_list(j) == output_file(i)%io_proc_id) THEN 
              wait_idx = j
              EXIT
            END IF
          END DO
          IF (wait_idx == -1) THEN
            nwait_list = nwait_list + 1
            wait_idx   = nwait_list
          END IF
          wait_list(wait_idx) = output_file(i)%io_proc_id
        END DO OUTFILE_LOOP
        DO i=1,nwait_list
          ! Blocking receive call:
          IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": wait for PE ",  wait_list(i)
          CALL p_recv(msg, wait_list(i), 0)
          ! Just for safety: Check if we got the correct tag
          IF(INT(msg) /= msg_io_done) CALL finish(modname, 'Compute PE: Got illegal I/O tag')
        END DO
      END IF
    ENDIF
    IF (jstep /= WAIT_UNTIL_FINISHED) THEN
      ! Wait in barrier until message is here
      IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": waiting in barrier (compute_wait_for_async_io)"
      CALL p_barrier(comm=p_comm_work)
      IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": barrier done (compute_wait_for_async_io)"
    END IF
  END SUBROUTINE compute_wait_for_async_io


  !-------------------------------------------------------------------------------------------------
  !> compute_start_async_io: Send a message to I/O PEs that they should start I/O
  !  The counterpart on the I/O side is async_io_wait_for_start
  !
  SUBROUTINE compute_start_async_io(jstep, output_pe_list, noutput_pe_list)
    INTEGER, INTENT(IN)          :: jstep
    INTEGER, INTENT(IN)          :: output_pe_list(:), noutput_pe_list
    ! local variables
    REAL(wp) :: msg(2)
    INTEGER  :: i

    IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": compute_start_async_io, jstep = ",jstep
    CALL p_barrier(comm=p_comm_work) ! make sure all are here
    msg(1) = REAL(msg_io_start,    wp)
    msg(2) = REAL(jstep,           wp)

    IF(p_pe_work==0) THEN

      ! When this subroutine is called, we have already proceeded to
      ! the next step. Send a "start message" to all I/O PEs which are
      ! due for output.
      DO i=1,noutput_pe_list
        IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": send signal to PE ",  output_pe_list(i)
        CALL p_isend(msg, output_pe_list(i), 0)
      END DO
      CALL p_wait()
    END IF
    IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": compute_start_async_io done."

  END SUBROUTINE compute_start_async_io


  !-------------------------------------------------------------------------------------------------
  !> compute_shutdown_async_io: Send a message to I/O PEs that they should shut down
  !  The counterpart on the I/O side is async_io_wait_for_start
  !
  SUBROUTINE compute_shutdown_async_io
    REAL(wp) :: msg(2)
    INTEGER  :: pe

    IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": compute_shutdown_async_io."
    CALL p_barrier(comm=p_comm_work) ! make sure all are here
    IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": compute_shutdown_async_io barrier done."
    msg(1) = REAL(msg_io_shutdown, wp)
    msg(2:) = 0._wp
    ! tell all I/O PEs about the shutdown
    IF(p_pe_work==0) THEN
      DO pe = p_io_pe0, (p_io_pe0+num_io_procs-1)
        CALL p_send(msg, pe, 0)
      END DO
    END IF
  END SUBROUTINE compute_shutdown_async_io

  !-------------------------------------------------------------------------------------------------
#endif

END MODULE mo_name_list_output
