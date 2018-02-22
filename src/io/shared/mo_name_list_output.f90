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
!! -------------------------------------------------------------------------
!!
!! The "namelist_output" module was originally written by Rainer
!! Johanni. Some data structures used therein are duplicates of those
!! created in the other parts of the model: In general, variable
!! fields are introduced in ICON through the "add_var" mechanism in
!! the module "shared/mo_var_list". This mechanism allocates "r_ptr"
!! POINTERs for REAL(wp) variable fields, see the data structure in
!! "t_var_list_element" (mo_var_list_element.f90). The "p_nh_state"
!! variables, for example, then point to the same location. In the
!! output, however, there exists a data structure "t_var_desc"
!! (variable descriptor) which also contains an "r_ptr" POINTER. This
!! also points to the original "r_ptr" location in memory.
!!
!! Exceptions and caveats for this described mechanism:
!!
!! - INTEGER fields are stored in "i_ptr" POINTERs.
!! - After gathering the output data, so-called "post-ops" are
!!   performed which modify the copied data (for example scaling from/to
!!   percent values).
!! - In asynchronous output mode, the "r_ptr" POINTERs are meaningless
!!   on those PEs which are dedicated for output. These are NULL
!!   pointers then.
!!
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
  USE mo_impl_constants,            ONLY: max_dom, SUCCESS, MAX_TIME_LEVELS, MAX_CHAR_LENGTH,       &
    &                                     ihs_ocean, BOUNDARY_MISSVAL
  USE mo_cdi_constants,             ONLY: GRID_REGULAR_LONLAT, GRID_UNSTRUCTURED_VERT,              &
    &                                     GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE
  USE mo_impl_constants_grf,        ONLY: grf_bdywidth_c
  USE mo_dynamics_config,           ONLY: iequations
  USE mo_cdi,                       ONLY: streamOpenWrite, FILETYPE_GRB2, streamDefTimestep, cdiEncodeTime, cdiEncodeDate, &
      &                                   CDI_UNDEFID, TSTEP_CONSTANT, FILETYPE_GRB, taxisDestroy, gridDestroy, &
      &                                   vlistDestroy, streamClose, streamWriteVarSlice, streamWriteVarSliceF, streamDefVlist, &
      &                                   streamSync, taxisDefVdate, taxisDefVtime, GRID_LONLAT, &
      &                                   streamOpenAppend, streamInqVlist, vlistInqTaxis, vlistNtsteps
  USE mo_util_cdi,                  ONLY: cdiGetStringError
  ! utility functions
  USE mo_io_units,                  ONLY: FILENAME_MAX, find_next_free_unit
  USE mo_exception,                 ONLY: finish, message, message_text
  USE mo_util_string,               ONLY: t_keyword_list, associate_keyword, with_keywords,         &
  &                                       int2string
  USE mo_dictionary,                ONLY: dict_finalize
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_write_output, ltimer,      &
    &                                     print_timer
  USE mo_name_list_output_gridinfo, ONLY: write_grid_info_grb2, GRID_INFO_NONE
  ! config
  USE mo_master_config,             ONLY: getModelBaseDir
  USE mo_grid_config,               ONLY: n_dom, l_limited_area
  USE mo_run_config,                ONLY: msg_level
  USE mo_io_config,                 ONLY: lkeep_in_sync,                   &
    &                                     config_lmask_boundary => lmask_boundary
  USE mo_gribout_config,            ONLY: gribout_config
  USE mo_parallel_config,           ONLY: p_test_run, use_dp_mpi2io, &
       num_io_procs, io_proc_chunk_size
  USE mo_name_list_output_config,   ONLY: use_async_name_list_io
  ! data types
  USE mo_var_metadata_types,        ONLY: t_var_metadata, POST_OP_SCALE, POST_OP_LUC
  USE mo_name_list_output_types,    ONLY: t_output_file, t_reorder_info, &
       msg_io_start, msg_io_done, msg_io_meteogram_flush, msg_io_shutdown, &
       all_events
  USE mo_output_event_types,        ONLY: t_sim_step_info, t_par_output_event
  ! parallelization
  USE mo_communication,             ONLY: exchange_data, t_comm_gather_pattern, idx_no, blk_no,     &
    &                                     idx_1d
  USE mo_mpi,                       ONLY: p_send, p_recv, p_barrier, stop_mpi,                      &
    &                                     p_mpi_wtime, p_irecv, p_wait, p_test, p_isend,            &
    &                                     p_comm_work, p_real_dp, p_real_sp, p_int,                 &
    &                                     my_process_is_stdio, my_process_is_mpi_test,              &
    &                                     my_process_is_mpi_workroot,                               &
    &                                     my_process_is_io, my_process_is_mpi_ioroot,               &
    &                                     process_mpi_all_test_id, process_mpi_all_workroot_id,     &
    &                                     num_work_procs, p_pe, p_pe_work, p_work_pe0, p_io_pe0,    &
    &                                     p_max
  ! calendar operations
  USE mtime,                        ONLY: datetime, newDatetime, deallocateDatetime, OPERATOR(-),   &
    &                                     timedelta, newTimedelta, deallocateTimedelta,             &
    &                                     MAX_DATETIME_STR_LEN
  ! output scheduling
  USE mo_output_event_handler,      ONLY: is_output_step, check_open_file, check_close_file,        &
    &                                     pass_output_step, get_current_filename,                   &
    &                                     get_current_date,                                         &
    &                                     is_output_step_complete, is_output_event_finished,        &
    &                                     check_write_readyfile, blocking_wait_for_irecvs
#ifndef NOMPI
  USE mo_output_event_handler,      ONLY: trigger_output_step_irecv
#endif
  USE mo_name_list_output_stats,    ONLY: set_reference_time, interval_start, interval_end,         &
    &                                     interval_write_psfile
  ! output initialization
  USE mo_name_list_output_init,     ONLY: init_name_list_output, setup_output_vlist,                &
    &                                     varnames_dict, out_varnames_dict,                         &
    &                                     output_file, patch_info, lonlat_info,                     &
    &                                     collect_requested_ipz_levels, create_vertical_axes
  USE mo_name_list_output_metadata, ONLY: metainfo_write_to_memwin, metainfo_get_from_memwin,       &
    &                                     metainfo_get_size, metainfo_get_timelevel
  USE mo_level_selection,           ONLY: create_mipz_level_selections
  USE mo_grib2_util,                ONLY: set_GRIB2_timedep_keys, set_GRIB2_timedep_local_keys
  ! model domain
  USE mo_model_domain,              ONLY: t_patch, p_patch
  USE mo_loopindices,               ONLY: get_indices_c
  ! post-ops

#ifndef __NO_ICON_ATMO__
  USE mo_post_op,                   ONLY: perform_post_op
  USE mo_meteogram_output,          ONLY: meteogram_init, meteogram_finalize, &
       meteogram_flush_file
  USE mo_meteogram_config,          ONLY: meteogram_output_config
  USE mo_intp_lonlat_types,         ONLY: lonlat_grids
#endif

  IMPLICIT NONE

  PRIVATE

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
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//"::open_output_file"
    CHARACTER(LEN=filename_max)    :: filename, filename_for_append
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText
    INTEGER                        :: tsID
    LOGICAL                        :: lexist, lappend
 
    ! open/append file: as this is a preliminary solution only, I do not try to
    ! streamline the conditionals
    filename            = TRIM(get_current_filename(of%out_event))
    filename_for_append = ''
    lappend             = .FALSE.

    ! check and reset filename, if data should be appended
    IF (split_output_filename(filename, filename_for_append, '_part_') > 0) THEN
      ! does the file to append to exist
      INQUIRE(file=TRIM(filename_for_append), exist=lexist)
      IF (lexist) THEN
        ! store the orginal allocated vlist (the handlers different) for later use with new files
        of%cdiVlistID_orig = of%cdiVlistID
        of%cdiTaxisID_orig = of%cdiTaxisID

        ! open for append
        of%cdiFileID       = streamOpenAppend(TRIM(filename_for_append))

        ! inquire the opened file for its associated vlist
        of%cdiVlistID      = streamInqVlist(of%cdiFileID)

        ! and time axis, the only components different to the previous model prepared vlist
        of%cdiTaxisID      = vlistInqTaxis(of%cdiVlistID)

        ! get the already stored number of time steps
        of%cdiTimeIndex    = vlistNtsteps(of%cdiVlistID)
        lappend            = .TRUE.
        of%appending       = .TRUE.
      ELSE
        IF (of%appending) THEN
          ! restore model internal vlist and time axis handler association
          of%cdiVlistID      = of%cdiVlistID_orig
          of%cdiVlistID_orig = CDI_UNDEFID
          of%cdiTaxisID      = of%cdiTaxisID_orig
          of%cdiTaxisID_orig = CDI_UNDEFID
        ENDIF
        ! file to append to does not exist that means we can use the name without part trailer
        filename = filename_for_append
        of%cdiFileID       = streamOpenWrite(TRIM(filename), of%output_type)
        of%appending       = .FALSE.
      ENDIF
    ELSE
      IF (of%appending) THEN
        ! restore model internal vlist and time axis handler association
        of%cdiVlistID      = of%cdiVlistID_orig
        of%cdiVlistID_orig = CDI_UNDEFID
        of%cdiTaxisID      = of%cdiTaxisID_orig
        of%cdiTaxisID_orig = CDI_UNDEFID
      ENDIF
        of%cdiFileID       = streamOpenWrite(TRIM(filename), of%output_type)
        of%appending       = .FALSE.
    ENDIF

    IF (of%cdiFileID < 0) THEN
      CALL cdiGetStringError(of%cdiFileID, cdiErrorText)
      WRITE(message_text,'(a)') TRIM(cdiErrorText)
      CALL message('',message_text,all_print=.TRUE.)
      CALL finish (routine, 'open failed on '//TRIM(filename))
    ELSE IF (msg_level >= 8) THEN
      IF (lappend) THEN
        CALL message (routine, 'to add more data, reopened '//TRIM(filename_for_append),all_print=.TRUE.)
      ELSE
        CALL message (routine, 'opened '//TRIM(filename),all_print=.TRUE.)
      END IF
    ENDIF

    IF (.NOT. lappend) THEN
      ! assign the vlist (which must have ben set before)
      CALL streamDefVlist(of%cdiFileID, of%cdiVlistID)
      ! set cdi internal time index to 0 for writing time slices in netCDF
      of%cdiTimeIndex = 0
    ENDIF

  CONTAINS

    FUNCTION split_output_filename(instring, string1, delim) RESULT (idx)
      INTEGER :: idx
      CHARACTER(len=*), INTENT(in)  :: instring ,delim
      CHARACTER(len=*), INTENT(out) :: string1
      CHARACTER(len=256) :: tmp_string

      tmp_string = TRIM(instring)

      idx = INDEX(tmp_string, delim)
      string1 = tmp_string(1:idx-1)

    END FUNCTION split_output_filename

  END SUBROUTINE open_output_file


  !------------------------------------------------------------------------------------------------
  !> Close all name_list files
  !
  SUBROUTINE close_name_list_output()
    ! local variables
    INTEGER :: i, ierror

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
#ifndef NOMPI
        IF(use_async_name_list_io .AND. .NOT. my_process_is_mpi_test()) THEN
          CALL mpi_win_free(output_file(i)%mem_win%mpi_win, ierror)
          IF (use_dp_mpi2io) THEN
            CALL mpi_free_mem(output_file(i)%mem_win%mem_ptr_dp, ierror)
          ELSE
            CALL mpi_free_mem(output_file(i)%mem_win%mem_ptr_sp, ierror)
          END IF
          CALL mpi_win_free(output_file(i)%mem_win%mpi_win_metainfo, ierror)
          CALL mpi_free_mem(output_file(i)%mem_win%mem_ptr_metainfo_pe0, ierror)
        END IF
#endif
        IF (output_file(i)%cdiFileID >= 0) THEN
          ! clean up level selection (if there is one):
          IF (ASSOCIATED(output_file(i)%level_selection)) THEN
            CALL output_file(i)%level_selection%finalize()
            DEALLOCATE(output_file(i)%level_selection)
            output_file(i)%level_selection => NULL()
          END IF
          CALL close_output_file(output_file(i))
          CALL destroy_output_vlist(output_file(i))
        END IF

        ! destroy vertical axes meta-data:
        CALL output_file(i)%verticalAxisList%finalize()
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
      CALL vlistDestroy(vlistID)
    ENDIF

    of%cdiVlistID      = CDI_UNDEFID
    of%cdiCellGridID   = CDI_UNDEFID
    of%cdiEdgeGridID   = CDI_UNDEFID
    of%cdiVertGridID   = CDI_UNDEFID
    of%cdiLonLatGridID = CDI_UNDEFID
    of%cdiTaxisID      = CDI_UNDEFID

  END SUBROUTINE destroy_output_vlist


  !------------------------------------------------------------------------------------------------
  !> Loop over all output_name_list's, write the ones for which output is due
  !  This routine also cares about opening the output files the first time
  !  and reopening the files after a certain number of steps.
  !
  SUBROUTINE write_name_list_output(jstep, opt_lhas_output)
    INTEGER,           INTENT(IN)   :: jstep             !< model step
    !> (Optional) Flag: .TRUE. if this async I/O PE has written during this step:
    LOGICAL, OPTIONAL, INTENT(OUT)  :: opt_lhas_output
    ! local variables
    CHARACTER(LEN=*), PARAMETER  :: routine = modname//"::write_name_list_output"
    INTEGER                           :: i, j, idate, itime, iret
    TYPE(datetime),           POINTER :: io_datetime => NULL()
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
      IF ((output_file(i)%io_proc_id == p_pe) .AND. (msg_level >= 8)) THEN
        WRITE(text,'(a,a,a,a,a,i0)') &
          & 'Output to ',TRIM(get_current_filename(output_file(i)%out_event)),        &
          & ' at simulation time ', TRIM(get_current_date(output_file(i)%out_event)), &
          & ' by PE ', p_pe
        CALL message(routine, text,all_print=.TRUE.)
      END IF

      IF (output_file(i)%io_proc_id == p_pe) THEN
        ! convert time stamp string into
        ! year/month/day/hour/minute/second values using the mtime
        ! library:
        io_datetime => newDatetime(TRIM(get_current_date(output_file(i)%out_event)))
        idate = cdiEncodeDate(INT(io_datetime%date%year),   &
          &                   INT(io_datetime%date%month),  &
          &                   INT(io_datetime%date%day))
        itime = cdiEncodeTime(INT(io_datetime%time%hour),   &
          &                   INT(io_datetime%time%minute), &
          &                   INT(io_datetime%time%second))
        CALL deallocateDatetime(io_datetime)
        CALL taxisDefVdate(output_file(i)%cdiTaxisID, idate)
        CALL taxisDefVtime(output_file(i)%cdiTaxisID, itime)
        iret = streamDefTimestep(output_file(i)%cdiFileId, output_file(i)%cdiTimeIndex)
        output_file(i)%cdiTimeIndex = output_file(i)%cdiTimeIndex + 1
      END IF

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
        IF (msg_level >= 8) THEN
          CALL message (routine, 'closed '//TRIM(get_current_filename(output_file(i)%out_event)),all_print=.TRUE.)
        END IF
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
        output_pe_list(list_idx) = output_file(i)%io_proc_id
      END IF

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
          HANDLE_COMPLETE_STEPS : DO WHILE (ASSOCIATED(ev))
            IF (is_output_step_complete(ev) .AND.  &
              & .NOT. is_output_event_finished(ev)) THEN
              !--- write ready file
              IF (check_write_readyfile(ev%output_event))  CALL write_ready_file(ev)
              ! launch a non-blocking request to all participating PEs to
              ! acknowledge the completion of the next output event
#ifndef NOMPI
              IF (use_async_name_list_io) THEN
                CALL trigger_output_step_irecv(ev)
              ELSE
#endif
                ev%output_event%i_event_step = ev%output_event%i_event_step + 1
#ifndef NOMPI
              END IF
#else
              ev => ev%next
#endif
            ELSE
              ev => ev%next
            END IF
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
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: dtime_string

    ! compute current forecast time (delta):
    mtime_date     => newDatetime(TRIM(get_current_date(ev)))
    mtime_begin    => newDatetime(TRIM(ev%output_event%event_data%sim_start))
    forecast_delta => newTimedelta("P01D")
    forecast_delta = mtime_date - mtime_begin

    WRITE (forecast_delta_str,'(4(i2.2))') forecast_delta%day, forecast_delta%hour, &
      &                                    forecast_delta%minute, forecast_delta%second
    WRITE (dtime_string,'(i4.4,2(i2.2),a,3(i2.2),a)')                                                 &
      &                      mtime_date%date%year, mtime_date%date%month, mtime_date%date%day, 'T',   &
      &                      mtime_date%time%hour, mtime_date%time%minute, mtime_date%time%second, 'Z'

    CALL deallocateDatetime(mtime_date)
    CALL deallocateDatetime(mtime_begin)
    CALL deallocateTimedelta(forecast_delta)

    ! substitute tokens in ready file name
    CALL associate_keyword("<path>",            TRIM(getModelBaseDir()),         keywords)
    CALL associate_keyword("<datetime>",        TRIM(get_current_date(ev)),      keywords)
    CALL associate_keyword("<ddhhmmss>",        TRIM(forecast_delta_str),        keywords)
    CALL associate_keyword("<datetime2>",       TRIM(dtime_string),              keywords)

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
    INTEGER,          PARAMETER                 :: iREAL_sp = 3

    INTEGER                                     :: tl, i_dom, i_log_dom, i, iv, jk, n_points, &
      &                                            nlevs, nindex, mpierr, lonlat_id,          &
      &                                            idata_type, lev_idx, lev
    INTEGER(i8)                                 :: ioff
    TYPE (t_var_metadata), POINTER              :: info
    TYPE(t_reorder_info),  POINTER              :: p_ri
    REAL(wp),          ALLOCATABLE              :: r_ptr(:,:,:)
    REAL(sp),          ALLOCATABLE              :: s_ptr(:,:,:)
    INTEGER,           ALLOCATABLE              :: i_ptr(:,:,:)
    REAL(wp),          ALLOCATABLE              :: r_out_recv(:)
    INTEGER,           ALLOCATABLE              :: r_out_int(:)
    REAL(sp),          ALLOCATABLE, TARGET      :: r_out_sp(:)
    REAL(dp),          ALLOCATABLE, TARGET      :: r_out_dp(:)
    TYPE(t_comm_gather_pattern), POINTER        :: p_pat
    LOGICAL                                     :: l_error
    LOGICAL                                     :: have_GRIB, lwrite_single_precision
    LOGICAL                                     :: var_ignore_level_selection
    INTEGER                                     :: nmiss    ! missing value indicator
    INTEGER                                     :: var_ref_pos, last_bdry_index
    REAL(wp)                                    :: missval
    INTEGER                                     :: rl_start, rl_end, i_nchdom, i_startblk, i_endblk, &
      &                                            i_startidx, i_endidx, local_idx, glb_idx, jc, jb
    TYPE(t_patch), POINTER                      :: ptr_patch

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

    DO iv = 1, of%num_vars
      ! Note that we provide the pointer "info_ptr" to the variable's
      ! info data object and not the modified copy "info".
      info => of%var_desc(iv)%info_ptr
      CALL metainfo_write_to_memwin(of%mem_win, iv, info)
    END DO

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

    ! "lmask_boundary": Some of the output fields are not updated with
    ! meaningful values in the vicinity of the lateral domain
    ! boundary. To avoid spurious data on these triangle cells (which
    ! could also spoil the GRIB range compression), the user may
    ! choose to set them to a "missing value". Implementation details:
    ! In the "synchronous" output mode, the implementation exploits
    ! the fact that all (global) indices for the lateral boundary
    ! region are ordered to the start of the data array. Therefore,
    ! only the computation of a limit index "last_bdry_index" is
    ! required to mask out the lateral points. In the asynchronous
    ! output mode, on the other hand, the compute processes possess
    ! only a portion of the output field and therefore need to loop
    ! over the lateral triangles block- and line-wise. This feature
    ! can be (de-)activated for specific variables through the
    ! "info%lmask_boundary" metadata flag. It also depends on a global
    ! namelist switch "io_nml/lmask_boundary" (LOGICAL, default:
    ! false).
    !
    ! Only for synchronous output mode: communicate the largest global
    ! index of the lateral boundary cells, if required:

    last_bdry_index = 0
    IF ( (.NOT.use_async_name_list_io .OR. my_process_is_mpi_test()) .AND.   &
      &  config_lmask_boundary )  THEN

      ptr_patch => p_patch(i_log_dom)
      rl_start   = 1
      rl_end     = grf_bdywidth_c
      i_nchdom   = MAX(1,ptr_patch%n_childdom)
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      glb_idx    = 0
      DO jb=i_startblk,i_endblk
        CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, rl_start, rl_end)
        DO jc=i_startidx, i_endidx
          local_idx = idx_1d(jc,jb)
          glb_idx   = MAX(glb_idx, ptr_patch%cells%decomp_info%glb_index(local_idx))
        END DO
      END DO
      last_bdry_index = p_max(glb_idx, p_comm_work)
    END IF

    ! ----------------------------------------------------
    ! Go over all name list variables for this output file
    ! ----------------------------------------------------
    DO iv = 1, of%num_vars

      info => of%var_desc(iv)%info

      ! inspect time-constant variables only if we are writing the
      ! first step in this file:
      IF ((info%isteptype == TSTEP_CONSTANT) .AND. .NOT. l_first_write) CYCLE

      ! Check if first dimension of array is nproma.
      ! Otherwise we got an array which is not suitable for this output scheme.
    ! IF(info%used_dimensions(1) /= nproma) &
    !   CALL finish(routine,'1st dim is not nproma: '//TRIM(info%name))

      idata_type = iUNKNOWN

      ! For time level dependent elements: set time level and check if
      ! time level is present:
        ! set a default time level (which is not used anyway, but must
        ! be a valid array subscript):
      tl = 1
#ifndef __NO_ICON_ATMO__
      IF (.NOT. ASSOCIATED(of%var_desc(iv)%r_ptr)  .AND. &
        & .NOT. ASSOCIATED(of%var_desc(iv)%s_ptr)  .AND. &
        & .NOT. ASSOCIATED(of%var_desc(iv)%i_ptr)) THEN
        tl = metainfo_get_timelevel(info,i_log_dom)
        IF(tl<=0 .OR. tl>max_time_levels) &
          CALL finish(routine, 'Illegal time level in nnow()/nnow_rcf()')
        ! Check if present
        IF (.NOT. ASSOCIATED(of%var_desc(iv)%tlev_rptr(tl)%p)   .AND.   &
          & .NOT. ASSOCIATED(of%var_desc(iv)%tlev_sptr(tl)%p)   .AND.   &
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
      ELSE IF (ASSOCIATED(of%var_desc(iv)%s_ptr) .OR.  &
        & ASSOCIATED(of%var_desc(iv)%tlev_sptr(tl)%p)) THEN
        idata_type = iREAL_sp
      ELSE IF (ASSOCIATED(of%var_desc(iv)%i_ptr) .OR.  &
        & ASSOCIATED(of%var_desc(iv)%tlev_iptr(tl)%p)) THEN
        idata_type = iINTEGER
      END IF

      SELECT CASE (info%ndims)
      CASE (1)
        IF (idata_type == iREAL)    ALLOCATE(r_ptr(info%used_dimensions(1),1,1))
        IF (idata_type == iREAL_sp) ALLOCATE(s_ptr(info%used_dimensions(1),1,1))
        IF (idata_type == iINTEGER) ALLOCATE(i_ptr(info%used_dimensions(1),1,1))

        IF (info%lcontained .AND. (info%var_ref_pos /= -1))  &
          & CALL finish(routine, "internal error")
        IF (ASSOCIATED(of%var_desc(iv)%r_ptr)) THEN
          r_ptr(:,1,1) = of%var_desc(iv)%r_ptr(:,1,1,1,1)
        ELSE IF (ASSOCIATED(of%var_desc(iv)%s_ptr)) THEN
          s_ptr(:,1,1) = of%var_desc(iv)%s_ptr(:,1,1,1,1)
        ELSE IF (ASSOCIATED(of%var_desc(iv)%i_ptr)) THEN
          i_ptr(:,1,1) = of%var_desc(iv)%i_ptr(:,1,1,1,1)
        ELSE
          CALL finish(routine, "Internal error!")
        ENDIF

      CASE (2)
        ! 2D fields: Make a 3D copy of the array
        IF (idata_type == iREAL)    ALLOCATE(r_ptr(info%used_dimensions(1),1,info%used_dimensions(2)))
        IF (idata_type == iREAL_sp) ALLOCATE(s_ptr(info%used_dimensions(1),1,info%used_dimensions(2)))
        IF (idata_type == iINTEGER) ALLOCATE(i_ptr(info%used_dimensions(1),1,info%used_dimensions(2)))

        var_ref_pos = 3
        IF (info%lcontained)  var_ref_pos = info%var_ref_pos

        IF (ASSOCIATED(of%var_desc(iv)%r_ptr)) THEN
          SELECT CASE(var_ref_pos)
          CASE (1)
            r_ptr(:,1,:) = of%var_desc(iv)%r_ptr(nindex,:,:,1,1)
          CASE (2)
            r_ptr(:,1,:) = of%var_desc(iv)%r_ptr(:,nindex,:,1,1)
          CASE (3)
            r_ptr(:,1,:) = of%var_desc(iv)%r_ptr(:,:,nindex,1,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
        ELSE IF (ASSOCIATED(of%var_desc(iv)%s_ptr)) THEN
          SELECT CASE(var_ref_pos)
          CASE (1)
            s_ptr(:,1,:) = of%var_desc(iv)%s_ptr(nindex,:,:,1,1)
          CASE (2)
            s_ptr(:,1,:) = of%var_desc(iv)%s_ptr(:,nindex,:,1,1)
          CASE (3)
            s_ptr(:,1,:) = of%var_desc(iv)%s_ptr(:,:,nindex,1,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
        ELSE IF (ASSOCIATED(of%var_desc(iv)%i_ptr)) THEN
          SELECT CASE(var_ref_pos)
          CASE (1)
            i_ptr(:,1,:) = of%var_desc(iv)%i_ptr(nindex,:,:,1,1)
          CASE (2)
            i_ptr(:,1,:) = of%var_desc(iv)%i_ptr(:,nindex,:,1,1)
          CASE (3)
            i_ptr(:,1,:) = of%var_desc(iv)%i_ptr(:,:,nindex,1,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
        ELSE IF (ASSOCIATED(of%var_desc(iv)%tlev_rptr(tl)%p)) THEN
          SELECT CASE(var_ref_pos)
          CASE (1)
            r_ptr(:,1,:) = of%var_desc(iv)%tlev_rptr(tl)%p(nindex,:,:,1,1)
          CASE (2)
            r_ptr(:,1,:) = of%var_desc(iv)%tlev_rptr(tl)%p(:,nindex,:,1,1)
          CASE (3)
            r_ptr(:,1,:) = of%var_desc(iv)%tlev_rptr(tl)%p(:,:,nindex,1,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
        ELSE IF (ASSOCIATED(of%var_desc(iv)%tlev_sptr(tl)%p)) THEN
          SELECT CASE(var_ref_pos)
          CASE (1)
            s_ptr(:,1,:) = of%var_desc(iv)%tlev_sptr(tl)%p(nindex,:,:,1,1)
          CASE (2)
            s_ptr(:,1,:) = of%var_desc(iv)%tlev_sptr(tl)%p(:,nindex,:,1,1)
          CASE (3)
            s_ptr(:,1,:) = of%var_desc(iv)%tlev_sptr(tl)%p(:,:,nindex,1,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
        ELSE IF (ASSOCIATED(of%var_desc(iv)%tlev_iptr(tl)%p)) THEN
          SELECT CASE(var_ref_pos)
          CASE (1)
            i_ptr(:,1,:) = of%var_desc(iv)%tlev_iptr(tl)%p(nindex,:,:,1,1)
          CASE (2)
            i_ptr(:,1,:) = of%var_desc(iv)%tlev_iptr(tl)%p(:,nindex,:,1,1)
          CASE (3)
            i_ptr(:,1,:) = of%var_desc(iv)%tlev_iptr(tl)%p(:,:,nindex,1,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
        ELSE
          CALL finish(routine, "Internal error!")
        ENDIF
      CASE (3)

        var_ref_pos = 4
        IF (info%lcontained)  var_ref_pos = info%var_ref_pos

        ! 3D fields: Here we could just set a pointer to the
        ! array... if there were no post-ops
        IF (idata_type == iREAL)    ALLOCATE(r_ptr(info%used_dimensions(1), &
          &                                        info%used_dimensions(2), &
          &                                        info%used_dimensions(3)))
        IF (idata_type == iREAL_sp) ALLOCATE(s_ptr(info%used_dimensions(1), &
          &                                        info%used_dimensions(2), &
          &                                        info%used_dimensions(3)))
        IF (idata_type == iINTEGER) ALLOCATE(i_ptr(info%used_dimensions(1), &
          &                                        info%used_dimensions(2), &
          &                                        info%used_dimensions(3)))

        IF(ASSOCIATED(of%var_desc(iv)%r_ptr)) THEN
          SELECT CASE(var_ref_pos)
          CASE (1)
            r_ptr = of%var_desc(iv)%r_ptr(nindex,:,:,:,1)
          CASE (2)
            r_ptr = of%var_desc(iv)%r_ptr(:,nindex,:,:,1)
          CASE (3)
            r_ptr = of%var_desc(iv)%r_ptr(:,:,nindex,:,1)
          CASE (4)
            r_ptr = of%var_desc(iv)%r_ptr(:,:,:,nindex,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
        ELSE IF (ASSOCIATED(of%var_desc(iv)%s_ptr)) THEN
          SELECT CASE(var_ref_pos)
          CASE (1)
            s_ptr = of%var_desc(iv)%s_ptr(nindex,:,:,:,1)
          CASE (2)
            s_ptr = of%var_desc(iv)%s_ptr(:,nindex,:,:,1)
          CASE (3)
            s_ptr = of%var_desc(iv)%s_ptr(:,:,nindex,:,1)
          CASE (4)
            s_ptr = of%var_desc(iv)%s_ptr(:,:,:,nindex,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
        ELSE IF (ASSOCIATED(of%var_desc(iv)%i_ptr)) THEN
          SELECT CASE(var_ref_pos)
          CASE (1)
            i_ptr = of%var_desc(iv)%i_ptr(nindex,:,:,:,1)
          CASE (2)
            i_ptr = of%var_desc(iv)%i_ptr(:,nindex,:,:,1)
          CASE (3)
            i_ptr = of%var_desc(iv)%i_ptr(:,:,nindex,:,1)
          CASE (4)
            i_ptr = of%var_desc(iv)%i_ptr(:,:,:,nindex,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
        ELSE IF (ASSOCIATED(of%var_desc(iv)%tlev_rptr(tl)%p)) THEN
          SELECT CASE(var_ref_pos)
          CASE (1)
            r_ptr = of%var_desc(iv)%tlev_rptr(tl)%p(nindex,:,:,:,1)
          CASE (2)
            r_ptr = of%var_desc(iv)%tlev_rptr(tl)%p(:,nindex,:,:,1)
          CASE (3)
            r_ptr = of%var_desc(iv)%tlev_rptr(tl)%p(:,:,nindex,:,1)
          CASE (4)
            r_ptr = of%var_desc(iv)%tlev_rptr(tl)%p(:,:,:,nindex,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
        ELSE IF (ASSOCIATED(of%var_desc(iv)%tlev_rptr(tl)%p)) THEN
          SELECT CASE(var_ref_pos)
          CASE (1)
            s_ptr = of%var_desc(iv)%tlev_sptr(tl)%p(nindex,:,:,:,1)
          CASE (2)
            s_ptr = of%var_desc(iv)%tlev_sptr(tl)%p(:,nindex,:,:,1)
          CASE (3)
            s_ptr = of%var_desc(iv)%tlev_sptr(tl)%p(:,:,nindex,:,1)
          CASE (4)
            s_ptr = of%var_desc(iv)%tlev_sptr(tl)%p(:,:,:,nindex,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
        ELSE IF (ASSOCIATED(of%var_desc(iv)%tlev_iptr(tl)%p)) THEN
          SELECT CASE(var_ref_pos)
          CASE (1)
            i_ptr = of%var_desc(iv)%tlev_iptr(tl)%p(nindex,:,:,:,1)
          CASE (2)
            i_ptr = of%var_desc(iv)%tlev_iptr(tl)%p(:,nindex,:,:,1)
          CASE (3)
            i_ptr = of%var_desc(iv)%tlev_iptr(tl)%p(:,:,nindex,:,1)
          CASE (4)
            i_ptr = of%var_desc(iv)%tlev_iptr(tl)%p(:,:,:,nindex,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
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

#ifndef __NO_ICON_ATMO__
      IF ( ANY((/POST_OP_SCALE, POST_OP_LUC/) == of%var_desc(iv)%info%post_op%ipost_op_type) ) THEN
        IF (idata_type == iREAL) THEN
          CALL perform_post_op(of%var_desc(iv)%info%post_op, r_ptr)
        ELSE IF (idata_type == iREAL_sp) THEN
          CALL perform_post_op(of%var_desc(iv)%info%post_op, s_ptr)
        ELSE IF (idata_type == iINTEGER) THEN
          CALL perform_post_op(of%var_desc(iv)%info%post_op, i_ptr)
        ENDIF
      END IF
#endif

      var_ignore_level_selection = .FALSE.
      IF(info%ndims < 3) THEN
        nlevs = 1
      ELSE
        ! handle the case that a few levels have been selected out of
        ! the total number of levels:
        IF (ASSOCIATED(of%level_selection)) THEN
          nlevs = 0
          ! Sometimes the user mixes level-selected variables with
          ! other fields on other z-axes (e.g. soil fields) in the
          ! output namelist. We try to catch this "wrong" user input
          ! here and handle it in the following way: if the current
          ! variable does not have one (or more) of the requested
          ! levels, then we completely ignore the level selection for
          ! this variable.
          !
          ! (... but note that we accept (nlevs+1) for an nlevs variable.)
          CHECK_LOOP : DO jk=1,MIN(of%level_selection%n_selected, info%used_dimensions(2))
            IF ((of%level_selection%global_idx(jk) < 1) .OR.  &
              & (of%level_selection%global_idx(jk) > (info%used_dimensions(2)+1))) THEN
              var_ignore_level_selection = .TRUE.
              IF (my_process_is_stdio() .AND. (msg_level >= 15)) &
                &   WRITE (0,*) "warning: ignoring level selection for variable ", TRIM(info%name)
              nlevs = info%used_dimensions(2)
              EXIT CHECK_LOOP
            ELSE
              IF ((of%level_selection%global_idx(jk) >= 1) .AND.  &
                & (of%level_selection%global_idx(jk) <= info%used_dimensions(2))) THEN
                nlevs = nlevs + 1
              END IF
            END IF
          END DO CHECK_LOOP
        ELSE
          nlevs = info%used_dimensions(2)
        END IF
      ENDIF

      ! Get pointer to appropriate reorder_info
      SELECT CASE (info%hgrid)
      CASE (GRID_UNSTRUCTURED_CELL)
        p_ri  => patch_info(i_dom)%cells
        p_pat => patch_info(i_dom)%p_pat_c
      CASE (GRID_LONLAT)
      CASE (GRID_UNSTRUCTURED_EDGE)
        p_ri  => patch_info(i_dom)%edges
        p_pat => patch_info(i_dom)%p_pat_e
      CASE (GRID_UNSTRUCTURED_VERT)
        p_ri  => patch_info(i_dom)%verts
        p_pat => patch_info(i_dom)%p_pat_v
#ifndef __NO_ICON_ATMO__
      CASE (GRID_REGULAR_LONLAT)
        lonlat_id = info%hor_interp%lonlat_id
        p_ri  => lonlat_info(lonlat_id, i_log_dom)%ri
        p_pat => lonlat_grids%list(lonlat_id)%p_pat(i_log_dom)
#endif
      CASE default
        CALL finish(routine,'unknown grid type')
      END SELECT

      IF (.NOT.use_async_name_list_io .OR. my_process_is_mpi_test()) THEN

        IF (info%hgrid == GRID_LONLAT) THEN
          n_points = 1
        ELSE
          n_points = p_ri%n_glb
        END IF

        have_GRIB = (of%output_type == FILETYPE_GRB) .OR. (of%output_type == FILETYPE_GRB2)
        lwrite_single_precision = (.NOT. use_dp_mpi2io) .AND. (.NOT. have_GRIB)

        IF (idata_type == iREAL) THEN
          ALLOCATE(r_out_dp(MERGE(n_points, 0, my_process_is_mpi_workroot())))
        END IF
        IF ((idata_type == iREAL_sp) .OR. lwrite_single_precision) THEN
          ALLOCATE(r_out_sp(MERGE(n_points, 0, my_process_is_mpi_workroot())))
        END IF
        IF (idata_type == iINTEGER) THEN
          ALLOCATE(r_out_int(MERGE(n_points, 0, my_process_is_mpi_workroot())))
        END IF

        IF(my_process_is_mpi_workroot()) THEN

          IF (my_process_is_mpi_test()) THEN

            IF (p_test_run .AND. use_dp_mpi2io) ALLOCATE(r_out_recv(n_points))

          ELSE

            ! Set some GRIB2 keys that may have changed during simulation.
            ! Note that (for synchronous output mode) we provide the
            ! pointer "info_ptr" to the variable's info data object and
            ! not the modified copy "info".
            IF  (of%output_type == FILETYPE_GRB2) THEN
              CALL set_GRIB2_timedep_keys( &
                & of%cdiFileID, info%cdiVarID, of%var_desc(iv)%info_ptr, &
                & of%out_event%output_event%event_data%sim_start,        &
                & get_current_date(of%out_event))
              CALL set_GRIB2_timedep_local_keys(of%cdiFileID, info%cdiVarID, &
                & gribout_config(of%phys_patch_id) )
            END IF
          END IF ! my_process_is_mpi_test()
        END IF ! my_process_is_mpi_workroot()


        ! set missval flag, if applicable
        !
        IF ( ( of%var_desc(iv)%info%lmiss .OR.                                            &
          &    ( of%var_desc(iv)%info%lmask_boundary .AND. config_lmask_boundary) ) .AND. &
          &  ( last_bdry_index >  0 ) ) THEN
          nmiss = 1
        ELSE
          nmiss = 0
        ENDIF

        ! For all levels (this needs to be done level-wise in order to reduce
        !                 memory consumption)
        DO lev = 1, nlevs
          ! -------------------
          ! No asynchronous I/O
          ! -------------------
          !
          ! gather the array on stdio PE and write it out there

          IF ( n_points == 1 ) THEN
            IF (my_process_is_mpi_workroot()) THEN
              !write(0,*)'#--- n_points:',n_points,'idata_type:',idata_type,'iREAL:',iREAL,'iINTEGER:',iINTEGER
              IF      (idata_type == iREAL ) THEN
                r_out_dp(:)  = r_ptr(:,1,1)
              ELSE IF (idata_type == iREAL_sp ) THEN
                r_out_sp(:)  = s_ptr(:,1,1)
              ELSE IF (idata_type == iINTEGER) THEN
                r_out_int(:) = i_ptr(:,1,1)
              END IF
            END IF
          ELSE ! n_points
            IF (idata_type == iREAL) THEN
              r_out_dp(:)  = 0._wp

              lev_idx = lev
              ! handle the case that a few levels have been selected out of
              ! the total number of levels:
              IF (      ASSOCIATED(of%level_selection)   .AND. &
                & (.NOT. var_ignore_level_selection)     .AND. &
                & (info%ndims > 2)) THEN
                lev_idx = of%level_selection%global_idx(lev_idx)
              END IF
              CALL exchange_data(in_array=r_ptr(:,lev_idx,:),                 &
                &                out_array=r_out_dp(:), gather_pattern=p_pat, &
                &                fill_value = BOUNDARY_MISSVAL)

            ELSE IF (idata_type == iREAL_sp) THEN
              r_out_sp(:)  = 0._wp

              lev_idx = lev
              ! handle the case that a few levels have been selected out of
              ! the total number of levels:
              IF (      ASSOCIATED(of%level_selection)   .AND. &
                & (.NOT. var_ignore_level_selection)     .AND. &
                & (info%ndims > 2)) THEN
                lev_idx = of%level_selection%global_idx(lev_idx)
              END IF
              CALL exchange_data(in_array=s_ptr(:,lev_idx,:),                 &
                &                out_array=r_out_sp(:), gather_pattern=p_pat)
              ! FIXME: Implement and use fill_value!

            ELSE IF (idata_type == iINTEGER) THEN
              r_out_int(:) = 0

              lev_idx = lev
              ! handle the case that a few levels have been selected out of
              ! the total number of levels:
              IF (      ASSOCIATED(of%level_selection)   .AND. &
                & (.NOT. var_ignore_level_selection)     .AND. &
                & (info%ndims > 2)) THEN
                lev_idx = of%level_selection%global_idx(lev_idx)
              END IF
              CALL exchange_data(in_array=i_ptr(:,lev_idx,:),                  &
                &                out_array=r_out_int(:), gather_pattern=p_pat)
              ! FIXME: Implement and use fill_value!

            END IF
          END IF ! n_points

          IF(my_process_is_mpi_workroot()) THEN

            SELECT CASE(idata_type)
            CASE(iREAL)
              !
              ! If single precision output is desired, we need to
              ! perform a type conversion:
              IF ( lwrite_single_precision ) THEN
                r_out_sp(:) = REAL(r_out_dp(:), sp)
              ENDIF

            CASE(iREAL_sp)
              !
              IF ( .NOT. lwrite_single_precision ) THEN
                r_out_dp(:) = REAL(r_out_sp(:), dp)
              ENDIF

            CASE(iINTEGER)
              !
              IF ( lwrite_single_precision ) THEN
                r_out_sp(:) = REAL(r_out_int(:), sp)
              ELSE
                r_out_dp(:) = REAL(r_out_int(:), dp)
              ENDIF
            END SELECT

            ! If required, set lateral boundary points to missing
            ! value. Note that this modifies only the output buffer!
            IF ( info%lmask_boundary                   .AND. &
              & (info%hgrid == GRID_UNSTRUCTURED_CELL) .AND. &
              & config_lmask_boundary ) THEN
              missval = BOUNDARY_MISSVAL
              IF (info%lmiss)  missval = info%missval%rval
              
              IF ( lwrite_single_precision ) THEN
                r_out_sp(1:last_bdry_index) = missval
              ELSE
                r_out_dp(1:last_bdry_index) = missval
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
                CALL p_recv(r_out_recv, process_mpi_all_workroot_id, 1)
                ! check for correctness
                l_error = .FALSE.
                DO i = 1, n_points
                  IF (r_out_recv(i) /= r_out_dp(i)) THEN
                    ! do detailed print-out only for "large" errors:
                    IF (ABS(r_out_recv(i) - r_out_dp(i)) > SYNC_ERROR_PRINT_TOL) THEN
                      WRITE (0,*) 'Sync error test PE/worker PEs for ', TRIM(info%name)
                      WRITE (0,*) "global pos (", idx_no(i), ",", blk_no(i),")"
                      WRITE (0,*) "level", lev, "/", nlevs
                      WRITE (0,*) "vals: ", r_out_recv(i), r_out_dp(i)
                      l_error = .TRUE.
                    END IF
                  END IF
                ENDDO
                IF (l_error)   CALL finish(routine,"Sync error!")
              ENDIF
            ENDIF

          ENDIF

          ! ----------
          ! write data
          ! ----------

          IF (my_process_is_mpi_workroot() .AND. .NOT. my_process_is_mpi_test()) THEN
            IF (.NOT. lwrite_single_precision) THEN
              CALL streamWriteVarSlice (of%cdiFileID, info%cdiVarID, lev-1, r_out_dp(:), nmiss)
            ELSE
              CALL streamWriteVarSliceF(of%cdiFileID, info%cdiVarID, lev-1, r_out_sp(:), nmiss)
            ENDIF
          END IF

        END DO ! lev = 1, nlevs

        IF (my_process_is_mpi_workroot() .AND. lkeep_in_sync .AND. &
          & .NOT. my_process_is_mpi_test()) CALL streamSync(of%cdiFileID)

        ! clean up
        IF (ALLOCATED(r_out_int)) DEALLOCATE(r_out_int)
        IF (ALLOCATED(r_out_sp))  DEALLOCATE(r_out_sp)
        IF (ALLOCATED(r_out_dp))  DEALLOCATE(r_out_dp)
        IF (ALLOCATED(r_out_recv)) DEALLOCATE(r_out_recv)

#ifndef NOMPI

      ELSE

        ! ------------------------
        ! Asynchronous I/O is used
        ! ------------------------

        ! just copy the OWN DATA points to the memory window
        DO jk = 1, nlevs
          ! handle the case that a few levels have been selected out of
          ! the total number of levels:
          IF (      ASSOCIATED(of%level_selection)   .AND. &
            & (.NOT. var_ignore_level_selection)     .AND. &
            & (info%ndims > 2)) THEN
            lev_idx = of%level_selection%global_idx(jk)
          ELSE
            lev_idx = jk
          END IF

          IF (use_dp_mpi2io) THEN
            IF (idata_type == iREAL) THEN
              DO i = 1, p_ri%n_own
                of%mem_win%mem_ptr_dp(ioff+INT(i,i8)) = &
                  & REAL(r_ptr(p_ri%own_idx(i),lev_idx,p_ri%own_blk(i)),dp)
              ENDDO
            END IF
            IF (idata_type == iREAL_sp) THEN
              DO i = 1, p_ri%n_own
                of%mem_win%mem_ptr_dp(ioff+INT(i,i8)) = &
                  & REAL(s_ptr(p_ri%own_idx(i),lev_idx,p_ri%own_blk(i)),dp)
              ENDDO
            END IF
            IF (idata_type == iINTEGER) THEN
              DO i = 1, p_ri%n_own
                of%mem_win%mem_ptr_dp(ioff+INT(i,i8)) = &
                  & REAL(i_ptr(p_ri%own_idx(i),lev_idx,p_ri%own_blk(i)),dp)
              ENDDO
            END IF

            ! If required, set lateral boundary points to missing
            ! value. Note that this modifies only the output buffer!
            IF ( info%lmask_boundary                    .AND. &
              &  (info%hgrid == GRID_UNSTRUCTURED_CELL) .AND. &
              &  config_lmask_boundary ) THEN
              missval = BOUNDARY_MISSVAL
              IF (info%lmiss) THEN
                IF (idata_type == iREAL) THEN
                  missval = info%missval%rval
                ELSE IF (idata_type == iINTEGER) THEN
                  missval = REAL(info%missval%ival,dp)
                END IF
              END IF
              ptr_patch => p_patch(i_log_dom)
              rl_start   = 1
              rl_end     = grf_bdywidth_c
              i_nchdom   = MAX(1,ptr_patch%n_childdom)
              i_startblk = ptr_patch%cells%start_blk(rl_start,1)
              i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
              CALL get_indices_c(ptr_patch, i_endblk, i_startblk, i_endblk, &
                i_startidx, i_endidx, rl_start, rl_end)
              DO i = 1, p_ri%n_own
                IF ( (p_ri%own_blk(i) < i_endblk) .OR. &
                  &  ((p_ri%own_blk(i) == i_endblk) .AND. &
                  &   (p_ri%own_idx(i) <= i_endidx)) ) THEN
                  of%mem_win%mem_ptr_dp(ioff+INT(i,i8)) = missval
                END IF
              END DO
            END IF

          ELSE
            IF (idata_type == iREAL) THEN
              DO i = 1, p_ri%n_own
                of%mem_win%mem_ptr_sp(ioff+INT(i,i8)) = &
                  & REAL(r_ptr(p_ri%own_idx(i),lev_idx,p_ri%own_blk(i)),sp)
              ENDDO
            END IF
            IF (idata_type == iREAL_sp) THEN
              DO i = 1, p_ri%n_own
                of%mem_win%mem_ptr_sp(ioff+INT(i,i8)) = &
                  & s_ptr(p_ri%own_idx(i),lev_idx,p_ri%own_blk(i))
              ENDDO
            END IF
            IF (idata_type == iINTEGER) THEN
              DO i = 1, p_ri%n_own
                of%mem_win%mem_ptr_sp(ioff+INT(i,i8)) = &
                  & REAL(i_ptr(p_ri%own_idx(i),lev_idx,p_ri%own_blk(i)),sp)
              ENDDO
            END IF

            ! If required, set lateral boundary points to missing
            ! value. Note that this modifies only the output buffer!
            IF ( info%lmask_boundary                    .AND. &
              &  (info%hgrid == GRID_UNSTRUCTURED_CELL) .AND. &
              &  config_lmask_boundary ) THEN
              missval = BOUNDARY_MISSVAL
              IF (info%lmiss) THEN
                IF (idata_type == iREAL) THEN
                  missval = info%missval%rval
                ELSE IF (idata_type == iINTEGER) THEN
                  missval = REAL(info%missval%ival,sp)
                END IF
              END IF
              ptr_patch => p_patch(i_log_dom)
              rl_start   = 1
              rl_end     = grf_bdywidth_c
              i_nchdom   = MAX(1,ptr_patch%n_childdom)
              i_startblk = ptr_patch%cells%start_blk(rl_start,1)
              i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
              CALL get_indices_c(ptr_patch, i_endblk, i_startblk, i_endblk, &
                i_startidx, i_endidx, rl_start, rl_end)
              DO i = 1, p_ri%n_own
                IF ( (p_ri%own_blk(i) < i_endblk) .OR. &
                  &  ((p_ri%own_blk(i) == i_endblk) .AND. &
                  &   (p_ri%own_idx(i) <= i_endidx)) ) THEN
                  of%mem_win%mem_ptr_sp(ioff+INT(i,i8)) = REAL(missval,sp)
                END IF
              END DO
            END IF

          END IF
          ioff = ioff + INT(p_ri%n_own,i8)
        END DO ! nlevs

#endif

      END IF

      ! clean up
      IF (ALLOCATED(r_ptr)) DEALLOCATE(r_ptr)
      IF (ALLOCATED(s_ptr)) DEALLOCATE(s_ptr)
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

  SUBROUTINE name_list_io_main_proc(sim_step_info)
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info
  END SUBROUTINE name_list_io_main_proc

#else


  !-------------------------------------------------------------------------------------------------
  !> Main routine for I/O PEs.
  !  Please note that this routine never returns.
  !
  SUBROUTINE name_list_io_main_proc(sim_step_info)
    !> Data structure containing all necessary data for mapping an
    !  output time stamp onto a corresponding simulation step index.
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info
    ! local variables:

#ifndef __NO_ICON_ATMO__
    LOGICAL             :: l_complete, lhas_output, &
      &                    lset_timers_for_idle_pe
    INTEGER             :: jg, jstep, i, action
    TYPE(t_par_output_event), POINTER :: ev

    ! define initial time stamp used as reference for output statistics
    CALL set_reference_time()

    ! Initialize name list output, this is a collective call for all PEs
    CALL init_name_list_output(sim_step_info)

    ! setup of meteogram output
    DO jg =1,n_dom
      IF (meteogram_output_config(jg)%lenabled) THEN
        CALL meteogram_init(meteogram_output_config(jg), jg,    &
          &                 grid_uuid=patch_info(jg)%grid_uuid, &
          &                 number_of_grid_used=patch_info(jg)%number_of_grid_used)
      END IF
    END DO


    ! Append the chosen p-levels, z-levels, i-levels to the levels
    ! sets for the corresponding domains:
    !
    ! Note that on pure I/O PEs we must call this *after* the
    ! "init_name_list_output", since some values (log_dom_id) are
    ! reuqired which are communicated there.
    CALL collect_requested_ipz_levels()
    IF (iequations/=ihs_ocean) THEN ! atm
      CALL create_mipz_level_selections(output_file)
    END IF
    CALL create_vertical_axes(output_file)

    ! Tell the compute PEs that we are ready to work
    IF (ANY(output_file(:)%io_proc_id == p_pe)) THEN
      CALL async_io_send_handshake(0)
    END IF

    ! Enter I/O loop
    ! skip loop, if this output PE is idle:
    IF (     ANY(                  output_file(:)%io_proc_id == p_pe) &
        .OR. ANY(meteogram_output_config(1:n_dom)%io_proc_id == p_pe)) THEN
      DO

        ! Wait for a message from the compute PEs to start
        CALL async_io_wait_for_start(action, jstep)

        IF(action == msg_io_shutdown) EXIT ! leave loop, we are done

        IF (action == msg_io_start) THEN
          ! perform I/O
          CALL write_name_list_output(jstep, opt_lhas_output=lhas_output)

          ! Inform compute PEs that we are done, if this I/O PE has
          ! written output:
          IF (lhas_output)  CALL async_io_send_handshake(jstep)

          ! Handle final pending "output step completed" messages: After
          ! all participating I/O PE's have acknowledged the completion of
          ! their write processes, we trigger a "ready file" on the first
          ! I/O PE.
          IF (.NOT. my_process_is_mpi_test()  .AND.  &
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
        ELSE IF (action == msg_io_meteogram_flush) THEN
          CALL meteogram_flush_file(jstep)
        END IF

      ENDDO
    ENDIF
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

    CALL interval_write_psfile("output_schedule.ps", "Output Timings", &
      &                        int2string(p_pe,'(i0)'), p_comm_work)

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
#ifndef NOMPI
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
    INTEGER                        :: nmiss    ! missing value indicator
    INTEGER                        :: ichunk, nchunks, chunk_start, chunk_end, &
      &                               this_chunk_nlevs, ilev

    !-- for timing
    CHARACTER(len=10)              :: ctime
    REAL(dp)                       :: t_get, t_write, t_copy, t_0, mb_get, mb_wr

  !------------------------------------------------------------------------------------------------

    CALL date_and_time(TIME=ctime)
    IF (msg_level >= 8) THEN
      WRITE (0, '(a,i0,a)') '#################### I/O PE ',p_pe,' starting I/O at '//ctime
    END IF
    CALL interval_start(TRIM(get_current_filename(of%out_event)))

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

    ! if no valid io_proc_chunk_size has been set by the parallel name list
    IF (io_proc_chunk_size <= 0) THEN
      io_proc_chunk_size = nlev_max
    ELSE
      io_proc_chunk_size = MIN(nlev_max, io_proc_chunk_size)
    END IF

    IF (use_dp_mpi2io) THEN
      ALLOCATE(var1_dp(nval*io_proc_chunk_size), STAT=ierrstat)
    ELSE
      ALLOCATE(var1_sp(nval*io_proc_chunk_size), STAT=ierrstat)
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

      ! WRITE (0,*) ">>>>>>>>>> ", info%name
      ! get also an update for this variable's meta-info (separate object)
      CALL metainfo_get_from_memwin(bufr_metainfo, iv, updated_info)

      IF ( of%output_type == FILETYPE_GRB2 ) THEN
        CALL set_GRIB2_timedep_keys(of%cdiFileID, info%cdiVarID, updated_info,       &
          &                         of%out_event%output_event%event_data%sim_start,  &
          &                         get_current_date(of%out_event))
        CALL set_GRIB2_timedep_local_keys(of%cdiFileID, info%cdiVarID,     &
          &                               gribout_config(of%phys_patch_id) )
      END IF

      ! Set missval flag, if applicable
      !
      ! Missing value masks are available in GRIB output format
      ! only. A missing value might be set by the user (info%lmiss) or
      ! automatically on nest boundary regions.
      !
      IF ( (of%output_type == FILETYPE_GRB) .OR. (of%output_type == FILETYPE_GRB2) ) THEN
        IF ( info%lmiss .OR.                                            &
          &  ( info%lmask_boundary    .AND. &
          &    config_lmask_boundary  .AND. &
          &    ((i_log_dom > 1) .OR. l_limited_area) ) ) THEN
          nmiss = 1
        ELSE
          nmiss = 0
        ENDIF
      ELSE  ! i.e. NETCDF
        nmiss = 0
      ENDIF

      ! inspect time-constant variables only if we are writing the
      ! first step in this file:
      IF ((info%isteptype == TSTEP_CONSTANT) .AND. .NOT. l_first_write) CYCLE

      IF(info%ndims == 2) THEN
        nlevs = 1
      ELSE
        ! handle the case that a few levels have been selected out of
        ! the total number of levels:
        IF (ASSOCIATED(of%level_selection)) THEN
          ! count the no. of selected levels for this variable:
          nlevs = 0
          CHECK_LOOP : DO jk=1,MIN(of%level_selection%n_selected, info%used_dimensions(2))
            IF ((of%level_selection%global_idx(jk) < 1) .OR.  &
              & (of%level_selection%global_idx(jk) > (info%used_dimensions(2)+1))) THEN
              nlevs = info%used_dimensions(2)
              EXIT CHECK_LOOP
            ELSE
              IF ((of%level_selection%global_idx(jk) >= 1) .AND.  &
                & (of%level_selection%global_idx(jk) <= info%used_dimensions(2))) THEN
                nlevs = nlevs + 1
              END IF
            END IF
          END DO CHECK_LOOP
        ELSE
          nlevs = info%used_dimensions(2)
        END IF
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

      ! no. of chunks of levels (each of size "io_proc_chunk_size"):
      nchunks = (nlevs-1)/io_proc_chunk_size + 1
      ! loop over all chunks (of levels)
      DO ichunk=1,nchunks

        chunk_start       = (ichunk-1)*io_proc_chunk_size + 1
        chunk_end         = MIN(chunk_start+io_proc_chunk_size-1, nlevs)
        this_chunk_nlevs  = (chunk_end - chunk_start + 1)

        ! Retrieve part of variable from every worker PE using MPI_Get
        nv_off  = 0
        DO np = 0, num_work_procs-1

          IF(p_ri%pe_own(np) == 0) CYCLE

          ! Number of words to transfer
          nval = p_ri%pe_own(np) * this_chunk_nlevs

          t_0 = p_mpi_wtime()
          CALL MPI_Win_lock(MPI_LOCK_SHARED, np, MPI_MODE_NOCHECK, &
            &               of%mem_win%mpi_win, mpierr)

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
          nval            = p_ri%pe_own(np)*this_chunk_nlevs
          nv_off          = nv_off + nval
          nv_off_np(np+1) = nv_off_np(np) + p_ri%pe_own(np)
        END DO

        DO ilev=chunk_start, chunk_end
          t_0 = p_mpi_wtime() ! performance measurement

          IF (use_dp_mpi2io) THEN
            var3_dp(:) = 0._wp

!$OMP PARALLEL
!$OMP DO PRIVATE(dst_start, dst_end, src_start, src_end)
            DO np = 0, num_work_procs-1
              dst_start = nv_off_np(np)+1
              dst_end   = nv_off_np(np+1)
              src_start = voff(np)+1
              src_end   = voff(np)+p_ri%pe_own(np)
              voff(np)  = src_end
              var3_dp(p_ri%reorder_index(dst_start:dst_end)) = &
                var1_dp(src_start:src_end)
            ENDDO
!$OMP END DO
!$OMP END PARALLEL
          ELSE

            IF (have_GRIB) THEN
              var3_dp(:) = 0._wp

              ! ECMWF GRIB-API/CDI has only a double precision interface at the
              ! date of coding this
!$OMP PARALLEL
!$OMP DO PRIVATE(dst_start, dst_end, src_start, src_end)
              DO np = 0, num_work_procs-1
                dst_start = nv_off_np(np)+1
                dst_end   = nv_off_np(np+1)
                src_start = voff(np)+1
                src_end   = voff(np)+p_ri%pe_own(np)
                voff(np)  = src_end
                var3_dp(p_ri%reorder_index(dst_start:dst_end)) = &
                  REAL(var1_sp(src_start:src_end), dp)
              ENDDO
!$OMP END DO
!$OMP END PARALLEL
            ELSE
              var3_sp(:) = 0._sp
!$OMP PARALLEL
!$OMP DO PRIVATE(dst_start, dst_end, src_start, src_end)
              DO np = 0, num_work_procs-1
                dst_start = nv_off_np(np)+1
                dst_end   = nv_off_np(np+1)
                src_start = voff(np)+1
                src_end   = voff(np)+p_ri%pe_own(np)
                voff(np)  = src_end
                var3_sp(p_ri%reorder_index(dst_start:dst_end)) = &
                  var1_sp(src_start:src_end)
              ENDDO
!$OMP END DO
!$OMP END PARALLEL
            END IF
          ENDIF
          t_copy = t_copy + p_mpi_wtime() - t_0 ! performance measurement
          ! Write calls (via CDIs) of the asynchronous I/O PEs:
          t_0 = p_mpi_wtime() ! performance measurement

          IF (use_dp_mpi2io .OR. have_GRIB) THEN
            CALL streamWriteVarSlice(of%cdiFileID, info%cdiVarID, ilev-1, var3_dp, nmiss)
            mb_wr = mb_wr + REAL(SIZE(var3_dp), wp)
          ELSE
            CALL streamWriteVarSliceF(of%cdiFileID, info%cdiVarID, ilev-1, var3_sp, nmiss)
            mb_wr = mb_wr + REAL(SIZE(var3_sp),wp)
          ENDIF
          t_write = t_write + p_mpi_wtime() - t_0 ! performance measurement

        ENDDO ! ilev

      ENDDO ! chunk loop

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
    IF (msg_level >= 8) THEN
      WRITE (0, '(a,i0,a)') '#################### I/O PE ',p_pe,' done at '//ctime
    END IF
    CALL interval_end(TRIM(get_current_filename(of%out_event)))

    ! Convert mb_get/mb_wr to MB
    IF (use_dp_mpi2io) THEN
      mb_get = mb_get*8*1.d-6
    ELSE
      mb_get = mb_get*4*1.d-6
    ENDIF
    mb_wr  = mb_wr*4*1.d-6 ! 4 byte since dp output is implicitly converted to sp
    ! writing this message causes a runtime error on the NEC because formatted output to stdio/stderr is limited to 132 chars
    IF (msg_level >= 12) THEN
      WRITE (0,'(10(a,f10.3))') &  ! remark: CALL message does not work here because it writes only on PE0
           & ' Got ',mb_get,' MB, time get: ',t_get,' s [',mb_get/MAX(1.e-6_wp,t_get), &
           & ' MB/s], time write: ',t_write,' s [',mb_wr/MAX(1.e-6_wp,t_write),        &
           & ' MB/s], times copy: ',t_copy,' s'
   !   CALL message('',message_text)
    ENDIF

    DEALLOCATE(bufr_metainfo, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  END SUBROUTINE io_proc_write_name_list
#endif

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
#ifndef NOMPI
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
#endif

  !-------------------------------------------------------------------------------------------------
  !> async_io_wait_for_start: Wait for a message from work PEs that we
  !  should start I/O or finish.  The counterpart on the compute side is
  !  compute_start_async_io/compute_shutdown_async_io
  !
#ifndef NOMPI
  SUBROUTINE async_io_wait_for_start(action, jstep)
    INTEGER, INTENT(OUT)          :: action ! pass on what to do
    INTEGER, INTENT(OUT)          :: jstep
    ! local variables
    INTEGER :: msg(2)
    TYPE(t_par_output_event), POINTER :: ev

    ! Set output parameters to default values
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

    SELECT CASE(msg(1))
    CASE(msg_io_start)
      jstep = msg(2)
    CASE(msg_io_meteogram_flush)
      jstep = msg(2)
    CASE(msg_io_shutdown)
    CASE DEFAULT
      ! Anything else is an error
      CALL finish(modname, 'I/O PE: Got illegal I/O tag')
    END SELECT
    action = msg(1)
  END SUBROUTINE async_io_wait_for_start
#endif

  !-------------------------------------------------------------------------------------------------
  ! ... called on compute procs:

  !-------------------------------------------------------------------------------------------------
  !> compute_wait_for_async_io: Wait for a message that the I/O is ready
  !  The counterpart on the I/O side is io_send_handshake
  !
#ifndef NOMPI
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
#endif

  !-------------------------------------------------------------------------------------------------
  !> compute_start_async_io: Send a message to I/O PEs that they should start I/O
  !  The counterpart on the I/O side is async_io_wait_for_start
  !
#ifndef NOMPI
  SUBROUTINE compute_start_async_io(jstep, output_pe_list, noutput_pe_list)
    INTEGER, INTENT(IN)          :: jstep
    INTEGER, INTENT(IN)          :: output_pe_list(:), noutput_pe_list
    ! local variables
    INTEGER :: msg(2)
    INTEGER  :: i

    IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": compute_start_async_io, jstep = ",jstep
    CALL p_barrier(comm=p_comm_work) ! make sure all are here
    msg(1) = msg_io_start
    msg(2) = jstep

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
#endif

  !-------------------------------------------------------------------------------------------------
  !> compute_shutdown_async_io: Send a message to I/O PEs that they should shut down
  !  The counterpart on the I/O side is async_io_wait_for_start
  !
#ifndef NOMPI
  SUBROUTINE compute_shutdown_async_io
    INTEGER :: msg(2)
    INTEGER  :: pe, i, ierror

    IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": compute_shutdown_async_io."
    CALL p_barrier(comm=p_comm_work) ! make sure all are here
    IF (ldebug)  WRITE (0,*) "pe ", p_pe, ": compute_shutdown_async_io barrier done."
    msg(1) = msg_io_shutdown
    msg(2) = 0
    ! tell all I/O PEs about the shutdown
    IF(p_pe_work==0) THEN
      DO pe = p_io_pe0, (p_io_pe0+num_io_procs-1)
        CALL p_send(msg, pe, 0)
      END DO
    END IF

    IF(use_async_name_list_io) THEN
      DO i = 1, SIZE(output_file)
        CALL mpi_win_free(output_file(i)%mem_win%mpi_win, ierror)
          IF (use_dp_mpi2io) THEN
            CALL mpi_free_mem(output_file(i)%mem_win%mem_ptr_dp, ierror)
          ELSE
            CALL mpi_free_mem(output_file(i)%mem_win%mem_ptr_sp, ierror)
          END IF
        CALL mpi_win_free(output_file(i)%mem_win%mpi_win_metainfo, ierror)
        CALL mpi_free_mem(output_file(i)%mem_win%mem_ptr_metainfo_pe0, ierror)
      END DO
    END IF
  END SUBROUTINE compute_shutdown_async_io
#endif

  !-------------------------------------------------------------------------------------------------
#endif

END MODULE mo_name_list_output
