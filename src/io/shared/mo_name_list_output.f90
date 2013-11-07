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
!!
!!   #define USE_CRAY_POINTER
!!
MODULE mo_name_list_output

#ifndef USE_CRAY_POINTER
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_intptr_t, c_f_pointer
#endif
! USE_CRAY_POINTER

  USE mo_kind,                      ONLY: wp, i8, dp, sp
  USE mo_impl_constants,            ONLY: zml_soil, max_dom, SUCCESS, MAX_TIME_LEVELS
  USE mo_grid_config,               ONLY: n_dom
  USE mo_cdi_constants              ! We need all
  USE mo_io_units,                  ONLY: filename_max, nnml, nnml_output, find_next_free_unit
  USE mo_io_config,                 ONLY: lkeep_in_sync
  USE mo_io_util,                   ONLY: get_file_extension
  USE mo_exception,                 ONLY: finish, message, message_text
  USE mo_var_metadata,              ONLY: t_var_metadata, POST_OP_NONE
  USE mo_var_list_element,          ONLY: level_type_ml, level_type_pl, level_type_hl,              &
    &                                     level_type_il, lev_type_str
  ! MPI Communication routines
  USE mo_mpi,                       ONLY: p_send, p_recv, p_bcast, p_barrier, p_stop,               &
    &                                     get_my_mpi_work_id, p_max,                                &
    &                                     get_my_mpi_work_communicator, p_mpi_wtime,                &
    &                                     p_irecv, p_wait, p_test, p_isend
  ! MPI Communicators
  USE mo_mpi,                       ONLY: p_comm_work
  ! MPI Data types
  USE mo_mpi,                       ONLY: p_int, p_int_i8, &
    &                                     p_real_dp, p_real_sp
  ! MPI Process type intrinsics
  USE mo_mpi,                       ONLY: my_process_is_stdio, my_process_is_mpi_test,              &
                                          my_process_is_mpi_workroot, my_process_is_mpi_seq,        &
                                          my_process_is_io, my_process_is_mpi_ioroot
  ! MPI Process IDs
  USE mo_mpi,                       ONLY: process_mpi_all_test_id, process_mpi_all_workroot_id,     &
                                          process_mpi_stdio_id
  ! MPI Process group sizes
  USE mo_mpi,                       ONLY: process_mpi_io_size, num_work_procs, p_n_work
  ! Processor numbers
  USE mo_mpi,                       ONLY: p_pe, p_pe_work, p_work_pe0, p_io_pe0

  USE mo_model_domain,              ONLY: t_patch, p_patch, p_phys_patch
  USE mo_parallel_config,           ONLY: nproma, p_test_run, use_dp_mpi2io

  USE mo_run_config,                ONLY: num_lev, num_levp1, dtime, ldump_states, ldump_dd,        &
    &                                     msg_level, output_mode, ltestcase
  USE mo_datetime,                  ONLY: t_datetime, cly360day_to_date
  USE mo_master_nml,                ONLY: model_base_dir
  USE mo_util_string,               ONLY: t_keyword_list, associate_keyword,                        &
    &                                     with_keywords, MAX_STRING_LEN,                            &
    &                                     tolower, int2string
  USE mo_communication,             ONLY: exchange_data, t_comm_pattern, idx_no, blk_no
  USE mo_math_constants,            ONLY: pi, pi_180
  USE mo_name_list_output_config,   ONLY: use_async_name_list_io, first_output_name_list
  USE mo_name_list_output_types,    ONLY:  l_output_phys_patch, t_output_name_list,                 &
  &                                        t_output_file, t_var_desc,                               &
  &                                        t_patch_info, t_reorder_info,                            &
  &                                        t_grid_info,                                             &
  &                                        REMAP_NONE, REMAP_REGULAR_LATLON,                        &
  &                                        msg_io_start, msg_io_done, msg_io_shutdown,              &
  &                                        IRLON, IRLAT, ILATLON, ICELL, IEDGE, IVERT,              &
  &                                        all_events
  USE mo_name_list_output_init,     ONLY:  init_name_list_output, setup_output_vlist,               &
    &                                      mem_ptr_sp, mem_ptr_dp, mpi_win,                         &
    &                                      varnames_dict, out_varnames_dict,                        &
    &                                      output_file, patch_info, lonlat_info,                    &
    &                                      i_sample
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_write_output, ltimer
  USE mo_dictionary,                ONLY: t_dictionary, dict_finalize
  USE mo_fortran_tools,             ONLY: assign_if_present
  ! post-ops
  USE mo_post_op,                   ONLY: perform_post_op
  USE mtime,                        ONLY: datetime, newDatetime, deallocateDatetime,                &
    &                                     PROLEPTIC_GREGORIAN, setCalendar
  USE mo_output_event_types,        ONLY: t_sim_step_info, t_par_output_event
  USE mo_output_event_handler,      ONLY: is_output_step, check_open_file, check_close_file,        &
    &                                     pass_output_step, get_current_filename,                   &
    &                                     get_current_date, trigger_output_step_irecv,              &
    &                                     is_output_step_complete, is_output_event_finished,        &
    &                                     check_write_readyfile, wait_for_final_irecvs

#ifdef __ICON_ATMO__
  USE mo_dynamics_config,           ONLY: nnow, nnow_rcf

! tool dependencies, maybe restructure
  USE mo_meteogram_output,          ONLY: meteogram_init, meteogram_finalize, meteogram_flush_file
  USE mo_meteogram_config,          ONLY: meteogram_output_config
  USE mo_lonlat_grid,               ONLY: t_lon_lat_grid, compute_lonlat_specs,                     &
    &                                     rotate_latlon_grid
  USE mo_intp_data_strc,            ONLY: lonlat_grid_list
#endif
! __ICON_ATMO__

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


CONTAINS

  !------------------------------------------------------------------------------------------------
  !> Writes the grid information in output file, GRIB2 format.
  !
#ifdef __ICON_ATMO__
  SUBROUTINE write_grid_info_grb2(of)
    TYPE (t_output_file), INTENT(INOUT) :: of
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::write_grid_info_grb2"

    TYPE t_grid_info_ptr
      TYPE (t_grid_info), POINTER :: ptr
    END TYPE t_grid_info_ptr

    INTEGER                        :: errstat, idom, igrid, idx(3), isize(3), idom_log
    TYPE (t_lon_lat_grid), POINTER :: grid
    REAL(wp), ALLOCATABLE          :: rotated_pts(:,:,:), r_out_dp(:,:), r_out_dp_1D(:)
    TYPE(t_grid_info_ptr)          :: gptr(3)

    ! skip this on test PE...
    IF (my_process_is_mpi_test()) RETURN

    SELECT CASE(of%name_list%remap)
    CASE (REMAP_NONE)
      idom     = of%phys_patch_id
      idom_log = patch_info(idom)%log_patch_id
      idx(:)    = (/ ICELL, IEDGE, IVERT /)
      isize(:)  = (/ patch_info(idom_log)%cells%n_glb, &
        &            patch_info(idom_log)%edges%n_glb, &
        &            patch_info(idom_log)%verts%n_glb /)
      gptr(1)%ptr => patch_info(idom_log)%grid_c
      gptr(2)%ptr => patch_info(idom_log)%grid_e
      gptr(3)%ptr => patch_info(idom_log)%grid_v
      DO igrid=1,3
        ! allocate data buffer:
        ALLOCATE(r_out_dp_1D(isize(igrid)), stat=errstat)
        IF (errstat /= SUCCESS) CALL finish(routine, 'ALLOCATE failed!')
        ! write RLON, RLAT
        r_out_dp_1D(:) = gptr(igrid)%ptr%lon(1:isize(igrid)) / pi_180
        CALL streamWriteVar(of%cdiFileID, of%cdi_grb2(idx(igrid),IRLON), r_out_dp_1D, 0)
        r_out_dp_1D(:) = gptr(igrid)%ptr%lat(1:isize(igrid)) / pi_180
        CALL streamWriteVar(of%cdiFileID, of%cdi_grb2(idx(igrid),IRLAT), r_out_dp_1D, 0)
        ! clean up
        DEALLOCATE(r_out_dp_1D, stat=errstat)
        IF (errstat /= SUCCESS) CALL finish(routine, 'DEALLOCATE failed!')
      END DO

    CASE (REMAP_REGULAR_LATLON)
      ! allocate data buffer:
      grid => lonlat_grid_list(of%name_list%lonlat_id)%grid
      ! compute some entries of lon-lat grid specification:
      CALL compute_lonlat_specs(grid)
      ALLOCATE(rotated_pts(grid%lon_dim, grid%lat_dim, 2), &
        &      r_out_dp(grid%lon_dim,grid%lat_dim), stat=errstat)
      IF (errstat /= SUCCESS) CALL finish(routine, 'ALLOCATE failed!')
      ! compute grid points of rotated lon/lat grid
      CALL rotate_latlon_grid(grid, rotated_pts)
      ! write RLON, RLAT
      r_out_dp(:,:) = rotated_pts(:,:,1) / pi_180
      CALL streamWriteVar(of%cdiFileID, of%cdi_grb2(ILATLON,IRLON), r_out_dp, 0)
      r_out_dp(:,:) = rotated_pts(:,:,2) / pi_180
      CALL streamWriteVar(of%cdiFileID, of%cdi_grb2(ILATLON,IRLAT), r_out_dp, 0)
      ! clean up
      DEALLOCATE(rotated_pts, r_out_dp, stat=errstat)
      IF (errstat /= SUCCESS) CALL finish(routine, 'DEALLOCATE failed!')

    CASE DEFAULT
      CALL finish(routine, "Unsupported grid type.")
    END SELECT

  END SUBROUTINE write_grid_info_grb2
#endif

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
    TYPE (t_datetime)                 :: rel_fct_time
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
    INTEGER :: i, jg

#ifndef NOMPI
#ifdef __ICON_ATMO__
    IF (use_async_name_list_io    .AND.  &
      & .NOT. my_process_is_io()  .AND.  &
      & .NOT. my_process_is_mpi_test()) THEN
      !-- compute PEs (senders):

#ifdef __ICON_ATMO__
      ! write recent samples of meteogram output
      DO jg = 1, n_dom
        IF (meteogram_output_config(jg)%lenabled) THEN
          CALL meteogram_flush_file(jg)
        END IF
      END DO
#endif
! __ICON_ATMO__

      CALL compute_wait_for_async_io()
      CALL compute_shutdown_async_io()

    ELSE
#endif
! __ICON_ATMO__
#endif
! NOMPI
      !-- asynchronous I/O PEs (receiver):
      DO i = 1, SIZE(output_file)
        IF (output_file(i)%cdiFileID >= 0) THEN
          CALL close_output_file(output_file(i))
          CALL destroy_output_vlist(output_file(i))
        END IF
      ENDDO
#ifndef NOMPI
#ifdef __ICON_ATMO__
    ENDIF
#endif
! __ICON_ATMO__
#endif
! NOMPI

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
      &                 ((.NOT. use_async_name_list_io) .AND. my_process_is_stdio())

    IF(of%cdiFileID /= CDI_UNDEFID) THEN
#ifdef __ICON_ATMO__
      IF (of%name_list%output_grid .AND. &
        & is_output_process        .AND. &
        & (of%name_list%filetype == FILETYPE_GRB2)) THEN
        CALL write_grid_info_grb2(of)
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
  SUBROUTINE write_name_list_output(jstep)
    INTEGER, INTENT(IN) :: jstep         !< model step
    ! local variables
    CHARACTER(LEN=*), PARAMETER  :: routine = modname//"::write_name_list_output"
    INTEGER                           :: i, idate, itime, iret, jg
    TYPE(t_output_name_list), POINTER :: p_onl
    TYPE(datetime),           POINTER :: mtime_datetime
    CHARACTER(LEN=filename_max+100)   :: text
    TYPE(t_par_output_event), POINTER :: ev

    IF (ltimer) CALL timer_start(timer_write_output)
#ifndef NOMPI
#ifdef __ICON_ATMO__
    IF(use_async_name_list_io) THEN
      IF(.NOT.my_process_is_io().AND..NOT.my_process_is_mpi_test()) THEN
        ! write recent samples of meteogram output
        DO jg = 1, n_dom
          IF (meteogram_output_config(jg)%lenabled) THEN
            CALL meteogram_flush_file(jg)
          END IF
        END DO
       
        ! If asynchronous I/O is enabled, the compute PEs have to make
        ! sure that the I/O PEs are ready with the last output step
        ! before writing data into the I/O memory window.  This
        ! routine (write_name_list_output) is also called from the I/O
        ! PEs, but in this case the calling routine cares about the
        ! flow control.
        CALL compute_wait_for_async_io()
      END IF
    ENDIF
#endif
! __ICON_ATMO__
#endif
! NOMPI

    ! Go over all output files
    OUTFILE_LOOP : DO i=1,SIZE(output_file)

      ! Skip this output file if it is not due for output!
      IF (.NOT. is_output_step(output_file(i)%out_event, jstep))  CYCLE OUTFILE_LOOP

      ! -------------------------------------------------
      ! Check if files have to be (re)opened
      ! -------------------------------------------------
      IF (check_open_file(output_file(i)%out_event) .AND.  &
        & (output_file(i)%io_proc_id == p_pe)) THEN 
        CALL setup_output_vlist(output_file(i))
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
      & .NOT.my_process_is_mpi_test()) CALL compute_start_async_io(jstep)
#endif

    ! Handle incoming "output step completed" messages: After all
    ! participating I/O PE's have acknowledged the completion of their
    ! write processes, we trigger a "ready file" on the first I/O PE.
    IF ((      use_async_name_list_io .AND. my_process_is_mpi_ioroot()) .OR.  &
      & (.NOT. use_async_name_list_io .AND. my_process_is_stdio())) THEN
      ev => all_events
      HANDLE_COMPLETE_STEPS : DO
        IF (.NOT. ASSOCIATED(ev)) EXIT HANDLE_COMPLETE_STEPS
        IF (.NOT. is_output_step_complete(ev) .OR.  &
          & is_output_event_finished(ev)) THEN 
          ev => ev%next
          CYCLE HANDLE_COMPLETE_STEPS
        END IF         
        !--- write ready file
        IF (check_write_readyfile(ev%output_event)) THEN
          CALL write_ready_file(TRIM(ev%output_event%event_data%name)//"_"//TRIM(get_current_date(ev))//".rdy")
        END IF
        ! launch a non-blocking request to all participating PEs to
        ! acknowledge the completion of the next output event
        CALL trigger_output_step_irecv(ev)
      END DO HANDLE_COMPLETE_STEPS
    END IF
    IF (ltimer) CALL timer_stop(timer_write_output)

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
  SUBROUTINE write_ready_file(rdy_filename)
    CHARACTER(LEN=*), INTENT(IN) :: rdy_filename
    ! local variables
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//"::write_ready_file"
    INTEGER                       :: iunit

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
  !> Write a output name list
  !
  SUBROUTINE write_name_list(of, l_first_write)

#ifndef NOMPI
#ifdef  __SUNPRO_F95
  INCLUDE "mpif.h"
#else
    USE mpi, ONLY: MPI_LOCK_EXCLUSIVE, MPI_MODE_NOCHECK
#endif
! __SUNPRO_F95
#endif
! NOMPI

    TYPE (t_output_file), INTENT(INOUT), TARGET :: of
    LOGICAL,              INTENT(IN)            :: l_first_write
    ! local variables:
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//"::write_name_list"
    REAL(wp),         PARAMETER    :: SYNC_ERROR_PRINT_TOL = 1e-13_wp
    INTEGER,          PARAMETER    :: iUNKNOWN = 0
    INTEGER,          PARAMETER    :: iINTEGER = 1
    INTEGER,          PARAMETER    :: iREAL    = 2

    INTEGER                        :: tl, i_dom, i_log_dom, i, iv, jk, n_points, &
      &                               nlevs, nblks, nindex, mpierr, lonlat_id,   &
      &                               ierrstat, idata_type
    INTEGER(i8)                    :: ioff
    TYPE (t_var_metadata), POINTER :: info
    TYPE(t_reorder_info),  POINTER :: p_ri
    REAL(wp),          ALLOCATABLE :: r_ptr(:,:,:)
    INTEGER,           ALLOCATABLE :: i_ptr(:,:,:)
    REAL(wp),          ALLOCATABLE :: r_tmp(:,:,:), r_out_recv(:,:)
    INTEGER,           ALLOCATABLE :: i_tmp(:,:,:)
    REAL(sp),          ALLOCATABLE :: r_out_sp(:,:)
    REAL(dp),          ALLOCATABLE :: r_out_dp(:,:)
    TYPE(t_comm_pattern),  POINTER :: p_pat
    LOGICAL                        :: l_error

    ! Offset in memory window for async I/O
    ioff = of%my_mem_win_off

    i_dom = of%phys_patch_id
    i_log_dom = of%log_patch_id

    tl = 0 ! to prevent warning

#ifndef NOMPI
    ! In case of async IO: Lock own window before writing to it
    IF(use_async_name_list_io .AND. .NOT.my_process_is_mpi_test()) &
      CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, p_pe_work, MPI_MODE_NOCHECK, mpi_win, mpierr)
#endif
! NOMPI

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
      IF(info%used_dimensions(1) /= nproma) &
        CALL finish(routine,'1st dim is not nproma: '//TRIM(info%name))

      idata_type = iUNKNOWN

      ! For time level dependent elements: set time level and check if
      ! time level is present:
        ! set a default time level (which is not used anyway, but must
        ! be a valid array subscript):
      tl = 1
#ifdef __ICON_ATMO__
      IF (.NOT. ASSOCIATED(of%var_desc(iv)%r_ptr)  .AND.    &
        & .NOT. ASSOCIATED(of%var_desc(iv)%i_ptr)) THEN
        SELECT CASE (info%tlev_source)
          CASE(0); tl = nnow(i_log_dom)
          CASE(1); tl = nnow_rcf(i_log_dom)
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
! __ICON_ATMO__

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

      IF (of%var_desc(iv)%info%post_op%ipost_op_type /= POST_OP_NONE) THEN
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
        IF(l_output_phys_patch) THEN
          p_pat => p_phys_patch(i_dom)%comm_pat_gather_c
        ELSE
          p_pat => p_patch(i_dom)%comm_pat_gather_c
        ENDIF
      CASE (GRID_UNSTRUCTURED_EDGE)
        p_ri => patch_info(i_dom)%edges
        IF(l_output_phys_patch) THEN
          p_pat => p_phys_patch(i_dom)%comm_pat_gather_e
        ELSE
          p_pat => p_patch(i_dom)%comm_pat_gather_e
        ENDIF
      CASE (GRID_UNSTRUCTURED_VERT)
        p_ri => patch_info(i_dom)%verts
        IF(l_output_phys_patch) THEN
          p_pat => p_phys_patch(i_dom)%comm_pat_gather_v
        ELSE
          p_pat => p_patch(i_dom)%comm_pat_gather_v
        ENDIF
#ifdef __ICON_ATMO__
      CASE (GRID_REGULAR_LONLAT)
        lonlat_id = info%hor_interp%lonlat_id
        p_ri  => lonlat_info(lonlat_id, i_log_dom)
        p_pat => lonlat_grid_list(lonlat_id)%p_pat(i_log_dom)
#endif
! __ICON_ATMO__
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
          IF(my_process_is_mpi_workroot()) THEN
            ALLOCATE(r_tmp(nproma,nlevs,nblks))
          ELSE
            ! Dimensions 1 and 2 of r_tmp must always be nproma and nlevs,
            ! otherwise exchange_data doesn't work!
            ALLOCATE(r_tmp(nproma,nlevs,1))
          ENDIF
          r_tmp(:,:,:) = 0._wp
          ! Gather data on root
          IF(my_process_is_mpi_seq()) THEN
            DO jk = 1, nlevs
              DO i = 1, p_ri%n_own
                r_tmp(p_ri%own_dst_idx(i),jk,p_ri%own_dst_blk(i)) = r_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i))
              ENDDO
            ENDDO
          ELSE
            CALL exchange_data(p_pat, RECV=r_tmp, SEND=r_ptr)
          ENDIF
        END IF
        IF (idata_type == iINTEGER) THEN
          IF(my_process_is_mpi_workroot()) THEN
            ALLOCATE(i_tmp(nproma,nlevs,nblks))
          ELSE
            ! Dimensions 1 and 2 of r_tmp must always be nproma and nlevs,
            ! otherwise exchange_data doesn't work!
            ALLOCATE(i_tmp(nproma,nlevs,1))
          ENDIF
          i_tmp(:,:,:) = 0
          ! Gather data on root
          IF(my_process_is_mpi_seq()) THEN
            DO jk = 1, nlevs
              DO i = 1, p_ri%n_own
                i_tmp(p_ri%own_dst_idx(i),jk,p_ri%own_dst_blk(i)) = i_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i))
              ENDDO
            ENDDO
          ELSE
            CALL exchange_data(p_pat, RECV=i_tmp, SEND=i_ptr)
          ENDIF
        END IF

        IF(my_process_is_mpi_workroot()) THEN

          ! De-block the array
          IF (use_dp_mpi2io) THEN
            ALLOCATE(r_out_dp(n_points, nlevs), STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
            IF (idata_type == iREAL) THEN
              DO jk = 1, nlevs
                r_out_dp(:,jk) = REAL(RESHAPE(r_tmp(:,jk,:), (/ n_points /)), dp)
              ENDDO
            END IF
            IF (idata_type == iINTEGER) THEN
              DO jk = 1, nlevs
                r_out_dp(:,jk) = REAL(RESHAPE(i_tmp(:,jk,:), (/ n_points /)), dp)
              ENDDO
            END IF
          ELSE
            ALLOCATE(r_out_sp(n_points, nlevs), STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
            IF (idata_type == iREAL) THEN
              DO jk = 1, nlevs
                r_out_sp(:,jk) = REAL(RESHAPE(r_tmp(:,jk,:), (/ n_points /)), sp)
              ENDDO
            END IF
            IF (idata_type == iINTEGER) THEN
              DO jk = 1, nlevs
                r_out_sp(:,jk) = REAL(RESHAPE(i_tmp(:,jk,:), (/ n_points /)), sp)
              ENDDO
            END IF
          ENDIF

          ! ------------------
          ! case of a test run
          ! ------------------
          !
          ! compare results on worker PEs and test PE
          IF (p_test_run  .AND.  .NOT. use_async_name_list_io  .AND.  use_dp_mpi2io) THEN
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

        ! ----------
        ! write data
        ! ----------

        IF (my_process_is_stdio() .AND. .NOT. my_process_is_mpi_test()) THEN
          IF (use_dp_mpi2io) THEN
            CALL streamWriteVar(of%cdiFileID, info%cdiVarID, r_out_dp, 0)
          ELSE
            CALL streamWriteVarF(of%cdiFileID, info%cdiVarID, r_out_sp, 0)
          ENDIF
          IF (lkeep_in_sync) CALL streamSync(of%cdiFileID)
        ENDIF

        ! clean up
        IF (ALLOCATED(r_tmp))  DEALLOCATE(r_tmp)
        IF (ALLOCATED(i_tmp))  DEALLOCATE(i_tmp)

        IF (my_process_is_mpi_workroot()) THEN
          IF (use_dp_mpi2io) THEN
            DEALLOCATE(r_out_dp, STAT=ierrstat)
          ELSE
            DEALLOCATE(r_out_sp, STAT=ierrstat)
          END IF
          IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
        ENDIF

      ELSE

        ! ------------------------
        ! Asynchronous I/O is used
        ! ------------------------
        !
        ! just copy the OWN DATA points to the memory window

        DO jk = 1, nlevs
          IF (use_dp_mpi2io) THEN
            IF (idata_type == iREAL) THEN
              DO i = 1, p_ri%n_own
                mem_ptr_dp(ioff+INT(i,i8)) = REAL(r_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i)),dp)
              ENDDO
            END IF
            IF (idata_type == iINTEGER) THEN
              DO i = 1, p_ri%n_own
                mem_ptr_dp(ioff+INT(i,i8)) = REAL(i_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i)),dp)
              ENDDO
            END IF
          ELSE
            IF (idata_type == iREAL) THEN
              DO i = 1, p_ri%n_own
                mem_ptr_sp(ioff+INT(i,i8)) = REAL(r_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i)),sp)
              ENDDO
            END IF
            IF (idata_type == iINTEGER) THEN
              DO i = 1, p_ri%n_own
                mem_ptr_sp(ioff+INT(i,i8)) = REAL(i_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i)),sp)
              ENDDO
            END IF
          END IF
          ioff = ioff + INT(p_ri%n_own,i8)
        END DO

      END IF

      ! clean up
      IF (ALLOCATED(r_ptr)) DEALLOCATE(r_ptr)
      IF (ALLOCATED(i_ptr)) DEALLOCATE(i_ptr)

    ENDDO

#ifndef NOMPI
    ! In case of async IO: Done writing to memory window, unlock it
    IF(use_async_name_list_io .AND. .NOT.my_process_is_mpi_test()) &
      CALL MPI_Win_unlock(p_pe_work, mpi_win, mpierr)
#endif
! NOMPI

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

#ifdef __ICON_ATMO__
    LOGICAL             :: done
    INTEGER             :: jg, jstep
    TYPE(t_par_output_event), POINTER :: ev

    ! If ldump_states or ldump_dd is set, the compute PEs will exit after dumping,
    ! there is nothing to do at all for I/O PEs

    IF(ldump_states .OR. ldump_dd) THEN
      CALL p_stop
      STOP
    ENDIF

    ! setup of meteogram output
    DO jg =1,n_dom
      IF (meteogram_output_config(jg)%lenabled) THEN
        CALL meteogram_init(meteogram_output_config(jg), jg)
      END IF
    END DO

    ! Initialize name list output, this is a collective call for all PEs
    CALL init_name_list_output(sim_step_info, opt_isample=isample)

    ! Tell the compute PEs that we are ready to work
    CALL async_io_send_handshake

    ! write recent samples of meteogram output
    DO jg = 1, n_dom
      IF (meteogram_output_config(jg)%lenabled) THEN
        CALL meteogram_flush_file(jg)
      END IF
    END DO

    ! Enter I/O loop
    DO
      ! Wait for a message from the compute PEs to start
      CALL async_io_wait_for_start(done, jstep)
      IF(done) EXIT ! leave loop, we are done

      ! perform I/O
      CALL write_name_list_output(jstep)

      ! Inform compute PEs that we are done
      CALL async_io_send_handshake

      ! write recent samples of meteogram output
      DO jg = 1, n_dom
        IF (meteogram_output_config(jg)%lenabled) THEN
          CALL meteogram_flush_file(jg)
        END IF
      END DO
    ENDDO

    ! Handle final pending "output step completed" messages: After all
    ! participating I/O PE's have acknowledged the completion of their
    ! write processes, we trigger a "ready file" on the first I/O PE.
    IF ((      use_async_name_list_io .AND. my_process_is_mpi_ioroot()) .OR.  &
      & (.NOT. use_async_name_list_io .AND. my_process_is_stdio())) THEN
      CALL wait_for_final_irecvs(all_events)
      ev => all_events
      HANDLE_COMPLETE_STEPS : DO
        IF (.NOT. ASSOCIATED(ev)) EXIT HANDLE_COMPLETE_STEPS

        !--- write ready file
        IF (check_write_readyfile(ev%output_event)) THEN
          CALL write_ready_file(TRIM(ev%output_event%event_data%name)//"_"//TRIM(get_current_date(ev))//".rdy")
        END IF
        ev => ev%next
      END DO HANDLE_COMPLETE_STEPS
    END IF

    ! Finalization sequence:
    CALL close_name_list_output

    ! finalize meteogram output
    DO jg = 1, n_dom
      IF (meteogram_output_config(jg)%lenabled) THEN
        CALL meteogram_finalize(jg)
      END IF
    END DO

    DO jg = 1, max_dom
      DEALLOCATE(meteogram_output_config(jg)%station_list)
    END DO
    ! Shut down MPI
    !
    CALL p_stop

    STOP
#endif
! __ICON_ATMO__

  END SUBROUTINE name_list_io_main_proc
  !------------------------------------------------------------------------------------------------


  !------------------------------------------------------------------------------------------------
  !> Output routine on the IO PE
  !
  SUBROUTINE io_proc_write_name_list(of, l_first_write)

#ifdef __SUNPRO_F95
    INCLUDE "mpif.h"
#else
    USE mpi, ONLY: MPI_ADDRESS_KIND, MPI_LOCK_SHARED, MPI_MODE_NOCHECK
#endif
! __SUNPRO_F95

    TYPE (t_output_file), TARGET, INTENT(IN) :: of
    LOGICAL                     , INTENT(IN) :: l_first_write

    CHARACTER(LEN=*), PARAMETER    :: routine = modname//"::io_proc_write_name_list"

    INTEGER                        :: nval, nlev_max, iv, jk, i, nlevs, mpierr, nv_off, np, i_dom, &
      &                               lonlat_id, i_log_dom, ierrstat
    INTEGER(KIND=MPI_ADDRESS_KIND) :: ioff(0:num_work_procs-1)
    INTEGER                        :: voff(0:num_work_procs-1)
    REAL(sp), ALLOCATABLE          :: var1_sp(:), var2_sp(:), var3_sp(:,:)
    REAL(dp), ALLOCATABLE          :: var1_dp(:), var2_dp(:), var3_dp(:,:)
    TYPE (t_var_metadata), POINTER :: info
    TYPE(t_reorder_info) , POINTER :: p_ri
    LOGICAL                        :: have_GRIB

    !-- for timing
    CHARACTER(len=10)              :: ctime
    REAL(dp)                       :: t_get, t_write, t_copy, t_intp, t_0, mb_get, mb_wr

  !------------------------------------------------------------------------------------------------
#if defined (__SX__) && !defined (NOMPI)
! It may be necessary that var1 is in global memory on NEC
! (Note: this is only allowed when we compile with MPI.)
!CDIR GM_ARRAY(var1)
#endif
! __SX__

    CALL date_and_time(TIME=ctime)
    WRITE (0, '(a,i0,a)') '#################### I/O PE ',p_pe,' starting I/O at '//ctime

    t_get   = 0.d0
    t_write = 0.d0
    t_copy  = 0.d0
    t_intp  = 0.d0
    mb_get  = 0.d0
    mb_wr   = 0.d0


    ! Get maximum number of data points in a slice and allocate tmp variables

    i_dom = of%phys_patch_id
    nval = MAX(patch_info(i_dom)%cells%n_glb, &
               patch_info(i_dom)%edges%n_glb, &
               patch_info(i_dom)%verts%n_glb)
#ifdef __ICON_ATMO__
    ! take also the lon-lat grids into account
    DO iv = 1, of%num_vars
      info => of%var_desc(iv)%info
      IF (info%hgrid == GRID_REGULAR_LONLAT) THEN
        lonlat_id = info%hor_interp%lonlat_id
        i_log_dom = of%log_patch_id
        p_ri  => lonlat_info(lonlat_id, i_log_dom)
        nval = MAX(nval, p_ri%n_glb)
      END IF
    END DO
#endif
! __ICON_ATMO__

    nlev_max = 1
    DO iv = 1, of%num_vars
      info => of%var_desc(iv)%info
      IF(info%ndims == 3) nlev_max = MAX(nlev_max, info%used_dimensions(2))
    ENDDO

    IF (use_dp_mpi2io) THEN
      ALLOCATE(var1_dp(nval*nlev_max), var2_dp(-1:nval), STAT=ierrstat)
    ELSE
      ALLOCATE(var1_sp(nval*nlev_max), var2_sp(-1:nval), STAT=ierrstat)
    ENDIF
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ioff(:) = of%mem_win_off(:)


    ! Go over all name list variables for this output file

    DO iv = 1, of%num_vars

      info => of%var_desc(iv)%info

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

#ifdef __ICON_ATMO__
        CASE (GRID_REGULAR_LONLAT)
          lonlat_id = info%hor_interp%lonlat_id
          i_log_dom = of%log_patch_id
          p_ri  => lonlat_info(lonlat_id, i_log_dom)
#endif
! __ICON_ATMO__

        CASE DEFAULT
          CALL finish(routine,'unknown grid type')
      END SELECT

      ! Retrieve part of variable from every worker PE using MPI_Get

      nv_off = 0
      DO np = 0, num_work_procs-1

        IF(p_ri%pe_own(np) == 0) CYCLE

        nval = p_ri%pe_own(np)*nlevs ! Number of words to transfer

        t_0 = p_mpi_wtime()
        CALL MPI_Win_lock(MPI_LOCK_SHARED, np, MPI_MODE_NOCHECK, mpi_win, mpierr)

        IF (use_dp_mpi2io) THEN
          CALL MPI_Get(var1_dp(nv_off+1), nval, p_real_dp, np, ioff(np), &
            &          nval, p_real_dp, mpi_win, mpierr)
        ELSE
          CALL MPI_Get(var1_sp(nv_off+1), nval, p_real_sp, np, ioff(np), &
            &          nval, p_real_sp, mpi_win, mpierr)
        ENDIF

        CALL MPI_Win_unlock(np, mpi_win, mpierr)
        t_get  = t_get  + p_mpi_wtime() - t_0
        mb_get = mb_get + nval

        ! Update the offset in var1
        nv_off = nv_off + nval

        ! Update the offset in the memory window on compute PEs
        ioff(np) = ioff(np) + INT(nval,i8)

      ENDDO

      ! compute the total offset for each PE
      nv_off = 0
      DO np = 0, num_work_procs-1
        voff(np) = nv_off
        nval     = p_ri%pe_own(np)*nlevs
        nv_off   = nv_off + nval
      END DO

      ! var1 is stored in the order in which the variable was stored on compute PEs,
      ! get it back into the global storage order

      t_0 = p_mpi_wtime()
      have_GRIB = of%output_type == FILETYPE_GRB .OR. of%output_type == FILETYPE_GRB2
      IF (use_dp_mpi2io .OR. have_GRIB) THEN
        ALLOCATE(var3_dp(p_ri%n_glb,nlevs), STAT=ierrstat) ! Must be allocated to exact size
      ELSE
        ALLOCATE(var3_sp(p_ri%n_glb,nlevs), STAT=ierrstat) ! Must be allocated to exact size
      ENDIF
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')


      LevelLoop: DO jk = 1, nlevs

        nv_off = 0
        IF (use_dp_mpi2io) THEN

          DO np = 0, num_work_procs-1
            var2_dp(-1) = 0._dp ! special value for lon-lat areas overlapping local patches
            var2_dp(nv_off+1:nv_off+p_ri%pe_own(np)) = var1_dp(voff(np)+1:voff(np)+p_ri%pe_own(np))
            nv_off   = nv_off+p_ri%pe_own(np)
            voff(np) = voff(np)+p_ri%pe_own(np)
          ENDDO

          var3_dp(1:p_ri%n_glb,jk) = var2_dp(p_ri%reorder_index(1:p_ri%n_glb))

        ELSE

          DO np = 0, num_work_procs-1
            var2_sp(-1) = 0._sp ! special value for lon-lat areas overlapping local patches
            var2_sp(nv_off+1:nv_off+p_ri%pe_own(np)) = var1_sp(voff(np)+1:voff(np)+p_ri%pe_own(np))
            nv_off   = nv_off+p_ri%pe_own(np)
            voff(np) = voff(np)+p_ri%pe_own(np)
          ENDDO

          IF (have_GRIB) THEN
            ! ECMWF GRIB-API/CDI has only a double precision interface at the date of coding this
            var3_dp(1:p_ri%n_glb,jk) = var2_sp(p_ri%reorder_index(1:p_ri%n_glb))
          ELSE
            var3_sp(1:p_ri%n_glb,jk) = var2_sp(p_ri%reorder_index(1:p_ri%n_glb))
          END IF

        ENDIF

      ENDDO LevelLoop
      t_copy = t_copy + p_mpi_wtime() - t_0


      t_0 = p_mpi_wtime()
      IF (use_dp_mpi2io .OR. have_GRIB) THEN
        CALL streamWriteVar(of%cdiFileID, info%cdiVarID, var3_dp, 0)
        mb_wr = mb_wr + REAL(SIZE(var3_dp), wp)
        t_write = t_write + p_mpi_wtime() - t_0
        DEALLOCATE(var3_dp, STAT=ierrstat)
      ELSE
        CALL streamWriteVarF(of%cdiFileID, info%cdiVarID, var3_sp, 0)
        mb_wr = mb_wr + REAL(SIZE(var3_sp),wp)
        t_write = t_write + p_mpi_wtime() - t_0
        DEALLOCATE(var3_sp, STAT=ierrstat)
      ENDIF
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    ENDDO ! Loop over output variables

    IF (use_dp_mpi2io) THEN
      DEALLOCATE(var1_dp, var2_dp, STAT=ierrstat)
    ELSE
      DEALLOCATE(var1_sp, var2_sp, STAT=ierrstat)
    ENDIF
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')


    !
    !-- timing report
    !
    CALL date_and_time(TIME=ctime)
    WRITE (0, '(a,i0,a)') '#################### I/O PE ',p_pe,' done at '//ctime
    ! Convert mb_get/mb_wr to MB
    mb_get = mb_get*8*1.d-6
    mb_wr  = mb_wr*4*1.d-6 ! 4 byte since dp output is implicitly converted to sp
    ! writing this message causes a runtime error on the NEC because formatted output to stdio/stderr is limited to 132 chars
#ifndef __SX__
    IF (msg_level >= 12) THEN
      WRITE (message_text,'(10(a,f10.3))') &
           & ' Got ',mb_get,' MB, time get: ',t_get,' s [',mb_get/MAX(1.e-6_wp,t_get), &
           & ' MB/s], time write: ',t_write,' s [',mb_wr/MAX(1.e-6_wp,t_write),        &
           & ' MB/s], times copy+intp: ',t_copy+t_intp,' s'
      CALL message('',message_text)
    ENDIF
#endif
! __SX__

    ! Convert mb_get/mb_wr to MB
    IF (use_dp_mpi2io) THEN
      mb_get = mb_get*8*1.d-6
    ELSE
      mb_get = mb_get*4*1.d-6
    ENDIF
    mb_wr = mb_wr*4*1.d-6 ! always 4 byte since dp output is implicitly converted to sp
    ! PRINT *,' Got ',mb_get,' MB, time get: ',t_get,' s [',mb_get/t_get,&
    !   ' MB/s], time write: ',t_write,' s [',mb_wr/t_write, &
    !   ' MB/s], times copy+intp: ',t_copy+t_intp,' s'

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
  SUBROUTINE async_io_send_handshake
    ! local variables
    REAL(wp) :: msg
    TYPE(t_par_output_event), POINTER :: ev

    ! make sure all are done
    CALL p_barrier(comm=p_comm_work)

    ! Simply send a message from I/O PE 0 to compute PE 0
    ! 
    ! Note: We have to do this in a non-blocking fashion in order to
    !       receive "ready file" messages.
    IF(p_pe_work == 0) THEN
      ! launch non-blocking receive request:
      CALL p_wait()
      msg = REAL(msg_io_done, wp)
      CALL p_isend(msg, p_work_pe0, 0)
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
          IF (check_write_readyfile(ev%output_event)) THEN
            CALL write_ready_file(TRIM(ev%output_event%event_data%name)//"_"//TRIM(get_current_date(ev))//".rdy")
          END IF
          ! launch a non-blocking request to all participating PEs to
          ! acknowledge the completion of the next output event
          CALL trigger_output_step_irecv(ev)
        END DO HANDLE_COMPLETE_STEPS

        IF (p_test()) EXIT
      END DO
      CALL p_wait()
    END IF

  END SUBROUTINE async_io_send_handshake


  !-------------------------------------------------------------------------------------------------
  !> async_io_wait_for_start: Wait for a message from work PEs that we
  !  should start I/O or finish.  The counterpart on the compute side is
  !  compute_start_io/compute_shutdown_io
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
    ! Note: We have to do this in a non-blocking fashion in order to
    !       receive "ready file" messages.
    IF(p_pe_work == 0) THEN
      ! launch non-blocking receive request:
      CALL p_wait()
      CALL p_irecv(msg, p_work_pe0, 0)
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
          IF (check_write_readyfile(ev%output_event)) THEN
            CALL write_ready_file(TRIM(ev%output_event%event_data%name)//"_"//TRIM(get_current_date(ev))//".rdy")
          END IF
          ! launch a non-blocking request to all participating PEs to
          ! acknowledge the completion of the next output event
          CALL trigger_output_step_irecv(ev)
        END DO HANDLE_COMPLETE_STEPS

        IF (p_test()) EXIT
      END DO
      CALL p_wait()
    END IF

    CALL p_bcast(msg, 0, comm=p_comm_work)
    
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
  SUBROUTINE compute_wait_for_async_io
    REAL(wp) :: msg

    ! First compute PE receives message from I/O leader
    IF(p_pe_work==0) THEN
      CALL p_recv(msg, p_io_pe0, 0)
      ! Just for safety: Check if we got the correct tag
      IF(INT(msg) /= msg_io_done) CALL finish(modname, 'Compute PE: Got illegal I/O tag')
    ENDIF
    ! Wait in barrier until message is here
    CALL p_barrier(comm=p_comm_work)
  END SUBROUTINE compute_wait_for_async_io


  !-------------------------------------------------------------------------------------------------
  !> compute_start_async_io: Send a message to I/O PEs that they should start I/O
  !  The counterpart on the I/O side is io_wait_for_start_message
  !
  SUBROUTINE compute_start_async_io(jstep)
    INTEGER, INTENT(IN)          :: jstep
    ! local variables
    REAL(wp) :: msg(2)

    CALL p_barrier(comm=p_comm_work) ! make sure all are here
    msg(1) = REAL(msg_io_start,    wp)
    msg(2) = REAL(jstep,           wp)
    IF(p_pe_work==0) CALL p_send(msg, p_io_pe0, 0)
  END SUBROUTINE compute_start_async_io


  !-------------------------------------------------------------------------------------------------
  !> compute_shutdown_async_io: Send a message to I/O PEs that they should shut down
  !  The counterpart on the I/O side is io_wait_for_start_message
  !
  SUBROUTINE compute_shutdown_async_io
    REAL(wp) :: msg(2)

    CALL p_barrier(comm=p_comm_work) ! make sure all are here
    msg(1) = REAL(msg_io_shutdown, wp)
    msg(2:) = 0._wp
    IF(p_pe_work==0) CALL p_send(msg, p_io_pe0, 0)
  END SUBROUTINE compute_shutdown_async_io

  !-------------------------------------------------------------------------------------------------
#endif
! NOMPI

END MODULE mo_name_list_output
