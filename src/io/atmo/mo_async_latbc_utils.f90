  !>
  !! This module contains the asynchronous I/O routine for lateral boundary nudging
  !!
  !! @author M. Pondkule (DWD)
  !!
  !!
  !! @par Revision History
  !! Initial version by M. Pondkule, DWD (2014-01-27)
  !! Allow boundary data from the ICON output by S. Brdar, DWD (2013-07-19)
  !!
  !! @par Copyright and License
  !!
  !! This code is subject to the DWD and MPI-M-Software-License-Agreement in
  !! its most recent form.
  !! Please see the file LICENSE in the root of the source tree for this code.
  !! Where software is supplied by third parties, it is indicated in the
  !! headers of the routines.
  !!
  !!
  !----------------------------
#include "omp_definitions.inc"
  !----------------------------

  MODULE mo_async_latbc_utils

#ifndef NOMPI
    USE mpi
    USE mo_mpi,                 ONLY: my_process_is_mpi_test, &
         &                            my_process_is_pref, my_process_is_work,   &
         &                            my_process_is_io, p_comm_work
    ! Processor numbers
    USE mo_mpi,                 ONLY: p_pref_pe0, p_pe_work, p_work_pe0, num_work_procs
    ! MPI Communication routines
    USE mo_mpi,                 ONLY: p_isend, p_irecv, p_barrier, p_wait, &
         &                            p_send, p_recv
    USE mo_latbc_read_recv,     ONLY: prefetch_cdi_2d, prefetch_cdi_3d, compute_data_receive
#endif

    USE mo_async_latbc_types,   ONLY: t_patch_data, t_reorder_data, t_latbc_data
    USE mo_kind,                ONLY: wp, sp, i8
    USE mo_parallel_config,     ONLY: nproma
    USE mo_model_domain,        ONLY: t_patch
    USE mo_grid_config,         ONLY: nroot
    USE mo_exception,           ONLY: message, message_text, finish
    USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, MODE_COSMODE, MODE_DWDANA, MODE_ICONVREMAP, &
                                      MODE_IAU_OLD, MODE_IAU, SUCCESS
    USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
    USE mo_io_units,            ONLY: filename_max
    USE mo_nonhydro_types,      ONLY: t_nh_state
    USE mo_intp_data_strc,      ONLY: t_int_state
    USE mo_nh_vert_interp,      ONLY: vert_interp
    USE mo_physical_constants,  ONLY: cpd, rd, cvd_o_rd, p0ref, vtmpc1
    USE mo_util_phys,           ONLY: virtual_temp
    USE mo_util_string,         ONLY: int2string
    USE mo_nh_init_utils,       ONLY: interp_uv_2_vn, convert_thdvars
    USE mo_sync,                ONLY: sync_patch_array, sync_patch_array_mult, SYNC_E, SYNC_C
    USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
    USE mtime,                  ONLY: timedelta, newTimedelta, deallocateTimedelta, &
         &                            event, newEvent, datetime, newDatetime,      &
         &                            isCurrentEventActive, deallocateDatetime,    &
         &                            MAX_DATETIME_STR_LEN, MAX_EVENTNAME_STR_LEN, &
         &                            MAX_TIMEDELTA_STR_LEN, getPTStringFromMS,    &
         &                            OPERATOR(>=), OPERATOR(-), OPERATOR(>),      &
         &                            OPERATOR(/=), datetimeToString,              &
         &                            timedelta, newTimedelta, deallocateTimedelta,&
         &                            OPERATOR(+), deallocateEvent, OPERATOR(*)
    USE mo_time_config,         ONLY: time_config
    USE mo_limarea_config,      ONLY: latbc_config, generate_filename_mtime, LATBC_TYPE_EXT
    USE mo_ext_data_types,      ONLY: t_external_data
    USE mo_run_config,          ONLY: iqv, iqc, iqi, iqr, iqs, ltransport, dtime, nsteps, msg_level
    USE mo_dynamics_config,     ONLY: nnow, nnow_rcf
    USE mo_initicon_config,     ONLY: init_mode
    USE mo_initicon_types,      ONLY: t_initicon_state
    USE mo_cdi,                 ONLY: streamOpenRead, streamClose, streamInqVlist, vlistInqTaxis, &
      &                               taxisInqVDate, taxisInqVTime, cdiDecodeTime, cdiDecodeDate
    USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE
    USE mo_util_cdi,            ONLY: cdiGetStringError
    USE mo_master_config,       ONLY: isRestart
    USE mo_fortran_tools,       ONLY: copy, init
    IMPLICIT NONE

    PRIVATE

    ! handshake subroutines
    PUBLIC :: async_pref_send_handshake
    PUBLIC :: compute_wait_for_async_pref
    PUBLIC :: async_pref_wait_for_start
    PUBLIC :: compute_start_async_pref
    PUBLIC :: compute_shutdown_async_pref
    PUBLIC :: allocate_pref_latbc_data

    PUBLIC ::  compute_init_latbc_data, async_init_latbc_data, read_latbc_data,     &
         &     update_lin_interpolation, recv_latbc_data



    !------------------------------------------------------------------------------------------------
    ! CONSTANTS
    !------------------------------------------------------------------------------------------------

    ! module name string
    CHARACTER(LEN=*), PARAMETER :: modname = 'mo_async_latbc_utils'

    ! Tags for communication between compute PEs and prefetching PEs
    INTEGER, PARAMETER :: msg_pref_start    = 56984
    INTEGER, PARAMETER :: msg_pref_done     = 26884
    INTEGER, PARAMETER :: msg_pref_shutdown = 48965
    INTEGER, PARAMETER :: msg_latbc_done    = 20883


  CONTAINS

    !-------------------------------------------------------------------------
    !>
    !! @par Revision History
    !! Initial version by S. Brdar, DWD (2013-06-13)
    !! Initial version by M.Pondkule, DWD (2014-01-27)
    !!
    SUBROUTINE allocate_pref_latbc_data(latbc, nlev_in, p_nh_state, ext_data, p_patch)

      TYPE(t_latbc_data),    INTENT(INOUT) :: latbc       !< variable buffer for latbc data
      INTEGER,               INTENT(IN)    :: nlev_in     !< no. of vertical input levels
      TYPE(t_nh_state),      INTENT(INOUT) :: p_nh_state  !< nonhydrostatic state on the global domain
      TYPE(t_external_data), INTENT(IN)    :: ext_data    !< external data on the global domain
      TYPE(t_patch),         INTENT(IN)    :: p_patch

#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::allocate_pref_latbc_data"
      INTEGER       :: tlev, nlev, nlevp1, nblks_c, nblks_e, ioper_mode, ierrstat

      ! Select operation mode determining the set of input fields to be allocated and read
      SELECT CASE (init_mode)
      CASE (MODE_COSMODE)
        ioper_mode = 2
      CASE (MODE_DWDANA, MODE_ICONVREMAP, MODE_IAU_OLD, MODE_IAU)
        ioper_mode = 3
      CASE DEFAULT
        ioper_mode = 1
      END SELECT

      ! Allocate memory for variables (3D and 2D) on work processors
      nlev    = p_patch%nlev
      nlevp1  = p_patch%nlevp1
      nblks_c = p_patch%nblks_c
      nblks_e = p_patch%nblks_e

      DO tlev = 1, 2
        ! For safety; tke is checked for being associated in vert_interp
        NULLIFY(latbc%latbc_data(tlev)%atm_in%tke)
        NULLIFY(latbc%latbc_data(tlev)%atm_in%tke_ifc)

         ! Basic icon_remap data
         ALLOCATE(latbc%latbc_data(tlev)%topography_c(nproma,nblks_c),     &
              latbc%latbc_data(tlev)%z_ifc       (nproma,nlevp1,nblks_c),  &
              latbc%latbc_data(tlev)%z_mc        (nproma,nlev,nblks_c), STAT=ierrstat)
         IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

         ! Allocate atmospheric input data
         ALLOCATE(latbc%latbc_data(tlev)%atm_in%psfc(nproma,     nblks_c), &
              latbc%latbc_data(tlev)%atm_in%phi_sfc(nproma,      nblks_c), &
              latbc%latbc_data(tlev)%atm_in%pres (nproma,nlev_in,nblks_c), &
              latbc%latbc_data(tlev)%atm_in%z3d  (nproma,nlev_in,nblks_c), &
              latbc%latbc_data(tlev)%atm_in%temp (nproma,nlev_in,nblks_c), &
              latbc%latbc_data(tlev)%atm_in%u    (nproma,nlev_in,nblks_c), &
              latbc%latbc_data(tlev)%atm_in%v    (nproma,nlev_in,nblks_c), &
              latbc%latbc_data(tlev)%atm_in%w    (nproma,nlev_in,nblks_c), &
              latbc%latbc_data(tlev)%atm_in%omega(nproma,nlev_in,nblks_c), &
              latbc%latbc_data(tlev)%atm_in%qv   (nproma,nlev_in,nblks_c), &
              latbc%latbc_data(tlev)%atm_in%qc   (nproma,nlev_in,nblks_c), &
              latbc%latbc_data(tlev)%atm_in%qi   (nproma,nlev_in,nblks_c), &
              latbc%latbc_data(tlev)%atm_in%qr   (nproma,nlev_in,nblks_c), &
              latbc%latbc_data(tlev)%atm_in%qs   (nproma,nlev_in,nblks_c), STAT=ierrstat)
         IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

         ! allocate also vn (though sometimes not needed)
         ALLOCATE(latbc%latbc_data(tlev)%atm_in%vn(nproma, nlev_in, nblks_e), STAT=ierrstat)
         IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

         IF (ioper_mode >= 2) THEN
           ALLOCATE(latbc%latbc_data(tlev)%atm_in%w_ifc(nproma,nlev_in+1,nblks_c), &
                    latbc%latbc_data(tlev)%atm_in%z3d_ifc(nproma,nlev_in+1,nblks_c), STAT=ierrstat)
           IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
         ENDIF
         IF (ioper_mode == 3) THEN
           ALLOCATE(latbc%latbc_data(tlev)%atm_in%rho(nproma,nlev_in,nblks_c),   &
                    latbc%latbc_data(tlev)%atm_in%theta_v(nproma,nlev_in,nblks_c), STAT=ierrstat)
           IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
         ENDIF

         ! Allocate atmospheric output data
         ALLOCATE(latbc%latbc_data(tlev)%atm%vn   (nproma,nlev,nblks_e), &
              latbc%latbc_data(tlev)%atm%u        (nproma,nlev,nblks_c), &
              latbc%latbc_data(tlev)%atm%v        (nproma,nlev,nblks_c), &
              latbc%latbc_data(tlev)%atm%w        (nproma,nlevp1,nblks_c), &
              latbc%latbc_data(tlev)%atm%temp     (nproma,nlev,nblks_c), &
              latbc%latbc_data(tlev)%atm%exner    (nproma,nlev,nblks_c), &
              latbc%latbc_data(tlev)%atm%pres     (nproma,nlev,nblks_c), &
              latbc%latbc_data(tlev)%atm%rho      (nproma,nlev,nblks_c), &
              latbc%latbc_data(tlev)%atm%theta_v  (nproma,nlev,nblks_c), &
              latbc%latbc_data(tlev)%atm%qv       (nproma,nlev,nblks_c), &
              latbc%latbc_data(tlev)%atm%qc       (nproma,nlev,nblks_c), &
              latbc%latbc_data(tlev)%atm%qi       (nproma,nlev,nblks_c), &
              latbc%latbc_data(tlev)%atm%qr       (nproma,nlev,nblks_c), &
              latbc%latbc_data(tlev)%atm%qs       (nproma,nlev,nblks_c), STAT=ierrstat)
         IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

         ! topography and metrics are time independent
!$OMP PARALLEL
         CALL copy(ext_data%atm%topography_c(:,:), &
              latbc%latbc_data(tlev)%topography_c(:,:))
         CALL copy(p_nh_state%metrics%z_ifc(:,:,:), &
              latbc%latbc_data(tlev)%z_ifc(:,:,:))
         CALL copy(p_nh_state%metrics%z_mc (:,:,:), &
              latbc%latbc_data(tlev)%z_mc (:,:,:))
!$OMP END PARALLEL

      END DO

#endif
    END SUBROUTINE allocate_pref_latbc_data


    !-------------------------------------------------------------------------
    !>
    !! @par Revision History
    !! Initial version by S. Brdar, DWD (2013-06-13)
    !! Modified version by M.Pondkule, DWD (2014-01-27)
    !!
    !! ** This subroutine is only called by the worker PEs. **
    !!
    SUBROUTINE compute_init_latbc_data(latbc, p_patch, p_int_state, p_nh_state, tlev)
      TYPE(t_latbc_data),     INTENT(INOUT) :: latbc
      TYPE(t_patch),          INTENT(IN)    :: p_patch
      TYPE(t_int_state),      INTENT(IN)    :: p_int_state
      TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state  !< nonhydrostatic state on the global domain
      INTEGER,                INTENT(OUT)   :: tlev

#ifndef NOMPI
      ! local variables
      REAL(wp)     :: delta_tstep_secs, delta_dtime_secs
      TYPE(timedelta), POINTER :: delta_tstep
      TYPE(datetime),  POINTER :: mtime_current
      INTEGER       :: i, add_delta, prev_latbc_tlev
      REAL(wp)      :: seconds
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::compute_init_latbc_data"
      CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  :: td_string

      IF (.NOT. my_process_is_work())  RETURN

      ! convert namelist parameter "limarea_nml/dtime_latbc" into
      ! mtime object:
      IF (latbc_config%dtime_latbc > 86400._wp) THEN
        CALL finish(routine, "Namelist setting of limarea_nml/dtime_latbc too large for mtime conversion!")
      END IF
      seconds = latbc_config%dtime_latbc*1000._wp
      CALL getPTStringFromMS(NINT(seconds,i8), td_string)
      latbc%delta_dtime => newTimedelta(td_string)

      ! create prefetching event:
      latbc%prefetchEvent => newEvent("Prefetch input", time_config%tc_startdate, &
           time_config%tc_current_date, time_config%tc_stopdate, latbc%delta_dtime)

      latbc%mtime_read  => newDatetime(time_config%tc_startdate)

      ! time interval delta_dtime_secs in seconds
      delta_dtime_secs = 86400 *INT(latbc%delta_dtime%day)    &
           &                  + 3600  *INT(latbc%delta_dtime%hour)   &
           &                  + 60    *INT(latbc%delta_dtime%minute) &
           &                  +        INT(latbc%delta_dtime%second)

      ! if there is lrestart flag then prefetch processor needs to start reading
      ! the data from the new date time of restart file. Below the time at which
      ! restart file starts is calculated and added to mtime_read which is the
      ! time step for reading the boundary data
      IF(isRestart()) THEN
         mtime_current => newDatetime(time_config%tc_current_date)
         delta_tstep => newTimedelta("PT0S")
         delta_tstep = latbc%mtime_read - mtime_current

         ! time interval delta_tstep_secs in seconds
         delta_tstep_secs = 86400 *INT(delta_tstep%day)    &
              &                  + 3600  *INT(delta_tstep%hour)   &
              &                  + 60    *INT(delta_tstep%minute) &
              &                  +        INT(delta_tstep%second)

         add_delta = FLOOR(delta_tstep_secs / delta_dtime_secs)
         ! no of times delta_dtime needs to be added to get the current time to
         ! read the boundary data file
         DO i = 1, add_delta
            latbc%mtime_read = latbc%mtime_read + latbc%delta_dtime
         ENDDO
         ! deallocating mtime and deltatime
         CALL deallocateTimedelta(delta_tstep)
         CALL deallocateDatetime(mtime_current)
      ENDIF

      ! prepare read/last indices
      tlev  = 1   ! read in the first time-level slot

      ! read first two time steps
      IF (latbc_config%init_latbc_from_fg) THEN

        prev_latbc_tlev = 3 - tlev
        tlev            = prev_latbc_tlev

        ! allocate input data for lateral boundary nudging
        CALL copy_fg_to_latbc(latbc%latbc_data, p_nh_state, tlev)

      ELSE
        CALL recv_latbc_data( latbc, p_patch, p_nh_state, p_int_state,            &
          &                   time_config%tc_current_date, lcheck_read=.FALSE., ltime_incr=.FALSE., &
          &                   tlev=tlev)

      ENDIF
      CALL recv_latbc_data( latbc, p_patch, p_nh_state, p_int_state,           &
        &                   time_config%tc_current_date, lcheck_read=.FALSE., ltime_incr=.TRUE., &
        &                   tlev=tlev)

      CALL message(routine,'done')
#endif
    END SUBROUTINE compute_init_latbc_data


    !-------------------------------------------------------------------------
    !>
    !! @par Revision History
    !! Initial version by S. Brdar, DWD (2013-06-13)
    !! Modified version by M.Pondkule, DWD (2014-01-27)
    !!
    !! ** This subroutine is only called by the prefetching PEs. **
    !!
    SUBROUTINE async_init_latbc_data(latbc)
      TYPE(t_latbc_data), INTENT(INOUT) :: latbc

#ifndef NOMPI
      ! local variables
      REAL(wp)     :: delta_tstep_secs, delta_dtime_secs
      TYPE(timedelta), POINTER :: delta_tstep
      TYPE(datetime),  POINTER :: mtime_current
      LOGICAL       :: done
      INTEGER       :: i, add_delta
      REAL(wp)      :: seconds
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::async_init_latbc_data"
      CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  :: td_string

      IF (.NOT. my_process_is_pref())  RETURN

      ! convert namelist parameter "limarea_nml/dtime_latbc" into
      ! mtime object:
      IF (latbc_config%dtime_latbc > 86400._wp) THEN
        CALL finish(routine, "Namelist setting of limarea_nml/dtime_latbc too large for mtime conversion!")
      END IF
      seconds = latbc_config%dtime_latbc*1000._wp
      CALL getPTStringFromMS(NINT(seconds,i8), td_string)
      latbc%delta_dtime => newTimedelta(td_string)

      ! create prefetching event:
      latbc%prefetchEvent => newEvent("Prefetch input", time_config%tc_exp_startdate, &
           time_config%tc_current_date, time_config%tc_stopdate, latbc%delta_dtime)

      latbc%mtime_read  => newDatetime(time_config%tc_exp_startdate)

      ! time interval delta_dtime_secs in seconds
      delta_dtime_secs = 86400 *INT(latbc%delta_dtime%day)    &
           &                  + 3600  *INT(latbc%delta_dtime%hour)   &
           &                  + 60    *INT(latbc%delta_dtime%minute) &
           &                  +        INT(latbc%delta_dtime%second)

      ! if there is lrestart flag then prefetch processor needs to start reading
      ! the data from the new date time of restart file. Below the time at which
      ! restart file starts is calculated and added to mtime_read which is the
      ! time step for reading the boundary data
      IF(isRestart()) THEN
         mtime_current => newDatetime(time_config%tc_current_date)
         delta_tstep => newTimedelta("PT0S")
         delta_tstep = latbc%mtime_read - mtime_current

         ! time interval delta_tstep_secs in seconds
         delta_tstep_secs = 86400 *INT(delta_tstep%day)    &
              &                  + 3600  *INT(delta_tstep%hour)   &
              &                  + 60    *INT(delta_tstep%minute) &
              &                  +        INT(delta_tstep%second)

         add_delta = FLOOR(delta_tstep_secs / delta_dtime_secs)
         ! no of times delta_dtime needs to be added to get the current time to
         ! read the boundary data file
         DO i = 1, add_delta
            latbc%mtime_read = latbc%mtime_read + latbc%delta_dtime
         ENDDO
         ! deallocating mtime and deltatime
         CALL deallocateTimedelta(delta_tstep)
         CALL deallocateDatetime(mtime_current)
      ENDIF

      ! read first two time steps
      IF (.NOT. latbc_config%init_latbc_from_fg) THEN

        CALL read_latbc_data(latbc, ltime_incr=.FALSE.)

        ! Inform compute PEs that we are done
        CALL async_pref_send_handshake()
        ! Wait for a message from the compute PEs to start
        CALL async_pref_wait_for_start(done)

      ENDIF

      CALL read_latbc_data(latbc, ltime_incr=.TRUE.)
      ! Inform compute PEs that we are done
      CALL async_pref_send_handshake()

      CALL message(routine,'done')
#endif
    END SUBROUTINE async_init_latbc_data


    !-------------------------------------------------------------------------
    !>
    !! Read horizontally interpolated atmospheric boundary data.
    !!
    !! The subroutine reads atmospheric boundary data and projects on
    !! the ICON global grid. 
    !!
    !! The following steps are performed:
    !! - Read atmospheric input data,
    !! - Write input data to memory window buffer. The offset for data
    !!   is set such that each of dataset belongs to the respective
    !!   compute processor,
    !!
    !! ** This subroutine is only called by the prefetching PE. **
    !!
    !! @par Revision History
    !! Initial version by S. Brdar, DWD (2013-07-19)
    !! Modified version by M. Pondkule, DWD (2014-02-11)
    !!

    SUBROUTINE read_latbc_data(latbc, ltime_incr)

      TYPE(t_latbc_data), INTENT(IN)    :: latbc

      ! flag to increment present datetime to next time level
      LOGICAL,            INTENT(IN)    :: ltime_incr

#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::read_latbc_data"

      INTEGER(KIND=MPI_ADDRESS_KIND)        :: ioff(0:num_work_procs-1)
      INTEGER                               :: jm, latbc_fileID, vlistID, taxisID, errno,  &
        &                                      idate, iyear, imonth, iday, itime, ihour,   &
        &                                      iminute, isecond, nlevs_read
      LOGICAL                               :: l_exist
      CHARACTER(LEN=filename_max)           :: latbc_filename, latbc_full_filename
      CHARACTER(LEN=132)                    :: message_text
      CHARACTER(LEN=MAX_CHAR_LENGTH)        :: cdiErrorText
      CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: dstringA, dstringB
      TYPE(datetime), POINTER               :: mtime_vdate


      ! Prepare the mtime_read for the next time level
      IF (ltime_incr) THEN
         latbc%mtime_read = latbc%mtime_read + latbc%delta_dtime
      ENDIF

      ! return if mtime_read is at least one full boundary data
      ! interval beyond the simulation end, implying that no further
      ! data are required for correct results
      IF (latbc%mtime_read >= time_config%tc_stopdate + latbc%delta_dtime) RETURN

      latbc_filename = generate_filename_mtime(nroot, latbc%patch_data%level, &
        &                                      latbc%mtime_read, time_config%tc_exp_startdate)
      latbc_full_filename = TRIM(latbc_config%latbc_path)//TRIM(latbc_filename)

      WRITE(0,*) 'reading boundary data: ', TRIM(latbc_filename)
      INQUIRE (FILE=TRIM(ADJUSTL(latbc_full_filename)), EXIST=l_exist)
      IF (.NOT. l_exist) THEN
         WRITE (message_text,'(a,a)') 'file not found:', TRIM(latbc_filename)
         CALL finish(routine, message_text)
      ENDIF

      ! opening and reading file
      latbc_fileID  = streamOpenRead(TRIM(latbc_full_filename))
      ! check if the file could be opened
      IF (latbc_fileID < 0) THEN
         CALL cdiGetStringError(latbc_fileID, cdiErrorText)
         WRITE(message_text,'(4a)') 'File ', TRIM(latbc_full_filename), &
              ' cannot be opened: ', TRIM(cdiErrorText)
         CALL finish(routine, TRIM(message_text))
      ENDIF

      ! consistency check: Make sure that the requested date is
      ! actually contained in the file.
      
      vlistID = streamInqVlist(latbc_fileID)
      taxisID = vlistInqTaxis(vlistID)
      idate   = taxisInqVDate(taxisID)
      itime   = taxisInqVTime(taxisID)
      CALL cdiDecodeDate(idate, iyear, imonth, iday)
      CALL cdiDecodeTime(itime, ihour, iminute, isecond)
      mtime_vdate => newDatetime(iyear, imonth, iday, ihour, iminute, isecond, 0, errno)
      IF (mtime_vdate /= latbc%mtime_read) THEN
        CALL finish(routine, "requested date does not match file '"//TRIM(latbc_full_filename))
      END IF
      IF (msg_level >= 10) THEN
        CALL datetimeToString(latbc%mtime_read, dstringA)
        CALL datetimeToString(mtime_vdate, dstringB)
        WRITE (0,*) "  read date: ", TRIM(dstringA)
        WRITE (0,*) "  file date: ", TRIM(dstringB)
      END IF
      CALL deallocateDatetime(mtime_vdate)

      ! initializing the displacement array for each compute processor
      ioff(:) = 0_MPI_ADDRESS_KIND
      !
      ! Perform CDI read operation
      !
      DO jm = 1, latbc%buffer%ngrp_vars
         ! Get pointer to appropriate reorder_info
         IF(latbc%buffer%nlev(jm) /= 1 ) THEN
            SELECT CASE (latbc%buffer%hgrid(jm))
            CASE(GRID_UNSTRUCTURED_CELL)
               ! Read 3d variables
               CALL prefetch_cdi_3d ( latbc_fileID, latbc%buffer%mapped_name(jm), latbc, &
                    &                 nlevs_read, latbc%buffer%hgrid(jm), ioff )

            CASE(GRID_UNSTRUCTURED_EDGE)
               CALL prefetch_cdi_3d ( latbc_fileID, latbc%buffer%mapped_name(jm), latbc, &
                    &                 nlevs_read, latbc%buffer%hgrid(jm), ioff )
            CASE default
               CALL finish(routine,'unknown grid type')
            END SELECT

            ! consistency check
            IF (nlevs_read /= latbc%buffer%nlev(jm)) THEN
              CALL finish(routine, "variable '"//TRIM(latbc%buffer%mapped_name(jm))//&
                &"': nlev="//TRIM(int2string(nlevs_read,'(i0)'))//" but expected "//&
                &TRIM(int2string(latbc%buffer%nlev(jm),'(i0)')))
            END IF

         ELSE
            SELECT CASE (latbc%buffer%hgrid(jm))
            CASE(GRID_UNSTRUCTURED_CELL)
               ! Read 2d variables
               CALL prefetch_cdi_2d ( latbc_fileID, latbc%buffer%mapped_name(jm), latbc, &
                    &                 latbc%buffer%hgrid(jm), ioff )

            CASE(GRID_UNSTRUCTURED_EDGE)
               CALL prefetch_cdi_2d ( latbc_fileID, latbc%buffer%mapped_name(jm), latbc, &
                    &                 latbc%buffer%hgrid(jm), ioff )

            CASE default
               CALL finish(routine,'unknown grid type')
            END SELECT

         ENDIF
      ENDDO

      ! close the open dataset file
      CALL streamClose(latbc_fileID)   
#endif
    END SUBROUTINE read_latbc_data


    !-------------------------------------------------------------------------
    !>
    !! Receive horizontally interpolated atmospheric boundary data
    !! from the prefetching PE. 
    !!
    !! ** This subroutine is only called by the worker PEs. **
    !!
    !! @par Revision History
    !! Initial version by S. Brdar, DWD (2013-07-19)
    !! Modified version by M. Pondkule, DWD (2014-02-11)
    !!

    SUBROUTINE recv_latbc_data(latbc, p_patch, p_nh_state, p_int, cur_datetime,&
      &                        lcheck_read, ltime_incr, tlev)

      TYPE(t_latbc_data),     INTENT(INOUT) :: latbc
      TYPE(t_patch),          INTENT(IN)    :: p_patch
      TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state   !< nonhydrostatic state on the global domain
      TYPE(t_int_state),      INTENT(IN)    :: p_int
      TYPE(datetime), POINTER, INTENT(IN)   :: cur_datetime !< current time
      INTEGER,                INTENT(INOUT) :: tlev

      ! flag to set wether compute processor need to read boundary
      ! data or need not and than they will return
      LOGICAL,      INTENT(IN)    :: lcheck_read

      ! flag to increment present datetime to next time level
      LOGICAL,      INTENT(IN)    :: ltime_incr

#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::recv_latbc_data"
      LOGICAL                               :: isactive
      CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: sim_cur
      REAL(wp)                              :: tdiff
      CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  :: tdiff_string
      TYPE(timedelta), POINTER              :: my_duration_slack
      INTEGER                               :: prev_latbc_tlev
      TYPE(datetime),  POINTER              :: mtime_date

      ! check for event been active
      my_duration_slack => newTimedelta("PT0S")
      my_duration_slack = time_config%tc_dt_model*0.5_wp
      isactive = isCurrentEventActive(latbc%prefetchEvent, cur_datetime, my_duration_slack)
      CALL deallocateTimedelta(my_duration_slack)

      ! do we need to read boundary data
      IF (lcheck_read .AND. (.NOT. isactive))  RETURN

      ! compute processors wait for msg from
      ! prefetch processor that they can start
      ! reading latbc data from memory window
      IF(.NOT. my_process_is_io() .AND. &
        & .NOT. my_process_is_mpi_test()) THEN
        CALL compute_wait_for_async_pref()
      END IF

      ! Prepare the mtime_read for the next time level
      IF (ltime_incr) THEN
         latbc%mtime_read = latbc%mtime_read + latbc%delta_dtime
      ENDIF

      ! Adjust read/last indices
      !
      ! New boundary data time-level is always read in latbc%latbc_data(tlev),
      ! whereas latbc%latbc_data(prev_latbc_tlev) always holds the last read boundary data
      !
      prev_latbc_tlev = 3 - tlev
      tlev  = prev_latbc_tlev

      ! start reading boundary data
      IF (latbc_config%itype_latbc == LATBC_TYPE_EXT) THEN
        CALL compute_latbc_intp_data( latbc, p_patch, p_nh_state, p_int, tlev )
      ELSE
        CALL compute_latbc_icon_data( latbc, p_patch, p_int, tlev )
      END IF

      ! Compute tendencies for nest boundary update
      IF (ltime_incr) CALL compute_boundary_tendencies(latbc%latbc_data, p_patch, p_nh_state, tlev)
      
#endif
    END SUBROUTINE recv_latbc_data


    !-------------------------------------------------------------------------
    !>
    !! Copy data read by the prefetch_latbc_intp_data routine
    !! (see comments above for the functionality of this processing mode)
    !!
    !! This subroutine is called by compute processors.
    !! The following steps are performed:
    !! - Copy from memory window buffer, atmospheric IFS analysis data,
    !! - compute the prognostic NH variable set (no vertical interpolation is performed!)
    !!
    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2014-05-19)
    !!
    SUBROUTINE compute_latbc_icon_data( latbc, p_patch, p_int, tlev )
      TYPE(t_latbc_data),     INTENT(INOUT), TARGET :: latbc
      TYPE(t_patch),          INTENT(IN), TARGET    :: p_patch
      TYPE(t_int_state),      INTENT(IN)            :: p_int
      INTEGER,                INTENT(IN)            :: tlev

#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::compute_latbc_icon_data"
      INTEGER(i8)                         :: eoff
      TYPE(t_reorder_data), POINTER       :: p_ri
      REAL(wp)                            :: temp_v(nproma,p_patch%nlev,p_patch%nblks_c)
      INTEGER                             :: jc, jk, jb, jm, j, jl, jv
      !-------------------------------------------------------------------------

      ! consistency check
      IF ((tlev <1) .OR. (SIZE(latbc%latbc_data) < tlev)) THEN
        CALL finish(routine, "Internal error!")
      END IF

      ! Offset in memory window for async prefetching
      eoff = 0_i8
      DO jv = 1, latbc%buffer%ngrp_vars
         ! Receive 2d and 3d variables
         CALL compute_data_receive (latbc%buffer%hgrid(jv), latbc%buffer%nlev(jv), &
                                    latbc%buffer%vars(jv)%buffer, eoff, latbc%patch_data)
      ENDDO

      ! Reading the next time step
#ifndef NOMPI
      IF((.NOT. my_process_is_io() .AND. &
           & .NOT. my_process_is_pref()) .AND. &
           & .NOT. my_process_is_mpi_test()) THEN
         CALL compute_start_async_pref()
      ENDIF
#endif

      ! Get patch ID
      p_ri => latbc%patch_data%cells

      !
      ! get prognostic 3d fields
      !
      ! copying tha variable values from prefetch buffer to the respective allocated variable

!$OMP PARALLEL PRIVATE(jm,jv,jc)
      jm = get_field_index(latbc, 'temp')
      IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'temp'!")

      ! consistency check
      IF ((SIZE(latbc%latbc_data(tlev)%atm%temp,2) < latbc%buffer%nlev(jm)) .OR.  &
        & (SIZE(latbc%buffer%vars(jm)%buffer,2)    < latbc%buffer%nlev(jm))) THEN
        CALL finish(routine, "Internal error!")
      END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc%buffer%nlev(jm)
         DO j = 1, p_ri%n_own ! p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            !       WRITE(0,*) 'compute_latbc_icon_data 01 ',latbc%latbc_data(tlev)%atm%temp(jl,jk,jb)
            latbc%latbc_data(tlev)%atm%temp(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO  !NOWAIT

      jm = get_field_index(latbc, 'u')
      IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'u'!")
      jv = get_field_index(latbc, 'v')
      IF (jv <= 0)  CALL finish(routine, "Internal error, invalid field index for 'v'!")

      ! consistency check
      IF ((SIZE(latbc%latbc_data(tlev)%atm_in%u,2) < latbc%buffer%nlev(jm)) .OR.  &
        & (SIZE(latbc%buffer%vars(jm)%buffer,2)    < latbc%buffer%nlev(jm))) THEN
        CALL finish(routine, "Internal error!")
      END IF
      IF ((SIZE(latbc%latbc_data(tlev)%atm_in%v,2) < latbc%buffer%nlev(jv)) .OR.  &
        & (SIZE(latbc%buffer%vars(jv)%buffer,2)    < latbc%buffer%nlev(jv))) THEN
        CALL finish(routine, "Internal error!")
      END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc%buffer%nlev(jm)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm_in%u(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
            latbc%latbc_data(tlev)%atm_in%v(jl,jk,jb) = REAL(latbc%buffer%vars(jv)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      jv = get_field_index(latbc, 'w')
      IF (jv <= 0)  CALL finish(routine, "Internal error, invalid field index for 'w'!")

      ! consistency check
      IF ((SIZE(latbc%latbc_data(tlev)%atm%w,2) < latbc%buffer%nlev(jv)) .OR.  &
        & (SIZE(latbc%buffer%vars(jv)%buffer,2) < latbc%buffer%nlev(jv))) THEN
        CALL finish(routine, "Internal error!")
      END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc%buffer%nlev(jv)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm%w(jl,jk,jb) = REAL(latbc%buffer%vars(jv)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      ! Read parameter Pressure
      jv = get_field_index(latbc, 'pres')
      IF (jv <= 0)  CALL finish(routine, "Internal error, invalid field index for 'pres'!")

      ! consistency check
      IF ((SIZE(latbc%latbc_data(tlev)%atm%pres,2) < latbc%buffer%nlev(jv)) .OR.  &
        & (SIZE(latbc%buffer%vars(jv)%buffer,2) < latbc%buffer%nlev(jv))) THEN
        CALL finish(routine, "Internal error!")
      END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc%buffer%nlev(jv)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm%pres(jl,jk,jb) = REAL(latbc%buffer%vars(jv)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      ! Read parameter qv
      jv = get_field_index(latbc, 'qv')
      IF (jv <= 0)  CALL finish(routine, "Internal error, invalid field index for 'qv'!")

      ! consistency check
      IF ((SIZE(latbc%latbc_data(tlev)%atm%qv,2) < latbc%buffer%nlev(jv)) .OR.  &
        & (SIZE(latbc%buffer%vars(jv)%buffer,2) < latbc%buffer%nlev(jv))) THEN
        CALL finish(routine, "Internal error!")
      END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc%buffer%nlev(jv)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm%qv(jl,jk,jb) = REAL(latbc%buffer%vars(jv)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      jc = get_field_index(latbc, 'qc')
      IF (jc <= 0)  CALL finish(routine, "Internal error, invalid field index for 'qc'!")

      ! consistency check
      IF ((SIZE(latbc%latbc_data(tlev)%atm%qc,2) < latbc%buffer%nlev(jc)) .OR.  &
        & (SIZE(latbc%buffer%vars(jc)%buffer,2) < latbc%buffer%nlev(jc))) THEN
        CALL finish(routine, "Internal error!")
      END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc%buffer%nlev(jc)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm%qc(jl,jk,jb) = REAL(latbc%buffer%vars(jc)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      ! Read parameter qi
      jm = get_field_index(latbc, 'qi')
      IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'qi'!")

      ! consistency check
      IF ((SIZE(latbc%latbc_data(tlev)%atm%qi,2) < latbc%buffer%nlev(jm)) .OR.  &
        & (SIZE(latbc%buffer%vars(jm)%buffer,2) < latbc%buffer%nlev(jm))) THEN
        CALL finish(routine, "Internal error!")
      END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc%buffer%nlev(jm)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm%qi(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      ! Read parameter qr
      jm = get_field_index(latbc, 'qr')
      IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'qr'!")

      ! consistency check
      IF ((SIZE(latbc%latbc_data(tlev)%atm%qr,2) < latbc%buffer%nlev(jm)) .OR.  &
        & (SIZE(latbc%buffer%vars(jm)%buffer,2) < latbc%buffer%nlev(jm))) THEN
        CALL finish(routine, "Internal error!")
      END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc%buffer%nlev(jm)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm%qr(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      ! Read parameter qs
      jv = get_field_index(latbc, 'qs')
      IF (jv <= 0)  CALL finish(routine, "Internal error, invalid field index for 'qs'!")

      ! consistency check
      IF ((SIZE(latbc%latbc_data(tlev)%atm%qs,2) < latbc%buffer%nlev(jv)) .OR.  &
        & (SIZE(latbc%buffer%vars(jv)%buffer,2) < latbc%buffer%nlev(jv))) THEN
        CALL finish(routine, "Internal error!")
      END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc%buffer%nlev(jv)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm%qs(jl,jk,jb) = REAL(latbc%buffer%vars(jv)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      ! boundary exchange for a 2-D and 3-D array to fill HALO region.
      ! This addition by M.Pondkule, DWD (11/06/2014)
      CALL sync_patch_array(SYNC_C, p_patch, latbc%latbc_data(tlev)%atm%temp)
      CALL sync_patch_array(SYNC_C, p_patch, latbc%latbc_data(tlev)%atm_in%u)
      CALL sync_patch_array(SYNC_C, p_patch, latbc%latbc_data(tlev)%atm_in%v)
      CALL sync_patch_array(SYNC_C, p_patch, latbc%latbc_data(tlev)%atm%w)
      CALL sync_patch_array(SYNC_C, p_patch, latbc%latbc_data(tlev)%atm%pres)
      CALL sync_patch_array(SYNC_C, p_patch, latbc%latbc_data(tlev)%atm%qv)
      CALL sync_patch_array(SYNC_C, p_patch, latbc%latbc_data(tlev)%atm%qc)
      CALL sync_patch_array(SYNC_C, p_patch, latbc%latbc_data(tlev)%atm%qi)
      CALL sync_patch_array(SYNC_C, p_patch, latbc%latbc_data(tlev)%atm%qr)
      CALL sync_patch_array(SYNC_C, p_patch, latbc%latbc_data(tlev)%atm%qs)
      CALL sync_patch_array(SYNC_C, p_patch, latbc%latbc_data(tlev)%atm%rho)
      CALL sync_patch_array(SYNC_C, p_patch, latbc%latbc_data(tlev)%atm%theta_v)
      !
      ! Convert u and v on cell points to vn at edge points
      !
      CALL interp_uv_2_vn( p_patch, p_int, latbc%latbc_data(tlev)%atm_in%u,              &
           &                  latbc%latbc_data(tlev)%atm_in%v, latbc%latbc_data(tlev)%atm%vn )

      !
      ! Compute virtual temperature
      !
      CALL virtual_temp( p_patch, latbc%latbc_data(tlev)%atm%temp, latbc%latbc_data(tlev)%atm%qv, &
           &                latbc%latbc_data(tlev)%atm%qc, latbc%latbc_data(tlev)%atm%qi,               &
           &                latbc%latbc_data(tlev)%atm%qr, latbc%latbc_data(tlev)%atm%qs,               &
           &                temp_v=temp_v )

      !
      ! Compute NH prognostic thermodynamical variables
      !
      CALL convert_thdvars( p_patch, latbc%latbc_data(tlev)%atm%pres, temp_v,              &
           &                   latbc%latbc_data(tlev)%atm%rho,                                   &
           &                   latbc%latbc_data(tlev)%atm%exner,                                 &
           &                   latbc%latbc_data(tlev)%atm%theta_v )

      CALL sync_patch_array(SYNC_E, p_patch, latbc%latbc_data(tlev)%atm%vn)

#endif
    END SUBROUTINE compute_latbc_icon_data
    !-------------------------------------------------------------------------


    !-------------------------------------------------------------------------
    !>
    !! Copy horizontally interpolated atmospheric analysis or forecast data
    !!
    !! This subroutine is called by compute processors.
    !! The following steps are performed:
    !! - Copy from memory window buffer, atmospheric analysis data,
    !! - interpolate vertically from intermediate remapicon grid to ICON grid
    !!   and compute the prognostic NH variable set,
    !!
    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2014-05-19)
    !!
    SUBROUTINE compute_latbc_intp_data( latbc, p_patch, p_nh_state, p_int, tlev)
      TYPE(t_latbc_data),     INTENT(INOUT), TARGET :: latbc  !< variable buffer for latbc data
      TYPE(t_patch),          INTENT(IN), TARGET    :: p_patch
      TYPE(t_nh_state),       INTENT(IN)            :: p_nh_state  !< nonhydrostatic state on the global domain
      TYPE(t_int_state),      INTENT(IN)            :: p_int
      INTEGER,                INTENT(IN)            :: tlev

#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::compute_latbc_intp_data"
      INTEGER(i8)                         :: eoff
      TYPE(t_reorder_data), POINTER       :: p_ri
      LOGICAL                             :: lconvert_omega2w
      INTEGER                             :: jc, jk, jb, jm, j, jl, jv
      INTEGER                             :: ioper_mode
      REAL(wp)                            :: log_exner, tempv
      INTEGER                             :: nlev     ! number of vertical levels in the boundary data

      ! consistency check
      IF ((tlev <1) .OR. (SIZE(latbc%latbc_data) < tlev)) THEN
        CALL finish(routine, "Internal error!")
      END IF

      ! Select operation mode determining the set of input fields to be allocated and read
      SELECT CASE (init_mode)
      CASE (MODE_COSMODE)
        ioper_mode = 2
      CASE (MODE_DWDANA, MODE_ICONVREMAP, MODE_IAU_OLD, MODE_IAU)
        ioper_mode = 3
      CASE DEFAULT
        ioper_mode = 1
      END SELECT

      ! Offset in memory window for async prefetching
      eoff = 0_i8

      DO jv = 1, latbc%buffer%ngrp_vars
         ! Receive 2d and 3d variables
         CALL compute_data_receive ( latbc%buffer%hgrid(jv), latbc%buffer%nlev(jv), &
                                     latbc%buffer%vars(jv)%buffer, eoff, latbc%patch_data)
      ENDDO

      IF (latbc%buffer%lthd_progvars .AND. ioper_mode /= 3) &
        CALL finish(routine, 'Lateral boundary variables not consistent with init_mode')

      ! Reading the next time step
      IF((.NOT. my_process_is_io() .AND. &
           & .NOT. my_process_is_pref()) .AND. &
           & .NOT. my_process_is_mpi_test()) THEN
         CALL compute_start_async_pref()
      END IF

      ! Get patch ID
      p_ri => latbc%patch_data%cells

      ! retrieve the no. of vertical levels from the temp/theta_v variable:
      IF (latbc%buffer%lthd_progvars) THEN
        jm = get_field_index(latbc, 'theta_v')
        IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'theta_v'!")
        nlev = latbc%buffer%nlev(jm)
      ELSE
        jm = get_field_index(latbc, 'temp')
        IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'temp'!")
        nlev = latbc%buffer%nlev(jm)
      END IF

      ! copying the variable values from prefetch buffer to the
      ! respective allocated variable

!$OMP PARALLEL PRIVATE(jm,jv,jc)

      IF (latbc%buffer%lthd_progvars) THEN
        jm = get_field_index(latbc, 'theta_v')
        IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'theta_v'!")

        ! consistency check
        IF ((SIZE(latbc%latbc_data(tlev)%atm_in%theta_v,2) < latbc%buffer%nlev(jm)) .OR.  &
          & (SIZE(latbc%buffer%vars(jm)%buffer,2)    < latbc%buffer%nlev(jm))) THEN
          CALL finish(routine, "Internal error!")
        END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
        DO jk=1, nlev
          DO j = 1, p_ri%n_own ! p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm_in%theta_v(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
          ENDDO
        ENDDO
!$OMP END DO

        jm = get_field_index(latbc, 'rho')
        IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'rho'!")

        ! consistency check
        IF ((SIZE(latbc%latbc_data(tlev)%atm_in%rho,2) < latbc%buffer%nlev(jm)) .OR.  &
          & (SIZE(latbc%buffer%vars(jm)%buffer,2)    < latbc%buffer%nlev(jm))) THEN
          CALL finish(routine, "Internal error!")
        END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
        DO jk=1, latbc%buffer%nlev(jm)
          DO j = 1, p_ri%n_own ! p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm_in%rho(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
          ENDDO
        ENDDO
!$OMP END DO

      ELSE

        jm = get_field_index(latbc, 'temp')
        IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'temp'!")

        ! consistency check
        IF ((SIZE(latbc%latbc_data(tlev)%atm_in%temp,2) < latbc%buffer%nlev(jm)) .OR.  &
          & (SIZE(latbc%buffer%vars(jm)%buffer,2)    < latbc%buffer%nlev(jm))) THEN
          CALL finish(routine, "Internal error!")
        END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
        DO jk=1, nlev
          DO j = 1, p_ri%n_own ! p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm_in%temp(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
          ENDDO
        ENDDO
!$OMP END DO
      ENDIF

      ! Read horizontal component of velocity (U and V)
      IF (latbc%buffer%lread_vn) THEN
         jm = get_field_index(latbc, 'vn')
         IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'vn'!")

        ! consistency check
        IF ((SIZE(latbc%latbc_data(tlev)%atm_in%vn,2) < latbc%buffer%nlev(jm)) .OR.  &
          & (SIZE(latbc%buffer%vars(jm)%buffer,2)    < latbc%buffer%nlev(jm))) THEN
          CALL finish(routine, "Internal error!")
        END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc%buffer%nlev(jm)
            DO j = 1, latbc%patch_data%edges%n_own !p_patch%n_patch_cells
               jb = latbc%patch_data%edges%own_blk(j) ! Block index in distributed patch
               jl = latbc%patch_data%edges%own_idx(j) ! Line  index in distributed patch
               latbc%latbc_data(tlev)%atm_in%vn(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO
      ELSE
         jm = get_field_index(latbc, 'u')
         IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'u'!")
         jv = get_field_index(latbc, 'v')
         IF (jv <= 0)  CALL finish(routine, "Internal error, invalid field index for 'v'!")

         ! consistency check
         IF ((SIZE(latbc%latbc_data(tlev)%atm_in%u,2) < latbc%buffer%nlev(jm)) .OR.  &
           & (SIZE(latbc%buffer%vars(jm)%buffer,2)    < latbc%buffer%nlev(jm))) THEN
           CALL finish(routine, "Internal error!")
         END IF
         IF ((SIZE(latbc%latbc_data(tlev)%atm_in%v,2) < latbc%buffer%nlev(jv)) .OR.  &
           & (SIZE(latbc%buffer%vars(jv)%buffer,2)    < latbc%buffer%nlev(jv))) THEN
           CALL finish(routine, "Internal error!")
         END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc%buffer%nlev(jm)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc%latbc_data(tlev)%atm_in%u(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
               latbc%latbc_data(tlev)%atm_in%v(jl,jk,jb) = REAL(latbc%buffer%vars(jv)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO
      ENDIF

      ! Read vertical component of velocity (W)
      IF (ioper_mode == 1) THEN
         lconvert_omega2w = .TRUE.
         jv = get_field_index(latbc, 'w')
         IF (jv <= 0)  CALL finish(routine, "Internal error, invalid field index for 'w'!")

         ! consistency check
         IF ((SIZE(latbc%latbc_data(tlev)%atm_in%omega,2) < latbc%buffer%nlev(jv)) .OR.  &
           & (SIZE(latbc%buffer%vars(jv)%buffer,2) < latbc%buffer%nlev(jv))) THEN
           CALL finish(routine, "Internal error!")
         END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc%buffer%nlev(jv)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc%latbc_data(tlev)%atm_in%omega(jl,jk,jb) = REAL(latbc%buffer%vars(jv)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO
      ELSE
         lconvert_omega2w = .FALSE.
         jv = get_field_index(latbc, 'w')
         IF (jv <= 0)  CALL finish(routine, "Internal error, invalid field index for 'w'!")

         ! consistency check
         IF ((SIZE(latbc%latbc_data(tlev)%atm_in%w_ifc,2) < latbc%buffer%nlev(jv)) .OR.  &
           & (SIZE(latbc%buffer%vars(jv)%buffer,2) < latbc%buffer%nlev(jv))) THEN
           CALL finish(routine, "Internal error!")
         END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc%buffer%nlev(jv)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc%latbc_data(tlev)%atm_in%w_ifc(jl,jk,jb) = REAL(latbc%buffer%vars(jv)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO
      ENDIF

      IF (ioper_mode >= 2) THEN
         ! Read parameter HHL
         jm = get_field_index(latbc, 'z_ifc')
         IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'z_ifc'!")

         ! consistency check
         IF ((SIZE(latbc%latbc_data(tlev)%atm_in%z3d_ifc,2) < latbc%buffer%nlev(jm)) .OR.  &
           & (SIZE(latbc%buffer%vars(jm)%buffer,2) < latbc%buffer%nlev(jm))) THEN
           CALL finish(routine, "Internal error!")
         END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc%buffer%nlev(jm)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc%latbc_data(tlev)%atm_in%z3d_ifc(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO

         ! Interpolate input 'z3d' and 'w' from the interface levels to the main levels
         !

         ! consistency check
         IF ((SIZE(latbc%latbc_data(tlev)%atm_in%w,2) < nlev) .OR.  &
           & (SIZE(latbc%latbc_data(tlev)%atm_in%w_ifc,2) < (nlev+1))) THEN
           CALL finish(routine, "Internal error!")
         END IF
         IF ((SIZE(latbc%latbc_data(tlev)%atm_in%z3d,2) < nlev) .OR.  &
           & (SIZE(latbc%latbc_data(tlev)%atm_in%z3d_ifc,2) < (nlev+1))) THEN
           CALL finish(routine, "Internal error!")
         END IF

!$OMP DO PRIVATE (jk,j,jb,jc) ICON_OMP_DEFAULT_SCHEDULE
         DO jk = 1, nlev    !!!!!!!!need to reset nlev_in to a new value from stored n_lev values
            DO j = 1, p_ri%n_own
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jc = p_ri%own_idx(j) ! Line  index in distributed patch

               ! Note: In future, we want to z3d from boundary data.
               !
               latbc%latbc_data(tlev)%atm_in%z3d(jc,jk,jb) = (latbc%latbc_data(tlev)%atm_in%z3d_ifc(jc,jk,jb) + &
                    &   latbc%latbc_data(tlev)%atm_in%z3d_ifc(jc,jk+1,jb) ) * 0.5_wp
               latbc%latbc_data(tlev)%atm_in%w(jc,jk,jb) = (latbc%latbc_data(tlev)%atm_in%w_ifc(jc,jk,jb) + &
                    &   latbc%latbc_data(tlev)%atm_in%w_ifc(jc,jk+1,jb)) * 0.5_wp
            ENDDO
         ENDDO
!$OMP END DO
      ENDIF

      ! Read parameter QV, QC and QI
      jv = get_field_index(latbc, 'qv')
      IF (jv <= 0)  CALL finish(routine, "Internal error, invalid field index for 'qv'!")

      ! consistency check
      IF ((SIZE(latbc%latbc_data(tlev)%atm_in%qv,2) < latbc%buffer%nlev(jv)) .OR.  &
        & (SIZE(latbc%buffer%vars(jv)%buffer,2) < latbc%buffer%nlev(jv))) THEN
        CALL finish(routine, "Internal error!")
      END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc%buffer%nlev(jv)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm_in%qv(jl,jk,jb) = REAL(latbc%buffer%vars(jv)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      jc = get_field_index(latbc, 'qc')
      IF (jc <= 0)  CALL finish(routine, "Internal error, invalid field index for 'qc'!")

      ! consistency check
      IF ((SIZE(latbc%latbc_data(tlev)%atm_in%qc,2) < latbc%buffer%nlev(jc)) .OR.  &
        & (SIZE(latbc%buffer%vars(jc)%buffer,2) < latbc%buffer%nlev(jc))) THEN
        CALL finish(routine, "Internal error!")
      END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc%buffer%nlev(jc)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm_in%qc(jl,jk,jb) = REAL(latbc%buffer%vars(jc)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      jm = get_field_index(latbc, 'qi')
      IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'qi'!")

      ! consistency check
      IF ((SIZE(latbc%latbc_data(tlev)%atm_in%qi,2) < latbc%buffer%nlev(jm)) .OR.  &
        & (SIZE(latbc%buffer%vars(jm)%buffer,2) < latbc%buffer%nlev(jm))) THEN
        CALL finish(routine, "Internal error!")
      END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc%buffer%nlev(jm)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc%latbc_data(tlev)%atm_in%qi(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      IF (latbc%buffer%lread_qr) THEN
         ! Read parameter QR
         jm = get_field_index(latbc, 'qr')
         IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for 'qr'!")

         ! consistency check
         IF ((SIZE(latbc%latbc_data(tlev)%atm_in%qr,2) < latbc%buffer%nlev(jm)) .OR.  &
           & (SIZE(latbc%buffer%vars(jm)%buffer,2) < latbc%buffer%nlev(jm))) THEN
           CALL finish(routine, "Internal error!")
         END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc%buffer%nlev(jm)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc%latbc_data(tlev)%atm_in%qr(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO
      ELSE
        CALL init(latbc%latbc_data(tlev)%atm_in%qr(:,:,:))
!$OMP BARRIER
      ENDIF

      IF (latbc%buffer%lread_qs) THEN
         ! Read parameter QS
         jv = get_field_index(latbc, 'qs')
         IF (jv <= 0)  CALL finish(routine, "Internal error, invalid field index for 'qs'!")

         ! consistency check
         IF ((SIZE(latbc%latbc_data(tlev)%atm_in%qs,2) < latbc%buffer%nlev(jv)) .OR.  &
           & (SIZE(latbc%buffer%vars(jv)%buffer,2) < latbc%buffer%nlev(jv))) THEN
           CALL finish(routine, "Internal error!")
         END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc%buffer%nlev(jv)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc%latbc_data(tlev)%atm_in%qs(jl,jk,jb) = REAL(latbc%buffer%vars(jv)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO
      ELSE
        CALL init(latbc%latbc_data(tlev)%atm_in%qs(:,:,:))
!$OMP BARRIER
      ENDIF

      IF (ioper_mode >= 2 .AND. .NOT. latbc%buffer%lthd_progvars) THEN
         ! Read parameter Pressure
         jv = get_field_index(latbc, 'pres')
         IF (jv <= 0)  CALL finish(routine, "Internal error, invalid field index for 'pres'!")

         ! consistency check
         IF ((SIZE(latbc%latbc_data(tlev)%atm_in%pres,2) < latbc%buffer%nlev(jv)) .OR.  &
           & (SIZE(latbc%buffer%vars(jv)%buffer,2) < latbc%buffer%nlev(jv))) THEN
           CALL finish(routine, "Internal error!")
         END IF

!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc%buffer%nlev(jv)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc%latbc_data(tlev)%atm_in%pres(jl,jk,jb) = REAL(latbc%buffer%vars(jv)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO
      ENDIF

      IF (ioper_mode /= 3) THEN
        ! Read parameter surface pressure (LNPS)
        jm = get_field_index(latbc, latbc%buffer%psvar,.TRUE.)
        IF (jm <= 0) THEN
          CALL finish(routine, "Internal error, invalid field index for "//TRIM(latbc%buffer%psvar)//"!")
        END IF

        ! consistency check
        IF ((SIZE(latbc%latbc_data(tlev)%atm_in%psfc,2) < latbc%buffer%nlev(jm)) .OR.  &
          & (SIZE(latbc%buffer%vars(jm)%buffer,2) < latbc%buffer%nlev(jm))) THEN
          CALL finish(routine, "Internal error!")
        END IF

!$OMP DO PRIVATE (j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
           jb = p_ri%own_blk(j) ! Block index in distributed patch
           jl = p_ri%own_idx(j) ! Line  index in distributed patch
           latbc%latbc_data(tlev)%atm_in%psfc(jl,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,1,jb), wp)
        ENDDO
!$OMP END DO

        ! Read parameter  surface Geopotential (GEOSP)
        jm = get_field_index(latbc, latbc%buffer%geop_ml_var,.TRUE.)
        IF (jm <= 0)  CALL finish(routine, "Internal error, invalid field index for geop_ml_var!")

        ! consistency check
        IF ((SIZE(latbc%latbc_data(tlev)%atm_in%phi_sfc,2) < latbc%buffer%nlev(jm)) .OR.  &
          & (SIZE(latbc%buffer%vars(jm)%buffer,2) < latbc%buffer%nlev(jm))) THEN
          CALL finish(routine, "Internal error!")
        END IF

!$OMP DO PRIVATE (j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
           jb = p_ri%own_blk(j) ! Block index in distributed patch
           jl = p_ri%own_idx(j) ! Line  index in distributed patch
           latbc%latbc_data(tlev)%atm_in%phi_sfc(jl,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,1,jb), wp)
         ENDDO
!$OMP END DO

      ELSE IF (latbc%buffer%lthd_progvars) THEN

        ! Diagnose pres and temp from prognostic ICON variables

!$OMP DO PRIVATE (jk,j,jb,jc,log_exner,tempv)
        DO jk = 1, nlev
          DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jc = p_ri%own_idx(j) ! Line  index in distributed patch

            log_exner = (1._wp/cvd_o_rd)*LOG(latbc%latbc_data(tlev)%atm_in%rho(jc,jk,jb)* &
              latbc%latbc_data(tlev)%atm_in%theta_v(jc,jk,jb)*rd/p0ref)
            tempv = latbc%latbc_data(tlev)%atm_in%theta_v(jc,jk,jb)*EXP(log_exner)

            latbc%latbc_data(tlev)%atm_in%pres(jc,jk,jb) = p0ref*EXP(cpd/rd*log_exner)
            latbc%latbc_data(tlev)%atm_in%temp(jc,jk,jb) = tempv / (1._wp + vtmpc1*latbc%latbc_data(tlev)%atm_in%qv(jc,jk,jb) - &
              (latbc%latbc_data(tlev)%atm_in%qc(jc,jk,jb) + latbc%latbc_data(tlev)%atm_in%qi(jc,jk,jb) +                        &
               latbc%latbc_data(tlev)%atm_in%qr(jc,jk,jb) + latbc%latbc_data(tlev)%atm_in%qs(jc,jk,jb)) )

          ENDDO
        ENDDO
!$OMP END DO
      ENDIF

!$OMP END PARALLEL

      ! boundary exchange for a 2-D and 3-D array, needed because the
      ! vertical interpolation includes the halo region (otherwise, the 
      ! syncs would have to be called after vert_interp)
      CALL sync_patch_array_mult(SYNC_C,p_patch,4,latbc%latbc_data(tlev)%atm_in%temp, &
        &                        latbc%latbc_data(tlev)%atm_in%z3d,                   &
        &                        latbc%latbc_data(tlev)%atm_in%w,                     &
        &                        latbc%latbc_data(tlev)%atm_in%pres)
      CALL sync_patch_array_mult(SYNC_C,p_patch,5,latbc%latbc_data(tlev)%atm_in%qv, &
        &                        latbc%latbc_data(tlev)%atm_in%qc,                  &
        &                        latbc%latbc_data(tlev)%atm_in%qi,                  &
        &                        latbc%latbc_data(tlev)%atm_in%qr,                  &
        &                        latbc%latbc_data(tlev)%atm_in%qs)
      IF (latbc%buffer%lread_vn) THEN
         CALL sync_patch_array(SYNC_E,p_patch,latbc%latbc_data(tlev)%atm_in%vn)
      ELSE
         CALL sync_patch_array_mult(SYNC_C,p_patch,2,latbc%latbc_data(tlev)%atm_in%u, &
           &                        latbc%latbc_data(tlev)%atm_in%v)
      ENDIF
      IF (lconvert_omega2w) CALL sync_patch_array(SYNC_C,p_patch,latbc%latbc_data(tlev)%atm_in%omega)
      IF (.NOT. latbc%buffer%lthd_progvars) THEN
        CALL sync_patch_array(SYNC_C,p_patch,latbc%latbc_data(tlev)%atm_in%phi_sfc)
        CALL sync_patch_array(SYNC_C,p_patch,latbc%latbc_data(tlev)%atm_in%psfc)
      ENDIF

      ! perform vertical interpolation of horizonally interpolated analysis data
      !
      CALL vert_interp(p_patch, p_int, p_nh_state%metrics, nlev, latbc%latbc_data(tlev),              &
           &    opt_convert_omega2w=lconvert_omega2w, opt_use_vn=latbc%buffer%lread_vn)

#endif
    END SUBROUTINE compute_latbc_intp_data


    ! Wrapper routine for copying prognostic variables from initial state to the 
    ! first time level of the lateral boundary data
    !
    SUBROUTINE copy_fg_to_latbc(latbc_data, p_nh, tlev)
      TYPE(t_initicon_state), INTENT(INOUT) :: latbc_data(:)
      TYPE(t_nh_state),       INTENT(IN)    :: p_nh
      INTEGER,                INTENT(IN)    :: tlev

      INTEGER :: jg

      jg = 1

!$OMP PARALLEL
      CALL copy(p_nh%prog(nnow(jg))%vn,      latbc_data(tlev)%atm%vn)
      CALL copy(p_nh%prog(nnow(jg))%w,       latbc_data(tlev)%atm%w)
      CALL copy(p_nh%prog(nnow(jg))%rho,     latbc_data(tlev)%atm%rho)
      CALL copy(p_nh%prog(nnow(jg))%theta_v, latbc_data(tlev)%atm%theta_v)
      CALL copy(p_nh%prog(nnow_rcf(jg))%tracer(:,:,:,iqv), latbc_data(tlev)%atm%qv)
      CALL copy(p_nh%prog(nnow_rcf(jg))%tracer(:,:,:,iqc), latbc_data(tlev)%atm%qc)
      CALL copy(p_nh%prog(nnow_rcf(jg))%tracer(:,:,:,iqi), latbc_data(tlev)%atm%qi)
      CALL copy(p_nh%prog(nnow_rcf(jg))%tracer(:,:,:,iqr), latbc_data(tlev)%atm%qr)
      CALL copy(p_nh%prog(nnow_rcf(jg))%tracer(:,:,:,iqs), latbc_data(tlev)%atm%qs)
!$OMP END PARALLEL

    END SUBROUTINE copy_fg_to_latbc

    !-------------------------------------------------------------------------
    !>
    !! @par Revision History
    !! Initial version by G. Zaengl, DWD (2013-10-22)
    !!
    SUBROUTINE compute_boundary_tendencies ( latbc_data, p_patch, p_nh, tlev )
      TYPE(t_initicon_state), INTENT(IN)    :: latbc_data(:)
      TYPE(t_patch),          INTENT(IN)    :: p_patch
      TYPE(t_nh_state),       INTENT(INOUT) :: p_nh
      INTEGER,                INTENT(IN)    :: tlev

#ifndef NOMPI
      ! Local variables
      INTEGER                         :: i_startblk, i_endblk, &
           &                                i_startidx, i_endidx
      INTEGER                         :: jc, jk, jb, je
      INTEGER                         :: nlev, nlevp1, prev_latbc_tlev
      REAL(wp)                        :: rdt

      prev_latbc_tlev = 3 - tlev

      nlev            = p_patch%nlev
      nlevp1          = p_patch%nlevp1

      ! Inverse value of boundary update frequency
      rdt = 1._wp/latbc_config%dtime_latbc

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

      ! a) Boundary tendency of horizontal velocity
      i_startblk = p_patch%edges%start_blk(1,1)
      i_endblk   = p_patch%edges%end_blk(grf_bdywidth_e,1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

         CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
              i_startidx, i_endidx, 1, grf_bdywidth_e)

         DO jk = 1, nlev
            DO je = i_startidx, i_endidx
               p_nh%diag%grf_tend_vn(je,jk,jb) = rdt * (            &
                    &   latbc_data(tlev)%atm%vn(je,jk,jb) &
                    & - latbc_data(prev_latbc_tlev)%atm%vn(je,jk,jb) )
            ENDDO
         ENDDO
      ENDDO
!$OMP END DO

      ! b) Boundary tendencies of variables at cell centers
      i_startblk = p_patch%cells%start_blk(1,1)
      i_endblk   = p_patch%cells%end_blk(grf_bdywidth_c,1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
              i_startidx, i_endidx, 1, grf_bdywidth_c)

         DO jk = 1, nlev
            DO jc = i_startidx, i_endidx

               p_nh%diag%grf_tend_rho(jc,jk,jb) = rdt * (            &
                    &   latbc_data(tlev)%atm%rho(jc,jk,jb) &
                    & - latbc_data(prev_latbc_tlev)%atm%rho(jc,jk,jb) )

               p_nh%diag%grf_tend_thv(jc,jk,jb) = rdt * (                &
                    &   latbc_data(tlev)%atm%theta_v(jc,jk,jb) &
                    & - latbc_data(prev_latbc_tlev)%atm%theta_v(jc,jk,jb) )

               p_nh%diag%grf_tend_w(jc,jk,jb) = rdt * (            &
                    &   latbc_data(tlev)%atm%w(jc,jk,jb) &
                    & - latbc_data(prev_latbc_tlev)%atm%w(jc,jk,jb) )

            ENDDO
         ENDDO

         DO jc = i_startidx, i_endidx
            p_nh%diag%grf_tend_w(jc,nlevp1,jb) = rdt * (            &
                 &   latbc_data(tlev)%atm%w(jc,nlevp1,jb) &
                 & - latbc_data(prev_latbc_tlev)%atm%w(jc,nlevp1,jb) )
         ENDDO

         IF (ltransport) THEN
            DO jk = 1, nlev
               DO jc = i_startidx, i_endidx
                  p_nh%diag%grf_tend_tracer(jc,jk,jb,iqv) =  rdt * (   &
                       &   latbc_data(tlev)%atm%qv(jc,jk,jb) &
                       & - latbc_data(prev_latbc_tlev)%atm%qv(jc,jk,jb) )

                  p_nh%diag%grf_tend_tracer(jc,jk,jb,iqc) =  rdt * (   &
                       &   latbc_data(tlev)%atm%qc(jc,jk,jb) &
                       & - latbc_data(prev_latbc_tlev)%atm%qc(jc,jk,jb) )

                  p_nh%diag%grf_tend_tracer(jc,jk,jb,iqi) =  rdt * (   &
                       &   latbc_data(tlev)%atm%qi(jc,jk,jb) &
                       & - latbc_data(prev_latbc_tlev)%atm%qi(jc,jk,jb) )

                  p_nh%diag%grf_tend_tracer(jc,jk,jb,iqr) =  rdt * (   &
                       &   latbc_data(tlev)%atm%qr(jc,jk,jb) &
                       & - latbc_data(prev_latbc_tlev)%atm%qr(jc,jk,jb) )

                  p_nh%diag%grf_tend_tracer(jc,jk,jb,iqs) =  rdt * (   &
                       &   latbc_data(tlev)%atm%qs(jc,jk,jb) &
                       & - latbc_data(prev_latbc_tlev)%atm%qs(jc,jk,jb) )

               ENDDO
            ENDDO
         ENDIF

      ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif

    END SUBROUTINE compute_boundary_tendencies

    !-------------------------------------------------------------------------------------------------
    !> Send a message to the work PEs that the input prefetching PEs is ready. The
    !  counterpart on the work PEs side is compute_wait_for_async_pref
    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2013-03-19)
    !
    SUBROUTINE async_pref_send_handshake()
#ifndef NOMPI
      INTEGER :: msg
      ! Simply send a message from Input prefetching PE 0 to work PE 0
      ! p_pe_work == 0 signifies processor 0 in Input prefetching PEs
      IF(p_pe_work == 0) THEN
         !     CALL p_wait()

        msg = msg_pref_done
        CALL p_isend(msg, p_work_pe0, 0)
      ENDIF
#endif
    END SUBROUTINE async_pref_send_handshake

    !-------------------------------------------------------------------------------------------------
    !> compute_wait_for_async_pref: Wait for a message that the prefetch PE is ready.
    !  The counterpart on the input Prefetching PEs side is async_pref_send_handshake
    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2013-03-19)
    !
    SUBROUTINE compute_wait_for_async_pref()
#ifndef NOMPI
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::compute_wait_for_async_pref"
      INTEGER :: msg, action_tag

      ! First compute PE receives message from input prefetching leader
      IF(p_pe_work==0) THEN
         CALL p_recv(msg, p_pref_pe0, 0)
         action_tag = msg
         ! Just for safety: Check if we got the correct tag
         IF(action_tag /= msg_pref_done) CALL finish(routine, 'Compute PE: Got illegal prefetching tag')
      ENDIF
      ! Wait in barrier until message is here
      CALL p_barrier(p_comm_work)
#endif
    END SUBROUTINE compute_wait_for_async_pref


    !-------------------------------------------------------------------------------------------------
    !> async_pref_wait_for_start: Wait for a message from compute PE that we should start
    !  tranferring the prefetch data or finish. The counterpart on the compute side is
    !  compute_start_async_pref/compute_shutdown_async_pref
    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2013-03-19)
    !
    SUBROUTINE async_pref_wait_for_start(done)
      LOGICAL, INTENT(INOUT) :: done ! flag if we should shut down

#ifndef NOMPI
      INTEGER :: msg
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::async_pref_wait_for_start"

      ! Set output parameters to default values
      done  = .FALSE.

      ! Receive message that we may start transferring the prefetching data (or should finish)
      IF(p_pe_work == 0) THEN
         ! launch non-blocking receive request
         CALL p_wait()
         CALL p_irecv(msg, p_work_pe0, 0)
         CALL p_wait()
      END IF

      SELECT CASE(msg)
      CASE(msg_pref_start)

      CASE(msg_pref_shutdown)
         done = .TRUE.

      CASE DEFAULT
         ! Anything else is an error
         CALL finish(routine, 'Prefetching PE: Got illegal prefetching tag')
      END SELECT
#endif
    END SUBROUTINE async_pref_wait_for_start

    !-------------------------------------------------------------------------------------------------
    !> compute_start_async_pref: Send a message to prefetching PEs that they should start
    !  prefetching input. The counterpart on the prefetching side is async_pref_wait_for_start
    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2013-03-19)
    !
    SUBROUTINE compute_start_async_pref()
#ifndef NOMPI
      ! local variables
      INTEGER :: msg
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::compute_start_async_pref"

      !   CALL p_barrier(comm=p_comm_work) ! make sure all are here

      msg = msg_pref_start

      IF(p_pe_work==0) CALL p_send(msg, p_pref_pe0, 0)
#endif
    END SUBROUTINE compute_start_async_pref

    !-------------------------------------------------------------------------------------------------
    !> compute_shutdown_async_pref: Send a message to prefetching PEs that they should shut down
    !  The counterpart on the prefetching side is async_pref_wait_for_start
    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2013-03-19)
    !
    SUBROUTINE compute_shutdown_async_pref
#ifndef NOMPI
      INTEGER :: msg

      !  CALL p_barrier(comm=p_comm_work) ! make sure all are here

      msg = msg_pref_shutdown

      IF(p_pe_work==0) CALL p_send(msg, p_pref_pe0, 0)
#endif
    END SUBROUTINE compute_shutdown_async_pref

    !-------------------------------------------------------------------------
    !>
    ! Return the index for a give variable in mapped variable list
    !
    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2013-05-19)
    !!
    FUNCTION get_field_index(latbc,name,opt_lmap) RESULT(result_varID)
      TYPE (t_latbc_data), INTENT(IN) :: latbc
      CHARACTER (LEN=*),   INTENT(IN) :: name !< variable name
      LOGICAL, OPTIONAL,   INTENT(IN) :: opt_lmap
      ! local variables
      LOGICAL, PARAMETER :: ldebug = .FALSE.
      INTEGER :: result_varID, varID
      LOGICAL :: lmap

      lmap = .FALSE.
      IF (PRESENT(opt_lmap)) THEN
        IF (opt_lmap) lmap = .TRUE.
      ENDIF

      IF (ldebug)  WRITE (0,*) "name : ", TRIM(name)

      result_varID = -1
      ! looping over variable list in internal or mapped name
      IF (lmap) THEN
        DO varID=1, latbc%buffer%ngrp_vars
          IF (ldebug)  WRITE (0,*) "mapped name : ", TRIM(latbc%buffer%mapped_name(varID))
          IF (TRIM(name) == TRIM(latbc%buffer%mapped_name(varID))) THEN
            result_varID = varID
          END IF
        END DO
      ELSE
        DO varID=1, latbc%buffer%ngrp_vars
          IF (ldebug)  WRITE (0,*) "internal name : ", TRIM(latbc%buffer%internal_name(varID))
          IF (TRIM(name) == TRIM(latbc%buffer%internal_name(varID))) THEN
            result_varID = varID
          END IF
        END DO
      ENDIF

    END FUNCTION get_field_index

    !-------------------------------------------------------------------------
    !>
    !  Update linear interpolation coefficients for a given time stamp
    !
    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2014-08-15)
    !!
    SUBROUTINE update_lin_interpolation( latbc, step_datetime )
      TYPE(t_latbc_data), INTENT(IN)    :: latbc
      TYPE(datetime), pointer  :: step_datetime
      TYPE(datetime), pointer  :: mtime_step
      TYPE(timedelta), pointer :: delta_tstep

      CHARACTER(LEN=*), PARAMETER         :: routine = modname//"::update_lin_interpolation"

#ifndef NOMPI
      delta_tstep => newTimedelta("PT0S")
      delta_tstep = latbc%mtime_read - step_datetime

      IF(delta_tstep%month /= 0) &
           CALL finish(routine, "time difference for reading boundary data cannot be more than a month.")

      ! compute the number of "dtime_latbc" intervals fitting into the time difference "delta_tstep":

      latbc_config%lc1 = (delta_tstep%day * 86400._wp + delta_tstep%hour * 3600._wp +  &
           delta_tstep%minute * 60._wp + delta_tstep%second) / latbc_config%dtime_latbc
      latbc_config%lc2 = 1._wp - latbc_config%lc1

      ! deallocating mtime and deltatime
      CALL deallocateTimedelta(delta_tstep)
#endif

    END SUBROUTINE update_lin_interpolation

    !-------------------------------------------------------------------------

  END MODULE mo_async_latbc_utils
