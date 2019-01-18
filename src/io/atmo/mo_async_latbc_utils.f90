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
    USE mo_mpi,                 ONLY: my_process_is_pref, my_process_is_work,   &
         &                            p_comm_work, my_process_is_io,            &
         &                            my_process_is_mpi_test
    ! Processor numbers
    USE mo_mpi,                 ONLY: p_pref_pe0, p_pe_work, p_work_pe0, num_work_procs
    ! MPI Communication routines
    USE mo_mpi,                 ONLY: p_isend, p_barrier, p_wait, &
         &                            p_send, p_recv, p_bcast
    USE mo_latbc_read_recv,     ONLY: prefetch_cdi_2d, prefetch_cdi_3d, compute_data_receive
#endif

    USE mo_async_latbc_types,   ONLY: t_latbc_data, t_buffer
    USE mo_reorder_info,        ONLY: t_reorder_info
    USE mo_kind,                ONLY: wp, i8
    USE mo_util_string,         ONLY: int2string
    USE mo_parallel_config,     ONLY: nproma
    USE mo_model_domain,        ONLY: t_patch
    USE mo_grid_config,         ONLY: nroot
    USE mo_exception,           ONLY: message, finish, message_text
    USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, SUCCESS
    USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE
    USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
    USE mo_io_units,            ONLY: filename_max
    USE mo_nonhydro_types,      ONLY: t_nh_state
    USE mo_intp_data_strc,      ONLY: t_int_state
    USE mo_nh_vert_interp,      ONLY: vert_interp
    USE mo_physical_constants,  ONLY: cpd, rd, cvd_o_rd, p0ref, vtmpc1
    USE mo_util_phys,           ONLY: virtual_temp
    USE mo_nh_init_utils,       ONLY: interp_uv_2_vn, convert_thdvars, convert_omega2w, &
      &                               compute_input_pressure_and_height
    USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
    USE mtime,                  ONLY: timedelta, newTimedelta, deallocateTimedelta, &
         &                            newEvent, datetime, newDatetime,             &
         &                            isCurrentEventActive, deallocateDatetime,    &
         &                            MAX_DATETIME_STR_LEN,                        &
         &                            MAX_TIMEDELTA_STR_LEN, getPTStringFromMS,    &
         &                            getTotalSecondsTimedelta,                    &
         &                            datetimeToString, timedeltaToString,         &
         &                            OPERATOR(>=), OPERATOR(-), OPERATOR(>),      &
         &                            OPERATOR(/=), OPERATOR(+), OPERATOR(*),      &
         &                            OPERATOR(<), OPERATOR(==),                   &
         &                            getTriggerNextEventAtDateTime
    USE mo_util_mtime,          ONLY: dummyDateTime
    USE mo_time_config,         ONLY: time_config
    USE mo_limarea_config,      ONLY: latbc_config, generate_filename, LATBC_TYPE_EXT
    USE mo_initicon_config,     ONLY: timeshift
    USE mo_ext_data_types,      ONLY: t_external_data
    USE mo_run_config,          ONLY: iqv, iqc, iqi, iqr, iqs, ltransport, msg_level
    USE mo_dynamics_config,     ONLY: nnow, nnow_rcf
    USE mo_initicon_types,      ONLY: t_init_state
    USE mo_cdi,                 ONLY: streamOpenRead, streamClose, streamInqVlist, vlistInqTaxis, &
      &                               taxisInqVDate, taxisInqVTime, cdiDecodeTime, cdiDecodeDate
    USE mo_util_cdi,            ONLY: cdiGetStringError
    USE mo_master_config,       ONLY: isRestart
    USE mo_fortran_tools,       ONLY: copy, init
    USE mo_util_string,         ONLY: tolower
    IMPLICIT NONE

    PRIVATE

    ! handshake subroutines
    PUBLIC :: async_pref_send_handshake
    PUBLIC :: compute_wait_for_async_pref
    PUBLIC :: async_pref_wait_for_start
    PUBLIC :: compute_shutdown_async_pref
    PUBLIC :: allocate_pref_latbc_data

    PUBLIC ::  compute_init_latbc_data, async_init_latbc_data, read_latbc_data,     &
         &     update_lin_interpolation, recv_latbc_data


    INTERFACE fetch_from_buffer
      MODULE PROCEDURE fetch_from_buffer_2D
      MODULE PROCEDURE fetch_from_buffer_3D_cells 
      MODULE PROCEDURE fetch_from_buffer_3D_generic
    END INTERFACE



    !------------------------------------------------------------------------------------------------
    ! CONSTANTS
    !------------------------------------------------------------------------------------------------

    ! module name string
    CHARACTER(LEN=*), PARAMETER :: modname = 'mo_async_latbc_utils'

    ! Tags for communication between compute PEs and prefetching PEs
    INTEGER, PARAMETER :: msg_pref_start    = 56984
    INTEGER, PARAMETER :: msg_pref_done     = 26884
    INTEGER, PARAMETER :: msg_pref_shutdown = 48965

    INTEGER, PARAMETER :: TAG_PREFETCH2WORK = 2001
    INTEGER, PARAMETER :: TAG_WORK2PREFETCH = 2002
    INTEGER, PARAMETER :: TAG_VDATETIME     = 2003

  CONTAINS


    !-------------------------------------------------------------------------
    !>
    !! @par Revision History
    !! Initial version by S. Brdar, DWD (2013-06-13)
    !! Initial version by M.Pondkule, DWD (2014-01-27)
    !!
    SUBROUTINE allocate_pref_latbc_data(latbc, nlev_in, p_nh_state, ext_data, p_patch)

      TYPE(t_latbc_data), TARGET, INTENT(INOUT) :: latbc       !< variable buffer for latbc data
      INTEGER,                    INTENT(IN)    :: nlev_in     !< no. of vertical input levels
      TYPE(t_nh_state),           INTENT(INOUT) :: p_nh_state  !< nonhydrostatic state on the global domain
      TYPE(t_external_data),      INTENT(IN)    :: ext_data    !< external data on the global domain
      TYPE(t_patch),              INTENT(IN)    :: p_patch

#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::allocate_pref_latbc_data"
      INTEGER               :: tlev, nlev, nlevp1, nblks_c, nblks_e, ierrstat

      ! Allocate memory for variables (3D and 2D) on work processors
      nlev    = p_patch%nlev
      nlevp1  = p_patch%nlevp1
      nblks_c = p_patch%nblks_c
      nblks_e = p_patch%nblks_e

      ALLOCATE(latbc%latbc_data_const)
      ALLOCATE(latbc%latbc_data_const%z_mc_in(nproma,nlev_in,nblks_c), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      
      ! topography and metrics are time independent
      latbc%latbc_data_const%topography_c => ext_data%atm%topography_c
      latbc%latbc_data_const%z_ifc        => p_nh_state%metrics%z_ifc
      latbc%latbc_data_const%z_mc         => p_nh_state%metrics%z_mc

      DO tlev = 1, 2
        ! For safety; tke is checked for being associated in vert_interp
        NULLIFY(latbc%latbc_data(tlev)%atm_in%tke)

         ! Allocate atmospheric input data
         latbc%latbc_data(tlev)%atm_in%nlev = nlev_in
         !
         ALLOCATE( &
           latbc%latbc_data(tlev)%atm_in%pres (nproma,nlev_in,nblks_c), &
           latbc%latbc_data(tlev)%atm_in%temp (nproma,nlev_in,nblks_c), &
           latbc%latbc_data(tlev)%atm_in%u    (nproma,nlev_in,nblks_c), &
           latbc%latbc_data(tlev)%atm_in%v    (nproma,nlev_in,nblks_c), &
           latbc%latbc_data(tlev)%atm_in%w    (nproma,nlev_in,nblks_c), &
           latbc%latbc_data(tlev)%atm_in%qv   (nproma,nlev_in,nblks_c), &
           latbc%latbc_data(tlev)%atm_in%qc   (nproma,nlev_in,nblks_c), &
           latbc%latbc_data(tlev)%atm_in%qi   (nproma,nlev_in,nblks_c), &
           latbc%latbc_data(tlev)%atm_in%qr   (nproma,nlev_in,nblks_c), &
           latbc%latbc_data(tlev)%atm_in%qs   (nproma,nlev_in,nblks_c), STAT=ierrstat)
         IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

         ! allocate also vn (though sometimes not needed)
         ALLOCATE(latbc%latbc_data(tlev)%atm_in%vn(nproma, nlev_in, nblks_e), STAT=ierrstat)
         IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

         IF (latbc%buffer%lread_theta_rho) THEN
           ALLOCATE(latbc%latbc_data(tlev)%atm_in%rho(nproma,nlev_in,nblks_c),   &
                    latbc%latbc_data(tlev)%atm_in%theta_v(nproma,nlev_in,nblks_c), STAT=ierrstat)
           IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
         ENDIF

        ! initialize validity dates with dymmy values
        latbc%latbc_data(tlev)%vDateTime = dummyDateTime()

        ! initialize with zero to simplify implementation of sparse lateral boundary condition mode
!$OMP PARALLEL
        CALL init(latbc%latbc_data(tlev)%atm_in%pres)
        CALL init(latbc%latbc_data(tlev)%atm_in%temp)
        CALL init(latbc%latbc_data(tlev)%atm_in%u)
        CALL init(latbc%latbc_data(tlev)%atm_in%v)
        CALL init(latbc%latbc_data(tlev)%atm_in%w)
        CALL init(latbc%latbc_data(tlev)%atm_in%qv)
        CALL init(latbc%latbc_data(tlev)%atm_in%qc)
        CALL init(latbc%latbc_data(tlev)%atm_in%qi)
        CALL init(latbc%latbc_data(tlev)%atm_in%qr)
        CALL init(latbc%latbc_data(tlev)%atm_in%qs)
        CALL init(latbc%latbc_data(tlev)%atm_in%vn)
        IF (latbc%buffer%lread_theta_rho) THEN
          CALL init(latbc%latbc_data(tlev)%atm_in%rho)
          CALL init(latbc%latbc_data(tlev)%atm_in%theta_v)
        ENDIF
!$OMP END PARALLEL

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

         latbc%latbc_data(tlev)%const => latbc%latbc_data_const

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
      TYPE(t_patch),          INTENT(INOUT) :: p_patch
      TYPE(t_int_state),      INTENT(IN)    :: p_int_state
      TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state  !< nonhydrostatic state on the global domain
      INTEGER,                INTENT(OUT)   :: tlev

#ifndef NOMPI
      ! local variables
      TYPE(datetime) :: nextActive          ! next trigger date for prefetch event
      TYPE(datetime), POINTER :: latbc_read_datetime ! next input date to be read
      INTEGER        :: ierr
      REAL(wp)       :: seconds
      INTEGER        :: prev_latbc_tlev
      CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  :: td_string
      CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: latbc_read_datetime_str
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::compute_init_latbc_data"

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
           time_config%tc_exp_startdate, time_config%tc_stopdate, latbc%delta_dtime)

      latbc%mtime_last_read  => newDatetime(time_config%tc_current_date)
      latbc_read_datetime    => newDatetime(time_config%tc_current_date)


      ! if there is lrestart flag then prefetch processor needs to start reading
      ! the data from the new date time of restart file. Below the time at which
      ! reading of the lbc files starts (latbc_read_datetime) is calculated. Basically, 
      ! we query the prefetch Event to deduce the last trigger datetime. This is exactly 
      ! the time at which reading of the lbc files has to (re-)start.
      IF(isRestart()) THEN
        !
        CALL getTriggerNextEventAtDateTime(latbc%prefetchEvent,time_config%tc_current_date,nextActive,ierr)
        IF (nextActive > time_config%tc_current_date) THEN
          latbc_read_datetime = nextActive + latbc%delta_dtime*(-1._wp)
        ELSE
          latbc_read_datetime = nextActive
        ENDIF
      ENDIF

      ! prepare read/last indices
      tlev  = 1   ! read in the first time-level slot


      ! read z_ifc from lbc file at vv=0 (if available)
      CALL recv_latbc_data( latbc              = latbc,                      &
        &                   p_patch            = p_patch,                    &
        &                   p_nh_state         = p_nh_state,                 &
        &                   p_int              = p_int_state,                &
        &                   cur_datetime       = time_config%tc_current_date,&
        &                   latbc_read_datetime= time_config%tc_exp_startdate, &
        &                   lcheck_read        = .FALSE.,                    &
        &                   tlev               = tlev,                       &
        &                   opt_lconst_data_only = .TRUE.                    )


      ! read first two time steps
      IF (latbc_config%init_latbc_from_fg) THEN

        CALL message('','take lbc for initial time from fg')
        CALL datetimeToString(latbc_read_datetime, latbc_read_datetime_str)
        WRITE (message_text, '(a,a)')  "  read date: ", TRIM(latbc_read_datetime_str)
        CALL message('', message_text)

        prev_latbc_tlev = 3 - tlev
        tlev            = prev_latbc_tlev

        ! allocate input data for lateral boundary nudging
        CALL copy_fg_to_latbc(latbc%latbc_data, p_nh_state, tlev)
        ! set validity Datetime
        latbc%latbc_data(tlev)%vDateTime = time_config%tc_current_date

      ELSE

        CALL recv_latbc_data( latbc              = latbc,                      &
          &                   p_patch            = p_patch,                    &
          &                   p_nh_state         = p_nh_state,                 &
          &                   p_int              = p_int_state,                &
          &                   cur_datetime       = time_config%tc_current_date,&
          &                   latbc_read_datetime= latbc_read_datetime,        &
          &                   lcheck_read        = .FALSE.,                    &
          &                   tlev               = tlev                        )

      ENDIF
      ! (note: if we initialize from FG, then the following call is
      ! the first file read-in!)
      IF(.NOT. isRestart()) THEN
        IF (timeshift%dt_shift < 0 ) THEN
          ! For IAU, the second frame is always taken at tc_exp_startdate
          ! which is equivalent to latbc_read_datetime + timeshift%mtime_absshift
          latbc_read_datetime = time_config%tc_exp_startdate
        ELSE
          latbc_read_datetime = latbc_read_datetime + latbc%delta_dtime
        ENDIF
      ELSE
        latbc_read_datetime = latbc_read_datetime + latbc%delta_dtime
      ENDIF


      CALL recv_latbc_data( latbc              = latbc,                           &
        &                   p_patch            = p_patch,                         &
        &                   p_nh_state         = p_nh_state,                      &
        &                   p_int              = p_int_state,                     &
        &                   cur_datetime       = time_config%tc_current_date,     &
        &                   latbc_read_datetime= latbc_read_datetime,             &
        &                   lcheck_read        = .FALSE.,                         &
        &                   tlev               = tlev                             )


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
      TYPE(datetime) :: nextActive             ! next trigger date for prefetch event
      INTEGER        :: ierr
      TYPE(datetime), POINTER :: latbc_read_datetime    ! next input date to be read
      LOGICAL        :: done
      REAL(wp)       :: seconds
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
           time_config%tc_exp_startdate, time_config%tc_stopdate, latbc%delta_dtime)

      latbc%mtime_last_read  => newDatetime(time_config%tc_current_date)
      latbc_read_datetime    => newDatetime(time_config%tc_current_date)


      ! if there is lrestart flag then prefetch processor needs to start reading
      ! the data from the new date time of restart file. Below the time at which
      ! reading of the lbc files starts (latbc_read_datetime) is calculated. Basically, 
      ! we query the prefetch Event to deduce the last trigger datetime. This is exactly 
      ! the time at which reading of the lbc files has to (re-)start.
      IF(isRestart()) THEN
        !
        CALL getTriggerNextEventAtDateTime(latbc%prefetchEvent,time_config%tc_current_date,nextActive,ierr)
        IF (nextActive > time_config%tc_current_date) THEN
          latbc_read_datetime = nextActive + latbc%delta_dtime*(-1._wp)
        ELSE
          latbc_read_datetime = nextActive
        ENDIF
      ENDIF


      ! read z_ifc from lbc file with validity time vv=0 (if available)
      !
      CALL read_latbc_data(latbc                = latbc,                        &
        &                  latbc_read_datetime  = time_config%tc_exp_startdate, &
        &                  opt_lconst_data_only = .TRUE.                        )

      ! Inform compute PEs that we are done
      CALL async_pref_send_handshake()
      ! Wait for a message from the compute PEs to start
      CALL async_pref_wait_for_start(done)


      ! read first two time steps
      IF (.NOT. latbc_config%init_latbc_from_fg) THEN

        CALL read_latbc_data(latbc              = latbc,              &
          &                  latbc_read_datetime= latbc_read_datetime )

        ! Inform compute PEs that we are done
        CALL async_pref_send_handshake()
        ! Wait for a message from the compute PEs to start
        CALL async_pref_wait_for_start(done)

      ENDIF
      !
      IF(.NOT. isRestart()) THEN
        IF (timeshift%dt_shift < 0 ) THEN
          ! For IAU, the second frame is always taken at tc_exp_startdate
          ! which is equivalent to latbc_read_datetime + timeshift%mtime_absshift
          latbc_read_datetime = time_config%tc_exp_startdate
        ELSE
          latbc_read_datetime = latbc_read_datetime + latbc%delta_dtime
        ENDIF
      ELSE
        latbc_read_datetime = latbc_read_datetime + latbc%delta_dtime
      ENDIF
      ! (note: if we initialize from FG, then the following call is
      ! the first file read-in!)
      CALL read_latbc_data(latbc               = latbc,                &
        &                  latbc_read_datetime = latbc_read_datetime   )


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

    SUBROUTINE read_latbc_data(latbc, latbc_read_datetime, opt_lconst_data_only)

      TYPE(t_latbc_data), TARGET, INTENT(INOUT)    :: latbc

      ! datetime of the next input time level
      TYPE(datetime), INTENT(IN)         :: latbc_read_datetime

      ! read only constant fields from file and store in buffer
      LOGICAL, OPTIONAL :: opt_lconst_data_only

#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::read_latbc_data"

      INTEGER(KIND=MPI_ADDRESS_KIND)        :: ioff(0:num_work_procs-1)
      INTEGER                               :: jm, latbc_fileID, vlistID, taxisID, errno,  &
        &                                      idate, iyear, imonth, iday, itime, ihour,   &
        &                                      iminute, isecond, nlevs_read, nlevs
      LOGICAL                               :: l_exist
      CHARACTER(LEN=filename_max)           :: latbc_filename, latbc_full_filename
      CHARACTER(LEN=MAX_CHAR_LENGTH)        :: cdiErrorText
      CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: dstringA, dstringB
      TYPE(datetime), POINTER               :: mtime_vdate
      LOGICAL                               :: lconst_data_only
      TYPE(datetime), POINTER               :: vDateTime_ptr  ! pointer to vDateTime
      character(len=max_datetime_str_len)   :: vDateTime_str  ! vDateTime in String format


      lconst_data_only = .FALSE.
      IF (PRESENT(opt_lconst_data_only)) lconst_data_only = opt_lconst_data_only



      ! return if latbc_read_datetime is at least one full boundary data
      ! interval beyond the simulation end, implying that no further
      ! data are required for correct results
      IF (latbc_read_datetime >= time_config%tc_stopdate + latbc%delta_dtime) RETURN

      latbc_filename = generate_filename(nroot, latbc%patch_data%level, &
        &                                latbc_read_datetime, time_config%tc_exp_startdate)
      latbc_full_filename = TRIM(latbc_config%latbc_path)//TRIM(latbc_filename)

      IF (lconst_data_only) THEN
        WRITE(0,*) 'reading const data: ', TRIM(latbc_filename)
      ELSE
        WRITE(0,*) 'reading boundary data: ', TRIM(latbc_filename)
      ENDIF
      INQUIRE (FILE=TRIM(ADJUSTL(latbc_full_filename)), EXIST=l_exist)
      IF (.NOT. l_exist) THEN
         CALL finish(routine, "File not found: "//TRIM(latbc_filename))
      ENDIF

      ! opening and reading file
      latbc_fileID  = streamOpenRead(TRIM(latbc_full_filename))
      ! check if the file could be opened
      IF (latbc_fileID < 0) THEN
         CALL cdiGetStringError(latbc_fileID, cdiErrorText)
         CALL finish(routine, "File "//TRIM(latbc_full_filename)//" cannot be opened: "//TRIM(cdiErrorText))
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
      IF (mtime_vdate /= latbc_read_datetime) THEN
        CALL finish(routine, "requested date does not match file '"//TRIM(latbc_full_filename))
      END IF
      IF (msg_level >= 10) THEN
        CALL datetimeToString(latbc_read_datetime, dstringA)
        CALL datetimeToString(mtime_vdate, dstringB)
        WRITE (0,*) "  read date: ", TRIM(dstringA)
        WRITE (0,*) "  file date: ", TRIM(dstringB)
      END IF


      ! store validity datetime of current boundary data timeslice
      latbc%buffer%vDateTime = mtime_vdate
      ! send validity datetime of current boundary data timeslice from prefetch PE0 
      ! to worker PE0
      IF(p_pe_work == 0) THEN
        vDateTime_ptr => latbc%buffer%vDateTime
        CALL datetimeToString(vDateTime_ptr, vDateTime_str, errno)
        CALL p_isend(vDateTime_str, p_work_pe0, TAG_VDATETIME)
      ENDIF

      CALL deallocateDatetime(mtime_vdate)

      ! initializing the displacement array for each compute processor
      ioff(:) = 0_MPI_ADDRESS_KIND
      !
      ! Perform CDI read operation
      !
      IF (lconst_data_only) THEN

        ZIFCLOOP: DO jm = 1, latbc%buffer%ngrp_vars

          IF (TRIM(tolower(latbc%buffer%internal_name(jm))) /= TRIM(tolower(latbc%buffer%hhl_var))) THEN
            ! when skipping a variable: be sure to move the buffer
            ! offset position, too
            SELECT CASE (latbc%buffer%hgrid(jm))
            CASE(GRID_UNSTRUCTURED_CELL)
              ioff(:) = ioff(:) + INT(latbc%buffer%nlev(jm) * latbc%patch_data%cells%pe_own(:),i8)
            CASE(GRID_UNSTRUCTURED_EDGE)
              ioff(:) = ioff(:) + INT(latbc%buffer%nlev(jm) * latbc%patch_data%edges%pe_own(:),i8)
            CASE default
              CALL finish(routine,'unknown grid type')
            END SELECT

            CYCLE ZIFCLOOP
          ENDIF

          !WRITE (0,*) routine,'fetch variable '//TRIM(latbc%buffer%mapped_name(jm))


          IF(latbc%buffer%nlev(jm) /= 1 ) THEN
            ! Read 3d variables
            CALL prefetch_cdi_3d ( latbc_fileID, latbc%buffer%mapped_name(jm), latbc, &
              &                    nlevs_read, latbc%buffer%hgrid(jm), ioff )
          ENDIF

          ! consistency check
          IF (nlevs_read /= latbc%buffer%nlev(jm)) THEN
            WRITE (message_text, *) &
              & "variable '", TRIM(latbc%buffer%mapped_name(jm)), &
              & "': nlev=", nlevs_read, " but expected ", latbc%buffer%nlev(jm)
            CALL finish(routine, message_text)
          END IF

        ENDDO ZIFCLOOP

      ELSE

        VARLOOP: DO jm = 1, latbc%buffer%ngrp_vars

          ! ------------------------------------------------------------
          ! skip constant variables
          ! ------------------------------------------------------------
          IF (TRIM(tolower(latbc%buffer%internal_name(jm))) == TRIM(tolower(latbc%buffer%hhl_var))) THEN
            nlevs = latbc%buffer%nlev(jm)
            ! when skipping a variable: be sure to move the buffer
            ! offset position, too
            ioff(:) = ioff(:) + INT(nlevs * latbc%patch_data%cells%pe_own(:),i8)
            CYCLE VARLOOP
          END IF


          !WRITE (0,*) routine,'fetch variable '//TRIM(latbc%buffer%mapped_name(jm))

          ! Get pointer to appropriate reorder_info
          IF(latbc%buffer%nlev(jm) /= 1 ) THEN

            SELECT CASE (latbc%buffer%hgrid(jm))
            CASE(GRID_UNSTRUCTURED_CELL)
           
              ! Read 3d variables
              CALL prefetch_cdi_3d ( latbc_fileID, latbc%buffer%mapped_name(jm), latbc, &
                &                    nlevs_read, latbc%buffer%hgrid(jm), ioff )
            
            CASE(GRID_UNSTRUCTURED_EDGE)
              CALL prefetch_cdi_3d ( latbc_fileID, latbc%buffer%mapped_name(jm), latbc, &
                &                    nlevs_read, latbc%buffer%hgrid(jm), ioff )
            CASE default
              CALL finish(routine,'unknown grid type')
            END SELECT
          
            ! consistency check
            IF (nlevs_read /= latbc%buffer%nlev(jm)) THEN
              WRITE (message_text, *) &
                & "variable '", TRIM(latbc%buffer%mapped_name(jm)), &
                & "': nlev=", nlevs_read, " but expected ", latbc%buffer%nlev(jm)
              CALL finish(routine, message_text)
            END IF
          
          ELSE
            SELECT CASE (latbc%buffer%hgrid(jm))
            CASE(GRID_UNSTRUCTURED_CELL)
              ! Read 2d variables
              CALL prefetch_cdi_2d ( latbc_fileID, latbc%buffer%mapped_name(jm), latbc, &
                &                    latbc%buffer%hgrid(jm), ioff )
            
            CASE(GRID_UNSTRUCTURED_EDGE)
              CALL prefetch_cdi_2d ( latbc_fileID, latbc%buffer%mapped_name(jm), latbc, &
                &                    latbc%buffer%hgrid(jm), ioff )
            
            CASE default
              CALL finish(routine,'unknown grid type')
            END SELECT
          
          ENDIF

        ENDDO VARLOOP
      ENDIF  ! lconst_data_only

      ! close the open dataset file
      CALL streamClose(latbc_fileID)

      ! Store mtime_last_read
      latbc%mtime_last_read = latbc_read_datetime

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
    SUBROUTINE recv_latbc_data(latbc, p_patch, p_nh_state, p_int, cur_datetime, &
      &                        latbc_read_datetime, lcheck_read, tlev, opt_lconst_data_only)

      TYPE(t_latbc_data),     INTENT(INOUT) :: latbc
      TYPE(t_patch),          INTENT(INOUT) :: p_patch
      TYPE(t_nh_state),       INTENT(INOUT) :: p_nh_state   !< nonhydrostatic state on the global domain
      TYPE(t_int_state),      INTENT(IN)    :: p_int
      TYPE(datetime), POINTER, INTENT(IN)   :: cur_datetime !< current time
      INTEGER,                INTENT(INOUT) :: tlev

      ! datetime of the next input time level
      TYPE(datetime),         INTENT(IN)    :: latbc_read_datetime

      ! flag to set wether compute processor need to read boundary
      ! data or need not and than they will return
      LOGICAL,      INTENT(IN)    :: lcheck_read

      ! fetch only constant fields from buffer
      LOGICAL, OPTIONAL :: opt_lconst_data_only

#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::recv_latbc_data"
      LOGICAL                               :: isactive
      TYPE(timedelta), POINTER              :: my_duration_slack
      INTEGER                               :: prev_latbc_tlev
      character(len=max_datetime_str_len)   :: vDateTime_str
      TYPE(datetime), POINTER               :: vDateTime_ptr
      LOGICAL                               :: lconst_data_only


      lconst_data_only = .FALSE.
      IF (PRESENT(opt_lconst_data_only)) lconst_data_only = opt_lconst_data_only

      ! check for event been active
      my_duration_slack => newTimedelta("PT0S")
      my_duration_slack = time_config%tc_dt_model*0.5_wp
      isactive = isCurrentEventActive(latbc%prefetchEvent, cur_datetime, my_duration_slack)
      CALL deallocateTimedelta(my_duration_slack)

      ! do we need to read boundary data
      IF (lcheck_read .AND. (.NOT. isactive))  RETURN


      ! return if latbc_read_datetime is at least one full boundary data
      ! interval beyond the simulation end, implying that no further
      ! data are required for correct results
      IF (latbc_read_datetime >= time_config%tc_stopdate + latbc%delta_dtime) RETURN


      ! compute processors wait for msg from
      ! prefetch processor that they can start
      ! reading latbc data from memory window
      IF(my_process_is_work()) CALL compute_wait_for_async_pref()


      ! receive validity dateTime of current boundary data timeslice 
      ! from prefetch PE0.
      IF(p_pe_work==0) THEN
         CALL p_recv(vDateTime_str, p_pref_pe0, TAG_VDATETIME)
         vDateTime_ptr => newDatetime(vDateTime_str)
         latbc%buffer%vDateTime = vDateTime_ptr
         CALL deallocateDatetime(vDateTime_ptr)
      ENDIF
      CALL p_bcast(latbc%buffer%vDateTime,0,p_comm_work)


      IF (lconst_data_only) THEN

        ! start reading time-constant data
        CALL recv_const_data( latbc, p_patch, tlev )

      ELSE
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
        CALL compute_boundary_tendencies(latbc%latbc_data, p_patch, p_nh_state, tlev)

      ENDIF

      ! Store mtime_last_read
      latbc%mtime_last_read = latbc_read_datetime
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
      TYPE(t_patch),          INTENT(INOUT), TARGET :: p_patch
      TYPE(t_int_state),      INTENT(IN)            :: p_int
      INTEGER,                INTENT(IN)            :: tlev

#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::compute_latbc_icon_data"
      INTEGER(i8)                         :: eoff
      REAL(wp)                            :: temp_v(nproma,p_patch%nlev,p_patch%nblks_c)
      INTEGER                             :: jv
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
      IF (my_process_is_work()) CALL compute_start_async_pref()

      !
      ! get prognostic 3d fields
      !
      ! copying tha variable values from prefetch buffer to the respective allocated variable

      CALL fetch_from_buffer(latbc, 'temp', latbc%latbc_data(tlev)%atm%temp)
      CALL fetch_from_buffer(latbc, 'u',    latbc%latbc_data(tlev)%atm_in%u)
      CALL fetch_from_buffer(latbc, 'v',    latbc%latbc_data(tlev)%atm_in%v)
      CALL fetch_from_buffer(latbc, 'w',    latbc%latbc_data(tlev)%atm%w)

      ! Read parameter Pressure
      CALL fetch_from_buffer(latbc, 'pres', latbc%latbc_data(tlev)%atm%pres)

      ! Read parameter qv
      CALL fetch_from_buffer(latbc, 'qv', latbc%latbc_data(tlev)%atm%qv)
      ! Read parameter qc
      CALL fetch_from_buffer(latbc, 'qc', latbc%latbc_data(tlev)%atm%qc)
      ! Read parameter qi
      CALL fetch_from_buffer(latbc, 'qi', latbc%latbc_data(tlev)%atm%qi)
      ! Read parameter qr
      CALL fetch_from_buffer(latbc, 'qr', latbc%latbc_data(tlev)%atm%qr)
      ! Read parameter qs
      CALL fetch_from_buffer(latbc, 'qs', latbc%latbc_data(tlev)%atm%qs)

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
    !! NOTE: Unfortunately, this subroutine is MPI-collective, since
    !!       it contains several synchronization calls. It must
    !!       therefore be passed by all worker PEs. However, there is
    !!       the common situation where vertical interpolation shall
    !!       performed on a subset of PEs only, while no valid data is
    !!       available on the remaining PE. For this situation we need
    !!       the optional "opt_lmask" parameter.
    !!
    !! TODO: NOT adjusted to the mask parameter are the following
    !!       subroutines:
    !!
    !!     convert_omega2w
    !!     (operating in the column only)
    !!
    !!
    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2014-05-19)
    !!
    SUBROUTINE compute_latbc_intp_data( latbc, p_patch, p_nh_state, p_int, tlev)
      TYPE(t_latbc_data),     INTENT(INOUT), TARGET :: latbc  !< variable buffer for latbc data
      TYPE(t_patch),          INTENT(INOUT)         :: p_patch
      TYPE(t_nh_state),       INTENT(IN)            :: p_nh_state  !< nonhydrostatic state on the global domain
      TYPE(t_int_state),      INTENT(IN)            :: p_int
      INTEGER,                INTENT(IN)            :: tlev

#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::compute_latbc_intp_data"
      INTEGER(i8)                         :: eoff
      TYPE(t_reorder_info), POINTER       :: p_ri
      INTEGER                             :: jc, jk, jb, j, jv, nlev_in, nblks_c, ierrstat
      REAL(wp)                            :: log_exner, tempv
      REAL(wp), ALLOCATABLE               :: psfc(:,:), phi_sfc(:,:),    &
        &                                    w_ifc(:,:,:), omega(:,:,:)


      nblks_c = p_patch%nblks_c
      nlev_in = latbc%latbc_data(tlev)%atm_in%nlev

      ! consistency check
      IF ((tlev <1) .OR. (SIZE(latbc%latbc_data) < tlev)) THEN
        CALL finish(routine, "Internal error!")
      END IF

      ! Offset in memory window for async prefetching
      eoff = 0_i8

      DO jv = 1, latbc%buffer%ngrp_vars
        ! Receive 2d and 3d variables
        CALL compute_data_receive ( latbc%buffer%hgrid(jv), latbc%buffer%nlev(jv), &
          latbc%buffer%vars(jv)%buffer, eoff, latbc%patch_data)
      ENDDO

      ! Reading the next time step
      IF (my_process_is_work()) CALL compute_start_async_pref()

      p_ri => latbc%patch_data%cells


      ! get validity DateTime of current boundary data timeslice
      latbc%latbc_data(tlev)%vDateTime = latbc%buffer%vDateTime

      ! copying the variable values from prefetch buffer to the
      ! respective allocated variable

      ! Read parameter QV, QC and QI
      CALL fetch_from_buffer(latbc, 'qv', latbc%latbc_data(tlev)%atm_in%qv)
      CALL fetch_from_buffer(latbc, 'qc', latbc%latbc_data(tlev)%atm_in%qc)
      CALL fetch_from_buffer(latbc, 'qi', latbc%latbc_data(tlev)%atm_in%qi)


      ! Read parameter QR
      IF (latbc%buffer%lread_qr) THEN
        CALL fetch_from_buffer(latbc, 'qr', latbc%latbc_data(tlev)%atm_in%qr)
      ELSE
!$OMP PARALLEL
        CALL init(latbc%latbc_data(tlev)%atm_in%qr(:,:,:))
!$OMP END PARALLEL
      ENDIF


      ! Read parameter QS
      IF (latbc%buffer%lread_qs) THEN
        CALL fetch_from_buffer(latbc, 'qs', latbc%latbc_data(tlev)%atm_in%qs)
      ELSE
!$OMP PARALLEL
        CALL init(latbc%latbc_data(tlev)%atm_in%qs(:,:,:))
!$OMP END PARALLEL
      ENDIF

      IF (latbc%buffer%lread_theta_rho) THEN
        
        CALL fetch_from_buffer(latbc, 'theta_v', latbc%latbc_data(tlev)%atm_in%theta_v)
        CALL fetch_from_buffer(latbc, 'rho',     latbc%latbc_data(tlev)%atm_in%rho)

        ! Diagnose pres and temp from prognostic ICON variables
!$OMP PARALLEL DO PRIVATE (jk,j,jb,jc,log_exner,tempv)
#ifdef __LOOP_EXCHANGE
        DO j = 1, p_ri%n_own !p_patch%n_patch_cells
          jb = p_ri%own_blk(j) ! Block index in distributed patch
          jc = p_ri%own_idx(j) ! Line  index in distributed patch
!DIR$ IVDEP
          DO jk = 1, nlev_in
#else
        DO jk = 1, nlev_in
          DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jc = p_ri%own_idx(j) ! Line  index in distributed patch
#endif

            log_exner = (1._wp/cvd_o_rd)*LOG(latbc%latbc_data(tlev)%atm_in%rho(jc,jk,jb)* &
              latbc%latbc_data(tlev)%atm_in%theta_v(jc,jk,jb)*rd/p0ref)
            tempv = latbc%latbc_data(tlev)%atm_in%theta_v(jc,jk,jb)*EXP(log_exner)

            latbc%latbc_data(tlev)%atm_in%pres(jc,jk,jb) = p0ref*EXP(cpd/rd*log_exner)
            latbc%latbc_data(tlev)%atm_in%temp(jc,jk,jb) = &
              &    tempv / (1._wp + vtmpc1*latbc%latbc_data(tlev)%atm_in%qv(jc,jk,jb) - &
              &             (latbc%latbc_data(tlev)%atm_in%qc(jc,jk,jb) + &
              &              latbc%latbc_data(tlev)%atm_in%qi(jc,jk,jb) + &
              &              latbc%latbc_data(tlev)%atm_in%qr(jc,jk,jb) + &
              &              latbc%latbc_data(tlev)%atm_in%qs(jc,jk,jb)) )

          ENDDO
        ENDDO
!$OMP END PARALLEL DO

      END IF


      IF (latbc%buffer%lread_pres) THEN
        CALL fetch_from_buffer(latbc, 'pres', latbc%latbc_data(tlev)%atm_in%pres)
      ENDIF
      IF (latbc%buffer%lread_temp) THEN
        CALL fetch_from_buffer(latbc, 'temp', latbc%latbc_data(tlev)%atm_in%temp)
      ENDIF
      IF (latbc%buffer%lread_vn) THEN
        CALL fetch_from_buffer(latbc, 'vn', latbc%latbc_data(tlev)%atm_in%vn, latbc%patch_data%edges)
      END IF
      IF (latbc%buffer%lread_u_v) THEN
        CALL fetch_from_buffer(latbc, 'u', latbc%latbc_data(tlev)%atm_in%u)
        CALL fetch_from_buffer(latbc, 'v', latbc%latbc_data(tlev)%atm_in%v)
      ENDIF



      ! Read vertical component of velocity (W) or OMEGA

      IF (latbc%buffer%lconvert_omega2w) THEN

         ! allocate temporary array:
         ALLOCATE(omega(nproma, (nlev_in+1), nblks_c), STAT=ierrstat)
         IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

         IF (latbc%buffer%lread_w) THEN
           CALL fetch_from_buffer(latbc, 'w', omega)
         ELSE
!$OMP PARALLEL
           CALL init(omega(:,:,:))
!$OMP END PARALLEL
         ENDIF

      ELSE

         ! allocate temporary array:
         ALLOCATE(w_ifc(nproma,    (nlev_in+1), nblks_c), STAT=ierrstat)
         IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

         IF (latbc%buffer%lread_w) CALL fetch_from_buffer(latbc, 'w', w_ifc)

         ! Interpolate input 'w' from the interface levels to the main levels:

         ! consistency check
         IF ((SIZE(latbc%latbc_data(tlev)%atm_in%w,2) < nlev_in) .OR. (SIZE(w_ifc,2) < (nlev_in+1))) THEN
           CALL finish(routine, "Internal error!")
         END IF

         IF (latbc%buffer%lread_w) THEN
!$OMP PARALLEL DO PRIVATE (jk,j,jb,jc) ICON_OMP_DEFAULT_SCHEDULE
           DO jk = 1, nlev_in
             DO j = 1, p_ri%n_own
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jc = p_ri%own_idx(j) ! Line  index in distributed patch

               latbc%latbc_data(tlev)%atm_in%w(jc,jk,jb) = (w_ifc(jc,jk,jb) + w_ifc(jc,jk+1,jb)) * 0.5_wp
             ENDDO
           ENDDO
!$OMP END PARALLEL DO
         ELSE
!$OMP PARALLEL
           CALL init(latbc%latbc_data(tlev)%atm_in%w(:,:,:))
!$OMP END PARALLEL
         ENDIF
      ENDIF


      IF (latbc%buffer%lread_ps_geop) THEN
        ! allocate local temporary arrays:
        ALLOCATE(psfc(nproma, nblks_c), phi_sfc(nproma, nblks_c), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

        ! Read parameter surface pressure (LNPS)
        CALL fetch_from_buffer(latbc, "pres_sfc", psfc)

        ! Read parameter  surface Geopotential (GEOSP)
        CALL fetch_from_buffer(latbc, latbc%buffer%geop_ml_var, phi_sfc)
      ENDIF

      ! boundary exchange for a 2-D and 3-D array, needed because the
      ! vertical interpolation includes the halo region (otherwise, the
      ! syncs would have to be called after vert_interp)

      IF (latbc%buffer%lcompute_hhl_pres) THEN ! i.e. atmospheric data from IFS
        IF (latbc%patch_data%cells%n_own > 0) THEN
          CALL compute_input_pressure_and_height(p_patch, psfc, phi_sfc, &
               latbc%latbc_data(tlev), latbc%patch_data%cell_mask)
        END IF
      END IF

      IF (latbc%buffer%lconvert_omega2w) THEN
        ! (note that "convert_omega2w" requires the pressure field
        ! which has been computed before)
        IF (latbc%patch_data%cells%n_own > 0) THEN
          CALL convert_omega2w(omega, &
            &                  latbc%latbc_data(tlev)%atm_in%w,     &
            &                  latbc%latbc_data(tlev)%atm_in%pres,  &
            &                  latbc%latbc_data(tlev)%atm_in%temp,  &
            &                  nblks_c, p_patch%npromz_c, nlev_in )
        END IF

      END IF

      ! cleanup
      IF (ALLOCATED(omega))    DEALLOCATE(omega)
      IF (ALLOCATED(psfc))     DEALLOCATE(psfc)
      IF (ALLOCATED(phi_sfc))  DEALLOCATE(phi_sfc)
      IF (ALLOCATED(w_ifc))    DEALLOCATE(w_ifc)


      ! perform vertical interpolation of horizontally interpolated analysis data.
      !
      ! - note that we compute RHO in this subroutine.
      ! - note that "vert_interp" is MPI-collective, we cannot skip
      !   this for single PEs
      !
      IF (latbc_config%lsparse_latbc) THEN
        IF (latbc%patch_data%cells%n_own > 0 .OR. latbc%patch_data%edges%n_own > 0) THEN
          CALL vert_interp(p_patch, p_int, p_nh_state%metrics, latbc%latbc_data(tlev),   &
            &    opt_use_vn=latbc%buffer%lread_vn,                                       &
            &    opt_lmask_c=latbc%patch_data%cell_mask,                           &
            &    opt_lmask_e=latbc%patch_data%edge_mask, opt_latbcmode=.TRUE.)
        ENDIF
      ELSE
        CALL vert_interp(p_patch, p_int, p_nh_state%metrics, latbc%latbc_data(tlev),   &
          &    opt_use_vn=latbc%buffer%lread_vn, opt_latbcmode=.TRUE.)
      ENDIF

#endif
    END SUBROUTINE compute_latbc_intp_data


    !-------------------------------------------------------------------------
    !>
    !! Receive time-constant data from prefetch PEs
    !!
    !! This subroutine is called by compute processors.
    !! The following steps are performed:
    !! - Copy from memory window buffer, atmospheric time-constant data
    !! - halo synchronization
    !!
    !! NOTE: Unfortunately, this subroutine is MPI-collective, since
    !!       it contains several synchronization calls. It must
    !!       therefore be passed by all worker PEs. However, there is
    !!       the common situation where vertical interpolation shall
    !!       performed on a subset of PEs only, while no valid data is
    !!       available on the remaining PE. For this situation we need
    !!       the optional "opt_lmask" parameter.
    !!
    !!
    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2014-05-19)
    !! Excerpt from compute_latbc_intp_data, Daniel Reinert, DWD (2018-02-02)
    !!
    SUBROUTINE recv_const_data( latbc, p_patch, tlev)
      TYPE(t_latbc_data),     INTENT(INOUT), TARGET :: latbc  !< variable buffer for latbc data
      TYPE(t_patch),          INTENT(INOUT)         :: p_patch
      INTEGER,                INTENT(IN)            :: tlev

#ifndef NOMPI
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::compute_latbc_intp_data"
      INTEGER(i8)                         :: eoff
      INTEGER                             :: jc, jk, jb, j, jv, nlev_in, nblks_c, ierrstat
      REAL(wp), ALLOCATABLE               :: z_ifc_in(:,:,:)


      nblks_c = p_patch%nblks_c
      nlev_in = latbc%latbc_data(tlev)%atm_in%nlev

      ! consistency check
      IF ((tlev <1) .OR. (SIZE(latbc%latbc_data) < tlev)) THEN
        CALL finish(routine, "Internal error!")
      END IF

      ! Offset in memory window for async prefetching
      eoff = 0_i8

      DO jv = 1, latbc%buffer%ngrp_vars
        ! Receive 2d and 3d variables
        CALL compute_data_receive ( latbc%buffer%hgrid(jv), latbc%buffer%nlev(jv), &
          latbc%buffer%vars(jv)%buffer, eoff, latbc%patch_data)
      ENDDO

      ! Reading the next time step
      IF((.NOT. my_process_is_io()       .AND. &
           & .NOT. my_process_is_pref()) .AND. &
           & .NOT. my_process_is_mpi_test()) THEN
         CALL compute_start_async_pref()
      END IF


      ! copying the variable values from prefetch buffer to the
      ! respective allocated variable

      IF (latbc%buffer%lread_hhl) THEN
        
        ! allocate temporary array:
        ALLOCATE(z_ifc_in(nproma, (nlev_in+1), nblks_c), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) THEN
          WRITE (0,*) "nlev = ", nlev_in, "; nblks_c=", nblks_c
          CALL finish(routine, "ALLOCATE failed!")
        END IF

        CALL fetch_from_buffer(latbc, latbc%buffer%hhl_var, z_ifc_in)


!$OMP PARALLEL DO PRIVATE (jk,j,jb,jc) ICON_OMP_DEFAULT_SCHEDULE
        DO jk = 1, nlev_in
          DO j = 1, latbc%patch_data%cells%n_own
            jb = latbc%patch_data%cells%own_blk(j) ! Block index in distributed patch
            jc = latbc%patch_data%cells%own_idx(j) ! Line  index in distributed patch

            latbc%latbc_data_const%z_mc_in(jc,jk,jb)  = (z_ifc_in(jc,jk,jb) + z_ifc_in(jc,jk+1,jb) ) * 0.5_wp
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

        ! cleanup
        IF (ALLOCATED(z_ifc_in)) DEALLOCATE(z_ifc_in)
      END IF

#endif
    END SUBROUTINE recv_const_data



    ! Wrapper routine for copying prognostic variables from initial state to the 
    ! first time level of the lateral boundary data
    !
    SUBROUTINE copy_fg_to_latbc(latbc_data, p_nh, tlev)
      TYPE(t_init_state),     INTENT(INOUT) :: latbc_data(:)
      TYPE(t_nh_state),       INTENT(IN)    :: p_nh
      INTEGER,                INTENT(IN)    :: tlev

      INTEGER, PARAMETER :: jg = 1

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
      TYPE(t_init_state),     INTENT(IN)    :: latbc_data(:)
      TYPE(t_patch),          INTENT(IN)    :: p_patch
      TYPE(t_nh_state),       INTENT(INOUT) :: p_nh
      INTEGER,                INTENT(IN)    :: tlev

#ifndef NOMPI
      ! Local variables
      INTEGER                         :: i_startblk, i_endblk, i_startidx, i_endidx,    &
        &                                jc, jk, jb, je, nlev, nlevp1, prev_latbc_tlev
      TYPE(timedelta)                 :: td
      REAL(wp)                        :: dt, rdt
      CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: vDateTime_str_cur, vDateTime_str_prev
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::compute_boundary_tendencies"


      prev_latbc_tlev = 3 - tlev

      ! check if boundary tendency computation should be conducted or skipped
      ! this is decided upon the validity dates of the time levels involved.
      IF ((latbc_data(prev_latbc_tlev)%vDateTime == dummyDateTime()) .OR. &
          (latbc_data(tlev)%vDateTime == latbc_data(prev_latbc_tlev)%vDateTime) ) RETURN

      nlev            = p_patch%nlev
      nlevp1          = p_patch%nlevp1

      ! compute timedelta from validity times of consecutive boundary data timeslices.
      td = latbc_data(tlev)%vDateTime - latbc_data(prev_latbc_tlev)%vDateTime
      ! update frequency in s
      dt = REAL(getTotalSecondsTimedelta(td, latbc_data(tlev)%vDateTime))


      IF (msg_level >= 15) THEN
        CALL message(routine, "")
        CALL datetimeToString(latbc_data(tlev)%vDateTime, vDateTime_str_cur)
        CALL datetimeToString(latbc_data(prev_latbc_tlev)%vDateTime, vDateTime_str_prev)
        WRITE (message_text, '(a,a)')  "lbc vdate current : ", TRIM(vDateTime_str_cur)
        CALL message("", message_text)
        WRITE (message_text, '(a,a)')  "lbc vdate previous: ", TRIM(vDateTime_str_prev)
        CALL message("", message_text)
        WRITE (message_text, '(a,f9.2)')  "boundary update frequency: ", dt
        CALL message("", message_text)
      ENDIF
      ! Inverse value of boundary update frequency
      rdt = 1._wp/dt      

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
!DIR$ IVDEP
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
        msg = msg_pref_done
        CALL p_isend(msg, p_work_pe0, TAG_PREFETCH2WORK)
        CALL p_wait()
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
      INTEGER :: msg

      ! First compute PE receives message from input prefetching leader
      IF(p_pe_work==0) THEN
         msg = 0
         CALL p_recv(msg, p_pref_pe0, TAG_PREFETCH2WORK)
         ! Just for safety: Check if we got the correct tag
         IF(msg /= msg_pref_done) THEN
           CALL finish(routine, 'Compute PE: Got illegal prefetching tag: '//int2string(msg,'(i0)'))
         END IF
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
         CALL p_recv(msg, p_work_pe0, TAG_WORK2PREFETCH)
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
      IF(p_pe_work==0) CALL p_send(msg, p_pref_pe0, TAG_WORK2PREFETCH)
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
      IF(p_pe_work==0) CALL p_send(msg, p_pref_pe0, TAG_WORK2PREFETCH)
#endif
    END SUBROUTINE compute_shutdown_async_pref

    !-------------------------------------------------------------------------
    !>
    ! Return the index for a given variable in mapped variable list
    !
    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2013-05-19)
    !!
    FUNCTION get_field_index(buffer,name) RESULT(result_varID)
      TYPE(t_buffer), INTENT(IN) :: buffer
      CHARACTER (LEN=*),   INTENT(IN) :: name !< variable name
      ! local variables
      LOGICAL, PARAMETER :: ldebug = .FALSE.
      INTEGER :: result_varID, varID, nvars
      CHARACTER(len=len(name)) :: name_lc

      result_varID = -1
      nvars = buffer%ngrp_vars
      if (nvars >= 1) name_lc = tolower(name)
      ! looping over variable list in internal name
      LOOP : DO varID=1, nvars
        IF (name_lc == tolower(buffer%internal_name(varID))) THEN
          result_varID = varID
          EXIT LOOP
        END IF
      END DO LOOP
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
      TYPE(datetime), INTENT(in) :: step_datetime
#ifndef NOMPI
      TYPE(timedelta)           :: delta_tstep
      REAL(wp)                  :: dtime_latbc      ! time delta between two consecutive 
                                                    ! boundary forcing time slices
      TYPE(timedelta)           :: td
      LOGICAL :: failure

      CHARACTER(LEN=*), PARAMETER  :: routine = modname//"::update_lin_interpolation"
      CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: vDateTime_str_cur, vDateTime_str_prv
      CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN):: delta_tstep_str
      INTEGER :: prv_tlev, cur_tlev
      ! compute boundary update timedelta
      cur_tlev = latbc%new_latbc_tlev
      prv_tlev = latbc%prev_latbc_tlev()
      td = latbc%latbc_data(cur_tlev)%vDateTime - latbc%latbc_data(prv_tlev)%vDateTime
      dtime_latbc =  REAL(getTotalSecondsTimedelta(td, latbc%latbc_data(cur_tlev)%vDateTime))

      delta_tstep = latbc%latbc_data(cur_tlev)%vDateTime - step_datetime

      failure = delta_tstep%month /= 0
      IF (msg_level >= 15 .OR. failure) THEN
        CALL message(routine, "", all_print=failure)
        CALL datetimeToString(latbc%latbc_data(prv_tlev)%vDateTime, vDateTime_str_prv)
        CALL datetimeToString(latbc%latbc_data(cur_tlev)%vDateTime, vDateTime_str_cur)
        CALL timedeltaToString(delta_tstep, delta_tstep_str)
        WRITE (message_text, '(a,a)')  "lbc vdate current : ", vDateTime_str_cur
        CALL message("", message_text, all_print=failure)
        WRITE (message_text, '(a,a)')  "lbc vdate previous: ", vDateTime_str_prv
        CALL message("", message_text, all_print=failure)
        WRITE (message_text, '(a,a)')  "delta_tstep_str: ", delta_tstep_str
        CALL message("", message_text, all_print=failure)
        IF (delta_tstep%month /= 0) &
          CALL finish(routine, "time difference for reading boundary&
          & data must not be more than a month.")
      ENDIF




      ! compute the number of "dtime_latbc" intervals fitting into the time difference "delta_tstep":

      latbc_config%lc1 = (delta_tstep%day * 86400._wp + delta_tstep%hour * 3600._wp +  &
           delta_tstep%minute * 60._wp + delta_tstep%second) / dtime_latbc
      latbc_config%lc2 = 1._wp - latbc_config%lc1

      ! deallocating mtime and deltatime
#endif

    END SUBROUTINE update_lin_interpolation


    ! ----------------------------------------------------------------------
    ! Auxiliary routine: fetches data from latbc buffer.
    !
    SUBROUTINE fetch_from_buffer_2D(latbc, name, target_buf, opt_p_ri)
      TYPE(t_latbc_data), INTENT(IN), TARGET :: latbc             !< contains read buffer
      CHARACTER(LEN=*),   INTENT(IN)         :: name              !< variable name
      REAL(wp),           INTENT(INOUT)      :: target_buf(:,:)   !< target buffer
      TYPE(t_reorder_info), POINTER, OPTIONAL, INTENT(IN) :: opt_p_ri         !< patch indices (cells, edges)
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::fetch_from_buffer_2D"
      INTEGER :: jm, j, jb, jl
      TYPE(t_reorder_info), POINTER :: p_ri !< patch indices (cells, edges)

      p_ri => latbc%patch_data%cells
      IF (PRESENT(opt_p_ri))  p_ri => opt_p_ri

      jm = get_field_index(latbc%buffer, TRIM(name))
      IF (jm <= 0)  CALL finish(routine//"_"//TRIM(name), "Internal error, invalid field index!")

!$OMP PARALLEL DO PRIVATE (j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO j = 1, p_ri%n_own ! p_patch%n_patch_cells
        jb = p_ri%own_blk(j) ! Block index in distributed patch
        jl = p_ri%own_idx(j) ! Line  index in distributed patch
        target_buf(jl,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,1,jb), wp)
      ENDDO
!$OMP END PARALLEL DO
    END SUBROUTINE fetch_from_buffer_2D


    SUBROUTINE fetch_from_buffer_3D_cells(latbc, name, target_buf)
      TYPE(t_latbc_data), INTENT(IN), TARGET :: latbc             !< contains read buffer
      CHARACTER(LEN=*),   INTENT(IN)         :: name              !< variable name
      REAL(wp),           INTENT(INOUT)      :: target_buf(:,:,:) !< target buffer

      CALL fetch_from_buffer_3D_generic(latbc, name, target_buf, latbc%patch_data%cells)
    END SUBROUTINE fetch_from_buffer_3D_cells


    SUBROUTINE fetch_from_buffer_3D_generic(latbc, name, target_buf, p_ri)
      TYPE(t_latbc_data),   INTENT(IN), TARGET :: latbc             !< contains read buffer
      CHARACTER(LEN=*),     INTENT(IN)         :: name              !< variable name
      REAL(wp),             INTENT(INOUT)      :: target_buf(:,:,:) !< target buffer
      TYPE(t_reorder_info), INTENT(IN)         :: p_ri              !< patch indices (cells, edges)
      ! local variables
      CHARACTER(LEN=*), PARAMETER :: routine = modname//"::fetch_from_buffer_3D_generic"
      INTEGER :: jm, jk, j, jb, jl

      jm = get_field_index(latbc%buffer, TRIM(name))
      IF (jm <= 0)  CALL finish(routine//"_"//TRIM(name), "Internal error, invalid field index!")

      ! consistency check
      IF ((SIZE(target_buf,2) < latbc%buffer%nlev(jm)) .OR.  &
        & (SIZE(latbc%buffer%vars(jm)%buffer,2) < latbc%buffer%nlev(jm))) THEN
        CALL finish(routine//"_"//TRIM(name), "Internal error!")
      END IF

!$OMP PARALLEL DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
      DO j = 1, p_ri%n_own ! p_patch%n_patch_cells
        jb = p_ri%own_blk(j) ! Block index in distributed patch
        jl = p_ri%own_idx(j) ! Line  index in distributed patch
        DO jk=1, latbc%buffer%nlev(jm)
#else
      DO jk=1, latbc%buffer%nlev(jm)
        DO j = 1, p_ri%n_own ! p_patch%n_patch_cells
          jb = p_ri%own_blk(j) ! Block index in distributed patch
          jl = p_ri%own_idx(j) ! Line  index in distributed patch
#endif
          target_buf(jl,jk,jb) = REAL(latbc%buffer%vars(jm)%buffer(jl,jk,jb), wp)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    END SUBROUTINE fetch_from_buffer_3D_generic

    !-------------------------------------------------------------------------

  END MODULE mo_async_latbc_utils
