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
  !! @par Copyright
  !! 2002-2013 by DWD and MPI-M
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

    USE mo_async_latbc_types,   ONLY: t_patch_data, t_reorder_data, latbc_buffer ! for testing win_put
    USE mo_kind,                ONLY: wp, sp, i8
    USE mo_parallel_config,     ONLY: nproma
    USE mo_model_domain,        ONLY: t_patch
    USE mo_grid_config,         ONLY: nroot
    USE mo_exception,           ONLY: message, message_text, finish
    USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, MODE_COSMODE, MODE_DWDANA, MODE_ICONVREMAP, &
                                      MODE_IAU_OLD, MODE_IAU
    USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
    USE mo_io_units,            ONLY: filename_max
    USE mo_nonhydro_types,      ONLY: t_nh_state
    USE mo_intp_data_strc,      ONLY: t_int_state
    USE mo_nh_vert_interp,      ONLY: vert_interp
    USE mo_physical_constants,  ONLY: cpd, rd, cvd_o_rd, p0ref, vtmpc1
    USE mo_util_phys,           ONLY: virtual_temp
    USE mo_nh_init_utils,       ONLY: interp_uv_2_vn, convert_thdvars
    USE mo_sync,                ONLY: sync_patch_array, sync_patch_array_mult, SYNC_E, SYNC_C
    USE mo_initicon_types,      ONLY: t_initicon_state
    USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
    USE mtime,                  ONLY: event, newEvent, datetime, newDatetime,      &
         &                            isCurrentEventActive, deallocateDatetime,    &
         &                            MAX_DATETIME_STR_LEN, MAX_EVENTNAME_STR_LEN, &
         &                            MAX_TIMEDELTA_STR_LEN, getPTStringFromMS,    &
         &                            OPERATOR(>=), OPERATOR(-), OPERATOR(>)
    USE mo_mtime_extensions,    ONLY: get_datetime_string, get_duration_string_real
    USE mo_datetime,            ONLY: t_datetime
    USE mo_time_config,         ONLY: time_config
    USE mo_limarea_config,      ONLY: latbc_config, generate_filename_mtime
    USE mo_ext_data_types,      ONLY: t_external_data
    USE mo_run_config,          ONLY: iqv, iqc, iqi, iqr, iqs, ltransport, dtime, nsteps
    USE mo_initicon_config,     ONLY: init_mode
    USE mtime_events,           ONLY: deallocateEvent
    USE mtime_timedelta,        ONLY: timedelta, newTimedelta, deallocateTimedelta, &
         &                            operator(+)
    USE mo_cdi,                 ONLY: streamOpenRead, streamClose
    USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE
    USE mo_util_cdi,            ONLY: cdiGetStringError
    USE mo_master_config,       ONLY: isRestart
    USE mo_fortran_tools,       ONLY: copy, init
    IMPLICIT NONE

    PRIVATE

    ! constants :
    PUBLIC :: msg_pref_start
    PUBLIC :: msg_pref_done
    PUBLIC :: msg_pref_shutdown
    PUBLIC :: msg_latbc_done
    ! handshake subroutines
    PUBLIC :: async_pref_send_handshake
    PUBLIC :: compute_wait_for_async_pref
    PUBLIC :: async_pref_wait_for_start
    PUBLIC :: compute_start_async_pref
    PUBLIC :: compute_shutdown_async_pref

    PUBLIC ::  prepare_pref_latbc_data, pref_latbc_data, &
         &     latbc_data, latbc_fileID, new_latbc_tlev, prev_latbc_tlev,  &
         &     update_lin_interpolation, deallocate_pref_latbc_data, mtime_read

    !------------------------------------------------------------------------------------------------
    ! CONSTANTS
    !------------------------------------------------------------------------------------------------

    ! Tags for communication between compute PEs and prefetching PEs

    INTEGER, PARAMETER :: msg_pref_start    = 56984
    INTEGER, PARAMETER :: msg_pref_done     = 26884
    INTEGER, PARAMETER :: msg_pref_shutdown = 48965
    INTEGER, PARAMETER :: msg_latbc_done    = 20883
    CHARACTER(len=*), PARAMETER :: version = '$Id$'
    CHARACTER(LEN=*), PARAMETER :: modname = 'mo_async_latbc_utils'
    INTEGER                :: latbc_fileID, &
         new_latbc_tlev, &  ! time level indices for  latbc_data. can be 1 or 2.
         prev_latbc_tlev    ! new_latbc_tlev is the time level index carrying the most recent data
    TYPE(t_initicon_state) :: latbc_data(2)     ! storage for two time-level boundary data
    INTEGER                :: nlev_in             ! number of vertical levels in the boundary data
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: sim_start, sim_end, sim_cur, sim_cur_read
    CHARACTER(LEN=MAX_EVENTNAME_STR_LEN) :: event_name
    TYPE(timedelta), pointer :: my_duration_slack
    TYPE(datetime), pointer :: mtime_date
    TYPE(datetime), pointer :: mtime_read
    TYPE(datetime), pointer :: mtime_end
    TYPE(event), pointer :: prefetchEvent
    TYPE(timedelta), pointer :: delta_dtime
    LOGICAL :: isactive


  CONTAINS

    !--------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !>
    !! @par Revision History
    !! Initial version by S. Brdar, DWD (2013-06-13)
    !! Initial version by M.Pondkule, DWD (2014-01-27)
    !!
    SUBROUTINE allocate_pref_latbc_data(opt_p_nh_state, opt_ext_data,opt_p_patch)

      TYPE(t_nh_state), OPTIONAL,      INTENT(INOUT):: opt_p_nh_state  !< nonhydrostatic state on the global domain
      TYPE(t_external_data), OPTIONAL, INTENT(IN)   :: opt_ext_data    !< external data on the global domain
      TYPE(t_patch),OPTIONAL,  INTENT(IN)   :: opt_p_patch

#ifndef NOMPI
      ! local variables
      INTEGER       :: tlev, nlev, nlevp1, nblks_c, nblks_e
      INTEGER       :: ioper_mode

      CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = &
           "mo_async_latbc_utils::allocate_pref_latbc_data"

      !  CALL message(TRIM(routine),'start')

      nlev_in = latbc_config%nlev_in
      IF(nlev_in == 0) THEN
         CALL finish(TRIM(routine), "Number of input levels <nlev_in> not yet initialized.")
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

      ! Allocate memory for variables (3D and 2D) on work processors
      nlev    = opt_p_patch%nlev
      nlevp1  = opt_p_patch%nlevp1
      nblks_c = opt_p_patch%nblks_c
      nblks_e = opt_p_patch%nblks_e

      DO tlev = 1, 2
         ! Basic icon_remap data
         ALLOCATE(latbc_data(tlev)%topography_c(nproma,nblks_c),     &
              latbc_data(tlev)%z_ifc       (nproma,nlevp1,nblks_c),  &
              latbc_data(tlev)%z_mc        (nproma,nlev,nblks_c)   )

         ! Allocate atmospheric input data
         ALLOCATE(latbc_data(tlev)%atm_in%psfc(nproma,     nblks_c), &
              latbc_data(tlev)%atm_in%phi_sfc(nproma,      nblks_c), &
              latbc_data(tlev)%atm_in%pres (nproma,nlev_in,nblks_c), &
              latbc_data(tlev)%atm_in%z3d  (nproma,nlev_in,nblks_c), &
              latbc_data(tlev)%atm_in%temp (nproma,nlev_in,nblks_c), &
              latbc_data(tlev)%atm_in%u    (nproma,nlev_in,nblks_c), &
              latbc_data(tlev)%atm_in%v    (nproma,nlev_in,nblks_c), &
              latbc_data(tlev)%atm_in%w    (nproma,nlev_in,nblks_c), &
              latbc_data(tlev)%atm_in%omega(nproma,nlev_in,nblks_c), &
              latbc_data(tlev)%atm_in%qv   (nproma,nlev_in,nblks_c), &
              latbc_data(tlev)%atm_in%qc   (nproma,nlev_in,nblks_c), &
              latbc_data(tlev)%atm_in%qi   (nproma,nlev_in,nblks_c), &
              latbc_data(tlev)%atm_in%qr   (nproma,nlev_in,nblks_c), &
              latbc_data(tlev)%atm_in%qs   (nproma,nlev_in,nblks_c)  )

         ! allocate also vn (though sometimes not needed)
         ALLOCATE(latbc_data(tlev)%atm_in%vn(nproma, nlev_in, nblks_e))

         IF (ioper_mode >= 2) THEN
           ALLOCATE(latbc_data(tlev)%atm_in%w_ifc(nproma,nlev_in+1,nblks_c), &
                    latbc_data(tlev)%atm_in%z3d_ifc(nproma,nlev_in+1,nblks_c))
         ENDIF
         IF (ioper_mode == 3) THEN
           ALLOCATE(latbc_data(tlev)%atm_in%rho(nproma,nlev_in,nblks_c),   &
                    latbc_data(tlev)%atm_in%theta_v(nproma,nlev_in,nblks_c))
         ENDIF

         ! For safety; tke is checked for being associated in vert_interp
         NULLIFY(latbc_data(tlev)%atm_in%tke)

         ! Allocate atmospheric output data
         ALLOCATE(latbc_data(tlev)%atm%vn   (nproma,nlev,nblks_e), &
              latbc_data(tlev)%atm%u        (nproma,nlev,nblks_c), &
              latbc_data(tlev)%atm%v        (nproma,nlev,nblks_c), &
              latbc_data(tlev)%atm%w        (nproma,nlevp1,nblks_c), &
              latbc_data(tlev)%atm%temp     (nproma,nlev,nblks_c), &
              latbc_data(tlev)%atm%exner    (nproma,nlev,nblks_c), &
              latbc_data(tlev)%atm%pres     (nproma,nlev,nblks_c), &
              latbc_data(tlev)%atm%rho      (nproma,nlev,nblks_c), &
              latbc_data(tlev)%atm%theta_v  (nproma,nlev,nblks_c), &
              latbc_data(tlev)%atm%qv       (nproma,nlev,nblks_c), &
              latbc_data(tlev)%atm%qc       (nproma,nlev,nblks_c), &
              latbc_data(tlev)%atm%qi       (nproma,nlev,nblks_c), &
              latbc_data(tlev)%atm%qr       (nproma,nlev,nblks_c), &
              latbc_data(tlev)%atm%qs       (nproma,nlev,nblks_c)  )

         ! topography and metrics are time independent
!$OMP PARALLEL
         CALL copy(opt_ext_data%atm%topography_c(:,:), &
              latbc_data(tlev)%topography_c(:,:))
         CALL copy(opt_p_nh_state%metrics%z_ifc(:,:,:), &
              latbc_data(tlev)%z_ifc(:,:,:))
         CALL copy(opt_p_nh_state%metrics%z_mc (:,:,:), &
              latbc_data(tlev)%z_mc (:,:,:))
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

    SUBROUTINE prepare_pref_latbc_data(patch_data, p_patch, p_int_state, p_nh_state, ext_data)
      TYPE(t_patch_data),     INTENT(IN)   :: patch_data
      TYPE(t_patch),          OPTIONAL,INTENT(IN)   :: p_patch
      TYPE(t_int_state),      OPTIONAL,INTENT(IN)   :: p_int_state
      TYPE(t_nh_state),       OPTIONAL,INTENT(INOUT):: p_nh_state  !< nonhydrostatic state on the global domain
      TYPE(t_external_data),  OPTIONAL,INTENT(IN)   :: ext_data    !< external data on the global domain

#ifndef NOMPI
      ! local variables
      REAL(wp)     :: delta_tstep_secs, delta_dtime_secs, end_date
      TYPE(timedelta), pointer :: delta_tstep
      TYPE(timedelta), pointer :: delta_tend
      TYPE(datetime), pointer :: mtime_current
      TYPE(datetime), pointer :: mtime_finish
      LOGICAL       :: done
      INTEGER       :: i, add_delta, end_delta, finish_delta
      REAL(wp)      :: tdiff, seconds
      CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  :: tdiff_string
      CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = &
           "mo_async_latbc_utils::prepare_pref_latbc_data"
      CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: td_string

      !   CALL message(TRIM(routine),'start')

      IF( my_process_is_work() ) THEN
         ! allocate input data for lateral boundary nudging
         CALL allocate_pref_latbc_data( opt_p_nh_state=p_nh_state, opt_ext_data=ext_data, opt_p_patch=p_patch)
      ENDIF

      ! compute sim_start, sim_end in a formate appropriate for mtime
      CALL get_datetime_string(sim_start, time_config%ini_datetime)
      CALL get_datetime_string(sim_end, time_config%end_datetime)
      CALL get_datetime_string(sim_cur_read, time_config%cur_datetime)

      event_name = 'Prefetch input'

      ! convert namelist parameter "limarea_nml/dtime_latbc" into
      ! mtime object:
      IF (latbc_config%dtime_latbc > 86400._wp) THEN
        CALL finish(routine, "Namelist setting of limarea_nml/dtime_latbc too large for mtime conversion!")
      END IF
      seconds = latbc_config%dtime_latbc*1000._wp
      CALL getPTStringFromMS(NINT(seconds,i8), td_string)

      ! create prefetching event:
      prefetchEvent => newEvent(TRIM(event_name), TRIM(sim_start), &
           TRIM(sim_cur_read), TRIM(sim_end), td_string)

      tdiff = (0.5_wp*dtime)
      CALL get_duration_string_real(tdiff, tdiff_string)
      my_duration_slack => newTimedelta(tdiff_string)

      delta_dtime => newTimedelta(td_string)

      mtime_read  => newDatetime(TRIM(sim_start))

      ! time interval delta_dtime_secs in seconds
      delta_dtime_secs = 86400 *INT(delta_dtime%day)    &
           &                  + 3600  *INT(delta_dtime%hour)   &
           &                  + 60    *INT(delta_dtime%minute) &
           &                  +        INT(delta_dtime%second)

      ! a check so that prefetch processor checks whether to
      ! read the boundary data from a file
      IF(nsteps /= 0) THEN
         mtime_end => newDatetime(TRIM(sim_start))

         end_date = nsteps * dtime

         add_delta = FLOOR(end_date / delta_dtime_secs)

         DO i = 1, add_delta+1
            mtime_end = mtime_end + delta_dtime
         ENDDO

      ELSE
         mtime_finish => newDatetime(TRIM(sim_end))
         delta_tend => newTimedelta("PT0S")
         delta_tend = mtime_finish - mtime_read

         finish_delta = 86400 *INT(delta_tend%day)    &
              &                  + 3600  *INT(delta_tend%hour)   &
              &                  + 60    *INT(delta_tend%minute) &
              &                  +        INT(delta_tend%second)

         end_delta = FLOOR(finish_delta/delta_dtime_secs) * delta_dtime_secs
         ! check if the difference is less than zero than
         ! point mtime_end to time sim_end
         IF((finish_delta - end_delta) < 1e-15) THEN
            mtime_end => newDatetime(TRIM(sim_end))
         ELSE
            ! deallocating mtime and deltatime
            mtime_end => newDatetime(TRIM(sim_end))
            mtime_end = mtime_end + delta_dtime
         ENDIF
         CALL deallocateDatetime(mtime_finish)
         CALL deallocateTimedelta(delta_tend)
      ENDIF

      ! if there is lrestart flag than prefetch processor needs to start reading
      ! the data from the new date time of restart file. Below the time at which
      ! restart file starts is calculated and added to mtime_read which is the
      ! time step for reading the boundary data
      IF(isRestart()) THEN
         mtime_current => newDatetime(TRIM(sim_cur_read))
         delta_tstep => newTimedelta("PT0S")
         delta_tstep = mtime_read - mtime_current

         ! time interval delta_tstep_secs in seconds
         delta_tstep_secs = 86400 *INT(delta_tstep%day)    &
              &                  + 3600  *INT(delta_tstep%hour)   &
              &                  + 60    *INT(delta_tstep%minute) &
              &                  +        INT(delta_tstep%second)

         add_delta = FLOOR(delta_tstep_secs / delta_dtime_secs)
         ! no of times delta_dtime needs to be added to get the current time to
         ! read the boundary data file
         DO i = 1, add_delta
            mtime_read = mtime_read + delta_dtime
         ENDDO
         ! deallocating mtime and deltatime
         CALL deallocateTimedelta(delta_tstep)
         CALL deallocateDatetime(mtime_current)
      ENDIF

      ! prepare read/last indices
      new_latbc_tlev = 1   ! read in the first time-level slot
      prev_latbc_tlev = 2

      ! read first two time steps
      IF( my_process_is_work()) THEN  ! IF (PRESENT(p_patch)) THEN
         CALL pref_latbc_data( patch_data, p_patch, p_nh_state, p_int_state, &
              &                   time_config%cur_datetime, lopt_check_read=.FALSE., lopt_time_incr=.FALSE.)
      ELSE IF( my_process_is_pref()) THEN
         CALL pref_latbc_data( patch_data, datetime=time_config%cur_datetime, &
              &                   lopt_check_read=.FALSE., lopt_time_incr=.FALSE. )
         ! Inform compute PEs that we are done
         CALL async_pref_send_handshake()
         ! Wait for a message from the compute PEs to start
         CALL async_pref_wait_for_start(done)
      END IF

      IF( my_process_is_work()) THEN  !IF (PRESENT(p_patch)) THEN
         CALL pref_latbc_data( patch_data, p_patch, p_nh_state, p_int_state, &
              &                   time_config%cur_datetime, lopt_check_read=.FALSE., lopt_time_incr=.TRUE.)
      ELSE IF( my_process_is_pref()) THEN
         CALL pref_latbc_data( patch_data,datetime=time_config%cur_datetime, &
              &                   lopt_check_read=.FALSE., lopt_time_incr=.TRUE.)
         ! Inform compute PEs that we are done
         CALL async_pref_send_handshake()
      END IF

      CALL message(TRIM(routine),'done')
#endif
    END SUBROUTINE prepare_pref_latbc_data

    !-------------------------------------------------------------------------
    !>
    !! Read horizontally interpolated atmospheric boundary data
    !!
    !! The subroutine reads atmospheric boundary data and projects on
    !! the ICON global grid
    !!
    !! @par Revision History
    !! Initial version by S. Brdar, DWD (2013-07-19)
    !! Modified version by M. Pondkule, DWD (2014-02-11)
    !!

    SUBROUTINE pref_latbc_data( patch_data, p_patch, p_nh_state, p_int, datetime,&
         &                   lopt_check_read, lopt_time_incr)
      TYPE(t_patch_data),     INTENT(IN), TARGET     :: patch_data
      TYPE(t_patch),          OPTIONAL,INTENT(IN)    :: p_patch
      TYPE(t_nh_state),       OPTIONAL,INTENT(INOUT) :: p_nh_state  !< nonhydrostatic state on the global domain
      TYPE(t_int_state),      OPTIONAL,INTENT(IN)    :: p_int
      TYPE(t_datetime),       OPTIONAL,INTENT(INOUT) :: datetime       !< current time
      LOGICAL,      INTENT(IN), OPTIONAL    :: lopt_check_read
      LOGICAL,      INTENT(IN), OPTIONAL    :: lopt_time_incr  !< increment latbc_datetime

#ifndef NOMPI
      ! local variables
      LOGICAL                               :: lcheck_read
      LOGICAL                               :: ltime_incr
      CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = "mo_async_latbc_utils::pref_latbc_data"

      IF( my_process_is_work()) THEN
         ! compute current datetime in a format appropriate for mtime
         CALL get_datetime_string(sim_cur, datetime) ! time_config%cur_datetime)
         mtime_date  => newDatetime(TRIM(sim_cur))

         ! check for event been active
         isactive = isCurrentEventActive(prefetchEvent, mtime_date, my_duration_slack)
      ENDIF

      ! flag to set wether compute processor need to read boundary
      ! data or need not and than they will return
      IF (PRESENT(lopt_check_read)) THEN
         lcheck_read = lopt_check_read
      ELSE
         lcheck_read = .TRUE.
      ENDIF

      ! flag to increment present datetime to next time level
      IF (PRESENT(lopt_time_incr)) THEN
         ltime_incr = lopt_time_incr
      ELSE
         ltime_incr = .TRUE.
      ENDIF

      !
      ! do we need to read boundary data
      !
      IF (.NOT. my_process_is_pref()) THEN
         IF ((lcheck_read .AND. (.NOT. isactive))) THEN
            RETURN
         ENDIF
      ENDIF

      ! compute processors wait for msg from
      ! prefetch processor that they can start
      ! reading latbc data from memory window
      IF((.NOT. my_process_is_io() .AND. &
           & .NOT. my_process_is_pref()) .AND. &
           & .NOT. my_process_is_mpi_test()) THEN
         CALL compute_wait_for_async_pref()
      END IF

      ! Prepare the mtime_read for the next time level
      IF(ltime_incr ) THEN
         mtime_read = mtime_read + delta_dtime
      ENDIF

      ! Adjust read/last indices
      !
      ! New boundary data time-level is always read in latbc_data(new_latbc_tlev),
      ! whereas latbc_data(prev_latbc_tlev) always holds the last read boundary data
      !
      new_latbc_tlev = prev_latbc_tlev
      prev_latbc_tlev = 3 - new_latbc_tlev

      !
      ! start reading boundary data
      !
      IF (latbc_config%itype_latbc == 1) THEN
         IF( my_process_is_pref()) THEN
            CALL prefetch_latbc_intp_data( patch_data )
         ELSE IF( my_process_is_work()) THEN
            CALL compute_latbc_intp_data( p_patch, patch_data, p_nh_state, p_int )
            ! NOMPI
            ! Compute tendencies for nest boundary update
            IF (ltime_incr) CALL compute_boundary_tendencies(p_patch, p_nh_state)
         ENDIF
      ELSE
         IF( my_process_is_pref()) THEN
            CALL prefetch_latbc_icon_data( patch_data )
         ELSE IF( my_process_is_work()) THEN
            CALL compute_latbc_icon_data( p_patch, patch_data, p_int )
            ! NOMPI
            ! Compute tendencies for nest boundary update
            IF (ltime_incr) CALL compute_boundary_tendencies(p_patch, p_nh_state)
         ENDIF
      ENDIF
#endif
    END SUBROUTINE pref_latbc_data

    !-------------------------------------------------------------------------
    !>
    !! Read atmospheric ICON output provided on the same vertical (and horizontal) grid
    !! as used in the limited-area run. This in practice means that the data need to come
    !! from a nested global or larger-scale ICON run whose nested domain is identical to
    !! the limited-area domain considered here. Thus, the functionality of this routine is
    !! restricted to a very specific test configuration.
    !!
    !! This subroutine is called by prefetch processor.
    !! The following steps are performed:
    !! - read atmospheric input data,
    !! - Write input data to memory window buffer. The offset for data
    !!   is set such that each of dataset belongs to the respective compute processor,

    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2014-05-07)
    !!
    SUBROUTINE prefetch_latbc_icon_data( patch_data )
      TYPE(t_patch_data), INTENT(IN)      :: patch_data

#ifndef NOMPI
      INTEGER(KIND=MPI_ADDRESS_KIND)      :: ioff(0:num_work_procs-1)
      ! local variables
      INTEGER                             :: jm, latbc_fileID
      LOGICAL                             :: l_exist
      CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = "mo_async_latbc_utils::prefetch_latbc_icon_data"
      CHARACTER(LEN=filename_max)           :: latbc_filename, latbc_full_filename
      CHARACTER(len=132)             :: message_text
      CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText

      ! if mtime_read is same as mtime_end the prefetch processor returns without further
      ! proceeding to generate filename and than looking for boundary data file
      IF(mtime_read >= mtime_end) &
           RETURN
      latbc_filename = generate_filename_mtime(nroot, patch_data%level, mtime_read)
      latbc_full_filename = TRIM(latbc_config%latbc_path)//TRIM(latbc_filename)
      WRITE(0,*) 'reading boundary data: ', TRIM(latbc_filename)
      INQUIRE (FILE=TRIM(ADJUSTL(latbc_full_filename)), EXIST=l_exist)
      IF (.NOT. l_exist) THEN
         WRITE (message_text,'(a,a)') 'file not found:', TRIM(latbc_filename)
         CALL finish(TRIM(routine), message_text)
      ENDIF
      !
      ! open file
      !
      latbc_fileID  = streamOpenRead(TRIM(latbc_full_filename))
      ! check if the file could be opened
      IF (latbc_fileID < 0) THEN
         CALL cdiGetStringError(latbc_fileID, cdiErrorText)
         WRITE(message_text,'(4a)') 'File ', TRIM(latbc_full_filename), &
              ' cannot be opened: ', TRIM(cdiErrorText)
         CALL finish(routine, TRIM(message_text))
      ENDIF

      ! initializing the displacement array for each compute processor
      ioff(:) = 0_MPI_ADDRESS_KIND

      !
      ! Prefetch ICON data using CDI read
      ! read prognostic 3d fields
      !
      DO jm = 1, latbc_buffer%ngrp_vars
         ! IF(jm == 1) THEN ! testing values for temperature
         IF(latbc_buffer%nlev(jm) /= 1 ) THEN
            SELECT CASE (latbc_buffer%hgrid(jm))
            CASE(GRID_UNSTRUCTURED_CELL)
               ! Read 3d variables
               CALL prefetch_cdi_3d ( latbc_fileID, latbc_buffer%mapped_name(jm), patch_data, &
                    &                 latbc_buffer%nlev(jm), latbc_buffer%hgrid(jm), ioff )
            CASE(GRID_UNSTRUCTURED_EDGE)
               CALL prefetch_cdi_3d ( latbc_fileID, latbc_buffer%mapped_name(jm), patch_data, &
                    &                 latbc_buffer%nlev(jm), latbc_buffer%hgrid(jm), ioff )
            CASE default
               CALL finish(routine,'unknown grid type')
            END SELECT
         ELSE
            SELECT CASE (latbc_buffer%hgrid(jm))
            CASE(GRID_UNSTRUCTURED_CELL)
               ! Read 2d variables
               CALL prefetch_cdi_2d ( latbc_fileID, latbc_buffer%mapped_name(jm), patch_data, &
                    &                 latbc_buffer%hgrid(jm), ioff )
            CASE(GRID_UNSTRUCTURED_EDGE)
               CALL prefetch_cdi_2d ( latbc_fileID, latbc_buffer%mapped_name(jm), patch_data, &
                    &                 latbc_buffer%hgrid(jm), ioff )
            CASE default
               CALL finish(routine,'unknown grid type')
            END SELECT
         ENDIF
      ENDDO

      !
      ! close the open dataset file
      !
      CALL streamClose(latbc_fileID)

#endif
    END SUBROUTINE prefetch_latbc_icon_data

    !-------------------------------------------------------------------------
    !>
    !! Copy data read by the prefetch_latbc_icon_data routine
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
    SUBROUTINE compute_latbc_icon_data( p_patch, patch_data, p_int )
      TYPE(t_patch),          TARGET      :: p_patch
      TYPE(t_patch_data),     TARGET, INTENT(IN) :: patch_data
      TYPE(t_int_state),      INTENT(IN)  :: p_int

#ifndef NOMPI
      ! local variables
      INTEGER(i8)                         :: eoff
      TYPE(t_reorder_data), POINTER       :: p_ri
      REAL(wp)                            :: temp_v(nproma,p_patch%nlev,p_patch%nblks_c)
      INTEGER                             :: jc, jk, jb, jm, tlev, j, jl, jv
      CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = "mo_async_latbc_utils::compute_latbc_icon_data"
      !-------------------------------------------------------------------------

      nlev_in   = latbc_config%nlev_in
      tlev      = new_latbc_tlev

      ! Offset in memory window for async prefetching
      eoff = 0_i8
      DO jv = 1, latbc_buffer%ngrp_vars
         ! Receive 2d and 3d variables
         CALL compute_data_receive (latbc_buffer%hgrid(jv), latbc_buffer%nlev(jv), &
                                    latbc_buffer%vars(jv)%buffer, eoff, patch_data)
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
      p_ri => patch_data%cells

      !
      ! get prognostic 3d fields
      !
      ! copying tha variable values from prefetch buffer to the respective allocated variable

!$OMP PARALLEL PRIVATE(jm,jv,jc)
      jm = get_field_index('temp')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc_buffer%nlev(jm)
         DO j = 1, p_ri%n_own ! p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            !       WRITE(0,*) 'compute_latbc_icon_data 01 ',latbc_data(tlev)%atm%temp(jl,jk,jb)
            latbc_data(tlev)%atm%temp(jl,jk,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO  !NOWAIT

      jm = get_field_index('u')
      jv = get_field_index('v')
      !   WRITE(0,*)'pref_latbc_cdi_data name ', latbc_buffer%mapped_name(jm), ' jm ', jm
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc_buffer%nlev(jm)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm_in%u(jl,jk,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,jk,jb), wp)
            latbc_data(tlev)%atm_in%v(jl,jk,jb) = REAL(latbc_buffer%vars(jv)%buffer(jl,jk,jb), wp)
            !          IF( p_pe_work == p_work_pe0 ) &
            !           &  WRITE(0,*) 'latbc_data(tlev)%atm_in%v value ', latbc_data(tlev)%atm_in%v(jl,jk,jb)
         ENDDO
      ENDDO
!$OMP END DO

      jv = get_field_index('w')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc_buffer%nlev(jv)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm%w(jl,jk,jb) = REAL(latbc_buffer%vars(jv)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      ! Read parameter Pressure
      jv = get_field_index('pres')
      !      WRITE(0,*)'pref_latbc_cdi_data jv ', latbc_buffer%mapped_name(jv), ' jv ', jv
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc_buffer%nlev(jv)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm%pres(jl,jk,jb) = REAL(latbc_buffer%vars(jv)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      ! Read parameter qv
      jv = get_field_index('qv')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc_buffer%nlev(jv)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm%qv(jl,jk,jb) = REAL(latbc_buffer%vars(jv)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      jc = get_field_index('qc')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc_buffer%nlev(jc)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm%qc(jl,jk,jb) = REAL(latbc_buffer%vars(jc)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      ! Read parameter qi
      jm = get_field_index('qi')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc_buffer%nlev(jm)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm%qi(jl,jk,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      ! Read parameter qr
      jm = get_field_index('qr')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc_buffer%nlev(jm)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm%qr(jl,jk,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      ! Read parameter qs
      jv = get_field_index('qs')
      !    WRITE(0,*)'pref_latbc_cdi_data jv ', latbc_buffer%mapped_name(jv), ' jv ', jv
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc_buffer%nlev(jv)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm%qs(jl,jk,jb) = REAL(latbc_buffer%vars(jv)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      ! boundary exchange for a 2-D and 3-D array to fill HALO region.
      ! This addition by M.Pondkule, DWD (11/06/2014)
      CALL sync_patch_array(SYNC_C, p_patch, latbc_data(tlev)%atm%temp)
      CALL sync_patch_array(SYNC_C, p_patch, latbc_data(tlev)%atm_in%u)
      CALL sync_patch_array(SYNC_C, p_patch, latbc_data(tlev)%atm_in%v)
      CALL sync_patch_array(SYNC_C, p_patch, latbc_data(tlev)%atm%w)
      CALL sync_patch_array(SYNC_C, p_patch, latbc_data(tlev)%atm%pres)
      CALL sync_patch_array(SYNC_C, p_patch, latbc_data(tlev)%atm%qv)
      CALL sync_patch_array(SYNC_C, p_patch, latbc_data(tlev)%atm%qc)
      CALL sync_patch_array(SYNC_C, p_patch, latbc_data(tlev)%atm%qi)
      CALL sync_patch_array(SYNC_C, p_patch, latbc_data(tlev)%atm%qr)
      CALL sync_patch_array(SYNC_C, p_patch, latbc_data(tlev)%atm%qs)
      CALL sync_patch_array(SYNC_C, p_patch, latbc_data(tlev)%atm%rho)
      CALL sync_patch_array(SYNC_C, p_patch, latbc_data(tlev)%atm%theta_v)
      !
      ! Convert u and v on cell points to vn at edge points
      !
      CALL interp_uv_2_vn( p_patch, p_int, latbc_data(tlev)%atm_in%u,              &
           &                  latbc_data(tlev)%atm_in%v, latbc_data(tlev)%atm%vn )

      !
      ! Compute virtual temperature
      !
      CALL virtual_temp( p_patch, latbc_data(tlev)%atm%temp, latbc_data(tlev)%atm%qv, &
           &                latbc_data(tlev)%atm%qc, latbc_data(tlev)%atm%qi,               &
           &                latbc_data(tlev)%atm%qr, latbc_data(tlev)%atm%qs,               &
           &                temp_v=temp_v )

      !
      ! Compute NH prognostic thermodynamical variables
      !
      CALL convert_thdvars( p_patch, latbc_data(tlev)%atm%pres, temp_v,              &
           &                   latbc_data(tlev)%atm%rho,                                   &
           &                   latbc_data(tlev)%atm%exner,                                 &
           &                   latbc_data(tlev)%atm%theta_v )

      CALL sync_patch_array(SYNC_E, p_patch, latbc_data(tlev)%atm%vn)

#endif
    END SUBROUTINE compute_latbc_icon_data
    !-------------------------------------------------------------------------


    !-------------------------------------------------------------------------
    !>
    !! Read horizontally interpolated atmospheric analysis or forecast data
    !!
    !! This subroutine is called by prefetch processor.
    !! The following steps are performed:
    !! - read atmospheric input data,
    !! - Write input data to memory window buffer. The offset for data
    !!   is set such that each of dataset belongs to the respective compute processor,
    !!
    !! @par Revision History
    !! Initial version by M. Pondkule, DWD (2014-04-25)
    !!
    SUBROUTINE prefetch_latbc_intp_data( patch_data )
      TYPE(t_patch_data),     INTENT(IN)  :: patch_data

#ifndef NOMPI
      ! local variables
      INTEGER(KIND=MPI_ADDRESS_KIND)      :: ioff(0:num_work_procs-1)
      INTEGER                             :: jm, latbc_fileID
      LOGICAL                             :: l_exist
      CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = "mo_async_latbc_utils::pref_latbc_intp_data"
      CHARACTER(LEN=filename_max)           :: latbc_filename, latbc_full_filename
      CHARACTER(LEN=132)             :: message_text
      CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText

      ! if mtime_read is same as mtime_end the prefetch processor returns without further
      ! proceeding to generate filename and than looking for boundary data file
      IF(mtime_read >= mtime_end) &
           RETURN
      latbc_filename = generate_filename_mtime(nroot, patch_data%level, mtime_read)
      latbc_full_filename = TRIM(latbc_config%latbc_path)//TRIM(latbc_filename)
      WRITE(0,*) 'reading boundary data: ', TRIM(latbc_filename)
      INQUIRE (FILE=TRIM(ADJUSTL(latbc_full_filename)), EXIST=l_exist)
      IF (.NOT. l_exist) THEN
         WRITE (message_text,'(a,a)') 'file not found:', TRIM(latbc_filename)
         CALL finish(TRIM(routine), message_text)
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

      ! initializing the displacement array for each compute processor
      ioff(:) = 0_MPI_ADDRESS_KIND

      !
      ! Perform CDI read operation
      !
      DO jm = 1, latbc_buffer%ngrp_vars
         ! Get pointer to appropriate reorder_info
         IF(latbc_buffer%nlev(jm) /= 1 ) THEN
            SELECT CASE (latbc_buffer%hgrid(jm))
            CASE(GRID_UNSTRUCTURED_CELL)
               ! Read 3d variables
               CALL prefetch_cdi_3d ( latbc_fileID, latbc_buffer%mapped_name(jm), patch_data, &
                    &                 latbc_buffer%nlev(jm), latbc_buffer%hgrid(jm), ioff )

            CASE(GRID_UNSTRUCTURED_EDGE)
               CALL prefetch_cdi_3d ( latbc_fileID, latbc_buffer%mapped_name(jm), patch_data, &
                    &                 latbc_buffer%nlev(jm), latbc_buffer%hgrid(jm), ioff )
            CASE default
               CALL finish(routine,'unknown grid type')
            END SELECT

         ELSE
            SELECT CASE (latbc_buffer%hgrid(jm))
            CASE(GRID_UNSTRUCTURED_CELL)
               ! Read 2d variables
               CALL prefetch_cdi_2d ( latbc_fileID, latbc_buffer%mapped_name(jm), patch_data, &
                    &                 latbc_buffer%hgrid(jm), ioff )

            CASE(GRID_UNSTRUCTURED_EDGE)
               CALL prefetch_cdi_2d ( latbc_fileID, latbc_buffer%mapped_name(jm), patch_data, &
                    &                 latbc_buffer%hgrid(jm), ioff )

            CASE default
               CALL finish(routine,'unknown grid type')
            END SELECT

         ENDIF
      ENDDO

      !
      ! close the open dataset file
      !
      CALL streamClose(latbc_fileID)
#endif
    END SUBROUTINE prefetch_latbc_intp_data

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
    SUBROUTINE compute_latbc_intp_data( p_patch, patch_data, p_nh_state, p_int)
      TYPE(t_patch),             TARGET   :: p_patch
      TYPE(t_patch_data), TARGET, INTENT(IN) :: patch_data
      TYPE(t_nh_state),       INTENT(IN)  :: p_nh_state  !< nonhydrostatic state on the global domain
      TYPE(t_int_state),      INTENT(IN)  :: p_int

#ifndef NOMPI
      ! local variables
      INTEGER(i8)                         :: eoff
      TYPE(t_reorder_data), POINTER       :: p_ri
      LOGICAL                             :: lconvert_omega2w
      INTEGER                             :: jc, jk, jb, jm, tlev, j, jl, jv
      INTEGER                             :: ioper_mode
      REAL(wp)                            :: log_exner, tempv
      CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = "mo_async_latbc_utils::compute_latbc_intp_data"

      nlev_in   = latbc_config%nlev_in
      tlev      = new_latbc_tlev

      ! Select operation mode determining the set of input fields to be allocated and read
      SELECT CASE (init_mode)
      CASE (MODE_COSMODE)
        ioper_mode = 2
      CASE (MODE_DWDANA, MODE_ICONVREMAP, MODE_IAU_OLD, MODE_IAU)
        ioper_mode = 3
        IF (.NOT. latbc_buffer%lthd_progvars) &
          CALL finish(TRIM(routine),'Limited area mode with ICON initial data requires prognostic&
            & thermodynamic variables in lateral boundary data')
      CASE DEFAULT
        ioper_mode = 1
      END SELECT

      ! Offset in memory window for async prefetching
      eoff = 0_i8

      DO jv = 1, latbc_buffer%ngrp_vars
         ! Receive 2d and 3d variables
         CALL compute_data_receive ( latbc_buffer%hgrid(jv), latbc_buffer%nlev(jv), &
                                     latbc_buffer%vars(jv)%buffer, eoff, patch_data)
      ENDDO

      IF (latbc_buffer%lthd_progvars .AND. ioper_mode /= 3) &
        CALL finish(TRIM(routine), 'Lateral boundary variables not consistent with init_mode')

      ! Reading the next time step
      IF((.NOT. my_process_is_io() .AND. &
           & .NOT. my_process_is_pref()) .AND. &
           & .NOT. my_process_is_mpi_test()) THEN
         CALL compute_start_async_pref()
      END IF

      ! Get patch ID
      p_ri => patch_data%cells

      ! copying the variable values from prefetch buffer to the respective allocated variable
!$OMP PARALLEL PRIVATE(jm,jv,jc)

      IF (latbc_buffer%lthd_progvars) THEN
        jm = get_field_index('theta_v')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
        DO jk=1, latbc_buffer%nlev(jm)
          DO j = 1, p_ri%n_own ! p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm_in%theta_v(jl,jk,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,jk,jb), wp)
          ENDDO
        ENDDO
!$OMP END DO
        jm = get_field_index('rho')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
        DO jk=1, latbc_buffer%nlev(jm)
          DO j = 1, p_ri%n_own ! p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm_in%rho(jl,jk,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,jk,jb), wp)
          ENDDO
        ENDDO
!$OMP END DO
      ELSE
        jm = get_field_index('temp')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
        DO jk=1, latbc_buffer%nlev(jm)
          DO j = 1, p_ri%n_own ! p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm_in%temp(jl,jk,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,jk,jb), wp)
          ENDDO
        ENDDO
!$OMP END DO
      ENDIF

      ! Read horizontal component of velocity (U and V)
      IF (latbc_buffer%lread_vn) THEN
         jm = get_field_index('vn')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc_buffer%nlev(jm)
            DO j = 1, patch_data%edges%n_own !p_patch%n_patch_cells
               jb = patch_data%edges%own_blk(j) ! Block index in distributed patch
               jl = patch_data%edges%own_idx(j) ! Line  index in distributed patch
               latbc_data(tlev)%atm_in%vn(jl,jk,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO
      ELSE
         jm = get_field_index('u')
         jv = get_field_index('v')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc_buffer%nlev(jm)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc_data(tlev)%atm_in%u(jl,jk,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,jk,jb), wp)
               latbc_data(tlev)%atm_in%v(jl,jk,jb) = REAL(latbc_buffer%vars(jv)%buffer(jl,jk,jb), wp)
               !          IF( p_pe_work == p_work_pe0 ) &
               !           &  WRITE(0,*) 'latbc_data(tlev)%atm_in%v value ', latbc_data(tlev)%atm_in%v(jl,jk,jb)
            ENDDO
         ENDDO
!$OMP END DO
      ENDIF

      ! Read vertical component of velocity (W)
      IF (ioper_mode == 1) THEN
         lconvert_omega2w = .TRUE.
         jv = get_field_index('w')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc_buffer%nlev(jv)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc_data(tlev)%atm_in%omega(jl,jk,jb) = REAL(latbc_buffer%vars(jv)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO
      ELSE
         lconvert_omega2w = .FALSE.
         jv = get_field_index('w')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc_buffer%nlev(jv)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc_data(tlev)%atm_in%w_ifc(jl,jk,jb) = REAL(latbc_buffer%vars(jv)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO
      ENDIF

      IF (ioper_mode >= 2) THEN
         ! Read parameter HHL
         jm = get_field_index('z_ifc')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc_buffer%nlev(jm)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc_data(tlev)%atm_in%z3d_ifc(jl,jk,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO

         ! Interpolate input 'z3d' and 'w' from the interface levels to the main levels
         !
!$OMP DO PRIVATE (jk,j,jb,jc) ICON_OMP_DEFAULT_SCHEDULE
         DO jk = 1, nlev_in    !!!!!!!!need to reset nlev_in to a new value from stored n_lev values
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jc = p_ri%own_idx(j) ! Line  index in distributed patch

               ! Note: In future, we want to z3d from boundary data.
               !
               latbc_data(tlev)%atm_in%z3d(jc,jk,jb) = (latbc_data(tlev)%atm_in%z3d_ifc(jc,jk,jb) + &
                    &   latbc_data(tlev)%atm_in%z3d_ifc(jc,jk+1,jb) ) * 0.5_wp
               latbc_data(tlev)%atm_in%w(jc,jk,jb) = (latbc_data(tlev)%atm_in%w_ifc(jc,jk,jb) + &
                    &   latbc_data(tlev)%atm_in%w_ifc(jc,jk+1,jb)) * 0.5_wp
            ENDDO
         ENDDO
!$OMP END DO
      ENDIF

      ! Read parameter QV, QC and QI
      jv = get_field_index('qv')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc_buffer%nlev(jv)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm_in%qv(jl,jk,jb) = REAL(latbc_buffer%vars(jv)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      jc = get_field_index('qc')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc_buffer%nlev(jc)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm_in%qc(jl,jk,jb) = REAL(latbc_buffer%vars(jc)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      jm = get_field_index('qi')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
      DO jk=1, latbc_buffer%nlev(jm)
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jl = p_ri%own_idx(j) ! Line  index in distributed patch
            latbc_data(tlev)%atm_in%qi(jl,jk,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,jk,jb), wp)
         ENDDO
      ENDDO
!$OMP END DO

      IF (latbc_buffer%lread_qr) THEN
         ! Read parameter QR
         jm = get_field_index('qr')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc_buffer%nlev(jm)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc_data(tlev)%atm_in%qr(jl,jk,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO
      ELSE
        CALL init(latbc_data(tlev)%atm_in%qr(:,:,:))
!$OMP BARRIER
      ENDIF

      IF (latbc_buffer%lread_qs) THEN
         ! Read parameter QS
         jv = get_field_index('qs')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc_buffer%nlev(jv)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc_data(tlev)%atm_in%qs(jl,jk,jb) = REAL(latbc_buffer%vars(jv)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO
      ELSE
        CALL init(latbc_data(tlev)%atm_in%qs(:,:,:))
!$OMP BARRIER
      ENDIF

      IF (ioper_mode == 2) THEN
         ! Read parameter Pressure
         jv = get_field_index('pres')
!$OMP DO PRIVATE (jk,j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO jk=1, latbc_buffer%nlev(jv)
            DO j = 1, p_ri%n_own !p_patch%n_patch_cells
               jb = p_ri%own_blk(j) ! Block index in distributed patch
               jl = p_ri%own_idx(j) ! Line  index in distributed patch
               latbc_data(tlev)%atm_in%pres(jl,jk,jb) = REAL(latbc_buffer%vars(jv)%buffer(jl,jk,jb), wp)
            ENDDO
         ENDDO
!$OMP END DO
      ENDIF

      IF (ioper_mode /= 3) THEN
        ! Read parameter surface pressure (LNPS)
        jm = get_field_index(latbc_buffer%psvar,.TRUE.)
!$OMP DO PRIVATE (j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
           jb = p_ri%own_blk(j) ! Block index in distributed patch
           jl = p_ri%own_idx(j) ! Line  index in distributed patch
           latbc_data(tlev)%atm_in%psfc(jl,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,1,jb), wp)
        ENDDO
!$OMP END DO
        ! Read parameter  surface Geopotential (GEOSP)
        jm = get_field_index(latbc_buffer%geop_ml_var,.TRUE.)
!$OMP DO PRIVATE (j,jb,jl) ICON_OMP_DEFAULT_SCHEDULE
         DO j = 1, p_ri%n_own !p_patch%n_patch_cells
           jb = p_ri%own_blk(j) ! Block index in distributed patch
           jl = p_ri%own_idx(j) ! Line  index in distributed patch
           latbc_data(tlev)%atm_in%phi_sfc(jl,jb) = REAL(latbc_buffer%vars(jm)%buffer(jl,1,jb), wp)
         ENDDO
!$OMP END DO
      ELSE
        ! Diagnose pres and temp from prognostic ICON variables
!$OMP DO PRIVATE (jk,j,jb,jc,log_exner,tempv)
        DO jk = 1, nlev_in
          DO j = 1, p_ri%n_own !p_patch%n_patch_cells
            jb = p_ri%own_blk(j) ! Block index in distributed patch
            jc = p_ri%own_idx(j) ! Line  index in distributed patch

            log_exner = (1._wp/cvd_o_rd)*LOG(latbc_data(tlev)%atm_in%rho(jc,jk,jb)* &
              latbc_data(tlev)%atm_in%theta_v(jc,jk,jb)*rd/p0ref)
            tempv = latbc_data(tlev)%atm_in%theta_v(jc,jk,jb)*EXP(log_exner)

            latbc_data(tlev)%atm_in%pres(jc,jk,jb) = p0ref*EXP(cpd/rd*log_exner)
            latbc_data(tlev)%atm_in%temp(jc,jk,jb) = tempv / (1._wp + vtmpc1*latbc_data(tlev)%atm_in%qv(jc,jk,jb) - &
              (latbc_data(tlev)%atm_in%qc(jc,jk,jb) + latbc_data(tlev)%atm_in%qi(jc,jk,jb) +                        &
               latbc_data(tlev)%atm_in%qr(jc,jk,jb) + latbc_data(tlev)%atm_in%qs(jc,jk,jb)) )

          ENDDO
        ENDDO
!$OMP END DO
      ENDIF

!$OMP END PARALLEL

      ! boundary exchange for a 2-D and 3-D array, needed because the
      ! vertical interpolation includes the halo region (otherwise, the 
      ! syncs would have to be called after vert_interp)
      CALL sync_patch_array_mult(SYNC_C,p_patch,4,latbc_data(tlev)%atm_in%temp,latbc_data(tlev)%atm_in%z3d, &
         latbc_data(tlev)%atm_in%w,latbc_data(tlev)%atm_in%pres)
      CALL sync_patch_array_mult(SYNC_C,p_patch,5,latbc_data(tlev)%atm_in%qv,latbc_data(tlev)%atm_in%qc, &
         latbc_data(tlev)%atm_in%qi,latbc_data(tlev)%atm_in%qr,latbc_data(tlev)%atm_in%qs)
      IF (latbc_buffer%lread_vn) THEN
         CALL sync_patch_array(SYNC_E,p_patch,latbc_data(tlev)%atm_in%vn)
      ELSE
         CALL sync_patch_array_mult(SYNC_C,p_patch,2,latbc_data(tlev)%atm_in%u,latbc_data(tlev)%atm_in%v)
      ENDIF
      IF (lconvert_omega2w) CALL sync_patch_array(SYNC_C,p_patch,latbc_data(tlev)%atm_in%omega)
      IF (.NOT. latbc_buffer%lthd_progvars) THEN
        CALL sync_patch_array(SYNC_C,p_patch,latbc_data(tlev)%atm_in%phi_sfc)
        CALL sync_patch_array(SYNC_C,p_patch,latbc_data(tlev)%atm_in%psfc)
      ENDIF

      ! perform vertical interpolation of horizonally interpolated analysis data
      !
      CALL vert_interp(p_patch, p_int, p_nh_state%metrics, nlev_in, latbc_data(tlev),              &
           &    opt_convert_omega2w=lconvert_omega2w, opt_use_vn=latbc_buffer%lread_vn)
#endif
    END SUBROUTINE compute_latbc_intp_data

    !-------------------------------------------------------------------------
    !>
    !! @par Revision History
    !! Initial version by S. Brdar, DWD (2013-06-13)
    !! Modified version by M. Pondkule, DWD (2013-04-17)
    !!
    SUBROUTINE deallocate_pref_latbc_data()

#ifndef NOMPI
      ! local variables
      INTEGER             :: tlev
      CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = &
           "mo_async_latbc_utils::deallocate_latbc_data"

      WRITE(message_text,'(a,a)') 'deallocating latbc data'
      CALL message(TRIM(routine), message_text)

      IF(my_process_is_work()) THEN
         !
         ! deallocate boundary data memory
         !
         DO tlev = 1, 2
            DEALLOCATE(latbc_data(tlev)%atm_in%psfc, &
                 latbc_data(tlev)%atm_in%phi_sfc, &
                 latbc_data(tlev)%atm_in%pres, &
                 latbc_data(tlev)%atm_in%temp, &
                 latbc_data(tlev)%atm_in%z3d, &
                 latbc_data(tlev)%atm_in%u, &
                 latbc_data(tlev)%atm_in%v, &
                 latbc_data(tlev)%atm_in%w, &
                 latbc_data(tlev)%atm_in%omega, &
                 latbc_data(tlev)%atm_in%qv, &
                 latbc_data(tlev)%atm_in%qc, &
                 latbc_data(tlev)%atm_in%qi, &
                 latbc_data(tlev)%atm_in%qr, &
                 latbc_data(tlev)%atm_in%qs )

            IF (init_mode == MODE_COSMODE) THEN
               DEALLOCATE(latbc_data(tlev)%atm_in%w_ifc)
               DEALLOCATE(latbc_data(tlev)%atm_in%z3d_ifc)
            ENDIF

            ! Allocate atmospheric output data
            DEALLOCATE(latbc_data(tlev)%atm%vn, &
                 latbc_data(tlev)%atm%u, &
                 latbc_data(tlev)%atm%v, &
                 latbc_data(tlev)%atm_in%vn, &
                 latbc_data(tlev)%atm%w, &
                 latbc_data(tlev)%atm%temp, &
                 latbc_data(tlev)%atm%exner, &
                 latbc_data(tlev)%atm%pres, &
                 latbc_data(tlev)%atm%rho, &
                 latbc_data(tlev)%atm%theta_v, &
                 latbc_data(tlev)%atm%qv, &
                 latbc_data(tlev)%atm%qc, &
                 latbc_data(tlev)%atm%qi, &
                 latbc_data(tlev)%atm%qr, &
                 latbc_data(tlev)%atm%qs )

         END DO
      ENDIF

      ! deallocating prefetch input event
      CALL deallocateEvent(prefetchEvent)
      ! deallocating Date and time
      CALL deallocateDatetime(mtime_date)
      CALL deallocateDatetime(mtime_read)
      CALL deallocateDatetime(mtime_end)
      CALL deallocateTimedelta(my_duration_slack)
      CALL deallocateTimedelta(delta_dtime)
#endif
    END SUBROUTINE deallocate_pref_latbc_data

    !-------------------------------------------------------------------------
    !>
    !! @par Revision History
    !! Initial version by G. Zaengl, DWD (2013-10-22)
    !!
    SUBROUTINE compute_boundary_tendencies ( p_patch, p_nh )
      TYPE(t_patch),    INTENT(IN)    :: p_patch
      TYPE(t_nh_state), INTENT(INOUT) :: p_nh

#ifndef NOMPI
      ! Local variables
      INTEGER                         :: i_startblk, i_endblk, &
           &                                i_startidx, i_endidx
      INTEGER                         :: jc, jk, jb, je
      INTEGER                         :: nlev, nlevp1
      REAL(wp)                        :: rdt

      nlev = p_patch%nlev
      nlevp1 = p_patch%nlevp1

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
                    &   latbc_data(new_latbc_tlev)%atm%vn(je,jk,jb) &
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
                    &   latbc_data(new_latbc_tlev)%atm%rho(jc,jk,jb) &
                    & - latbc_data(prev_latbc_tlev)%atm%rho(jc,jk,jb) )

               p_nh%diag%grf_tend_thv(jc,jk,jb) = rdt * (                &
                    &   latbc_data(new_latbc_tlev)%atm%theta_v(jc,jk,jb) &
                    & - latbc_data(prev_latbc_tlev)%atm%theta_v(jc,jk,jb) )

               p_nh%diag%grf_tend_w(jc,jk,jb) = rdt * (            &
                    &   latbc_data(new_latbc_tlev)%atm%w(jc,jk,jb) &
                    & - latbc_data(prev_latbc_tlev)%atm%w(jc,jk,jb) )

            ENDDO
         ENDDO

         DO jc = i_startidx, i_endidx
            p_nh%diag%grf_tend_w(jc,nlevp1,jb) = rdt * (            &
                 &   latbc_data(new_latbc_tlev)%atm%w(jc,nlevp1,jb) &
                 & - latbc_data(prev_latbc_tlev)%atm%w(jc,nlevp1,jb) )
         ENDDO

         IF (ltransport) THEN
            DO jk = 1, nlev
               DO jc = i_startidx, i_endidx
                  p_nh%diag%grf_tend_tracer(jc,jk,jb,iqv) =  rdt * (   &
                       &   latbc_data(new_latbc_tlev)%atm%qv(jc,jk,jb) &
                       & - latbc_data(prev_latbc_tlev)%atm%qv(jc,jk,jb) )

                  p_nh%diag%grf_tend_tracer(jc,jk,jb,iqc) =  rdt * (   &
                       &   latbc_data(new_latbc_tlev)%atm%qc(jc,jk,jb) &
                       & - latbc_data(prev_latbc_tlev)%atm%qc(jc,jk,jb) )

                  p_nh%diag%grf_tend_tracer(jc,jk,jb,iqi) =  rdt * (   &
                       &   latbc_data(new_latbc_tlev)%atm%qi(jc,jk,jb) &
                       & - latbc_data(prev_latbc_tlev)%atm%qi(jc,jk,jb) )

                  p_nh%diag%grf_tend_tracer(jc,jk,jb,iqr) =  rdt * (   &
                       &   latbc_data(new_latbc_tlev)%atm%qr(jc,jk,jb) &
                       & - latbc_data(prev_latbc_tlev)%atm%qr(jc,jk,jb) )

                  p_nh%diag%grf_tend_tracer(jc,jk,jb,iqs) =  rdt * (   &
                       &   latbc_data(new_latbc_tlev)%atm%qs(jc,jk,jb) &
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
      REAL(wp) :: msg
      ! Simply send a message from Input prefetching PE 0 to work PE 0
      ! p_pe_work == 0 signifies processor 0 in Input prefetching PEs
      ! Note: We have to do this in a non-blocking fashion in order to
      !       receive "ready file" messages
      IF(p_pe_work == 0) THEN
         !     CALL p_wait()
         msg = REAL(msg_pref_done, wp)
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
      REAL(wp) :: msg
      ! First compute PE receives message from input prefetching leader
      IF(p_pe_work==0) THEN
         CALL p_recv(msg, p_pref_pe0, 0)
         ! Just for safety: Check if we got the correct tag
         IF(INT(msg) /= msg_pref_done) CALL finish(modname, 'Compute PE: Got illegal prefetching tag')
      ENDIF
      ! Wait in barrier until message is here
      CALL p_barrier(comm=p_comm_work)
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
      LOGICAL, INTENT(INOUT)          :: done ! flag if we should shut down

#ifndef NOMPI
      REAL(wp) :: msg(2)
      CHARACTER(*), PARAMETER :: method_name = "async_pref_wait_for_start"

      ! Set output parameters to default values
      done  = .FALSE.

      ! Receive message that we may start transferring the prefetching data (or should finish)
      ! Note: We have to do this in a non-blocking fashion in order to
      !       receive "ready file" messages.
      IF(p_pe_work == 0) THEN
         ! launch non-blocking receive request
         CALL p_wait()
         CALL p_irecv(msg, p_work_pe0, 0)
         CALL p_wait()
      END IF

      SELECT CASE(INT(msg(1)))
      CASE(msg_pref_start)

      CASE(msg_pref_shutdown)
         done = .TRUE.

      CASE DEFAULT
         ! Anything else is an error
         CALL finish(modname, 'Prefetching PE: Got illegal prefetching tag')
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
      REAL(wp) :: msg(2)
      CHARACTER(*), PARAMETER :: method_name = "compute_start_async_pref"

      !   CALL p_barrier(comm=p_comm_work) ! make sure all are here

      msg(1) = REAL(msg_pref_start,  wp)

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
      REAL(wp) :: msg(2)

      !  CALL p_barrier(comm=p_comm_work) ! make sure all are here

      msg(1) = REAL(msg_pref_shutdown, wp)
      msg(2) = 0._wp

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
    FUNCTION get_field_index(name,opt_lmap) RESULT(result_varID)
      CHARACTER (LEN=*), INTENT(IN) :: name !< variable name
      LOGICAL, OPTIONAL, INTENT(IN) :: opt_lmap
      INTEGER :: result_varID, varID
      LOGICAL :: lmap

      lmap = .FALSE.
      IF (PRESENT(opt_lmap)) THEN
        IF (opt_lmap) lmap = .TRUE.
      ENDIF

      result_varID = -1
      ! looping over variable list in internal or mapped name
      IF (lmap) THEN
        DO varID=1, latbc_buffer%ngrp_vars
          IF (TRIM(name) == TRIM(latbc_buffer%mapped_name(varID))) THEN
            result_varID = varID
          END IF
        END DO
      ELSE
        DO varID=1, latbc_buffer%ngrp_vars
          IF (TRIM(name) == TRIM(latbc_buffer%internal_name(varID))) THEN
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
    SUBROUTINE update_lin_interpolation( step_datetime )
      TYPE(t_datetime), INTENT(INOUT) :: step_datetime
      TYPE(datetime), pointer :: mtime_step
      TYPE(timedelta), pointer :: delta_tstep
      CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: sim_step
      CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: routine = &
           "mo_async_latbc_utils::update_lin_interpolation"

      ! compute current datetime in a format appropriate for mtime
      CALL get_datetime_string(sim_step, step_datetime)
      mtime_step  => newDatetime(TRIM(sim_step))

      delta_tstep => newTimedelta("PT0S")
      delta_tstep = mtime_read - mtime_step

      IF(delta_tstep%month /= 0) &
           CALL finish(TRIM(routine), "time difference for reading boundary data cannot be more than a month.")

      ! compute the number of "dtime_latbc" intervals fitting into the time difference "delta_tstep":

      latbc_config%lc1 = (delta_tstep%day * 86400._wp + delta_tstep%hour * 3600._wp +  &
           delta_tstep%minute * 60._wp + delta_tstep%second) / latbc_config%dtime_latbc
      latbc_config%lc2 = 1._wp - latbc_config%lc1

      ! deallocating mtime and deltatime
      CALL deallocateTimedelta(delta_tstep)
      CALL deallocateDatetime(mtime_step)

    END SUBROUTINE update_lin_interpolation

    !-------------------------------------------------------------------------

  END MODULE mo_async_latbc_utils
