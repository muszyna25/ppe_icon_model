!> Routines that describe when regular output steps and ready file
!> events shall be triggered. The definitions in this module control
!> the behaviour of "mo_output_event_handler".
!!
!! See "mo_output_event_handler" for a detailed description.
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
MODULE mo_output_event_control

  USE mo_mpi,                ONLY: my_process_is_mpi_test, my_process_is_mpi_workroot
  USE mo_impl_constants,     ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_exception,          ONLY: finish
  USE mo_kind,               ONLY: wp, i4, i8
  USE mo_master_config,      ONLY: getModelBaseDir
  USE mtime,                 ONLY: MAX_DATETIME_STR_LEN, MAX_DATETIME_STR_LEN,          &
    &                              MAX_TIMEDELTA_STR_LEN, datetime,                     &
    &                              deallocateDatetime, datetimeToString,                &
    &                              newDatetime, OPERATOR(>=), OPERATOR(*),              &
    &                              OPERATOR(+), OPERATOR(-), timedelta, newTimedelta,   &
    &                              deallocateTimedelta, OPERATOR(<=), OPERATOR(>),      &
    &                              OPERATOR(<), OPERATOR(==),                           &
    &                              divisionquotienttimespan,                            &
    &                              dividedatetimedifferenceinseconds,                   &
    &                              getPTStringFromMS, getPTStringFromSeconds,           &
    &                              timedeltaToString
  USE mo_var_list_element,   ONLY: lev_type_str
  USE mo_output_event_types, ONLY: t_sim_step_info, t_event_step_data
  USE mo_util_string,        ONLY: t_keyword_list, associate_keyword, with_keywords,    &
    &                              int2string, tolower
  USE mo_name_list_output_types, ONLY: t_fname_metadata


  IMPLICIT NONE

  ! public subroutines + functions:
  PUBLIC :: compute_matching_sim_steps
  PUBLIC :: generate_output_filenames

  !---------------------------------------------------------------
  ! constants
  
  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_output_event_control'

  !> Internal switch for debugging output
  LOGICAL,          PARAMETER :: ldebug  = .FALSE.

CONTAINS

  ! --------------------------------------------------------------------------------------------------
  !> Compute the matching simulation step for each date-time stamp.
  !
  !  The list of date-time stamp string is given in @p date_string. We
  !  assume that the date_string list is ordered in time. We choose
  !  the advection time step matching or following the time stamp
  !  string. The "true" date string is also returned.
  !
  !  @note The mtime library offers no arithmetic subroutines to
  !        compute dates. Therefore we must actually step through all
  !        dynamic timesteps.
  !
  !  @author F. Prill, DWD
  !
  ! --------------------------------------------------------------------------------------------------
  SUBROUTINE compute_matching_sim_steps(nstrings, date_string, sim_step_info, &
    &                                   result_steps, result_exactdate)

    INTEGER,              INTENT(IN)    :: nstrings                             !< no. of string to convert
    CHARACTER(len=*),     INTENT(IN)    :: date_string(:)                       !< array of ISO 8601 time stamp strings
    TYPE(t_sim_step_info),INTENT(IN)    :: sim_step_info                        !< definitions: time step size, etc.
    INTEGER,              INTENT(INOUT) :: result_steps(:)                      !< resulting step indices (+last sim step)
    CHARACTER(LEN=*),     INTENT(INOUT) :: result_exactdate(:)                  !< resulting (exact) time step strings  (+last sim step)

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::compute_matching_sim_steps"
    INTEGER                  :: idtime_ms, ilist
    TYPE(datetime),  POINTER :: mtime_begin, mtime_end, mtime_date1, &
         &                      mtime_dom_start, mtime_dom_end
    TYPE(timedelta), POINTER :: delta
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: dtime_string

    ! debugging output
    IF (ldebug)  WRITE (0,*) "date_string: ", date_string(1:nstrings)

    ! build an ISO 8601 duration string from the given "dtime" value:
    idtime_ms = NINT(sim_step_info%dtime*1000._wp)
    CALL getPTStringFromMS(INT(idtime_ms,i8), dtime_string)
    ! create a time delta of "dtime" seconds length
    delta => newTimedelta(TRIM(dtime_string))

    ! Domains (and their output) can be activated and deactivated
    ! during the simulation. This is determined by the parameters
    ! "dom_start_time" and "dom_end_time". Therefore, we must create
    ! a corresponding event.
    mtime_dom_start => newDatetime(TRIM(sim_step_info%dom_start_time))
    mtime_dom_end   => newDatetime(TRIM(sim_step_info%dom_end_time))
    mtime_begin     => newDatetime(TRIM(sim_step_info%run_start))
    mtime_end       => newDatetime(TRIM(sim_step_info%sim_end  ))

    result_steps(:)     = -1
    result_exactdate(:) = ""

    DO ilist = 1,nstrings
      mtime_date1 => newDatetime(TRIM(date_string(ilist)))

      ! check if domain is inactive:
      IF (((mtime_date1 >= mtime_dom_start) .AND. (mtime_date1 >= mtime_begin)) .OR.  &
        & ((mtime_date1 <= mtime_end) .AND. (mtime_date1 <= mtime_dom_end))) THEN 
        CALL compute_step(mtime_date1, mtime_begin, mtime_end,                &
          &               sim_step_info%dtime, delta,                         &
          &               sim_step_info%jstep0,                               &
          &               result_steps(ilist), result_exactdate(ilist))
        IF (ldebug) THEN
          WRITE (0,*) ilist, ": ", result_steps(ilist), " -> ", TRIM(result_exactdate(ilist))
        END IF
      END IF
      CALL deallocateDatetime(mtime_date1)
    END DO
    ! clean up
    CALL deallocateDatetime(mtime_dom_start)
    CALL deallocateDatetime(mtime_dom_end)
    CALL deallocateDatetime(mtime_begin)
    CALL deallocateDatetime(mtime_end)
    CALL deallocateTimedelta(delta)

  END SUBROUTINE compute_matching_sim_steps


  ! --------------------------------------------------------------------------------------------------
  !> Utility function: Compute "(date1-sim_start)/(dtime)"
  !
  !  @author F. Prill, DWD
  ! --------------------------------------------------------------------------------------------------
  SUBROUTINE compute_step(mtime_current, mtime_begin, mtime_end, dtime,  &
    &                     delta, step_offset, step, exact_date)
    TYPE(datetime),  POINTER                         :: mtime_current       !< input date to translated into step
    TYPE(datetime),  POINTER                         :: mtime_begin         !< begin of run (note: restart cases!)
    TYPE(datetime),  POINTER                         :: mtime_end           !< end of run
    REAL(wp),                            INTENT(IN)  :: dtime               !< [s] length of a time step
    TYPE(timedelta), POINTER                         :: delta
    INTEGER,                             INTENT(IN)  :: step_offset
    INTEGER,                             INTENT(OUT) :: step                !< result: corresponding simulations step
    CHARACTER(len=MAX_DATETIME_STR_LEN), INTENT(OUT) :: exact_date          !< result: corresponding simulation date
    ! local variables
    REAL                                 :: intvlmillisec
    TYPE(datetime),  POINTER             :: mtime_step
    CHARACTER(len=max_timedelta_str_len) :: td_string
    TYPE(divisionquotienttimespan)       :: tq     
    TYPE(timedelta), POINTER             :: vlsec => NULL()

    ! first, we compute the dynamic time step which is equal or larger than
    ! the desired date "mtime_current"
    ! intvlsec    = REAL(dtime)
    ! step        = CEILING(datetimedividebyseconds(mtime_begin, mtime_date1, intvlsec))

    intvlmillisec = NINT(dtime*1000._wp)
    CALL getPTStringFromMS(INT(intvlmillisec,i8), td_string)
    !CALL getptstringfromseconds(INT(intvlsec,i8), td_string)
    vlsec => newtimedelta(td_string)

    CALL divideDatetimeDifferenceInSeconds(mtime_current, mtime_begin, vlsec, tq)

    step = INT(tq%quotient,i4)
    
    mtime_step => newDatetime('0001-01-01T00:00:00')
    IF (step >= 0) THEN
      mtime_step = mtime_begin + step * vlsec
      CALL datetimeToString(mtime_step, exact_date)
    END IF
    CALL deallocateDatetime(mtime_step)

    ! then we add the offset "jstep0" (nonzero for restart cases):
    step        = step + step_offset

  END SUBROUTINE compute_step


  ! --------------------------------------------------------------------------------------------------
  !> Function for generating output file names.
  !
  !  This routine also prescribes, in which steps the output file is
  !  opened or closed.
  !
  !  @author F. Prill, DWD
  !
  ! --------------------------------------------------------------------------------------------------
  FUNCTION generate_output_filenames(nstrings, date_string, sim_steps, &
    &                                sim_step_info, fname_metadata, skipped_dates)  RESULT(result_fnames)
    INTEGER,                INTENT(IN)    :: nstrings           !< no. of string to convert
    CHARACTER(len=*),       INTENT(IN)    :: date_string(:)     !< array of ISO 8601 time stamp strings
    INTEGER,                INTENT(IN)    :: sim_steps(:)       !< array of corresponding simulation steps
    TYPE(t_sim_step_info),  INTENT(IN)    :: sim_step_info      !< definitions: time step size, etc.
    TYPE(t_fname_metadata), INTENT(IN)    :: fname_metadata     !< additional meta-data for generating output filename
    INTEGER,                INTENT(IN)    :: skipped_dates
    TYPE(t_event_step_data) :: result_fnames(SIZE(date_string))
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::generate_output_filenames"
    INTEGER                             :: i, j, ifile, ipart, total_index, this_jfile
    CHARACTER(len=MAX_CHAR_LENGTH)      :: cfilename 
    TYPE (t_keyword_list), POINTER      :: keywords     => NULL()
    CHARACTER(len=MAX_CHAR_LENGTH)      :: fname(nstrings)        ! list for duplicate check
    INTEGER                             :: ifname                 ! current length of "ifname"
    TYPE(datetime),  POINTER            :: file_end, step_date, run_start, sim_start, &
      &                                    mtime_first, mtime_date
    TYPE(timedelta), POINTER            :: delta, forecast_delta
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: dtime_string, forecast_delta_str

    run_start  => newDatetime(TRIM(sim_step_info%run_start))
    sim_start  => newDatetime(TRIM(sim_step_info%sim_start))

    ! ---------------------------------------------------
    ! prescribe, which file is used in which output step.
    ! ---------------------------------------------------
    ! a) user has specified "steps_per_file"
    ! b) user has specified "file_interval"
    IF (TRIM(fname_metadata%file_interval) == "") THEN
      ! case a): steps_per_file
      mtime_first  => newDatetime(TRIM(date_string(1)))
      ! special treatment for the initial time step written at the
      ! begin of the simulation: The first file gets one extra entry
      IF ( (run_start == mtime_first)              .AND.  &
        &  fname_metadata%steps_per_file_inclfirst   .AND.  &
        &  (fname_metadata%jfile_offset == 0) ) THEN
        result_fnames(1)%jfile = 1
        result_fnames(1)%jpart = 1
        DO i=2,nstrings
          result_fnames(i)%jfile = (i-2)/fname_metadata%steps_per_file + 1 ! INTEGER division!
          result_fnames(i)%jpart = i - 1 - (result_fnames(i)%jfile-1)*fname_metadata%steps_per_file
          IF (result_fnames(i)%jfile == 1)  result_fnames(i)%jpart = result_fnames(i)%jpart + 1
        END DO
      ELSE
        DO i=1,nstrings
          result_fnames(i)%jfile = (i-1)/fname_metadata%steps_per_file + 1 ! INTEGER division!
          result_fnames(i)%jpart = i - (result_fnames(i)%jfile-1)*fname_metadata%steps_per_file
        END DO
      END IF
      CALL deallocateDatetime(mtime_first)
    ELSE
      ! case b): file_interval
      ifile      =  1
      ipart      =  0
      delta     =>  newTimedelta(TRIM(fname_metadata%file_interval))
      file_end  =>  newDatetime(date_string(1))

      file_end   =  file_end + delta
      OUTSTEP_LOOP : DO i=1,nstrings
        step_date => newDatetime(date_string(i))
        IF (ldebug) THEN
          CALL datetimeToString(file_end, dtime_string)
          WRITE (0,*) "step_date = ", date_string(i), "; file_end = ", dtime_string
        END IF
        IF (step_date >= file_end) THEN
          ifile = ifile + 1
          ipart = 0
          DO
            file_end =  file_end + delta
            IF (file_end > step_date) EXIT
          END DO
        END IF
        ipart = ipart + 1
        result_fnames(i)%jfile = ifile
        result_fnames(i)%jpart = ipart
        CALL deallocateDatetime(step_date)
      END DO OUTSTEP_LOOP
      CALL deallocateDatetime(file_end)
      CALL deallocateTimedelta(delta)
    END IF
    ! add offset to file numbers:
    IF (fname_metadata%jfile_offset /= 0) THEN
      result_fnames(1:nstrings)%jfile = result_fnames(1:nstrings)%jfile + fname_metadata%jfile_offset - 1
    END IF

    ! --------------------------------------------------------------
    ! prescribe, in which steps the output file is opened or closed:
    ! --------------------------------------------------------------

    DO i=1,nstrings
      IF (i == 1) THEN
        result_fnames(i)%l_open_file = .TRUE.
      ELSE
        result_fnames(i)%l_open_file = (result_fnames(i-1)%jfile /= result_fnames(i)%jfile)
      END IF
      IF (i==nstrings) THEN
        result_fnames(i)%l_close_file = .TRUE.
      ELSE
        result_fnames(i)%l_close_file = (result_fnames(i+1)%jfile /= result_fnames(i)%jfile)
      END IF
    END DO

    ! ----------------------------------------------
    ! Set actual output file name (insert keywords):
    ! ----------------------------------------------
    forecast_delta => newTimedelta("P01D")
    ifname = 0
    DO i=1,nstrings
      IF (.NOT. result_fnames(i)%l_open_file) THEN
        ! if no file is opened: filename is identical to last step
        result_fnames(i)%filename_string = result_fnames(i-1)%filename_string
        CYCLE
      END IF
      ! otherwise, generate filename
      !
      ! define keywords:
      CALL associate_keyword("<path>",            TRIM(getModelBaseDir()),                                  keywords)
      CALL associate_keyword("<output_filename>", TRIM(fname_metadata%filename_pref),                       keywords)
      CALL associate_keyword("<physdom>",         TRIM(int2string(fname_metadata%phys_patch_id, "(i2.2)")), keywords)
      CALL associate_keyword("<levtype>",         TRIM(lev_type_str(fname_metadata%ilev_type)),             keywords)
      CALL associate_keyword("<levtype_l>",       TRIM(tolower(lev_type_str(fname_metadata%ilev_type))),    keywords)
      CALL associate_keyword("<npartitions>",     TRIM(int2string(fname_metadata%npartitions)),             keywords)
      CALL associate_keyword("<ifile_partition>", TRIM(int2string(fname_metadata%ifile_partition)),         keywords)
      CALL associate_keyword("<datetime>",        TRIM(date_string(i)),                                     keywords)
      ! keywords: compute current forecast time (delta):
      mtime_date => newDatetime(TRIM(date_string(i)))
      forecast_delta = mtime_date - sim_start
      WRITE (forecast_delta_str,'(4(i2.2))') forecast_delta%day, forecast_delta%hour, &
        &                                    forecast_delta%minute, forecast_delta%second 
      CALL associate_keyword("<ddhhmmss>",        TRIM(forecast_delta_str),                                 keywords)
      forecast_delta_str = ""
      WRITE (forecast_delta_str,'(i3.3,2(i2.2))') forecast_delta%day*24 + forecast_delta%hour, &
        &                                         forecast_delta%minute, forecast_delta%second 
      CALL associate_keyword("<hhhmmss>",        TRIM(forecast_delta_str),                                 keywords)

      ! keywords: compose other variants of the absolute date-time
      !
      ! "YYYYMMDDThhmmssZ"     for the basic format of ISO8601 without the
      !                        separators "-" and ":" of the extended date-time format
      WRITE (dtime_string,'(i4.4,2(i2.2),a,3(i2.2),a)')                                                 &
        &                      mtime_date%date%year, mtime_date%date%month, mtime_date%date%day, 'T',   &
        &                      mtime_date%time%hour, mtime_date%time%minute, mtime_date%time%second, 'Z'
      CALL associate_keyword("<datetime2>",       TRIM(dtime_string),                                       keywords)

      ! "YYYYMMDDThhmmss.sssZ" for the basic format of ISO8601 with 3-digit 
      !                        fractions of seconds
      WRITE (dtime_string,'(i4.4,2(i2.2),a,3(i2.2),a,i3.3,a)')                                                 &
        &                      mtime_date%date%year, mtime_date%date%month, mtime_date%date%day, 'T',          &
        &                      mtime_date%time%hour, mtime_date%time%minute, mtime_date%time%second, '.',      &
        &                      mtime_date%time%ms, 'Z'
      CALL associate_keyword("<datetime3>",       TRIM(dtime_string),                                       keywords)
      CALL deallocateDatetime(mtime_date)

      ! "total_index": If we have split one namelist into concurrent,
      !                alternating files, compute the file index,
      !                which an "unsplit" namelist would have
      !                produced:
      IF (fname_metadata%npartitions > 1) THEN
        total_index = fname_metadata%npartitions*(result_fnames(i)%jfile+skipped_dates-1) + &
          &           fname_metadata%ifile_partition 
        this_jfile  = total_index
      ELSE
        total_index = result_fnames(i)%jfile
        this_jfile  = result_fnames(i)%jfile
      END IF
      CALL associate_keyword("<total_index>", TRIM(int2string(total_index,"(i4.4)")),  keywords)
      CALL associate_keyword("<jfile>",       TRIM(int2string(this_jfile, "(i4.4)")), keywords)

      cfilename = TRIM(with_keywords(keywords, fname_metadata%filename_format))
      IF(my_process_is_mpi_test()) THEN
         ! (a) use filename with extension "_TEST" (then any post-processing needs to be adapted)
         WRITE(result_fnames(i)%filename_string,'(a,"_TEST",a)') TRIM(cfilename),TRIM(fname_metadata%extn)
!        ! (b) use standard filename
!        WRITE(result_fnames(i)%filename_string,'(a,a)')         TRIM(cfilename),TRIM(fname_metadata%extn)
      ELSE
        WRITE(result_fnames(i)%filename_string,'(a,a)')         TRIM(cfilename),TRIM(fname_metadata%extn)
      ENDIF

      ! consistency check: test if the user has accidentally chosen a
      ! file name syntax which yields duplicate file names:
      DO j=1,ifname
        IF (result_fnames(i)%filename_string == fname(j)) THEN
          ! found a duplicate filename:
          CALL finish(routine, "Ambiguous output file name: '"//TRIM(fname(j))//"'")
        END IF
      END DO
      ifname = ifname + 1
      fname(ifname) = result_fnames(i)%filename_string
    END DO
    CALL deallocateTimedelta(forecast_delta)

    ! clean up
    CALL deallocateDatetime(run_start)
    CALL deallocateDatetime(sim_start)

  END FUNCTION generate_output_filenames

END MODULE mo_output_event_control
