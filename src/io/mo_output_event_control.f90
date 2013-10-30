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
MODULE mo_output_event_control

  USE mo_mpi,                ONLY: my_process_is_mpi_test
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_exception,          ONLY: finish
  USE mo_kind,               ONLY: wp
  USE mo_master_nml,         ONLY: model_base_dir
  USE mtime,                 ONLY: MAX_DATETIME_STR_LEN, MAX_DATETIME_STR_LEN,          &
    &                              MAX_TIMEDELTA_STR_LEN, PROLEPTIC_GREGORIAN,          &
    &                              event, datetime, newEvent,                           &
    &                              setCalendar, resetCalendar,                          &
    &                              deallocateDatetime, datetimeToString,                &
    &                              deallocateEvent, newDatetime, OPERATOR(>=),          &
    &                              OPERATOR(+), timedelta, newTimedelta,                &
    &                              deallocateTimedelta, OPERATOR(<=), OPERATOR(>),      &
    &                              OPERATOR(<), OPERATOR(==)
  USE mo_mtime_extensions,   ONLY: isCurrentEventActive, getPTStringFromMS
  USE mo_var_list_element,   ONLY: lev_type_str
  USE mo_output_event_types, ONLY: t_sim_step_info, t_event_step_data
  USE mo_util_string,        ONLY: t_keyword_list, associate_keyword, with_keywords,    &
    &                              int2string, tolower, MAX_STRING_LEN
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
  !  the dynamic time step matching or following the time stamp
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
    INTEGER                  :: idtime_ms, sim_step, ilist, i
    LOGICAL                  :: l_dyn_active, l_adv_active
    TYPE(event),     POINTER :: mtime_dyn, mtime_adv
    TYPE(datetime),  POINTER :: mtime_dyn_date, mtime_cur_date, mtime_begin, mtime_end, &
      &                         mtime_dom_start, mtime_dom_end
    TYPE(timedelta), POINTER :: delta
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: dtime_string

    CALL setCalendar(PROLEPTIC_GREGORIAN)

    ! debugging output
    IF (ldebug)  WRITE (0,*) "date_string: ", date_string(1:nstrings)

    ! build an ISO 8601 duration string from the given "dtime" value:
    idtime_ms = NINT(sim_step_info%dtime*1000._wp)
    CALL getPTStringFromMS(idtime_ms, dtime_string)
    ! create a time delta of "dtime" seconds length
    delta => newTimedelta(TRIM(dtime_string))

    ! Domains (and their output) can be activated and deactivated
    ! during the simulation. This is determined by the parameters
    ! "dom_start_time" and "dom_end_time". Therefore, we must create
    ! a corresponding event.
    mtime_dom_start => newDatetime(TRIM(sim_step_info%dom_start_time))
    mtime_dom_end   => newDatetime(TRIM(sim_step_info%dom_end_time))

    mtime_begin     => newDatetime(TRIM(sim_step_info%sim_start))
    mtime_end       => newDatetime(TRIM(sim_step_info%sim_end  ))

    ! create an event for the dynamic time stepping with the "mtime" library:
    mtime_dyn => newEvent("dynamic",                                 &  ! name
      &                   TRIM(sim_step_info%sim_start),             &  ! ref.date
      &                   TRIM(sim_step_info%sim_start),             &  ! first
      &                   TRIM(sim_step_info%sim_end),               &  ! last
      &                   TRIM(dtime_string))                           ! interval

    ! create an event for the advection time stepping with the "mtime" library:
    idtime_ms = sim_step_info%iadv_rcf * idtime_ms
    CALL getPTStringFromMS(idtime_ms, dtime_string)
    mtime_adv => newEvent("advection",                               &  ! name
      &                   TRIM(sim_step_info%sim_start),             &  ! ref.date
      &                   TRIM(sim_step_info%sim_start),             &  ! first
      &                   TRIM(sim_step_info%sim_end),               &  ! last
      &                   TRIM(dtime_string))                           ! interval

    sim_step            =  0
    ilist               =  0
    result_steps(:)     = -1
    result_exactdate(:) = ""

    mtime_dyn_date => newDatetime(TRIM(sim_step_info%sim_start))
    mtime_cur_date => newDatetime(date_string(1))
    EVENT_LOOP: DO
      ! consider only time steps which are also advection steps and
      ! choose the dynamic time step matching or following the time
      ! stamp string:

      l_adv_active = isCurrentEventActive(mtime_adv, mtime_dyn_date) .OR. &
        &            (mtime_dyn_date == mtime_end)

      IF ((mtime_dyn_date >= mtime_dom_start) .AND.  &
        & (mtime_dyn_date >= mtime_cur_date ) .AND.  &
        & l_adv_active ) THEN

        CALL datetimeToString(mtime_dyn_date, dtime_string)
        ilist = ilist + 1
        result_exactdate(ilist) = dtime_string
        result_steps(ilist)     = sim_step
        ! debugging output
        IF (ldebug)  WRITE (0,*) ilist, ": ", sim_step, " -> ", TRIM(dtime_string)

        IF (ilist >= nstrings) EXIT EVENT_LOOP
        CALL deallocateDatetime(mtime_cur_date)
        mtime_cur_date => newDatetime(date_string(ilist+1))

        ! consistency check: make sure that the output event
        ! intervals are not too small:
        IF (mtime_dyn_date >= mtime_cur_date) THEN
          WRITE (0,*) "ilist    = ", ilist
          WRITE (0,*) "dyn_date = ", TRIM(dtime_string)
          WRITE (0,*) "cur_date = ", TRIM(date_string(ilist+1))
          CALL finish(routine, "Output interval chosen too small!")
        END IF
      END IF
      
      ! increase simulation step
      sim_step = sim_step + 1
      mtime_dyn_date = mtime_dyn_date + delta

      ! check if domain is inactive:
      IF (mtime_dyn_date > mtime_dom_end)  EXIT EVENT_LOOP

      ! check if dynamic time step has reached the end of the simulation:
      l_dyn_active = isCurrentEventActive(mtime_dyn, mtime_dyn_date)
      IF (.NOT. l_dyn_active) EXIT EVENT_LOOP
    END DO EVENT_LOOP

    ! debugging output
    IF (ldebug)  WRITE (0,*) "result_exactdate: ", result_exactdate(1:ilist)

    ! clean up
    CALL deallocateEvent(mtime_dyn)
    CALL deallocateEvent(mtime_adv)
    CALL deallocateDatetime(mtime_dom_start)
    CALL deallocateDatetime(mtime_dom_end)
    CALL deallocateDatetime(mtime_cur_date)
    CALL deallocateDatetime(mtime_dyn_date)
    CALL deallocateDatetime(mtime_begin)
    CALL deallocateDatetime(mtime_end)
    CALL deallocateTimedelta(delta)
    CALL resetCalendar()
  END SUBROUTINE compute_matching_sim_steps


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
    &                                sim_step_info, fname_metadata)  RESULT(result_fnames)
    INTEGER,                INTENT(IN)    :: nstrings           !< no. of string to convert
    CHARACTER(len=*),       INTENT(IN)    :: date_string(:)     !< array of ISO 8601 time stamp strings
    INTEGER,                INTENT(IN)    :: sim_steps(:)       !< array of corresponding simulation steps
    TYPE(t_sim_step_info),  INTENT(IN)    :: sim_step_info      !< definitions: time step size, etc.
    TYPE(t_fname_metadata), INTENT(IN)    :: fname_metadata     !< additional meta-data for generating output filename
    TYPE(t_event_step_data) :: result_fnames(SIZE(date_string))
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::generate_output_filenames"
    INTEGER                             :: i, j, ifile, ipart
    CHARACTER(len=MAX_STRING_LEN)       :: cfilename
    TYPE (t_keyword_list), POINTER      :: keywords     => NULL()
    CHARACTER(len=MAX_STRING_LEN)       :: fname(nstrings)        ! list for duplicate check
    INTEGER                             :: ifname                 ! current length of "ifname"
    TYPE(datetime),  POINTER            :: file_end, step_date
    TYPE(timedelta), POINTER            :: delta
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: dtime_string

    ! ---------------------------------------------------
    ! prescribe, which file is used in which output step.
    ! ---------------------------------------------------
    ! a) user has specified "steps_per_file"
    ! b) user has specified "file_interval"
    IF (TRIM(fname_metadata%file_interval) == "") THEN
      ! case a): steps_per_file
      DO i=1,nstrings
        result_fnames(i)%jfile = (i-1)/fname_metadata%steps_per_file + 1 ! INTEGER division!
        result_fnames(i)%jpart = i - (result_fnames(i)%jfile-1)*fname_metadata%steps_per_file
      END DO
    ELSE
      ! case b): file_interval
      CALL setCalendar(PROLEPTIC_GREGORIAN)
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
    DO i=1,nstrings
      IF (.NOT. result_fnames(i)%l_open_file) THEN
        ! if no file is opened: filename is identical to last step
        result_fnames(i)%filename_string = result_fnames(i-1)%filename_string
        CYCLE
      END IF
      ! otherwise, generate filename
      CALL associate_keyword("<path>",            TRIM(model_base_dir),                                     keywords)
      CALL associate_keyword("<output_filename>", TRIM(fname_metadata%filename_pref),                       keywords)
      CALL associate_keyword("<physdom>",         TRIM(int2string(fname_metadata%phys_patch_id, "(i2.2)")), keywords)
      CALL associate_keyword("<levtype>",         TRIM(lev_type_str(fname_metadata%ilev_type)),             keywords)
      CALL associate_keyword("<levtype_l>",       TRIM(tolower(lev_type_str(fname_metadata%ilev_type))),    keywords)
      CALL associate_keyword("<jfile>",           TRIM(int2string(result_fnames(i)%jfile, "(i4.4)")),       keywords)
      CALL associate_keyword("<ddhhmmss>",        TRIM(date_string(i)),                                     keywords)
      cfilename = TRIM(with_keywords(keywords, fname_metadata%filename_format))
      IF(my_process_is_mpi_test()) THEN
        WRITE(result_fnames(i)%filename_string,'(a,"_TEST",a)') TRIM(cfilename),TRIM(fname_metadata%extn)
      ELSE
        WRITE(result_fnames(i)%filename_string,'(a,a)')         TRIM(cfilename),TRIM(fname_metadata%extn)
      ENDIF

      ! consistency check: test if the user has accidentally chosen a
      ! file name syntax which yields duplicate file names:
      ifname = 0
      DO j=1,ifname
        IF (result_fnames(i)%filename_string == fname(j)) THEN
          ! found a duplicate filename:
          CALL finish(routine, "Ambiguous output file name: '"//TRIM(fname(j))//"'")
        END IF
        ifname = ifname + 1
        fname(ifname) = result_fnames(i)%filename_string
      END DO
    END DO
  END FUNCTION generate_output_filenames

END MODULE mo_output_event_control
