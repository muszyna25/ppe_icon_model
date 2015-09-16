!>
!! This is the master program of the ICON model.
!!
!!
!! @par Revision History
!!   
!! 
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
PROGRAM icon

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  USE, INTRINSIC :: ieee_features
  USE, INTRINSIC :: ieee_arithmetic
  USE, INTRINSIC :: ieee_exceptions

  USE mo_kind, ONLY: wp
#endif
#endif
#if defined (__INTEL_COMPILER) && ! defined (VARLIST_INITIZIALIZE_WITH_NAN)
  USE, INTRINSIC :: ieee_arithmetic
#endif
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_io_units,            ONLY: filename_max
  USE mo_mpi,                 ONLY: start_mpi , stop_mpi, my_process_is_global_root
  USE mo_master_control,      ONLY: init_master_control,                                &
    &                               get_my_namelist_filename, get_my_process_type,      &
    &                               atmo_process, ocean_process, testbed_process
  USE mo_time_config,         ONLY: restart_experiment
  USE mo_util_signal
  USE mo_util_sysinfo,        ONLY: util_user_name, util_os_system, util_node_name
  USE mo_util_vcs,            ONLY: util_repository_url,                                &
    &                               util_branch_name,                                   &
    &                               util_revision_key

#ifndef __NO_ICON_OCEAN__
  USE mo_ocean_model,         ONLY: ocean_model
#endif

#ifndef __NO_ICON_TESTBED__
  USE mo_icon_testbed,        ONLY: icon_testbed
#endif

#ifndef __NO_ICON_ATMO__
  USE mo_atmo_model,          ONLY: atmo_model
#endif

  USE mo_cdi,                 ONLY: gribapiLibraryVersion
  USE mo_cf_convention          ! We need all ? 
#ifdef _MTIME_DEBUG
  USE mo_kind, ONLY: wp
  USE mtime,   ONLY: timedelta, newTimedelta, deallocateTimedelta,  &
    &                timedeltaToString, ASSIGNMENT(=), OPERATOR(*), &
    &                MAX_TIMEDELTA_STR_LEN, setCalendar,            &
    &                PROLEPTIC_GREGORIAN, datetime,                 &
    &                datetimetostring, newDatetime, OPERATOR(+),    &
    &                deallocateDatetime, MAX_DATETIME_STR_LEN
#endif
  
  IMPLICIT NONE

  INTEGER                     :: master_control_status, my_process_component
  CHARACTER(len=filename_max) :: my_namelist_filename
  CHARACTER(len=filename_max) :: master_namelist_filename="icon_master.namelist"

  CHARACTER(len=256) :: repository  = ''
  CHARACTER(len=256) :: branch      = ''
  CHARACTER(len=256) :: revision    = ''
  CHARACTER(len=256) :: executable  = ''
  CHARACTER(len=256) :: user_name   = ''
  CHARACTER(len=256) :: os_name     = ''
  CHARACTER(len=256) :: host_name   = ''
  CHARACTER(len=256) :: tmp_string  = ''
  CHARACTER(len=  8) :: date_string = ''
  CHARACTER(len= 10) :: time_string = ''

  INTEGER :: nlen, nlena, nlenb, nlenc, nlend
  INTEGER :: grb_major_version, grb_minor_version, grb_revision_version

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  TYPE(ieee_status_type)      :: saved_fpscr
  LOGICAL                     :: halting_mode,  current_flags(size(ieee_all))
  REAL(wp)                    :: r
#endif
#endif

#if defined (__xlC__)
  INTEGER                     :: core_dump_flag
  INTEGER                     :: signals(1)
  INTEGER                     :: iret
#endif

#ifdef _MTIME_DEBUG
  TYPE(timedelta),  POINTER             :: mtime_td
  TYPE(datetime),   POINTER             :: mtime_date
  CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  :: td_string
  CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: dstring
#endif

!--------------------------------------------------------------------

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  CALL ieee_get_status(saved_fpscr)
  CALL ieee_set_halting_mode(ieee_all, .TRUE.)
#endif
#endif

#if defined (__INTEL_COMPILER) && ! defined (VARLIST_INITIZIALIZE_WITH_NAN)
  ! Important on Intel: disable underflow exceptions:
  CALL ieee_set_halting_mode(ieee_underflow, .FALSE.)
#endif

#if defined (__SX__)
  ! sxf90 is not Fortran standard compliant, use vendor extension:

  ! export environment variable F_NORCW=65535
  ! (this SX environment variable specifies that a control record is
  !  not added/expected, s.t. the file content can be treated like a
  !  stream of characters)
  CALL putenv ("F_NORCW=65535")
#endif

  !-------------------------------------------------------------------
  ! Initialize MPI, this should always be the first call
  CALL start_mpi('ICON')
  
  !-------------------------------------------------------------------
  !set up signal trapping on IBM: export USE_SIGNAL_HANDLING=yes

#if defined (__xlC__) 
  core_dump_flag = 0
  signals(1)     = 0

  iret = signal_trap(core_dump_flag, signals)

  IF (iret == -2) THEN
    CALL message('', 'Signal trapping disabled by environment')
  ELSE IF (iret == -1) THEN
     WRITE(message_text,'(a,i0)') 'Error: ', iret
     CALL message('', message_text)
  ELSE IF (iret == 0) THEN
    CALL message('', 'FPE trapping is not set')
  ELSE
    WRITE(message_text,'(a,i0)') 'FPE trapping mode =', iret
    CALL message('', message_text)
  END IF
#endif
#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  WRITE(message_text,'(a,l1)') ' IEEE standard supported: ', ieee_support_standard(r)
  CALL message('', message_text)
#endif
#endif

  ! --- VCS version
  !
  ! When using SVN, the revision number is printed as it is
  ! generated by the "config/pvcs.pl" command-line tool during the
  ! make process. This revision number is stored in the file
  ! build/<config.guess triplet>/src/version.c.
  !
  ! Comments on the revision number format: Consider, for example
  !
  !     21000:21099M
  !
  ! "M"           : The working copy has local modifications.
  !
  ! 21000:21099 :  This is a mixed-revision working copy. Not all
  !                files have the same revision number and revision
  !                numbers are within the given range. This usually
  !                happens when a developer makes a checkout of
  !                revision 21000, then changes a single file and
  !                commits, resulting in revision 21099 for this
  !                file. As long as the user does not "svn update",
  !                his working copy has mixed revisions.

  nlen = 256
  call util_repository_url(repository, nlen)
  nlen = 256
  call util_branch_name(branch, nlen)
  nlen = 256
  call util_revision_key(revision, nlen)
  
  CALL get_command_argument(0, executable, nlend)
  CALL date_and_time(date_string, time_string)

  tmp_string = ''
  CALL util_os_system (tmp_string, nlena)
  os_name = tmp_string(1:nlena)

  tmp_string = ''
  CALL util_user_name (tmp_string, nlenb)
  user_name = tmp_string(1:nlenb)

  tmp_string = ''
  CALL util_node_name (tmp_string, nlenc)
  host_name = tmp_string(1:nlenc)

  CALL set_cf_global(title       = 'ICON simulation',                                                    &
       &             institution = 'Max Planck Institute for Meteorology/Deutscher Wetterdienst',        &
       &             source      = TRIM(repository)//'@'//TRIM(revision),                                &
       &             history     = executable(1:nlend)//' at '//date_string(1:8)//' '//time_string(1:6), &
       &             references  = 'see MPIM/DWD publications',                                          &
       &             comment     = TRIM(user_name)//' on '//TRIM(host_name)//' ('//TRIM(os_name)//')')
 
  !-------------------------------------------------------------------
  IF (my_process_is_global_root()) THEN

    WRITE(message_text,'(a,a)') 'Repository: ', TRIM(repository)
    CALL message('',message_text)
    WRITE(message_text,'(a,a)') 'Branch    : ', TRIM(branch)
    CALL message('',message_text)
    WRITE(message_text,'(a,a)') 'Revision  : ', TRIM(revision)
    CALL message('',message_text)
    WRITE(message_text,'(a,a)') 'Executable: ', TRIM(cf_global_info%history)
    CALL message('',message_text)
    WRITE(message_text,'(a,a)') 'User      : ', TRIM(cf_global_info%comment)
    CALL message('',message_text)

    ! --- Ask CDI for the GRIB version in use
    CALL gribapiLibraryVersion(grb_major_version, grb_minor_version, grb_revision_version)
    IF (grb_major_version > 0) THEN
      WRITE(message_text,'(a,i0,a,i0,a,i0)') 'GRIB-API  : ', grb_major_version, ".", grb_minor_version, ".",&
        &                                    grb_revision_version
      CALL message('',message_text)
    END IF

    ! --- Cray compiler: print compiler version in use:
#ifdef _CRAYFTN
    WRITE(message_text,'(a,a)') 'Compiler  : ', _RELEASE_STRING
    CALL message('',message_text)
#endif

    CALL message('','')
  END IF

#ifdef _MTIME_DEBUG
  WRITE(message_text,'(a,a)') 'executable running with MTIME_DEBUG enabled.'
  CALL message('',message_text)

  !-------------------------------------------------------------------
  ! a small test of mtime functionality
  !-------------------------------------------------------------------
  IF (my_process_is_global_root()) THEN
    CALL setCalendar(PROLEPTIC_GREGORIAN)
    
    mtime_td => newTimedelta("PT1H1M1S")
    mtime_td = mtime_td * 0.3_wp
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0,*) "PT1H1M1S * 0.3 = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT1H1M1S")
    mtime_td = mtime_td * 0.5_wp
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0,*) "PT1H1M1S * 0.5 = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT1H1M1S")
    mtime_td = mtime_td * 1.5_wp
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0,*) "PT1H1M1S * 1.5 = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT1H1M1S")
    mtime_td = mtime_td * 2.0_wp
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0,*) "PT1H1M1S * 2.0 = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT98765S")
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0,*) "PT98765S = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT987654321S")
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0,*) "PT987654321S = ", TRIM(td_string)

    mtime_date => newDatetime("1970-01-01T00:00:00")
    CALL datetimetostring(mtime_date, dstring)
    WRITE (0,*) "1970-01-01T00:00:00 = ", TRIM(dstring)
    mtime_date = mtime_date + mtime_td
    CALL datetimetostring(mtime_date, dstring)
    WRITE (0,*) "1970-01-01T00:00:00 + PT987654321S = ", TRIM(dstring)
    CALL deallocateDatetime(mtime_date)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT10000000S")
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0,*) "PT10000000S = ", TRIM(td_string)

    mtime_date => newDatetime("2014-06-01T00:00:00")
    CALL datetimetostring(mtime_date, dstring)
    WRITE (0,*) "2014-06-01T00:00:00 = ", TRIM(dstring)
    mtime_date = mtime_date + mtime_td
    CALL datetimetostring(mtime_date, dstring)
    WRITE (0,*) "2014-06-01T00:00:00 + PT10000000S = ", TRIM(dstring)
    CALL deallocateDatetime(mtime_date)
    CALL deallocateTimedelta(mtime_td)
  END IF
#endif

  !-------------------------------------------------------------------
  ! Initialize the master control

  master_control_status = init_master_control(TRIM(master_namelist_filename))
  
  my_namelist_filename = get_my_namelist_filename()
  my_process_component = get_my_process_type()

  SELECT CASE (my_process_component)

#ifndef __NO_ICON_ATMO__
  CASE (atmo_process)
    CALL atmo_model  (my_namelist_filename, TRIM(master_namelist_filename))
#endif

#ifndef __NO_ICON_OCEAN__
  CASE (ocean_process)
    CALL ocean_model (my_namelist_filename, TRIM(master_namelist_filename))
#endif

#ifndef __NO_ICON_TESTBED__
  CASE (testbed_process)
    CALL icon_testbed(my_namelist_filename, TRIM(master_namelist_filename))
#endif

  CASE default
    CALL finish("icon","my_process_component is unkown")
    
  END SELECT
      
  ! write the control.status file
  IF (my_process_is_global_root()) THEN
    OPEN (500, FILE="finish.status")
    IF (restart_experiment) THEN
      WRITE(500,*) "RESTART"
    ELSE
      WRITE(500,*) "OK"
    ENDIF
    CLOSE(500)    
  END IF

  ! Shut down MPI
  !
  CALL stop_mpi

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  CALL ieee_set_status(saved_fpscr)
#endif
#endif

END PROGRAM icon
