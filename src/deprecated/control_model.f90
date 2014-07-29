!>
!! @page pagecontrolmodelf901 ICON Master Program
!!
!!
!! @author
!!     Leonidas Linardakis, Luis KOrnblueh
!!     (MPI-M)
!!
!! @date 2013-06-04
!!

!>
!! This is the master progam of the ICON model.
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
PROGRAM control_model

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
  USE mo_master_control,      ONLY: init_master_control,  &
    & get_my_namelist_filename, get_my_process_type,      &
    & testbed_process,  atmo_process, ocean_process!, radiation_process
  USE mo_time_config,         ONLY: restart_experiment
  USE mo_util_signal
  USE mo_util_vcs,            ONLY: util_repository_url, &
       &                            util_branch_name,    &
       &                            util_revision_key

#ifndef __NO_ICON_OCEAN__
  USE mo_ocean_model,         ONLY: ocean_model
#endif

#ifndef __NO_ICON_TESTBED__
  USE mo_icon_testbed,        ONLY: icon_testbed
#endif

#ifndef __NO_ICON_ATMO__
  USE mo_atmo_model,          ONLY: atmo_model
!   USE mo_radiation_model,     ONLY: radiation_model
#endif

#ifdef MESSY
  USE messy_main_constants_mem, ONLY: modver
  USE messy_main_compilerinfo_mem, ONLY: compiler_version, compiler_call &
    &   , compiler_flags                                                 &
    &   , compiler_cppdefs, compiler_includes
#endif

  IMPLICIT NONE

  INTEGER    :: master_control_status
  
  INTEGER    :: my_process_component

  CHARACTER(len=filename_max) :: my_namelist_filename
  CHARACTER(len=filename_max) :: master_namelist_filename="icon_master.namelist"

  CHARACTER(len=256) :: repository = ''
  CHARACTER(len=256) :: branch = ''
  CHARACTER(len=256) :: revision = ''
  
  INTEGER :: nlen

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  TYPE(ieee_status_type) :: saved_fpscr
  LOGICAL :: halting_mode,  current_flags(size(ieee_all))
  REAL(wp) :: r
#endif
#endif

#if defined (__xlC__)
  INTEGER :: core_dump_flag
  INTEGER :: signals(1)
  INTEGER :: iret
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
  ! Initialize MPI, this should aleays be the first call
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

  !-------------------------------------------------------------------
  ! PRINT "DEPRECTATED" WARNING
  !-------------------------------------------------------------------
  IF (my_process_is_global_root()) THEN

    WRITE (0,*) "  -----------------------------------------------"
    WRITE (0,*) "    ___  ___ ___ ___ ___ ___   _ _____ ___ ___   "
    WRITE (0,*) "   |   \| __| _ \ _ \ __/ __| /_\_   _| __|   \  "
    WRITE (0,*) "   | |) | _||  _/   / _| (__ / _ \| | | _|| |) | "
    WRITE (0,*) "   |___/|___|_| |_|_\___\___/_/ \_\_| |___|___/  "
    WRITE (0,*) "                                                 "
    WRITE (0,*) "  -----------------------------------------------"
    WRITE (0,*) "                                                 "
    WRITE (0,*) "       ICON'S MAIN PROGRAM HAS BEEN RENAMED!     "
    WRITE (0,*) "                                                 "
    WRITE (0,*) "             control_model  ->  icon             "
    WRITE (0,*) "                                                 "
    WRITE (0,*) "     PLEASE CHANGE YOUR RUN SCRIPT ACCORDINGLY.  "
    WRITE (0,*) "                                                 "
    WRITE (0,*) "     In future versions the control_model        "
    WRITE (0,*) "     executable will be REMOVED!                 "
    WRITE (0,*) "                      	                          "
    WRITE (0,*) "  -----------------------------------------------"

  END IF
  !-------------------------------------------------------------------


  !-------------------------------------------------------------------
  ! Print VCS version
  IF (my_process_is_global_root()) THEN

    nlen = 256
    call util_repository_url(repository, nlen)
    nlen = 256
    call util_branch_name(branch, nlen)
    nlen = 256
    call util_revision_key(revision, nlen)

    WRITE(message_text,'(a,a)') 'Repository: ', TRIM(repository)
    CALL message('',message_text)
    WRITE(message_text,'(a,a)') 'Branch    : ', TRIM(branch)
    CALL message('',message_text)
    WRITE(message_text,'(a,a)') 'Revision  : ', TRIM(revision)
    CALL message('',message_text)
    CALL message('','')

  END IF

  !-------------------------------------------------------------------
  ! Initialize the master control

  master_control_status = init_master_control(TRIM(master_namelist_filename))
  
  my_namelist_filename = get_my_namelist_filename()
  my_process_component = get_my_process_type()

!   write(0,*) 'control model. my_process_component=',my_process_component
  
  SELECT CASE (my_process_component)

#ifndef __NO_ICON_ATMO__
  CASE (atmo_process)
    CALL atmo_model(my_namelist_filename,TRIM(master_namelist_filename))
#endif

!   CASE (radiation_process)
!     CALL radiation_model(my_namelist_filename, TRIM(master_namelist_filename))

#ifndef __NO_ICON_TESTBED__
  CASE (testbed_process)
    CALL icon_testbed(my_namelist_filename, TRIM(master_namelist_filename))
#endif

#ifndef __NO_ICON_OCEAN__
  CASE (ocean_process)
    CALL ocean_model(my_namelist_filename, TRIM(master_namelist_filename))
#endif

  CASE default
    CALL finish("control_model","my_process_component is unkown")
    
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

END PROGRAM control_model
