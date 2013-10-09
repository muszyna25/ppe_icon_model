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
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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
  USE mo_mpi,                 ONLY: start_mpi , p_stop, my_process_is_global_root
  USE mo_master_control,      ONLY: init_master_control,  &
    & get_my_namelist_filename, get_my_process_type,      &
    & testbed_process,  atmo_process, ocean_process!, radiation_process
  USE mo_time_config,        ONLY: restart_experiment
  USE mo_util_signal
  USE mo_util_svn,           ONLY: printSVNVersion

  USE mo_ocean_model,         ONLY: ocean_model
  USE mo_icon_testbed,        ONLY: icon_testbed

#ifndef __ICON_OCEAN_ONLY__
  USE mo_atmo_model,          ONLY: atmo_model
!   USE mo_radiation_model,     ONLY: radiation_model
#endif

  IMPLICIT NONE

  INTEGER    :: master_control_status
  
  INTEGER    :: my_process_component

  CHARACTER(len=filename_max) :: my_namelist_filename
  CHARACTER(len=filename_max) :: master_namelist_filename="icon_master.namelist"

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


  !declaration of OpenMP Runtime Library Routines:
! !$  INTEGER omp_get_max_threads
! !!$  INTEGER omp_get_num_threads
! !$  INTEGER omp_get_num_procs
! !!$  INTEGER omp_get_thread_num
! !$  LOGICAL omp_get_dynamic
! 
! !$  INTEGER :: max_threads_omp, num_procs_omp
! !$  LOGICAL :: ldynamic_omp
! 
!   !--------------------------------------------------------------------
!   !BOC
! 
!   !print out some information about OpenMP parallelization
! !$  max_threads_omp  = omp_get_max_threads()
! !$  num_procs_omp    = omp_get_num_procs()
! !$  ldynamic_omp     = omp_get_dynamic()
! !$  WRITE(message_text,'(A,I3,A,I3)')                &
! !$    & "OpenMP:  MAX_THREADS = ", max_threads_omp,  &
! !$    & ",  NUM_PROCS = ", num_procs_omp
! !$  CALL message('control_model',message_text)
! !$  WRITE(message_text,'(A,L3)')  &
! !$    & "OpenMP:  DYNAMIC = ", ldynamic_omp
! !$  CALL message('control_model',message_text)

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
  ! Print SVN version
  IF (my_process_is_global_root()) THEN
    CALL printSVNVersion() ! (... if defined)
  END IF

  !-------------------------------------------------------------------
  ! Initialize the master control

  master_control_status = init_master_control(TRIM(master_namelist_filename))
  
  my_namelist_filename = get_my_namelist_filename()
  my_process_component = get_my_process_type()

!   write(0,*) 'control model. my_process_component=',my_process_component
  
  SELECT CASE (my_process_component)

#ifndef __ICON_OCEAN_ONLY__
  CASE (atmo_process)
    CALL atmo_model(my_namelist_filename,TRIM(master_namelist_filename))


!   CASE (radiation_process)
!     CALL radiation_model(my_namelist_filename, TRIM(master_namelist_filename))
#endif

  CASE (testbed_process)
    CALL icon_testbed(my_namelist_filename, TRIM(master_namelist_filename))

  CASE (ocean_process)
    CALL ocean_model(my_namelist_filename, TRIM(master_namelist_filename))

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
  CALL p_stop

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  CALL ieee_set_status(saved_fpscr)
#endif
#endif
END PROGRAM control_model

