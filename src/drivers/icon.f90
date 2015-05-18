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
    &                               testbed_process,  atmo_process, ocean_process
  USE mo_time_config,         ONLY: restart_experiment
  USE mo_util_signal
  USE mo_util_vcs,            ONLY: util_repository_url,                                &
       &                            util_branch_name,                                   &
       &                            util_revision_key

#ifndef __NO_ICON_OCEAN__
  USE mo_ocean_model,         ONLY: ocean_model
#endif

#ifndef __NO_ICON_TESTBED__
  USE mo_icon_testbed,        ONLY: icon_testbed
#endif

#ifndef __NO_ICON_ATMO__
  USE mo_atmo_model,          ONLY: atmo_model
#endif

  USE mo_cdi_constants          ! We need all  

  USE mo_expression, ONLY: expression, real_kind ! DEVELOPMENT

  IMPLICIT NONE

  INTEGER                     :: master_control_status, my_process_component, nlen
  CHARACTER(len=filename_max) :: my_namelist_filename
  CHARACTER(len=filename_max) :: master_namelist_filename="icon_master.namelist"
  CHARACTER(len=256)          :: repository = ''
  CHARACTER(len=256)          :: branch = ''
  CHARACTER(len=256)          :: revision = ''
  integer                     :: grb_major_version, grb_minor_version, grb_revision_version

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

  TYPE(expression)          :: formula                   ! DEVELOPMENT
  CLASS(*), POINTER         :: val_2D(:,:)   => NULL()   ! DEVELOPMENT
  REAL(real_kind), TARGET   :: z_sfc(2,2)                ! DEVELOPMENT

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
  IF (my_process_is_global_root()) THEN

    ! --- Print VCS version
    !
    ! When using SVN, the revision number is printed as it is
    ! generated by the "svnversion" command-line tool during the make
    ! process. This revision number is stored in the file
    ! build/.../src/version.c.
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
    !
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

  !-------------------------------------------------------------------
  ! Initialize the master control

  master_control_status = init_master_control(TRIM(master_namelist_filename))
  
  my_namelist_filename = get_my_namelist_filename()
  my_process_component = get_my_process_type()

  SELECT CASE (my_process_component)

#ifndef __NO_ICON_ATMO__
  CASE (atmo_process)
    IF (my_process_is_global_root()) THEN

      ! DEVELOPMENT ---------------------------------------------

      WRITE (0,*) "DEVELOPMENT ---------------------------------------------"
      ! create some dummy data
      z_sfc(1,:) = (/ 1, 2 /)
      z_sfc(2,:) = (/ 3, 4 /)
      
      ! --- scalar example:
      formula = expression("[z_sfc] + sin(45*pi/180.) * 10 + 5")
      CALL formula%link("z_sfc", z_sfc)
      CALL formula%evaluate(val_2D)
      SELECT TYPE(val_2D)
      TYPE is (REAL(real_kind))
        WRITE (0,*) "formula: result = ", val_2D
      END SELECT
      DEALLOCATE(val_2D)
      CALL formula%finalize()
      WRITE (0,*) "END DEVELOPMENT -----------------------------------------"
    END IF

    ! END DEVELOPMENT -----------------------------------------

    CALL atmo_model(my_namelist_filename,TRIM(master_namelist_filename))
#endif

#ifndef __NO_ICON_TESTBED__
  CASE (testbed_process)
    CALL icon_testbed(my_namelist_filename, TRIM(master_namelist_filename))
#endif

#ifndef __NO_ICON_OCEAN__
  CASE (ocean_process)
    CALL ocean_model(my_namelist_filename, TRIM(master_namelist_filename))
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
