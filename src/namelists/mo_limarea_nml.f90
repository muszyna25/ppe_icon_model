!>
!! Contains the setup of variables related to the boundary relaxation of limited
!! area models from an external dataset
!!
!! @author S. Brdar (DWD)
!!
!!
!! @par Revision History
!! Initial release by S. Brdar (2013-06-13)
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
MODULE mo_limarea_nml

  USE mo_kind,                ONLY: wp, i8
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_impl_constants,      ONLY: max_dom, MAX_CHAR_LENGTH
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist     , &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_limarea_config,      ONLY: latbc_config
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mtime,                  ONLY: MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN,            &
    &                               newTimedelta, deallocateTimedelta, OPERATOR(==),        &
    &                               getPTStringFromMS, timedelta,                           &
    &                               getTotalMilliSecondsTimeDelta, datetime, newDatetime,   &
    &                               deallocateDatetime

  IMPLICIT NONE
  PRIVATE
  PUBLIC read_limarea_namelist

  !------------------------------------------------------------------------
  ! Namelist variables
  !------------------------------------------------------------------------
  INTEGER                         :: itype_latbc    ! type of limited area boundary nudging
  REAL(wp)                        :: dtime_latbc    ! dt between two consequtive external latbc files
  CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: dt_latbc       ! dt between two consequtive external latbc files    
  INTEGER                         :: nlev_latbc     ! number of vertical levels in boundary data
  CHARACTER(LEN=filename_max)     :: latbc_filename ! prefix of latbc files
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: latbc_path     ! directory containing external latbc files

  NAMELIST /limarea_nml/ itype_latbc, dtime_latbc, nlev_latbc, latbc_filename, latbc_path, dt_latbc

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_limarea_nml'

CONTAINS
  !>
  !!
  SUBROUTINE read_limarea_namelist( filename )
    CHARACTER(LEN=*), INTENT(IN) :: filename
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::configure_latbc"

    INTEGER                              :: istat, funit
    INTEGER                              :: iunit, errno
    REAL(wp)                             :: dtime_latbc_in_ms
    TYPE(timedelta), POINTER             :: td
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: dt_latbc_str
    TYPE(datetime), POINTER              :: base_dt

    !------------------------------------------------------------
    ! Default settings
    !------------------------------------------------------------
    itype_latbc      = 0
    dtime_latbc      = -1._wp
    dt_latbc         = ''
    nlev_latbc       = 0
    latbc_filename   = "prepiconR<nroot>B<jlev>_<y><m><d><h>.nc"
    latbc_path       = "./"

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('limarea_nml')
      READ(funit,NML=limarea_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('limarea_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, limarea_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, limarea_nml)
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, limarea_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! Fill the configuration state
    !----------------------------------------------------
    latbc_config% itype_latbc     = itype_latbc
    latbc_config% nlev_in         = nlev_latbc
    latbc_config% latbc_filename  = latbc_filename
    latbc_config% latbc_path      = TRIM(latbc_path)//'/'

    ! There exist to alternative ways to set the update interval for
    ! lateral bc data. If both parameters are used, we test for
    ! inconsistency. Otherwise we translate the value of one setting
    ! to the other one.

    IF (LEN_TRIM(dt_latbc) == 0) THEN
      dtime_latbc_in_ms = 1000._wp * dtime_latbc
      CALL getPTStringFromMS(NINT(dtime_latbc_in_ms,i8), latbc_config%dt_latbc)
      latbc_config%dtime_latbc_mtime => newTimedelta(latbc_config%dt_latbc, errno)
      IF (errno /= 0)  CALL finish(routine, "Error in initialization of dtime_latbc time delta.")
    ELSE
      latbc_config%dtime_latbc_mtime => newTimedelta(latbc_config%dt_latbc, errno)
      IF (errno /= 0)  CALL finish(routine, "Error in initialization of dtime_latbc time delta.")

      IF (dtime_latbc > 0.) THEN
        ! test for inconsistency
        dtime_latbc_in_ms = 1000._wp * dtime_latbc
        CALL getPTStringFromMS(NINT(dtime_latbc_in_ms,i8), dt_latbc_str)
        td => newTimedelta(dt_latbc_str, errno)
        IF (errno /= 0)  CALL finish(routine, "Error in initialization of dtime_latbc time delta.")
        IF (.NOT. (td == latbc_config%dtime_latbc_mtime)) THEN
          CALL finish(routine, "Inconsistent setting of dtime_latbc time delta.")
        END IF
        CALL deallocateTimedelta(td)        
      ELSE
        base_dt => newDatetime("2011-01-01T00:00:00Z", errno)
        IF (errno /= 0)  CALL finish(routine, "Error in setup of reference datetime.")
        dtime_latbc = getTotalMilliSecondsTimeDelta(latbc_config%dtime_latbc_mtime, base_dt)/1000._wp
        CALL deallocateDatetime(base_dt)
      END IF
    END IF
    latbc_config%dtime_latbc       = dtime_latbc
    latbc_config%dt_latbc          = TRIM(dt_latbc)

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio()) THEN
      funit = open_tmpfile()
      WRITE(funit,NML=limarea_nml)
      CALL store_and_close_namelist(funit, 'limarea_nml')
    ENDIF

    !-----------------------------------------------------
    !write the contents of the namelist to an ASCII file
    !-----------------------------------------------------
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=limarea_nml)

  END SUBROUTINE read_limarea_namelist

END MODULE mo_limarea_nml
