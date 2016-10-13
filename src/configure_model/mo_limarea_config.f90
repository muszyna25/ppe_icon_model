!>
!! Contains the setup of variables related to the lateral boundary
!! condition for limited area models
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
MODULE mo_limarea_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_dom, MAX_CHAR_LENGTH
  USE mo_io_units,           ONLY: filename_max
  USE mo_util_string,        ONLY: t_keyword_list,                   &
                                   associate_keyword, with_keywords, &
                                   int2string
  USE mo_exception,          ONLY: message, message_text, finish
  USE mtime,                 ONLY: MAX_TIMEDELTA_STR_LEN, datetime,  &
    &                              timedelta

  IMPLICIT NONE

  PUBLIC :: t_latbc_config, latbc_config, configure_latbc, generate_filename, &
    &       generate_filename_mtime 
  !>
  !!----------------------------------------------------------------------------
  !! Derived type containing control variables specific to the nonhydrostatic 
  !! atm model

  !------------------------------------------------------------------------
  TYPE t_latbc_config

    ! variables from namelist
    INTEGER                              :: itype_latbc       ! type of limited area boundary nudging
    INTEGER                              :: nlev_in           
    CHARACTER(LEN=filename_max)          :: latbc_filename    ! prefix of latbc files
    CHARACTER(LEN=MAX_CHAR_LENGTH)       :: latbc_path        ! directory containing external latbc files
    REAL(wp)                             :: lc1, lc2          ! linear interpolation coefficients
                                                             
    REAL(wp)                             :: dtime_latbc       ! dt between two consequtive external latbc files
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: dt_latbc          ! ISO string (dt_latbc)
    TYPE(timedelta), POINTER             :: dtime_latbc_mtime ! dt between two consequtive external latbc files
  END TYPE t_latbc_config
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  TYPE(t_latbc_config), TARGET :: latbc_config
  !!----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE configure_latbc()
  !--------------------------------------------------------------------------------------
  !  Set up parameters 
  !--------------------------------------------------------------------------------------
    CHARACTER(*), PARAMETER :: routine = &
      "mo_limarea_config::configure_latbc"

    !----------------------------------------------------
    ! Sanity check and Prints
    !----------------------------------------------------

    IF (latbc_config%itype_latbc == 0) THEN

       WRITE(message_text,'(a)')'Lateral boundary nudging using the initial boundary data.'
       CALL message(TRIM(routine),message_text)

    ELSE IF (latbc_config%itype_latbc == 1) THEN

       WRITE(message_text,'(a)')'Lateral boundary condition using the IFS or COSMO-DE boundary data.'
       CALL message(TRIM(routine),message_text)

    ELSE IF (latbc_config%itype_latbc == 2) THEN

       WRITE(message_text,'(a)')'Lateral boundary condition using the ICON global boundary data.'
       CALL message(TRIM(routine),message_text)

    ELSE

       WRITE(message_text,'(a,i8)') 'Wrong lateral boundary condition mode:', latbc_config%itype_latbc
       CALL finish(TRIM(routine),message_text)

    END IF
  
  END SUBROUTINE configure_latbc
  !--------------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------------
  FUNCTION generate_filename(nroot, jlev, latbc_mtime) RESULT(result_str)
    INTEGER,          INTENT(IN)                :: nroot, jlev
    TYPE(datetime),   INTENT(IN)                :: latbc_mtime
    CHARACTER(MAX_CHAR_LENGTH )                 :: result_str

    ! Local variables
    TYPE (t_keyword_list), POINTER              :: keywords => NULL()
    CHARACTER(MAX_CHAR_LENGTH)                  :: str
    CHARACTER(MAX_CHAR_LENGTH), PARAMETER       :: &
      &  routine = 'mo_limarea_config::generate_filename_mtime:'
    
    WRITE(str,'(i4)')   latbc_mtime%date%year
    CALL associate_keyword("<y>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%date%month
    CALL associate_keyword("<m>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%date%day
    CALL associate_keyword("<d>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%hour
    CALL associate_keyword("<h>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%minute
    CALL associate_keyword("<min>",       TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%second !FLOOR(latbc_mtime%time%second)
    CALL associate_keyword("<sec>",       TRIM(str),                        keywords)
      
    CALL associate_keyword("<nroot>",     TRIM(int2string(nroot,'(i1)')),   keywords)
    CALL associate_keyword("<nroot0>",    TRIM(int2string(nroot,'(i2.2)')), keywords)
    CALL associate_keyword("<jlev>",      TRIM(int2string(jlev, '(i2.2)')), keywords)
    CALL associate_keyword("<dom>",       TRIM(int2string(1,'(i2.2)')),     keywords)

    ! replace keywords in latbc_filename
    result_str = TRIM(with_keywords(keywords, TRIM(latbc_config%latbc_filename)))
  END FUNCTION generate_filename
  !--------------------------------------------------------------------------------------
  FUNCTION generate_filename_mtime(nroot, jlev, latbc_mtime) RESULT(result_str)
    INTEGER,          INTENT(IN)                :: nroot, jlev
    TYPE(datetime),   INTENT(IN)                :: latbc_mtime
    CHARACTER(MAX_CHAR_LENGTH )                 :: result_str

    ! Local variables
    TYPE (t_keyword_list), POINTER              :: keywords => NULL()
    CHARACTER(MAX_CHAR_LENGTH)                  :: str
    CHARACTER(MAX_CHAR_LENGTH), PARAMETER       :: &
      &  routine = 'mo_limarea_config::generate_filename_mtime:'
    
    WRITE(str,'(i4)')   latbc_mtime%date%year
    CALL associate_keyword("<y>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%date%month
    CALL associate_keyword("<m>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%date%day
    CALL associate_keyword("<d>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%hour
    CALL associate_keyword("<h>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%minute
    CALL associate_keyword("<min>",       TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%second !FLOOR(latbc_mtime%time%second)
    CALL associate_keyword("<sec>",       TRIM(str),                        keywords)
      
    CALL associate_keyword("<nroot>",     TRIM(int2string(nroot,'(i1)')),   keywords)
    CALL associate_keyword("<nroot0>",    TRIM(int2string(nroot,'(i2.2)')), keywords)
    CALL associate_keyword("<jlev>",      TRIM(int2string(jlev, '(i2.2)')), keywords)
    CALL associate_keyword("<dom>",       TRIM(int2string(1,'(i2.2)')),     keywords)

    ! replace keywords in latbc_filename
    result_str = TRIM(with_keywords(keywords, TRIM(latbc_config%latbc_filename)))
  END FUNCTION generate_filename_mtime
  !--------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
END MODULE mo_limarea_config
