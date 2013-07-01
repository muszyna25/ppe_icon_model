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
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_latbc_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_dom, MAX_CHAR_LENGTH
  USE mo_io_units,           ONLY: filename_max
  USE mo_util_string,        ONLY: t_keyword_list, MAX_STRING_LEN,   &
                                   associate_keyword, with_keywords, &
                                   int2string
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_datetime,           ONLY: t_datetime


  IMPLICIT NONE

  PUBLIC :: t_latbc_config, latbc_config, configure_latbc, generate_filename

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !>
  !!----------------------------------------------------------------------------
  !! Derived type containing control variables specific to the nonhydrostatic 
  !! atm model

  !------------------------------------------------------------------------
  TYPE t_latbc_config

    ! variables from namelist
    INTEGER                       :: itype_latbc      ! type of limited area boundary nudging
    REAL(wp)                      :: dtime_latbc      ! dt between two consequtive external latbc files
    CHARACTER(LEN=filename_max)   :: latbc_filename   ! prefix of latbc files
    CHARACTER(LEN=MAX_STRING_LEN) :: latbc_path       ! directory containing external latbc files

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
      "mo_latbc_config::configure_latbc"

    !----------------------------------------------------
    ! Sanity check and Prints
    !----------------------------------------------------

    IF(latbc_config%itype_latbc == 0) THEN

       WRITE(message_text,'(a)')'Internal lateral boundary condition (Using the boundary data from the previous time step).'
       CALL message(TRIM(routine),message_text)

    ELSE IF(latbc_config%itype_latbc == 1) THEN

       WRITE(message_text,'(a)')'External lateral boundary condition (Currently only IFS files supported).'
       CALL message(TRIM(routine),message_text)

    ELSE
       WRITE(message_text,'(a,i8)') 'Wrong lateral boundary condition mode:', latbc_config%itype_latbc
       CALL finish(TRIM(routine),message_text)
    END IF
  
  END SUBROUTINE configure_latbc
  !--------------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------------
  FUNCTION generate_filename(nroot, jlev, latbc_datetime) RESULT(result_str)
    INTEGER,          INTENT(IN)                :: nroot, jlev
    TYPE(t_datetime), INTENT(IN)                :: latbc_datetime
    CHARACTER(MAX_STRING_LEN)                   :: result_str

    ! Local variables
    CHARACTER(MAX_STRING_LEN)                   :: snroot, sjlev, stlev, slatbc_filename
    CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: &
      &  routine = 'mo_latbc_config::generate_filename:'
    
    WRITE(slatbc_filename,*) TRIM(latbc_config%latbc_path)//TRIM(latbc_config%latbc_filename)
    WRITE(snroot,'(i1)') nroot
    WRITE(sjlev,'(i2.2)') jlev
    WRITE(stlev,'(i4,4i2.2)') latbc_datetime%year, latbc_datetime%month, latbc_datetime%day, latbc_datetime%hour

    result_str = TRIM(slatbc_filename)//'R'//TRIM(snroot)//'B'//TRIM(sjlev)//'_DOM01_'//TRIM(stlev)//'.nc'
  END FUNCTION generate_filename
  !--------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
END MODULE mo_latbc_config
