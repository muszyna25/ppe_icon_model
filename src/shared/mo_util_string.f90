!>
!!    This module holds string conversion utilities.
!!
!!
!! @par Revision History
!!    Original version from the ECHAM5 model (version 5.3.01).
!!    Modified to the ProTeX-style by Hui Wan, MPI-M (2007-01-17)
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_util_string
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

! !PUBLIC ENTITIES
  PUBLIC :: tolower        ! Conversion   : 'ABCXYZ' -> 'abcxyz'
  PUBLIC :: toupper        ! Conversion   : 'abcxyz' -> 'ABCXYZ'
  PUBLIC :: char2          ! Conversion   : INTEGER  -> CHAR (LEN=2)
  PUBLIC :: separator      ! Format string: (/"-----...-----"/)

! !DEFINED PARAMETERS
  CHARACTER(len=*), PARAMETER :: separator = REPEAT('-',78)

  CONTAINS

!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
!

  !>
  !!   Conversion: Uppercase -> Lowercase.
  !!
  !!
  FUNCTION tolower (upper)


    CHARACTER(LEN=*)              ,INTENT(in) :: upper
    CHARACTER(LEN=LEN_TRIM(upper))            :: tolower

    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('a')-ICHAR('A')

!-------------------------------------------------------------------------

    DO i=1,LEN_TRIM(upper)
      IF (ICHAR(upper(i:i)) >= ICHAR('A') .AND. &
          ICHAR(upper(i:i)) <= ICHAR('Z')) THEN
        tolower(i:i) = CHAR( ICHAR(upper(i:i)) + idel )
      ELSE
        tolower(i:i) = upper(i:i)
      END IF
    END DO

  END FUNCTION tolower
!-------------------------------------------------------------------------
!

  !>
  !!   Conversion: Lowercase -> Uppercase.
  !!
  !!
  FUNCTION toupper (lower)


    CHARACTER(LEN=*)              ,INTENT(in) :: lower
    CHARACTER(LEN=LEN_TRIM(lower))            :: toupper

    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('A')-ICHAR('a')

!-------------------------------------------------------------------------

    DO i=1,LEN_TRIM(lower)
      IF (ICHAR(lower(i:i)) >= ICHAR('a') .AND. &
          ICHAR(lower(i:i)) <= ICHAR('z')) THEN
        toupper(i:i) = CHAR( ICHAR(lower(i:i)) + idel )
      ELSE
        toupper(i:i) = lower(i:i)
      END IF
    END DO

  END FUNCTION toupper

!-------------------------------------------------------------------------
!

  !>
  !!    Conversion: INTEGER -> CHARACTER(LEN=2).
  !!
  !!
  FUNCTION char2 (i, zero)


    CHARACTER(LEN=2)                       :: char2 ! result
    INTEGER          ,INTENT(in)           :: i     ! argument
    CHARACTER        ,INTENT(in) ,OPTIONAL :: zero  ! padding instead of '0'

    INTEGER ,PARAMETER :: i0 = ICHAR ('0')
!-------------------------------------------------------------------------

    IF (i>99 .OR. i<0) THEN
      char2 = '**'
    ELSE
      char2(1:1) = CHAR(    i/10  + i0)
      char2(2:2) = CHAR(MOD(i,10) + i0)
    ENDIF

    IF(PRESENT(zero)) THEN
      IF(char2(1:1) == '0') char2(1:1) = zero
      IF(char2(2:2) == '0') char2(2:2) = zero
    ENDIF
  END FUNCTION char2

END MODULE mo_util_string
