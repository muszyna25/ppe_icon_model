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
  !
  ! String conversion utilities
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: tolower        ! Conversion   : 'ABCXYZ' -> 'abcxyz'   
  PUBLIC :: toupper        ! Conversion   : 'abcxyz' -> 'ABCXYZ'
  PUBLIC :: char2          ! Conversion   : INTEGER  -> CHAR (LEN=2)
  PUBLIC :: separator      ! Format string: (/"-----...-----"/)
  PUBLIC :: int2string     ! returns integer n as a string
  PUBLIC :: real2string    ! returns real n as a string
  PUBLIC :: logical2string ! returns logical n as a string
  !
  INTERFACE real2string
    MODULE PROCEDURE float2string
    MODULE PROCEDURE double2string
  END INTERFACE real2string
  !
  CHARACTER(len=*), PARAMETER :: separator = REPEAT('-',100)
  !
CONTAINS
  !
  !------------------------------------------------------------------------------------------------
  !
  ! Conversion: Uppercase -> Lowercase
  !
  FUNCTION tolower (uppercase)
    CHARACTER(len=*), INTENT(in) :: uppercase
    CHARACTER(len=LEN_TRIM(uppercase)) :: tolower
    !
    INTEGER, PARAMETER :: idel = ICHAR('a')-ICHAR('A')
    INTEGER :: i
    !
    DO i = 1, LEN_TRIM(uppercase)
      IF (ICHAR(uppercase(i:i)) >= ICHAR('A') .AND. ICHAR(uppercase(i:i)) <= ICHAR('Z')) THEN
        tolower(i:i) = CHAR( ICHAR(uppercase(i:i)) + idel )
      ELSE
        tolower(i:i) = uppercase(i:i)
      ENDIF
    ENDDO
    !
  END FUNCTION tolower
  !------------------------------------------------------------------------------------------------
  !
  ! Conversion: Lowercase -> Uppercase
  !
  FUNCTION toupper (lowercase)
    CHARACTER(len=*), INTENT(in) :: lowercase
    CHARACTER(len=LEN_TRIM(lowercase)) :: toupper
    !
    INTEGER, PARAMETER :: idel = ICHAR('A')-ICHAR('a')
    INTEGER :: i
    !
    DO i = 1, LEN_TRIM(lowercase)
      IF (ICHAR(lowercase(i:i)) >= ICHAR('a') .AND. ICHAR(lowercase(i:i)) <= ICHAR('z')) THEN
        toupper(i:i) = CHAR( ICHAR(lowercase(i:i)) + idel )
      ELSE
        toupper(i:i) = lowercase(i:i)
      ENDIF
    ENDDO
    !
  END FUNCTION toupper
  !------------------------------------------------------------------------------------------------
  !
  ! Conversion: INTEGER -> CHARACTER(LEN=2)
  !
  FUNCTION char2 (i, zero)
    CHARACTER(len=2) :: char2 ! result
    INTEGER,   INTENT(in)           :: i     ! argument
    CHARACTER, INTENT(in), OPTIONAL :: zero  ! padding instead of '0'
    !
    INTEGER, PARAMETER :: i0 = ICHAR ('0')
    !
    IF (i > 99 .OR. i < 0) THEN
      char2 = '**'
    ELSE
      char2(1:1) = CHAR(    i/10  + i0)
      char2(2:2) = CHAR(MOD(i,10) + i0)
    ENDIF
    !
    IF (PRESENT(zero)) THEN
      IF(char2(1:1) == '0') char2(1:1) = zero
      IF(char2(2:2) == '0') char2(2:2) = zero
    ENDIF
    !
  END FUNCTION char2
  !------------------------------------------------------------------------------------------------
  !
  ! returns integer n as a string (often needed in printing messages)
  !
  FUNCTION int2string(n)
    CHARACTER(len=10) :: int2string ! result
    INTEGER, INTENT(in) :: n
    !
    WRITE(int2string,'(I10)') n
    int2string = ADJUSTL(int2string)
    !
  END FUNCTION int2string
  !------------------------------------------------------------------------------------------------
  !
  ! returns real n as a string (often needed in printing messages)
  !
  FUNCTION float2string(n) 
    CHARACTER(len=32) :: float2string ! result
    REAL, INTENT(in) :: n
    !
    WRITE(float2string,'(g32.5)') n
    float2string = ADJUSTL(float2string)
    !
  END FUNCTION float2string
  !
  FUNCTION double2string(n) 
    CHARACTER(len=32) :: double2string ! result
    DOUBLE PRECISION, INTENT(in) :: n
    !
    WRITE(double2string,'(g32.5)') n
    double2string = ADJUSTL(double2string)
    !
  END FUNCTION double2string
  !------------------------------------------------------------------------------------------------
  !
  ! returns integer n as a string (often needed in printing messages)
  !
  FUNCTION logical2string(n)
    CHARACTER(len=10) :: logical2string ! result
    LOGICAL, INTENT(in) :: n
    !
    WRITE(logical2string,'(l10)') n
    logical2string = ADJUSTL(logical2string)
    !
  END FUNCTION logical2string
  !
END MODULE mo_util_string
