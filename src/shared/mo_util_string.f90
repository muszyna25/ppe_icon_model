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
  PUBLIC :: split_string         ! splits string into words
  PUBLIC :: string_contains_word ! searches in a string list

  !
  PUBLIC :: normal, bold
  PUBLIC :: fg_black, fg_red, fg_green, fg_blue, fg_magenta, fg_cyan, fg_white, fg_default
  PUBLIC :: bg_black, bg_red, bg_green, bg_blue, bg_magenta, bg_cyan, bg_white, bg_default
  !
  INTERFACE real2string
    MODULE PROCEDURE float2string
    MODULE PROCEDURE double2string
  END INTERFACE real2string
  !
  ! ANSI color sequences 
  !
  CHARACTER(len=1), PARAMETER :: esc = ACHAR(27)
  CHARACTER(len=1), PARAMETER :: orb = ACHAR(91)
  !  
  CHARACTER(len=5) :: normal = esc//orb//'22m'
  CHARACTER(len=4) :: bold   = esc//orb//'1m'
  !
  CHARACTER(len=5) :: fg_black   = esc//orb//'30m'
  CHARACTER(len=5) :: fg_red     = esc//orb//'31m'
  CHARACTER(len=5) :: fg_green   = esc//orb//'32m'
  CHARACTER(len=5) :: fg_yellow  = esc//orb//'33m'
  CHARACTER(len=5) :: fg_blue    = esc//orb//'34m'
  CHARACTER(len=5) :: fg_magenta = esc//orb//'35m'
  CHARACTER(len=5) :: fg_cyan    = esc//orb//'36m'
  CHARACTER(len=5) :: fg_white   = esc//orb//'37m'
  CHARACTER(len=5) :: fg_default = esc//orb//'39m'
  !
  CHARACTER(len=5) :: bg_black   = esc//orb//'40m'
  CHARACTER(len=5) :: bg_red     = esc//orb//'41m'
  CHARACTER(len=5) :: bg_green   = esc//orb//'42m'
  CHARACTER(len=5) :: bg_yellow  = esc//orb//'43m'
  CHARACTER(len=5) :: bg_blue    = esc//orb//'44m'
  CHARACTER(len=5) :: bg_magenta = esc//orb//'45m'
  CHARACTER(len=5) :: bg_cyan    = esc//orb//'46m'
  CHARACTER(len=5) :: bg_white   = esc//orb//'47m'
  CHARACTER(len=5) :: bg_default = esc//orb//'49m'
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
  ! Converts multiple spaces and tabs to single spaces and removes leading spaces.
  !
  SUBROUTINE tocompact(string)
    CHARACTER(len=*), INTENT(inout) :: string
    !
    CHARACTER(len=LEN_TRIM(string)) :: tmp_string
    CHARACTER(len=1):: char
    !
    INTEGER :: i, k, spaces
    !
    string = ADJUSTL(string)
    tmp_string = ' '
    spaces = 0
    k = 0
    !
    DO i = 1, LEN_TRIM(string)
      char = string(i:i)
      SELECT CASE(IACHAR(char))
      CASE (9,32)     ! SPACE and TAB
        IF (spaces == 0) THEN
          k = k+1
          tmp_string(k:k) = ' '
        ENDIF
        spaces = 1
      CASE (33:)      ! everything else
        k = k+1
        tmp_string(k:k) = char
        spaces = 0
      END SELECT
    END DO
    !
    string = ADJUSTL(tmp_string)
    !
  END SUBROUTINE tocompact
  !
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


  !> parses a character string, splits string into words.
  !  This routine takes a comma-separated string like
  !  str = "iconR2B02-grid_DOM01-grid.nc , iconR2B02-grid.nc"
  !  as input and splits it into the components, returning the
  !  number of parts, the start indices and the respective
  !  lengths.
  !  Whitespace is ignored.
  
  SUBROUTINE split_string(zline, n, pos, ilength)

    CHARACTER, PARAMETER :: delim = ',' ! delimiter

    CHARACTER(len=*), INTENT(IN)      :: zline              ! string containing list
    INTEGER,          INTENT(OUT)     :: n                  ! number of parts
    INTEGER,          INTENT(INOUT)   :: pos(:), ilength(:) ! position, lengths of parts
    ! local variables
    INTEGER       :: i           ! index position
    LOGICAL       :: l_word_open ! flag. if true, index "i" is part of a word
    INTEGER       :: istart

    l_word_open = .FALSE.
    n           = 0
    istart      = 1
    DO i=1,LEN(zline)
      IF (.NOT. ((IACHAR(zline(i:i)) ==  9) .OR.  &
        &        (IACHAR(zline(i:i)) == 32) .OR.  &
        &        (zline(i:i) == "'")        .OR.  &
        &        (zline(i:i) == '"')        .OR.  &
        &        (zline(i:i) == delim) )) THEN
        l_word_open = .TRUE.
      ELSE
        IF (l_word_open) THEN
          n = n + 1
          pos(n)  = istart
          ilength(n) = LEN(TRIM(zline(istart:(i-1))))
        END IF
        istart = i+1
        l_word_open = .FALSE.
      END IF
    END DO

  END SUBROUTINE split_string


  !> searches in a string list that has been previously parsed by
  !> "split_string"

  FUNCTION string_contains_word(zword, zline, n, pos, ilength) RESULT(lflag)

    LOGICAL                       :: lflag              ! result
    CHARACTER(len=*), INTENT(IN)  :: zword              ! search word
    CHARACTER(len=*), INTENT(IN)  :: zline              ! string containing list
    INTEGER,          INTENT(IN)  :: n                  ! number of parts
    INTEGER,          INTENT(IN)  :: pos(:), ilength(:) ! position, lengths of parts
    ! local variables
    INTEGER :: i, iwordlen

    lflag     = .FALSE.
    iwordlen  = LEN_TRIM(ADJUSTL(zword))

    DO i=1,n
      IF (ilength(i) /= iwordlen) CYCLE
      lflag = lflag .OR.   &
        &     (zline(pos(i):(pos(i)+ilength(i)-1)) == TRIM(ADJUSTL(zword)))
      IF (lflag) EXIT
    END DO

  END FUNCTION string_contains_word
  !
END MODULE mo_util_string
