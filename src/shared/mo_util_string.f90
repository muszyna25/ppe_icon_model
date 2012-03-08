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
  PUBLIC :: tocompact      ! remove gaps in string
  PUBLIC :: str_replace       ! replace any occurrence of keyword by substring
  PUBLIC :: t_keyword_list
  PUBLIC :: associate_keyword ! add a pair (keyword -> substitution) to a keyword list
  PUBLIC :: with_keywords     ! subroutine for keyword substitution
  PUBLIC :: MAX_STRING_LEN
  PUBLIC :: remove_duplicates

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
  ! String length constants used for keyword substitution
  INTEGER, PARAMETER :: MAX_STRING_LEN  = 256
  INTEGER, PARAMETER :: MAX_KEYWORD_LEN = 128

  ! Linked list used for keyword substitution in strings
  TYPE t_keyword_list
    CHARACTER(len=MAX_KEYWORD_LEN) :: keyword    !< keyword string ...
    CHARACTER(len=MAX_KEYWORD_LEN) :: subst      !< ... will be substituted by "subst"
    TYPE(t_keyword_list), POINTER  :: next
  END TYPE t_keyword_list

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
  ! Converts multiple spaces and tabs to single spaces and removes
  ! leading spaces.
  !
  SUBROUTINE tocompact(string)
    CHARACTER(len=*), INTENT(inout) :: string
    ! local variables
    INTEGER   :: offset, i, i_max
    CHARACTER :: char
    LOGICAL   :: lspaces

    offset = 0
    i      = 0
    i_max  = LEN_TRIM(string)
    LOOP : DO 
      i = i + 1       ! current write pos
      IF ((i+offset) > i_max) EXIT LOOP
      lspaces = .FALSE.
      LOOKAHEAD : DO
        char = string((i+offset):(i+offset))
        SELECT CASE(IACHAR(char))
        CASE (9,32)     ! SPACE and TAB
          offset  = offset + 1
          IF ((i+offset) > i_max) EXIT LOOP
          lspaces = (i>1)
        CASE default
          IF (lspaces) THEN
            string(i:i) = ' '
            i      = i      + 1
            offset = offset - 1
          END IF
          string(i:i) = char
          EXIT LOOKAHEAD
        END SELECT
      END DO LOOKAHEAD
    END DO LOOP
    string = string(1:(i-1))
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
  FUNCTION int2string(n, opt_fmt)
    CHARACTER(len=10) :: int2string ! result
    INTEGER, INTENT(in) :: n
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt
    !
    CHARACTER(len=10) :: fmt

    IF (PRESENT(opt_fmt)) THEN
      fmt = opt_fmt
    ELSE
      fmt = '(I10)'
    END IF
    WRITE(int2string,fmt) n
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

  !==============================================================================
  !+ Utility function: Insert (keyword, substitution) pair into keyword list
  !------------------------------------------------------------------------------
  SUBROUTINE keyword_list_push(keyword, subst, list_head)
    ! Parameters
    CHARACTER(len=*),        INTENT(IN) :: keyword, subst
    TYPE(t_keyword_list),    POINTER    :: list_head
    ! Local parameters
    TYPE (t_keyword_list),   POINTER    :: tmp
    INTEGER                             :: errstat

    ! throw error if keyword, subst are too long
    ! note: we don't call "finish" to avoid circular dep
    IF ((LEN_TRIM(keyword) > MAX_KEYWORD_LEN) .OR.  &
      & (LEN_TRIM(subst)   > MAX_KEYWORD_LEN))      &
      &  WRITE (0,*) "ERROR: keyword_list_push: keyword too long"

    ! insert element into linked list
    tmp => list_head
    ALLOCATE(list_head, stat=errstat)
    IF (errstat /= 0) &
      & WRITE (0,*) "ERROR: keyword_list_push: ALLOCATE"
    list_head%keyword = keyword
    list_head%subst   = subst
    list_head%next    => tmp

  END SUBROUTINE keyword_list_push
  

  !==============================================================================
  !+ Utility function: Get (keyword, substitution) pair from keyword list
  !------------------------------------------------------------------------------
  SUBROUTINE keyword_list_pop(list_head, keyword, subst)
    ! Parameters
    CHARACTER(len=MAX_KEYWORD_LEN),  INTENT(OUT) :: keyword, subst
    TYPE(t_keyword_list),            POINTER     :: list_head
    ! Local parameters
    TYPE (t_keyword_list),   POINTER    :: tmp
    INTEGER                             :: errstat
    
    IF (.NOT. ASSOCIATED(list_head)) THEN
      keyword = ""
      subst   = ""
    ELSE
      ! remove list head
      keyword =  list_head%keyword
      subst   =  list_head%subst
      tmp     => list_head%next
      DEALLOCATE(list_head, STAT=errstat)
      ! note: we don't call "finish" to avoid circular dep
      IF (errstat /= 0) &
        & WRITE (0,*) "ERROR: keyword_list_pop: DEALLOCATE"
      list_head => tmp
    END IF
    
  END SUBROUTINE keyword_list_pop


  !==============================================================================
  !+ Utility function: replace any occurrence of keyword by substring
  !------------------------------------------------------------------------------
  FUNCTION str_replace(in_str, keyword, subst) RESULT(out_str)
    ! Parameters
    CHARACTER(len=MAX_KEYWORD_LEN), INTENT(IN) :: keyword, subst
    CHARACTER(len=*),               INTENT(IN) :: in_str
    CHARACTER(len=MAX_STRING_LEN)              :: out_str
    ! Local parameters
    INTEGER :: kw_len, in_len, subs_len, pos, out_pos
    
    out_str = ""
    kw_len   = LEN_TRIM(keyword)
    subs_len = LEN_TRIM(subst)
    in_len   = LEN_TRIM(in_str)
    pos      = 1
    out_pos  = 1
    DO
      IF (pos > in_len) EXIT
      IF (in_str(pos:(pos+kw_len-1)) == keyword) THEN
        pos     = pos + kw_len
        ! note: we don't call "finish" to avoid circular dep
        IF ((out_pos + subs_len) > MAX_STRING_LEN) &
          & WRITE (0,*) "ERROR: str_replace: string too long"
        out_str(out_pos:(out_pos+subs_len-1)) = subst(1:subs_len)
        out_pos = out_pos + subs_len
      ELSE
        IF (out_pos > MAX_STRING_LEN) &
          & WRITE (0,*) "ERROR: str_replace: string too long"
        out_str(out_pos:out_pos) = in_str(pos:pos)
        pos     = pos + 1
        out_pos = out_pos + 1
      END IF
    END DO
    
  END FUNCTION str_replace


  !==============================================================================
  !+ Remove duplicate entries from a list of strings.
  !
  ! @note This is a very crude implementation, quadratic complexity.
  !
  SUBROUTINE remove_duplicates(str_list, nitems)
    CHARACTER(len=*),          INTENT(INOUT) :: str_list(:)
    INTEGER,                   INTENT(INOUT) :: nitems
    ! local variables
    INTEGER :: iwrite, iread, nitems_old, i
    LOGICAL :: l_duplicate
    
    nitems_old = nitems
    
    iwrite = 1
    DO iread=1,nitems
      ! check if item already in string list (1:iwrite-1):
      l_duplicate = .FALSE.
      CHECK_LOOP : DO i=1,(iwrite-1)
        IF (TRIM(str_list(i)) == TRIM(str_list(iread))) THEN
          l_duplicate = .TRUE.
          EXIT CHECK_LOOP
        END IF
      END DO CHECK_LOOP
      IF (.NOT. l_duplicate) THEN
        str_list(iwrite) = str_list(iread)
        iwrite = iwrite + 1
      END IF
    END DO
    nitems = iwrite-1
    
    ! clear the rest of the list
    DO iwrite=(nitems+1),nitems_old
      str_list(iwrite) = ' '
    END DO
  END SUBROUTINE remove_duplicates


  !==============================================================================
  !+ Add a pair (keyword -> substitution) to a keyword list
  ! see FUNCTION with_keywords for further documentation.
  !------------------------------------------------------------------------------
  SUBROUTINE associate_keyword(keyword, subst, keyword_list)
    CHARACTER(len=*), INTENT(IN)   :: keyword, subst
    TYPE(t_keyword_list), POINTER  :: keyword_list
    
    CALL keyword_list_push(keyword, subst, keyword_list)
  END SUBROUTINE associate_keyword


  !==============================================================================
  !+ Subroutine for keyword substitution
  ! Usage example: Consider the following code snippet:
  ! \code
  !    CHARACTER(len=*), PARAMETER    :: filename = "<path>/bin/<prefix>grid.nc"
  !    TYPE (t_keyword_list), POINTER :: keywords => NULL()
  !    CALL associate_keyword("<path>",   "/usr/local", keywords)
  !    CALL associate_keyword("<prefix>", "exp01_",     keywords)
  ! \endcode
  ! Then, by calling 'with_keywords(keywords, filename)',
  ! the filename is transformed into '/usr/local/bin/exp01_grid.nc'.
  ! 
  !------------------------------------------------------------------------------
  FUNCTION with_keywords(keyword_list, in_str) RESULT(result_str)
    TYPE(t_keyword_list), POINTER  :: keyword_list
    CHARACTER(len=*), INTENT(IN)   :: in_str
    CHARACTER(len=MAX_STRING_LEN)  :: result_str
    CHARACTER(len=MAX_KEYWORD_LEN) :: keyword, subst
    
    ! note: we don't call "finish" to avoid circular dep
    IF (LEN_TRIM(in_str) > MAX_STRING_LEN) &
      & WRITE (0,*) "ERROR: with_keywords: string too long"
    result_str = in_str
    IF (.NOT. ASSOCIATED(keyword_list)) RETURN
    DO
      CALL keyword_list_pop(keyword_list, keyword, subst)
      result_str = str_replace(result_str, keyword, subst)
      IF (.NOT. ASSOCIATED(keyword_list)) RETURN
    END DO
  END FUNCTION with_keywords
  
END MODULE mo_util_string
