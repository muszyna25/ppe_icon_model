!>
!!    This module holds string conversion utilities.
!!
!!
!! @par Revision History
!!    Original version from the ECHAM5 model (version 5.3.01).
!!    Modified to the ProTeX-style by Hui Wan, MPI-M (2007-01-17)
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

! String conversion utilities
MODULE mo_util_string

  USE ISO_C_BINDING,     ONLY: C_INT8_T, C_CHAR
  ! Note: This file must not use mo_exception:finish() to avoid a circular dependency.
  USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH
  USE mo_util_sort,      ONLY: quicksort
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: tolower              ! Conversion   : 'ABCXYZ' -> 'abcxyz'   
  PUBLIC :: toupper              ! Conversion   : 'abcxyz' -> 'ABCXYZ'
  PUBLIC :: separator            ! Format string: (/"-----...-----"/)
  PUBLIC :: int2string           ! returns integer n as a string
  PUBLIC :: real2string          ! returns real n as a string
  PUBLIC :: logical2string       ! returns logical n as a string
  PUBLIC :: split_string         ! splits string into words
  PUBLIC :: string_contains_word ! searches in a string list
  PUBLIC :: tocompact            ! remove gaps in string
  PUBLIC :: str_replace          ! replace any occurrence of keyword by substring
  PUBLIC :: t_keyword_list
  PUBLIC :: associate_keyword    ! add a pair (keyword -> substitution) to a keyword list
  PUBLIC :: with_keywords        ! subroutine for keyword substitution
  PUBLIC :: remove_duplicates
  PUBLIC :: difference
  PUBLIC :: add_to_list
  PUBLIC :: one_of
  PUBLIC :: insert_group
  PUBLIC :: delete_keyword_list
  PUBLIC :: sort_and_compress_list
  PUBLIC :: tohex                ! For debugging: Produce a hex dump of the given string, revealing any unprintable characters.
  PUBLIC :: remove_whitespace
  PUBLIC :: pretty_print_string_list
  PUBLIC :: find_trailing_number

  !functions to handle character arrays as strings
  PUBLIC :: toCharArray     ! convert a fortran string to a character array of kind = c_char
  PUBLIC :: toCharacter     ! convert a character array of kind = c_char back to a fortran string
  PUBLIC :: charArray_dup   ! make a copy of a character array
  PUBLIC :: charArray_equal ! compare two character arrays for equality
  PUBLIC :: charArray_toLower   ! canonicalize to lower case

  !
  PUBLIC :: normal, bold
  PUBLIC :: fg_black, fg_red, fg_green, fg_blue, fg_magenta, fg_cyan, fg_white, fg_default
  PUBLIC :: bg_black, bg_red, bg_green, bg_blue, bg_magenta, bg_cyan, bg_white, bg_default
  !
  INTERFACE real2string
    MODULE PROCEDURE float2string
    MODULE PROCEDURE double2string
  END INTERFACE real2string

  INTERFACE charArray_equal
    MODULE PROCEDURE charArray_equal_array
    MODULE PROCEDURE charArray_equal_char
  END INTERFACE charArray_equal

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

  ! Linked list used for keyword substitution in strings
  TYPE t_keyword_list
    CHARACTER(len=MAX_CHAR_LENGTH) :: keyword    !< keyword string ...
    CHARACTER(len=MAX_CHAR_LENGTH) :: subst      !< ... will be substituted by "subst"
    TYPE(t_keyword_list), POINTER  :: next
  END TYPE t_keyword_list

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_util_string"

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
      ! To eliminate LF and CR generates error reading namelists in restart file by gfortran and NAG!
      ! CASE (9,32,10,13)     ! SPACE and TAB, LF and CR
        CASE (9,32)           ! SPACE and TAB
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


  !------------------------------------------------------------------------------------------------
  !
  ! returns integer n as a string (often needed in printing messages)
  !
  FUNCTION int2string(n, opt_fmt)
    CHARACTER(len=11) :: int2string ! result
    INTEGER, INTENT(in) :: n
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt
    !
    CHARACTER(len=11) :: fmt

    IF (PRESENT(opt_fmt)) THEN
      fmt = opt_fmt
    ELSE
      fmt = '(I11)'
    END IF
    WRITE(int2string,fmt) n
    int2string = ADJUSTL(int2string)
    !
  END FUNCTION int2string
  !------------------------------------------------------------------------------------------------
  !
  ! returns real n as a string (often needed in printing messages)
  !
  FUNCTION float2string(n, opt_fmt) 
    CHARACTER(len=32) :: float2string ! result
    REAL, INTENT(in) :: n
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt
    !
    CHARACTER(len=10) :: fmt
    !
    IF (PRESENT(opt_fmt)) THEN
      fmt = opt_fmt
    ELSE
      fmt = '(g32.5)'
    END IF
    WRITE(float2string,fmt) n
    float2string = ADJUSTL(float2string)
    !
  END FUNCTION float2string
  !
  FUNCTION double2string(n, opt_fmt) 
    CHARACTER(len=32) :: double2string ! result
    DOUBLE PRECISION, INTENT(in) :: n
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt
    !
    CHARACTER(len=10) :: fmt
    !
    IF (PRESENT(opt_fmt)) THEN
      fmt = opt_fmt
    ELSE
      fmt = '(g32.5)'
    END IF
    WRITE(double2string,fmt) n
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


  !> Function for convenience
  !
  !  If "in_str" is matching one of the arguments "arg(i)" return the
  !  index "i". Returns "-1" if none of the strings matches.
  !
  FUNCTION one_of(in_str, arg)
    INTEGER :: one_of
    CHARACTER(len=*), INTENT(IN)           :: in_str    ! input string
    CHARACTER(len=*), INTENT(IN)           :: arg(:)
    ! local variables:
    INTEGER :: i, n, in_str_tlen, arg_tlen
    CHARACTER(len=len_trim(in_str)) :: in_str_upper

#ifndef _CRAYFTN
    one_of = -1
    n = SIZE(arg)
    IF (n > 0) THEN
      in_str_upper = toupper(in_str)
      in_str_tlen = LEN_TRIM(in_str)
      DO i=1,n
        arg_tlen = LEN_TRIM(arg(i))
        IF (arg_tlen == in_str_tlen) THEN
          IF (in_str_upper == toupper(arg(i))) THEN
            one_of=i
            EXIT
          ENDIF
        END IF
      END DO
    END IF
#else ! The crap compiler is to brain-dead to compile and USE the above code.
    one_of = -1
    IF (SIZE(arg) > 0) THEN
      in_str_tlen = LEN_TRIM(in_str)
      DO i=1,SIZE(arg)
        arg_tlen = LEN_TRIM(arg(i))
        IF (arg_tlen == in_str_tlen) THEN
          IF (toupper(in_str(1:in_str_tlen)) == toupper(arg(i)(1:arg_tlen))) THEN
            one_of=i
            EXIT
          ENDIF
        END IF
      END DO
    END IF
#endif
  END FUNCTION one_of


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
    IF ((LEN_TRIM(keyword) > MAX_CHAR_LENGTH) .OR.  &
      & (LEN_TRIM(subst)   > MAX_CHAR_LENGTH))      &
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
    CHARACTER(len=*),   INTENT(OUT) :: keyword
    CHARACTER(len=*),   INTENT(OUT) :: subst
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
    CHARACTER(len=*), INTENT(IN)           :: keyword, subst
    CHARACTER(len=*), INTENT(IN)           :: in_str
    CHARACTER(len=MAX_CHAR_LENGTH)          :: out_str
    ! Local parameters
    INTEGER :: kw_len, in_len, subs_len, pos, out_pos, upper

    out_str = ""
    kw_len   = LEN_TRIM(keyword)
    subs_len = LEN_TRIM(subst)
    in_len   = LEN_TRIM(in_str)
    pos      = 1
    out_pos  = 1
    DO
      IF (pos > in_len) EXIT
      upper = MIN((pos+kw_len-1), in_len)
      IF (in_str(pos:upper) == keyword) THEN
        pos     = pos + kw_len
        ! note: we don't call "finish" to avoid circular dep
        IF ((out_pos+subs_len+in_len-pos-kw_len) > MAX_CHAR_LENGTH) THEN
          WRITE (0,*) "ERROR: str_replace: string too long"
        END IF
        out_str(out_pos:(out_pos+subs_len-1)) = subst(1:subs_len)
        out_pos = out_pos + subs_len
      ELSE
        IF ((out_pos+in_len-pos) > MAX_CHAR_LENGTH) THEN
          WRITE (0,*) "ERROR: str_replace: string too long"
        END IF
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

    nitems_old = nitems

    iwrite = 0
    ITEM_LOOP: DO iread=1,nitems
      ! check if item already in string list (1:iwrite-1):
      DO i=1,iwrite
        IF (str_list(i) == str_list(iread)) CYCLE item_loop
      END DO
      iwrite = iwrite + 1
      IF (iwrite /= iread) str_list(iwrite) = str_list(iread)
    END DO ITEM_LOOP
    nitems = iwrite

    ! clear the rest of the list
    DO iwrite = iwrite+1, nitems_old
      str_list(iwrite) = ' '
    END DO
  END SUBROUTINE remove_duplicates


  !==============================================================================
  !+ Remove entries from a list of strings which occur in a second list.
  !
  ! @note This is a very crude implementation, quadratic complexity.
  !
  SUBROUTINE difference(str_list1, nitems1, str_list2, nitems2)
    CHARACTER(len=*),          INTENT(INOUT) :: str_list1(:)
    INTEGER,                   INTENT(INOUT) :: nitems1
    CHARACTER(len=*),          INTENT(IN)    :: str_list2(:)
    INTEGER,                   INTENT(IN)    :: nitems2
    ! local variables
    INTEGER :: iwrite, iread, nitems_old, i

    nitems_old = nitems1

    iwrite = 1
    ITEM1_LOOP: DO iread=1,nitems_old
      ! check if item is in string list 2:
      DO i=1,nitems2
        IF (str_list2(i) == str_list1(iread)) CYCLE ITEM1_LOOP
      END DO
      ! can only be reached for non-duplicate entries
      IF (iwrite /= iread) str_list1(iwrite) = str_list1(iread)
      iwrite = iwrite + 1
    END DO ITEM1_LOOP
    nitems1 = iwrite-1

    ! clear the rest of the list
    DO iwrite=iwrite,nitems_old
      str_list1(iwrite) = ' '
    END DO
  END SUBROUTINE difference




  !==============================================================================
  !+ Add entries from list 2 to list 1, if they are not already present
  !+ in list 1.
  !
  ! @note This is a very crude implementation, quadratic complexity.
  !
  SUBROUTINE add_to_list(str_list1, nitems1, str_list2, nitems2)
    CHARACTER(len=*),          INTENT(INOUT) :: str_list1(:)
    INTEGER,                   INTENT(INOUT) :: nitems1
    CHARACTER(len=*),          INTENT(IN)    :: str_list2(:)
    INTEGER,                   INTENT(IN)    :: nitems2
    ! local variables
    INTEGER :: iread, i
    LOGICAL :: l_duplicate


    ! Loop over all items that should potentially be added
    DO iread=1,nitems2
      ! check if item is already in string list 1:
      l_duplicate = .FALSE.
      ! Loop over all items in the target list (list 1)
      CHECK_LOOP : DO i=1,nitems1
        IF (TRIM(str_list1(i)) == TRIM(str_list2(iread))) THEN
          l_duplicate = .TRUE.
          EXIT CHECK_LOOP
        END IF
      END DO CHECK_LOOP
      IF (.NOT. l_duplicate) THEN
        str_list1(nitems1+1) = str_list2(iread)
        nitems1 = nitems1+1
      END IF
    END DO

  END SUBROUTINE add_to_list


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
    CHARACTER(len=MAX_CHAR_LENGTH) :: result_str, subst, keyword

    ! note: we don't call "finish" to avoid circular dep
    IF (LEN_TRIM(in_str) > MAX_CHAR_LENGTH) &
      & WRITE (0,*) "ERROR: with_keywords: string too long"
    result_str = in_str
    IF (.NOT. ASSOCIATED(keyword_list)) RETURN
    DO
      CALL keyword_list_pop(keyword_list, keyword, subst)
      result_str = str_replace(result_str, keyword, subst)
      IF (.NOT. ASSOCIATED(keyword_list)) RETURN
    END DO
  END FUNCTION with_keywords


  !==============================================================================
  !+ Subroutine for keyword substitution
  !
  ! If we have a list of strings, for example "u", "v", "tracers",
  ! then we can use this function to replace a keyword that denotes a
  ! whole group of variables (like "tracers"), for example by
  ! group_list="Q1", "Q2", etc.
  !
  ! @param[in]  varlist    original array of strings (variable names)
  ! @param[in]  vname_len  length of each string
  ! @param[in]  n          length of list
  ! @param[in]  group_name substitution keyword (i.e. variable group name)
  ! @param[in]  group_list array of strings that will be inserted
  !
  ! @return     contents of @p varlist where @p group_name has been replaced.
  !------------------------------------------------------------------------------
  SUBROUTINE insert_group(varlist, vname_len, n, group_name, group_list, result_list)
    INTEGER,                        INTENT(IN)    :: vname_len, n
    CHARACTER(LEN=vname_len),       INTENT(INOUT) :: result_list(n)
    CHARACTER(LEN=*),               INTENT(IN)    :: varlist(:), group_list(:)
    CHARACTER(LEN=*),               INTENT(IN)    :: group_name
    ! local variables
    INTEGER :: i,j,k,m,ngroups
    CHARACTER(len=LEN_TRIM(group_name)) :: group_name_uc

    k=0
    ngroups = SIZE(group_list)
    m = SIZE(varlist)
    IF (m > 0) THEN
      group_name_uc = toupper(group_name)
      DO i=1,m
        IF (varlist(i) == ' ') EXIT
        IF (toupper(varlist(i)) == group_name_uc) THEN
          DO j=1,ngroups
            k = k+1
            result_list(k) = group_list(j)
          END DO
        ELSE
          k = k+1
          result_list(k) = varlist(i)
        END IF
      END DO
      CALL remove_duplicates(result_list, k )
    END IF
    DO i=k+1,n
      result_list(i) = " "
    END DO

  END SUBROUTINE insert_group

  !==============================================================================
  SUBROUTINE delete_keyword_list(list_head)
    ! Parameters
    TYPE(t_keyword_list),    POINTER    :: list_head, next

    DO WHILE (ASSOCIATED(list_head))
      next => list_head%next
      DEALLOCATE(list_head)
      list_head => next
    END DO

  END SUBROUTINE delete_keyword_list
  !==============================================================================


  !> Utility function: Takes a list of integer values as an input
  !> (without duplicates) and returns a string containing this list in
  !> an ordered, compressed fashion.
  !
  !  E.g., the list
  !     ( 1, 10, 9, 8, 3, 5, 6 )
  !  is transformed into
  !      1,3,5,6,8-10
  !
  !  Initial implementation: 2014-02-14   F. Prill (DWD)
  !
  SUBROUTINE sort_and_compress_list(idx_list, dst)
    INTEGER,          INTENT(IN)  :: idx_list(:)
    CHARACTER(LEN=*), INTENT(OUT) :: dst
    ! local variables
    INTEGER :: list(SIZE(idx_list)),  &  ! sorted copy
      &        nnext(SIZE(idx_list))
    INTEGER :: i, j, N

    dst = " "
    N = SIZE(idx_list)
    list(:) = idx_list(:)
    ! sort the list
    CALL quicksort(list)
    ! find out, how many direct successors follow:
    j        = 1
    nnext(:) = 0
    nnext(1) = 1
    DO i=2,N
      IF (list(i) == list(i-1)) THEN
        WRITE (0,*) "ERROR: sort_and_compress_list operates on non-unique list entries"
      END IF
      IF (list(i) == (list(i-1)+1)) THEN
        nnext(j) = nnext(j) + 1
      ELSE
        j = i
        nnext(j) = 1
      END IF
    END DO
    ! build the result string:
    i = 1
    DO WHILE (i <= n)
      IF (nnext(i) > 1) THEN
        IF (nnext(i) == 2) THEN
          dst = TRIM(dst)//TRIM(int2string(list(i)))//","//TRIM(int2string(list(i+1)))
        ELSE
          dst = TRIM(dst)//TRIM(int2string(list(i)))//"-"//TRIM(int2string(list(i+nnext(i)-1)))
        END IF
        i = i + nnext(i)
      ELSE
        dst = TRIM(dst)//TRIM(int2string(list(i)))
        i = i + 1
      END IF
      IF (i <= N)  dst = TRIM(dst)//", "
    END DO
  END SUBROUTINE sort_and_compress_list

  FUNCTION tohex_internal(inData) RESULT(resultVar)
    INTEGER(KIND = C_INT8_T), INTENT(IN) :: inData(:)
    CHARACTER(LEN = 3*SIZE(inData, 1) - 1) :: resultVar

    CHARACTER(LEN = 16), PARAMETER :: nibbles = "0123456789abcdef"
    INTEGER :: inputIndex, outputIndex, curChar, nibble1, nibble2

    outputIndex = 1
    DO inputIndex = 1, SIZE(inData, 1)
        IF(inputIndex /= 1) THEN
            resultVar(outputIndex:outputIndex) = " "
            outputIndex = outputIndex + 1
        END IF
        curChar = inData(inputIndex)
        IF(curChar < 0) curChar = curChar + 256
        nibble1 = ISHFT(curChar, -4) + 1
        nibble2 = IAND(curChar, 15) + 1
        resultVar(outputIndex:outputIndex) = nibbles(nibble1:nibble1)
        resultVar(outputIndex+1:outputIndex+1) = nibbles(nibble2:nibble2)
        outputIndex = outputIndex + 2
    END DO
  END FUNCTION tohex_internal

  FUNCTION tohex(string) RESULT(resultVar)
    CHARACTER(LEN = *), INTENT(IN) :: string
    CHARACTER(LEN = 3*LEN(string) - 1) :: resultVar

    INTEGER(KIND = C_INT8_T) :: mold(1)
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":tohex"

    IF(LEN(resultVar) /= LEN(tohex_internal(TRANSFER(string, mold)))) THEN
        ! throw error if the returned SIZE is wrong
        ! note: we don't call "finish" to avoid circular dep
        WRITE(0,*) "fatal error: "//modname//":tohex_internal() returned string of unexpected length"
        resultVar = "fatal error: "//modname//":tohex_internal() returned string of unexpected length"
    ELSE
        resultVar = tohex_internal(TRANSFER(string, mold))
    END IF
  END FUNCTION tohex

  !------------------------------------------------------------------------------------------------
  !> Remove all white space from a string (also between "words").
  !
  FUNCTION remove_whitespace(in_str)
    CHARACTER(len=*), INTENT(in)    :: in_str
    CHARACTER(len=LEN_TRIM(in_str)) :: remove_whitespace
    ! local variables
    INTEGER   :: i,j, ichar

    remove_whitespace = " "
    j = 0
    DO i = 1,LEN(in_str)
      ichar = IACHAR(in_str(i:i))
      IF ((ichar /= 9) .AND. (ichar /= 32)) THEN
        j = j + 1
        remove_whitespace(j:j) = in_str(i:i)
      END IF
    END DO
  END FUNCTION remove_whitespace

  FUNCTION toCharArray(string) RESULT(resultVar)
    CHARACTER(LEN = *), INTENT(IN) :: string
    CHARACTER(KIND = C_CHAR), POINTER :: resultVar(:)
    INTEGER :: i, error

    CHARACTER(LEN = *), PARAMETER :: routine = modName//":toCharArray"

    ALLOCATE(resultVar(LEN(string)), STAT = error)
    ! note: we don't call "finish" to avoid circular dependency
    IF(error /= 0) WRITE(0,*) "memory allocation error"
    DO i = 1, LEN(string)
        resultVar(i) = string(i:i)
    END DO
  END FUNCTION toCharArray

  FUNCTION toCharacter(charArray) RESULT(resultVar)
    CHARACTER(KIND = C_CHAR), INTENT(IN) :: charArray(:)
    CHARACTER(LEN = :), POINTER :: resultVar
    INTEGER :: i, error, stringSize

    CHARACTER(LEN = *), PARAMETER :: routine = modName//":toCharacter"

    stringSize = SIZE(charArray, 1) !XXX: This may not be merged into the next line, because that triggers a bug in gfortran
    ALLOCATE(CHARACTER(LEN = stringSize) :: resultVar, STAT = error)
    ! note: we don't call "finish" to avoid circular dependency
    IF(error /= 0) WRITE(0,*) "memory allocation error"
    DO i = 1, SIZE(charArray, 1)
        resultVar(i:i) = charArray(i)
    END DO
  END FUNCTION toCharacter

  FUNCTION charArray_dup(charArray) RESULT(resultVar)
    CHARACTER(KIND = C_CHAR), INTENT(IN) :: charArray(:)
    CHARACTER(KIND = C_CHAR), POINTER :: resultVar(:)
    INTEGER :: error

    CHARACTER(LEN = *), PARAMETER :: routine = modName//":charArray_dup"

    ALLOCATE(resultVar(SIZE(charArray, 1)), STAT = error)
    ! note: we don't call "finish" to avoid circular dependency
    IF(error /= 0) WRITE(0,*) "memory allocation error"
    resultVar(:) = charArray(:)
  END FUNCTION charArray_dup

  LOGICAL FUNCTION charArray_equal_array(stringA, stringB) RESULT(resultVar)
    CHARACTER(KIND = C_CHAR), INTENT(IN) :: stringA(:), stringB(:)
    INTEGER :: i

    resultVar = .FALSE.
    IF(SIZE(stringA, 1) /= SIZE(stringB, 1)) RETURN
    DO i = 1, SIZE(stringA, 1)
        IF(stringA(i) /= stringB(i)) RETURN
    END DO
    resultVar = .TRUE.
  END FUNCTION charArray_equal_array

  LOGICAL FUNCTION charArray_equal_char(stringA, stringB) RESULT(resultVar)
    CHARACTER(KIND = C_CHAR), INTENT(IN) :: stringA(:)
    CHARACTER(*), INTENT(IN) :: stringB
    INTEGER :: i

    resultVar = .FALSE.
    IF(SIZE(stringA, 1) /= LEN_TRIM(stringB)) RETURN
    DO i = 1, SIZE(stringA, 1)
        IF(stringA(i) /= stringB(i:i)) RETURN
    END DO
    resultVar = .TRUE.
  END FUNCTION charArray_equal_char

  SUBROUTINE charArray_toLower(string)
    CHARACTER(KIND = C_CHAR), INTENT(INOUT) :: string(:)
    INTEGER :: i, curChar

    DO i = 1, SIZE(string, 1)
        curChar = IACHAR(string(i))
        IF(curChar >= IACHAR('A') .AND. curChar <= IACHAR('Z')) THEN
            curChar = curChar - IACHAR('A') + IACHAR('a')
            string(i) = ACHAR(curChar)
        END IF
    END DO
  END SUBROUTINE charArray_toLower

  !> "pretty-print" a list of strings (comma-separated).
  !
  !  After at most "max_ll" characters-per-line a new line is inserted.
  !  Each line is indented by an (optional) prefix string.
  !
  SUBROUTINE pretty_print_string_list(list, opt_max_ll, opt_dst, opt_prefix)
    CHARACTER(LEN=*),           INTENT(IN) :: list(:)     ! string list for print-out
    INTEGER, OPTIONAL,          INTENT(IN) :: opt_max_ll  ! max. line length
    INTEGER, OPTIONAL,          INTENT(IN) :: opt_dst     ! (optional:) WRITE destination
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: opt_prefix  ! (optional:) line prefix

    INTEGER :: dst, max_ll, ccnt, i, len_i
    CHARACTER(len=:), ALLOCATABLE :: prefix

    dst = 0
    IF (PRESENT(opt_dst))     dst = opt_dst
    max_ll = 80
    IF (PRESENT(opt_max_ll))  max_ll = opt_max_ll

    IF (PRESENT(opt_prefix)) THEN
      prefix = opt_prefix
    ELSE
      ALLOCATE(CHARACTER(1) :: prefix)
      prefix = " "
    END IF

    ccnt = max_ll
    DO i=1,SIZE(list)
      len_i = LEN_TRIM(list(i))
      ! start new line (if necessary)
      IF (ccnt + len_i + 2 > max_ll) THEN
        IF (i>1)  WRITE (dst,"(a)") " "
        WRITE (dst,"(a)", advance='no') prefix
        ccnt = LEN(prefix)
      END IF
      WRITE (dst,"(a)", advance='no') TRIM(list(i))
      IF (i < SIZE(list))  WRITE (dst,"(a)", advance='no') ", "
      ccnt = ccnt + len_i + 2
    END DO
    WRITE (dst,"(a)") " "

    DEALLOCATE(prefix)
  END SUBROUTINE pretty_print_string_list


  !> find position of numeric suffix in the character string "str",
  !  return "-1" if no such suffix is found.
  !
  RECURSIVE FUNCTION find_trailing_number(str) RESULT(pos)
    INTEGER                      :: pos
    CHARACTER(LEN=*), INTENT(IN) :: str  !< input string
    INTEGER :: len

    pos   = -1
    len   = LEN_TRIM(str)
    IF (len == 0) RETURN

    IF (is_number(str(len:len))) THEN
      IF (len > 1)   pos = find_trailing_number(str(1:(len-1)))
      IF (pos == -1) pos = len
    END IF

  CONTAINS
    LOGICAL FUNCTION is_number(char)
      CHARACTER, INTENT(IN) :: char
      is_number = (IACHAR(char) - IACHAR('0')) <= 9
    END FUNCTION is_number
  END FUNCTION find_trailing_number

END MODULE mo_util_string
