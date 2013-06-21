!>
!! Utility module: print namelist to file, annotating all changed values.
!!
!! The implementation is more or less a "poor man's approach". We write
!! and re-read the namelists in text form, since Fortran is missing the
!! language features for looping over namelist elements.
!!
!! @par Revision History
!!  by F. Prill, DWD (2013-06-13)
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_nml_annotate

#ifdef __ICON__
  USE mo_namelist,    ONLY: position_nml, POSITIONED
  USE mo_util_string, ONLY: int2string, str_replace, tolower,    &
    &                       tocompact
  USE mo_io_units,    ONLY: find_next_free_unit, nnml
  USE mo_exception,   ONLY: finish
#else
  USE mo_utilities, ONLY: nnml, POSITIONED,                      &
    &                     int2string, str_replace, tolower,      &
    &                     finish, tocompact, find_next_free_unit
#endif
  USE mo_mpi,       ONLY: get_my_mpi_work_id

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: MAX_STRBUF_LENGTH  = 32767 !< max. size of a single namelist
  INTEGER, PARAMETER :: MAX_STRBUF_LINES   =  1024 !< max. no. of lines in a namelist

  INTEGER            :: tmpnml1            =  3    !< file handle for temporary text file with defaults
  INTEGER            :: tmpnml2            =  4    !< file handle for temporary text file with settings

  CHARACTER (LEN=*), PARAMETER :: TMP_FILENAME1 = "tmpfile1.mo_nml_annotate" !< temporary file with defaults
  CHARACTER (LEN=*), PARAMETER :: TMP_FILENAME2 = "tmpfile2.mo_nml_annotate" !< temporary file with settings

  !> internal buffer variables
  CHARACTER (len=MAX_STRBUF_LENGTH)  :: nml_log_buffer1, nml_log_buffer2

  !> start positions of the several namelist lines in the buffer
  INTEGER :: startpos1(MAX_STRBUF_LINES), startpos2(MAX_STRBUF_LINES)

  LOGICAL :: linitialized = .FALSE.

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_nml_annotate')

  PUBLIC :: temp_defaults
  PUBLIC :: temp_settings
  PUBLIC :: log_nml_settings 

CONTAINS


  !> Decode name of the namelist from buffer
  !
  SUBROUTINE parse_nml_name(buffer, ipos, name)
    CHARACTER (len=*),               INTENT(IN)    :: buffer  !< character buffer containing namelist
    INTEGER,                         INTENT(INOUT) :: ipos    !< current/new position in buffer
    CHARACTER (len=*), INTENT(OUT)   :: name    !< result: name of the namelist
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::parse_nml_name'
    INTEGER :: ipos_end

    ! starting token is an '&' character
    ipos = INDEX(buffer(ipos:), '&')
    IF (ipos == 0) THEN
      WRITE (0,*) "buffer = ", buffer(ipos:)
      CALL finish(routine,"Invalid start for namelist block!")
    END IF
    ! determine position of next ' ' character
    ipos_end = ipos + INDEX(buffer(ipos:), ' ')
    ! return substring
    name = tolower(buffer((ipos+1):(ipos_end-1)))
    CALL tocompact(name)
    ipos = ipos_end-1
  END SUBROUTINE parse_nml_name


  !> Decode the next pair of "key = value" from the buffer.
  !
  !  Note: Leading, trailing and internal blanks are removed from the value
  !        string.
  !
  FUNCTION parse_key_value_pair(buffer, ipos, iline, startpos, key, value)
    LOGICAL :: parse_key_value_pair
    CHARACTER (len=*),                 INTENT(IN)    :: buffer       !< character buffer containing namelist
    INTEGER,                           INTENT(INOUT) :: ipos, iline  !< current/new position in buffer
    INTEGER,                           INTENT(IN)    :: startpos(:)  !< array of line start positions
    CHARACTER (len=*),                 INTENT(OUT)   :: key          !< result: key (string)
    CHARACTER (len=MAX_STRBUF_LENGTH), INTENT(OUT)   :: value        !< result: value (string)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::parse_key_value_pair'
    INTEGER :: ipos_end, ipos_end2, lt

    parse_key_value_pair = .TRUE.
    IF (ipos == 0) THEN
      parse_key_value_pair = .FALSE.
      RETURN
    END IF

    ! determine position of next '=' character
    ipos_end = INDEX(buffer(ipos:), '=')
    IF (ipos_end == 0) THEN
      ! not found
      parse_key_value_pair = .FALSE.
      RETURN
    END IF
    ipos_end = ipos_end + ipos

    key = tolower(buffer(ipos:(ipos_end-2)))
    key = ADJUSTL(key)
    ! determine position of next ',' character
    ipos_end2 = startpos(iline+2)
    IF (ipos_end2 == 0) ipos_end2 = LEN_TRIM(buffer)+2
    value = buffer((ipos_end):(ipos_end2-2))
    ! strip "value" from trailing commas and slashes:
    lt = LEN_TRIM(value)
    IF (value(lt:lt) == "/") value(lt:lt) = " "
    lt = LEN_TRIM(value)
    IF (value(lt:lt) == ",") value(lt:lt) = " "
    ! if "value" is a string (delimiter: "'") then trim leading, trailing
    ! blanks:
    CALL tocompact(value)
    value = str_replace(value, "' ", "'")
    value = str_replace(value, " '", "'")
    value = ADJUSTL(value)

    ipos  = ipos_end2
    iline = iline + 1
  END FUNCTION parse_key_value_pair


  SUBROUTINE parse_namelist(buffer1, buffer2, col_width, gap_width, opt_file)
    CHARACTER(LEN=*), INTENT(IN) :: buffer1, buffer2       !< character buffers containing defaults and settings
    INTEGER,          INTENT(IN) :: col_width, gap_width   !< column width (in characters), width of gap between columns
    INTEGER,          INTENT(IN), OPTIONAL     :: opt_file !< log file name (optional)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::parse_namelist'
    CHARACTER (len=MAX_STRBUF_LENGTH)  :: name, key, def_key
    CHARACTER (len=MAX_STRBUF_LENGTH)  :: def_value, value
    INTEGER                            :: ipos1,ipos2,iline1,iline2, dst_file

    dst_file = 0
    IF (PRESENT(opt_file))  dst_file = opt_file

    ! parse namelist
    iline1 = 1
    iline2 = 1
    ipos1  = 1
    CALL parse_nml_name(buffer1, ipos1, name)
    ipos2  = ipos1
    WRITE (dst_file,*) " "
    WRITE (dst_file,'(a,a)') "NAMELIST: ", TRIM(name)
    WRITE (dst_file,*) " "
    DO 
      IF (.NOT. parse_key_value_pair(buffer1, ipos1, iline1, startpos1, def_key, def_value)) EXIT
      IF (.NOT. parse_key_value_pair(buffer2, ipos2, iline2, startpos2, key,     value))     EXIT
      IF (TRIM(def_key) /= TRIM(key)) THEN
        WRITE (0,*) "key = ", TRIM(key), "; default key = ", TRIM(def_key)
        CALL finish(routine, "Internal error!")
      END IF
      IF (def_value == value) THEN
        WRITE (dst_file,'(a,a, a'//int2string(gap_width)//',a)') "  ", def_key(1:col_width), " ", TRIM(def_value(1:col_width))
      ELSE
        WRITE (dst_file,'(a,a, a'//int2string(gap_width)//',a,a,a,a)') "  ", def_key(1:col_width), " ", &
          &      value(1:col_width),  &
          &      '     (default was: ', TRIM(def_value), ')'
      END IF
    END DO
    WRITE (dst_file,*) " "
  END SUBROUTINE parse_namelist


  !> Opens a new temporary text file.
  !
  FUNCTION temp_defaults()
    INTEGER :: temp_defaults
    tmpnml1 = find_next_free_unit(10,100)
    OPEN (tmpnml1, file=TMP_FILENAME1//int2string(get_my_mpi_work_id()), &
      &   status='replace', action='write', delim='apostrophe')
    temp_defaults = tmpnml1
  END FUNCTION temp_defaults


  !> Opens a new temporary text file.
  !
  FUNCTION temp_settings()
    INTEGER :: temp_settings
    tmpnml2 = find_next_free_unit(10,100)
    OPEN (tmpnml2, file=TMP_FILENAME2//int2string(get_my_mpi_work_id()), &
      &   status='replace', action='write', delim='apostrophe')
    temp_settings = tmpnml2
  END FUNCTION temp_settings


  !> Loop over buffer, determine start indices of new namelist key/value
  !  pairs.
  !
  SUBROUTINE compute_startpos(buffer, startpos)
    CHARACTER (len=*),   INTENT(IN)    :: buffer       !< character buffer containing namelist
    INTEGER,             INTENT(INOUT) :: startpos(:)  !< array of key/value start positions
    ! local variables
    INTEGER   :: i, ilevel, last_pos, iline, i_max, gap_start
    CHARACTER :: str_delim(10) ! if level>0 : character delimiting string (' or ")
    LOGICAL   :: lgap

    i           = 0
    i_max       = LEN_TRIM(buffer)
    ilevel      = 0 ! ("0" means: outside string clauses)
    last_pos    = 1
    startpos(:) = 0
    startpos(1) = 1
    iline       = 2
    lgap        = .TRUE.
    gap_start   = 1
    ! loop over all characters of the string
    DO
      i = i + 1
      IF (i > i_max) EXIT

      ! check, if we open or close a string (which means, that we do not parse
      ! the current character any further
      IF ((buffer(i:i) == '"') .OR. (buffer(i:i) == "'")) THEN
        IF ((ilevel > 0) .AND. (str_delim(ilevel+1) == buffer(i:i))) THEN
          ilevel = ilevel - 1
        ELSE
          ilevel = ilevel + 1
          str_delim(ilevel+1) = buffer(i:i)
        END IF
      END IF
      IF (ilevel > 0) CYCLE

      ! ignore blanks (SPACE and TAB); ignore ","
      IF ((IACHAR(buffer(i:i)) == 9) .OR. (IACHAR(buffer(i:i)) == 32)  .OR.  &
        & (buffer(i:i) == ",")) THEN
        IF (.NOT. lgap) gap_start = i
        lgap = .TRUE.
        CYCLE
      END IF
        
      ! set "startpos" to the beginning of the last keyword
      IF (buffer(i:i) == "=") THEN
        startpos(iline) = last_pos
        iline = iline + 1
      ELSE
        IF (lgap) THEN
          last_pos = i
          lgap     = .FALSE.
        END IF
      END IF
    END DO
  END SUBROUTINE compute_startpos


  !> Read namelist back into string buffer
  !
  SUBROUTINE nml_read_defaults()
    ! local variables
    INTEGER :: ipos,iline, ios

    ! close previously opened temporary text file with namelist
    CLOSE (tmpnml1)
    tmpnml1 = find_next_free_unit(10,100)
    OPEN (tmpnml1, file=TMP_FILENAME1//int2string(get_my_mpi_work_id()), status='old', action='read')
    ipos  = 1
    iline = 1
    DO
      READ(tmpnml1, '(A)', iostat=ios, END=888) nml_log_buffer1(ipos:)
      IF (ios /= 0) EXIT
      nml_log_buffer1 = TRIM(nml_log_buffer1)
      startpos1(iline) = ipos
      ipos  = LEN(TRIM(nml_log_buffer1)) + 1
      iline = iline + 1
    END DO
888 CONTINUE
    CLOSE (tmpnml1, status='delete')
    CALL compute_startpos(nml_log_buffer1, startpos1)
  END SUBROUTINE nml_read_defaults


  !> Read namelist with settings back into string buffer.
  !
  SUBROUTINE nml_read_settings()
    ! local variables
    INTEGER :: ipos,iline

    ! close previously opened temporary text file with namelist
    CLOSE (tmpnml2)
    tmpnml2 = find_next_free_unit(10,100)
    OPEN (tmpnml2, file=TMP_FILENAME2//int2string(get_my_mpi_work_id()), status='old', action='read')
    ipos  = 1
    iline = 1
    DO
      READ(tmpnml2, '(A)', END=888) nml_log_buffer2(ipos:)
      nml_log_buffer2 = TRIM(nml_log_buffer2)
      ipos  = LEN(TRIM(nml_log_buffer2)) + 1
    END DO
888 CONTINUE
    CLOSE (tmpnml2, status='delete')
    CALL compute_startpos(nml_log_buffer2, startpos2)
  END SUBROUTINE nml_read_settings


  !> Read defaults and settings into string buffer, compare them and do the
  !  print-out.
  !
  SUBROUTINE log_nml_settings(opt_filename)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: opt_filename !< log file name (optional)
    ! local variables
    INTEGER :: dst_file

    IF (get_my_mpi_work_id() /= 0) THEN
      CLOSE (tmpnml1, status='delete')
      CLOSE (tmpnml2, status='delete')
      RETURN
    END IF

    CALL nml_read_defaults()
    CALL nml_read_settings()

    dst_file = 0
    IF (PRESENT(opt_filename)) THEN
      dst_file = find_next_free_unit(10,100)
      IF (.NOT. linitialized) THEN
        OPEN (dst_file, file=opt_filename, status='replace', action='write')
        linitialized = .TRUE.
      ELSE
        OPEN (dst_file, file=opt_filename, position='append', action='write')
      END IF
    END IF
    ! parse namelist
    CALL parse_namelist(nml_log_buffer1, nml_log_buffer2, &
      &                 col_width=64, gap_width=3,        &
      &                 opt_file=dst_file)
    IF (PRESENT(opt_filename)) THEN
      CLOSE(dst_file)
    END IF
  END SUBROUTINE log_nml_settings

END MODULE mo_nml_annotate
