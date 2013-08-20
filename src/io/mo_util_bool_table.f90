!>
!! Contains utility routines for a TRUE/FALSE table.
!!
!! The table has the following layout: Each table row is determined by
!! a row name. Then, in each column the table entry is marked by an
!! "x" where the column contains this row name.
!!
!! @author F. Prill, DWD
!!
!!
!! @par Revision History
!! Initial revision: 2013-08-20 : F. Prill, DWD
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
!! @todo Implement a destruction routine.
!!
MODULE mo_util_bool_table

  USE mo_impl_constants, ONLY : MAX_CHAR_LENGTH, SUCCESS
  USE mo_var_metadata,   ONLY : VARNAME_LEN
  USE mo_exception,      ONLY : finish
  USE mo_util_string,    ONLY : remove_duplicates, tolower, int2string


  IMPLICIT NONE

  PRIVATE
  PUBLIC  :: init_bool_table
  PUBLIC  :: add_column
  PUBLIC  :: print_bool_table
  PUBLIC  :: t_column
  PUBLIC  :: t_bool_table

  ! ---------------------------------------------------------------------
  ! CONSTANTS
  ! ---------------------------------------------------------------------

  INTEGER, PARAMETER :: MAX_TABLE_COLUMNS =  100          !< max num. of table columns
  INTEGER, PARAMETER :: MAX_TABLE_ROWS    = 1000          !< max num. of table rows
  INTEGER, PARAMETER :: MAX_TITLE_LEN     =   64          !< max length of column title
  INTEGER, PARAMETER :: MAX_ROWNAME_LEN   = VARNAME_LEN   !< max length of row name

  CHARACTER(LEN=*), PARAMETER :: DELIMITER     = ' | '    !< vertical line
  CHARACTER(LEN=*), PARAMETER :: CHAR_TRUE     = ' x '    !< TRUE  sign
  CHARACTER(LEN=*), PARAMETER :: CHAR_FALSE    = '   '    !< FALSE sign

  CHARACTER(LEN=*), PARAMETER :: modname   = 'mo_util_bool_table'

  ! ---------------------------------------------------------------------
  ! TYPE DEFINITIONS
  ! ---------------------------------------------------------------------

  !> Type definition for a single table column
  TYPE t_column
    CHARACTER(LEN=MAX_TITLE_LEN) :: title
    INTEGER                      :: width
    CHARACTER(LEN=MAX_ROWNAME_LEN), ALLOCATABLE :: true_entries(:)
  END TYPE t_column

  !> Type definition for a complete table
  TYPE t_bool_table
    INTEGER                         :: n_columns
    TYPE (t_column)                 :: column(MAX_TABLE_COLUMNS) 
    CHARACTER(LEN=MAX_ROWNAME_LEN)  :: rowname(MAX_TABLE_ROWS) 
    INTEGER                         :: n_rows
    INTEGER                         :: rowname_len                !< max. length of a row name
  END TYPE t_bool_table


CONTAINS

  !> Initialize the table data structure.
  !
  SUBROUTINE init_bool_table(table)
    TYPE (t_bool_table), INTENT(INOUT) :: table

    table%n_columns   = 0
    table%n_rows      = 0
    table%rowname_len = 0
  END SUBROUTINE init_bool_table


  !> Add a column to the table.
  !
  SUBROUTINE add_column(table, colname, str_list, nitems)
    TYPE (t_bool_table),        INTENT(INOUT) :: table       !< true/false table object
    CHARACTER(LEN=*),      INTENT(IN)    :: colname     !< column name
    CHARACTER(LEN=*),      INTENT(IN)    :: str_list(:) !< list of "true" entries
    INTEGER,               INTENT(IN)    :: nitems      !< length of given list
    ! local variables    
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::add_column'
    INTEGER :: ierrstat, i, rowname_len

    table%n_columns = table%n_columns + 1
    ALLOCATE(table%column(table%n_columns)%true_entries(nitems), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! set column title string
    IF (LEN_TRIM(colname) > MAX_TITLE_LEN) THEN
      CALL finish(routine, "title name length too long!")
    ELSE
      table%column(table%n_columns)%title = TRIM(colname)
      table%column(table%n_columns)%width = LEN_TRIM(colname)
    END IF
    ! fill list of row entries
    DO i=1,nitems
      rowname_len = LEN_TRIM(str_list(i))
      IF (rowname_len > MAX_ROWNAME_LEN) THEN
        CALL finish(routine, "row name length too long!")
      ELSE
        table%rowname_len = MAX(table%rowname_len, rowname_len)
        table%column(table%n_columns)%true_entries(i) = TRIM(tolower(TRIM(str_list(i))))

        ! add entry also to the list of rows (duplicates will be removed later)
        table%n_rows = table%n_rows + 1
        table%rowname(table%n_rows) = TRIM(tolower(TRIM(str_list(i))))
      END IF
    END DO
    ! throw out duplicates:
    CALL remove_duplicates(table%rowname, table%n_rows)
  END SUBROUTINE add_column


  !> @return string corresponding to TRUE, if @p rowname is set in
  !          column @p column, otherwise return the string
  !          corresponding to FALSE.
  !
  FUNCTION get_bool_table_entry(rowname, column)
    CHARACTER(LEN=MAX_ROWNAME_LEN) :: get_bool_table_entry
    CHARACTER(len=*), INTENT(IN) :: rowname
    TYPE (t_column),  INTENT(IN) :: column
    ! local variables
    INTEGER :: i
    
    get_bool_table_entry = CHAR_FALSE
    FIND_LOOP : DO i=1,SIZE(column%true_entries)
      IF (TRIM(tolower(TRIM(rowname))) == TRIM(tolower(TRIM(column%true_entries(i))))) THEN
        get_bool_table_entry = CHAR_TRUE
        EXIT FIND_LOOP
      END IF
    END DO FIND_LOOP
  END FUNCTION get_bool_table_entry


  !> Print the TRUE/FALSE table
  !
  SUBROUTINE print_bool_table(table)
    TYPE (t_bool_table),        INTENT(IN) :: table       !< true/false table object
    ! local variables
    INTEGER :: dst, i, width, icol
    CHARACTER(LEN=MAX_TITLE_LEN)   :: title_str
    CHARACTER(LEN=MAX_ROWNAME_LEN) :: entry_str
    CHARACTER(LEN=64)              :: format_str

    dst = 0 ! output file / output stream

    ! construct and print the table header
    WRITE (dst,*) "" ! new line
    ! first column:
    width      = table%rowname_len
    format_str = "(a"//TRIM(int2string(table%rowname_len))//',a)'
    WRITE (dst,TRIM(format_str), advance='no') " ", "    "
    ! other, TRUE/FALSE columns:
    DO icol = 1, table%n_columns
      title_str = table%column(icol)%title
      width     = table%column(icol)%width
      format_str = "(a"//TRIM(int2string(width))//',a)'
      WRITE (dst,TRIM(format_str), advance='no') ADJUSTL(title_str(1:width)), DELIMITER
    END DO
    WRITE (dst,*) "" ! new line
    WRITE (dst,*) "" ! new line

    ! write the table contents
    LOOP : DO i=1,table%n_rows
      ! first column:
      width      = table%rowname_len
      format_str = "(a"//TRIM(int2string(table%rowname_len))//',a)'
      WRITE (dst,TRIM(format_str), advance='no') table%rowname(i), "    "
      ! other, TRUE/FALSE columns:
      DO icol = 1, table%n_columns
        entry_str = get_bool_table_entry(table%rowname(i), table%column(icol))
        format_str = "(a"//TRIM(int2string(table%column(icol)%width))//',a'//")"
        WRITE (dst,TRIM(format_str), advance='no') entry_str, DELIMITER
      END DO
      WRITE (dst,*) "" ! new line
    END DO LOOP
    WRITE (dst,*) "" ! new line
  END SUBROUTINE print_bool_table

END MODULE mo_util_bool_table
