!>
!! Contains utility routines for simple table output.
!!
!! @author F. Prill, DWD
!!
!!
!! @par Revision History
!! Initial revision: 2013-07-03 : F. Prill, DWD
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
MODULE mo_util_table

  USE mo_impl_constants, ONLY : SUCCESS
  USE mo_exception,      ONLY : finish
  USE mo_util_string,    ONLY : int2string

  IMPLICIT NONE

  PRIVATE
  PUBLIC  :: initialize_table
  PUBLIC  :: finalize_table
  PUBLIC  :: add_table_column
  PUBLIC  :: set_table_entry
  PUBLIC  :: print_table
  ! types
  PUBLIC  :: t_table, t_column


  ! ---------------------------------------------------------------------
  ! CONSTANTS
  ! ---------------------------------------------------------------------

  INTEGER, PARAMETER :: MAX_TABLE_COLUMNS = 100
  INTEGER, PARAMETER :: MAX_TITLE_LEN     =  64
  INTEGER, PARAMETER :: MAX_COLUMN_LEN    =  64

  CHARACTER(LEN=*), PARAMETER :: modname   = 'mo_util_table'

  ! ---------------------------------------------------------------------
  ! TYPE DEFINITIONS
  ! ---------------------------------------------------------------------

  !> Type definition for a single table column
  TYPE t_column
    CHARACTER(LEN=MAX_TITLE_LEN)               :: title
    INTEGER                                    :: width
    INTEGER                                    :: n_rows
    CHARACTER(LEN=MAX_COLUMN_LEN), ALLOCATABLE :: row(:)
  END TYPE t_column

  !> Type definition for a complete table
  TYPE t_table
    INTEGER                                    :: n_columns, n_rows
    TYPE (t_column)                            :: column(MAX_TABLE_COLUMNS) 
  END TYPE t_table


CONTAINS

  !> Initialize the table data structure.
  !
  SUBROUTINE initialize_table(table)
    TYPE (t_table),   INTENT(INOUT) :: table
    CALL finalize_table(table)
  END SUBROUTINE initialize_table


  !> Deallocate table data structure.
  !
  SUBROUTINE finalize_table(table)
    TYPE (t_table),   INTENT(INOUT) :: table
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::finalize_table"
    INTEGER :: i, ierrstat
    
    DO i=1,SIZE(table%column)
      IF (ALLOCATED(table%column(i)%row)) THEN
        DEALLOCATE(table%column(i)%row, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')    
      END IF
    END DO

    table%n_columns = 0
    table%n_rows    = 0
  END SUBROUTINE finalize_table


  !> Initialize the table layout (title, columns).
  !
  SUBROUTINE add_table_column(table, column_title)
    TYPE (t_table),   INTENT(INOUT) :: table
    CHARACTER(LEN=*), INTENT(IN)    :: column_title

    table%n_columns = table%n_columns + 1
    table%column(table%n_columns)%title  = TRIM(column_title)
    table%column(table%n_columns)%width  = LEN_TRIM(column_title)
    table%column(table%n_columns)%n_rows = 0
    table%n_rows = 0
  END SUBROUTINE add_table_column


  !> Set column row size to a given size.
  !
  SUBROUTINE resize_column(column, n_rows)
    TYPE (t_column), INTENT(INOUT) :: column
    INTEGER,         INTENT(IN)    :: n_rows
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::resize_column"
    INTEGER :: new_size, errstat
    CHARACTER(LEN=MAX_COLUMN_LEN), TARGET, ALLOCATABLE :: tmp(:)

    ! triangle copy
    new_size = n_rows
    IF (ALLOCATED(column%row)) THEN
      IF (new_size <= column%n_rows) RETURN
      IF (column%n_rows > SIZE(column%row)) &
        & CALL finish (routine, 'Unexpected array size!')
      ALLOCATE(tmp(column%n_rows), STAT=errstat)
      IF (errstat /= 0) CALL finish (routine, 'Error in ALLOCATE operation!')
      tmp(1:column%n_rows) = column%row(1:column%n_rows)
      DEALLOCATE(column%row, STAT=errstat)
      IF (errstat /= 0) CALL finish (routine, 'Error in DEALLOCATE operation!')
    END IF
    ALLOCATE(column%row(new_size), STAT=errstat)
    IF (errstat /= 0)  CALL finish (routine, 'Error in ALLOCATE operation!')
    column%row((column%n_rows+1):new_size) = " "
    IF (ALLOCATED(tmp)) THEN
      column%row(1:column%n_rows) = tmp(1:column%n_rows)
      DEALLOCATE(tmp, STAT=errstat)
      IF (errstat /= 0)  CALL finish (routine, 'Error in DEALLOCATE operation!')
    END IF
    column%n_rows = new_size
  END SUBROUTINE resize_column


  FUNCTION get_column_index(table, column_title) RESULT(icol)
    INTEGER :: icol
    TYPE (t_table),    INTENT(IN) :: table
    CHARACTER(LEN=*),  INTENT(IN)    :: column_title
    ! local variables
    INTEGER :: i

    ! find the (first) matching column index
    i    =  1
    icol = -1
    COLFIND_LOOP : DO
      IF (i > table%n_columns)  EXIT COLFIND_LOOP
      IF (TRIM(column_title) == TRIM(table%column(i)%title)) THEN
        icol = i
        EXIT COLFIND_LOOP
      END IF
      i = i + 1
    END DO COLFIND_LOOP
  END FUNCTION get_column_index
  

  !> @return Set table entry in the given row and for a given column.
  !
  SUBROUTINE set_table_entry(table, irow, column_title, entry_str)
    TYPE (t_table),    INTENT(INOUT) :: table
    INTEGER,           INTENT(IN)    :: irow
    CHARACTER(LEN=*),  INTENT(IN)    :: column_title
    CHARACTER(LEN=*),  INTENT(IN)    :: entry_str
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::set_table_entry"
    INTEGER :: icol

    icol =  get_column_index(table, column_title)

    ! for convenience: create colums when called for the first row
    IF (icol == -1) THEN
      IF (irow == 1) THEN
        CALL add_table_column(table, column_title)
        icol =  get_column_index(table, column_title)
      ELSE
        CALL finish(routine, "Column title not found!")
      END IF
    END IF
    CALL resize_column(table%column(icol), irow)
    table%column(icol)%row(irow) = TRIM(entry_str)
    table%column(icol)%width     = MAX(table%column(icol)%width, LEN_TRIM(entry_str))
    table%n_rows = MAX(table%n_rows, irow)
  END SUBROUTINE set_table_entry


  !> Print a table.
  !
  SUBROUTINE print_table(table, opt_delimiter, opt_dstfile, opt_hline)
    TYPE (t_table),   INTENT(IN)           :: table
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: opt_delimiter   !< (optional) delimiter character
    INTEGER,          INTENT(IN), OPTIONAL :: opt_dstfile     !< (optional) output file
    LOGICAL,          INTENT(IN), OPTIONAL :: opt_hline       !< (optional) draw head/foot horizontal line

    ! local variables
    INTEGER                       :: dst, irow, icol, width
    CHARACTER(LEN=MAX_TITLE_LEN)  :: title_str
    CHARACTER(LEN=64)             :: format_str
    CHARACTER(LEN=MAX_COLUMN_LEN) :: entry_str
    CHARACTER(LEN=3)              :: delimiter
    LOGICAL                       :: hline

    hline = .FALSE.
    IF (PRESENT(opt_hline))  hline = opt_hline

    delimiter = ' | '
    IF (PRESENT(opt_delimiter))  delimiter = opt_delimiter
    
    dst = 0
    IF (PRESENT(opt_dstfile)) dst = opt_dstfile

    ! construct and print the table header
    WRITE (dst,*) "" ! new line

    ! optional head line
    IF (hline) THEN
      WRITE (dst,"(a1)", advance='no') " " ! indent
      DO icol = 1, table%n_columns
        width = table%column(icol)%width
        format_str = "(a"//TRIM(int2string(width))//',a)'
        WRITE (dst,TRIM(format_str), advance='no') REPEAT('-',width), delimiter
      END DO
      WRITE (dst,*) "" ! new line
    END IF

    WRITE (dst,"(a1)", advance='no') " " ! indent
    DO icol = 1, table%n_columns
      title_str = table%column(icol)%title
      width = MIN(table%column(icol)%width,MAX_TITLE_LEN)
      format_str = "(a"//TRIM(int2string(width))//',a)'
      WRITE (dst,TRIM(format_str), advance='no') ADJUSTL(title_str(1:width)), delimiter
    END DO
    WRITE (dst,*) "" ! new line
    WRITE (dst,"(a1)", advance='no') " " ! indent
    DO icol = 1, table%n_columns
      width = table%column(icol)%width
      format_str = "(a"//TRIM(int2string(width))//',a)'
      WRITE (dst,TRIM(format_str), advance='no') REPEAT('-',width), delimiter
    END DO
    WRITE (dst,*) "" ! new line
    WRITE (dst,*) "" ! new line

    ! write the table contents
    LOOP : DO irow=1,table%n_rows
      WRITE (dst,"(a1)", advance='no') " " ! indent
      DO icol = 1, table%n_columns
        entry_str  = " "
        IF (table%column(icol)%n_rows >= irow) THEN
          entry_str  = table%column(icol)%row(irow)
        END IF
        format_str = "(a"//TRIM(int2string(table%column(icol)%width))//',a'//")"
        WRITE (dst,TRIM(format_str), advance='no') entry_str, delimiter
      END DO
      WRITE (dst,*) "" ! new line
    END DO LOOP

    ! optional foot line
    IF (hline) THEN
      WRITE (dst,"(a1)", advance='no') " " ! indent
      DO icol = 1, table%n_columns
        width = table%column(icol)%width
        format_str = "(a"//TRIM(int2string(width))//',a)'
        IF (icol < table%n_columns) THEN
          WRITE (dst,TRIM(format_str), advance='no') REPEAT('-',width), "---"
        ELSE
          WRITE (dst,TRIM(format_str), advance='no') REPEAT('-',width), "   "
        END IF
      END DO
      WRITE (dst,*) "" ! new line
    END IF
  END SUBROUTINE print_table

END MODULE mo_util_table
