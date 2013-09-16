!>
!! Contains utility routines for CDI I/O (table summary).
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
MODULE mo_util_cdi_table

  USE mo_impl_constants, ONLY : MAX_CHAR_LENGTH
  USE mo_util_string,    ONLY : int2string
  USE mo_exception,      ONLY : finish


  IMPLICIT NONE
  INCLUDE 'cdi.inc'

  PRIVATE
  PUBLIC  :: print_cdi_summary

  ! ---------------------------------------------------------------------
  ! CONSTANTS
  ! ---------------------------------------------------------------------

  INTEGER, PARAMETER :: MAX_TABLE_COLUMNS = 100
  INTEGER, PARAMETER :: MAX_TITLE_LEN     =  64
  INTEGER, PARAMETER :: MAX_COLUMN_LEN    =  64

  INTEGER, PARAMETER :: GRIB2_SHORTNAME   =  1
  INTEGER, PARAMETER :: GRIB2_CATEGORY    =  2
  INTEGER, PARAMETER :: GRIB2_DISCIPLINE  =  3
  INTEGER, PARAMETER :: GRIB2_NUMBER      =  4
  INTEGER, PARAMETER :: GRIB2_TRIPLE      =  5
  INTEGER, PARAMETER :: VALIDITY_DATE     =  7
  INTEGER, PARAMETER :: VALIDITY_TIME     =  8
  INTEGER, PARAMETER :: LEVEL_TYPE        =  9
  INTEGER, PARAMETER :: RUN_TYPE          = 10
  INTEGER, PARAMETER :: TIME_VVMM         = 11
  INTEGER, PARAMETER :: NUM_LEVELS        = 13 
  INTEGER, PARAMETER :: RUN_CLASS         = 14 
  INTEGER, PARAMETER :: EXP_ID            = 15 
  INTEGER, PARAMETER :: GRID_ID           = 16
  INTEGER, PARAMETER :: NGRIDREF          = 17

  CHARACTER(LEN=*), PARAMETER :: DELIMITER     = ' | '

  CHARACTER(LEN=*), PARAMETER :: modname   = 'mo_util_cdi_table'

  ! ---------------------------------------------------------------------
  ! TYPE DEFINITIONS
  ! ---------------------------------------------------------------------

  !> Type definition for a single table column
  TYPE t_column
    CHARACTER(LEN=MAX_TITLE_LEN) :: title
    INTEGER                      :: id, width
  END TYPE t_column

  !> Type definition for a complete table
  TYPE t_table
    INTEGER                      :: n_columns
    TYPE (t_column)              :: column(MAX_TABLE_COLUMNS) 
  END TYPE t_table


CONTAINS

  !> Initialize the table layout (title, columns).
  !
  SUBROUTINE setup_table_output(table)
    TYPE (t_table), INTENT(INOUT) :: table

    table%n_columns = 12
    !                           title      column ID        width
    table%column( 1) = t_column("name",    GRIB2_SHORTNAME,  10)
    table%column( 2) = t_column("triple",  GRIB2_TRIPLE,     11)
    table%column( 3) = t_column("date",    VALIDITY_DATE,    10)
    table%column( 4) = t_column("time",    VALIDITY_TIME,     8)
    table%column( 5) = t_column("lvt",     LEVEL_TYPE,        3)
    table%column( 6) = t_column("runtype", RUN_TYPE,          7)
    table%column( 7) = t_column("vvmm",    TIME_VVMM,         4)
    table%column( 8) = t_column("nlv",     NUM_LEVELS,        3)
    table%column( 9) = t_column("clas",    RUN_CLASS,         4)
    table%column(10) = t_column("expid",   EXP_ID,            5)
    table%column(11) = t_column("grid",    GRID_ID,           5)
    table%column(12) = t_column("rgrid",   NGRIDREF,          5)
  END SUBROUTINE setup_table_output


  !> @return Table entry for a given variable and a given column.
  !
  FUNCTION get_table_entry(vlistID, ivar, column)
    CHARACTER(LEN=MAX_COLUMN_LEN) :: get_table_entry
    INTEGER,         INTENT(IN) :: vlistID
    TYPE (t_column), INTENT(IN) :: column
    INTEGER,         INTENT(IN) :: ivar
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::get_table_entry'
    CHARACTER (len=MAX_CHAR_LENGTH) :: name
    CHARACTER (len=32)              :: wdth
    INTEGER :: param, number, category, discipline, &
      &        zaxisID, nlev, iexp_id, ilevtyp,     &
      &        iruntype, irunclass, ingridused,     &
      &        ingridref, gridID

    if (column%width > MAX_COLUMN_LEN) &
      &  CALL finish(routine, "Internal error: Column width exceeds maximum length!")

    get_table_entry = " "
    wdth = TRIM(int2string(column%width))
    SELECT CASE (column%id)
    CASE(GRIB2_SHORTNAME)
      CALL vlistInqVarName(vlistID, ivar, name)
      WRITE (get_table_entry, "(a"//TRIM(wdth)//")") ADJUSTL(name(1:column%width))
      !
    CASE(GRIB2_CATEGORY)
      ! get GRIB2 information on variable:
      param = vlistInqVarParam(vlistID, ivar)
      CALL cdiDecodeParam(param, number, category, discipline)
      WRITE (get_table_entry, "(i"//TRIM(wdth)//")") category
      !
    CASE(GRIB2_DISCIPLINE)
      ! get GRIB2 information on variable:
      param = vlistInqVarParam(vlistID, ivar)
      CALL cdiDecodeParam(param, number, category, discipline)
      WRITE (get_table_entry, "(i"//TRIM(wdth)//")") discipline
      !
    CASE(GRIB2_NUMBER)
      ! get GRIB2 information on variable:
      param = vlistInqVarParam(vlistID, ivar)
      CALL cdiDecodeParam(param, number, category, discipline)
      WRITE (get_table_entry, "(i"//trim(wdth)//")") number
      !
    CASE(GRIB2_TRIPLE)
      ! get GRIB2 information on variable:
      param = vlistInqVarParam(vlistID, ivar)
      CALL cdiDecodeParam(param, number, category, discipline)
      WRITE (get_table_entry, "(i3.1,a1,i3.1,a1,i3.1)") discipline,".",category,".",number
      !
    CASE(VALIDITY_DATE)
      ! NOT YET IMPLEMENTED!
      get_table_entry = "n/a"
      !
    CASE(VALIDITY_TIME)
      ! NOT YET IMPLEMENTED!
      get_table_entry = "n/a"
      !
    CASE(LEVEL_TYPE)
      ilevtyp = vlistInqVarIntKey(vlistID, ivar, "typeOfFirstFixedSurface")
      WRITE (get_table_entry, "(i"//trim(wdth)//")") ilevtyp
      !
    CASE(RUN_TYPE)
      iruntype = vlistInqVarIntKey(vlistID, ivar, "typeOfGeneratingProcess")
      WRITE (get_table_entry, "(i"//trim(wdth)//")") iruntype
      !
    CASE(TIME_VVMM)
      ! NOT YET IMPLEMENTED!
      get_table_entry = "n/a"
      !
    CASE(NUM_LEVELS)
      zaxisID = vlistInqVarZaxis(vlistID, ivar)
      nlev    = zaxisInqSize(zaxisID)
      WRITE (get_table_entry, "(i"//trim(wdth)//")") nlev
      !
    CASE(RUN_CLASS)
      irunclass = vlistInqVarIntKey(vlistID, ivar, "backgroundProcess")
      WRITE (get_table_entry, "(i"//trim(wdth)//")") irunclass
      !
    CASE(EXP_ID)
      iexp_id = vlistInqVarIntKey(vlistID, ivar, "localNumberOfExperiment")
      WRITE (get_table_entry, "(i"//trim(wdth)//")") iexp_id
      !
    CASE(GRID_ID)
      gridID = vlistInqVarGrid(vlistID, ivar)
      ingridused = gridInqNumber(gridID)
      WRITE (get_table_entry, "(i"//trim(wdth)//")") ingridused
      !
    CASE(NGRIDREF)
      gridID = vlistInqVarGrid(vlistID, ivar)
      ingridref = gridInqPosition(gridID)
      WRITE (get_table_entry, "(i"//trim(wdth)//")") ingridref
      !
    CASE DEFAULT
      CALL finish(routine, "Internal error: Unknown table column!")
    END SELECT

  END FUNCTION get_table_entry


  !> Print a summary of a file opened with the CDI
  !
  !  @note We assume that the file has already been opened by the CDI.
  !
  SUBROUTINE print_cdi_summary(vlistID, opt_dstfile)
    INTEGER, INTENT(IN)           :: vlistID
    INTEGER, INTENT(IN), OPTIONAL :: opt_dstfile !< (optional) output file
    ! local variables
    INTEGER                       :: dst, nvars, varID, line_len, icol, width
    TYPE (t_table)                :: table
    CHARACTER(LEN=MAX_TITLE_LEN)  :: title_str
    CHARACTER(LEN=MAX_COLUMN_LEN) :: entry_str
    CHARACTER(LEN=64)             :: format_str

    dst = 0
    IF (PRESENT(opt_dstfile)) dst = opt_dstfile

    ! initialize the table layout (title, columns)  
    CALL setup_table_output(table)

    ! construct and print the table header
    WRITE (dst,*) "" ! new line
    line_len = 0
    DO icol = 1, table%n_columns
      title_str = table%column(icol)%title
      width = table%column(icol)%width
      format_str = "(a"//TRIM(int2string(width))//',a)'
      WRITE (dst,TRIM(format_str), advance='no') ADJUSTL(title_str(1:width)), DELIMITER
      line_len = line_len + table%column(icol)%width + LEN(DELIMITER)
    END DO
    WRITE (dst,*) "" ! new line
    WRITE (dst,*) "" ! new line

    ! write the table contents
    nvars = vlistNvars(vlistID)
    LOOP : DO varID=0,(nvars-1)
      DO icol = 1, table%n_columns
        entry_str = get_table_entry(vlistID, varID, table%column(icol))
        format_str = "(a"//TRIM(int2string(table%column(icol)%width))//',a'//")"
        WRITE (dst,TRIM(format_str), advance='no') entry_str, DELIMITER
      END DO
      WRITE (dst,*) "" ! new line
    END DO LOOP
    WRITE (dst,*) "" ! new line
  END SUBROUTINE print_cdi_summary

END MODULE mo_util_cdi_table
