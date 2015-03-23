!>
!! Contains utility routines for CDI I/O (table summary).
!!
!! @author F. Prill, DWD
!!
!!
!! @par Revision History
!! Initial revision: 2013-07-03 : F. Prill, DWD
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_util_cdi_table

  USE mo_impl_constants,   ONLY : MAX_CHAR_LENGTH
  USE mo_util_string,      ONLY : int2string
  USE mo_exception,        ONLY : finish
  USE mo_util_uuid,        ONLY : t_uuid, char2uuid, uuid_unparse, uuid_string_length
  USE mo_util_hash,        ONLY : util_hashword
  USE mtime,               ONLY : newDatetime, newTimedelta, datetime, timedelta, &
    &                             timedeltaToString, max_timedelta_str_len,       &
    &                             max_datetime_str_len, deallocateDatetime,       &
    &                             deallocateTimedelta
  USE mo_mtime_extensions, ONLY : getTimeDeltaFromDateTime
  USE mo_dictionary,       ONLY : t_dictionary, dict_get


  IMPLICIT NONE
  INCLUDE 'cdi.inc'

  PRIVATE
  PUBLIC  :: print_cdi_summary
  PUBLIC  :: new_inventory_list
  PUBLIC  :: complete_inventory_list
  PUBLIC  :: delete_inventory_list
  PUBLIC  :: find_inventory_list_element

  PUBLIC  :: t_inventory_list
  PUBLIC  :: t_inventory_element

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
  INTEGER, PARAMETER :: H_UUID            = 18
  INTEGER, PARAMETER :: V_UUID            = 19

  CHARACTER(LEN=*), PARAMETER :: DELIMITER     = ' | '

  CHARACTER(LEN=*), PARAMETER :: modname   = 'mo_util_cdi_table'

  ! ---------------------------------------------------------------------
  ! TYPE DEFINITIONS
  ! ---------------------------------------------------------------------

  !> Type definition for a single table column
  TYPE t_column
    CHARACTER(LEN=MAX_TITLE_LEN) :: title
    INTEGER                      :: id, width
    LOGICAL                      :: lprint 
  END TYPE t_column

  !> Type definition for a complete table
  TYPE t_table
    INTEGER                      :: n_columns
    TYPE (t_column)              :: column(MAX_TABLE_COLUMNS)
  END TYPE t_table


  TYPE t_var_inventory_element
    CHARACTER(len=128)            :: name                     ! variable name
    INTEGER                       :: key                      ! hash value of name (for fast search)
    INTEGER                       :: discipline               ! GRIB2 discipline
    INTEGER                       :: category                 ! GRIB2 category
    INTEGER                       :: number                   ! GRIB2 number
    TYPE(datetime)                :: vdatetime                ! Validity date/time
    TYPE(timedelta)               :: forecast_time            ! forecast_time
    INTEGER                       :: levelType                ! level type
    INTEGER                       :: num_levels               ! number of levels
    INTEGER                       :: backgroundProc           ! GRIB2 background process
    INTEGER                       :: numberOfHGrid            ! GRIB2 numberOfGridUsed (horizontal)
    INTEGER                       :: gridInReference          ! GRIB2 numberOfGridInReference
    INTEGER                       :: typeOfGeneratingProcess  ! GRIB2 typeOfGeneratingProcess
    INTEGER                       :: typeOfFirstFixedSurface  ! typeOfFirstFixedSurface
    TYPE(t_uuid)                  :: uuidOfHGrid              ! horizontal grid UUID
    TYPE(t_uuid)                  :: uuidOfVGrid              ! vertical grid UUID
  END TYPE t_var_inventory_element

  TYPE t_inventory_element
    TYPE(t_var_inventory_element)       :: field
    TYPE(t_inventory_element), POINTER  :: next_list_element
  END TYPE t_inventory_element

  TYPE t_inventory_list_intrinsic
    INTEGER                            :: list_elements      ! allocated elements
    INTEGER                            :: filetype           ! CDI filetypeID
    TYPE(t_inventory_element), POINTER :: first_list_element ! reference to first 
  END TYPE t_inventory_list_intrinsic
  
  TYPE t_inventory_list
    TYPE(t_inventory_list_intrinsic), POINTER :: p
  END TYPE t_inventory_list

CONTAINS

  !> Initialize the table layout (title, columns).
  !
  SUBROUTINE setup_table_output(table)
    TYPE (t_table), INTENT(INOUT) :: table

    table%n_columns = 14
    !                           title            column ID        width   print
    table%column( 1) = t_column("name",          GRIB2_SHORTNAME,  10,    .TRUE. )
    table%column( 2) = t_column("triple",        GRIB2_TRIPLE,     11,    .TRUE. )
    table%column( 3) = t_column("Vdate",         VALIDITY_DATE,    10,    .TRUE. )
    table%column( 4) = t_column("Vtime",         VALIDITY_TIME,     8,    .TRUE. )
    table%column( 5) = t_column("lvt",           LEVEL_TYPE,        3,    .TRUE. )
    table%column( 6) = t_column("runtyp",        RUN_TYPE,          6,    .TRUE. )
    table%column( 7) = t_column("vvmm",          TIME_VVMM,         6,    .TRUE. )
    table%column( 8) = t_column("nlv",           NUM_LEVELS,        3,    .TRUE. )
    table%column( 9) = t_column("clas",          RUN_CLASS,         4,    .TRUE. )
    table%column(10) = t_column("expid",         EXP_ID,            5,    .TRUE. )
    table%column(11) = t_column("grid",          GRID_ID,           5,    .TRUE. )
    table%column(12) = t_column("rgrid",         NGRIDREF,          5,    .TRUE. )
    table%column(13) = t_column("uuidOfHGrid",   H_UUID,           11,    .FALSE.)
    table%column(14) = t_column("uuidOfVGrid",   V_UUID,           11,    .FALSE.)
  END SUBROUTINE setup_table_output


  !> @return Table entry for a given variable and a given column.
  !  @return: Store table entry in linked list (optional)
  !
  SUBROUTINE get_table_entry(vlistID, ivar, column, entry_str, opt_list_element)
    INTEGER,         INTENT(IN) :: vlistID
    TYPE (t_column), INTENT(IN) :: column
    INTEGER,         INTENT(IN) :: ivar
    CHARACTER(LEN=MAX_COLUMN_LEN),         INTENT(OUT) :: entry_str
    TYPE(t_inventory_element), OPTIONAL, INTENT(INOUT) :: opt_list_element
    !
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::get_table_entry'
    CHARACTER (len=MAX_CHAR_LENGTH) :: name
    CHARACTER (len=32)              :: wdth
    INTEGER :: param, number, category, discipline, &
      &        zaxisID, nlev, iexp_id, ilevtyp,     &
      &        iruntype, irunclass, ingridused,     &
      &        ingridref, gridID, taxisID
    INTEGER :: vdate, vtime, rdate, rtime,          &
      &        hour, minute, second,                &
      &        year, month, day
    TYPE(t_uuid)             :: uuidOfHGrid, uuidOfVGrid
    CHARACTER(len=1)         :: uuid_string(16)
    CHARACTER(len=uuid_string_length) :: uuid_unparsed
    TYPE(datetime),  POINTER :: mtime_vdatetime, mtime_rdatetime
    TYPE(timedelta), POINTER :: forecast_time
    CHARACTER(len=max_timedelta_str_len) :: forecast_time_string

    if (column%width > MAX_COLUMN_LEN) &
      &  CALL finish(routine, "Internal error: Column width exceeds maximum length!")

    entry_str = " "
    wdth = TRIM(int2string(column%width))
    SELECT CASE (column%id)
    CASE(GRIB2_SHORTNAME)
      CALL vlistInqVarName(vlistID, ivar, name)
      WRITE (entry_str, "(a"//TRIM(wdth)//")") ADJUSTL(name(1:column%width))
      IF (PRESENT(opt_list_element)) opt_list_element%field%name = TRIM(name)
      !
    CASE(GRIB2_CATEGORY)
      ! get GRIB2 information on variable:
      param = vlistInqVarParam(vlistID, ivar)
      CALL cdiDecodeParam(param, number, category, discipline)
      WRITE (entry_str, "(i"//TRIM(wdth)//")") category
      IF (PRESENT(opt_list_element)) opt_list_element%field%category = category
      !
    CASE(GRIB2_DISCIPLINE)
      ! get GRIB2 information on variable:
      param = vlistInqVarParam(vlistID, ivar)
      CALL cdiDecodeParam(param, number, category, discipline)
      WRITE (entry_str, "(i"//TRIM(wdth)//")") discipline
      IF (PRESENT(opt_list_element)) opt_list_element%field%discipline = discipline
      !
    CASE(GRIB2_NUMBER)
      ! get GRIB2 information on variable:
      param = vlistInqVarParam(vlistID, ivar)
      CALL cdiDecodeParam(param, number, category, discipline)
      WRITE (entry_str, "(i"//trim(wdth)//")") number
      IF (PRESENT(opt_list_element)) opt_list_element%field%number = number
      !
    CASE(GRIB2_TRIPLE)
      ! get GRIB2 information on variable:
      param = vlistInqVarParam(vlistID, ivar)
      CALL cdiDecodeParam(param, number, category, discipline)
      WRITE (entry_str, "(i3.1,a1,i3.1,a1,i3.1)") discipline,".",category,".",number
      !
    CASE(VALIDITY_DATE)
      taxisID = vlistInqTaxis(vlistID)
      vdate   = taxisInqVdate(taxisID)
      CALL cdiDecodeDate(vdate, year, month, day)
      WRITE(entry_str, "(i4,a1,i2.2,a1,i2.2)") year,"-",month,"-",day
      !
    CASE(VALIDITY_TIME)
      taxisID = vlistInqTaxis(vlistID)
      vtime   = taxisInqVtime(taxisID)
      CALL cdiDecodeTime(vtime, hour, minute, second)
      WRITE(entry_str, "(i2,a1,i2.2,a1,i2.2)") hour,":",minute,":",second
      !
    CASE(LEVEL_TYPE)
      ilevtyp = vlistInqVarIntKey(vlistID, ivar, "typeOfFirstFixedSurface")
      WRITE (entry_str, "(i"//trim(wdth)//")") ilevtyp
      IF (PRESENT(opt_list_element)) opt_list_element%field%levelType = ilevtyp
      !
    CASE(RUN_TYPE)
      iruntype = vlistInqVarTypeOfGeneratingProcess(vlistID, ivar)
      WRITE (entry_str, "(i"//trim(wdth)//")") iruntype
      IF (PRESENT(opt_list_element)) opt_list_element%field%typeOfGeneratingProcess = iruntype
      !
    CASE(TIME_VVMM)
      taxisID = vlistInqTaxis(vlistID)
      ! get vdate and vtime
      vdate   = taxisInqVdate(taxisID)
      vtime   = taxisInqVtime(taxisID)
      CALL cdiDecodeTime(vdate, year, month, day)
      CALL cdiDecodeTime(vtime, hour, minute, second)
      mtime_vdatetime => newDatetime(year, month , day,  &
        &                            hour, minute, second, ims=0)
      ! get rdate and rtime
      rdate   = taxisInqRdate(taxisID)
      rtime   = taxisInqRtime(taxisID)
      CALL cdiDecodeTime(rdate, year, month, day)
      CALL cdiDecodeTime(rtime, hour, minute, second)
      mtime_rdatetime => newDatetime(year, month , day,  &
        &                            hour, minute, second, ims=0)

      ! this 'initialization' is necessary, in order to correctly deal with 
      ! timedelta=0.
      forecast_time =>newTimedelta("PT00H")
      ! compute forecastTime
      CALL getTimeDeltaFromDateTime(mtime_vdatetime, mtime_rdatetime, forecast_time)
      CALL timedeltaToString(forecast_time, forecast_time_string)

      entry_str = TRIM(forecast_time_string)
      !
      IF (PRESENT(opt_list_element)) THEN
        opt_list_element%field%vdatetime     = mtime_vdatetime
        opt_list_element%field%forecast_time = forecast_time
      ENDIF
      !
      CALL deallocateDatetime(mtime_vdatetime)
      CALL deallocateDatetime(mtime_rdatetime)
      CALL deallocateTimedelta(forecast_time)
      !
    CASE(NUM_LEVELS)
      zaxisID = vlistInqVarZaxis(vlistID, ivar)
      nlev    = zaxisInqSize(zaxisID)
      WRITE (entry_str, "(i"//trim(wdth)//")") nlev
      IF (PRESENT(opt_list_element)) opt_list_element%field%num_levels = nlev
      !
    CASE(RUN_CLASS)
      irunclass = vlistInqVarIntKey(vlistID, ivar, "backgroundProcess")
      WRITE (entry_str, "(i"//trim(wdth)//")") irunclass
      IF (PRESENT(opt_list_element)) opt_list_element%field%backgroundProc = irunclass
      !
    CASE(EXP_ID)
      iexp_id = vlistInqVarIntKey(vlistID, ivar, "localNumberOfExperiment")
      WRITE (entry_str, "(i"//trim(wdth)//")") iexp_id
      !
    CASE(GRID_ID)
      gridID = vlistInqVarGrid(vlistID, ivar)
      ingridused = gridInqNumber(gridID)
      WRITE (entry_str, "(i"//trim(wdth)//")") ingridused
      IF (PRESENT(opt_list_element)) opt_list_element%field%numberOfHGrid = ingridused
      !
    CASE(NGRIDREF)
      gridID = vlistInqVarGrid(vlistID, ivar)
      ingridref = gridInqPosition(gridID)
      WRITE (entry_str, "(i"//trim(wdth)//")") ingridref
      !
    CASE(H_UUID)
      gridID = vlistInqVarGrid(vlistID, ivar)
      CALL gridInqUUID(gridID, uuid_string)
      CALL char2uuid(uuid_string, uuidOfHGrid)      ! Char array to uuid of type t_uuid
      CALL uuid_unparse(uuidOfHGrid, uuid_unparsed) ! uuid of type t_uuid to character string
      entry_str = TRIM(uuid_unparsed)
      IF (PRESENT(opt_list_element)) opt_list_element%field%uuidOfHGrid = uuidOfHGrid
      !
    CASE(V_UUID)
      !
      zaxisID = vlistInqVarZaxis(vlistID, ivar)
      CALL zaxisInqUUID(zaxisID, uuid_string)
      CALL char2uuid(uuid_string, uuidOfVGrid)      ! Char array to uuid of type t_uuid
      CALL uuid_unparse(uuidOfVGrid, uuid_unparsed) ! uuid of type t_uuid to character string
      entry_str = TRIM(uuid_unparsed)
      IF (PRESENT(opt_list_element)) opt_list_element%field%uuidOfVGrid = uuidOfVGrid
      !
    CASE DEFAULT
      CALL finish(routine, "Internal error: Unknown table column!")
    END SELECT

  END SUBROUTINE get_table_entry


  !> Print a summary of a file opened with the CDI
  !
  !  @note We assume that the file has already been opened by the CDI.
  !
  SUBROUTINE print_cdi_summary(vlistID, opt_dstlist, opt_dstfile)
    INTEGER, INTENT(IN)           :: vlistID
    INTEGER, INTENT(IN), OPTIONAL :: opt_dstfile !< (optional) output file
    TYPE(t_inventory_list), OPTIONAL, INTENT(INOUT) :: &
      &    opt_dstlist                           !< optional inventory list (linked list)
                                                 !< for storing the retrieved meta-info
    ! local variables
    INTEGER                       :: dst, nvars, varID, line_len, icol, width
    TYPE (t_table)                :: table
    CHARACTER(LEN=MAX_TITLE_LEN)  :: title_str
    CHARACTER(LEN=MAX_COLUMN_LEN) :: entry_str
    CHARACTER(LEN=64)             :: format_str

    TYPE(t_inventory_element), POINTER  :: current_list_element =>NULL()
    TYPE(t_inventory_element), POINTER  :: new_list_element =>NULL()

    dst = 0
    IF (PRESENT(opt_dstfile)) dst = opt_dstfile

    ! initialize the table layout (title, columns)  
    CALL setup_table_output(table)

    ! construct and print the table header
    WRITE (dst,*) "" ! new line
    line_len = 0
    DO icol = 1, table%n_columns
      IF (.NOT. table%column(icol)%lprint) CYCLE
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

      ! create new list element
      !
      IF (PRESENT(opt_dstlist)) THEN
        ! create new list element and append
        CALL append_inventory_list_element(opt_dstlist, new_list_element)
        current_list_element => new_list_element
      ENDIF

      DO icol = 1, table%n_columns

        CALL get_table_entry(vlistID, varID, table%column(icol), & ! in
          &                  entry_str,                          & ! in
          &                  opt_list_element=current_list_element ) ! inout

        format_str = "(a"//TRIM(int2string(table%column(icol)%width))//',a'//")"
        IF (table%column(icol)%lprint) THEN
          WRITE (dst,TRIM(format_str), advance='no') entry_str, DELIMITER
        ENDIF
      END DO

      WRITE (dst,*) "" ! new line
    END DO LOOP
    WRITE (dst,*) "" ! new line
  END SUBROUTINE print_cdi_summary



  !-------------
  !>
  !! SUBROUTINE new_inventory_list
  !! Initialize a variable of type t_inventory_list with default values 
  !! and nullify anchor to linked list
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2014-07-25)
  !!
  !!
  SUBROUTINE new_inventory_list(this_list)
    TYPE(t_inventory_list), INTENT(INOUT) :: this_list

    !------------------------------------------------------------

    ALLOCATE(this_list%p)

    this_list%p%list_elements       = 0
    this_list%p%filetype            = -999
    this_list%p%first_list_element  => NULL()
  END SUBROUTINE new_inventory_list



  !-------------
  !>
  !! SUBROUTINE create_inventory_list_element
  !! Allocate additional element for inventory list
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2014-07-25)
  !!
  !!
  SUBROUTINE create_inventory_list_element(this_list, current_list_element)
    !
    TYPE(t_inventory_list), INTENT(INOUT)    :: this_list
    TYPE(t_inventory_element), POINTER, INTENT(INOUT) :: current_list_element

    INTEGER :: ist
    !------------------------------------------------------------

    ALLOCATE (current_list_element, STAT=ist)
    IF (ist /= 0) THEN
      CALL finish('create_inventory_list_element','Cannot add element to linked list ...')
    ENDIF
    this_list%p%list_elements = this_list%p%list_elements + 1 
    !
    current_list_element%next_list_element => NULL()
    
  END SUBROUTINE create_inventory_list_element



  !-------------
  !>
  !! SUBROUTINE append_inventory_list_element
  !! Append additional element to inventory list
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2014-07-25)
  !!
  !!
  SUBROUTINE append_inventory_list_element(this_list, new_list_element)
    !
    TYPE(t_inventory_list),             INTENT(INOUT) :: this_list
    TYPE(t_inventory_element), POINTER, INTENT(INOUT) :: new_list_element

    TYPE(t_inventory_element), POINTER :: current_list_element
    !------------------------------------------------------------

    !
    ! insert as first element if list is empty
    !
    IF (.NOT. ASSOCIATED(this_list%p%first_list_element)) THEN
      CALL create_inventory_list_element(this_list, this_list%p%first_list_element)
      new_list_element => this_list%p%first_list_element
      RETURN
    ENDIF
    !
    ! loop over list elements to find position
    !
    current_list_element => this_list%p%first_list_element
    DO WHILE (ASSOCIATED(current_list_element%next_list_element)) 
      current_list_element => current_list_element%next_list_element
    ENDDO
    !
    ! insert element
    !
    CALL create_inventory_list_element (this_list, new_list_element)
    new_list_element%next_list_element => current_list_element%next_list_element
    current_list_element%next_list_element => new_list_element

  END SUBROUTINE append_inventory_list_element



  !-------------
  !>
  !! SUBROUTINE complete_inventory_list
  !! Names are translated into internal names, in case that GRIB2 fields 
  !! are read in. In addition, a hash value is created for those names 
  !! for search speed up.
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2014-07-28)
  !!
  !!
  SUBROUTINE complete_inventory_list(filetype, dict, this_list)
    INTEGER,                INTENT(IN)    :: filetype    ! CDI filetype ID
    TYPE(t_dictionary),     INTENT(IN)    :: dict        ! dictionary
    TYPE(t_inventory_list), INTENT(INOUT) :: this_list
    TYPE(t_inventory_element), POINTER    :: current_element

    ! store information on filetype
    this_list%p%filetype = filetype

    current_element => this_list%p%first_list_element

    DO WHILE (ASSOCIATED(current_element))
      ! In case of GRIB2 translate GRIB2 varname to netcdf varname and add to group
      IF (this_list%p%filetype == FILETYPE_GRB2) THEN
        ! GRIB2 -> Netcdf
        !
        current_element%field%name = TRIM(dict_get(dict,                       &
          &                          TRIM(current_element%field%name),         &
          &                          default=TRIM(current_element%field%name), &
          &                          linverse=.TRUE.))
      ENDIF
      !
      ! create hash value
      current_element%field%key = util_hashword(current_element%field%name, &
        &                         LEN_TRIM(current_element%field%name), 0)

      current_element => current_element%next_list_element
    ENDDO

  END SUBROUTINE complete_inventory_list



  !-------------
  !>
  !! SUBROUTINE delete_inventory_list
  !! Delete all elements of a linked list of type t_inventory_list 
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2014-07-25)
  !!
  !!
  SUBROUTINE delete_inventory_list(this_list)
    !
    TYPE(t_inventory_list), INTENT(INOUT) :: this_list

    !------------------------------------------------------------

    IF (ASSOCIATED(this_list%p)) THEN
      CALL delete_inventory_list_elements(this_list, this_list%p%first_list_element)

      this_list%p%first_list_element  => NULL()
      DEALLOCATE(this_list%p)
    ENDIF
  END SUBROUTINE delete_inventory_list



  !-------------
  !>
  !! SUBROUTINE delete_inventory_list_elements
  !! Deallocate a list element and all its successors 
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2014-07-25)
  !!
  !!
  SUBROUTINE delete_inventory_list_elements(this_list, this_list_element)
    TYPE(t_inventory_list), INTENT(INOUT)   :: this_list
    TYPE(t_inventory_element), POINTER      :: this_list_element
    !
    TYPE(t_inventory_element), POINTER      :: this, next
    !------------------------------------------------------------
    next => this_list_element
    this_list_element => NULL()
    !
    DO
      IF (.NOT. ASSOCIATED(next)) EXIT
      this => next
      next => this%next_list_element
      !
      this_list%p%list_elements = this_list%p%list_elements-1
      DEALLOCATE (this)
     END DO

  END SUBROUTINE delete_inventory_list_elements



  !-------------
  !>
  !! SUBROUTINE find_inventory_list_element
  !! Find element with matching name 
  !!
  !! @par Revision History
  !! Initial version by Luis Kornblueh, MPI (2011-??-??)
  !! - Modification by Daniel Reinert, DWD  (2014-07-28)
  !!   Adapted to inventory lists
  !!
  !!
  FUNCTION find_inventory_list_element (this_list, name) RESULT(this_list_element)
    !
    TYPE(t_inventory_list),   INTENT(in) :: this_list
    CHARACTER(len=*),         INTENT(in) :: name
    !  
    TYPE(t_inventory_element), POINTER :: this_list_element
    INTEGER :: key
    !
    key = util_hashword(name, LEN_TRIM(name), 0)
    !
    this_list_element => this_list%p%first_list_element
    DO WHILE (ASSOCIATED(this_list_element))
      IF (key == this_list_element%field%key) THEN
        RETURN
      ENDIF
      this_list_element => this_list_element%next_list_element
    ENDDO
    !
    NULLIFY (this_list_element)
    !
  END FUNCTION find_inventory_list_element


END MODULE mo_util_cdi_table
