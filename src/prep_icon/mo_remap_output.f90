!! This module contains subroutines for file output of re-mapped data
!! using the CDI.
!!
!! @author F. Prill, DWD
!!
MODULE mo_remap_output

#ifdef __ICON__
  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_exception,          ONLY: finish
  USE mo_util_string,        ONLY: tolower
  USE mo_util_uuid,          ONLY: t_uuid
#else
  USE mo_utilities,     ONLY: wp, SUCCESS, tolower, finish, t_uuid
#endif

  USE mo_mpi,           ONLY: get_my_mpi_work_id
  USE mo_remap_config,  ONLY: dbg_level, MAX_INPUT_FIELDS, MAX_NZAXIS, &
    &                         MAX_NAME_LENGTH
  USE mo_remap_shared,  ONLY: t_grid, GRID_TYPE_REGULAR, GRID_TYPE_ICON
  USE mo_remap_input,   ONLY: t_field_metadata, t_global_metadata, t_zaxis_metadata, &
    &                         CELL_GRID, EDGE_GRID
  USE mo_remap_io,      ONLY: t_file_metadata

  IMPLICIT NONE
  INCLUDE 'cdi.inc'

  PRIVATE
  PUBLIC :: load_metadata_output
  PUBLIC :: open_output, close_output
  PUBLIC :: store_field
  ! output variables:
  PUBLIC :: varID
  PUBLIC :: t_output_grid

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_output')

  ! name and version of this tool
  ! (will be written as meta-data to output data file):
  CHARACTER (LEN=*), PARAMETER :: tool_name    = "prepicon_conversion_tool"
  CHARACTER (LEN=*), PARAMETER :: tool_version = "version $Rev$"

  ! generic interface for storing data to output file
  ! (may later be expanded to 3D fields)
  INTERFACE store_field
    MODULE PROCEDURE store_field2D
  END INTERFACE

  !> Data structure containing output grid description.
  !
  !  This meta-data is read from the provided output grid file.
  !
  TYPE t_output_grid
    TYPE(t_uuid)                  :: grid_uuid               !< uuid of grid
    ! unstructured cell grid description:
    INTEGER                       :: nglb_c                  !< global no. of cells
    INTEGER                       :: nglb_e                  !< global no. of edges
    INTEGER                       :: cell_type               !< 3: triangles
    REAL(wp), ALLOCATABLE         :: clon(:),    clat(:),   &
      &                              elon(:),    elat(:),   &
      &                              clonv(:,:), clatv(:,:),&
      &                              elonv(:,:), elatv(:,:)
    ! structured grid description:
    INTEGER                       :: structure               !< LONLAT or UNSTRUCTURED
    INTEGER                       :: nx, ny                  !< hor. grid dimensions
    REAL(wp), ALLOCATABLE         :: xvals(:), yvals(:)      !< regular horizontal grid
  END TYPE t_output_grid

  ! CDI axis ID's for output file
  INTEGER  :: cdiTaxisID, varID(MAX_INPUT_FIELDS), cdiZaxisID(MAX_NZAXIS)

CONTAINS

  !> Read all meta-data from the output grid file which will be later
  !  needed for the output data file.
  !
  SUBROUTINE load_metadata_output(file_metadata, grid_metadata, rank0)
    TYPE (t_file_metadata), INTENT(IN)  :: file_metadata
    TYPE (t_output_grid),   INTENT(OUT) :: grid_metadata
    INTEGER,                INTENT(IN)  :: rank0
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::load_metadata_output')
    INTEGER :: ierrstat, xsize, ysize
    ! CHARACTER :: dummy

    ! read data only on I/O process:
    IF (get_my_mpi_work_id() /= rank0) RETURN

    SELECT CASE(file_metadata%structure)
    CASE (GRID_TYPE_ICON)
      IF (dbg_level >= 10)  WRITE (0,*) "# load meta-data for output, ICON grid"
      ! read unstructured grid description from file:
      grid_metadata%cell_type = 3
      grid_metadata%structure = GRID_UNSTRUCTURED
      grid_metadata%nglb_c    = gridInqSize(file_metadata%c_gridID)
      grid_metadata%nglb_e    = gridInqSize(file_metadata%e_gridID)
      ALLOCATE(grid_metadata%clon(grid_metadata%nglb_c),     &
        &      grid_metadata%clat(grid_metadata%nglb_c),     &
        &      grid_metadata%elon(grid_metadata%nglb_e),     &
        &      grid_metadata%elat(grid_metadata%nglb_e),     &
        &      grid_metadata%clonv(grid_metadata%cell_type, grid_metadata%nglb_c), &
        &      grid_metadata%clatv(grid_metadata%cell_type, grid_metadata%nglb_c), &
        &      grid_metadata%elonv(4, grid_metadata%nglb_e), &
        &      grid_metadata%elatv(4, grid_metadata%nglb_e), &
        &      STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      xsize = gridInqXvals    (file_metadata%c_gridID, grid_metadata%clon)
      ysize = gridInqYvals    (file_metadata%c_gridID, grid_metadata%clat)
      xsize = gridInqXvals    (file_metadata%e_gridID, grid_metadata%elon)
      ysize = gridInqYvals    (file_metadata%e_gridID, grid_metadata%elat)
      xsize = gridInqXbounds  (file_metadata%c_gridID, grid_metadata%clonv)
      ysize = gridInqYbounds  (file_metadata%c_gridID, grid_metadata%clatv)
      xsize = gridInqXbounds  (file_metadata%e_gridID, grid_metadata%elonv)
      ysize = gridInqYbounds  (file_metadata%e_gridID, grid_metadata%elatv)
      ! read grid UUID:
      ! DISABLED for the time being (segfault problem):
      ! dummy = gridInqUUID(file_metadata%c_gridID, grid_metadata%grid_uuid%data)
    CASE (GRID_TYPE_REGULAR)
      IF (dbg_level >= 10)  WRITE (0,*) "# load meta-data for output, regular grid."
      ! read regular grid structure from file:
      grid_metadata%structure = GRID_LONLAT
      grid_metadata%cell_type = 4
      grid_metadata%nx        = gridInqXsize(file_metadata%c_gridID)
      grid_metadata%ny        = gridInqYsize(file_metadata%c_gridID)
      ALLOCATE(grid_metadata%xvals(grid_metadata%nx), &
        &      grid_metadata%yvals(grid_metadata%ny), &
        &      STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      xsize = gridInqXvals(file_metadata%c_gridID, grid_metadata%xvals)
      ysize = gridInqYvals(file_metadata%c_gridID, grid_metadata%yvals)
      ! read grid UUID:
      ! DISABLED for the time being (segfault problem):
      ! dummy = gridInqUUID(file_metadata%gridID, grid_metadata%grid_uuid%data)
    CASE DEFAULT
      CALL finish(routine, "Unknown grid type")
    END SELECT
  END SUBROUTINE load_metadata_output


  !> Open output data file.
  !
  SUBROUTINE open_output(out_filename, global_metadata, field_metadata, nfields, &
    &                    zaxis_metadata, n_zaxis, grid_metadata, file_metadata, &
    &                    infile_metadata)
    CHARACTER(LEN=*),         INTENT(IN)    :: out_filename
    TYPE (t_global_metadata), INTENT(IN)    :: global_metadata
    TYPE (t_field_metadata),  INTENT(IN)    :: field_metadata(:)
    TYPE (t_zaxis_metadata),  INTENT(IN)    :: zaxis_metadata(:)
    TYPE (t_output_grid),     INTENT(INOUT) :: grid_metadata
    INTEGER,                  INTENT(IN)    :: nfields, n_zaxis
    TYPE (t_file_metadata),   INTENT(OUT)   :: file_metadata
    TYPE (t_file_metadata),   INTENT(IN)    :: infile_metadata
    ! local constants
    CHARACTER(LEN=*), PARAMETER :: routine  = TRIM(modname)//'::open_output'
    CHARACTER(LEN=*), PARAMETER :: att_name = TRIM(tool_name)
    CHARACTER(LEN=*), PARAMETER :: att_txt  = TRIM(tool_version)
    LOGICAL,          PARAMETER :: ldefine = .TRUE.
    ! local variables
    INTEGER :: i, j, iret, cdiTimeIndex, iaxis, output_type, gridID, ngrids, tmp
    CHARACTER(len=MAX_NAME_LENGTH) :: zname

    cdiZaxisID(:) = CDI_UNDEFID
    output_type   = FILETYPE_NC2
    file_metadata%streamID  = streamOpenWrite(TRIM(out_filename), output_type)
    file_metadata%vlistID   = vlistCreate() ! create cdi vlist

    file_metadata%structure = infile_metadata%structure

    IF (ldefine) THEN
      ! define output grid (using CDI)
      CALL vlistDefInstitut(file_metadata%vlistID, global_metadata%cdiInstID)
      tmp = vlistDefAttTxt(file_metadata%vlistID, CDI_GLOBAL, &
        &    TRIM(att_name),LEN(TRIM(att_txt)),TRIM(att_txt))
      CALL define_output_grid(grid_metadata, file_metadata)
    ELSE
      ! copy output grid (using CDI)
      CALL vlistCopy(file_metadata%vlistID, infile_metadata%vlistID)

      SELECT CASE(file_metadata%structure)
      CASE (GRID_TYPE_ICON)
        ! set file gridID to cell-grid:
        ngrids = vlistNgrids(file_metadata%vlistID)
        LOOP : DO i=1,ngrids
          gridID = vlistGrid(file_metadata%vlistID, 0)
          CALL gridInqXname(gridID, zname)
          IF (tolower(TRIM(zname)) == "clon") file_metadata%c_gridID = gridID
          IF (tolower(TRIM(zname)) == "elon") file_metadata%e_gridID = gridID
          IF (tolower(TRIM(zname)) == "vlon") file_metadata%v_gridID = gridID
        END DO LOOP
      CASE (GRID_TYPE_REGULAR)
        ! set file gridID to regular grid:
        file_metadata%c_gridID = vlistGrid(file_metadata%vlistID, 0)
      CASE DEFAULT
        CALL finish(routine, "Unknown grid type")
      END SELECT
    END IF

    ! define vertical levels:
    ! -----------------------

    IF (dbg_level >= 5) WRITE (0,*) "# define levels"
    DO i=1,n_zaxis
      cdiZaxisID(i) = zaxisCreate(zaxis_metadata(i)%zaxis_type, zaxis_metadata(i)%zaxis_size)
      CALL zaxisDefLevels(cdiZaxisID(i), zaxis_metadata(i)%levels)
      IF (zaxis_metadata(i)%vct_size > 0) &
        &  CALL zaxisDefVct(cdiZaxisID(i), zaxis_metadata(i)%vct_size, zaxis_metadata(i)%vct)
    END DO

    ! define time axis:
    ! -----------------

    IF (dbg_level >= 5) WRITE (0,*) "# define time axis"
    cdiTaxisID = taxisCreate(global_metadata%timetype)
    CALL taxisDefTunit   (cdiTaxisID, global_metadata%tunit   )
    CALL taxisDefCalendar(cdiTaxisID, global_metadata%calendar)
    CALL taxisDefRdate   (cdiTaxisID, global_metadata%rdate   )
    CALL taxisDefRtime   (cdiTaxisID, global_metadata%rtime   )
    CALL taxisDefVdate   (cdiTaxisID, global_metadata%vdate   )
    CALL taxisDefVtime   (cdiTaxisID, global_metadata%vtime   )

    CALL vlistDefTaxis(file_metadata%vlistID, cdiTaxisID)

    ! define variables:
    ! -----------------

    IF (dbg_level >= 5) WRITE (0,*) "# define variables"
    DO i=1,nfields
      ! find z-axis for this variable
      iaxis = -1
      DO j=1,n_zaxis
        IF (zaxis_metadata(j)%zaxisID == field_metadata(i)%zaxisID) iaxis=j
      END DO
      IF (iaxis == -1)  CALL finish(routine, "Internal error!")
      ! create variable
      SELECT CASE (field_metadata(i)%grid_type)
      CASE(CELL_GRID)
        gridID = file_metadata%c_gridID
      CASE(EDGE_GRID)
        gridID = file_metadata%e_gridID
      CASE DEFAULT
        CALL finish(routine, "Unknown grid type")
      END SELECT
        
      varID(i) = vlistDefVar(file_metadata%vlistID, gridID, cdizaxisID(iaxis), &
        &                    field_metadata(i)%steptype)
      CALL vlistDefVarName(file_metadata%vlistID, varID(i), TRIM(field_metadata(i)%outputname))

      ! set GRIB2 triplet
      CALL vlistDefVarParam(file_metadata%vlistID, varID(i), &
        &  cdiEncodeParam(field_metadata(i)%grib2%number,    &
        &                 field_metadata(i)%grib2%category,  &
        &                 field_metadata(i)%grib2%discipline) )
      CALL vlistDefVarDatatype(file_metadata%vlistID, varID(i), field_metadata(i)%grib2%bits)

      CALL vlistDefVarStdName  (file_metadata%vlistID, varID(i), TRIM(field_metadata(i)%cf%standard_name))
      CALL vlistDefVarLongname (file_metadata%vlistID, varID(i), TRIM(field_metadata(i)%cf%long_name)    )
      CALL vlistDefVarUnits    (file_metadata%vlistID, varID(i), TRIM(field_metadata(i)%cf%units)        )
      CALL vlistDefVarTsteptype(file_metadata%vlistID, varID(i), field_metadata(i)%steptype              )
      CALL vlistDefVarDatatype (file_metadata%vlistID, varID(i), field_metadata(i)%cf%datatype           )

      ! set missing value information on variable:
      CALL vlistDefVarMissval(file_metadata%vlistID, varID(i), field_metadata(i)%missval)
    END DO

    CALL streamDefVlist(file_metadata%streamID, file_metadata%vlistID)

    cdiTimeIndex = 0
    iret = streamDefTimestep(file_metadata%streamID, cdiTimeIndex)
  END SUBROUTINE open_output


  !> Define output grid in data file (using CDI)
  !
  SUBROUTINE define_output_grid(grid_metadata, file_metadata)
    TYPE (t_output_grid),   INTENT(INOUT) :: grid_metadata
    TYPE (t_file_metadata), INTENT(INOUT) :: file_metadata
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::define_output_grid'
    INTEGER :: nx, ny

    SELECT CASE(grid_metadata%structure)
    CASE (GRID_LONLAT)
      ! define regular horizontal grid:
      nx = grid_metadata%nx
      ny = grid_metadata%ny
      IF (dbg_level >= 10) WRITE (0,*) "# define horizontal grid with size ", nx*ny

      file_metadata%c_gridID = gridCreate(grid_metadata%structure, nx * ny)
      CALL gridDefXsize (file_metadata%c_gridID, nx)
      CALL gridDefXname (file_metadata%c_gridID, 'lon')
      CALL gridDefXunits(file_metadata%c_gridID, 'degrees_east')
      CALL gridDefYsize (file_metadata%c_gridID, ny)
      CALL gridDefYname (file_metadata%c_gridID, 'lat')
      CALL gridDefYunits(file_metadata%c_gridID, 'degrees_north')

      CALL gridDefXvals (file_metadata%c_gridID, grid_metadata%xvals)
      CALL gridDefYvals (file_metadata%c_gridID, grid_metadata%yvals)
      ! define grid UUID:
      ! DISABLED for the time being (segfault problem):
      ! CALL gridDefUUID(file_metadata%c_gridID, grid_metadata%grid_uuid%data)
    CASE (GRID_UNSTRUCTURED)
      ! define unstructured horizontal cell grid:
      IF (dbg_level >= 10) WRITE (0,*) "# define horizontal grid with size ", grid_metadata%nglb_c

      ! cell grid:
      file_metadata%c_gridID = gridCreate(grid_metadata%structure, grid_metadata%nglb_c)
      CALL gridDefNvertex  (file_metadata%c_gridID, grid_metadata%cell_type)
      CALL gridDefXname    (file_metadata%c_gridID, 'clon')
      CALL gridDefXlongname(file_metadata%c_gridID, 'center longitude')
      CALL gridDefXunits   (file_metadata%c_gridID, 'degrees east')
      CALL gridDefYname    (file_metadata%c_gridID, 'clat')
      CALL gridDefYlongname(file_metadata%c_gridID, 'degrees north')
      CALL gridDefYunits   (file_metadata%c_gridID, 'radian')

      CALL gridDefXvals    (file_metadata%c_gridID, grid_metadata%clon)
      CALL gridDefYvals    (file_metadata%c_gridID, grid_metadata%clat)
      CALL gridDefXbounds  (file_metadata%c_gridID, grid_metadata%clonv)
      CALL gridDefYbounds  (file_metadata%c_gridID, grid_metadata%clatv)

      ! edge grid:
      file_metadata%e_gridID = gridCreate(grid_metadata%structure, grid_metadata%nglb_e)
      CALL gridDefNvertex  (file_metadata%e_gridID, 4)
      CALL gridDefXname    (file_metadata%e_gridID, 'elon')
      CALL gridDefXlongname(file_metadata%e_gridID, 'edge midpoint longitude')
      CALL gridDefXunits   (file_metadata%e_gridID, 'radian')
      CALL gridDefYname    (file_metadata%e_gridID, 'elat')
      CALL gridDefYlongname(file_metadata%e_gridID, 'edge midpoint latitude')
      CALL gridDefYunits   (file_metadata%e_gridID, 'radian')

      CALL gridDefXvals    (file_metadata%e_gridID, grid_metadata%elon)
      CALL gridDefYvals    (file_metadata%e_gridID, grid_metadata%elat)
      CALL gridDefXbounds  (file_metadata%e_gridID, grid_metadata%elonv)
      CALL gridDefYbounds  (file_metadata%e_gridID, grid_metadata%elatv)

      ! define grid UUID:
      ! DISABLED for the time being (segfault problem):
      ! CALL gridDefUUID(file_metadata%gridID, grid_metadata%grid_uuid%data)
    CASE DEFAULT
      CALL finish(routine, "Unknown grid type")
    END SELECT

    IF (dbg_level >= 10) WRITE (0,*) "# done."
  END SUBROUTINE define_output_grid


  !> Close output data file
  !
  SUBROUTINE close_output(file_metadata, grid_metadata)
    TYPE (t_file_metadata), INTENT(INOUT) :: file_metadata
    TYPE (t_output_grid),   INTENT(INOUT) :: grid_metadata
    ! local variables
    INTEGER :: j
    TYPE (t_file_metadata)      :: empty

    IF (dbg_level >= 5) WRITE (0,*) "# clean up grid"

    IF(file_metadata%c_gridID /= CDI_UNDEFID) CALL gridDestroy(file_metadata%c_gridID)
    IF(file_metadata%e_gridID /= CDI_UNDEFID) CALL gridDestroy(file_metadata%e_gridID)
    file_metadata%c_gridID = CDI_UNDEFID
    file_metadata%e_gridID = CDI_UNDEFID
    CALL taxisDestroy(cdiTaxisID)
    DO j = 1, SIZE(cdiZaxisID)
      IF(cdiZaxisID(j) /= CDI_UNDEFID) CALL zaxisDestroy(cdiZaxisID(j))
      cdiZaxisID(j) = CDI_UNDEFID
    ENDDO
    if (ALLOCATED (grid_metadata%clon )) DEALLOCATE (grid_metadata%clon)
    if (ALLOCATED (grid_metadata%clat )) DEALLOCATE (grid_metadata%clat)
    if (ALLOCATED (grid_metadata%elon )) DEALLOCATE (grid_metadata%elon)
    if (ALLOCATED (grid_metadata%elat )) DEALLOCATE (grid_metadata%elat)
    if (ALLOCATED (grid_metadata%clonv)) DEALLOCATE (grid_metadata%clonv)
    if (ALLOCATED (grid_metadata%clatv)) DEALLOCATE (grid_metadata%clatv)
    if (ALLOCATED (grid_metadata%elonv)) DEALLOCATE (grid_metadata%elonv)
    if (ALLOCATED (grid_metadata%elatv)) DEALLOCATE (grid_metadata%elatv)
    if (ALLOCATED (grid_metadata%xvals)) DEALLOCATE (grid_metadata%xvals)
    if (ALLOCATED (grid_metadata%yvals)) DEALLOCATE (grid_metadata%yvals)

    CALL streamClose(file_metadata%streamID)
    file_metadata = empty
  END SUBROUTINE close_output


  !> Writes data field to output data file (using the CDI).
  !
  !  @note The output file is assumed to be already open.
  !
  SUBROUTINE store_field2D(file_metadata, varID, ilevel, field)
    TYPE (t_file_metadata),  INTENT(IN) :: file_metadata
    INTEGER,                 INTENT(IN) :: varID, ilevel
    REAL(wp),                INTENT(IN) :: field(:,:)

    ! get levelID:
    CALL streamWriteVarSlice(file_metadata%streamID, varID, ilevel-1, field, 0)
  END SUBROUTINE store_field2D

END MODULE mo_remap_output
