MODULE mo_remap_io

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: nnml
  USE mo_mpi,                ONLY: get_my_mpi_work_id
  USE mo_namelist,           ONLY: POSITIONED, position_nml, open_nml, close_nml
  USE mo_exception,          ONLY: finish
  USE mo_util_string,        ONLY: tolower
  USE mo_util_netcdf,        ONLY: nf
  USE mo_remap_config,       ONLY: dbg_level, MAX_NAME_LENGTH, MAX_NSTENCIL
  USE mo_remap_sync,         ONLY: t_gather_c, scatter_field2D_c
  USE mo_remap_shared,       ONLY: GRID_TYPE_REGULAR, GRID_TYPE_ICON    
  IMPLICIT NONE

  INCLUDE 'cdi.inc'
  INCLUDE 'netcdf.inc'

  PRIVATE
  PUBLIC :: read_remap_namelist
  PUBLIC :: open_file, close_file
  PUBLIC :: read_field2D, read_field3D
  PUBLIC :: get_varID
  PUBLIC :: in_filename,      out_filename,      &
    &       in_grid_filename, out_grid_filename, &
    &       in_type,          out_type
  PUBLIC :: t_file_metadata
  PUBLIC :: l_have3dbuffer
  PUBLIC :: s_maxsize
  PUBLIC :: in_file_gribedition

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_io')  

  CHARACTER (len=MAX_NAME_LENGTH)    :: in_filename         !< input file name.
  CHARACTER (len=MAX_NAME_LENGTH)    :: out_filename        !< output file name.
  CHARACTER (len=MAX_NAME_LENGTH)    :: in_grid_filename    !< input grid file name.
  CHARACTER (len=MAX_NAME_LENGTH)    :: out_grid_filename   !< output grid file name.

  !> work-around: GRIB edition number of input data file (as long as we cannot
  !  determine this via CDI)
  INTEGER                            :: in_file_gribedition

  INTEGER :: in_type, out_type !< structure of input/output grid

  ! Flag. True if a 3D buffer is used for output:
  LOGICAL  :: l_have3dbuffer

  ! Maximum size of sequential list for very large stencils
  INTEGER  :: s_maxsize

  ! namelist definition: main namelist
  NAMELIST/remap_nml/   in_grid_filename,  in_filename,  in_type,  &
    &                   out_grid_filename, out_filename, out_type, &
    &                   l_have3dbuffer, s_maxsize,                 &
    &                   in_file_gribedition

  ! Derived type containing IDs for opened grid/data files.
  !
  ! This is necessary because structured grid files are opened using
  ! the CDI, unstructured grid files are opened using the NetCDF
  ! library.
  !
  TYPE t_file_metadata
    INTEGER :: structure, streamID, vlistID, ncfileID

    ! CDI internal ID for the gaussian grid or the cells/edges/verts
    ! grids (unstructured case):
    INTEGER :: gridID, c_gridID, e_gridID, v_gridID
  END TYPE t_file_metadata

CONTAINS

  !> Opens namelist file, reads interpolation config.
  !
  SUBROUTINE read_remap_namelist(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename !< main namelist file
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::read_remap_namelist')
    INTEGER :: istat

    IF (dbg_level >= 5) WRITE (0,*) "# read io namelist"
    ! default settings
    in_grid_filename  = ""
    in_filename       = ""
    in_type           =  1
    out_grid_filename = ""
    out_filename      = ""
    out_type          =  2
    l_have3dbuffer    = .FALSE.
    s_maxsize         = 500000
    in_file_gribedition = 2

    ! read user's (new) specifications
    CALL open_nml(TRIM(filename))
    CALL position_nml ('remap_nml', status=istat)
    IF (istat == POSITIONED) THEN
      READ (nnml, remap_nml)
    ELSE
      CALL finish(routine, "Internal error!")
    END IF
    CALL close_nml
    ! status output
    IF (dbg_level >= 2)  WRITE (0,*) "# stencil size: ", MAX_NSTENCIL, "/", s_maxsize
  END SUBROUTINE read_remap_namelist


  !> Open grid/data file for input
  !
  !  @note structured grid files are opened using the CDI,
  !        unstructured grid files are also opened using the NetCDF
  !        library.
  !
  SUBROUTINE open_file(filename, istructure, file_metadata, rank0)
    CHARACTER(LEN=*),       INTENT(IN)  :: filename
    INTEGER,                INTENT(IN)  :: istructure
    TYPE (t_file_metadata), INTENT(OUT) :: file_metadata
    INTEGER,                INTENT(IN)  :: rank0
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM('mo_grid_interface::open_file')
    INTEGER                     :: ngrids, i, gridID
    CHARACTER(len=MAX_NAME_LENGTH) :: zname

    file_metadata%structure = istructure

    IF (get_my_mpi_work_id() == rank0) THEN
      SELECT CASE(file_metadata%structure) 
      CASE (GRID_TYPE_ICON)
        IF (dbg_level >= 2)  WRITE (0,*) "# open file '",TRIM(filename),"' (unstructured NetCDF)"
        CALL nf(nf_open(TRIM(filename), NF_NOWRITE, file_metadata%ncfileID))
        file_metadata%structure = GRID_TYPE_ICON
        file_metadata%streamID  = streamOpenRead(TRIM(filename))
        file_metadata%vlistID   = streamInqVlist(file_metadata%streamID)
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
        IF (dbg_level >= 2)  WRITE (0,*) "# open file '",TRIM(filename),"' (structured)"
        file_metadata%structure = GRID_TYPE_REGULAR
        file_metadata%streamID  = streamOpenRead(TRIM(filename))
        file_metadata%vlistID   = streamInqVlist(file_metadata%streamID)
        ! set file gridID to regular grid:
        file_metadata%gridID    = vlistGrid(file_metadata%vlistID, 0)
      CASE DEFAULT
        CALL finish(routine, "Unknown grid type")
      END SELECT
    END IF
  END SUBROUTINE open_file


  !> Open grid/data file for input
  !
  SUBROUTINE close_file(file_metadata)
    TYPE (t_file_metadata), INTENT(INOUT) :: file_metadata
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM('mo_grid_interface::close_file')
    
    SELECT CASE(file_metadata%structure) 
    CASE (GRID_TYPE_ICON)
      CALL nf(nf_close(file_metadata%ncfileID))
      CALL streamClose(file_metadata%streamID)
    CASE (GRID_TYPE_REGULAR)
      CALL streamClose(file_metadata%streamID)
    CASE DEFAULT
      CALL finish(routine, "Unknown grid type")
    END SELECT

    file_metadata%structure = -1
    file_metadata%streamID  = -1
    file_metadata%vlistID   = -1
    file_metadata%ncfileID  = -1
  END SUBROUTINE close_file


  !> Utility routine: Read and communicate 2D field.
  !
  SUBROUTINE read_field3D(file_metadata, var_name, var_code, nlev, gather_c, rfield1D, &
    &                     rfield2D, rfield3D_loc, opt_zaxisID)
    TYPE (t_file_metadata), INTENT(IN)    :: file_metadata      !< input file meta-data
    CHARACTER(LEN=*),       INTENT(IN)    :: var_name           !< variable name string
    INTEGER,                INTENT(IN)    :: var_code           !< GRIB variable code
    INTEGER,                INTENT(IN)    :: nlev               !< no. of levels to read
    TYPE(t_gather_c),       INTENT(IN)    :: gather_c           !< communication pattern
    REAL(wp),               INTENT(INOUT) :: rfield1D(:),   &   !< global 1D field (on rank0)
      &                                      rfield2D(:,:), &   !< global 2D field (on rank0)
      &                                      rfield3D_loc(:,:,:)!< local 2D field
    INTEGER, INTENT(OUT), OPTIONAL   :: opt_zaxisID
    ! local variables:
    INTEGER :: streamID, vlistID, varID, ilev, nmiss
    INTEGER :: shape2d(2)

    shape2d = (/ UBOUND(rfield2D,1), UBOUND(rfield2D,2) /)

    IF (get_my_mpi_work_id() == gather_c%rank0) THEN
      streamID = file_metadata%streamID
      vlistID  = file_metadata%vlistID
      IF (in_file_gribedition == 1) THEN
        varID    = get_varID(vlistID, code=var_code)
      ELSE
        varID    = get_varID(vlistID, name=TRIM(var_name))
      END IF
      IF (PRESENT(opt_zaxisID)) THEN
        opt_zaxisID = vlistInqVarZaxis(vlistID, varID)
      END IF
    END IF

    DO ilev=1,nlev
      IF (get_my_mpi_work_id() == gather_c%rank0) THEN
        ! read record as 1D field
        CALL streamReadVarSlice(streamID, varID, ilev-1, rfield1D, nmiss)
        ! reshape record into (nproma, nblks) field
        rfield2D(:,:) = RESHAPE(rfield1D(:), shape2d, (/ 0._wp /))
      END IF ! if rank0
      ! for distributed computation: scatter interpolation input to work PEs
      CALL scatter_field2D_c(gather_c, rfield2D, rfield3D_loc(:,ilev,:))
    END DO ! ilev
  END SUBROUTINE read_field3D


  !> Utility routine: Read and communicate 2D field.
  !
  SUBROUTINE read_field2D(file_metadata, var_name, var_code, gather_c, rfield1D, &
    &                     rfield2D, rfield2D_loc, opt_ilev)
    TYPE (t_file_metadata), INTENT(IN)    :: file_metadata      !< input file meta-data
    CHARACTER(LEN=*),       INTENT(IN)    :: var_name           !< variable name string
    INTEGER,                INTENT(IN)    :: var_code           !< GRIB variable code
    TYPE(t_gather_c),       INTENT(IN)    :: gather_c           !< communication pattern
    REAL(wp),               INTENT(INOUT) :: rfield1D(:),   &   !< global 1D field (on rank0)
      &                                      rfield2D(:,:), &   !< global 2D field (on rank0)
      &                                      rfield2D_loc(:,:)  !< local 2D field
    INTEGER, INTENT(IN), OPTIONAL :: opt_ilev !< (Optional) level index
    ! local variables:
    INTEGER :: streamID, vlistID, varID, ilev, nmiss
    INTEGER :: shape2d(2)

    shape2d = (/ UBOUND(rfield2D,1), UBOUND(rfield2D,2) /)

    ilev = 1
    IF (PRESENT(opt_ilev)) ilev = opt_ilev

    IF (get_my_mpi_work_id() == gather_c%rank0) THEN
      streamID = file_metadata%streamID
      vlistID  = file_metadata%vlistID
      IF (in_file_gribedition == 1) THEN
        varID    = get_varID(vlistID, code=var_code)
      ELSE
        varID    = get_varID(vlistID, name=TRIM(var_name))
      END IF

      ! read record as 1D field
      CALL streamReadVarSlice(streamID, varID, ilev-1, rfield1D, nmiss)
      ! reshape record into (nproma, nblks) field
      rfield2D(:,:) = RESHAPE(rfield1D(:), shape2d, (/ 0._wp /))
    END IF ! if rank0
    ! for distributed computation: scatter interpolation input to work PEs
    CALL scatter_field2D_c(gather_c, rfield2D, rfield2D_loc)
  END SUBROUTINE read_field2D


  !> @return vlist variable ID for a given variable name
  !
  !  Uses cdilib for file access. 
  ! 
  FUNCTION get_varID(vlistID, name, code) RESULT(result_varID)
    INTEGER                                 :: result_varID
    INTEGER,           INTENT(IN)           :: vlistID             !< link to GRIB file vlist
    CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: name                !< variable name
    INTEGER,           INTENT(IN), OPTIONAL :: code                !< GRIB1 variable code
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::get_varID')
    CHARACTER(len=MAX_NAME_LENGTH) :: zname
    LOGICAL                        :: l_found
    INTEGER                        :: nvars, varID, icode

    IF (.NOT. (PRESENT(name) .OR. PRESENT(code))) &
      & CALL finish(routine, "Internal error!")

    zname = ""
    icode =  0

    result_varID = -1
    ! total number of available fields:
    nvars = vlistNvars(vlistID)
    ! loop over vlist, find the corresponding varID
    l_found = .FALSE.
    LOOP : DO varID=0,(nvars-1)
      IF (PRESENT(name)) CALL vlistInqVarName(vlistID, varID, zname)
      IF (PRESENT(code)) icode = vlistInqVarCode(vlistID, varID)
      IF (dbg_level >= 10) THEN
        WRITE (0,*) "# scanning [", tolower(TRIM(zname)), "], code: ", icode
      END IF

      IF (PRESENT(name)) THEN
        IF (tolower(TRIM(zname)) == tolower(TRIM(name))) THEN
          result_varID = varID
          l_found = .TRUE.
          EXIT LOOP
        END IF
      END IF
      IF (PRESENT(code)) THEN
        IF (code == icode) THEN
          result_varID = varID
          l_found = .TRUE.
          EXIT LOOP
        END IF
      END IF
    END DO LOOP
    IF (.NOT. l_found) THEN
      IF (PRESENT(name)) THEN
        CALL finish(routine, "Variable "//TRIM(name)//" not found!")
      END IF
      IF (PRESENT(code)) THEN
        WRITE (0,*) "Variable code: ", code
        CALL finish(routine, "Variable not found!")
      END IF
    END IF
  END FUNCTION get_varID

END MODULE mo_remap_io
