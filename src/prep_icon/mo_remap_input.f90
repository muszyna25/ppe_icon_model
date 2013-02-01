!> Defines, loads input data.
!
!  @author F. Prill, DWD
!
MODULE mo_remap_input

  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: finish, message
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_communication,      ONLY: blk_no
  USE mo_timer,              ONLY: tic, toc
  USE mo_cf_convention,      ONLY: t_cf_var
  USE mo_grib2,              ONLY: t_grib2_var
  USE mo_io_units,           ONLY: nnml
  USE mo_namelist,           ONLY: POSITIONED, position_nml, open_nml, close_nml
  USE mo_util_sort,          ONLY: quicksort
  USE mo_util_string,        ONLY: tolower
  USE mo_mpi,                ONLY: get_my_mpi_work_id, p_comm_work, p_int, p_bcast
  USE mo_remap_config,       ONLY: dbg_level, MAX_NAME_LENGTH, MAX_NZAXIS, MAX_INPUT_FIELDS
  USE mo_remap_shared,       ONLY: t_grid, GRID_TYPE_REGULAR
  USE mo_remap_sync,         ONLY: t_gather, scatter_field2D
  USE mo_remap_io,           ONLY: t_file_metadata, get_varID

  IMPLICIT NONE
  INCLUDE 'cdi.inc'

  PRIVATE
  ! subroutines:
  PUBLIC :: read_input_namelist
  PUBLIC :: load_metadata_input
  PUBLIC :: close_input
  PUBLIC :: input_import_data
  PUBLIC :: generate_missval_mask
  ! types and variables:
  PUBLIC :: n_input_fields, input_field
  PUBLIC :: global_metadata
  PUBLIC :: n_zaxis, zaxis_metadata
  PUBLIC :: field_id_u, field_id_v
  PUBLIC :: t_field_metadata
  PUBLIC :: t_global_metadata
  PUBLIC :: t_zaxis_metadata
  PUBLIC :: t_field_adjustment
  PUBLIC :: CONST_UNINITIALIZED, CONST_FALSE, CONST_TRUE
  PUBLIC :: INTP_CONS, INTP_RBF, INTP_NONE
  PUBLIC :: CELL_GRID, EDGE_GRID

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_input')

  ! some constants (for better readability):
  INTEGER, PARAMETER :: CONST_UNINITIALIZED = -1
  INTEGER, PARAMETER :: CONST_FALSE         =  0
  INTEGER, PARAMETER :: CONST_TRUE          =  1

  ! interpolation method:
  INTEGER, PARAMETER :: INTP_NONE           =  1 !< no interpolation
  INTEGER, PARAMETER :: INTP_CONS           =  2 !< conservative interpolation weights
  INTEGER, PARAMETER :: INTP_RBF            =  3 !< RBF interpolation weights

  ! important for unstructured grid: grid type CELL/EDGE
  INTEGER, PARAMETER :: CELL_GRID           =  1 !< unstructured cell grid or structured grid
  INTEGER, PARAMETER :: EDGE_GRID           =  2 !< unstructured cell grid or structured grid
  
  CHARACTER(LEN=13), PARAMETER :: INTP_NAME_STR(3) = &
    & (/ "none/indirect", "conservative ", "rbf          " /)


  !> Data structure containing global meta data
  TYPE t_global_metadata
    ! meta-data from input file:
    INTEGER                       :: cdiInstID             !< output generating institute
    INTEGER                       :: calendar              !< calendar ID
    INTEGER                       :: tunit                 !< time axis unit
    INTEGER                       :: timetype              !< time axis type
    INTEGER                       :: rdate, rtime          !< time axis reference data,time
    INTEGER                       :: vdate, vtime          !< time axis verification data,time
    ! internal IDs for CDI
    INTEGER                       :: gridID
  END TYPE t_global_metadata

  !> Meta-data for field adjustments (e.g. hydrostatic correction)
  !  applied prior to horizontal interpolation.
  TYPE t_field_adjustment
    LOGICAL                         :: lhydrostatic_correction
    CHARACTER (LEN=MAX_NAME_LENGTH) :: var_temp                  !< field name: "temperature"
    INTEGER                         :: code_temp
    CHARACTER (LEN=MAX_NAME_LENGTH) :: var_geosp                 !< field name: "surface geopotential"
    INTEGER                         :: code_geosp
    CHARACTER (LEN=MAX_NAME_LENGTH) :: var_qv                    !< field name: "specific humidity
    INTEGER                         :: code_qv
    REAL(wp)                        :: hpbl1                     !< height above ground of surface inversion top
    REAL(wp)                        :: hpbl2                     !< top of layer used to estimate the vertical
                                                                 !  temperature gradient above the inversion
  END TYPE t_field_adjustment

  !> Data structure containing meta data for a single input field
  !
  !  @note Not all metadata is available on all working PEs in
  !        distributed computations.
  TYPE t_field_metadata
    ! user-defined (by namelist):
    CHARACTER (len=MAX_NAME_LENGTH)  :: inputname, outputname !< variable name
    INTEGER                          :: code
    CHARACTER (len=MAX_NAME_LENGTH)  :: type_of_layer         !< level type
    ! important for unstructured grid: grid type CELL/EDGE
    INTEGER                          :: grid_type             !< CELL / EDGE
    ! meta-data (from input file):
    TYPE (t_cf_var)                  :: cf                    !< CF convention information
    TYPE(t_grib2_var)                :: grib2                 !< GRIB2 related information
    INTEGER                          :: steptype
    INTEGER                          :: nlev                  !< no. of vertical levels
    ! meta-data on missing values
    INTEGER                          :: has_missvals          !< Flag: miss vals present? (-1: uninitialized,0:no,1:yes)
    REAL(wp)                         :: missval               !< missing value
    ! meta-data for field adjustment
    TYPE (t_field_adjustment)        :: fa
    ! internal IDs for CDI (IDs correspond to input data file)
    INTEGER                          :: varID, zaxisID
    ! meta-data for interpolation method
    INTEGER                          :: intp_method           !< INTP_CONS / INTP_RBF / INTP_NONE
  END TYPE t_field_metadata

  !> Data structure containing meta data for a single z-axis object
  TYPE t_zaxis_metadata
    ! meta-data (from input file):
    INTEGER                       :: zaxis_type            !< ZAXIS_SURFACE, ...
    INTEGER                       :: zaxis_size            !< no. of levels
    REAL(wp), ALLOCATABLE         :: levels(:)             !< vertical levels
    INTEGER                       :: vct_size
    REAL(wp), ALLOCATABLE         :: vct(:)
    ! internal IDs for CDI
    INTEGER                       :: zaxisID
  END TYPE t_zaxis_metadata

  TYPE(t_global_metadata)         :: global_metadata       !< global meta-data
  INTEGER                         :: n_input_fields        !< no. of INPUT fields
  TYPE (t_field_metadata), ALLOCATABLE :: input_field(:)   !< list of INPUT fields
  INTEGER                         :: n_zaxis               !< no. of z-axis objects
  TYPE(t_zaxis_metadata)          :: zaxis_metadata(MAX_NZAXIS) !< z-axis meta-data
  INTEGER                         :: field_id_u, field_id_v !< field indices for "u", "v" component (if available)

  ! input variables for the *current* field namelist (there's more than one!)
  CHARACTER (len=MAX_NAME_LENGTH) :: inputname, outputname
  INTEGER                         :: code
  CHARACTER (len=MAX_NAME_LENGTH) :: type_of_layer

  LOGICAL                         :: lhydrostatic_correction
  CHARACTER (LEN=MAX_NAME_LENGTH) :: var_temp                  !< field name: "temperature"
  INTEGER                         :: code_temp
  CHARACTER (LEN=MAX_NAME_LENGTH) :: var_geosp                 !< field name: "surface geopotential"
  INTEGER                         :: code_geosp
  CHARACTER (LEN=MAX_NAME_LENGTH) :: var_qv                    !< field name: "specific humidity
  INTEGER                         :: code_qv
  REAL(wp)                        :: hpbl1                     !< height above ground of surface inversion top
  REAL(wp)                        :: hpbl2                     !< top of layer used to estimate the vertical
                                                               !  temperature gradient above the inversion


  ! namelist definition: namelist for a single field
  NAMELIST/input_field_nml/ inputname, outputname, code, type_of_layer,   &
    &                       lhydrostatic_correction, var_temp, code_temp, &
    &                       var_geosp, code_geosp, var_qv, code_qv,       &
    &                       hpbl1, hpbl2


CONTAINS

  !> Opens input namelist file and loads meta data.
  !
  SUBROUTINE read_input_namelist(input_cfg_filename, rank0, lcompute_vn, lcompute_vt)
    CHARACTER(LEN=*), INTENT(IN) :: input_cfg_filename !< field config namelist file
    INTEGER,          INTENT(IN) :: rank0
    LOGICAL,          INTENT(IN) :: lcompute_vn   !< Flag: .TRUE., if normal wind shall be computed (instead of u,v)
    LOGICAL,          INTENT(IN) :: lcompute_vt   !< Flag: .TRUE., if tangential wind shall be computed
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::read_input_namelist')
    INTEGER :: istat, ierrstat
    LOGICAL :: lrewind

    n_input_fields = 0
    field_id_u = CONST_UNINITIALIZED
    field_id_v = CONST_UNINITIALIZED

    ALLOCATE(input_field(MAX_INPUT_FIELDS), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! -- second read field namelists (there's more than one)
    IF (dbg_level >= 5) &
      &  WRITE (0,*) "# read field namelists from '"//TRIM(input_cfg_filename)//"'"
    lrewind = .TRUE.
    CALL open_nml(TRIM(input_cfg_filename))
    CFG_LOOP : DO
      ! move to the next field namelist in the config file:
      CALL position_nml ("input_field_nml", lrewind=lrewind, status=istat)
      IF(istat /= POSITIONED) THEN
        CALL close_nml
        EXIT CFG_LOOP
      ENDIF
      lrewind = .FALSE.

      ! default settings
      inputname     = " "
      outputname    = " "
      code          =   0
      type_of_layer = " "

      lhydrostatic_correction = .FALSE.
      var_temp      = "T"
      code_temp     = 130
      var_geosp     = "FI"
      code_geosp    = 129
      var_qv        = "QV"
      code_qv       = 133
      hpbl1         = 500._wp
      hpbl2         = 1000._wp

      ! read user's (new) specifications
      READ (nnml, input_field_nml)
      IF(istat > 0)  CALL finish(routine, "Internal error!")

      n_input_fields = n_input_fields + 1
      input_field(n_input_fields)%inputname       = TRIM(inputname)
      input_field(n_input_fields)%outputname      = TRIM(outputname)
      input_field(n_input_fields)%code            = code
      input_field(n_input_fields)%type_of_layer   = TRIM(type_of_layer)
      input_field(n_input_fields)%grid_type       = CELL_GRID

      ! default interpolation method: INTP_CONS, only "u", "v" are
      ! treated with RBF interpolation weights:
      IF ((lcompute_vt) .OR. (lcompute_vn)) THEN
        ! NOTE: For development reasons, we do not write "vn"/"vt"
        ! *instead* of "u", "v", but in addition to these fields.
        IF (TRIM(tolower(inputname)) == "u") THEN
          ! input_field(n_input_fields)%intp_method = INTP_NONE  
          input_field(n_input_fields)%intp_method = INTP_CONS        
          field_id_u = n_input_fields
        ELSE IF (TRIM(tolower(inputname)) == "v") THEN
          ! input_field(n_input_fields)%intp_method = INTP_NONE 
          input_field(n_input_fields)%intp_method = INTP_CONS        
          field_id_v = n_input_fields
        ELSE
          input_field(n_input_fields)%intp_method = INTP_CONS        
        END IF
      ELSE
        input_field(n_input_fields)%intp_method = INTP_CONS        
      END IF

      input_field(n_input_fields)%fa%lhydrostatic_correction  = lhydrostatic_correction
      input_field(n_input_fields)%fa%var_temp   = var_temp
      input_field(n_input_fields)%fa%code_temp  = code_temp
      input_field(n_input_fields)%fa%var_geosp  = var_geosp
      input_field(n_input_fields)%fa%code_geosp = code_geosp
      input_field(n_input_fields)%fa%var_qv     = var_qv
      input_field(n_input_fields)%fa%code_qv    = code_qv
      input_field(n_input_fields)%fa%hpbl1      = hpbl1
      input_field(n_input_fields)%fa%hpbl2      = hpbl2

      ! not yet initialized:
      input_field(n_input_fields)%steptype      = -1
      input_field(n_input_fields)%varID         = CDI_UNDEFID
      input_field(n_input_fields)%zaxisID       = CDI_UNDEFID

      IF (get_my_mpi_work_id() == rank0) &
        &  CALL input_print_metadata(input_field(n_input_fields))
    END DO CFG_LOOP

    ! ------------------
    ! consistency checks
    ! ------------------

    ! error message if 0 input fields were configured:
    IF (n_input_fields == 0) &
      &  CALL finish(routine, "No input fields found in configuration namelist!")

    IF (ANY((input_field(:)%fa%lhydrostatic_correction) .AND.   &
      &     (input_field(:)%intp_method  == INTP_RBF))) THEN
      CALL finish(routine, "No hydrostatic correction with RBF interpolated fields!")
    END IF

  END SUBROUTINE read_input_namelist


  !> Formatted print-out of input field configuration.
  !
  SUBROUTINE input_print_metadata(field)
    TYPE (t_field_metadata), INTENT(IN) :: field

    WRITE (0,*)        "# input field '"//TRIM(field%inputname)//"' -> '"//TRIM(field%outputname)//"'"
    WRITE (0,'(a,i3)') " #     code          : ", field%code
    WRITE (0,'(a,a)')  " #     type_of_layer : ", TRIM(field%type_of_layer)
    WRITE (0,'(a,a)')  " #     interpolation : ", TRIM(INTP_NAME_STR(field%intp_method))
  END SUBROUTINE input_print_metadata


  !> @return vlist variable ID for a given input field.
  !
  !  Uses cdilib for file access. Checks for consistency of meta data.
  !
  !  @todo Consistency checks of variable metadata still missing!
  !
  FUNCTION get_field_varID(vlistID, gribedition, field_info) RESULT(result_varID)
    INTEGER                             :: result_varID
    INTEGER,                 INTENT(IN) :: vlistID             !< link to GRIB file vlist
    INTEGER,                 INTENT(IN) :: gribedition         !< grib edition 1/2
    TYPE (t_field_metadata), INTENT(IN) :: field_info          !< field meta-data

    IF (gribedition == 1) THEN
      result_varID = get_varID(vlistID, code=field_info%code)
    ELSE
      result_varID = get_varID(vlistID, name=field_info%inputname)
    END IF
  END FUNCTION get_field_varID


  !> loads variable meta-data
  !
  SUBROUTINE load_metadata_input(file, rank0)
    TYPE (t_file_metadata), INTENT(INOUT) :: file
    INTEGER,                INTENT(IN)    :: rank0
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//':load_metadata_input')
    INTEGER :: i, param, taxisID, ierrstat, zaxisID, vlistID, ntsteps
    INTEGER, ALLOCATABLE :: zaxisIDlist(:)

    ! load global meta-data
    IF (get_my_mpi_work_id() == rank0) THEN
      vlistID = file%vlistID
      global_metadata%cdiInstID = vlistInqInstitut(vlistID)

      taxisID = vlistInqTaxis(vlistID)
      global_metadata%calendar  = taxisInqCalendar(taxisID)
      global_metadata%tunit     = taxisInqTunit(taxisID)
      global_metadata%timetype  = taxisInqType(taxisID)
      global_metadata%rdate     = taxisInqRdate(taxisID)
      global_metadata%rtime     = taxisInqRtime(taxisID)
      global_metadata%vdate     = taxisInqVdate(taxisID)
      global_metadata%vtime     = taxisInqVtime(taxisID)
      ntsteps = vlistNtsteps(vlistID)
      IF (ntsteps > 1) CALL finish(routine, "File contains more than one time step!")
    END IF

    input_field(:)%nlev = 0
    ! Loop over all fields from input config file, read file only on PE rank0:
    DO i=1,n_input_fields
      IF (get_my_mpi_work_id() == rank0) THEN
        input_field(i)%varID = get_field_varID(vlistID, file%gribedition, input_field(i))
      END IF
      ! distributed computation: broadcast to other worker PEs:
      CALL p_bcast(input_field(i)%varID, rank0, p_comm_work)
      IF (input_field(i)%varID /= -1) THEN
        IF (get_my_mpi_work_id() == rank0) THEN
          ! get the corresponding grid ID
          input_field(i)%zaxisID = vlistInqVarZaxis(vlistID, input_field(i)%varID)
          ! get the field dimensions:
          input_field(i)%nlev     = zaxisInqSize(input_field(i)%zaxisID)
          IF ((input_field(i)%type_of_layer == "surface") .AND.  &
            & (input_field(i)%nlev /= 1)) THEN
            CALL finish(routine, TRIM(input_field(i)%inputname)//": Inconsistent input data (level no.)")
          END IF
          input_field(i)%steptype = vlistInqVarTsteptype(vlistID, input_field(i)%varID)
          ! get CF information on variable:
          CALL vlistInqVarLongname(vlistID, input_field(i)%varID, input_field(i)%cf%long_name    )
          CALL vlistInqVarStdName( vlistID, input_field(i)%varID, input_field(i)%cf%standard_name)
          CALL vlistInqVarUnits(   vlistID, input_field(i)%varID, input_field(i)%cf%units        )
          input_field(i)%cf%datatype = vlistInqVarDatatype(vlistID, input_field(i)%varID)
          ! get GRIB2 information on variable:
          param = vlistInqVarParam(vlistID, input_field(i)%varID)
          CALL cdiDecodeParam(param, input_field(i)%grib2%number,    &
            &                        input_field(i)%grib2%category,  &
            &                        input_field(i)%grib2%discipline)
          input_field(i)%grib2%bits = vlistInqVarDatatype(vlistID, input_field(i)%varID)
          ! get missing value information on variable:
          input_field(i)%missval = vlistInqVarMissval(vlistID, input_field(i)%varID)
        END IF ! rank0
        ! distributed computation: broadcast to other worker PEs:
        CALL p_bcast(input_field(i)%nlev,rank0,p_comm_work)
        CALL p_bcast(input_field(i)%missval,rank0,p_comm_work)
        input_field(i)%has_missvals = CONST_UNINITIALIZED ! unitialized
      END IF ! (varID /= -1)
    END DO

    ! create a unique list of z-axes which are in use:
    IF (dbg_level >= 10)  WRITE (0,*) "# read z-axis meta-data"
    IF (get_my_mpi_work_id() == rank0) THEN
      n_zaxis = 0
      ALLOCATE(zaxisIDlist(n_input_fields), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      zaxisIDlist(:) = input_field(1:n_input_fields)%zaxisID
      CALL quicksort(zaxisIDlist)
      DO i=1,n_input_fields
        IF (i>1) THEN
          IF (zaxisIDlist(i) == zaxisIDlist(i-1)) CYCLE
        END IF
        n_zaxis = n_zaxis + 1
        zaxisID = zaxisIDlist(i)
        zaxis_metadata(n_zaxis)%zaxisID     = zaxisID
        zaxis_metadata(n_zaxis)%zaxis_type  = zaxisInqType(zaxisID)
        zaxis_metadata(n_zaxis)%zaxis_size  = zaxisInqSize(zaxisID)
        IF (zaxis_metadata(n_zaxis)%zaxis_size > 0) THEN
          ALLOCATE(zaxis_metadata(n_zaxis)%levels(zaxis_metadata(n_zaxis)%zaxis_size), &
            &      STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
          CALL zaxisInqLevels(zaxisID, zaxis_metadata(n_zaxis)%levels)
        END IF

        zaxis_metadata(n_zaxis)%vct_size = zaxisInqVctSize(zaxisID)
        IF (zaxis_metadata(n_zaxis)%vct_size > 0) THEN
          ALLOCATE(zaxis_metadata(n_zaxis)%vct(zaxis_metadata(n_zaxis)%vct_size), &
            &      STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
          CALL zaxisInqVct(zaxisID, zaxis_metadata(n_zaxis)%vct);
        END IF
      END DO
      DEALLOCATE(zaxisIDlist, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
  END SUBROUTINE load_metadata_input


  !> Utility routine: Generate LOGICAL-array mask from a given field
  !  variable, corresponding to the missing value pattern.  The value
  !  .TRUE. means: missing value exists.
  !
  SUBROUTINE generate_missval_mask(in_field, missval, out_mask2D, nblks, npromz)
    REAL(wp),           INTENT(IN)    :: in_field(:,:)
    REAL(wp),           INTENT(IN)    :: missval
    LOGICAL,            INTENT(INOUT) :: out_mask2D(:,:)
    INTEGER,            INTENT(IN)    :: nblks, npromz
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = &
      &  TRIM(TRIM(modname)//':generate_missval_mask')
    INTEGER :: i,j, size1, size2, nmiss, i_endidx

    IF (dbg_level >= 1) WRITE (0,*) "# generate mask for missing values."

    size1 = SIZE(in_field,1)
    size2 = SIZE(in_field,2)
    IF ((size1 /= SIZE(out_mask2D,1)) .OR. &
      & (size2 /= SIZE(out_mask2D,2)) .OR. &
      & (size1 /= nproma)             .OR. &
      & (size2 > nblks)) THEN
      CALL finish(routine, "Invalid array dimension!")
    END IF

!$OMP PARALLEL 
!$OMP DO PRIVATE(i,j,i_endidx)
    DO j=1,nblks
      i_endidx = nproma
      IF (j==nblks) i_endidx=npromz
      DO i=1,i_endidx
        out_mask2D(i,j) = (in_field(i,j) == missval) ! Attention: Comparison of REALs
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    IF (dbg_level >= 1) THEN
      nmiss = COUNT(out_mask2D)
      WRITE (0,*) "# masking ", nmiss, " values."
    END IF

  END SUBROUTINE generate_missval_mask


  !> Main subroutine: Reads data fields.
  !
  !  Requires a pre-allocated variable list ("varlist") as output.
  !
  ! @todo Check if grid ID corresponds to the grid for which the
  !       interpolation weights are computed.
  !
  SUBROUTINE input_import_data(file_metadata, ivar, gather_c, &
    &                          dst_field, src_grid,           &
    &                          opt_time_comm, opt_time_read   )

    TYPE (t_file_metadata), INTENT(IN) :: file_metadata
    INTEGER,           INTENT(IN)    :: ivar
    TYPE(t_gather),    INTENT(IN)    :: gather_c         !< communication pattern
    REAL(wp),          INTENT(INOUT) :: dst_field(:,:,:)
    TYPE (t_grid),     INTENT(IN)    :: src_grid
    REAL, INTENT(INOUT), OPTIONAL    :: opt_time_comm, opt_time_read
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//':input_import_data'
    INTEGER :: streamID, ierrstat, nmiss, nsize,  &
      &        shape2D_glb(2), shape2D_loc(2), ilev, maxblk, has_missvals
    REAL(wp), ALLOCATABLE :: rfield1D(:), rfield2D(:,:), rfield2D_loc(:,:)
    REAL    :: time_comm, time_read

    has_missvals = CONST_UNINITIALIZED

    streamID = file_metadata%streamID
    IF (dbg_level >= 2)  CALL message(routine, "Start")
    ! allocate global input field only on PE rank0
    IF (get_my_mpi_work_id() == gather_c%rank0) THEN
      shape2D_glb = (/ nproma,blk_no(src_grid%p_patch%n_patch_cells_g) /)
    ELSE
      shape2D_glb = (/ 0,0 /)
    END IF
    shape2D_loc = (/ nproma,src_grid%p_patch%nblks_c /)
    ALLOCATE(rfield2D(    shape2D_glb(1), shape2D_glb(2)), &
      &      rfield2D_loc(shape2D_loc(1), shape2D_loc(2)), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! Loop over all levels, reading each variable in 2D slices.
    IF (input_field(ivar)%varID /= -1) THEN
      IF (get_my_mpi_work_id() == gather_c%rank0) THEN
        ! get the 1D input field size
        nsize = gridInqSize(file_metadata%c_gridID)
        IF (.NOT. ALLOCATED(rfield1D)) THEN
          ! allocate a temporary field to hold the GRIB data
          ALLOCATE(rfield1D(nsize), STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        ELSE
          ! paranoia, again: check if the field size has changed
          IF (nsize /= SIZE(rfield1D)) CALL finish(routine, "Internal error!")
        END IF
      END IF ! rank0

      ! read field data:
      ! ----------------
      DO ilev=1,input_field(ivar)%nlev
        IF (get_my_mpi_work_id() == gather_c%rank0) THEN
          ! read record as 1D field
          IF (dbg_level >= 10) WRITE (0,*) "# ",TRIM(input_field(ivar)%inputname),": read record ", ilev
          CALL tic(time_read)  ! performance measurement: start
          CALL streamReadVarSlice(streamID, input_field(ivar)%varID, ilev-1, rfield1D, nmiss)
          IF (PRESENT(opt_time_read)) opt_time_read = opt_time_read + toc(time_read)
          ! reshape record into (nproma, nblks) field
          rfield2D(:,:) = RESHAPE(rfield1D(:), shape2D_glb, (/ 0._wp /))

          ! check if missing values are present, set flag "has_missvals"
          !
          ! note: probably we could use a cdi internal function for this, but
          !       which one?
          IF (has_missvals == CONST_UNINITIALIZED) THEN
            IF (dbg_level >= 2) &
              &  WRITE (0,*) "# testing for missing values (missval = ", input_field(ivar)%missval, ")"
            has_missvals = CONST_FALSE
            ! Attention: Comparison of REALs:
            IF (ANY(rfield2D(:,:) ==  input_field(ivar)%missval)) THEN
              has_missvals = CONST_TRUE
            END IF
          END IF
        END IF ! if rank0
        IF (input_field(ivar)%has_missvals == CONST_UNINITIALIZED) THEN
          IF (get_my_mpi_work_id() == gather_c%rank0) input_field(ivar)%has_missvals = has_missvals
          CALL p_bcast(input_field(ivar)%has_missvals,gather_c%rank0,p_comm_work)
        END IF

        ! for distributed computation: scatter interpolation input to work PEs
        CALL tic(time_comm)  ! performance measurement: start
        CALL scatter_field2D(gather_c, rfield2D, rfield2D_loc)
        IF (PRESENT(opt_time_comm)) opt_time_comm = opt_time_comm + toc(time_comm)

        !-----------------------
        ! Copy level to 3D field
        ! (horizontal interpolation moved to caller!)
        !-----------------------
        maxblk = min (size (dst_field, dim=3), size (rfield2D_loc, dim=2))
        dst_field(:,ilev,1:maxblk)  = rfield2D_loc(:,1:maxblk)
        dst_field(:,ilev,maxblk+1:) = 0._wp
      END DO ! ilev
    END IF ! (varID /= -1)

    ! clean up
    DEALLOCATE(rfield2D, rfield2D_loc, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    IF (ALLOCATED(rfield1D)) THEN
      DEALLOCATE(rfield1D, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
    IF (dbg_level >= 2)  CALL message(routine, "Done")
  END SUBROUTINE input_import_data


  !> Close regular grid file
  !
  SUBROUTINE close_input()
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//':close_input')
    INTEGER :: ierrstat, i

    DO i=1,n_zaxis
      IF (ALLOCATED(zaxis_metadata(i)%vct)) THEN
        DEALLOCATE(zaxis_metadata(i)%vct,    STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed 1!")
      END IF
      IF (ALLOCATED(zaxis_metadata(i)%levels)) THEN
        DEALLOCATE(zaxis_metadata(i)%levels, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed 2!")
      END IF
    END DO
    DEALLOCATE(input_field, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed 3!")
  END SUBROUTINE close_input

END MODULE mo_remap_input
