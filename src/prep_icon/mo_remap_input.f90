!> Defines, loads and remaps input data.
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
  USE mo_mpi,                ONLY: get_my_mpi_work_id, p_comm_work, p_int, p_bcast
  USE mo_remap_config,       ONLY: dbg_level, MAX_NAME_LENGTH
  USE mo_remap_shared,       ONLY: t_grid, GRID_TYPE_REGULAR
  USE mo_remap_intp,         ONLY: t_intp_data, interpolate_c
  USE mo_remap_sync,         ONLY: t_gather_c, scatter_field2D_c
  USE mo_remap_io,           ONLY: t_file_metadata, get_varID
  USE mo_remap_hydcorr,      ONLY: t_field_adjustment, hydrostatic_correction

  IMPLICIT NONE
  INCLUDE 'cdi.inc'

  PRIVATE
  ! subroutines:
  PUBLIC :: read_input_namelist
  PUBLIC :: load_metadata_input
  PUBLIC :: close_input
  PUBLIC :: input_import_data
  ! types and variables:
  PUBLIC :: n_input_fields, input_field
  PUBLIC :: global_metadata
  PUBLIC :: n_zaxis, zaxis_metadata
  PUBLIC :: t_field_metadata
  PUBLIC :: t_global_metadata
  PUBLIC :: t_zaxis_metadata
  PUBLIC :: MAX_INPUT_FIELDS, MAX_NZAXIS

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_input')  

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

  !> Data structure containing meta data for a single input field
  !
  !  @note Not all metadata is available on all working PEs in
  !        distributed computations.
  TYPE t_field_metadata
    ! user-defined (by namelist):
    CHARACTER (len=MAX_NAME_LENGTH)  :: inputname, outputname !< variable name
    INTEGER                          :: code, table
    CHARACTER (len=MAX_NAME_LENGTH)  :: type_of_layer         !< level type
    ! meta-data (from input file):
    TYPE (t_cf_var)                  :: cf                    !< CF convention information 
    TYPE(t_grib2_var)                :: grib2                 !< GRIB2 related information
    INTEGER                          :: steptype
    INTEGER                          :: idim(2)               !< horizontal field dimensions
    INTEGER                          :: nlev                  !< no. of vertical levels
    ! Meta-data for field adjustment
    TYPE (t_field_adjustment)        :: fa
    ! internal IDs for CDI
    INTEGER :: varID, gridID, zaxisID
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

  ! meta data ("config state")
  INTEGER, PARAMETER              :: MAX_INPUT_FIELDS = 50 !< maximum number of INPUT fields
  INTEGER, PARAMETER              :: MAX_NZAXIS     = 15   !< maximum number of z-axis objects

  TYPE(t_global_metadata)         :: global_metadata       !< global meta-data
  INTEGER                         :: n_input_fields        !< no. of INPUT fields
  TYPE (t_field_metadata), ALLOCATABLE :: input_field(:)   !< list of INPUT fields
  INTEGER                         :: n_zaxis               !< no. of z-axis objects
  TYPE(t_zaxis_metadata)          :: zaxis_metadata(MAX_NZAXIS) !< z-axis meta-data

  ! input variables for the *current* field namelist (there's more than one!)
  CHARACTER (len=MAX_NAME_LENGTH) :: inputname, outputname
  INTEGER                         :: code
  CHARACTER (len=MAX_NAME_LENGTH) :: type_of_layer

  LOGICAL                         :: lhydrostatic_correction
  CHARACTER (LEN=MAX_NAME_LENGTH) :: var_temp                  !< field name: "temperature"
  CHARACTER (LEN=MAX_NAME_LENGTH) :: var_z                     !< field name: "geopotential"
  CHARACTER (LEN=MAX_NAME_LENGTH) :: var_geosp                 !< field name: "surface geopotential"
  CHARACTER (LEN=MAX_NAME_LENGTH) :: var_qv                    !< field name: "specific humidity
  REAL(wp)                        :: hpbl1                     !< height above ground of surface inversion top
  REAL(wp)                        :: hpbl2                     !< top of layer used to estimate the vertical 

  ! namelist definition: namelist for a single field
  NAMELIST/input_field_nml/ inputname, outputname, code, type_of_layer,  &
    &                       lhydrostatic_correction, var_temp, var_z,    &
    &                       var_geosp, var_qv, hpbl1, hpbl2                    


CONTAINS

  !> Opens input namelist file and loads meta data.
  !
  SUBROUTINE read_input_namelist(input_cfg_filename, rank0)
    CHARACTER(LEN=*), INTENT(IN) :: input_cfg_filename !< field config namelist file
    INTEGER,                INTENT(IN)    :: rank0
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::read_input_namelist')
    INTEGER :: istat, ierrstat
    LOGICAL :: lrewind

    n_input_fields = 0    
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
      var_temp      = ""            
      var_z         = ""            
      var_geosp     = ""            
      var_qv        = ""            
      hpbl1         = 500
      hpbl2         = 1000

      ! read user's (new) specifications
      READ (nnml, input_field_nml)
      IF(istat > 0)  CALL finish(routine, "Internal error!")

      n_input_fields = n_input_fields + 1
      input_field(n_input_fields)%inputname       = TRIM(inputname)
      input_field(n_input_fields)%outputname      = TRIM(outputname)    
      input_field(n_input_fields)%code            = code                     
      input_field(n_input_fields)%type_of_layer   = TRIM(type_of_layer)

      input_field(n_input_fields)%fa%lhydrostatic_correction  = lhydrostatic_correction     
      input_field(n_input_fields)%fa%var_temp                 = var_temp                    
      input_field(n_input_fields)%fa%var_z                    = var_z                       
      input_field(n_input_fields)%fa%var_geosp                = var_geosp                   
      input_field(n_input_fields)%fa%var_qv                   = var_qv                      
      input_field(n_input_fields)%fa%hpbl1                    = hpbl1                       
      input_field(n_input_fields)%fa%hpbl2                    = hpbl2                       

      IF (get_my_mpi_work_id() == rank0) &
        &  CALL input_print_metadata(input_field(n_input_fields))
    END DO CFG_LOOP

    ! error message if 0 input fields were configured:
    IF (n_input_fields == 0) &
      &  CALL finish(routine, "No input fields found in configuration namelist!")

  END SUBROUTINE read_input_namelist


  !> Formatted print-out of input field configuration.
  !
  SUBROUTINE input_print_metadata(field)
    TYPE (t_field_metadata), INTENT(IN) :: field

    WRITE (0,*)        "# input field '"//TRIM(field%inputname)//"' -> '"//TRIM(field%outputname)//"'"
    WRITE (0,'(a,i3)') " #     code          : ", field%code
    WRITE (0,'(a,a)')  " #     type_of_layer : ", TRIM(field%type_of_layer)
  END SUBROUTINE input_print_metadata


  !> @return vlist variable ID for a given input field.
  !
  !  Uses cdilib for file access. Checks for consistency of meta data.
  !
  ! @todo Consistency checks of variable metadata still missing!
  ! 
  FUNCTION get_field_varID(vlistID, field_info) RESULT(result_varID)
    INTEGER                             :: result_varID
    INTEGER,                 INTENT(IN) :: vlistID             !< link to GRIB file vlist
    TYPE (t_field_metadata), INTENT(IN) :: field_info          !< field meta-data

    result_varID = get_varID(vlistID, field_info%inputname)
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
        input_field(i)%varID = get_field_varID(vlistID, input_field(i))
      END IF
      ! distributed computation: broadcast to other worker PEs:
      CALL p_bcast(input_field(i)%varID, rank0, p_comm_work)
      IF (input_field(i)%varID /= -1) THEN
        IF (get_my_mpi_work_id() == rank0) THEN
          ! get the corresponding grid ID
          input_field(i)%gridID  = vlistInqVarGrid( vlistID, input_field(i)%varID)
          input_field(i)%zaxisID = vlistInqVarZaxis(vlistID, input_field(i)%varID)
          ! get the field dimensions:
          input_field(i)%idim(:)  = (/ gridInqXsize(input_field(i)%gridID), &
            &                          gridInqYsize(input_field(i)%gridID) /)
          input_field(i)%nlev     = zaxisInqSize(input_field(i)%zaxisID)
          IF ((input_field(i)%type_of_layer == "surface") .AND.  &
            & (input_field(i)%nlev /= 1)) THEN
            CALL finish(routine, "Inconsistent input data (level no.)")
          END IF
          input_field(i)%steptype = vlistInqVarTsteptype(vlistID, input_field(i)%varID)
          input_field(i)%table    = vlistInqVarTable(vlistID, input_field(i)%varID)
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
        END IF ! rank0
        ! distributed computation: broadcast to other worker PEs:
        CALL p_bcast(input_field(i)%nlev,rank0,p_comm_work)
      END IF ! (varID /= -1)
    END DO

    IF (get_my_mpi_work_id() == rank0) THEN
      ! make sure that there is only one grid definition for regular
      ! grids:
      IF (file%structure == GRID_TYPE_REGULAR) THEN
        DO i=2,n_input_fields    
          IF (input_field(1)%gridID /= input_field(2)%gridID) &
            &   CALL finish(routine, "Only single-grid source file supported!")          
        END DO
      END IF
    END IF

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


  !> Performs horizontal remapping of a single data field.
  !
  !  Requires pre-computed interpolation weights.
  !
  SUBROUTINE input_remap_horz(file_metadata, field_src, field_dst, src_grid, dst_grid, &
    &                         fa, dst_intp_data, gather_c)
    TYPE (t_file_metadata),    INTENT(IN)    :: file_metadata    !< input file meta-data
    REAL(wp),                  INTENT(INOUT) :: field_src(:,:)   !< input field (may be modified)
    REAL(wp),                  INTENT(INOUT) :: field_dst(:,:)   !< output field
    TYPE (t_grid),             INTENT(IN)    :: src_grid         !< source grid     
    TYPE (t_grid),             INTENT(IN)    :: dst_grid         !< destination grid
    TYPE (t_field_adjustment), INTENT(IN)    :: fa               !< meta-data for hydrostatic correction
    TYPE(t_intp_data),         INTENT(IN)    :: dst_intp_data    !< input file meta-data
    TYPE(t_gather_c),          INTENT(IN)    :: gather_c         !< communication pattern

    IF (dbg_level >= 5) &
      & WRITE (0,*) "horizontal remapping..."
    IF (.NOT. fa%lhydrostatic_correction) THEN
      CALL interpolate_c(field_src, field_dst, dst_grid, dst_intp_data)
    ELSE
      CALL hydrostatic_correction(dst_intp_data, src_grid, dst_grid, file_metadata, &
        &                         gather_c, fa, field_src, field_dst)
    END IF

  END SUBROUTINE input_remap_horz


  !> Main subroutine: Reads and remaps data fields.
  !
  !  Requires a pre-allocated variable list ("varlist") as output.
  ! 
  ! @todo Check if grid ID corresponds to the grid for which the
  !       interpolation weights are computed.
  !
  SUBROUTINE input_import_data(file_metadata, ivar, gather_c,                &
    &                          dst_field, src_grid, dst_grid, dst_intp_data, &
    &                          opt_time_comm, opt_time_read, opt_time_intp)
    
    TYPE (t_file_metadata), INTENT(IN) :: file_metadata
    INTEGER,           INTENT(IN)    :: ivar
    TYPE(t_gather_c),  INTENT(IN)    :: gather_c         !< communication pattern
    REAL(wp),          INTENT(INOUT) :: dst_field(:,:,:)
    TYPE (t_grid),     INTENT(IN)    :: src_grid, dst_grid
    TYPE(t_intp_data), INTENT(IN)    :: dst_intp_data
    REAL, INTENT(INOUT), OPTIONAL    :: opt_time_comm, opt_time_read, opt_time_intp
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//':input_import_data')
    INTEGER :: streamID, ierrstat, nmiss, nsize, &
      &        shape2D_glb(2), shape2D_loc(2), ilev
    REAL(wp), ALLOCATABLE :: rfield1D(:), rfield2D(:,:), rfield2D_loc(:,:)
    REAL    :: time_comm, time_read, time_intp

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
        nsize = gridInqSize(input_field(ivar)%gridID)
        IF (file_metadata%structure == GRID_TYPE_REGULAR) THEN
          IF (nsize /= (input_field(ivar)%idim(1)*input_field(ivar)%idim(2))) &
            &   CALL finish(routine, "Internal error!") ! paranoia
        END IF

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
        END IF ! if rank0
        ! for distributed computation: scatter interpolation input to work PEs
        CALL tic(time_comm)  ! performance measurement: start
        CALL scatter_field2D_c(gather_c, rfield2D, rfield2D_loc)
        IF (PRESENT(opt_time_comm)) opt_time_comm = opt_time_comm + toc(time_comm)

        ! perform horizontal interpolation
        CALL tic(time_intp)  ! performance measurement: start
        CALL input_remap_horz(file_metadata, rfield2D_loc, dst_field(:,ilev,:), &
          &                   src_grid, dst_grid, input_field(ivar)%fa,         &
          &                   dst_intp_data, gather_c)
        IF (PRESENT(opt_time_intp)) opt_time_intp = opt_time_intp + toc(time_intp)
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
      IF (ALLOCATED(zaxis_metadata(i)%vct))    DEALLOCATE(zaxis_metadata(i)%vct,    STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
      IF (ALLOCATED(zaxis_metadata(i)%levels)) DEALLOCATE(zaxis_metadata(i)%levels, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END DO
    DEALLOCATE(input_field, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END SUBROUTINE close_input

END MODULE mo_remap_input
