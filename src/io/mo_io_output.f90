MODULE mo_io_output
  !
  USE mo_kind,          ONLY: wp
  USE mo_exception,     ONLY: finish, message, message_text
  USE mo_var_metadata,  ONLY: t_var_metadata
  USE mo_linked_list,   ONLY: t_list_element
  USE mo_var_list,      ONLY: t_var_list, nvar_lists, var_lists
  USE mo_cdi_constants
  USE mo_util_string,   ONLY: separator
  USE mo_util_sysinfo,  ONLY: util_user_name, util_os_system, util_node_name 
  USE mo_io_distribute, ONLY: gather_cells, gather_edges, gather_vertices
  USE mo_mpi,           ONLY: my_process_is_stdio
#ifndef NOMPI
  USE mo_model_domain,  ONLY: t_patch
#endif
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: set_output_time
  PUBLIC :: set_output_vct, set_output_depth, set_output_height
  PUBLIC :: init_output
  PUBLIC :: open_output_files
  PUBLIC :: write_output
  PUBLIC :: close_output_files
  PUBLIC :: finish_output
  !
  TYPE t_output_files
    CHARACTER(len=64) :: functionality
    CHARACTER(len=64) :: filename
  END type t_output_files
  INTEGER, PARAMETER :: max_output_files = 257
  INTEGER, SAVE :: noutput_files = 0 
  TYPE(t_output_files), ALLOCATABLE :: output_files(:)
  !
  TYPE t_h_grid
    INTEGER :: type
    INTEGER :: nelements
    INTEGER :: nvertices
  END type t_h_grid
  !
  INTEGER, SAVE :: nh_grids = 0 
  TYPE(t_h_grid) :: hgrid_def(3)
  !
  TYPE t_v_grid
    INTEGER :: type
    INTEGER :: nlevels
  END type t_v_grid
  !
  INTEGER, SAVE :: nv_grids = 0 
  TYPE(t_v_grid) :: vgrid_def(7)
  !
  TYPE t_t_axis
    INTEGER :: type
  END type t_t_axis
  !
  INTEGER, SAVE :: nt_axis = 0 
  TYPE(t_t_axis) :: taxis_def(2)
  !
  CHARACTER(len=32) :: private_initial_time = '' 
  CHARACTER(len=32) :: private_output_time = '' 
  REAL(wp), ALLOCATABLE :: private_vct(:)
  REAL(wp), ALLOCATABLE :: private_depth_full(:),  private_depth_half(:)
  REAL(wp), ALLOCATABLE :: private_height_full(:),  private_height_half(:) 
  !
  LOGICAL, SAVE :: lvct_initialised = .FALSE. 
  LOGICAL, SAVE :: ldepth_initialised = .FALSE. 
  LOGICAL, SAVE :: lheight_initialised = .FALSE. 
  !
  LOGICAL, SAVE :: loutput_initialised = .FALSE. 
  !
  INTEGER :: private_nc  = -1
  INTEGER :: private_ncv = -1
  INTEGER :: private_nv  = -1
  INTEGER :: private_nvv = -1
  INTEGER :: private_ne  = -1
  INTEGER :: private_nev = -1
  !
  CHARACTER(len=1024) :: private_base_filename
  !
  !
  !------------------------------------------------------------------------------------------------
CONTAINS
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE set_output_filenames(functionality, filename)
    CHARACTER(len=*), INTENT(in) :: functionality
    CHARACTER(len=*), INTENT(in) :: filename
    IF (.NOT. ALLOCATED(output_files)) THEN
      ALLOCATE(output_files(max_output_files))
    ENDIF
    noutput_files = noutput_files+1
    IF (noutput_files > max_output_files) THEN
      CALL finish('set_output_filenames','too many output filenames')
    ELSE
      output_files(noutput_files)%functionality = functionality
      output_files(noutput_files)%filename      = filename
    ENDIF    
  END SUBROUTINE set_output_filenames
  !
  SUBROUTINE get_output_filenames()
  END SUBROUTINE get_output_filenames
  !
  !------------------------------------------------------------------------------------------------
  !
  ! YYYYMMDDThhmmssZ (T is a separator and Z means UTC as timezone)
  !
  SUBROUTINE set_output_time(iso8601)
    CHARACTER(len=*), INTENT(in) :: iso8601
    IF (private_initial_time == '') THEN
      private_initial_time = iso8601
    ELSE
      private_output_time = iso8601
    ENDIF
  END SUBROUTINE set_output_time
  !------------------------------------------------------------------------------------------------
  !
  !  VCT as in echam (first half of vector contains a and second half b
  !
  SUBROUTINE set_output_vct(vct)
    REAL(wp), INTENT(in) :: vct(:)
    IF (lvct_initialised) RETURN
    ALLOCATE(private_vct(SIZE(vct)))
    private_vct(:) = vct(:)
    lvct_initialised = .TRUE.
  END SUBROUTINE set_output_vct
  !------------------------------------------------------------------------------------------------
  !
  !  height based vertical coordinates
  !
  SUBROUTINE set_output_height(zh, zf)
    REAL(wp), INTENT(in) :: zh(:), zf(:)
    IF (lheight_initialised) RETURN
    ALLOCATE(private_height_half(SIZE(zh)), private_height_full(SIZE(zf)))
    private_height_half(:) = zh(:)
    private_height_full(:) = zf(:)
    lheight_initialised = .TRUE.
  END SUBROUTINE set_output_height
  !------------------------------------------------------------------------------------------------
  !
  !  depth based vertical coordinates
  !
  SUBROUTINE set_output_depth(zh, zf)
    REAL(wp), INTENT(in) :: zh(:), zf(:)
    IF (ldepth_initialised) RETURN
    ALLOCATE(private_depth_half(SIZE(zh)), private_depth_full(SIZE(zf)))
    private_depth_half(:) = zh(:)
    private_depth_full(:) = zf(:)
    ldepth_initialised = .TRUE.
  END SUBROUTINE set_output_depth
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE set_horizontal_grid(grid_type, nelements, nvertices)
    INTEGER, INTENT(in) :: grid_type
    INTEGER, INTENT(in) :: nelements
    INTEGER, INTENT(in) :: nvertices
    !    
    nh_grids = nh_grids+1
    !
    hgrid_def(nh_grids)%type      = grid_type
    hgrid_def(nh_grids)%nelements = nelements
    hgrid_def(nh_grids)%nvertices = nvertices
    !
  END SUBROUTINE set_horizontal_grid
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE set_vertical_grid(type, nlevels)
    INTEGER, INTENT(in) :: type
    INTEGER, INTENT(in) :: nlevels
    !    
    nv_grids = nv_grids+1
    !
    vgrid_def(nv_grids)%type    = type
    vgrid_def(nv_grids)%nlevels = nlevels
    !
  END SUBROUTINE set_vertical_grid
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE set_time_axis(type)
    INTEGER, INTENT(in) :: type
    !    
    nt_axis = nt_axis+1
    !
    taxis_def(nt_axis)%type = type
    !
  END SUBROUTINE set_time_axis
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE init_output(model_name, model_version, &
       &                 experiment_id,             &
       &                 domain, root, refinement,  & 
       &                 nc, ncv, nv, nvv, ne, nev, &
       &                 nlev, ndepth, nheight)
    CHARACTER(len=*), INTENT(in) :: model_name
    CHARACTER(len=*), INTENT(in) :: model_version
    CHARACTER(len=*), INTENT(in) :: experiment_id
    INTEGER,          INTENT(in) :: domain
    INTEGER,          INTENT(in) :: root
    INTEGER,          INTENT(in) :: refinement
    INTEGER,          INTENT(in) :: nc
    INTEGER,          INTENT(in) :: ncv
    INTEGER,          INTENT(in) :: nv
    INTEGER,          INTENT(in) :: nvv
    INTEGER,          INTENT(in) :: ne
    INTEGER,          INTENT(in) :: nev
    INTEGER,          INTENT(in) :: nlev
    INTEGER,          INTENT(in) :: ndepth
    INTEGER,          INTENT(in) :: nheight
    !
    CHARACTER(len= 256) :: executable
    CHARACTER(len= 256) :: user_name
    CHARACTER(len= 256) :: os_name
    CHARACTER(len= 256) :: host_name
    CHARACTER(len= 256) :: tmp_string
    CHARACTER(len=   8) :: date_string
    CHARACTER(len=  10) :: time_string
    !
    INTEGER :: nlena, nlenb, nlenc, nlend
    !
    IF (loutput_initialised) RETURN
    !
    IF (private_initial_time == '') THEN
      CALL finish('open_output_files','initial time not set')
    ENDIF

    ! set base output filename  
    !
    WRITE(private_base_filename,'(a,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,a)') &
         experiment_id,                                                &
         '_D', domain, 'R', root, 'B', refinement, 'L', nlev,          &
         '_', private_initial_time
write (0,'(/,a,a,a,a,a,a/)') &
     'LK: ', model_name, '-', model_version, ': ', TRIM(private_base_filename)
    !
    executable = ''
    user_name  = ''
    os_name    = ''
    host_name  = ''
    !
    CALL get_command_argument(0, executable, nlend)
    CALL date_and_time(date_string, time_string)
    !
    tmp_string = ''
    CALL util_os_system (tmp_string, nlena)
    os_name = tmp_string(1:nlena)

    tmp_string = ''
    CALL util_user_name (tmp_string, nlenb)
    user_name = tmp_string(1:nlenb)
    
    tmp_string = ''
    CALL util_node_name (tmp_string, nlenc)
    host_name = tmp_string(1:nlenc)    
    !
    ! set CD-Convention required output attributes
    !
!!$    CALL set_output_attribute('title',       &
!!$         'ICON simulation')
!!$    CALL set_output_attribute('institution', &
!!$         'Max Planck Institute for Meteorology/Deutscher Wetterdienst')
!!$    CALL set_output_attribute('source',      &
!!$         model_name//'-'//model_version)
!!$    CALL set_output_attribute('history',     &
!!$         executable(1:nlend)//' at '//date_string(1:8)//' '//time_string(1:6))
!!$    CALL set_output_attribute('references',  &
!!$         'see MPIM/DWD publications')
!!$    CALL set_output_attribute('comment',     &
!!$         TRIM(user_name)//' on '//TRIM(host_name)//' ('//TRIM(os_name)//')')
    !
    ! define horizontal grids
    !
    CALL set_horizontal_grid(GRID_UNSTRUCTURED_CELL, nc, ncv)
    CALL set_horizontal_grid(GRID_UNSTRUCTURED_VERT, nv, nvv)
    CALL set_horizontal_grid(GRID_UNSTRUCTURED_EDGE, ne, nev)
    !
    ! define vertical grids 
    !
    CALL set_vertical_grid(ZAXIS_SURFACE,     1)
    CALL set_vertical_grid(ZAXIS_HYBRID,      nlev)
    CALL set_vertical_grid(ZAXIS_HYBRID_HALF, nlev+1)
    CALL set_vertical_grid(ZAXIS_DEPTH_BELOW_SEA, ndepth)
    CALL set_vertical_grid(ZAXIS_DEPTH_BELOW_SEA, ndepth+1)
    CALL set_vertical_grid(ZAXIS_HEIGHT, nheight)
    CALL set_vertical_grid(ZAXIS_HEIGHT, nheight+1)    !
    !
    ! define time axis
    !
    CALL set_time_axis(TAXIS_ABSOLUTE)
    CALL set_time_axis(TAXIS_RELATIVE)
    !
    private_nc  = nc
    private_ncv = ncv
    private_nv  = nv
    private_nvv = nvv
    private_ne  = ne
    private_nev = nev
    !
    loutput_initialised = .TRUE.
    !
  END SUBROUTINE init_output
  !------------------------------------------------------------------------------------------------
  !
  ! Loop over all the output streams and open the associated files. Set
  ! unit numbers (file IDs) for all streams associated with a file.
  !
  SUBROUTINE open_output_files
    !
    CHARACTER(len=1024) :: output_filename
    INTEGER :: i ,j, k, ihg, ivg, nlevp1
    REAL(wp), ALLOCATABLE :: levels(:)
    !
    ! first check for output time
    !
    IF (private_initial_time == '') THEN
      CALL finish('open_output_files','initial time not set')
    ENDIF
    !
    ! print header for output var_lists
    !
    CALL message('',separator)
    CALL message('','')
    CALL message('','Open output files:')
    CALL message('','')
    WRITE(message_text,'(t1,a,t50,a,t84,a,t94,a)') 'file', 'var list', 'file ID', 'output'
    CALL message('',message_text)
    CALL message('','')
    !
    ! Loop over all var_lists and open the associated files. Set
    ! file IDs if necessary.
    !
    DO i = 1, nvar_lists
      var_lists(i)%p%first = .FALSE.
    ENDDO
    !
    DO i = 1, nvar_lists
      !
      ! skip, if file is already opened
      !
      IF (var_lists(i)%p%output_opened) CYCLE
      !
      ! skip, if var_list is not required for outputing
      !
      IF (.NOT. var_lists(i)%p%loutput) CYCLE
      !
      ! check output file type
      !
      SELECT CASE (var_lists(i)%p%output_type)       
      CASE (FILETYPE_NC)
        CALL finish('open_output_files','netCDF classic not supported')
      CASE (FILETYPE_NC2, FILETYPE_NC4)
        ! this is ok, both formats can write more than 2GB files
      CASE (FILETYPE_GRB)
        CALL finish('open_output_files','GRIB1 not supported')
      CASE (FILETYPE_GRB2)
        CALL message('open_output_files','GRIB2 support experimental')
      CASE default
        CALL finish('open_output_files','unknown output_type')
      END SELECT
      !
      var_lists(i)%p%first = .TRUE.
      !
      output_filename = TRIM(private_base_filename)      &
           &           //'_'//TRIM(var_lists(i)%p%model_type)//'.nc'
      !
      IF (my_process_is_stdio()) THEN
        SELECT CASE (var_lists(i)%p%output_type)
        CASE (FILETYPE_NC2)
          var_lists(i)%p%cdiFileID_output = streamOpenWrite(output_filename, FILETYPE_NC2)
        CASE (FILETYPE_NC4)
          var_lists(i)%p%cdiFileID_output = streamOpenWrite(output_filename, FILETYPE_NC4)
        CASE (FILETYPE_GRB2)
          var_lists(i)%p%cdiFileID_output = streamOpenWrite(output_filename, FILETYPE_GRB2)
        END SELECT
        !
        var_lists(i)%p%filename = TRIM(output_filename)
        !
        IF (var_lists(i)%p%cdiFileID_output < 0) THEN
          WRITE(message_text,'(a)') cdiStringError(var_lists(i)%p%cdiFileID_output)
          CALL message('',message_text)
          CALL finish ('open_output_files', 'open failed on '//TRIM(output_filename))
        ELSE
          CALL message ('LK: ', 'open succedded on '//TRIM(output_filename))
write (0,*) 'LK: fileid - ', var_lists(i)%p%cdiFileID_output           
          var_lists(i)%p%output_opened = .TRUE.
        END IF
        !
        ! The following sections add the file global properties collected in init_output
        !
        ! 1. create cdi vlist 
        !
        var_lists(i)%p%cdiVlistID = vlistCreate()
        !
        !    set cdi internal time index to 0 for writing time slices in netCDF
        !
        var_lists(i)%p%cdiTimeIndex = 0
        !
        ! 2. add global attributes for netCDF
        !
        ! 3. add horizontal grid descriptions
        !
        DO ihg = 1, nh_grids
          SELECT CASE (hgrid_def(ihg)%type)
          CASE (GRID_UNSTRUCTURED_CELL)
            var_lists(i)%p%cdiCellGridID = gridCreate(GRID_UNSTRUCTURED, hgrid_def(ihg)%nelements)
            CALL gridDefNvertex(var_lists(i)%p%cdiCellGridID, hgrid_def(ihg)%nvertices)
            !
            CALL gridDefXname(var_lists(i)%p%cdiCellGridID, 'clon')
            CALL gridDefXlongname(var_lists(i)%p%cdiCellGridID, 'center longitude')
            CALL gridDefXunits(var_lists(i)%p%cdiCellGridID, 'radians')
            !
            CALL gridDefYname(var_lists(i)%p%cdiCellGridID, 'clat')
            CALL gridDefYlongname(var_lists(i)%p%cdiCellGridID, 'center latitude')
            CALL gridDefYunits(var_lists(i)%p%cdiCellGridID, 'radians')
            !
          CASE (GRID_UNSTRUCTURED_VERT)
            var_lists(i)%p%cdiVertGridID = gridCreate(GRID_UNSTRUCTURED, hgrid_def(ihg)%nelements)
            CALL gridDefNvertex(var_lists(i)%p%cdiVertGridID, hgrid_def(ihg)%nvertices)
            !
            CALL gridDefXname(var_lists(i)%p%cdiVertGridID, 'vlon')
            CALL gridDefXlongname(var_lists(i)%p%cdiVertGridID, 'vertex longitude')
            CALL gridDefXunits(var_lists(i)%p%cdiVertGridID, 'radians')
            !
            CALL gridDefYname(var_lists(i)%p%cdiVertGridID, 'vlat')
            CALL gridDefYlongname(var_lists(i)%p%cdiVertGridID, 'vertex latitude')
            CALL gridDefYunits(var_lists(i)%p%cdiVertGridID, 'radians')
            !
          CASE (GRID_UNSTRUCTURED_EDGE)
            var_lists(i)%p%cdiEdgeGridID = gridCreate(GRID_UNSTRUCTURED, hgrid_def(ihg)%nelements)
            CALL gridDefNvertex(var_lists(i)%p%cdiEdgeGridID, hgrid_def(ihg)%nvertices)
            !
            CALL gridDefXname(var_lists(i)%p%cdiEdgeGridID, 'elon')
            CALL gridDefXlongname(var_lists(i)%p%cdiEdgeGridID, 'edge longitude')
            CALL gridDefXunits(var_lists(i)%p%cdiEdgeGridID, 'radians')
            !
            CALL gridDefYname(var_lists(i)%p%cdiEdgeGridID, 'elat')
            CALL gridDefYlongname(var_lists(i)%p%cdiEdgeGridID, 'edge latitude')
            CALL gridDefYunits(var_lists(i)%p%cdiEdgeGridID, 'radians')
            !
          END SELECT
        ENDDO
        !
        ! 4. add vertical grid descriptions
        !
        DO ivg = 1, nv_grids
          SELECT CASE (vgrid_def(ivg)%type)
          CASE (ZAXIS_SURFACE)
            var_lists(i)%p%cdiSurfZaxisID = zaxisCreate(ZAXIS_SURFACE, vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(1))
            levels(1) = 0.0_wp
            CALL zaxisDefLevels(var_lists(i)%p%cdiSurfZaxisID, levels)
            DEALLOCATE(levels)
          CASE (ZAXIS_HYBRID)
            IF (.NOT. lvct_initialised) CYCLE
            var_lists(i)%p%cdiFullZaxisID = zaxisCreate(ZAXIS_HYBRID, vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(vgrid_def(ivg)%nlevels))
            DO k = 1, vgrid_def(ivg)%nlevels
              levels(k) = REAL(k,wp)
            END DO
            CALL zaxisDefLevels(var_lists(i)%p%cdiFullZaxisID, levels)
            DEALLOCATE(levels)
            nlevp1 = vgrid_def(ivg)%nlevels+1
            CALL zaxisDefVct(var_lists(i)%p%cdiFullZaxisID, 2*nlevp1, private_vct(1:2*nlevp1))
          CASE (ZAXIS_HYBRID_HALF)
            IF (.NOT. lvct_initialised) CYCLE
            var_lists(i)%p%cdiHalfZaxisID  = zaxisCreate(ZAXIS_HYBRID_HALF, vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(vgrid_def(ivg)%nlevels))
            DO k = 1, vgrid_def(ivg)%nlevels
              levels(k) = REAL(k,wp)
            END DO
            CALL zaxisDefLevels(var_lists(i)%p%cdiHalfZaxisID, levels)
            DEALLOCATE(levels)
            nlevp1 = vgrid_def(ivg)%nlevels
            CALL zaxisDefVct(var_lists(i)%p%cdiHalfZaxisID, 2*nlevp1, private_vct(1:2*nlevp1))
          CASE (ZAXIS_DEPTH_BELOW_SEA)
            IF (.NOT. ldepth_initialised) CYCLE
            IF (SIZE(private_depth_full) == vgrid_def(ivg)%nlevels) THEN
              var_lists(i)%p%cdiDepthFullZaxisID = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, &
                   &                                           vgrid_def(ivg)%nlevels)
              CALL zaxisDefLevels(var_lists(i)%p%cdiDepthFullZaxisID, &
                   &              private_depth_full)
            ELSE IF (SIZE(private_depth_half) == vgrid_def(ivg)%nlevels) THEN
              var_lists(i)%p%cdiDepthHalfZaxisID = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, &
                   &                                           vgrid_def(ivg)%nlevels)
              CALL zaxisDefLevels(var_lists(i)%p%cdiDepthHalfZaxisID, &
                   &              private_depth_half)
            ELSE
              CALL finish('open_output_files','Number of depth levels not available.')
            ENDIF
          CASE (ZAXIS_HEIGHT)
            IF (.NOT. lheight_initialised) CYCLE
            IF (SIZE(private_height_full) == vgrid_def(ivg)%nlevels) THEN
              var_lists(i)%p%cdiHeightFullZaxisID = zaxisCreate(ZAXIS_HEIGHT, &
                   &                                            vgrid_def(ivg)%nlevels)
              CALL zaxisDefLevels(var_lists(i)%p%cdiHeightFullZaxisID, &
                   &              private_height_full)
            ELSE IF (SIZE(private_height_half) == vgrid_def(ivg)%nlevels) THEN
              var_lists(i)%p%cdiHeightHalfZaxisID = zaxisCreate(ZAXIS_HEIGHT, &
                   &                                            vgrid_def(ivg)%nlevels)
              CALL zaxisDefLevels(var_lists(i)%p%cdiHeightHalfZaxisID, &
                   &              private_height_half)
            ELSE
              CALL finish('open_output_files','Number of height levels not available.')
            ENDIF
          END SELECT
        ENDDO
        !
        ! 5. output does contain absolute time 
        !
        var_lists(i)%p%cdiTaxisID = taxisCreate(TAXIS_ABSOLUTE)
        CALL vlistDefTaxis(var_lists(i)%p%cdiVlistID, var_lists(i)%p%cdiTaxisID)
      ENDIF
      !
      ! add variables
      !
      IF (my_process_is_stdio()) THEN
        !
        CALL addVarListToVlist(var_lists(i), var_lists(i)%p%cdiVlistID)
        !
        WRITE(message_text,'(t1,a49,t50,a31,t84,i6,t94,l5)')        &
             output_filename, var_lists(i)%p%name,             &
             var_lists(i)%p%cdiFileID_output, var_lists(i)%p%loutput
        CALL message('',message_text)
      ENDIF
      !
      ! loop over all other output var_lists eventually corresponding to the same file
      !
      DO j = 1, nvar_lists
        !
        IF (var_lists(j)%p%output_opened) CYCLE
        IF (.NOT. var_lists(j)%p%loutput) CYCLE
        !
        IF (var_lists(j)%p%output_type /= var_lists(i)%p%output_type &
             .AND. var_lists(j)%p%cdiFileID_output == var_lists(i)%p%cdiFileID_output) THEN
          CALL finish('open_output_streams', 'different file types for the same output file')
        ENDIF
        !
        IF (var_lists(i)%p%model_type == var_lists(j)%p%model_type) THEN 
          var_lists(j)%p%output_opened = .TRUE.
          var_lists(j)%p%filename = var_lists(i)%p%filename          
          !
          ! set file IDs of all associated output files
          !
          var_lists(j)%p%cdiFileID_output     = var_lists(i)%p%cdiFileID_output
          var_lists(j)%p%cdiVlistID           = var_lists(i)%p%cdiVlistID
          var_lists(j)%p%cdiCellGridID        = var_lists(i)%p%cdiCellGridID
          var_lists(j)%p%cdiVertGridID        = var_lists(i)%p%cdiVertGridID
          var_lists(j)%p%cdiEdgeGridID        = var_lists(i)%p%cdiEdgeGridID
          var_lists(j)%p%cdiSurfZaxisID       = var_lists(i)%p%cdiSurfZaxisID 
          var_lists(j)%p%cdiFullZaxisID       = var_lists(i)%p%cdiFullZaxisID 
          var_lists(j)%p%cdiHalfZaxisID       = var_lists(i)%p%cdiHalfZaxisID 
          var_lists(j)%p%cdiDepthFullZaxisID  = var_lists(i)%p%cdiDepthFullZaxisID 
          var_lists(j)%p%cdiDepthHalfZaxisID  = var_lists(i)%p%cdiDepthHalfZaxisID 
          var_lists(j)%p%cdiHeightFullZaxisID = var_lists(i)%p%cdiHeightFullZaxisID 
          var_lists(j)%p%cdiHeightHalfZaxisID = var_lists(i)%p%cdiHeightHalfZaxisID
          var_lists(j)%p%cdiTaxisID           = var_lists(i)%p%cdiTaxisID
          !
          ! add variables to already existing cdi vlists
          !
          IF (my_process_is_stdio()) THEN
            !
            CALL addVarListToVlist(var_lists(j), var_lists(j)%p%cdiVlistID)
            !
            WRITE(message_text,'(t1,a49,t50,a31,t84,i6,t94,l5)')        &
                 output_filename, var_lists(j)%p%name,             &
                 var_lists(j)%p%cdiFileID_output, var_lists(j)%p%loutput
            CALL message('',message_text)
          ENDIF
        ENDIF
      ENDDO
      !      
      IF (my_process_is_stdio() .AND. var_lists(i)%p%first) THEN
        CALL streamDefVlist(var_lists(i)%p%cdiFileID_output, var_lists(i)%p%cdiVlistID)
      ENDIF
      !
    END DO
    !    
    CALL message('','')
    !
  END SUBROUTINE open_output_files
  !------------------------------------------------------------------------------------------------
  !
  ! define variables and attributes
  !
  SUBROUTINE addVarListToVlist(this_list, vlistID)
    TYPE (t_var_list), INTENT(inout) :: this_list
    INTEGER,           INTENT(inout) :: vlistID
    !      
    TYPE (t_var_metadata), POINTER :: info
    TYPE (t_list_element), POINTER :: element
    TYPE (t_list_element), TARGET  :: start_with
    !
    INTEGER :: varID, gridID, zaxisID
    !
    REAL(wp) :: casted_missval
    !
    element => start_with
    element%next_list_element => this_list%p%first_list_element
    !
    for_all_list_elements: DO
      !
      element => element%next_list_element
      IF (.NOT.ASSOCIATED(element)) EXIT
      !
      ! retrieve information from actual linked list element
      !
      info => element%field%info
      !
      ! skip this field ?
      !
      IF (.NOT. info%loutput) CYCLE
      !
      ! set grid ID
      !
      SELECT CASE (info%hgrid)
      CASE(GRID_UNSTRUCTURED_CELL)
        info%cdiGridID = this_list%p%cdiCellGridID
        gridID = info%cdiGridID
      CASE(GRID_UNSTRUCTURED_VERT)
        info%cdiGridID = this_list%p%cdiVertGridID
        gridID = info%cdiGridID
      CASE(GRID_UNSTRUCTURED_EDGE)
        info%cdiGridID = this_list%p%cdiEdgeGridID
        gridID = info%cdiGridID
      END SELECT
      !
      ! set z axis ID
      !
      SELECT CASE (info%vgrid)
      CASE (ZAXIS_SURFACE)
        info%cdiZaxisID =  this_list%p%cdiSurfZaxisID
        zaxisID = info%cdiZaxisID
      CASE (ZAXIS_HYBRID)
        info%cdiZaxisID =  this_list%p%cdiFullZaxisID
        zaxisID = info%cdiZaxisID
      CASE (ZAXIS_HYBRID_HALF)
        info%cdiZaxisID =  this_list%p%cdiHalfZaxisID
        zaxisID = info%cdiZaxisID
      CASE (ZAXIS_DEPTH_BELOW_SEA)
        IF (info%used_dimensions(2) == SIZE(private_depth_half)) THEN
          info%cdiZaxisID =  this_list%p%cdiDepthHalfZaxisID
          zaxisID = info%cdiZaxisID
        ELSE IF (info%used_dimensions(2) == SIZE(private_depth_full)) THEN
          info%cdiZaxisID =  this_list%p%cdiDepthFullZaxisID
          zaxisID = info%cdiZaxisID
        ENDIF
      CASE (ZAXIS_HEIGHT)
        IF (info%used_dimensions(2) == SIZE(private_height_half)) THEN
          info%cdiZaxisID =  this_list%p%cdiHeightHalfZaxisID
          zaxisID = info%cdiZaxisID
        ELSE IF (info%used_dimensions(2) == SIZE(private_height_full)) THEN
          info%cdiZaxisID =  this_list%p%cdiHeightFullZaxisID
          zaxisID = info%cdiZaxisID
        ENDIF
      END SELECT
      !
      IF ( gridID  == -1 ) THEN
        CALL finish('addStreamToVlist', 'GRID definition missing for '//TRIM(info%name))
      END IF
      IF ( zaxisID == -1 ) THEN
        CALL finish('addStreamToVlist', 'ZAXIS definition missing for '//TRIM(info%name))
      END IF
      !
      info%cdiVarID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE)
      varID = info%cdiVarID 
      !
      CALL vlistDefVarDatatype(vlistID, varID, DATATYPE_FLT64)
      CALL vlistDefVarName(vlistID, varID, info%name)
      !
      IF (info%cf%long_name /= '') CALL vlistDefVarLongname(vlistID, varID, info%cf%long_name)
      IF (info%cf%units /= '') CALL vlistDefVarUnits(vlistID, varID, info%cf%units)

      IF (info%lmiss) THEN
        IF (ASSOCIATED(element%field%r_ptr)) THEN
          casted_missval = info%missval%rval
        ELSE IF (ASSOCIATED(element%field%i_ptr)) THEN
          casted_missval = REAL(info%missval%ival,wp)
        ELSE
          IF (info%missval%lval) THEN
            casted_missval = 1.0_wp
          ELSE
            casted_missval = 0.0_wp
          ENDIF
        ENDIF
        CALL vlistDefVarMissval(vlistID, varID, casted_missval)
      ENDIF
      !
    ENDDO for_all_list_elements
    !
  END SUBROUTINE addVarListToVlist
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE close_output_files
    !
    ! Loop over all the output streams and close the associated files, set
    ! opened to false
    !
    INTEGER :: i, j, fileID, vlistID
    !
    CALL message('',separator)
    CALL message('','')
    CALL message('','Close output files:')
    CALL message('','')
    WRITE(message_text,'(t1,a,t84,a)') 'file', 'file ID'
    CALL message('',message_text)
    CALL message('','')
    !
    close_all_lists: DO i = 1, nvar_lists
      !
      IF (var_lists(i)%p%output_opened) THEN
        IF (my_process_is_stdio() .AND. var_lists(i)%p%first) THEN
          !
          fileID = var_lists(i)%p%cdiFileID_output
          !
          IF (fileID /= CDI_UNDEFID) THEN
            CALL streamClose(var_lists(i)%p%cdiFileID_output)
            DO j = 1, nvar_lists
              IF (fileID == var_lists(j)%p%cdiFileID_output) THEN
                var_lists(j)%p%cdiFileID_output = CDI_UNDEFID
              ENDIF
            ENDDO
          ENDIF
          !
          WRITE(message_text,'(t1,a,t84,i6)') &
               TRIM(var_lists(i)%p%filename), fileID
          CALL message('',message_text)
          !
        ENDIF
        var_lists(i)%p%filename   = ''
      ENDIF
    ENDDO close_all_lists
    !
    for_all_vlists: DO i = 1, nvar_lists
      vlistID = var_lists(i)%p%cdiVlistID
      IF (vlistID /= CDI_UNDEFID) THEN
        CALL vlistDestroy(var_lists(i)%p%cdiVlistID)
        DO j = 1, nvar_lists
          IF (vlistID == var_lists(j)%p%cdiVlistID) THEN
            var_lists(j)%p%cdiVlistID = CDI_UNDEFID
          ENDIF
        ENDDO
      ENDIF
    ENDDO for_all_vlists
    CALL message('','')
    !
    ! reset all var list properties related to cdi files
    !
    reset_all_lists: DO i = 1, nvar_lists
      var_lists(i)%p%output_opened = .FALSE.
      var_lists(i)%p%first  = .FALSE.
    ENDDO reset_all_lists
    !
    private_output_time = ''
    !
  END SUBROUTINE close_output_files
  !------------------------------------------------------------------------------------------------
  !
  ! loop over all var_lists for output
  !
  SUBROUTINE write_output(p_patch)
#ifndef NOMPI
    TYPE(t_patch), OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,       OPTIONAL, INTENT(in) :: p_patch
#endif
    INTEGER :: i,j
    LOGICAL :: write_info
    !
    write_info   = .TRUE.
    !
    ! pick up first stream associated with each file
    !
    DO i = 1, nvar_lists
      IF (var_lists(i)%p%first) THEN
        IF (write_info) THEN
          SELECT CASE (var_lists(i)%p%output_type)
          CASE (FILETYPE_NC2)
            CALL message('','Write netCDF2 output for : '//TRIM(private_output_time))
          CASE (FILETYPE_NC4)
            IF (var_lists(i)%p%compression_type == COMPRESS_ZIP) THEN
              CALL message('', &
                   'Write compressed netCDF4 output for : '//TRIM(private_output_time))
            ELSE
              CALL message('','Write netCDF4 output for : '//TRIM(private_output_time))
            END IF
          CASE (FILETYPE_GRB2)
            IF (var_lists(i)%p%compression_type == COMPRESS_SZIP) THEN
              CALL message('','Write compressed GRIB2 output for : '//TRIM(private_output_time))
            ELSE
              CALL message('','Write GRIB2 output for : '//TRIM(private_output_time))
            ENDIF
          END SELECT
        ENDIF
        write_info = .FALSE.
        !
        ! write time information to netCDF file
        !
        IF (my_process_is_stdio()) THEN
          CALL write_time_to_output(var_lists(i))
        ENDIF
        !
        ! loop over all streams associated with the file
        !
        DO j = i, nvar_lists

          IF (var_lists(j)%p%cdiFileID_output == var_lists(i)%p%cdiFileID_output) THEN 
            !
            !
            ! write variables
            !
            CALL write_output_var_list(var_lists(j), p_patch=p_patch)
            !
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    !
  END SUBROUTINE write_output
  !------------------------------------------------------------------------------------------------
  !
  ! set time for output in cdi format
  !
  SUBROUTINE write_time_to_output (this_list)
    TYPE (t_var_list), INTENT(inout) :: this_list
    !
    INTEGER :: fileID, idate, itime, iret
    !
    fileID = this_list%p%cdiFileID_output
    !
    CALL get_date_components(private_output_time, idate, itime)
    !
    CALL taxisDefVdate(this_list%p%cdiTaxisID, idate)
    CALL taxisDefVtime(this_list%p%cdiTaxisID, itime)
    !
    iret = streamDefTimestep(fileID, this_list%p%cdiTimeIndex)
    this_list%p%cdiTimeIndex = this_list%p%cdiTimeIndex + 1
    ! 
  CONTAINS
    !
    SUBROUTINE get_date_components(iso8601, idate, itime)
      CHARACTER(len=*), INTENT(in)  :: iso8601
      INTEGER,          INTENT(out) :: idate, itime
      !
      INTEGER :: it, iz 
      !
      it = INDEX(iso8601, 'T')
      iz = INDEX(iso8601, 'Z')
      READ(iso8601(1:it-1), '(i10)') idate 
      READ(iso8601(it+1:iz-1), '(i10)') itime 
      !
    END SUBROUTINE get_date_components
    !
  END SUBROUTINE write_time_to_output
  !------------------------------------------------------------------------------------------------
  !
  ! write variables of a list for output
  !
  SUBROUTINE write_output_var_list(this_list, p_patch)
    TYPE (t_var_list) ,INTENT(in) :: this_list
#ifndef NOMPI
    TYPE(t_patch), OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,       OPTIONAL, INTENT(in) :: p_patch
#endif
    !
    INTEGER           :: gridtype
    !
    ! variables of derived type used in linked list
    !
    TYPE (t_var_metadata), POINTER :: info
    TYPE (t_list_element), POINTER :: element
    TYPE (t_list_element), TARGET  :: start_with
    !
    REAL(wp), POINTER :: rptr2d(:,:)   ! 2d field distributed over processors
    REAL(wp), POINTER :: rptr3d(:,:,:) ! 3d field distributed over processors
    !
    REAL(wp), POINTER :: r5d(:,:,:,:,:) ! field gathered on I/O processor
    !
    INTEGER :: gdims(5), nindex
    !
    ! Loop over all fields in linked list
    !
    element => start_with
    element%next_list_element => this_list%p%first_list_element    
    !
    for_all_list_elements: DO
      !
      element => element%next_list_element
      IF (.NOT.ASSOCIATED(element)) EXIT
      !
      rptr2d => NULL()
      rptr3d => NULL()
      !
      ! retrieve information from actual linked list element
      !
      info => element%field%info 
      !
      ! skip this field ?
      !
      IF (.NOT. info%loutput) CYCLE
      !
      IF (info%lcontained) THEN 
        nindex = info%ncontained
      ELSE
        nindex = 1
      ENDIF
      !
      SELECT CASE (info%ndims)
      CASE (1)
        CALL finish('write_output_var_list','1d arrays not handled yet.')
      CASE (2)
        rptr2d => element%field%r_ptr(:,:,nindex,1,1)
      CASE (3)
        rptr3d => element%field%r_ptr(:,:,:,nindex,1)
      CASE (4)
        CALL finish('write_output_var_list','4d arrays not handled yet.')
      CASE (5)
        CALL finish('write_output_var_list','5d arrays not handled yet.')
      CASE DEFAULT 
        CALL finish('write_output_var_list','dimension not set.')        
      END SELECT
      !
      gridtype = info%hgrid
      !
      ! allocate temporary global array on output processor
      ! and gather field from other processors
      !
      r5d => NULL()
      !
      SELECT CASE (gridtype)
      CASE (GRID_UNSTRUCTURED_CELL)
        IF (info%ndims == 2) THEN
          gdims(:) = (/ private_nc, 1, 1, 1, 1 /)
          ALLOCATE(r5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_cells(rptr2d, r5d, p_patch=p_patch)
        ELSE
          gdims(:) = (/ private_nc, info%used_dimensions(2), 1, 1, 1 /)
          ALLOCATE(r5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_cells(rptr3d, r5d, p_patch=p_patch)
        ENDIF
      CASE (GRID_UNSTRUCTURED_VERT)
        IF (info%ndims == 2) THEN
          gdims(:) = (/ private_nv, 1, 1, 1, 1 /)
          ALLOCATE(r5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_vertices(rptr2d, r5d, p_patch=p_patch)
        ELSE
          gdims(:) = (/ private_nv, info%used_dimensions(2), 1, 1, 1 /)
          ALLOCATE(r5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_vertices(rptr3d, r5d, p_patch=p_patch)
        ENDIF
      CASE (GRID_UNSTRUCTURED_EDGE)
        IF (info%ndims == 2) THEN
          gdims(:) = (/ private_ne, 1, 1, 1, 1 /)
          ALLOCATE(r5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_edges(rptr2d, r5d, p_patch=p_patch)
        ELSE
          gdims(:) = (/ private_ne, info%used_dimensions(2), 1, 1, 1 /)
          ALLOCATE(r5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_edges(rptr3d, r5d, p_patch=p_patch)
        ENDIF
      CASE default
        CALL finish('out_stream','unknown grid type')
      END SELECT
      !
      ! write data
      !
      IF (my_process_is_stdio()) THEN
        CALL write_var (this_list, info, r5d)
      END IF
      !
      ! deallocate temporary global arrays
      !
      IF (ASSOCIATED (r5d)) DEALLOCATE (r5d)
      !
    END DO for_all_list_elements
    !
  END SUBROUTINE write_output_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! finally write data ...
  !
  SUBROUTINE write_var (this_list, info, array)
    TYPE (t_var_list),     INTENT(in) :: this_list
    TYPE (t_var_metadata), INTENT(in) :: info             ! field description
    REAL(wp),              INTENT(in) :: array(:,:,:,:,:) ! output field
    !
    INTEGER :: fileID                       ! File ID
    INTEGER :: varID                        ! Variable ID
    INTEGER :: nmiss = 0
    !
    fileID  = this_list%p%cdiFileID_output
    varID   = info%cdiVarID
    !
    CALL streamWriteVar(fileID, varID, array, nmiss)
    !
  END SUBROUTINE write_var
  !------------------------------------------------------------------------------------------------
  !
  ! deallocate module variables
  !
  SUBROUTINE finish_output

    INTEGER :: i

    for_all_var_lists: DO i = 1, nvar_lists    
      IF (var_lists(i)%p%cdiFileID_output >= 0) THEN
        CALL gridDestroy(var_lists(i)%p%cdiCellGridID)
        CALL gridDestroy(var_lists(i)%p%cdiVertGridID)
        CALL gridDestroy(var_lists(i)%p%cdiEdgeGridID)
        IF (var_lists(i)%p%cdiSurfZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiSurfZaxisID)
        IF (var_lists(i)%p%cdiFullZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiFullZaxisID)
        IF (var_lists(i)%p%cdiHalfZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiHalfZaxisID)
        IF (var_lists(i)%p%cdiDepthFullZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiDepthFullZaxisID)
        IF (var_lists(i)%p%cdiDepthHalfZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiDepthHalfZaxisID)
        IF (var_lists(i)%p%cdiHeightFullZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiHeightFullZaxisID)
        IF (var_lists(i)%p%cdiHeightHalfZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiHeightHalfZaxisID)
        var_lists(i)%p%cdiFileId_output     = CDI_UNDEFID
        var_lists(i)%p%cdiVlistId           = CDI_UNDEFID
        var_lists(i)%p%cdiCellGridID        = CDI_UNDEFID
        var_lists(i)%p%cdiVertGridID        = CDI_UNDEFID
        var_lists(i)%p%cdiEdgeGridID        = CDI_UNDEFID
        var_lists(i)%p%cdiSurfZaxisID       = CDI_UNDEFID
        var_lists(i)%p%cdiHalfZaxisID       = CDI_UNDEFID
        var_lists(i)%p%cdiFullZaxisID       = CDI_UNDEFID
        var_lists(i)%p%cdiDepthHalfZaxisID  = CDI_UNDEFID
        var_lists(i)%p%cdiDepthFullZaxisID  = CDI_UNDEFID
        var_lists(i)%p%cdiHeightHalfZaxisID = CDI_UNDEFID
        var_lists(i)%p%cdiHeightFullZaxisID = CDI_UNDEFID
        var_lists(i)%p%cdiTaxisID           = CDI_UNDEFID
        var_lists(i)%p%cdiTimeIndex         = CDI_UNDEFID
      ENDIF
    ENDDO for_all_var_lists
    !
    IF (ALLOCATED(private_vct)) DEALLOCATE(private_vct)
    lvct_initialised = .FALSE.
    IF (ALLOCATED(private_depth_full)) DEALLOCATE(private_depth_full)
    IF (ALLOCATED(private_depth_half)) DEALLOCATE(private_depth_half)
    ldepth_initialised = .FALSE.
    IF (ALLOCATED(private_height_full)) DEALLOCATE(private_height_full)
    IF (ALLOCATED(private_height_half)) DEALLOCATE(private_height_half)
    lheight_initialised = .FALSE.
    !
    nh_grids   = 0
    nv_grids   = 0
    nt_axis    = 0
    loutput_initialised = .FALSE.
    !
  END SUBROUTINE finish_output
  !
END MODULE mo_io_output
