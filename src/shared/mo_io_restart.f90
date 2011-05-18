MODULE mo_io_restart
  !
  USE mo_kind,          ONLY: wp
  USE mo_exception,     ONLY: finish, message, message_text
  USE mo_util_string,   ONLY: separator
  USE mo_util_sysinfo,  ONLY: util_user_name, util_os_system, util_node_name 
  USE mo_util_symlink,  ONLY: util_symlink, util_rename, util_islink, util_unlink
  USE mo_var_metadata,  ONLY: t_var_metadata
  USE mo_linked_list,   ONLY: t_list_element
  USE mo_var_list,      ONLY: t_var_list, nvar_lists, var_lists
  USE mo_cdi_constants
#ifndef NOMPI
  USE mo_mpi,           ONLY: p_parallel_io
#endif
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: set_restart_time
  PUBLIC :: set_restart_vct
  PUBLIC :: init_restart
  PUBLIC :: open_writing_restart_files
  PUBLIC :: write_restart
  PUBLIC :: close_writing_restart_files
  PUBLIC :: open_reading_restart_files
!  PUBLIC :: read_restart
!  PUBLIC :: close_reading_restart_files
  PUBLIC :: cleanup_restart
  !
  TYPE t_gl_att
    CHARACTER(len=64) :: name
    CHARACTER(len=64) :: val
  END TYPE t_gl_att
  !
  INTEGER, PARAMETER :: max_gl_atts = 128
  INTEGER, SAVE :: nglob_atts = 0 
  TYPE(t_gl_att), POINTER :: global_restart_attributes(:) => NULL()
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
  TYPE(t_v_grid) :: vgrid_def(3)
  !
  TYPE t_t_axis
    INTEGER :: type
  END type t_t_axis
  !
  INTEGER, SAVE :: nt_axis = 0 
  TYPE(t_t_axis) :: taxis_def(2)
  !
  CHARACTER(len=32) :: private_restart_time = '' 
  REAL(wp), ALLOCATABLE :: private_vct(:)
  !
  LOGICAL, SAVE :: lvct_initialised = .FALSE. 
  LOGICAL, SAVE :: lrestart_initialised = .FALSE. 
  !
  INTEGER :: private_nc  = -1
  INTEGER :: private_ncv = -1
  INTEGER :: private_nv  = -1
  INTEGER :: private_nvv = -1
  INTEGER :: private_ne  = -1
  INTEGER :: private_nev = -1
  !
#ifdef NOMPI
  LOGICAL :: p_parallel_io = .TRUE.
#endif
  !
  INTERFACE gather_cells
    MODULE PROCEDURE gather_cells_2d
    MODULE PROCEDURE gather_cells_3d
  END INTERFACE gather_cells
  !
  INTERFACE gather_edges
    MODULE PROCEDURE gather_edges_2d
    MODULE PROCEDURE gather_edges_3d
  END INTERFACE gather_edges
  !
  INTERFACE gather_vertices
    MODULE PROCEDURE gather_vertices_2d
    MODULE PROCEDURE gather_vertices_3d
  END INTERFACE gather_vertices
  !
  !------------------------------------------------------------------------------------------------
CONTAINS
  !------------------------------------------------------------------------------------------------
  !
  ! YYYYMMDDThhmmssZ (T is a separator and Z means UTC as timezone)
  !
  SUBROUTINE set_restart_time(iso8601)
    CHARACTER(len=*), INTENT(in) :: iso8601
    private_restart_time = iso8601
  END SUBROUTINE set_restart_time
  !------------------------------------------------------------------------------------------------
  !
  !  VCT as in echam (first half of vector contains a and second half b
  !
  SUBROUTINE set_restart_vct(vct)
    REAL(wp), INTENT(in) :: vct(:)
    IF (lvct_initialised) RETURN
    ALLOCATE(private_vct(SIZE(vct)))
    private_vct(:) = vct(:)
    lvct_initialised = .TRUE.
  END SUBROUTINE set_restart_vct
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE set_restart_attribute(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in) :: attribute_name
    CHARACTER(len=*), INTENT(in) :: attribute_value
    IF (.NOT. ASSOCIATED(global_restart_attributes)) THEN
      ALLOCATE(global_restart_attributes(max_gl_atts))
    ENDIF
    nglob_atts = nglob_atts+1
    IF (nglob_atts > max_gl_atts) THEN
      CALL finish('set_restart_attributes','too many global attributes for restart file')
    ELSE
      global_restart_attributes(nglob_atts)%name = attribute_name
      global_restart_attributes(nglob_atts)%val = attribute_value
    ENDIF    
  END SUBROUTINE set_restart_attribute
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
  SUBROUTINE init_restart(model_name, model_version, nc, ncv, nv, nvv, ne, nev, nlev)
    CHARACTER(len=*), INTENT(in) :: model_name
    CHARACTER(len=*), INTENT(in) :: model_version
    INTEGER,          INTENT(in) :: nc
    INTEGER,          INTENT(in) :: ncv
    INTEGER,          INTENT(in) :: nv
    INTEGER,          INTENT(in) :: nvv
    INTEGER,          INTENT(in) :: ne
    INTEGER,          INTENT(in) :: nev
    INTEGER,          INTENT(in) :: nlev
    !
    CHARACTER(len=256) :: executable
    CHARACTER(len=256) :: user_name
    CHARACTER(len=256) :: os_name
    CHARACTER(len=256) :: host_name
    CHARACTER(len=256) :: tmp_string
    CHARACTER(len=  8) :: date_string
    CHARACTER(len= 10) :: time_string

    INTEGER :: nlena, nlenb, nlenc, nlend
    !
    IF (lrestart_initialised) RETURN
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
    ! set CD-Convention required global attributes
    !
    CALL set_restart_attribute('title', 'ICON simulation')
    CALL set_restart_attribute('institution', &
         &                     'Max Planck Institute for Meteorology/Deutscher Wetterdienst')
    CALL set_restart_attribute('source', model_name//'-'//model_version)
    CALL set_restart_attribute('history', &
         &              executable(1:nlend)//' at '//date_string(1:8)//' '//time_string(1:6))
    CALL set_restart_attribute('references', 'see MPIM/DWD publications')
    CALL set_restart_attribute('comment', &
         &                TRIM(user_name)//' on '//TRIM(host_name)//' ('//TRIM(os_name)//')')
    !
    ! define horizontal grids
    !
    CALL set_horizontal_grid(GRID_UNSTRUCTURED_CELL, nc, ncv)
    CALL set_horizontal_grid(GRID_UNSTRUCTURED_VERT, nv, nvv)
    CALL set_horizontal_grid(GRID_UNSTRUCTURED_EDGE, ne, nev)
    !
    ! define vertical grids 
    !LK test for vct ...
    !
    CALL set_vertical_grid(ZAXIS_SURFACE,     1)
    CALL set_vertical_grid(ZAXIS_HYBRID,      nlev)
    CALL set_vertical_grid(ZAXIS_HYBRID_HALF, nlev+1)
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
    lrestart_initialised = .TRUE.
    !
  END SUBROUTINE init_restart
  !------------------------------------------------------------------------------------------------
  !
  ! Loop over all the output streams and open the associated files. Set
  ! unit numbers (file IDs) for all streams associated with a file.
  !
  SUBROUTINE open_writing_restart_files(basename)
    CHARACTER(len=*), INTENT(in) :: basename
    !
    CHARACTER(len=1024) :: restart_filename
    INTEGER :: status, i ,j, k, ia, ihg, ivg, nlevp1
    REAL(wp), ALLOCATABLE :: levels(:)
    !
    ! first set restart file name
    !
    IF (private_restart_time == '') THEN
      CALL finish('open_restart_files','restart time not set')
    ENDIF
    !
    ! print header for restart var_lists
    !
    CALL message('',separator)
    CALL message('','')
    CALL message('','Open restart files:')
    CALL message('','')
    WRITE(message_text,'(t1,a,t50,a,t84,a,t94,a)') 'file', 'var list', 'file ID', 'restart'
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
      IF (var_lists(i)%p%opened) CYCLE
      !
      ! skip, if var_list is not required for restarting
      !
      IF (.NOT. var_lists(i)%p%lrestart) CYCLE
      !
      ! check restart file type
      !
      SELECT CASE (var_lists(i)%p%restart_type)       
      CASE (FILETYPE_NC)
        CALL finish('open_restart_files','netCDF classic not supported')
      CASE (FILETYPE_NC2, FILETYPE_NC4)
        ! this is ok, both formats can write more than 2GB files
      CASE default
        CALL finish('open_restart_files','unknown restart_type')
      END SELECT
      !
      var_lists(i)%p%first = .TRUE.
      !
      restart_filename = basename//'_'//TRIM(private_restart_time) &
           &                     //'_'//TRIM(var_lists(i)%p%model_type)//'.nc'
      !
      IF (p_parallel_io) THEN
        SELECT CASE (var_lists(i)%p%restart_type)
        CASE (FILETYPE_NC2)
          var_lists(i)%p%cdiFileID = streamOpenWrite(restart_filename, FILETYPE_NC2)
        CASE (FILETYPE_NC4)
          var_lists(i)%p%cdiFileID = streamOpenWrite(restart_filename, FILETYPE_NC4)
        END SELECT
        !
        var_lists(i)%p%filename = TRIM(restart_filename)
        !
        IF (var_lists(i)%p%cdiFileID < 0) THEN
          WRITE(message_text,'(a)') cdiStringError(var_lists(i)%p%cdiFileID)
          CALL message('',message_text)
          CALL finish ('open_restart_files', 'open failed on '//TRIM(restart_filename))
        ELSE
          var_lists(i)%p%opened = .TRUE.
        END IF
        !
        ! The following sections add the file global properties collected in init_restart
        !
        ! 1. create cdi vlist 
        !
        var_lists(i)%p%cdiVlistID = vlistCreate()
        !
        !    set cdi internal time index to 0 for writing time slices in netCDF
        !
        var_lists(i)%p%cdiTimeIndex = 0
        !
        ! 2. add global variables
        !
        DO ia = 1, nglob_atts
          status = vlistDefAttTxt(var_lists(i)%p%cdiVlistID, CDI_GLOBAL,       &
               &                  TRIM(global_restart_attributes(ia)%name),    &
               &                  LEN_TRIM(global_restart_attributes(ia)%val), &
               &                  TRIM(global_restart_attributes(ia)%val))
        END DO
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
            var_lists(i)%p%cdiHalfZaxisID  = zaxisCreate(ZAXIS_HYBRID_HALF, vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(vgrid_def(ivg)%nlevels))
            DO k = 1, vgrid_def(ivg)%nlevels
              levels(k) = REAL(k,wp)
            END DO
            CALL zaxisDefLevels(var_lists(i)%p%cdiHalfZaxisID, levels)
            DEALLOCATE(levels)
            nlevp1 = vgrid_def(ivg)%nlevels
            CALL zaxisDefVct(var_lists(i)%p%cdiHalfZaxisID, 2*nlevp1, private_vct(1:2*nlevp1))
          END SELECT
        ENDDO
        !
        ! 5. restart does contain absolute time 
        !
        var_lists(i)%p%cdiTaxisID = taxisCreate(TAXIS_ABSOLUTE)
        CALL vlistDefTaxis(var_lists(i)%p%cdiVlistID, var_lists(i)%p%cdiTaxisID)
      ENDIF
      !
      ! add variables
      !
      IF (p_parallel_io) THEN
        !
        CALL addVarListToVlist(var_lists(i), var_lists(i)%p%cdiVlistID)
        !
        WRITE(message_text,'(t1,a49,t50,a31,t84,i6,t94,l5)')        &
             restart_filename, var_lists(i)%p%name,             &
             var_lists(i)%p%cdiFileID, var_lists(i)%p%lrestart
        CALL message('',message_text)
      ENDIF
      !
      ! loop over all other output var_lists eventually corresponding to the same file
      !
      DO j = 1, nvar_lists
        !
        IF (var_lists(j)%p%opened) CYCLE
        IF (.NOT. var_lists(j)%p%lrestart) CYCLE
        !
        IF (var_lists(j)%p%restart_type /= var_lists(i)%p%restart_type) THEN
          CALL finish('open_output_streams', 'different file types for the same restart file')
        ENDIF
        !
        IF (var_lists(i)%p%model_type == var_lists(j)%p%model_type) THEN 
          var_lists(j)%p%opened = .TRUE.
          var_lists(j)%p%filename = var_lists(i)%p%filename          
          !
          ! set file IDs of all associated restart files
          !
          var_lists(j)%p%cdiFileID      = var_lists(i)%p%cdiFileID
          var_lists(j)%p%cdiVlistID     = var_lists(i)%p%cdiVlistID
          var_lists(j)%p%cdiCellGridID  = var_lists(i)%p%cdiCellGridID
          var_lists(j)%p%cdiVertGridID  = var_lists(i)%p%cdiVertGridID
          var_lists(j)%p%cdiEdgeGridID  = var_lists(i)%p%cdiEdgeGridID
          var_lists(j)%p%cdiSurfZaxisID = var_lists(i)%p%cdiSurfZaxisID 
          var_lists(j)%p%cdiFullZaxisID = var_lists(i)%p%cdiFullZaxisID 
          var_lists(j)%p%cdiHalfZaxisID = var_lists(i)%p%cdiHalfZaxisID 
          var_lists(j)%p%cdiTaxisID     = var_lists(i)%p%cdiTaxisID
          !
          ! add variables to already existing cdi vlists
          !
          IF (p_parallel_io) THEN
            !
            CALL addVarListToVlist(var_lists(j), var_lists(j)%p%cdiVlistID)
            !
            WRITE(message_text,'(t1,a49,t50,a31,t84,i6,t94,l5)')        &
                 restart_filename, var_lists(j)%p%name,             &
                 var_lists(j)%p%cdiFileID, var_lists(j)%p%lrestart
            CALL message('',message_text)
          ENDIF
        ENDIF
      ENDDO
      !      
      IF (p_parallel_io .AND. var_lists(i)%p%first) THEN
        CALL streamDefVlist(var_lists(i)%p%cdiFileID, var_lists(i)%p%cdiVlistID)
      ENDIF
      !
    END DO
    !    
    CALL message('','')
    !
  END SUBROUTINE open_writing_restart_files
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
      IF (.NOT. info%lrestart) CYCLE
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
          casted_missval = info%missval%ival
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
  SUBROUTINE close_writing_restart_files
    !
    ! Loop over all the output streams and close the associated files, set
    ! opened to false
    !
    CHARACTER(len=80) :: linkname
    INTEGER :: i, iret, fileID
    !
    CALL message('',separator)
    CALL message('','')
    CALL message('','Close restart files:')
    CALL message('','')
    WRITE(message_text,'(t1,a,t50,a,t84,a)') 'file', 'link target', 'file ID'
    CALL message('',message_text)
    CALL message('','')
    !
    close_all_lists: DO i = 1, nvar_lists
      !
      IF (var_lists(i)%p%opened) THEN
        IF (p_parallel_io .AND. var_lists(i)%p%first) THEN
          !
          fileID = var_lists(i)%p%cdiFileID
          !
          CALL streamClose(var_lists(i)%p%cdiFileID)
          CALL vlistDestroy(var_lists(i)%p%cdiVlistID)
          !
          linkname = 'restart_'//TRIM(var_lists(i)%p%model_type)//'.nc'
          IF (util_islink(TRIM(linkname))) THEN
            iret = util_unlink(TRIM(linkname))
          ENDIF
          iret = util_symlink(TRIM(var_lists(i)%p%filename),TRIM(linkname))
          !
          WRITE(message_text,'(t1,a,t50,a,t84,i6)') &
               TRIM(var_lists(i)%p%filename), TRIM(linkname), fileID
          CALL message('',message_text)
          !
        ENDIF
        var_lists(i)%p%filename   = ''
      ENDIF
    ENDDO close_all_lists
    CALL message('','')
    !
    ! reset all var list properties related to cdi files
    !
    reset_all_lists: DO i = 1, nvar_lists
      var_lists(i)%p%cdiFileID = -1
      var_lists(i)%p%opened = .FALSE.
      var_lists(i)%p%first  = .FALSE.
    ENDDO reset_all_lists
    !
    private_restart_time = ''
    !
  END SUBROUTINE close_writing_restart_files
  !------------------------------------------------------------------------------------------------
  !
  ! loop over all var_lists for restart
  !
  SUBROUTINE write_restart
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
          SELECT CASE (var_lists(i)%p%restart_type)
          CASE (FILETYPE_NC2)
            CALL message('','Write netCDF2 restart for : '//TRIM(private_restart_time))
          CASE (FILETYPE_NC4)
            IF (var_lists(i)%p%compression_type == COMPRESS_ZIP) THEN
              CALL message('','Write compressed netCDF4 restart for : '//TRIM(private_restart_time))
            ELSE
              CALL message('','Write netCDF4 restart for : '//TRIM(private_restart_time))
            END IF
          END SELECT
        ENDIF
        write_info = .FALSE.
        !
        ! write time information to netCDF file
        !
        IF (p_parallel_io) THEN
          CALL write_time_to_restart(var_lists(i))
        ENDIF
        !
        ! loop over all streams associated with the file
        !
        DO j = i, nvar_lists
          IF (var_lists(j)%p%cdiFileID == var_lists(i)%p%cdiFileID) THEN 
            !
            !
            ! write variables
            !
            CALL write_restart_var_list(var_lists(j))
            !
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    !
  END SUBROUTINE write_restart
  !------------------------------------------------------------------------------------------------
  !
  ! set time for restart in cdi format
  !
  SUBROUTINE write_time_to_restart (this_list)
    TYPE (t_var_list), INTENT(in) :: this_list
    !
    INTEGER :: fileID, idate, itime, iret
    !
    fileID = this_list%p%cdiFileID
    !
    CALL get_date_components(private_restart_time, idate, itime)
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
  END SUBROUTINE write_time_to_restart
  !------------------------------------------------------------------------------------------------
  !
  ! write variables of a list for restart
  !
  SUBROUTINE write_restart_var_list(this_list)
    TYPE (t_var_list) ,INTENT(in) :: this_list
    !
    INTEGER           :: gridtype
    !
    ! variables of derived type used in linked list
    !
    TYPE (t_var_metadata), POINTER :: info
    TYPE (t_list_element), POINTER :: element
    TYPE (t_list_element), TARGET  :: start_with
    !
    REAL(wp), POINTER :: ptr2d(:,:)   ! 2d field distributed over processors
    REAL(wp), POINTER :: ptr3d(:,:,:) ! 3d field distributed over processors
    !
    REAL(wp), POINTER :: z5d(:,:,:,:,:) ! field gathered on I/O processor
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
      ptr2d => NULL()
      ptr3d => NULL()
      !
      ! retrieve information from actual linked list element
      !
      info => element%field%info 
      !
      ! skip this field ?
      !
      IF (.NOT. info%lrestart) CYCLE
      !
      IF (info%lcontained) THEN 
        nindex = info%ncontained
      ELSE
        nindex = 1
      ENDIF
      !
      SELECT CASE (info%ndims)
      CASE (1)
        CALL finish('write_restart_var_list','1d arrays not handled yet.')
      CASE (2)
        ptr2d => element%field%r_ptr(:,:,nindex,1,1)
      CASE (3)
        ptr3d => element%field%r_ptr(:,:,:,nindex,1)
      CASE (4)
        CALL finish('write_restart_var_list','4d arrays not handled yet.')
      CASE (5)
        CALL finish('write_restart_var_list','5d arrays not handled yet.')
      CASE DEFAULT 
        CALL finish('write_restart_var_list','dimension not set.')        
      END SELECT
      !
      gridtype = info%hgrid
      !
      ! allocate temporary global array on output processor
      ! and gather field from other processors
      !
      NULLIFY(z5d)
      !
      SELECT CASE (gridtype)
      CASE (GRID_UNSTRUCTURED_CELL)
        IF (info%ndims == 2) THEN
          gdims(:) = (/ private_nc, 1, 1, 1, 1 /)
          ALLOCATE(z5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_cells(ptr2d, z5d, info%name)
        ELSE
          gdims(:) = (/ private_nc, info%used_dimensions(2), 1, 1, 1 /)
          ALLOCATE(z5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_cells(ptr3d, z5d, info%name)
        ENDIF
      CASE (GRID_UNSTRUCTURED_VERT)
        IF (info%ndims == 2) THEN
          gdims(:) = (/ private_nv, 1, 1, 1, 1 /)
          ALLOCATE(z5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_vertices(ptr2d, z5d, info%name)
        ELSE
          gdims(:) = (/ private_nv, info%used_dimensions(2), 1, 1, 1 /)
          ALLOCATE(z5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_vertices(ptr3d, z5d, info%name)
        ENDIF
      CASE (GRID_UNSTRUCTURED_EDGE)
        IF (info%ndims == 2) THEN
          gdims(:) = (/ private_ne, 1, 1, 1, 1 /)
          ALLOCATE(z5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_edges(ptr2d, z5d, info%name)
        ELSE
          gdims(:) = (/ private_ne, info%used_dimensions(2), 1, 1, 1 /)
          ALLOCATE(z5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_edges(ptr3d, z5d, info%name)
        ENDIF
      CASE default
        CALL finish('out_stream','unknown grid type')
      END SELECT
      !
      ! write data
      !
      IF (p_parallel_io) THEN
        CALL write_var (this_list, info, z5d)
      END IF
      !
      ! deallocate temporary global arrays
      !
      IF (ASSOCIATED (z5d)) DEALLOCATE (z5d)
      !
    END DO for_all_list_elements
    !
  END SUBROUTINE write_restart_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! finally write data ...
  !
  SUBROUTINE write_var (this_list, info, array)
    TYPE (t_var_list),     INTENT(in) :: this_list
    TYPE (t_var_metadata), INTENT(in) :: info             ! field description
    REAL(wp),              INTENT(in) :: array(:,:,:,:,:) ! restart field
    !
    INTEGER :: fileID                       ! File ID
    INTEGER :: varID                        ! Variable ID
    !
    fileID  = this_list%p%cdiFileID
    varID   = info%cdiVarID
    !
    CALL streamWriteVar(fileID, varID, array, 0)
    !
  END SUBROUTINE write_var
  !------------------------------------------------------------------------------------------------
  !
  ! deallocate module variables
  !
  SUBROUTINE cleanup_restart

    INTEGER :: i

    DO i = 1, nvar_lists    
      IF (var_lists(i)%p%cdiFileID >= 0) THEN
        CALL gridDestroy(var_lists(i)%p%cdiCellGridID)
        CALL gridDestroy(var_lists(i)%p%cdiVertGridID)
        CALL gridDestroy(var_lists(i)%p%cdiEdgeGridID)
        CALL zaxisDestroy(var_lists(i)%p%cdiSurfZaxisID)
        CALL zaxisDestroy(var_lists(i)%p%cdiFullZaxisID)
        CALL zaxisDestroy(var_lists(i)%p%cdiHalfZaxisID)
        var_lists(i)%p%cdiFileId      = CDI_UNDEFID
        var_lists(i)%p%cdiVlistId     = CDI_UNDEFID
        var_lists(i)%p%cdiCellGridID  = CDI_UNDEFID
        var_lists(i)%p%cdiVertGridID  = CDI_UNDEFID
        var_lists(i)%p%cdiEdgeGridID  = CDI_UNDEFID
        var_lists(i)%p%cdiSurfZaxisID = CDI_UNDEFID
        var_lists(i)%p%cdiHalfZaxisID = CDI_UNDEFID
        var_lists(i)%p%cdiFullZaxisID = CDI_UNDEFID
        var_lists(i)%p%cdiTaxisID     = CDI_UNDEFID
        var_lists(i)%p%cdiTimeIndex   = CDI_UNDEFID
      ENDIF
    ENDDO
    !
    DEALLOCATE(private_vct)
    lvct_initialised = .FALSE.
    !
    nglob_atts = 0
    nh_grids   = 0
    nv_grids   = 0
    nt_axis    = 0
    lrestart_initialised = .FALSE.
    !
  END SUBROUTINE cleanup_restart
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE gather_cells_2d(in_array, out_array, name)
    REAL(wp),                   INTENT(in) :: in_array(:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z1d(:)
    z1d => out_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder_2d(in_array, z1d)
#endif
  END SUBROUTINE gather_cells_2d
  !
  SUBROUTINE gather_cells_3d(in_array, out_array, name)
    REAL(wp),                   INTENT(in) :: in_array(:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z2d(:,:)
    z2d => out_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder_3d(in_array, z2d)
#endif
  END SUBROUTINE gather_cells_3d
  !
  SUBROUTINE gather_vertices_2d(in_array, out_array, name)
    REAL(wp),                   INTENT(in) :: in_array(:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z1d(:)
    z1d => out_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder_2d(in_array, z1d)
#endif
  END SUBROUTINE gather_vertices_2d
  !
  SUBROUTINE gather_vertices_3d(in_array, out_array, name)
    REAL(wp),                   INTENT(in) :: in_array(:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z2d(:,:)
    z2d => out_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder_3d(in_array, z2d)
#endif
  END SUBROUTINE gather_vertices_3d
  !
  SUBROUTINE gather_edges_2d(in_array, out_array, name)
    REAL(wp),                   INTENT(in) :: in_array(:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z1d(:)
    z1d => out_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder_2d(in_array, z1d)
#endif
  END SUBROUTINE gather_edges_2d
  !
  SUBROUTINE gather_edges_3d(in_array, out_array, name)
    REAL(wp),                   INTENT(in) :: in_array(:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z2d(:,:)
    z2d => out_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder_3d(in_array, z2d)
#endif
  END SUBROUTINE gather_edges_3d
  !
  SUBROUTINE reorder_2d(in, out)
    REAL(wp), INTENT(in)    :: in(:,:)
    REAL(wp), INTENT(inout) :: out(:)
    !
    LOGICAL, ALLOCATABLE    :: lmask(:)
    INTEGER ::  isize_in, isize_out
    INTEGER :: idiscrep
    !
    isize_in  = SIZE(in)
    isize_out = SIZE(out)
    idiscrep = isize_in-isize_out

    IF(idiscrep == 0 )THEN
      out = RESHAPE(in,(/ isize_out /))
    ELSE
      ALLOCATE (lmask(isize_in))
      lmask(1:isize_out) = .TRUE.
      lmask(isize_out+1:isize_in) = .FALSE.
      out = PACK(RESHAPE(in,(/isize_in/)),lmask)
      DEALLOCATE (lmask)
    ENDIF

  END SUBROUTINE reorder_2d
  !
  SUBROUTINE reorder_3d(in, out)
    REAL(wp), INTENT(in)    :: in(:,:,:)
    REAL(wp), INTENT(inout) :: out(:,:)
    !
    LOGICAL, ALLOCATABLE    :: lmask(:)
    INTEGER ::isize_in, isize_out, isize_lev
    INTEGER :: idiscrep, k
    !
    isize_in  = SIZE(in,1)*SIZE(in,3)
    isize_out = SIZE(out,1)
    isize_lev = SIZE(in,2)
    idiscrep = isize_in-isize_out
    !
    IF (idiscrep /= 0 )THEN
      ALLOCATE (lmask(isize_in))
      lmask(1:isize_out) = .TRUE.
      lmask(isize_out+1:isize_in) = .FALSE.
    ENDIF
    !
    DO k = 1, isize_lev
      IF (idiscrep /= 0 )THEN
        out(:,k) = PACK(RESHAPE(in(:,k,:),(/isize_in/)),lmask)
      ELSE
        out(:,k) =      RESHAPE(in(:,k,:),(/isize_out/))
      ENDIF
    ENDDO
    !   
    IF (idiscrep /= 0 )THEN
      DEALLOCATE (lmask)
    ENDIF
    !
  END SUBROUTINE reorder_3d
  !
  SUBROUTINE open_reading_restart_files(model_type)
    CHARACTER(len=*), INTENT(in) :: model_type
    !
!!$    TYPE (t_list_element), POINTER :: element
!!$    TYPE (t_list_element), TARGET  :: start_with
    !
    CHARACTER(len=80) :: restart_filename, name
    !
    INTEGER :: fileID, vlistID, taxisID, varID
    INTEGER :: idate, itime
    !
    INTEGER :: iret
    !
    restart_filename = 'restart_'//TRIM(model_type)//'.nc'
    !
    IF (.NOT. util_islink(TRIM(restart_filename))) THEN
      iret = util_rename(TRIM(restart_filename), TRIM(restart_filename)//'.bak')
      iret = util_symlink(TRIM(restart_filename)//'.bak', TRIM(restart_filename))
    ENDIF
    !
    fileID  = streamOpenRead(restart_filename)
    vlistID = streamInqVlist(fileID)
    taxisID = vlistInqTaxis(vlistID)
    !
    idate = taxisInqVdate(taxisID)
    itime = taxisInqVtime(taxisID)
    !
!    write(0,*) idate, itime
    DO varID = 0, vlistNvars(vlistID)-1
      CALL vlistInqVarName(vlistID, varID, name)
!      write (0,*) ' ... read name '//TRIM(name)
    ENDDO
    !
!!$    for_all_lists: DO i = 1, nvar_lists
!!$
!!$      element => start_with
!!$      element%next_list_element => var_lists(i)%p%first_list_element
!!$      !
!!$      for_all_list_elements: DO
!!$        element => element%next_list_element
!!$        IF (.NOT.ASSOCIATED(element)) EXIT
!!$        !
!!$
!!$
!!$        !    vlistInqVarName
!!$
!!$        !    streamReadVar
!!$
!!$      ENDDO for_all_list_elements
!!$    ENDDO for_all_lists

    CALL streamClose(fileID)
!    CALL vlistDestroy(vlistID)

  END SUBROUTINE open_reading_restart_files
  !
END MODULE mo_io_restart
