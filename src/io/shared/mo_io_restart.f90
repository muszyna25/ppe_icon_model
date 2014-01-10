!OPTION! -pvctl conflict
MODULE mo_io_restart
  !
  USE mo_kind,                  ONLY: wp
  USE mo_mpi,                   ONLY: p_barrier,p_comm_work
  USE mo_exception,             ONLY: finish, message, message_text, get_filename_noext
  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH
  USE mo_var_metadata,          ONLY: t_var_metadata
  USE mo_linked_list,           ONLY: t_var_list, t_list_element, find_list_element
  USE mo_var_list,              ONLY: nvar_lists, var_lists, get_var_timelevel
  USE mo_cdi_constants
  USE mo_util_string,           ONLY: separator
  USE mo_util_sysinfo,          ONLY: util_user_name, util_os_system, util_node_name
  USE mo_util_file,             ONLY: util_symlink, util_rename, util_islink, util_unlink
  USE mo_util_hash,             ONLY: util_hashword
  USE mo_util_uuid,             ONLY: t_uuid
  USE mo_io_restart_namelist,   ONLY: nmls, restart_namelist
  USE mo_io_restart_attributes, ONLY: set_restart_attribute, get_restart_attribute, &
       &                              read_attributes, delete_attributes,           &
       &                              restart_attributes_count_text,                &
       &                              restart_attributes_count_real,                &
       &                              restart_attributes_count_int,                 &
       &                              restart_attributes_count_bool
  USE mo_scatter,               ONLY: scatter_cells, scatter_edges, scatter_vertices
  USE mo_io_units,              ONLY: find_next_free_unit, filename_max
  USE mo_datetime,              ONLY: t_datetime,iso8601
  USE mo_run_config,            ONLY: ltimer, output_mode
  USE mo_timer,                 ONLY: timer_start, timer_stop,                      &
    &                                 timer_write_restart_file, timer_write_output
  USE mo_math_utilities,        ONLY: set_zlev

  USE mo_dynamics_config,       ONLY: iequations, nold, nnow, nnew, nnew_rcf, nnow_rcf

#ifndef __NO_ICON_ATMO__
!LK comment: should not be here !!!!!! polution of namespace !!!!!!
  USE mo_impl_constants,        ONLY: IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER, &
    &                                 LEAPFROG_EXPL, LEAPFROG_SI
  USE mo_ha_dyn_config,         ONLY: ha_dyn_config
#endif
#ifndef __NO_ICON_OCEAN__
  USE mo_ocean_nml,             ONLY: n_zlev, dzlev_m
#endif

  USE mo_model_domain,          ONLY: t_patch
  USE mo_mpi,                   ONLY: my_process_is_mpi_workroot, &
    & my_process_is_mpi_test
  USE mo_communication,         ONLY: t_comm_gather_pattern, exchange_data

  !
  IMPLICIT NONE
  !
  PRIVATE

  INCLUDE 'netcdf.inc'
  !
  PUBLIC :: set_restart_time
  PUBLIC :: set_restart_vct, set_restart_depth
  PUBLIC :: set_restart_depth_lnd
  PUBLIC :: set_restart_height_snow
  PUBLIC :: init_restart
  PUBLIC :: open_writing_restart_files
  PUBLIC :: create_restart_file
  PUBLIC :: close_writing_restart_files
  PUBLIC :: read_restart_files
  PUBLIC :: finish_restart
  PUBLIC :: write_restart_info_file
  PUBLIC :: read_restart_info_file
  !
  TYPE t_restart_files
    CHARACTER(len=64) :: functionality
    CHARACTER(len=64) :: filename
    CHARACTER(len=64) :: linkname
  END type t_restart_files
  INTEGER, PARAMETER :: max_restart_files = 257
  INTEGER, PARAMETER :: max_vertical_axes = 20
  INTEGER, SAVE :: nrestart_files = 0
  TYPE(t_restart_files), ALLOCATABLE :: restart_files(:)
  !
  TYPE t_h_grid
    INTEGER :: type
    INTEGER :: nelements
    INTEGER :: nvertices
    TYPE(t_uuid) :: uuid
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
  TYPE(t_v_grid) :: vgrid_def(max_vertical_axes)
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
  REAL(wp), ALLOCATABLE :: private_depth_full(:),  private_depth_half(:)
  REAL(wp), ALLOCATABLE :: private_depth_lnd_full(:),  private_depth_lnd_half(:)
  REAL(wp), ALLOCATABLE :: private_generic_level(:)
  REAL(wp), ALLOCATABLE :: private_height_snow_half(:), private_height_snow_full(:)
  !
  LOGICAL, SAVE :: lvct_initialised         = .FALSE.
  LOGICAL, SAVE :: ldepth_initialised       = .FALSE.
  LOGICAL, SAVE :: ldepth_lnd_initialised   = .FALSE.
  LOGICAL, SAVE :: lheight_snow_initialised = .FALSE.
  !
  LOGICAL, SAVE :: lrestart_initialised     = .FALSE.
  !
  INTEGER :: private_nc  = -1
  INTEGER :: private_ncv = -1
  INTEGER :: private_nv  = -1
  INTEGER :: private_nvv = -1
  INTEGER :: private_ne  = -1
  INTEGER :: private_nev = -1
  !
  !
  CHARACTER(len=12), PARAMETER :: restart_info_file = 'restart.info'
  !
  !------------------------------------------------------------------------------------------------
CONTAINS
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE set_restart_filenames(functionality, filename)
    CHARACTER(len=*), INTENT(in) :: functionality
    CHARACTER(len=*), INTENT(in) :: filename
    IF (.NOT. ALLOCATED(restart_files)) THEN
      ALLOCATE(restart_files(max_restart_files))
    ENDIF
    nrestart_files = nrestart_files+1
    IF (nrestart_files > max_restart_files) THEN
      CALL finish('set_restart_filenames','too many restart filenames')
    ELSE
      restart_files(nrestart_files)%functionality = functionality
      restart_files(nrestart_files)%filename      = filename
      ! need to get the links target name
    ENDIF
  END SUBROUTINE set_restart_filenames
  !
  SUBROUTINE get_restart_filenames()
  END SUBROUTINE get_restart_filenames
  !
  SUBROUTINE write_restart_info_file
    INTEGER :: nrf
    nrf = find_next_free_unit(10,100)
    OPEN(nrf,file=restart_info_file)
    WRITE(nrf, '(a,a)')                &
         'gridspec: grid_D01_atm.nc', &
         ' ! here should be the physical filename including path'
    CLOSE(nrf)
    CALL message('',restart_info_file//' written')
  END SUBROUTINE write_restart_info_file
  !
  SUBROUTINE read_restart_info_file(gridspecfile, status)
    LOGICAL,          INTENT(out) :: status
    CHARACTER(len=*), INTENT(out) :: gridspecfile
    !
    LOGICAL :: lexist
    INTEGER :: delimiter, comment
    INTEGER :: nrf, ios
    CHARACTER(len=256) :: iomsg, functionality
    CHARACTER(len=filename_max) :: buffer, line
    !
    CALL message('',restart_info_file//' to be read')
    status = .FALSE. ! initially
    gridspecfile = ''
    !
    INQUIRE(FILE=restart_info_file,exist=lexist)
    IF (.NOT. lexist) THEN
      CALL finish('read_restart_info_file',restart_info_file//' is missing')
    ENDIF
    !
    nrf = find_next_free_unit(10,100)
    ! just to be pedantic: status='old'
    OPEN(nrf,file=restart_info_file,STATUS='OLD')
    for_all_lines: DO
#ifdef __SX__
      READ(nrf, '(a)', IOSTAT=ios) buffer
      iomsg = ' unknown - sxf90 does not support IOMSG'
#else
      READ(nrf, '(a)', IOSTAT=ios, IOMSG=iomsg) buffer
#endif
      IF (ios < 0) THEN
        EXIT ! information missing
      ELSE IF (ios > 0) THEN
        CALL finish('read_restart_info_file',restart_info_file//' read error '//TRIM(iomsg))
      END IF
      line = ADJUSTL(buffer)
      delimiter = INDEX(line,':')
      IF (delimiter == 0) CYCLE
      functionality = TRIM(line(1:delimiter))
      IF (functionality(1:8) == 'gridspec') THEN
        comment = INDEX(line,'!')
        gridspecfile = TRIM(ADJUSTL(line(delimiter+1:comment-1)))
        status = .TRUE.
        EXIT
      ENDIF
    ENDDO for_all_lines
    CLOSE(nrf)
    !
  END SUBROUTINE read_restart_info_file
  !
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
  !  depth based vertical coordinates
  !
  SUBROUTINE set_restart_depth(zh, zf)
    REAL(wp), INTENT(in) :: zh(:), zf(:)
    IF (ldepth_initialised) RETURN
    ALLOCATE(private_depth_half(SIZE(zh)), private_depth_full(SIZE(zf)))
    private_depth_half(:) = zh(:)
    private_depth_full(:) = zf(:)
    ldepth_initialised = .TRUE.
  END SUBROUTINE set_restart_depth
  !------------------------------------------------------------------------------------------------
  !
  !  depth based vertical coordinates
  !
  SUBROUTINE set_restart_depth_lnd(zh, zf)
    REAL(wp), INTENT(in) :: zh(:), zf(:)
    IF (ldepth_lnd_initialised) RETURN
    ALLOCATE(private_depth_lnd_half(SIZE(zh)), private_depth_lnd_full(SIZE(zf)))
    private_depth_lnd_half(:) = zh(:)
    private_depth_lnd_full(:) = zf(:)
    ldepth_lnd_initialised = .TRUE.
  END SUBROUTINE set_restart_depth_lnd
  !------------------------------------------------------------------------------------------------
  !
  !  height based vertical coordinates for multi layer snow model (TERRA)
  !
  SUBROUTINE set_restart_height_snow(zh, zf)
    REAL(wp), INTENT(in) :: zh(:), zf(:)
    IF (lheight_snow_initialised) RETURN
    ALLOCATE(private_height_snow_half(SIZE(zh)), private_height_snow_full(SIZE(zf)))
    private_height_snow_half(:) = zh(:)
    private_height_snow_full(:) = zf(:)
    lheight_snow_initialised = .TRUE.
  END SUBROUTINE set_restart_height_snow
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE set_horizontal_grid(grid_type, nelements, nvertices, grid_uuid)
    INTEGER,      INTENT(in)           :: grid_type
    INTEGER,      INTENT(in)           :: nelements
    INTEGER,      INTENT(in)           :: nvertices
    TYPE(t_uuid), INTENT(in), OPTIONAL :: grid_uuid
    !
    nh_grids = nh_grids+1
    !
    hgrid_def(nh_grids)%type      = grid_type
    hgrid_def(nh_grids)%nelements = nelements
    hgrid_def(nh_grids)%nvertices = nvertices
    !
    IF (PRESENT(grid_uuid)) THEN
      hgrid_def(nh_grids)%uuid = grid_uuid
    ENDIF
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
  SUBROUTINE init_restart(model_name, model_version, &
       &                  nc, ncv, nv, nvv, ne, nev, &
       &                  nlev, ndepth, nlev_soil,   &
       &                  nlev_snow, nice_class)
    CHARACTER(len=*), INTENT(in) :: model_name
    CHARACTER(len=*), INTENT(in) :: model_version
    INTEGER,          INTENT(in) :: nc
    INTEGER,          INTENT(in) :: ncv
    INTEGER,          INTENT(in) :: nv
    INTEGER,          INTENT(in) :: nvv
    INTEGER,          INTENT(in) :: ne
    INTEGER,          INTENT(in) :: nev
    INTEGER,          INTENT(in) :: nlev
    INTEGER,          INTENT(in) :: ndepth
    INTEGER,          INTENT(in) :: nlev_soil
    INTEGER,          INTENT(in) :: nlev_snow
    INTEGER,          INTENT(in) :: nice_class
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
    ! set CD-Convention required restart attributes
    !
    CALL set_restart_attribute('title',       &
         'ICON simulation')
    CALL set_restart_attribute('institution', &
         'Max Planck Institute for Meteorology/Deutscher Wetterdienst')
    CALL set_restart_attribute('source',      &
         model_name//'-'//model_version)
    CALL set_restart_attribute('history',     &
         executable(1:nlend)//' at '//date_string(1:8)//' '//time_string(1:6))
    CALL set_restart_attribute('references',  &
         'see MPIM/DWD publications')
    CALL set_restart_attribute('comment',     &
         TRIM(user_name)//' on '//TRIM(host_name)//' ('//TRIM(os_name)//')')
    !
    ! define horizontal grids
    !
    CALL set_horizontal_grid(GRID_UNSTRUCTURED_CELL, nc, ncv)
    CALL set_horizontal_grid(GRID_UNSTRUCTURED_VERT, nv, nvv)
    CALL set_horizontal_grid(GRID_UNSTRUCTURED_EDGE, ne, nev)
    !
    ! define vertical grids
    !
    CALL set_vertical_grid(ZA_SURFACE             , 1          )
    CALL set_vertical_grid(ZA_HYBRID              , nlev       )
    CALL set_vertical_grid(ZA_HYBRID_HALF         , nlev+1     )
    CALL set_vertical_grid(ZA_DEPTH_BELOW_LAND    , nlev_soil  )
    CALL set_vertical_grid(ZA_DEPTH_BELOW_LAND_P1 , nlev_soil+1)
    CALL set_vertical_grid(ZA_SNOW                , nlev_snow  )
    CALL set_vertical_grid(ZA_SNOW_HALF           , nlev_snow+1)
    CALL set_vertical_grid(ZA_HEIGHT_2M           , 1          )
    CALL set_vertical_grid(ZA_HEIGHT_10M          , 1          )
    CALL set_vertical_grid(ZA_LAKE_BOTTOM         , 1          )
    CALL set_vertical_grid(ZA_LAKE_BOTTOM_HALF    , 1          )
    CALL set_vertical_grid(ZA_MIX_LAYER           , 1          )
    CALL set_vertical_grid(ZA_SEDIMENT_BOTTOM_TW_HALF, 1       )
    CALL set_vertical_grid(ZA_TOA                 , 1          )
    CALL set_vertical_grid(ZA_DEPTH_BELOW_SEA     , ndepth     )
    CALL set_vertical_grid(ZA_DEPTH_BELOW_SEA_HALF, ndepth+1   )
    CALL set_vertical_grid(ZA_GENERIC_ICE         , nice_class )
    CALL set_vertical_grid(ZA_DEPTH_RUNOFF_S      , 1          )
    CALL set_vertical_grid(ZA_DEPTH_RUNOFF_G      , 1          )
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

!AD(9July-2013): The following condition seemed unnecessary. So after
!  discussing with DR we decided to get rid of it for time being

!    IF (.NOT. (lvct_initialised .OR. ldepth_initialised            &
!      & .OR. ldepth_lnd_initialised )) THEN
!      CALL finish('init_restart','none of the vertical grids is initialised')
!      ! more consistency checks need to follow
!    ENDIF

    !
    lrestart_initialised = .TRUE.
    !
  END SUBROUTINE init_restart
  !------------------------------------------------------------------------------------------------


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
    REAL(wp), ALLOCATABLE :: levels(:), ubounds(:), lbounds(:)
    !
    CHARACTER(len=64) :: attribute_name, text_attribute
    REAL(wp) :: real_attribute
    INTEGER :: int_attribute
    LOGICAL :: bool_attribute

    IF (my_process_is_mpi_test()) RETURN
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
      IF (var_lists(i)%p%restart_opened) CYCLE
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
      IF (my_process_is_mpi_workroot()) THEN
        SELECT CASE (var_lists(i)%p%restart_type)
        CASE (FILETYPE_NC2)
          var_lists(i)%p%cdiFileID_restart = streamOpenWrite(restart_filename, FILETYPE_NC2)
        CASE (FILETYPE_NC4)
          var_lists(i)%p%cdiFileID_restart = streamOpenWrite(restart_filename, FILETYPE_NC4)
        END SELECT
        !
        var_lists(i)%p%filename = TRIM(restart_filename)
        !
        IF (var_lists(i)%p%cdiFileID_restart < 0) THEN
          WRITE(message_text,'(a)') cdiStringError(var_lists(i)%p%cdiFileID_restart)
          CALL message('',message_text)
          CALL finish ('open_restart_files', 'open failed on '//TRIM(restart_filename))
        ELSE
          var_lists(i)%p%restart_opened = .TRUE.
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
        ! 2. add global attributes
        !
        ! 2.1 namelists as text attributes
        !
        DO ia = 1, nmls
          status = vlistDefAttTxt(var_lists(i)%p%cdiVlistID, CDI_GLOBAL, &
               &                  TRIM(restart_namelist(ia)%name),       &
               &                  LEN_TRIM(restart_namelist(ia)%text),   &
               &                  TRIM(restart_namelist(ia)%text))
        ENDDO
        !
        ! 2.2 text attributes
        !
        DO ia = 1, restart_attributes_count_text()
          CALL get_restart_attribute(ia, attribute_name, text_attribute)
          status = vlistDefAttTxt(var_lists(i)%p%cdiVlistID, CDI_GLOBAL, &
               &                  TRIM(attribute_name),                  &
               &                  LEN_TRIM(text_attribute),              &
               &                  TRIM(text_attribute))
        END DO
        !
        ! 2.3 real attributes
        !
        DO ia = 1, restart_attributes_count_real()
          CALL get_restart_attribute(ia, attribute_name, real_attribute)
          status = vlistDefAttFlt(var_lists(i)%p%cdiVlistID, CDI_GLOBAL, &
               &                  TRIM(attribute_name),                  &
               &                  DATATYPE_FLT64,                        &
               &                  1,                                     &
               &                  real_attribute)
        ENDDO
        !
        ! 2.4 integer attributes
        !
        DO ia = 1, restart_attributes_count_int()
          CALL get_restart_attribute(ia, attribute_name, int_attribute)
          status = vlistDefAttInt(var_lists(i)%p%cdiVlistID, CDI_GLOBAL, &
               &                  TRIM(attribute_name),                  &
               &                  DATATYPE_INT32,                        &
               &                  1,                                     &
               &                  int_attribute)
        ENDDO
        !
        ! 2.5 logical attributes
        !
        DO ia = 1, restart_attributes_count_bool()
          CALL get_restart_attribute(ia, attribute_name, bool_attribute)
          IF (bool_attribute) THEN
            int_attribute = 1
          ELSE
            int_attribute = 0
          ENDIF
          status = vlistDefAttInt(var_lists(i)%p%cdiVlistID, CDI_GLOBAL, &
               &                  TRIM(attribute_name),                  &
               &                  DATATYPE_INT32,                        &
               &                  1,                                     &
               &                  int_attribute)
        ENDDO
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
            CALL gridDefXlongname(var_lists(i)%p%cdiEdgeGridID, 'edge midpoint longitude')
            CALL gridDefXunits(var_lists(i)%p%cdiEdgeGridID, 'radians')
            !
            CALL gridDefYname(var_lists(i)%p%cdiEdgeGridID, 'elat')
            CALL gridDefYlongname(var_lists(i)%p%cdiEdgeGridID, 'edge midpoint latitude')
            CALL gridDefYunits(var_lists(i)%p%cdiEdgeGridID, 'radians')
            !
          END SELECT
        ENDDO

        !
        ! 4. add vertical grid descriptions
        !
        DO ivg = 1, nv_grids

          SELECT CASE (vgrid_def(ivg)%type)

          CASE (ZA_SURFACE)
            var_lists(i)%p%cdiSurfZaxisID = zaxisCreate(ZAXIS_SURFACE, vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(1))
            levels(1) = 0.0_wp
            CALL zaxisDefLevels(var_lists(i)%p%cdiSurfZaxisID, levels)
            DEALLOCATE(levels)

          CASE (ZA_HYBRID)
            var_lists(i)%p%cdiFullZaxisID = zaxisCreate(ZAXIS_HYBRID, vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(vgrid_def(ivg)%nlevels))
            DO k = 1, vgrid_def(ivg)%nlevels
              levels(k) = REAL(k,wp)
            END DO
            CALL zaxisDefLevels(var_lists(i)%p%cdiFullZaxisID, levels)
            DEALLOCATE(levels)
            IF (.NOT. lvct_initialised) CYCLE
            nlevp1 = vgrid_def(ivg)%nlevels+1
            CALL zaxisDefVct(var_lists(i)%p%cdiFullZaxisID, 2*nlevp1, private_vct(1:2*nlevp1))

          CASE (ZA_HYBRID_HALF)
            var_lists(i)%p%cdiHalfZaxisID  = zaxisCreate(ZAXIS_HYBRID_HALF, vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(vgrid_def(ivg)%nlevels))
            DO k = 1, vgrid_def(ivg)%nlevels
              levels(k) = REAL(k,wp)
            END DO
            CALL zaxisDefLevels(var_lists(i)%p%cdiHalfZaxisID, levels)
            DEALLOCATE(levels)
            IF (.NOT. lvct_initialised) CYCLE
            nlevp1 = vgrid_def(ivg)%nlevels
            CALL zaxisDefVct(var_lists(i)%p%cdiHalfZaxisID, 2*nlevp1, private_vct(1:2*nlevp1))

          CASE (ZA_DEPTH_BELOW_LAND)
            IF (.NOT. ldepth_lnd_initialised) CYCLE
            var_lists(i)%p%cdiDepthFullZaxisID = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, &
                 &                                            vgrid_def(ivg)%nlevels)
            CALL zaxisDefLevels(var_lists(i)%p%cdiDepthFullZaxisID, &
                 &              private_depth_lnd_full)

          CASE (ZA_DEPTH_BELOW_LAND_P1)
            IF (.NOT. ldepth_lnd_initialised) CYCLE
            var_lists(i)%p%cdiDepthHalfZaxisID = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, &
                 &                                            vgrid_def(ivg)%nlevels)
            CALL zaxisDefLevels(var_lists(i)%p%cdiDepthHalfZaxisID, &
                 &              private_depth_lnd_half)

          CASE (ZA_DEPTH_RUNOFF_S)
            var_lists(i)%p%cdiDepthRunoff_sZaxisID = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(1))
            levels(1) = 1.0_wp
            CALL zaxisDefLevels(var_lists(i)%p%cdiDepthRunoff_sZaxisID, levels)
            DEALLOCATE(levels)

          CASE (ZA_DEPTH_RUNOFF_G)
            var_lists(i)%p%cdiDepthRunoff_gZaxisID = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(1))
            levels(1) = 1.0_wp
            CALL zaxisDefLevels(var_lists(i)%p%cdiDepthRunoff_gZaxisID, levels)
            DEALLOCATE(levels)

          CASE (ZA_SNOW)
            IF (.NOT. lheight_snow_initialised) CYCLE
            var_lists(i)%p%cdiSnowGenericZaxisID = zaxisCreate(ZAXIS_GENERIC, &
              &                                            vgrid_def(ivg)%nlevels)
            CALL zaxisDefLevels(var_lists(i)%p%cdiSnowGenericZaxisID, &
              &              private_height_snow_full)

          CASE (ZA_SNOW_HALF)
            IF (.NOT. lheight_snow_initialised) CYCLE
            var_lists(i)%p%cdiSnowHalfGenericZaxisID = zaxisCreate(ZAXIS_GENERIC, &
              &                                            vgrid_def(ivg)%nlevels)
            CALL zaxisDefLevels(var_lists(i)%p%cdiSnowHalfGenericZaxisID, &
              &              private_height_snow_half)

          CASE (ZA_TOA)
            var_lists(i)%p%cdiToaZaxisID = zaxisCreate(ZAXIS_TOA, vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(1))
            levels(1) = 1.0_wp
            CALL zaxisDefLevels(var_lists(i)%p%cdiToaZaxisID, levels)
            DEALLOCATE(levels)

          CASE (ZA_HEIGHT_2M)
            var_lists(i)%p%cdiH2mZaxisID = zaxisCreate(ZAXIS_HEIGHT, vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(1))
            levels(1) = 2.0_wp
            CALL zaxisDefLevels(var_lists(i)%p%cdiH2mZaxisID, levels)
            DEALLOCATE(levels)

          CASE (ZA_HEIGHT_10M)
            var_lists(i)%p%cdiH10mZaxisID = zaxisCreate(ZAXIS_HEIGHT, vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(1))
            levels(1) = 10.0_wp
            CALL zaxisDefLevels(var_lists(i)%p%cdiH10mZaxisID, levels)
            DEALLOCATE(levels)

          CASE (ZA_LAKE_BOTTOM)
            var_lists(i)%p%cdiLakeBottomZaxisID = zaxisCreate(ZAXIS_LAKE_BOTTOM, vgrid_def(ivg)%nlevels)
            ALLOCATE(lbounds(1), ubounds(1), levels(1))
            lbounds(1)= 1._wp
            ubounds(1)= 0._wp
            levels(1) = 1._wp
            CALL zaxisDefLbounds(var_lists(i)%p%cdiLakeBottomZaxisID, lbounds) !necessary for GRIB2
            CALL zaxisDefUbounds(var_lists(i)%p%cdiLakeBottomZaxisID, ubounds) !necessary for GRIB2
            CALL zaxisDefLevels (var_lists(i)%p%cdiLakeBottomZaxisID, levels)
            CALL zaxisDefUnits  (var_lists(i)%p%cdiLakeBottomZaxisID, "m")
            DEALLOCATE(lbounds, ubounds, levels)

          CASE (ZA_LAKE_BOTTOM_HALF)
            var_lists(i)%p%cdiLakeHalfBottomZaxisID = zaxisCreate(ZAXIS_LAKE_BOTTOM, vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(1))
            levels(1) = 1._wp
            CALL zaxisDefLevels(var_lists(i)%p%cdiLakeHalfBottomZaxisID, levels)
            CALL zaxisDefUnits(var_lists(i)%p%cdiLakeHalfBottomZaxisID, "m")
            DEALLOCATE(levels)

          CASE (ZA_MIX_LAYER)
            var_lists(i)%p%cdiLakeMixLayerZaxisID = zaxisCreate(ZAXIS_MIX_LAYER, vgrid_def(ivg)%nlevels)
            ALLOCATE(lbounds(1), ubounds(1), levels(1))
            lbounds(1)= 1._wp
            ubounds(1)= 0._wp
            levels(1) = 1._wp
            CALL zaxisDefLbounds(var_lists(i)%p%cdiLakeMixLayerZaxisID, lbounds) !necessary for GRIB2
            CALL zaxisDefUbounds(var_lists(i)%p%cdiLakeMixLayerZaxisID, ubounds) !necessary for GRIB2
            CALL zaxisDefLevels (var_lists(i)%p%cdiLakeMixLayerZaxisID, levels)
            CALL zaxisDefUnits  (var_lists(i)%p%cdiLakeMixLayerZaxisID, "m")
            DEALLOCATE(lbounds, ubounds, levels)

          CASE (ZA_SEDIMENT_BOTTOM_TW_HALF)
            var_lists(i)%p%cdiLakeHalfSedBottomTwZaxisID = zaxisCreate(ZAXIS_SEDIMENT_BOTTOM_TW, &
              &                                            vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(1))
            levels(1) = 0._wp
            CALL zaxisDefLevels(var_lists(i)%p%cdiLakeHalfSedBottomTwZaxisID, levels)
            CALL zaxisDefUnits(var_lists(i)%p%cdiLakeHalfSedBottomTwZaxisID, "m")
            DEALLOCATE(levels)

          !
          ! Ocean
          !
          CASE (ZA_DEPTH_BELOW_SEA)
            IF (.NOT. ldepth_initialised) CYCLE
            var_lists(i)%p%cdiDepthFullZaxisID = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, &
              &                                           vgrid_def(ivg)%nlevels)
            CALL zaxisDefLevels(var_lists(i)%p%cdiDepthFullZaxisID, &
              &              private_depth_full)

          CASE (ZA_DEPTH_BELOW_SEA_HALF)
            IF (.NOT. ldepth_initialised) CYCLE
            var_lists(i)%p%cdiDepthHalfZaxisID = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, &
              &                                           vgrid_def(ivg)%nlevels)
            CALL zaxisDefLevels(var_lists(i)%p%cdiDepthHalfZaxisID, &
              &              private_depth_half)

          CASE (ZA_GENERIC_ICE)
            !!!!!!!!! ATTENTION: !!!!!!!!!!!
            ! As soon as i_no_ice_thick_class is set to i_no_ice_thick_class>1 this no longer works
            var_lists(i)%p%cdiIceGenericZaxisID = zaxisCreate(ZAXIS_GENERIC, &
              &                                           vgrid_def(ivg)%nlevels)
            ALLOCATE(levels(1))
            levels(1) = 1.0_wp
            CALL zaxisDefLevels(var_lists(i)%p%cdiIceGenericZaxisID, levels)
            DEALLOCATE(levels)

          CASE DEFAULT
            CALL finish('open_writing_restart_files','Vertical grid description not found.')
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
      IF (my_process_is_mpi_workroot()) THEN
        !
        CALL addVarListToVlist(var_lists(i), var_lists(i)%p%cdiVlistID)
        !
        WRITE(message_text,'(t1,a49,t50,a31,t84,i6,t94,l5)')        &
             restart_filename, var_lists(i)%p%name,             &
             var_lists(i)%p%cdiFileID_restart, var_lists(i)%p%lrestart
        CALL message('',message_text)
      ENDIF


      !
      ! loop over all other output var_lists eventually corresponding to the same file
      !
      DO j = 1, nvar_lists
        !
        IF (var_lists(j)%p%restart_opened) CYCLE
        IF (.NOT. var_lists(j)%p%lrestart) CYCLE
        !
        IF (var_lists(j)%p%restart_type /= var_lists(i)%p%restart_type) THEN
          CALL finish('open_output_streams', 'different file types for the same restart file')
        ENDIF
        !
        IF (var_lists(i)%p%model_type == var_lists(j)%p%model_type) THEN
          var_lists(j)%p%restart_opened = .TRUE.
          var_lists(j)%p%filename = var_lists(i)%p%filename
          !
          ! set file IDs of all associated restart files
          !
          var_lists(j)%p%cdiFileID_restart       = var_lists(i)%p%cdiFileID_restart
          var_lists(j)%p%cdiVlistID              = var_lists(i)%p%cdiVlistID
          var_lists(j)%p%cdiCellGridID           = var_lists(i)%p%cdiCellGridID
          var_lists(j)%p%cdiVertGridID           = var_lists(i)%p%cdiVertGridID
          var_lists(j)%p%cdiEdgeGridID           = var_lists(i)%p%cdiEdgeGridID
          var_lists(j)%p%cdiSurfZaxisID          = var_lists(i)%p%cdiSurfZaxisID
          var_lists(j)%p%cdiGenericZaxisID       = var_lists(i)%p%cdiGenericZaxisID
          var_lists(j)%p%cdiFullZaxisID          = var_lists(i)%p%cdiFullZaxisID
          var_lists(j)%p%cdiHalfZaxisID          = var_lists(i)%p%cdiHalfZaxisID
          var_lists(j)%p%cdiDepthFullZaxisID     = var_lists(i)%p%cdiDepthFullZaxisID
          var_lists(j)%p%cdiDepthHalfZaxisID     = var_lists(i)%p%cdiDepthHalfZaxisID
          var_lists(j)%p%cdiDepthRunoff_sZaxisID = var_lists(i)%p%cdiDepthRunoff_sZaxisID
          var_lists(j)%p%cdiDepthRunoff_gZaxisID = var_lists(i)%p%cdiDepthRunoff_gZaxisID
          var_lists(j)%p%cdiIceGenericZaxisID    = var_lists(i)%p%cdiIceGenericZaxisID
          var_lists(j)%p%cdiSnowGenericZaxisID   = var_lists(i)%p%cdiSnowGenericZaxisID
          var_lists(j)%p%cdiSnowHalfGenericZaxisID = var_lists(i)%p%cdiSnowHalfGenericZaxisID
          var_lists(j)%p%cdiToaZaxisID           = var_lists(i)%p%cdiToaZaxisID
          var_lists(j)%p%cdiH2mZaxisID           = var_lists(i)%p%cdiH2mZaxisID
          var_lists(j)%p%cdiH10mZaxisID          = var_lists(i)%p%cdiH10mZaxisID
          var_lists(j)%p%cdiLakeBottomZaxisID    = var_lists(i)%p%cdiLakeBottomZaxisID
          var_lists(j)%p%cdiLakeHalfBottomZaxisID= var_lists(i)%p%cdiLakeHalfBottomZaxisID
          var_lists(j)%p%cdiLakeMixLayerZaxisID  = var_lists(i)%p%cdiLakeMixLayerZaxisID
          var_lists(j)%p%cdiLakeHalfSedBottomTwZaxisID  = var_lists(i)%p%cdiLakeHalfSedBottomTwZaxisID
          var_lists(j)%p%cdiTaxisID              = var_lists(i)%p%cdiTaxisID
          !
          ! add variables to already existing cdi vlists
          !
          IF (my_process_is_mpi_workroot()) THEN
            !
            CALL addVarListToVlist(var_lists(j), var_lists(j)%p%cdiVlistID)
            !
            WRITE(message_text,'(t1,a49,t50,a31,t84,i6,t94,l5)')        &
                 restart_filename, var_lists(j)%p%name,             &
                 var_lists(j)%p%cdiFileID_restart, var_lists(j)%p%lrestart
            CALL message('',message_text)
          ENDIF
        ENDIF
      ENDDO

      !
      IF (my_process_is_mpi_workroot() .AND. var_lists(i)%p%first) THEN
        CALL streamDefVlist(var_lists(i)%p%cdiFileID_restart, var_lists(i)%p%cdiVlistID)
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
    INTEGER :: time_level
    INTEGER :: tlev_skip

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
      ! skip this field because of wrong time index ?
      !
      ! get time index of current field
      time_level = get_var_timelevel(element%field)

#ifndef __NO_ICON_ATMO__
! this should be removed !
      ! get information about timelevel to be skipped for current field
      ! for the time being this will work with the global patch only
      IF (element%field%info%tlev_source == 0) THEN
        tlev_skip = nnew(1)          ! ATTENTION: 1 (global patch) hardcoded
      ELSE IF (element%field%info%tlev_source == 1) THEN
        tlev_skip = nnew_rcf(1)      ! ATTENTION: 1 (global patch) hardcoded
      ELSE
        tlev_skip = -99
      ENDIF


      SELECT CASE (iequations)
      CASE(IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER)

        IF ( time_level == tlev_skip                        &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_EXPL &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_SI   ) CYCLE   ! skip field
      CASE default
        IF ( time_level == tlev_skip ) CYCLE   ! skip field
      END SELECT
!  #else
      ! this is not clear if needed
      IF ( time_level == tlev_skip ) CYCLE   ! skip field
#endif

      !
      ! set grid ID
      !
      SELECT CASE (info%hgrid)
      CASE(GRID_UNSTRUCTURED_CELL)
        info%cdiGridID = this_list%p%cdiCellGridID
      CASE(GRID_UNSTRUCTURED_VERT)
        info%cdiGridID = this_list%p%cdiVertGridID
      CASE(GRID_UNSTRUCTURED_EDGE)
        info%cdiGridID = this_list%p%cdiEdgeGridID
      END SELECT

      !
      ! set z axis ID
      !
      SELECT CASE (info%vgrid)
      CASE (ZA_SURFACE)
        info%cdiZaxisID =  this_list%p%cdiSurfZaxisID
      CASE (ZA_HYBRID)
        info%cdiZaxisID =  this_list%p%cdiFullZaxisID
      CASE (ZA_HYBRID_HALF)
        info%cdiZaxisID =  this_list%p%cdiHalfZaxisID
      CASE (ZA_DEPTH_BELOW_LAND)
        info%cdiZaxisID =  this_list%p%cdiDepthFullZaxisID
      CASE (ZA_DEPTH_BELOW_LAND_P1)
        info%cdiZaxisID =  this_list%p%cdiDepthHalfZaxisID
      CASE (ZA_DEPTH_RUNOFF_S)
        info%cdiZaxisID =  this_list%p%cdiDepthRunoff_sZaxisID
      CASE (ZA_DEPTH_RUNOFF_G)
        info%cdiZaxisID =  this_list%p%cdiDepthRunoff_gZaxisID
      CASE (ZA_SNOW)
        info%cdiZaxisID =  this_list%p%cdiSnowGenericZaxisID
      CASE (ZA_SNOW_HALF)
        info%cdiZaxisID =  this_list%p%cdiSnowHalfGenericZaxisID
      CASE (ZA_TOA)
        info%cdiZaxisID =  this_list%p%cdiToaZaxisID
      CASE (ZA_HEIGHT_2M)
        info%cdiZaxisID =  this_list%p%cdiH2mZaxisID
      CASE (ZA_HEIGHT_10M)
        info%cdiZaxisID =  this_list%p%cdiH10mZaxisID
      CASE (ZA_LAKE_BOTTOM)
        info%cdiZaxisID =  this_list%p%cdiLakeBottomZaxisID
      CASE (ZA_LAKE_BOTTOM_HALF)
        info%cdiZaxisID =  this_list%p%cdiLakeHalfBottomZaxisID
      CASE (ZA_MIX_LAYER)
        info%cdiZaxisID =  this_list%p%cdiLakeMixLayerZaxisID
      CASE (ZA_SEDIMENT_BOTTOM_TW_HALF)
        info%cdiZaxisID =  this_list%p%cdiLakeHalfSedBottomTwZaxisID
      !
      ! ocean
      !
      CASE (ZA_DEPTH_BELOW_SEA)
        info%cdiZaxisID =  this_list%p%cdiDepthFullZaxisID
      CASE (ZA_DEPTH_BELOW_SEA_HALF)
        info%cdiZaxisID =  this_list%p%cdiDepthHalfZaxisID
      CASE (ZA_GENERIC_ICE)
        info%cdiZaxisID =  this_list%p%cdiIceGenericZaxisID
      END SELECT

      gridID  = info%cdiGridID
      zaxisID = info%cdiZaxisID


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

  !-------------
  !>
  !!
  !! Hui Wan (MPI-M, 2011-05)
  !!
  SUBROUTINE create_restart_file( patch, datetime,             &
                                & jstep,                       &
                                & opt_pvct,                    &
                                & opt_t_elapsed_phy,           &
                                & opt_lcall_phy, opt_sim_time, &
                                & opt_jstep_adv_ntsteps,       &
                                & opt_jstep_adv_marchuk_order, &
                                & opt_depth, opt_depth_lnd,    &
                                & opt_nlev_snow,               &
                                & opt_nice_class)

    TYPE(t_patch),   INTENT(IN) :: patch
    TYPE(t_datetime),INTENT(IN) :: datetime
    INTEGER, INTENT(IN) :: jstep                ! simulation step

    REAL(wp), INTENT(IN), OPTIONAL :: opt_pvct(:)
    INTEGER,  INTENT(IN), OPTIONAL :: opt_depth
    INTEGER,  INTENT(IN), OPTIONAL :: opt_depth_lnd   ! vertical levels soil model
    REAL(wp), INTENT(IN), OPTIONAL :: opt_t_elapsed_phy(:,:)
    LOGICAL , INTENT(IN), OPTIONAL :: opt_lcall_phy(:,:)
    REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time
    INTEGER,  INTENT(IN), OPTIONAL :: opt_jstep_adv_ntsteps
    INTEGER,  INTENT(IN), OPTIONAL :: opt_jstep_adv_marchuk_order
    INTEGER,  INTENT(IN), OPTIONAL :: opt_nlev_snow
    INTEGER,  INTENT(IN), OPTIONAL :: opt_nice_class

    INTEGER :: klev, jg, kcell, kvert, kedge, icelltype
    INTEGER :: izlev, inlev_soil, inlev_snow, i, nice_class
    REAL(wp), ALLOCATABLE :: zlevels_full(:), zlevels_half(:)


    CHARACTER(LEN=132) :: string
    CHARACTER(len=MAX_CHAR_LENGTH) :: attname   ! attribute name
    INTEGER :: jp, jp_end   ! loop index and array size


    IF (ltimer) CALL timer_start(timer_write_restart_file)
    !----------------
    ! Initialization
    klev      = patch%nlev
    jg        = patch%id
    kcell     = patch%n_patch_cells_g
    kvert     = patch%n_patch_verts_g
    kedge     = patch%n_patch_edges_g
    icelltype = patch%cell_type

    CALL set_restart_attribute( 'current_caltime', datetime%caltime )
    CALL set_restart_attribute( 'current_calday' , datetime%calday )

    CALL set_restart_attribute( 'current_daysec' , datetime%daysec )

    CALL set_restart_attribute( 'nold'    , nold    (jg))
    CALL set_restart_attribute( 'nnow'    , nnow    (jg))
    CALL set_restart_attribute( 'nnew'    , nnew    (jg))
    CALL set_restart_attribute( 'nnow_rcf', nnow_rcf(jg))
    CALL set_restart_attribute( 'nnew_rcf', nnew_rcf(jg))

    ! set simulation step
    CALL set_restart_attribute( 'jstep', jstep )

    !----------------
    ! additional restart-output for nonhydrostatic model
    IF (PRESENT(opt_sim_time)) THEN
      WRITE(attname,'(a,i2.2)') 'sim_time_DOM',jg
      CALL set_restart_attribute( TRIM(attname), opt_sim_time )
    ENDIF

    !-------------------------------------------------------------
    ! DR
    ! WORKAROUND FOR FIELDS WHICH NEED TO GO INTO THE RESTART FILE,
    ! BUT SO FAR CANNOT BE HANDELED CORRECTLY BY ADD_VAR OR
    ! SET_RESTART_ATTRIBUTE
    !-------------------------------------------------------------

    IF (PRESENT(opt_jstep_adv_ntsteps)) THEN
        WRITE(attname,'(a,i2.2)') 'jstep_adv_ntsteps_DOM',jg
        CALL set_restart_attribute( TRIM(attname), opt_jstep_adv_ntsteps )
    ENDIF

    IF (PRESENT(opt_jstep_adv_marchuk_order)) THEN
        WRITE(attname,'(a,i2.2)') 'jstep_adv_marchuk_order_DOM',jg
        CALL set_restart_attribute( TRIM(attname), opt_jstep_adv_marchuk_order )
    ENDIF

    IF (PRESENT(opt_t_elapsed_phy) .AND. PRESENT(opt_lcall_phy)) THEN
      ! Inquire array size
      jp_end = SIZE(opt_t_elapsed_phy,2)
      DO jp = 1, jp_end
        WRITE(attname,'(a,i2.2,a,i2.2)') 't_elapsed_phy_DOM',jg,'_PHY',jp
        CALL set_restart_attribute( TRIM(attname), opt_t_elapsed_phy(jg,jp) )
      ENDDO
      ! Inquire array size
      jp_end = SIZE(opt_lcall_phy,2)
      DO jp = 1, jp_end
        WRITE(attname,'(a,i2.2,a,i2.2)') 'lcall_phy_DOM',jg,'_PHY',jp
        CALL set_restart_attribute( TRIM(attname), opt_lcall_phy(jg,jp) )
      ENDDO
    ENDIF

    IF (PRESENT(opt_pvct)) CALL set_restart_vct( opt_pvct )  ! Vertical coordinate (A's and B's)
    IF (PRESENT(opt_depth_lnd)) THEN            ! geometrical depth for land module
      !This part is only called if opt_depth_lnd > 0
      IF (opt_depth_lnd > 0) THEN
        inlev_soil = opt_depth_lnd
        ALLOCATE(zlevels_full(inlev_soil))
        ALLOCATE(zlevels_half(inlev_soil+1))
        DO i = 1, inlev_soil
          zlevels_full(i) = REAL(i,wp)
        END DO
        DO i = 1, inlev_soil+1
          zlevels_half(i) = REAL(i,wp)
        END DO
        CALL set_restart_depth_lnd(zlevels_half, zlevels_full)
        DEALLOCATE(zlevels_full)
        DEALLOCATE(zlevels_half)
      ELSE
       inlev_soil = 0
      END IF
    ELSE
      inlev_soil = 0
    ENDIF
    IF (PRESENT(opt_nlev_snow)) THEN  ! number of snow levels (multi layer snow model)
      !This part is only called if opt_nlev_snow > 0
      IF (opt_nlev_snow > 0) THEN
        inlev_snow = opt_nlev_snow
        ALLOCATE(zlevels_full(inlev_snow))
        ALLOCATE(zlevels_half(inlev_snow+1))
        DO i = 1, inlev_snow
          zlevels_full(i) = REAL(i,wp)
        END DO
        DO i = 1, inlev_snow+1
          zlevels_half(i) = REAL(i,wp)
        END DO
        CALL set_restart_height_snow(zlevels_half, zlevels_full)
        DEALLOCATE(zlevels_full)
        DEALLOCATE(zlevels_half)
      ELSE
        inlev_snow = 0
      ENDIF
    ELSE
      inlev_snow = 0
    ENDIF
!DR end preliminary fix
    izlev = 0
#ifndef __NO_ICON_OCEAN__
    IF (PRESENT(opt_depth)) THEN                              ! Ocean depth
      !This part is only called if opt_depth > 0
      IF(opt_depth>0)THEN
        izlev = opt_depth
        ALLOCATE(zlevels_full(izlev))
        ALLOCATE(zlevels_half(izlev+1))
        CALL set_zlev(zlevels_half, zlevels_full, n_zlev, dzlev_m)
        CALL set_restart_depth(zlevels_half, zlevels_full)
        DEALLOCATE(zlevels_full)
        DEALLOCATE(zlevels_half)
!      ELSE
!        izlev = 0
      END IF
!    ELSE
!      izlev = 0
    END IF
#endif

    IF (.NOT.PRESENT(opt_nice_class)) THEN
      nice_class = 1
    ELSE
      nice_class = opt_nice_class
    END IF

    CALL init_restart( 'ICON',            &! model name
                     & '1.4.00',          &! model version
                     & kcell, icelltype,  &! total # of cells, # of vertices per cell
                     & kvert, 9-icelltype,&! total # of vertices, # of vertices per dual cell
                     & kedge, 4,          &! total # of cells, shape of control volume for edge
                     & klev,              &! total # of vertical layers
                     & izlev,             &! total # of depths below sea
                     & inlev_soil,        &! total # of depths below land (TERRA or JSBACH)
                     & inlev_snow,        &! total # of vertical snow layers (TERRA)
                     & nice_class         )! total # of ice classes (sea ice)

    CALL set_restart_time( iso8601(datetime) )  ! Time tag

    ! Open new file, write data, close and then clean-up.
    message_text = get_filename_noext(patch%grid_filename)
    WRITE(string,'(a,a)') 'restart.',TRIM(message_text)

    CALL open_writing_restart_files( TRIM(string) )

    CALL write_restart( patch )

    CALL close_writing_restart_files
    CALL finish_restart

    IF (ltimer) CALL timer_stop(timer_write_restart_file)

  END SUBROUTINE create_restart_file

  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE close_writing_restart_files
    !
    ! Loop over all the output streams and close the associated files, set
    ! opened to false
    !
    CHARACTER(len=80) :: linkname
    INTEGER :: i, j, iret, fileID, vlistID

    IF (my_process_is_mpi_test()) RETURN
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
      IF (var_lists(i)%p%restart_opened) THEN
        IF (my_process_is_mpi_workroot() .AND. var_lists(i)%p%first) THEN
          !
          fileID = var_lists(i)%p%cdiFileID_restart
          !
          IF (fileID /= CDI_UNDEFID) THEN
            CALL streamClose(var_lists(i)%p%cdiFileID_restart)
            DO j = 1, nvar_lists
              IF (fileID == var_lists(j)%p%cdiFileID_restart) THEN
                var_lists(j)%p%cdiFileID_restart = CDI_UNDEFID
              ENDIF
            ENDDO
          ENDIF
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
      var_lists(i)%p%restart_opened = .FALSE.
      var_lists(i)%p%first          = .FALSE.
    ENDDO reset_all_lists
    !
    private_restart_time = ''
    !
  END SUBROUTINE close_writing_restart_files

  !------------------------------------------------------------------------------------------------
  !
  ! loop over all var_lists for restart
  !
  SUBROUTINE write_restart(p_patch)
    TYPE(t_patch), INTENT(in) :: p_patch

    INTEGER :: i,j
    LOGICAL :: write_info
    !
    IF (my_process_is_mpi_test()) RETURN

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
              CALL message('', &
                   'Write compressed netCDF4 restart for : '//TRIM(private_restart_time))
            ELSE
              CALL message('','Write netCDF4 restart for : '//TRIM(private_restart_time))
            END IF
          END SELECT
        ENDIF
        write_info = .FALSE.
        !
        ! write time information to netCDF file
        !
        IF (my_process_is_mpi_workroot()) THEN
          CALL write_time_to_restart(var_lists(i))
        ENDIF
        !
        ! loop over all streams associated with the file
        !
        DO j = i, nvar_lists

          IF (var_lists(j)%p%cdiFileID_restart == var_lists(i)%p%cdiFileID_restart) THEN
            !
            !
            ! write variables
            !
            CALL write_restart_var_list(var_lists(j), p_patch=p_patch)
            !
          ENDIF
        ENDDO
      ENDIF
    ENDDO
!PR
  CALL message('','Finished Write netCDF2 restart for : '//TRIM(private_restart_time))
    !
  END SUBROUTINE write_restart

  !------------------------------------------------------------------------------------------------
  !
  ! set time for restart in cdi format
  !
  SUBROUTINE write_time_to_restart (this_list)
    TYPE (t_var_list), INTENT(inout) :: this_list
    !
    INTEGER :: fileID, idate, itime, iret
    !
    fileID = this_list%p%cdiFileID_restart
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
  SUBROUTINE write_restart_var_list(this_list, p_patch)
    TYPE (t_var_list) ,INTENT(in) :: this_list
    TYPE(t_patch), TARGET, INTENT(in) :: p_patch
    !
    ! variables of derived type used in linked list
    !
    TYPE (t_var_metadata), POINTER :: info
    TYPE (t_list_element), POINTER :: element
    TYPE (t_list_element), TARGET  :: start_with
    !
    REAL(wp), ALLOCATABLE :: r2d(:,:) ! field gathered on I/O processor
    REAL(wp), ALLOCATABLE :: r1d(:) ! field gathered on I/O processor
    !
    INTEGER :: nindex, private_n
    TYPE(t_comm_gather_pattern), POINTER :: gather_pattern
    !
    INTEGER :: time_level
    INTEGER :: tlev_skip

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
      ! retrieve information from actual linked list element
      !
      info => element%field%info

      !
      ! skip this field ?
      !
      IF (.NOT. info%lrestart) CYCLE

      !
      ! skip this field because of wrong time index ?
      !
      ! get time index of current field
      time_level = get_var_timelevel(element%field)

#ifndef __NO_ICON_ATMO__
! this should be removed !
      ! get information about timelevel to be skipped for current field
      ! for the time being this will work with the global patch only
      IF (element%field%info%tlev_source == 0) THEN
        tlev_skip = nnew(1)          ! ATTENTION: 1 (global patch) hardcoded
      ELSE IF (element%field%info%tlev_source == 1) THEN
        tlev_skip = nnew_rcf(1)      ! ATTENTION: 1 (global patch) hardcoded
      ELSE
        tlev_skip = -99
      ENDIF

      SELECT CASE (iequations)
      CASE(IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER)

        IF ( time_level == tlev_skip                        &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_EXPL &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_SI   ) CYCLE   ! skip field
      CASE default
        IF ( time_level == tlev_skip ) CYCLE   ! skip field
      END SELECT
! #else
     IF ( time_level == tlev_skip ) CYCLE   ! skip field
#endif

      nindex = MERGE(info%ncontained, 1, info%lcontained)

      SELECT CASE (info%hgrid)
      CASE (GRID_UNSTRUCTURED_CELL)
        private_n = private_nc
        gather_pattern => p_patch%comm_pat_gather_c
      CASE (GRID_UNSTRUCTURED_VERT)
        private_n = private_nv
        gather_pattern => p_patch%comm_pat_gather_v
      CASE (GRID_UNSTRUCTURED_EDGE)
        private_n = private_ne
        gather_pattern => p_patch%comm_pat_gather_e
      CASE default
        CALL finish('out_stream','unknown grid type')
      END SELECT

      SELECT CASE (info%ndims)
      CASE (1)
        CALL finish('write_restart_var_list','1d arrays not handled yet.')
      CASE (2)
        ALLOCATE(r1d(MERGE(private_n, 0, my_process_is_mpi_workroot())))
        CALL exchange_data(element%field%r_ptr(:,:,nindex,1,1), r1d, &
          &                gather_pattern)
      CASE (3)
        ALLOCATE(r2d(MERGE(private_n, 0, my_process_is_mpi_workroot()), &
          &          info%used_dimensions(2)))
        CALL exchange_data(element%field%r_ptr(:,:,:,nindex,1), r2d, &
          &                gather_pattern)
      CASE (4)
        CALL finish('write_restart_var_list','4d arrays not handled yet.')
      CASE (5)
        CALL finish('write_restart_var_list','5d arrays not handled yet.')
      CASE DEFAULT
        CALL finish('write_restart_var_list','dimension not set.')
      END SELECT
      !
      ! write data
      !
      IF (my_process_is_mpi_workroot()) THEN
        write (0,*)' ... write ',info%name

        IF (info%ndims == 2) THEN
          CALL streamWriteVar(this_list%p%cdiFileID_restart, info%cdiVarID, &
            &                 r1d(:), 0)
        ELSE IF (info%ndims == 3) THEN
          CALL streamWriteVar(this_list%p%cdiFileID_restart, info%cdiVarID, &
            &                 r2d(:,:), 0)
        END IF
      END IF
      !
      ! deallocate temporary global arrays
      !
      IF (ALLOCATED(r1d)) DEALLOCATE(r1d)
      IF (ALLOCATED(r2d)) DEALLOCATE(r2d)
      !
    END DO for_all_list_elements
    !
  END SUBROUTINE write_restart_var_list

  !------------------------------------------------------------------------------------------------
  !
  ! deallocate module variables
  !
  SUBROUTINE finish_restart

    INTEGER :: i

    IF (my_process_is_mpi_test()) RETURN

    for_all_var_lists: DO i = 1, nvar_lists
      IF (var_lists(i)%p%cdiFileID_restart >= 0) THEN
        CALL gridDestroy(var_lists(i)%p%cdiCellGridID)
        CALL gridDestroy(var_lists(i)%p%cdiVertGridID)
        CALL gridDestroy(var_lists(i)%p%cdiEdgeGridID)
        IF (var_lists(i)%p%cdiSurfZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiSurfZaxisID)
        IF (var_lists(i)%p%cdiGenericZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiGenericZaxisID)
        IF (var_lists(i)%p%cdiFullZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiFullZaxisID)
        IF (var_lists(i)%p%cdiHalfZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiHalfZaxisID)
        IF (var_lists(i)%p%cdiDepthFullZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiDepthFullZaxisID)
        IF (var_lists(i)%p%cdiDepthHalfZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiDepthHalfZaxisID)
        IF (var_lists(i)%p%cdiDepthRunoff_sZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiDepthRunoff_sZaxisID)
        IF (var_lists(i)%p%cdiDepthRunoff_gZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiDepthRunoff_gZaxisID)
        IF (var_lists(i)%p%cdiIceGenericZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiIceGenericZaxisID)
        IF (var_lists(i)%p%cdiH2mZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiH2mZaxisID)
        IF (var_lists(i)%p%cdiH10mZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiH10mZaxisID)
        IF (var_lists(i)%p%cdiLakeBottomZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiLakeBottomZaxisID)
        IF (var_lists(i)%p%cdiLakeHalfBottomZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiLakeHalfBottomZaxisID)
        IF (var_lists(i)%p%cdiLakeMixLayerZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiLakeMixLayerZaxisID)
        IF (var_lists(i)%p%cdiLakeHalfSedBottomTwZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiLakeHalfSedBottomTwZaxisID)
        IF (var_lists(i)%p%cdiSnowGenericZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiSnowGenericZaxisID)
        IF (var_lists(i)%p%cdiSnowHalfGenericZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiSnowHalfGenericZaxisID)
        IF (var_lists(i)%p%cdiToaZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiToaZaxisID)
        var_lists(i)%p%cdiFileId_restart       = CDI_UNDEFID
        var_lists(i)%p%cdiVlistId              = CDI_UNDEFID
        var_lists(i)%p%cdiCellGridID           = CDI_UNDEFID
        var_lists(i)%p%cdiVertGridID           = CDI_UNDEFID
        var_lists(i)%p%cdiEdgeGridID           = CDI_UNDEFID
        var_lists(i)%p%cdiSurfZaxisID          = CDI_UNDEFID
        var_lists(i)%p%cdiGenericZaxisID       = CDI_UNDEFID
        var_lists(i)%p%cdiHalfZaxisID          = CDI_UNDEFID
        var_lists(i)%p%cdiFullZaxisID          = CDI_UNDEFID
        var_lists(i)%p%cdiDepthHalfZaxisID     = CDI_UNDEFID
        var_lists(i)%p%cdiDepthFullZaxisID     = CDI_UNDEFID
        var_lists(i)%p%cdiDepthRunoff_sZaxisID = CDI_UNDEFID
        var_lists(i)%p%cdiDepthRunoff_gZaxisID = CDI_UNDEFID
        var_lists(i)%p%cdiIceGenericZaxisID    = CDI_UNDEFID
        var_lists(i)%p%cdiH2mZaxisID           = CDI_UNDEFID
        var_lists(i)%p%cdiH10mZaxisID          = CDI_UNDEFID
        var_lists(i)%p%cdiLakeBottomZaxisID    = CDI_UNDEFID
        var_lists(i)%p%cdiLakeHalfBottomZaxisID= CDI_UNDEFID
        var_lists(i)%p%cdiLakeMixLayerZaxisID  = CDI_UNDEFID
        var_lists(i)%p%cdiLakeHalfSedBottomTwZaxisID = CDI_UNDEFID
        var_lists(i)%p%cdiToaZaxisID           = CDI_UNDEFID
        var_lists(i)%p%cdiTaxisID              = CDI_UNDEFID
        var_lists(i)%p%cdiTimeIndex            = CDI_UNDEFID
        var_lists(i)%p%cdiSnowGenericZaxisID   = CDI_UNDEFID
        var_lists(i)%p%cdiSnowHalfGenericZaxisID = CDI_UNDEFID
      ENDIF
    ENDDO for_all_var_lists
    !
    IF (ALLOCATED(private_vct)) DEALLOCATE(private_vct)
    lvct_initialised = .FALSE.
    IF (ALLOCATED(private_depth_full)) DEALLOCATE(private_depth_full)
    IF (ALLOCATED(private_depth_half)) DEALLOCATE(private_depth_half)
    ldepth_initialised = .FALSE.
    IF (ALLOCATED(private_depth_lnd_full)) DEALLOCATE(private_depth_lnd_full)
    IF (ALLOCATED(private_depth_lnd_half)) DEALLOCATE(private_depth_lnd_half)
    ldepth_lnd_initialised = .FALSE.
    IF (ALLOCATED(private_height_snow_full)) DEALLOCATE(private_height_snow_full)
    IF (ALLOCATED(private_height_snow_half)) DEALLOCATE(private_height_snow_half)
    lheight_snow_initialised = .FALSE.
    !
    CALL delete_attributes
    !
    nh_grids   = 0
    nv_grids   = 0
    nt_axis    = 0
    lrestart_initialised = .FALSE.
    !
  END SUBROUTINE finish_restart

  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE read_restart_files(p_patch)
    !
    TYPE(t_patch), OPTIONAL, INTENT(in) :: p_patch
    !
    CHARACTER(len=80) :: restart_filename, name
    !
    INTEGER :: fileID, vlistID, gridID, zaxisID, taxisID, varID
    INTEGER :: idate, itime
    INTEGER :: ic, il
    !
    TYPE model_search
      CHARACTER(len=8) :: abbreviation
      INTEGER          :: key
    END type model_search
    TYPE(model_search) :: abbreviations(nvar_lists)
    !
    TYPE (t_list_element), POINTER :: element
    TYPE (t_var_metadata), POINTER :: info
    !
    CHARACTER(len=8) :: model_type
    INTEGER :: n, nfiles, i, iret, istat, key, vgrid
    INTEGER :: gdims(5), nindex, nmiss
    !
    REAL(wp), POINTER :: r5d(:,:,:,:,:)
    !
    REAL(wp), POINTER :: rptr2d(:,:)
    REAL(wp), POINTER :: rptr3d(:,:,:)
    !
    INTEGER :: string_length  !, ncid

    IF (my_process_is_mpi_workroot()) write(0,*) "read_restart_files, nvar_lists=", nvar_lists
    abbreviations(1:nvar_lists)%key = 0
    abbreviations(1:nvar_lists)%abbreviation = ""
    key = 0
    n = 1
    for_all_model_types: DO i = 1, nvar_lists
! --------------------------------------------------------------
     key = util_hashword(TRIM(var_lists(i)%p%model_type), LEN_TRIM(var_lists(i)%p%model_type), 0)
      IF (.NOT. ANY(abbreviations(1:n)%key == key)) THEN
        abbreviations(n)%abbreviation = var_lists(i)%p%model_type
        abbreviations(n)%key = key
        n = n+1
      ENDIF
    ENDDO for_all_model_types
! --------------------------------------------------------------
    nfiles = n-1
!     write(0,*) 'nfiles=', nfiles
    !
!    CALL message('--','--')
!     CALL message('--',separator)
!     CALL message('','')
    !
    for_all_files: DO n = 1, nfiles
      model_type =TRIM(abbreviations(n)%abbreviation)
      restart_filename = 'restart_'//TRIM(model_type)//'.nc'
!       write(0,*) "n=", n
!       write(0,*) "model_type=", model_type
!       write(0,*) "restart_filename=", restart_filename
!       write(0,*) "util_islink(TRIM(restart_filename)=", &
!        util_islink(TRIM(restart_filename))
      !
      IF (.NOT. util_islink(TRIM(restart_filename))) THEN
        iret = util_rename(TRIM(restart_filename), TRIM(restart_filename)//'.bak')
        write(0,*) "util_rename returned:", iret
        iret = util_symlink(TRIM(restart_filename)//'.bak', TRIM(restart_filename))
        write(0,*) "util_symlink:", iret
      ENDIF
      !

      string_length=LEN_TRIM(restart_filename)
      name = TRIM(restart_filename)//CHAR(0)
      ! check if the netcdf open works
!       write(0,*) "nf_open ", TRIM(restart_filename)
!       CALL nf(nf_open(TRIM(restart_filename), nf_nowrite, ncid))
!       CALL nf(nf_close(ncid))

      IF (my_process_is_mpi_workroot()) write(0,*) "streamOpenRead ", TRIM(restart_filename)

      fileID  = streamOpenRead(name)
      IF (my_process_is_mpi_workroot()) write(0,*) "fileID=",fileID
      vlistID = streamInqVlist(fileID)
!       write(0,*) "vlistID=",vlistID

      taxisID = vlistInqTaxis(vlistID)
!       write(0,*) "taxisID=",taxisID
      !
      idate = taxisInqVdate(taxisID)
!       write(0,*) "idate=",idate
      itime = taxisInqVtime(taxisID)
!       write(0,*) "itime=",itime
      !
      WRITE(message_text,'(a,i8.8,a,i6.6,a,a)') &
           'Read restart for : ', idate, 'T', itime, 'Z from ',TRIM(restart_filename)
      CALL message('read_restart_files',message_text)
      !
      CALL read_attributes(vlistID)
      !
      for_all_vars: DO varID = 0, vlistNvars(vlistID)-1
        !
        CALL vlistInqVarName(vlistID, varID, name)
        !
        for_all_lists: DO i = 1, nvar_lists
          IF (var_lists(i)%p%model_type == model_type) THEN
            element => find_list_element(var_lists(i), TRIM(name))
            IF (ASSOCIATED(element)) THEN
              IF (element%field%info%lrestart) THEN
                !
                info => element%field%info
                !
                ! allocate temporary global array on output processor
                ! and gather field from other processors
                !
                !
                NULLIFY(r5d)
                !
                NULLIFY(rptr2d)
                NULLIFY(rptr3d)
                !
                gridID = vlistInqVarGrid(vlistID, varID)
                ic = gridInqSize(gridID)
                zaxisID = vlistInqVarZaxis(vlistID, varID)
                vgrid = zaxisInqType(zaxisID)
                IF (vgrid == ZAXIS_SURFACE) THEN
                  il = 1
                ELSE
                  il = zaxisInqSize(zaxisID)
                ENDIF
                !
                gdims(:) = (/ ic, il, 1, 1, 1 /)
                !
                IF (my_process_is_mpi_workroot()) THEN
                  ALLOCATE(r5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)),STAT=istat)
                  IF (istat /= 0) THEN
                    CALL finish('','allocation of r5d failed ...')
                  ENDIF
                  CALL streamReadVar(fileID, varID, r5d, nmiss)
                ELSE
                  ALLOCATE(r5d(gdims(1),1,gdims(3),gdims(4),gdims(5)),STAT=istat)
                ENDIF
                CALL p_barrier(comm=p_comm_work)
                !
                IF (info%lcontained) THEN
                  nindex = info%ncontained
                ELSE
                  nindex = 1
                ENDIF
                !
                SELECT CASE (info%hgrid)
                CASE (GRID_UNSTRUCTURED_CELL)
                  IF (info%ndims == 2) THEN
                    rptr2d => element%field%r_ptr(:,:,nindex,1,1)
                    CALL scatter_cells(r5d, rptr2d, p_patch=p_patch)
                  ELSE
                    rptr3d => element%field%r_ptr(:,:,:,nindex,1)
                    CALL scatter_cells(r5d, rptr3d, p_patch=p_patch)
                  ENDIF
                CASE (GRID_UNSTRUCTURED_VERT)
                  IF (info%ndims == 2) THEN
                    rptr2d => element%field%r_ptr(:,:,nindex,1,1)
                    CALL scatter_vertices(r5d, rptr2d, p_patch=p_patch)
                  ELSE
                    rptr3d => element%field%r_ptr(:,:,:,nindex,1)
                    CALL scatter_vertices(r5d, rptr3d, p_patch=p_patch)
                  ENDIF
                CASE (GRID_UNSTRUCTURED_EDGE)
                  IF (info%ndims == 2) THEN
                    rptr2d => element%field%r_ptr(:,:,nindex,1,1)
                    CALL scatter_edges(r5d, rptr2d, p_patch=p_patch)
                  ELSE
                    rptr3d => element%field%r_ptr(:,:,:,nindex,1)
                    CALL scatter_edges(r5d, rptr3d, p_patch=p_patch)
                  ENDIF
                CASE default
                  CALL finish('out_stream','unknown grid type')
                END SELECT
                !
                ! deallocate temporary global arrays
                !
                IF (ASSOCIATED (r5d)) DEALLOCATE (r5d)
                !
                IF (my_process_is_mpi_workroot()) THEN
                  write (0,*) ' ... read ',TRIM(element%field%info%name)
                ENDIF
                CYCLE for_all_vars
              ENDIF
            ENDIF
          ENDIF
        ENDDO for_all_lists
        CALL message('reading_restart_file','Variable '//TRIM(name)//' not defined.')
      ENDDO for_all_vars
      !
      CALL streamClose(fileID)
      !
    ENDDO for_all_files
    !
    CALL message('','')
    CALL message('',separator)
    CALL message('','')
    !
  END SUBROUTINE read_restart_files
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE nf(STATUS)
    INTEGER, INTENT(in) :: STATUS

    IF (STATUS /= nf_noerr) THEN
      CALL finish('mo_io_grid netCDF error', nf_strerror(STATUS))
    ENDIF

  END SUBROUTINE nf
  !-------------------------------------------------------------------------
END MODULE mo_io_restart
