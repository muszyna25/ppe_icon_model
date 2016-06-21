!> Module for writing restart files (synchronously) and for reading restart files.
!!
!! Note: The asynchronous implementation of the restart output can be
!!       found in the module "mo_io_restart_async"
!!
!! Initial implementation: L. Kornblueh
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! --------------------------------------------------------------------------------
!!
!! Generated files:
!! ================
!!
!! 1. For each domain 1, ..., n_dom, and for each restart output time step:
!!
!!    Restart data file
!!    -----------------
!!      "<gridfile>_restart_<modeltype>_<timestamp>.nc",     e.g.
!!      "iconR2B06_DOM01_restart_atm_20110101T001200Z.nc"    (NetCDF format)
!!
!!    This filename can be customized using the namelist parameter
!!    "mo_run_nml/restart_filename".
!!
!!    This file contains
!!    -  data
!!    -  namelists
!!    -  several attributes
!!
!!    Note:
!!    -  We read the namelists only once and assume that these
!!       are identical for all domains.
!!    -  Since we do not know about the total number of domains at startup,
!!       we have to ask the current restart file for the attribute "n_dom"
!!
!! 2. For each domain 1, ..., n_dom, and for the LAST restart output time step:
!!
!!    Symbolic link to data file
!!    --------------------------
!!      "restart_<modeltype>_DOMxx.nc"
!!
!!    Note:
!!    -  The domain-dependent suffix "...DOMxx" is also required for non-nested setups.
!!
!! --------------------------------------------------------------------------------
!!
!OPTION! -pvctl conflict
MODULE mo_io_restart
  !
  USE mo_kind,                  ONLY: wp
  USE mo_mpi,                   ONLY: p_barrier,p_comm_work,p_bcast
  USE mo_exception,             ONLY: finish, message, message_text, get_filename_noext
  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH, TLEV_NNOW, TLEV_NNOW_RCF
  USE mo_var_metadata_types,    ONLY: t_var_metadata
  USE mo_linked_list,           ONLY: t_var_list, t_list_element
  USE mo_var_list,              ONLY: nvar_lists, var_lists, get_var_timelevel, find_list_element
  USE mo_cdi,                   ONLY: FILETYPE_NC, FILETYPE_NC2, FILETYPE_NC4, ZAXIS_SURFACE, CDI_UNDEFID, COMPRESS_ZIP, &
                                    & DATATYPE_FLT64, TIME_VARIABLE, CDI_GLOBAL, DATATYPE_INT32, GRID_UNSTRUCTURED, &
                                    & TAXIS_ABSOLUTE, ZAXIS_DEPTH_BELOW_LAND, ZAXIS_GENERIC, ZAXIS_DEPTH_BELOW_SEA, ZAXIS_HEIGHT, &
                                    & ZAXIS_HYBRID, ZAXIS_HYBRID_HALF, ZAXIS_LAKE_BOTTOM, ZAXIS_MIX_LAYER, &
                                    & ZAXIS_SEDIMENT_BOTTOM_TW, ZAXIS_TOA, TAXIS_RELATIVE, streamOpenRead, streamInqVlist, &
                                    & vlistInqTaxis, taxisInqVdate, taxisInqVtime, vlistNvars, vlistInqVarGrid, gridInqSize, &
                                    & vlistInqVarZaxis, zaxisInqType, zaxisInqSize, streamDefTimestep, vlistDefVar, zaxisCreate, &
                                    & taxisCreate, gridCreate, vlistDefAttInt, vlistDefAttFlt, vlistDefAttTxt, vlistCreate, &
                                    & streamOpenWrite, zaxisDestroy, gridDestroy, vlistDestroy, streamClose, streamWriteVarSlice, &
                                    & streamWriteVar, streamDefVlist, vlistDefVarDatatype, vlistDefVarName, &
                                    & vlistInqVarName, zaxisDefLevels, gridDefNvertex, streamReadVar, zaxisDefLbounds, &
                                    & zaxisDefUbounds, zaxisDefVct, zaxisDefUnits, vlistDefVarLongname, vlistDefVarUnits, &
                                    & vlistDefVarMissval, gridDefXlongname, gridDefYlongname, vlistDefTaxis, taxisDefVdate, &
                                    & taxisDefVtime, gridDefXname, gridDefYname, gridDefXunits, gridDefYunits
  USE mo_util_cdi,              ONLY: cdiGetStringError
  USE mo_cdi_constants,         ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, &
                                    & ZA_HYBRID, ZA_HYBRID_HALF, ZA_DEPTH_BELOW_LAND, ZA_DEPTH_BELOW_LAND_P1, ZA_DEPTH_RUNOFF_S, &
                                    & ZA_DEPTH_RUNOFF_G, ZA_SNOW, ZA_SNOW_HALF, ZA_TOA, ZA_HEIGHT_2M, ZA_HEIGHT_10M, &
                                    & ZA_LAKE_BOTTOM, ZA_LAKE_BOTTOM_HALF, ZA_MIX_LAYER, ZA_SEDIMENT_BOTTOM_TW_HALF, &
                                    & ZA_DEPTH_BELOW_SEA, ZA_DEPTH_BELOW_SEA_HALF, ZA_GENERIC_ICE, ZA_OCEAN_SEDIMENT
  USE mo_cf_convention
  USE mo_util_string,           ONLY: t_keyword_list, associate_keyword, with_keywords, &
    &                                 int2string, separator
  USE mo_util_file,             ONLY: util_symlink, util_rename, util_islink, util_unlink
  USE mo_util_hash,             ONLY: util_hashword
  USE mo_util_uuid,             ONLY: t_uuid
  USE mo_io_restart_namelist,   ONLY: nmls, restart_namelist,                       &
    &                                 read_and_bcast_restart_namelists
  USE mo_io_restart_attributes, ONLY: set_restart_attribute, get_restart_attribute, &
    &                                 delete_attributes,                            &
    &                                 restart_attributes_count_text,                &
    &                                 restart_attributes_count_real,                &
    &                                 restart_attributes_count_int,                 &
    &                                 restart_attributes_count_bool,                &
    &                                 read_and_bcast_attributes
  USE mo_datetime,              ONLY: t_datetime,iso8601,iso8601extended
  USE mo_run_config,            ONLY: ltimer, restart_filename
  USE mo_timer,                 ONLY: timer_start, timer_stop,                      &
    &                                 timer_write_restart_file

  USE mo_dynamics_config,       ONLY: iequations, nold, nnow, nnew, nnew_rcf, nnow_rcf
  USE mo_grid_config,           ONLY: l_limited_area

#ifndef __NO_ICON_ATMO__
!LK comment: should not be here !!!!!! polution of namespace !!!!!!
!GZ: but then we need an alternative method of skipping unnecessary time levels!
  USE mo_impl_constants,        ONLY: IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER, &
    &                                 LEAPFROG_EXPL, LEAPFROG_SI, INH_ATMOSPHERE
  USE mo_ha_dyn_config,         ONLY: ha_dyn_config
#endif

  USE mo_model_domain,          ONLY: t_patch
  USE mo_mpi,                   ONLY: my_process_is_mpi_workroot, &
    & my_process_is_mpi_test
  USE mo_communication,         ONLY: t_comm_gather_pattern, exchange_data, &
    &                                 t_scatterPattern

#ifndef __NO_ICON__OCEAN
  USE mo_ocean_nml,                ONLY: lhamocc
  USE mo_sedmnt,                ONLY: ks, ksp, dzsed
  USE mo_math_utilities,        ONLY: set_zlev
#endif
  !
  IMPLICIT NONE
  !
  PRIVATE

  INCLUDE 'netcdf.inc'
  !
  PUBLIC :: create_restart_file
  PUBLIC :: read_restart_files
  PUBLIC :: read_restart_header
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
  REAL(wp), ALLOCATABLE :: private_height_snow_half(:), private_height_snow_full(:)
  !
  LOGICAL, SAVE :: lvct_initialised         = .FALSE.
  LOGICAL, SAVE :: use_ocean_levels       = .FALSE.
  LOGICAL, SAVE :: ldepth_lnd_initialised   = .FALSE.
  LOGICAL, SAVE :: lheight_snow_initialised = .FALSE.
  !
  LOGICAL, SAVE :: lrestart_initialised     = .FALSE.
  !
  INTEGER :: private_nc  = -1
  INTEGER :: private_nv  = -1
  INTEGER :: private_ne  = -1
  !
  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_io_restart'
  !
  !------------------------------------------------------------------------------------------------
CONTAINS
  !------------------------------------------------------------------------------------------------
  !> Reads attributes and namelists for all available domains from restart file.
  SUBROUTINE read_restart_header(modeltype_str)
    CHARACTER(LEN=*), INTENT(IN)  :: modeltype_str
    ! local variables
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//"::read_restart_header"
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: rst_filename
    CHARACTER(len=132) :: message_text
    LOGICAL                        :: lsuccess, lexists
    INTEGER                        :: idom, total_dom, fileID, vlistID, root_pe
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText

    ! rank of broadcast root PE
    root_pe = 0

    idom = 1
    rst_filename = "restart_"//TRIM(modeltype_str)//"_DOM"//TRIM(int2string(idom, "(i2.2)"))//".nc"

    ! test if the domain-dependent restart file exists:
    INQUIRE(file=TRIM(rst_filename), exist=lexists)
    ! otherwise, give a warning and resort to the old naming scheme:
    IF (.NOT. lexists) THEN
      CALL finish(routine, "Restart file not found! Expected name: restart_<model_type>_DOM<2 digit domain number>.nc")
    END IF

    ! Read all namelists used in the previous run
    ! and store them in a buffer. These values will overwrite the
    ! model default, and will later be overwritten if the user has
    ! specified something different for this integration.
    !
    ! Note: We read the namelists only once and assume that these
    !       are identical for all domains.
    !
    IF (my_process_is_mpi_workroot()) THEN
      fileID  = streamOpenRead(TRIM(rst_filename))
      ! check if the file could be opened
      IF (fileID < 0) THEN
        CALL cdiGetStringError(fileID, cdiErrorText)
        WRITE(message_text,'(4a)') 'File ', TRIM(rst_filename), &
             ' cannot be opened: ', TRIM(cdiErrorText)
        CALL finish(routine, TRIM(message_text))
      ENDIF

      vlistID = streamInqVlist(fileID)
    END IF
    CALL read_and_bcast_restart_namelists(vlistID, my_process_is_mpi_workroot(), &
      &                                   root_pe, p_comm_work)
    ! note: we don't close the restart file for domain 1 yet...
    CALL message(TRIM(routine), 'read namelists from restart file')

    total_dom = 1
    DO
      IF (idom > total_dom) EXIT

      ! Read all global attributs in the restart file and store them in a buffer.

      IF (idom > 1) INQUIRE(file=TRIM(rst_filename), exist=lexists)
      IF (my_process_is_mpi_workroot()) THEN
        WRITE(0,*) "streamOpenRead ", TRIM(rst_filename)
        ! note: we have opened the file for domain 1 already
        IF (idom > 1) THEN
          IF (lexists) THEN
            fileID  = streamOpenRead(TRIM(rst_filename))
            vlistID = streamInqVlist(fileID)
          ENDIF
        END IF
      END IF
      IF (lexists) THEN
        CALL read_and_bcast_attributes(vlistID, my_process_is_mpi_workroot(), root_pe, p_comm_work)
        IF (my_process_is_mpi_workroot()) THEN
          CALL streamClose(fileID)
        END IF

        CALL message(TRIM(routine), 'read global attributes from restart file')
      ELSE
        CALL message(TRIM(routine), 'Warning: domain not active at restart time')
      ENDIF

      ! since we do not know about the total number of domains yet,
      ! we have to ask the current restart file for this
      ! information:
      CALL get_restart_attribute( 'n_dom', total_dom )

      idom = idom + 1
      rst_filename = "restart_"//TRIM(modeltype_str)//"_DOM"//TRIM(int2string(idom, "(i2.2)"))//".nc"
    END DO

  END SUBROUTINE read_restart_header
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
    IF (ALLOCATED(private_vct))  DEALLOCATE(private_vct)
    ALLOCATE(private_vct(SIZE(vct)))
    private_vct(:) = vct(:)
    lvct_initialised = .TRUE.
  END SUBROUTINE set_restart_vct
  !------------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------
  !  depth based vertical coordinates
  SUBROUTINE set_restart_depth_lnd(zh, zf)
    REAL(wp), INTENT(in) :: zh(:), zf(:)
    IF (ldepth_lnd_initialised) RETURN
    IF (ALLOCATED(private_depth_lnd_half)) &
      &  DEALLOCATE(private_depth_lnd_half, private_depth_lnd_full)
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
    IF (ALLOCATED(private_height_snow_half))  &
      &  DEALLOCATE(private_height_snow_half, private_height_snow_full)
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
    IF (lrestart_initialised) RETURN
    !
    ! set CF-Convention required restart attributes
    !
    CALL set_restart_attribute('title',       TRIM(cf_global_info%title))
    CALL set_restart_attribute('institution', TRIM(cf_global_info%institution))
    CALL set_restart_attribute('source',      TRIM(cf_global_info%source))
    CALL set_restart_attribute('history',     TRIM(cf_global_info%history))
    CALL set_restart_attribute('references',  TRIM(cf_global_info%references))
    CALL set_restart_attribute('comment',     TRIM(cf_global_info%comment))
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
    if(lhamocc)CALL set_vertical_grid(ZA_OCEAN_SEDIMENT      , ks         )
    !
    ! define time axis
    !
    CALL set_time_axis(TAXIS_ABSOLUTE)
    CALL set_time_axis(TAXIS_RELATIVE)
    !
    private_nc  = nc
    private_nv  = nv
    private_ne  = ne
    !

!AD(9July-2013): The following condition seemed unnecessary. So after
!  discussing with DR we decided to get rid of it for time being

!    IF (.NOT. (lvct_initialised .OR. use_ocean_levels            &
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
  SUBROUTINE open_writing_restart_files(jg, restart_filename)
    INTEGER,          INTENT(IN) :: jg                   !< patch ID
    CHARACTER(LEN=*), INTENT(IN) :: restart_filename
    !
    INTEGER :: status, i ,j, k, ia, ihg, ivg, nlevp1
    REAL(wp), ALLOCATABLE :: levels(:), ubounds(:), lbounds(:), levels_sp(:)
    !
    CHARACTER(len=64) :: attribute_name
    CHARACTER(len=256) :: text_attribute
    REAL(wp) :: real_attribute(1)
    INTEGER :: int_attribute(1)
    LOGICAL :: bool_attribute
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText

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
    WRITE(message_text,'(t1,a,t70,a,t84,a,t94,a)') 'file', 'var list', 'file ID', 'restart'
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
      ! skip, if var_list does not match the current patch ID
      IF (var_lists(i)%p%patch_id /= jg) CYCLE
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
      IF (my_process_is_mpi_workroot()) THEN
        var_lists(i)%p%cdiFileID_restart = streamOpenWrite(TRIM(restart_filename), var_lists(i)%p%restart_type)
        var_lists(i)%p%filename          = TRIM(restart_filename)
        !
        IF (var_lists(i)%p%cdiFileID_restart < 0) THEN
          CALL cdiGetStringError(var_lists(i)%p%cdiFileID_restart, cdiErrorText)
          WRITE(message_text,'(a)') TRIM(cdiErrorText)
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
          CALL get_restart_attribute(ia, attribute_name, real_attribute(1))
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
          CALL get_restart_attribute(ia, attribute_name, int_attribute(1))
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
            int_attribute(1) = 1
          ELSE
            int_attribute(1) = 0
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
            IF (.NOT. use_ocean_levels) CYCLE
            var_lists(i)%p%cdiDepthFullZaxisID = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, &
              &                                           vgrid_def(ivg)%nlevels)
            CALL zaxisDefLevels(var_lists(i)%p%cdiDepthFullZaxisID, &
              &              private_depth_full)

          CASE (ZA_DEPTH_BELOW_SEA_HALF)
            IF (.NOT. use_ocean_levels) CYCLE
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
    
          CASE (ZA_OCEAN_SEDIMENT)
           if(lhamocc)then
           ! HAMOCC sediment
            var_lists(i)%p%cdiOceanSedGenericZaxisID = zaxisCreate(ZAXIS_GENERIC, &
              &                                           vgrid_def(ivg)%nlevels)
           ALLOCATE(levels(ks))
           ALLOCATE(levels_sp(ksp))
           CALL set_zlev(levels_sp, levels, ks, dzsed*1000._wp)
           CALL zaxisDefLevels(var_lists(i)%p%cdiOceanSedGenericZaxisID, REAL(levels,wp))
           DEALLOCATE(levels)
           DEALLOCATE(levels_sp)
          endif  
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
        CALL addVarListToVlist(var_lists(i), var_lists(i)%p%cdiVlistID, jg)
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
        IF (i == j)                        CYCLE
        IF (var_lists(j)%p%restart_opened) CYCLE
        IF (.NOT. var_lists(j)%p%lrestart) CYCLE
        IF (var_lists(j)%p%patch_id /= jg) CYCLE
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
          var_lists(j)%p%cdiOceanSedGenericZaxisID    = var_lists(i)%p%cdiOceanSedGenericZaxisID
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
            CALL addVarListToVlist(var_lists(j), var_lists(j)%p%cdiVlistID, jg)
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
  SUBROUTINE addVarListToVlist(this_list, vlistID, jg)
    TYPE (t_var_list), INTENT(inout) :: this_list
    INTEGER,           INTENT(inout) :: vlistID
    INTEGER,           INTENT(in)    :: jg
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
    LOGICAL :: lskip_timelev, lskip_extra_timelevs

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

      lskip_timelev = .FALSE.
      lskip_extra_timelevs = .FALSE.
#ifndef __NO_ICON_ATMO__
      IF (iequations == INH_ATMOSPHERE .AND. .NOT. (l_limited_area .AND. jg == 1)) THEN
        lskip_extra_timelevs = .TRUE.
      ENDIF

      ! get information about timelevel to be skipped for current field
      IF (element%field%info%tlev_source == TLEV_NNOW ) THEN
        IF (time_level == nnew(jg))                    lskip_timelev = .TRUE.
        ! this is needed to skip the extra time levels allocated for nesting
        IF (lskip_extra_timelevs .AND. time_level > 2) lskip_timelev = .TRUE.
      ELSE IF (element%field%info%tlev_source == TLEV_NNOW_RCF) THEN
        IF (time_level == nnew_rcf(jg))  lskip_timelev = .TRUE.
      ENDIF


      SELECT CASE (iequations)
      CASE(IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER)

        IF ( lskip_timelev                                  &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_EXPL &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_SI   ) CYCLE   ! skip field
      CASE default
        IF ( lskip_timelev ) CYCLE   ! skip field
      END SELECT
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
      ! set z axis ID (only if cdiZaxisID is not already set)
      !
      IF (info%cdiZaxisID < 0) THEN
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
        CASE (ZA_OCEAN_SEDIMENT)
          info%cdiZaxisID =  this_list%p%cdiOceanSedGenericZaxisID
        END SELECT
      END IF

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
      CALL vlistDefVarName(vlistID, varID, TRIM(info%name))
      !
      IF (info%cf%long_name /= '') CALL vlistDefVarLongname(vlistID, varID, TRIM(info%cf%long_name))
      IF (info%cf%units /= '') CALL vlistDefVarUnits(vlistID, varID, TRIM(info%cf%units))

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
                                & model_type,                  &
                                & opt_pvct,                    &
                                & opt_t_elapsed_phy,           &
                                & opt_lcall_phy, opt_sim_time, &
                                & opt_ndyn_substeps,           &
                                & opt_jstep_adv_marchuk_order, &
                                & opt_depth_lnd,               &
                                & opt_nlev_snow,               &
                                & opt_nice_class,              &
                                & opt_ndom,                    &
                                & opt_output_jfile,            &
                                & ocean_Zlevels,               &
                                & ocean_Zheight_CellMiddle,    &
                                & ocean_Zheight_CellInterfaces)

    TYPE(t_patch),       INTENT(IN) :: patch
    TYPE(t_datetime),    INTENT(IN) :: datetime
    INTEGER,             INTENT(IN) :: jstep                ! simulation step
    CHARACTER(len=*),    INTENT(IN) :: model_type           ! store model type

    REAL(wp), INTENT(IN), OPTIONAL :: opt_pvct(:)
    INTEGER,  INTENT(IN), OPTIONAL :: opt_depth_lnd               ! vertical levels soil model
    REAL(wp), INTENT(IN), OPTIONAL :: opt_t_elapsed_phy(:,:)
    LOGICAL , INTENT(IN), OPTIONAL :: opt_lcall_phy(:,:)
    REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time
    INTEGER,  INTENT(IN), OPTIONAL :: opt_ndyn_substeps
    INTEGER,  INTENT(IN), OPTIONAL :: opt_jstep_adv_marchuk_order
    INTEGER,  INTENT(IN), OPTIONAL :: opt_nlev_snow
    INTEGER,  INTENT(IN), OPTIONAL :: opt_nice_class
    INTEGER,  INTENT(IN), OPTIONAL :: opt_ndom                    !< no. of domains (appended to symlink name)
    INTEGER,  INTENT(IN), OPTIONAL :: opt_output_jfile(:)
    INTEGER,  INTENT(IN), OPTIONAL :: ocean_Zlevels
    REAL(wp), INTENT(IN), OPTIONAL :: ocean_Zheight_CellMiddle(:)
    REAL(wp), INTENT(IN), OPTIONAL :: ocean_Zheight_CellInterfaces(:)


    INTEGER :: klev, jg, kcell, kvert, kedge, max_cell_connectivity, max_vertex_connectivity
    INTEGER :: inlev_soil, inlev_snow, i, nice_class
    INTEGER :: ndepth    ! depth of n
    REAL(wp), ALLOCATABLE :: zlevels_full(:), zlevels_half(:)
    CHARACTER(len=MAX_CHAR_LENGTH)  :: string

    CHARACTER(len=MAX_CHAR_LENGTH) :: attname   ! attribute name
    INTEGER :: jp, jp_end   ! loop index and array size

    TYPE (t_keyword_list), POINTER :: keywords => NULL()

    IF (ltimer) CALL timer_start(timer_write_restart_file)

    !AD (14 Oct 2014):Necessary to avoid writing of attributes of
    !multiple domains in the restart file of the parent domain when
    !started from a restart file.
    !Funny: This call is already there in async restart but not here!
    !We should avoid maintaing two pieces of code doing same thing
    CALL delete_attributes()

    !----------------
    ! Initialization
    klev      = patch%nlev
    jg        = patch%id
    kcell     = patch%n_patch_cells_g
    kvert     = patch%n_patch_verts_g
    kedge     = patch%n_patch_edges_g
    max_cell_connectivity   = patch%cells%max_connectivity
    max_vertex_connectivity = patch%verts%max_connectivity

    CALL set_restart_attribute( 'current_caltime', datetime%caltime )
    CALL set_restart_attribute( 'current_calday' , datetime%calday )

    CALL set_restart_attribute( 'current_daysec' , datetime%daysec )

    CALL set_restart_attribute( 'nold_DOM'//TRIM(int2string(jg, "(i2.2)"))    , nold    (jg))
    CALL set_restart_attribute( 'nnow_DOM'//TRIM(int2string(jg, "(i2.2)"))    , nnow    (jg))
    CALL set_restart_attribute( 'nnew_DOM'//TRIM(int2string(jg, "(i2.2)"))    , nnew    (jg))
    CALL set_restart_attribute( 'nnow_rcf_DOM'//TRIM(int2string(jg, "(i2.2)")), nnow_rcf(jg))
    CALL set_restart_attribute( 'nnew_rcf_DOM'//TRIM(int2string(jg, "(i2.2)")), nnew_rcf(jg))

    ! set no. of domains
    IF (PRESENT(opt_ndom)) THEN
      CALL set_restart_attribute( 'n_dom', opt_ndom)
    ELSE
      CALL set_restart_attribute( 'n_dom', 1)
    END IF

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

    IF (PRESENT(opt_ndyn_substeps)) THEN
        WRITE(attname,'(a,i2.2)') 'ndyn_substeps_DOM',jg
        CALL set_restart_attribute( TRIM(attname), opt_ndyn_substeps )
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

    ndepth = 0
    IF (PRESENT(ocean_Zheight_CellMiddle) ) THEN
      IF (.not. PRESENT(ocean_Zheight_CellInterfaces) .or. .not. PRESENT(ocean_Zlevels) ) &
        CALL finish('create_restart_file','Ocean level parameteres not complete')
      IF (.not. use_ocean_levels) THEN
        ! initialize the ocean levels
        IF (ALLOCATED(private_depth_half))  &
          &  DEALLOCATE(private_depth_half, private_depth_full)
        ALLOCATE(private_depth_half(ocean_Zlevels+1), private_depth_full(ocean_Zlevels))
        private_depth_half(:) = ocean_Zheight_CellInterfaces(:)
        private_depth_full(:) = ocean_Zheight_CellMiddle(:)
      END IF
      ndepth = ocean_Zlevels
      use_ocean_levels = .TRUE.
     END IF


    IF (PRESENT(opt_output_jfile)) THEN
      DO i=1,SIZE(opt_output_jfile)
        WRITE(attname,'(a,i2.2)') 'output_jfile_',i
        CALL set_restart_attribute( TRIM(attname), opt_output_jfile(i) )
      END DO
    END IF

    IF (.NOT.PRESENT(opt_nice_class)) THEN
      nice_class = 1
    ELSE
      nice_class = opt_nice_class
    END IF

    CALL init_restart( 'ICON',            &! model name
                     & '1.4.00',          &! model version
                     & kcell, max_cell_connectivity,  &! total # of cells, # of vertices per cell
                     & kvert, max_vertex_connectivity,&! total # of vertices, # of vertices per dual cell
                     & kedge, 4,          &! total # of cells, shape of control volume for edge
                     & klev,              &! total # of vertical layers
                     & ndepth,            &! total # of depths below sea
                     & inlev_soil,        &! total # of depths below land (TERRA or JSBACH)
                     & inlev_snow,        &! total # of vertical snow layers (TERRA)
                     & nice_class         )! total # of ice classes (sea ice)

    CALL set_restart_time( iso8601(datetime) )  ! Time tag
    ! in preparation for move to mtime
    CALL set_restart_attribute('tc_startdate', iso8601extended(datetime))

    ! Open new file, write data, close and then clean-up.
    CALL associate_keyword("<gridfile>",   TRIM(get_filename_noext(patch%grid_filename)),  keywords)
    CALL associate_keyword("<idom>",       TRIM(int2string(jg, "(i2.2)")),                 keywords)
    CALL associate_keyword("<rsttime>",    TRIM(private_restart_time),                     keywords)
    CALL associate_keyword("<mtype>",      TRIM(model_type),                               keywords)
    ! replace keywords in file name
    string = TRIM(with_keywords(keywords, TRIM(restart_filename)))

    CALL open_writing_restart_files( jg, TRIM(string) )

    !CALL write_restart( patch )
    !CALL write_restart( patch, TRIM(restart_filename) )
    CALL write_restart( patch, TRIM(string) )

    CALL close_writing_restart_files(jg, opt_ndom)
    CALL finish_restart

    IF (ltimer) CALL timer_stop(timer_write_restart_file)

  END SUBROUTINE create_restart_file

  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE close_writing_restart_files(jg, opt_ndom)
    INTEGER,  INTENT(IN)           :: jg           !< patch ID
    INTEGER,  INTENT(IN), OPTIONAL :: opt_ndom     !< no. of domains (appended to symlink name)
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
          IF (PRESENT(opt_ndom)) THEN
            IF (opt_ndom > 1) THEN
              linkname = 'restart_'//TRIM(var_lists(i)%p%model_type)//"_DOM"//TRIM(int2string(jg, "(i2.2)"))//'.nc'
            ELSE
              linkname = 'restart_'//TRIM(var_lists(i)%p%model_type)//'_DOM01.nc'
            END IF
          ELSE
            linkname = 'restart_'//TRIM(var_lists(i)%p%model_type)//'_DOM01.nc'
          END IF
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
  SUBROUTINE write_restart(p_patch, restart_filename)
    TYPE(t_patch), INTENT(in) :: p_patch
    CHARACTER(len=*)          :: restart_filename

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
            CALL message('write_restart: ', &
              &  'Write netCDF2 restart for : '//TRIM(private_restart_time)//' restart filename: '//TRIM(restart_filename))
          CASE (FILETYPE_NC4)
            IF (var_lists(i)%p%compression_type == COMPRESS_ZIP) THEN
              CALL message('write_restart: ', &
                &  'Write compressed netCDF4 restart for : '//&
                &TRIM(private_restart_time)//' restart filename: '//restart_filename)
            ELSE
              CALL message('write_restart: ', &
                &  'Write netCDF4 restart for : '//TRIM(private_restart_time)//' restart filename: '//TRIM(restart_filename))
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

          ! skip var_list if it does not match the current patch ID
          IF (var_lists(j)%p%patch_id /= p_patch%id) CYCLE

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
  CALL message('write_restart: ', &
    &  'Finished Write netCDF2 restart for : '//TRIM(private_restart_time)//' restart filename: '//TRIM(restart_filename))
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
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//"::write_restart_var_list"
    TYPE (t_var_metadata), POINTER :: info
    TYPE (t_list_element), POINTER :: element
    TYPE (t_list_element), TARGET  :: start_with
    !
    REAL(wp), ALLOCATABLE :: r_out_wp(:) ! field gathered on I/O processor
    !
    INTEGER :: nindex, private_n, lev
    TYPE(t_comm_gather_pattern), POINTER :: gather_pattern
    !
    INTEGER :: time_level
    INTEGER :: jg, var_ref_pos
    LOGICAL :: lskip_timelev, lskip_extra_timelevs

    jg = p_patch%id

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

      lskip_timelev = .FALSE.
      lskip_extra_timelevs = .FALSE.
#ifndef __NO_ICON_ATMO__
      IF (iequations == INH_ATMOSPHERE .AND. .NOT. (l_limited_area .AND. jg == 1)) THEN
        lskip_extra_timelevs = .TRUE.
      ENDIF

      ! get information about timelevel to be skipped for current field
      IF (element%field%info%tlev_source == TLEV_NNOW ) THEN
        IF (time_level == nnew(jg))                    lskip_timelev = .TRUE.
        ! this is needed to skip the extra time levels allocated for nesting
        IF (lskip_extra_timelevs .AND. time_level > 2) lskip_timelev = .TRUE.
      ELSE IF (element%field%info%tlev_source == TLEV_NNOW_RCF) THEN
        IF (time_level == nnew_rcf(jg))  lskip_timelev = .TRUE.
      ENDIF

      SELECT CASE (iequations)
      CASE(IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER)

        IF ( lskip_timelev                                  &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_EXPL &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_SI   ) CYCLE   ! skip field
      CASE default
        IF ( lskip_timelev ) CYCLE   ! skip field
      END SELECT
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
        CALL finish(routine,'unknown grid type')
      END SELECT

      ALLOCATE(r_out_wp(MERGE(private_n, 0, my_process_is_mpi_workroot())))

      IF (info%ndims == 1) THEN
        CALL finish(routine, '1d arrays not handled yet.')
      ELSE IF (info%ndims == 2) THEN
        var_ref_pos = 3
        IF (info%lcontained)  var_ref_pos = info%var_ref_pos
        SELECT CASE(var_ref_pos)
        CASE (1)
          CALL exchange_data(in_array=element%field%r_ptr(nindex,:,:,1,1), &
            &                out_array=r_out_wp, gather_pattern=gather_pattern)
        CASE (2)
          CALL exchange_data(in_array=element%field%r_ptr(:,nindex,:,1,1), &
            &                out_array=r_out_wp, gather_pattern=gather_pattern)
        CASE (3)
          CALL exchange_data(in_array=element%field%r_ptr(:,:,nindex,1,1), &
            &                out_array=r_out_wp, gather_pattern=gather_pattern)
        CASE default
          CALL finish(routine, "internal error!")
        END SELECT
        !
        ! write data
        !
        IF (my_process_is_mpi_workroot()) THEN
          WRITE (0,*)' ... write ',info%name
          CALL streamWriteVar(this_list%p%cdiFileID_restart, info%cdiVarID, &
            &                 r_out_wp(:), 0)
        END IF
      ELSE IF (info%ndims == 3) THEN
        var_ref_pos = 4
        IF (info%lcontained)  var_ref_pos = info%var_ref_pos

        IF (my_process_is_mpi_workroot()) &
          write (0,*)' ... write ',info%name
        DO lev = 1, info%used_dimensions(2)

          SELECT CASE(var_ref_pos)
          CASE (1)
            CALL exchange_data( in_array=element%field%r_ptr(nindex,:,lev,:,1), &
              &                 out_array=r_out_wp, gather_pattern=gather_pattern)
          CASE (2)
            CALL exchange_data( in_array=element%field%r_ptr(:,nindex,lev,:,1), &
              &                 out_array=r_out_wp, gather_pattern=gather_pattern)
          CASE (3)
            CALL exchange_data( in_array=element%field%r_ptr(:,lev,nindex,:,1), &
              &                 out_array=r_out_wp, gather_pattern=gather_pattern)
          CASE (4)
            CALL exchange_data( in_array=element%field%r_ptr(:,lev,:,nindex,1), &
              &                 out_array=r_out_wp, gather_pattern=gather_pattern)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT

          IF (my_process_is_mpi_workroot()) &
            CALL streamWriteVarSlice(this_list%p%cdiFileID_restart, &
              &                      info%cdiVarID, lev-1, r_out_wp(:), 0)
        END DO
      ELSE IF (info%ndims > 3) THEN
        CALL finish(routine, 'arrays with more then three dimensions not handled yet.')
      ELSE
        CALL finish(routine,'dimension not set.')
      END IF

      !
      ! deallocate temporary global arrays
      !
      DEALLOCATE(r_out_wp)
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
        IF (var_lists(i)%p%cdiOceanSedGenericZaxisID /= CDI_UNDEFID) &
             CALL zaxisDestroy(var_lists(i)%p%cdiOceanSedGenericZaxisID)
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
        var_lists(i)%p%cdiOceanSedGenericZaxisID    = CDI_UNDEFID
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
    use_ocean_levels = .FALSE.
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
  SUBROUTINE read_restart_files(p_patch, opt_ndom)

    TYPE(t_patch), INTENT(in) :: p_patch
    INTEGER,       OPTIONAL, INTENT(in) :: opt_ndom
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::read_restart_files"
    !
    TYPE model_search
      CHARACTER(len=8) :: abbreviation
      INTEGER          :: key
    END type model_search
    TYPE(model_search) :: abbreviations(nvar_lists)
    !
    TYPE (t_list_element), POINTER   :: element
    TYPE (t_var_metadata), POINTER   :: info
    CHARACTER(len=80)                :: restart_filename, name
    INTEGER                          :: fileID, vlistID, gridID, zaxisID, taxisID,    &
      &                                 varID, idate, itime, ic, il, n, nfiles, i,   &
      &                                 iret, istat, key, vgrid, gdims(5), nindex,   &
      &                                 nmiss, nvars, root_pe, var_ref_pos, lev
    CHARACTER(len=8)                 :: model_type
    REAL(wp), POINTER                :: r1d(:), rptr2d(:,:), rptr3d(:,:,:)
    INTEGER, POINTER                 :: glb_index(:)
    CLASS(t_scatterPattern), POINTER :: scatter_pattern

    ! rank of broadcast root PE
    root_pe = 0

    IF (my_process_is_mpi_workroot()) THEN
      WRITE(0,*) "read_restart_files, nvar_lists=", nvar_lists
    END IF
    abbreviations(1:nvar_lists)%key          = 0
    abbreviations(1:nvar_lists)%abbreviation = ""
    key = 0
    n   = 1
    for_all_model_types: DO i = 1, nvar_lists
      ! skip var_list if it does not match the current patch ID
      IF (var_lists(i)%p%patch_id /= p_patch%id) CYCLE

      key = util_hashword(TRIM(var_lists(i)%p%model_type), LEN_TRIM(var_lists(i)%p%model_type), 0)
      IF (.NOT. ANY(abbreviations(1:n)%key == key)) THEN
        abbreviations(n)%abbreviation = var_lists(i)%p%model_type
        abbreviations(n)%key = key
        n = n+1
      ENDIF
    ENDDO for_all_model_types

    ALLOCATE(r1d(1),STAT=istat)
    IF (istat /= 0) THEN
      CALL finish('','allocation of r1d failed ...')
    ENDIF

    nfiles = n-1
    !
    for_all_files: DO n = 1, nfiles
      model_type =TRIM(abbreviations(n)%abbreviation)
      IF (PRESENT(opt_ndom)) THEN
        IF (opt_ndom > 1) THEN
          restart_filename = 'restart_'//TRIM(model_type)//"_DOM"//TRIM(int2string(p_patch%id, "(i2.2)"))//'.nc'
        ELSE
          restart_filename = 'restart_'//TRIM(model_type)//'_DOM01.nc'
        END IF
      ELSE
        restart_filename = 'restart_'//TRIM(model_type)//'_DOM01.nc'
      END IF
      !
      IF (.NOT. util_islink(TRIM(restart_filename))) THEN
        iret = util_rename(TRIM(restart_filename), TRIM(restart_filename)//'.bak')
        WRITE(0,*) "util_rename returned:", iret
        iret = util_symlink(TRIM(restart_filename)//'.bak', TRIM(restart_filename))
        WRITE(0,*) "util_symlink:", iret
      ENDIF
      !
      name          = TRIM(restart_filename)//CHAR(0)

      IF (my_process_is_mpi_workroot()) THEN
        WRITE(0,*) "streamOpenRead ", TRIM(restart_filename)

        fileID  = streamOpenRead(TRIM(name))
        vlistID = streamInqVlist(fileID)
        taxisID = vlistInqTaxis(vlistID)
        !
        idate   = taxisInqVdate(taxisID)
        itime   = taxisInqVtime(taxisID)
        !
        WRITE(message_text,'(a,i8.8,a,i6.6,a,a)') &
          'Read restart for : ', idate, 'T', itime, 'Z from ',TRIM(restart_filename)
      END IF

      CALL message('read_restart_files',message_text)
      !
      IF (my_process_is_mpi_workroot()) THEN
        nvars = vlistNvars(vlistID)
      END IF
      CALL p_bcast(nvars, root_pe, comm=p_comm_work)

      !AD: This call is already made in read_restart_header. Why do it twice?
      CALL read_and_bcast_attributes(vlistID, my_process_is_mpi_workroot(), root_pe, p_comm_work)
      !
      for_all_vars: DO varID = 0, (nvars-1)
        !
        IF (my_process_is_mpi_workroot()) THEN
          CALL vlistInqVarName(vlistID, varID, name)
        END IF
        CALL p_bcast(name, root_pe, comm=p_comm_work)
        !
        for_all_lists: DO i = 1, nvar_lists
          ! skip var_list if it does not match the current patch ID
          IF (var_lists(i)%p%patch_id   /= p_patch%id) CYCLE
          IF (var_lists(i)%p%model_type /= model_type) CYCLE

          element => find_list_element(var_lists(i), TRIM(name))
          IF (ASSOCIATED(element)) THEN
            IF (element%field%info%lrestart) THEN
              !
              info => element%field%info
              !
              ! allocate temporary global array on output processor
              ! and gather field from other processors
              !
              NULLIFY(rptr2d)
              NULLIFY(rptr3d)
              !
              IF (my_process_is_mpi_workroot()) THEN
                gridID  = vlistInqVarGrid(vlistID, varID)
                zaxisID = vlistInqVarZaxis(vlistID, varID)
                vgrid   = zaxisInqType(zaxisID)
                ic = gridInqSize(gridID)
                IF (SIZE(r1d, 1) /= ic) THEN
                  DEALLOCATE(r1d)
                  ALLOCATE(r1d(ic),STAT=istat)
                  IF (istat /= 0) THEN
                    CALL finish('','allocation of r1d failed ...')
                  ENDIF
                END IF
                IF (vgrid == ZAXIS_SURFACE) THEN
                  il = 1
                ELSE
                  il = zaxisInqSize(zaxisID)
                ENDIF
              END IF

              CALL p_bcast(il, root_pe, comm=p_comm_work)
              !
              IF (info%lcontained) THEN
                nindex = info%ncontained
              ELSE
                nindex = 1
              ENDIF
              !

              SELECT CASE(info%ndims)
              CASE (2)
                var_ref_pos = 3
                IF (info%lcontained)  var_ref_pos = info%var_ref_pos
                SELECT CASE(var_ref_pos)
                CASE (1)
                  rptr2d => element%field%r_ptr(nindex,:,:,1,1)
                CASE (2)
                  rptr2d => element%field%r_ptr(:,nindex,:,1,1)
                CASE (3)
                  rptr2d => element%field%r_ptr(:,:,nindex,1,1)
                CASE default
                  CALL finish(routine, "internal error!")
                END SELECT
              CASE (3)
                var_ref_pos = 4
                IF (info%lcontained)  var_ref_pos = info%var_ref_pos
                SELECT CASE(var_ref_pos)
                CASE (1)
                  rptr3d => element%field%r_ptr(nindex,:,:,:,1)
                CASE (2)
                  rptr3d => element%field%r_ptr(:,nindex,:,:,1)
                CASE (3)
                  rptr3d => element%field%r_ptr(:,:,nindex,:,1)
                CASE (4)
                  rptr3d => element%field%r_ptr(:,:,:,nindex,1)
                CASE default
                  CALL finish(routine, "internal error!")
                END SELECT
              CASE DEFAULT
                CALL finish(routine, "internal error!")
              END SELECT

              SELECT CASE (info%hgrid)
              CASE (GRID_UNSTRUCTURED_CELL)
                scatter_pattern => p_patch%comm_pat_scatter_c
              CASE (GRID_UNSTRUCTURED_VERT)
                scatter_pattern => p_patch%comm_pat_scatter_v
              CASE (GRID_UNSTRUCTURED_EDGE)
                scatter_pattern => p_patch%comm_pat_scatter_e
              CASE default
                CALL finish('out_stream','unknown grid type')
              END SELECT

              DO lev = 1, il
                IF (my_process_is_mpi_workroot()) THEN
                  CALL streamReadVarSlice(fileID, varID, lev-1, r1d, nmiss)
                ENDIF
                IF (info%ndims == 2) THEN
                  CALL scatter_pattern%distribute(r1d, rptr2d, .FALSE.)
                ELSE
                  CALL scatter_pattern%distribute(r1d, rptr3d(:,lev,:), .FALSE.)
                ENDIF
              END DO

              !
              ! deallocate temporary global arrays
              !
              !
              IF (my_process_is_mpi_workroot()) THEN
                WRITE (0,*) ' ... read ',TRIM(element%field%info%name)
              ENDIF
              CYCLE for_all_vars
            ENDIF
          ENDIF
        ENDDO for_all_lists
        CALL message('reading_restart_file','Variable '//TRIM(name)//' not defined.')
      ENDDO for_all_vars
      !
      IF (my_process_is_mpi_workroot())  CALL streamClose(fileID)
      !
    ENDDO for_all_files

    DEALLOCATE (r1d)
    !
    CALL message('','')
    CALL message('',separator)
    CALL message('','')
    !
  END SUBROUTINE read_restart_files
  !-------------------------------------------------------------------------
END MODULE mo_io_restart
