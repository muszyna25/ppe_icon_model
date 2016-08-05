!> Module for writing restart files (synchronously) and for reading restart files.
!!
!! Note: The asynchronous implementation of the restart output can be
!!       found in the module "mo_restart_async"
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
MODULE mo_restart

  USE mo_kind,                  ONLY: wp
  USE mo_mpi,                   ONLY: p_comm_work,p_bcast
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_fortran_tools,         ONLY: t_ptr_2d
  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH, TLEV_NNOW, TLEV_NNOW_RCF, SUCCESS
  USE mo_var_metadata_types,    ONLY: t_var_metadata
  USE mo_linked_list,           ONLY: t_var_list, t_list_element
  USE mo_var_list,              ONLY: nvar_lists, var_lists, get_var_timelevel, find_list_element
  USE mo_cdi,                   ONLY: FILETYPE_NC, FILETYPE_NC2, FILETYPE_NC4, ZAXIS_SURFACE, CDI_UNDEFID, COMPRESS_ZIP, &
                                    & streamOpenRead, streamInqVlist, vlistInqTaxis, taxisInqVdate, taxisInqVtime, vlistNvars, &
                                    & vlistInqVarGrid, gridInqSize, vlistInqVarZaxis, zaxisInqType, zaxisInqSize, &
                                    & streamClose, streamWriteVarSlice, streamWriteVar, vlistInqVarName
  USE mo_util_cdi,              ONLY: cdiGetStringError
  USE mo_cdi_constants,         ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, &
                                    & ZA_HYBRID, ZA_HYBRID_HALF, ZA_DEPTH_BELOW_LAND, ZA_DEPTH_BELOW_LAND_P1, ZA_DEPTH_RUNOFF_S, &
                                    & ZA_DEPTH_RUNOFF_G, ZA_SNOW, ZA_SNOW_HALF, ZA_TOA, ZA_HEIGHT_2M, ZA_HEIGHT_10M, &
                                    & ZA_LAKE_BOTTOM, ZA_LAKE_BOTTOM_HALF, ZA_MIX_LAYER, ZA_SEDIMENT_BOTTOM_TW_HALF, &
                                    & ZA_DEPTH_BELOW_SEA, ZA_DEPTH_BELOW_SEA_HALF, ZA_GENERIC_ICE, ZA_OCEAN_SEDIMENT, ZA_COUNT
  USE mo_util_restart,          ONLY: t_v_grid, t_restart_cdi_ids, set_vertical_grid, setGeneralRestartAttributes, &
                                    & setDynamicPatchRestartAttributes, setPhysicsRestartAttributes, create_restart_file_link, &
                                    & getRestartFilename, getLevelPointers
  USE mo_util_string,           ONLY: int2string, separator
  USE mo_util_file,             ONLY: util_symlink, util_rename, util_islink
  USE mo_util_hash,             ONLY: util_hashword
  USE mo_restart_namelist,   ONLY: read_and_bcast_restart_namelists
  USE mo_restart_attributes, ONLY: t_RestartAttributeList, RestartAttributeList_make, setAttributesForRestarting
  USE mo_datetime,              ONLY: t_datetime, iso8601
  USE mo_run_config,            ONLY: ltimer, restart_filename
  USE mo_timer,                 ONLY: timer_start, timer_stop, timer_write_restart_file

  USE mo_dynamics_config,       ONLY: iequations, nold, nnow, nnew, nnew_rcf, nnow_rcf
  USE mo_grid_config,           ONLY: l_limited_area, n_dom

#ifndef __NO_ICON_ATMO__
!LK comment: should not be here !!!!!! polution of namespace !!!!!!
!GZ: but then we need an alternative method of skipping unnecessary time levels!
  USE mo_impl_constants,        ONLY: IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER, LEAPFROG_EXPL, LEAPFROG_SI, INH_ATMOSPHERE
  USE mo_ha_dyn_config,         ONLY: ha_dyn_config
#endif

  USE mo_model_domain,          ONLY: t_patch
  USE mo_mpi,                   ONLY: my_process_is_mpi_workroot, my_process_is_mpi_test, p_comm_rank
  USE mo_communication,         ONLY: t_comm_gather_pattern, exchange_data, t_scatterPattern

#ifndef __NO_ICON__OCEAN
  USE mo_ocean_nml,                ONLY: lhamocc
  USE mo_sedmnt,                ONLY: ks, ksp, dzsed
  USE mo_math_utilities,        ONLY: set_zlev
#endif

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'

  PUBLIC :: create_restart_file
  PUBLIC :: read_restart_files
  PUBLIC :: read_restart_header

  INTEGER, SAVE :: nv_grids = 0
  TYPE(t_v_grid) :: vgrid_def(ZA_COUNT)

  CHARACTER(len=32) :: private_restart_time = ''
  REAL(wp), ALLOCATABLE :: private_vct(:)
  REAL(wp), ALLOCATABLE :: private_depth_full(:),  private_depth_half(:)
  REAL(wp), ALLOCATABLE :: private_depth_lnd_full(:),  private_depth_lnd_half(:)
  REAL(wp), ALLOCATABLE :: private_height_snow_half(:), private_height_snow_full(:)

  LOGICAL, SAVE :: lvct_initialised         = .FALSE.
  LOGICAL, SAVE :: use_ocean_levels       = .FALSE.
  LOGICAL, SAVE :: ldepth_lnd_initialised   = .FALSE.
  LOGICAL, SAVE :: lheight_snow_initialised = .FALSE.

  LOGICAL, SAVE :: lrestart_initialised     = .FALSE.

  INTEGER :: private_nc  = -1
  INTEGER :: private_nv  = -1
  INTEGER :: private_ne  = -1

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_restart'

  !------------------------------------------------------------------------------------------------
CONTAINS
  !------------------------------------------------------------------------------------------------
  !> Reads attributes and namelists for all available domains from restart file.
  !>
  !> XXX: The code that I found READ the restart attributes from all the files,
  !>      *appending* them to the list of attributes. Since all restart files contain the same attributes,
  !>      this ONLY served to duplicate them. The current code just ignores the attributes IN all but the first file,
  !>      just like the original code ignored the namelists IN all but the first file.
  !>      However, it might be a good idea to add some consistency checking on the attributes/namelists of the other files.
  SUBROUTINE read_restart_header(modeltype_str)
    CHARACTER(LEN=*), INTENT(IN)  :: modeltype_str
    ! local variables
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//"::read_restart_header"
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: rst_filename
    CHARACTER(len=132) :: message_text
    LOGICAL                        :: lsuccess, lexists
    INTEGER                        :: idom, total_dom, fileID, vlistID, myRank
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText
    INTEGER, PARAMETER :: root_pe = 0
    TYPE(t_RestartAttributeList), SAVE, POINTER :: restartAttributes

    ! rank of broadcast root PE
    myRank = p_comm_rank(p_comm_work)

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

    ! Note: We read the namelists AND attributes only once and assume that these
    !       are identical for all domains (which IS guaranteed by the way the restart files are written).

    IF (myRank == root_pe) THEN
        fileID  = streamOpenRead(TRIM(rst_filename))
        ! check if the file could be opened
        IF (fileID < 0) THEN
            CALL cdiGetStringError(fileID, cdiErrorText)
            WRITE(message_text,'(4a)') 'File ', TRIM(rst_filename), ' cannot be opened: ', TRIM(cdiErrorText)
            CALL finish(routine, TRIM(message_text))
        ENDIF

        vlistID = streamInqVlist(fileID)
    END IF
    CALL read_and_bcast_restart_namelists(vlistID, root_pe, p_comm_work)
    restartAttributes => RestartAttributeList_make(vlistID, root_pe, p_comm_work)
    CALL setAttributesForRestarting(restartAttributes)
    IF (myRank == root_pe) CALL streamClose(fileID)

    CALL message(TRIM(routine), 'read namelists AND attributes from restart file')

    ! since we do not know about the total number of domains yet,
    ! we have to ask the restart file for this information:
    total_dom = restartAttributes%getInteger( 'n_dom' )

    ! check whether we have all the restart files we expect
    IF (myRank == root_pe) THEN
        DO idom = 2, total_dom
            rst_filename = "restart_"//TRIM(modeltype_str)//"_DOM"//TRIM(int2string(idom, "(i2.2)"))//".nc"
            IF (idom > 1) INQUIRE(file=TRIM(rst_filename), exist=lexists)
            IF (lexists) THEN
                CALL message(TRIM(routine), 'read global attributes from restart file')
            ELSE
                CALL message(TRIM(routine), 'Warning: domain not active at restart time')
            ENDIF
        END DO
    END IF

  END SUBROUTINE read_restart_header

  !  VCT as in echam (first half of vector contains a and second half b
  SUBROUTINE set_restart_vct(vct)
    REAL(wp), INTENT(in) :: vct(:)
    IF (lvct_initialised) RETURN
    IF (ALLOCATED(private_vct))  DEALLOCATE(private_vct)
    ALLOCATE(private_vct(SIZE(vct)))
    private_vct(:) = vct(:)
    lvct_initialised = .TRUE.
  END SUBROUTINE set_restart_vct

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

  !  height based vertical coordinates for multi layer snow model (TERRA)
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

  SUBROUTINE defineRestartAttributes(restartAttributes, datetime, jstep, opt_ndom, opt_ndyn_substeps, &
                                    &opt_jstep_adv_marchuk_order, opt_output_jfile, opt_sim_time, opt_t_elapsed_phy, opt_lcall_phy)
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    TYPE(t_datetime), INTENT(IN) :: datetime
    INTEGER, VALUE :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_ndom, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_output_jfile(:)
    REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time, opt_t_elapsed_phy(:,:)
    LOGICAL , INTENT(IN), OPTIONAL :: opt_lcall_phy(:,:)

    INTEGER :: i, jg, effectiveDomainCount
    CHARACTER(LEN = 2) :: jgString
    CHARACTER(LEN = :), ALLOCATABLE :: prefix
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":defineRestartAttributes"

    ! first the attributes that are independent of the domain
    effectiveDomainCount = 1
    IF(PRESENT(opt_ndom)) effectiveDomainCount = opt_ndom
    CALL setGeneralRestartAttributes(restartAttributes, datetime, effectiveDomainCount, jstep, opt_output_jfile)

    ! now the stuff that depends on the domain
    DO jg = 1, n_dom
        jgString = TRIM(int2string(jg, "(i2.2)"))
        CALL setDynamicPatchRestartAttributes(restartAttributes, jg, nold(jg), nnow(jg), nnew(jg), nnow_rcf(jg), nnew_rcf(jg))

        !----------------
        ! additional restart-output for nonhydrostatic model
        IF (PRESENT(opt_sim_time)) CALL restartAttributes%setReal('sim_time_DOM'//jgString, opt_sim_time )

        !-------------------------------------------------------------
        ! DR
        ! WORKAROUND FOR FIELDS WHICH NEED TO GO INTO THE RESTART FILE,
        ! BUT SO FAR CANNOT BE HANDELED CORRECTLY BY ADD_VAR OR
        ! SET_RESTART_ATTRIBUTE
        !-------------------------------------------------------------

        IF (PRESENT(opt_ndyn_substeps)) CALL restartAttributes%setInteger('ndyn_substeps_DOM'//jgString, opt_ndyn_substeps)
        IF (PRESENT(opt_jstep_adv_marchuk_order)) CALL restartAttributes%setInteger('jstep_adv_marchuk_order_DOM'//jgString, &
                                                                                   &opt_jstep_adv_marchuk_order)

        IF (PRESENT(opt_t_elapsed_phy) .AND. PRESENT(opt_lcall_phy)) THEN
            CALL setPhysicsRestartAttributes(restartAttributes, jg, opt_t_elapsed_phy(jg,:), opt_lcall_phy(jg,:))
        ENDIF
    END DO
  END SUBROUTINE defineRestartAttributes

  SUBROUTINE init_restart(nc, nv, ne, &
       &                  nlev, ndepth, nlev_soil,   &
       &                  nlev_snow, nice_class, nlev_sediment)
    INTEGER,          INTENT(in) :: nc
    INTEGER,          INTENT(in) :: nv
    INTEGER,          INTENT(in) :: ne
    INTEGER,          INTENT(in) :: nlev
    INTEGER,          INTENT(in) :: ndepth
    INTEGER,          INTENT(in) :: nlev_soil
    INTEGER,          INTENT(in) :: nlev_snow
    INTEGER,          INTENT(in) :: nice_class
    INTEGER,          INTENT(IN) :: nlev_sediment

    INTEGER :: error
    REAL(wp), ALLOCATABLE :: levels(:), levels_sp(:)
    CHARACTER(*), PARAMETER :: routine = modname//":init_restart"

    IF (lrestart_initialised) RETURN

    ! define vertical grids

    CALL set_vertical_grid(vgrid_def, nv_grids, ZA_SURFACE, 0._wp)
    CALL set_vertical_grid(vgrid_def, nv_grids, ZA_HYBRID, nlev)
    CALL set_vertical_grid(vgrid_def, nv_grids, ZA_HYBRID_HALF, nlev+1)
    CALL set_vertical_grid(vgrid_def, nv_grids, ZA_HEIGHT_2M, 2._wp)
    CALL set_vertical_grid(vgrid_def, nv_grids, ZA_HEIGHT_10M, 10._wp)
    CALL set_vertical_grid(vgrid_def, nv_grids, ZA_LAKE_BOTTOM, 1._wp)
    CALL set_vertical_grid(vgrid_def, nv_grids, ZA_LAKE_BOTTOM_HALF, 1._wp)
    CALL set_vertical_grid(vgrid_def, nv_grids, ZA_MIX_LAYER, 1._wp)
    CALL set_vertical_grid(vgrid_def, nv_grids, ZA_SEDIMENT_BOTTOM_TW_HALF, 0._wp)
    CALL set_vertical_grid(vgrid_def, nv_grids, ZA_TOA, 1._wp)
    CALL set_vertical_grid(vgrid_def, nv_grids, ZA_GENERIC_ICE, 1._wp)
    CALL set_vertical_grid(vgrid_def, nv_grids, ZA_DEPTH_RUNOFF_S, 1._wp)
    CALL set_vertical_grid(vgrid_def, nv_grids, ZA_DEPTH_RUNOFF_G, 1._wp)
    IF(ldepth_lnd_initialised) CALL set_vertical_grid(vgrid_def, nv_grids, ZA_DEPTH_BELOW_LAND, private_depth_lnd_full)
    IF(ldepth_lnd_initialised) CALL set_vertical_grid(vgrid_def, nv_grids, ZA_DEPTH_BELOW_LAND_P1, private_depth_lnd_half)
    IF(lheight_snow_initialised) CALL set_vertical_grid(vgrid_def, nv_grids, ZA_SNOW, private_height_snow_full)
    IF(lheight_snow_initialised) CALL set_vertical_grid(vgrid_def, nv_grids, ZA_SNOW_HALF, private_height_snow_half)
    IF(use_ocean_levels) CALL set_vertical_grid(vgrid_def, nv_grids, ZA_DEPTH_BELOW_SEA, private_depth_full)
    IF(use_ocean_levels) CALL set_vertical_grid(vgrid_def, nv_grids, ZA_DEPTH_BELOW_SEA_HALF, private_depth_half)
    IF(lhamocc) THEN
        ! HAMOCC sediment
        ALLOCATE(levels(nlev_sediment), levels_sp(nlev_sediment + 1), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

        CALL set_zlev(levels_sp, levels, nlev_sediment, dzsed*1000._wp)
        CALL set_vertical_grid(vgrid_def, nv_grids, ZA_OCEAN_SEDIMENT, REAL(levels,wp))

        DEALLOCATE(levels, levels_sp)
    END IF

    private_nc  = nc
    private_nv  = nv
    private_ne  = ne

!AD(9July-2013): The following condition seemed unnecessary. So after
!  discussing with DR we decided to get rid of it for time being

!    IF (.NOT. (lvct_initialised .OR. use_ocean_levels            &
!      & .OR. ldepth_lnd_initialised )) THEN
!      CALL finish('init_restart','none of the vertical grids is initialised')
!      ! more consistency checks need to follow
!    ENDIF

    lrestart_initialised = .TRUE.
  END SUBROUTINE init_restart

  ! Loop over all the output streams and open the associated files. Set
  ! unit numbers (file IDs) for all streams associated with a file.
  SUBROUTINE open_writing_restart_files(patch, restart_filename, restartAttributes, cdiIds, datetime)
    TYPE(t_patch), INTENT(IN) :: patch
    CHARACTER(LEN=*), INTENT(IN) :: restart_filename
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    TYPE(t_restart_cdi_ids), INTENT(INOUT) :: cdiIds(:)
    TYPE(t_datetime), INTENT(IN) :: datetime

    INTEGER :: status, i ,j, jg, ia, ihg, ivg

    CHARACTER(len=64) :: attribute_name
    CHARACTER(len=256) :: text_attribute
    REAL(wp) :: real_attribute(1)
    INTEGER :: int_attribute(1), axisIds(ZA_COUNT)
    LOGICAL :: bool_attribute
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":open_writing_restart_files"

    IF (my_process_is_mpi_test()) RETURN

    jg = patch%id

    ! first set restart file name

    IF (private_restart_time == '') THEN
      CALL finish('open_restart_files','restart time not set')
    ENDIF

    ! print header for restart var_lists

    CALL message('',separator)
    CALL message('','')
    CALL message('','Open restart files:')
    CALL message('','')
    WRITE(message_text,'(t1,a,t70,a,t84,a,t94,a)') 'file', 'var list', 'file ID', 'restart'
    CALL message('',message_text)
    CALL message('','')

    ! Loop over all var_lists and open the associated files. Set
    ! file IDs if necessary.

    DO i = 1, nvar_lists
      var_lists(i)%p%first = .FALSE.
    ENDDO

    DO i = 1, nvar_lists

      ! skip, if file is already opened

      IF (var_lists(i)%p%restart_opened) CYCLE

      ! skip, if var_list is not required for restarting

      IF (.NOT. var_lists(i)%p%lrestart) CYCLE

      ! skip, if var_list does not match the current patch ID
      IF (var_lists(i)%p%patch_id /= jg) CYCLE

      ! check restart file type

      SELECT CASE (var_lists(i)%p%restart_type)
      CASE (FILETYPE_NC)
        CALL finish('open_restart_files','netCDF classic not supported')
      CASE (FILETYPE_NC2, FILETYPE_NC4)
        ! this is ok, both formats can write more than 2GB files
      CASE default
        CALL finish('open_restart_files','unknown restart_type')
      END SELECT

      var_lists(i)%p%first = .TRUE.

      IF (my_process_is_mpi_workroot()) THEN
        IF(lvct_initialised) THEN
            CALL cdiIds(i)%openRestartAndCreateIds(TRIM(restart_filename), var_lists(i)%p%restart_type, restartAttributes, &
                                        &patch%n_patch_cells_g, patch%n_patch_verts_g, patch%n_patch_edges_g, &
                                        &patch%geometry_info%cell_type, vgrid_def(1:nv_grids), private_vct)
        ELSE
            CALL cdiIds(i)%openRestartAndCreateIds(TRIM(restart_filename), var_lists(i)%p%restart_type, restartAttributes, &
                                        &patch%n_patch_cells_g, patch%n_patch_verts_g, patch%n_patch_edges_g, &
                                        &patch%geometry_info%cell_type, vgrid_def(1:nv_grids))
        END IF
        ! set the related fields IN the var_lists

        var_lists(i)%p%filename = TRIM(restart_filename)
        var_lists(i)%p%restart_opened = .TRUE.

        ! 6. add variables

        CALL addVarListToVlist(var_lists(i), cdiIds(i), jg)

        WRITE(message_text,'(t1,a49,t50,a31,t84,i6,t94,l5)')        &
             restart_filename, var_lists(i)%p%name,             &
             cdiIds(i)%file, var_lists(i)%p%lrestart
        CALL message('',message_text)
      ENDIF



      ! loop over all other output var_lists eventually corresponding to the same file

      DO j = 1, nvar_lists

        IF (i == j)                        CYCLE
        IF (var_lists(j)%p%restart_opened) CYCLE
        IF (.NOT. var_lists(j)%p%lrestart) CYCLE
        IF (var_lists(j)%p%patch_id /= jg) CYCLE

        IF (var_lists(j)%p%restart_type /= var_lists(i)%p%restart_type) THEN
          CALL finish('open_output_streams', 'different file types for the same restart file')
        ENDIF

        IF (var_lists(i)%p%model_type == var_lists(j)%p%model_type) THEN
          var_lists(j)%p%restart_opened = .TRUE.
          var_lists(j)%p%filename = var_lists(i)%p%filename

          ! set file IDs of all associated restart files

          cdiIds(j)%file = cdiIds(i)%file   ! IN write_restart() we USE this field to identify the var_lists that USE the same file, which IS why we need to copy this over.

          ! add variables to already existing cdi vlists

          IF (my_process_is_mpi_workroot()) THEN

            CALL addVarListToVlist(var_lists(j), cdiIds(i), jg)

            WRITE(message_text,'(t1,a49,t50,a31,t84,i6,t94,l5)')        &
                 restart_filename, var_lists(j)%p%name,             &
                 cdiIds(j)%file, var_lists(j)%p%lrestart
            CALL message('',message_text)
          ENDIF
        ENDIF
      ENDDO


      IF (my_process_is_mpi_workroot() .AND. var_lists(i)%p%first) THEN
        CALL cdiIds(i)%finalizeVlist(datetime)
      ENDIF

    END DO

    CALL message('','')

  END SUBROUTINE open_writing_restart_files

  !------------------------------------------------------------------------------------------------

  ! define variables and attributes

  SUBROUTINE addVarListToVlist(this_list, cdiIds, jg)
    TYPE (t_var_list), INTENT(inout) :: this_list
    TYPE(t_restart_cdi_ids), INTENT(INOUT) :: cdiIds
    INTEGER, VALUE :: jg

    TYPE (t_var_metadata), POINTER :: info
    TYPE (t_list_element), POINTER :: element
    TYPE (t_list_element), TARGET  :: start_with

    INTEGER :: time_level
    LOGICAL :: lskip_timelev, lskip_extra_timelevs

    element => start_with
    element%next_list_element => this_list%p%first_list_element

    for_all_list_elements: DO

      element => element%next_list_element
      IF (.NOT.ASSOCIATED(element)) EXIT

      ! retrieve information from actual linked list element

      info => element%field%info

      ! skip this field ?

      IF (.NOT. info%lrestart) CYCLE

      ! skip this field because of wrong time index ?

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

      CALL cdiIds%defineVariable(info)
    ENDDO for_all_list_elements

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


    INTEGER :: klev, jg, kcell, kvert, kedge
    INTEGER :: inlev_soil, inlev_snow, i, nice_class, error
    INTEGER :: ndepth    ! depth of n
    REAL(wp), ALLOCATABLE :: zlevels_full(:), zlevels_half(:)
    CHARACTER(len=MAX_CHAR_LENGTH)  :: string
    TYPE(t_restart_cdi_ids), ALLOCATABLE :: cdiIds(:)

    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":create_restart_file"

    IF (ltimer) CALL timer_start(timer_write_restart_file)

    !----------------
    ! Initialization
    klev      = patch%nlev
    jg        = patch%id
    kcell     = patch%n_patch_cells_g
    kvert     = patch%n_patch_verts_g
    kedge     = patch%n_patch_edges_g

    restartAttributes => RestartAttributeList_make()
    CALL defineRestartAttributes(restartAttributes, datetime, jstep, opt_ndom, opt_ndyn_substeps, &
                                &opt_jstep_adv_marchuk_order, opt_output_jfile, opt_sim_time, opt_t_elapsed_phy, opt_lcall_phy)

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


    IF (.NOT.PRESENT(opt_nice_class)) THEN
      nice_class = 1
    ELSE
      nice_class = opt_nice_class
    END IF

    CALL init_restart( kcell,             &! total # of cells
                     & kvert,             &! total # of vertices
                     & kedge,             &! total # of cells
                     & klev,              &! total # of vertical layers
                     & ndepth,            &! total # of depths below sea
                     & inlev_soil,        &! total # of depths below land (TERRA or JSBACH)
                     & inlev_snow,        &! total # of vertical snow layers (TERRA)
                     & nice_class,        &! total # of ice classes (sea ice)
                     & ks)                 ! total # of sediment layers (HAMOCC)
    ALLOCATE(cdiIds(nvar_lists), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    DO i = 1, SIZE(cdiIds, 1)
        CALL cdiIds(i)%init()
    END DO

    private_restart_time = iso8601(datetime)
    string = getRestartFilename(patch%grid_filename, jg, datetime, model_type)

    CALL open_writing_restart_files(patch, TRIM(string), restartAttributes, cdiIds, datetime)

    CALL write_restart(patch, cdiIds)

    CALL close_writing_restart_files(jg, cdiIds, opt_ndom)
    CALL finish_restart()

    CALL restartAttributes%destruct()
    DEALLOCATE(cdiIds)
    DEALLOCATE(restartAttributes)

    IF (ltimer) CALL timer_stop(timer_write_restart_file)

  END SUBROUTINE create_restart_file

  !------------------------------------------------------------------------------------------------

  SUBROUTINE close_writing_restart_files(jg, cdiIds, opt_ndom)
    INTEGER,  INTENT(IN)           :: jg           !< patch ID
    TYPE(t_restart_cdi_ids), INTENT(INOUT) :: cdiIds(:)
    INTEGER,  INTENT(IN), OPTIONAL :: opt_ndom     !< no. of domains (appended to symlink name)

    ! Loop over all the output streams and close the associated files, set
    ! opened to false

    INTEGER :: i, j, iret, fileID, vlistID

    IF (my_process_is_mpi_test()) RETURN

    CALL message('',separator)
    CALL message('','')
    CALL message('','Close restart files:')
    CALL message('','')
    WRITE(message_text,'(t1,a,t50,a,t84,a)') 'file', 'link target', 'file ID'
    CALL message('',message_text)
    CALL message('','')

    close_all_lists: DO i = 1, nvar_lists

      IF (var_lists(i)%p%restart_opened) THEN
        IF (my_process_is_mpi_workroot() .AND. var_lists(i)%p%first) THEN

          fileID = cdiIds(i)%file

          ! reset the copies of this file ID
          DO j = 1, nvar_lists
            IF(i /= j .AND. fileID == cdiIds(j)%file) cdiIds(j)%file = CDI_UNDEFID
          END DO

          ! close the file
          CALL cdiIds(i)%closeAndDestroyIds()

          CALL create_restart_file_link(TRIM(var_lists(i)%p%filename), TRIM(var_lists(i)%p%model_type), 0, jg, opt_ndom = opt_ndom)
        ENDIF
        var_lists(i)%p%filename   = ''
      ENDIF
    ENDDO close_all_lists

    CALL message('','')

    ! reset all var list properties related to cdi files

    reset_all_lists: DO i = 1, nvar_lists
      var_lists(i)%p%restart_opened = .FALSE.
      var_lists(i)%p%first          = .FALSE.
    ENDDO reset_all_lists

    private_restart_time = ''

  END SUBROUTINE close_writing_restart_files

  !------------------------------------------------------------------------------------------------

  ! loop over all var_lists for restart

  SUBROUTINE write_restart(p_patch, cdiIds)
    TYPE(t_patch), INTENT(in) :: p_patch
    TYPE(t_restart_cdi_ids), INTENT(INOUT) :: cdiIds(:)

    INTEGER :: i,j
    LOGICAL :: write_info

    IF (my_process_is_mpi_test()) RETURN

    write_info   = .TRUE.

    ! pick up first stream associated with each file

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

        ! loop over all streams associated with the file

        DO j = i, nvar_lists

          ! skip var_list if it does not match the current patch ID
          IF (var_lists(j)%p%patch_id /= p_patch%id) CYCLE

          IF (cdiIds(j)%file == cdiIds(i)%file) THEN


            ! write variables

            CALL write_restart_var_list(var_lists(j), p_patch, cdiIds(i)%file)

          ENDIF
        ENDDO
      ENDIF
    ENDDO
!PR
  CALL message('','Finished Write netCDF2 restart for : '//TRIM(private_restart_time))

  END SUBROUTINE write_restart

  !------------------------------------------------------------------------------------------------

  ! write variables of a list for restart

  SUBROUTINE write_restart_var_list(this_list, p_patch, fileId)
    TYPE (t_var_list) ,INTENT(in) :: this_list
    TYPE(t_patch), TARGET, INTENT(in) :: p_patch
    INTEGER, VALUE :: fileId

    ! variables of derived type used in linked list

    CHARACTER(LEN=*), PARAMETER    :: routine = modname//"::write_restart_var_list"
    TYPE (t_var_metadata), POINTER :: info
    TYPE (t_list_element), POINTER :: element
    TYPE (t_list_element), TARGET  :: start_with

    REAL(wp), ALLOCATABLE :: r_out_wp(:) ! field gathered on I/O processor

    INTEGER :: private_n, lev
    TYPE(t_ptr_2d), ALLOCATABLE :: levelPointers(:)
    TYPE(t_comm_gather_pattern), POINTER :: gather_pattern

    INTEGER :: time_level
    INTEGER :: jg
    LOGICAL :: lskip_timelev, lskip_extra_timelevs

    jg = p_patch%id


    ! Loop over all fields in linked list

    element => start_with
    element%next_list_element => this_list%p%first_list_element

    for_all_list_elements: DO

      element => element%next_list_element
      IF (.NOT.ASSOCIATED(element)) EXIT

      ! retrieve information from actual linked list element

      info => element%field%info


      ! skip this field ?

      IF (.NOT. info%lrestart) CYCLE


      ! skip this field because of wrong time index ?

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

      ! we are committed to writing now
      IF (my_process_is_mpi_workroot()) write (0,*)' ... write ',info%name

      ! ALLOCATE the global array to gather the DATA on the master process
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

      ! get pointers to the local DATA
      CALL getLevelPointers(element%field%info, element%field%r_ptr, levelPointers)

      ! gather the DATA IN the master process AND WRITE it to disk
      DO lev = 1, SIZE(levelPointers)
        CALL exchange_data(in_array = levelPointers(lev)%p, out_array = r_out_wp, gather_pattern = gather_pattern)
        IF(my_process_is_mpi_workroot()) THEN
            SELECT CASE(info%ndims)
                CASE(2)
                    CALL streamWriteVar(fileId, info%cdiVarID, r_out_wp(:), 0)
                CASE(3)
                    CALL streamWriteVarSlice(fileId, info%cdiVarID, lev-1, r_out_wp(:), 0)
                CASE DEFAULT
                    CALL finish(routine, TRIM(int2string(info%ndims))//"d arrays not handled yet.")
            END SELECT
        END IF
      END DO

      ! deallocate temporary global arrays
      DEALLOCATE(r_out_wp)
      ! no deallocation of levelPointers so that the next invocation of getLevelPointers() may reuse the allocation
    END DO for_all_list_elements
  END SUBROUTINE write_restart_var_list

  !------------------------------------------------------------------------------------------------

  ! deallocate module variables

  SUBROUTINE finish_restart()

    IF (my_process_is_mpi_test()) RETURN

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

    nv_grids   = 0
    lrestart_initialised = .FALSE.

  END SUBROUTINE finish_restart

  !------------------------------------------------------------------------------------------------

  SUBROUTINE read_restart_files(p_patch, opt_ndom)

    TYPE(t_patch), INTENT(in) :: p_patch
    INTEGER,       OPTIONAL, INTENT(in) :: opt_ndom
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::read_restart_files"

    TYPE model_search
      CHARACTER(len=8) :: abbreviation
      INTEGER          :: key
    END type model_search
    TYPE(model_search) :: abbreviations(nvar_lists)

    TYPE (t_list_element), POINTER   :: element
    TYPE (t_var_metadata), POINTER   :: info
    CHARACTER(len=80)                :: restart_filename, name
    INTEGER                          :: fileID, vlistID, gridID, zaxisID, taxisID,    &
      &                                 varID, idate, itime, ic, il, n, nfiles, i,   &
      &                                 iret, istat, key, vgrid, gdims(5), nindex,   &
      &                                 nmiss, nvars, root_pe, var_ref_pos, lev
    CHARACTER(len=8)                 :: model_type
    REAL(wp), POINTER                :: r1d(:), rptr2d(:,:), rptr3d(:,:,:)
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
      name = TRIM(restart_filename)//CHAR(0)

      IF (my_process_is_mpi_workroot()) THEN
        WRITE(0,*) "streamOpenRead ", TRIM(restart_filename)

        fileID  = streamOpenRead(TRIM(name))
        IF(fileID < 0) CALL finish(routine, "could not open file '"//TRIM(NAME)//"'")
        vlistID = streamInqVlist(fileID)
        taxisID = vlistInqTaxis(vlistID)

        idate   = taxisInqVdate(taxisID)
        itime   = taxisInqVtime(taxisID)

        WRITE(message_text,'(a,i8.8,a,i6.6,a,a)') &
          'Read restart for : ', idate, 'T', itime, 'Z from ',TRIM(restart_filename)
      END IF

      CALL message('read_restart_files',message_text)

      IF (my_process_is_mpi_workroot()) THEN
        nvars = vlistNvars(vlistID)
      END IF
      CALL p_bcast(nvars, root_pe, comm=p_comm_work)

      for_all_vars: DO varID = 0, (nvars-1)

        IF (my_process_is_mpi_workroot()) THEN
          CALL vlistInqVarName(vlistID, varID, name)
        END IF
        CALL p_bcast(name, root_pe, comm=p_comm_work)

        for_all_lists: DO i = 1, nvar_lists
          ! skip var_list if it does not match the current patch ID
          IF (var_lists(i)%p%patch_id   /= p_patch%id) CYCLE
          IF (var_lists(i)%p%model_type /= model_type) CYCLE

          element => find_list_element(var_lists(i), TRIM(name))
          IF (ASSOCIATED(element)) THEN
            IF (element%field%info%lrestart) THEN

              info => element%field%info

              ! allocate temporary global array on output processor
              ! and gather field from other processors

              NULLIFY(rptr2d)
              NULLIFY(rptr3d)

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

              IF (info%lcontained) THEN
                nindex = info%ncontained
              ELSE
                nindex = 1
              ENDIF


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


              ! deallocate temporary global arrays


              IF (my_process_is_mpi_workroot()) THEN
                WRITE (0,*) ' ... read ',TRIM(element%field%info%name)
              ENDIF
              CYCLE for_all_vars
            ENDIF
          ENDIF
        ENDDO for_all_lists
        CALL message('reading_restart_file','Variable '//TRIM(name)//' not defined.')
      ENDDO for_all_vars

      IF (my_process_is_mpi_workroot())  CALL streamClose(fileID)

    ENDDO for_all_files

    DEALLOCATE (r1d)

    CALL message('','')
    CALL message('',separator)
    CALL message('','')

  END SUBROUTINE read_restart_files
  !-------------------------------------------------------------------------
END MODULE mo_restart
