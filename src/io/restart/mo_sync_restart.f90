!> Module for writing restart files (synchronously)
!!
!! Note: The asynchronous implementation of the restart output can be
!!       found in the module "mo_async_restart"
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
MODULE mo_sync_restart
  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_fortran_tools,         ONLY: t_ptr_2d
  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH, TLEV_NNOW, TLEV_NNOW_RCF, SUCCESS
  USE mo_var_metadata_types,    ONLY: t_var_metadata
  USE mo_linked_list,           ONLY: t_var_list, t_list_element
  USE mo_var_list,              ONLY: nvar_lists, var_lists, get_var_timelevel
  USE mo_cdi,                   ONLY: FILETYPE_NC, FILETYPE_NC2, FILETYPE_NC4, CDI_UNDEFID, COMPRESS_ZIP, &
                                    & streamWriteVarSlice, streamWriteVar
  USE mo_cdi_constants,         ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, &
                                    & ZA_HYBRID, ZA_HYBRID_HALF, ZA_DEPTH_BELOW_LAND, ZA_DEPTH_BELOW_LAND_P1, ZA_DEPTH_RUNOFF_S, &
                                    & ZA_DEPTH_RUNOFF_G, ZA_SNOW, ZA_SNOW_HALF, ZA_TOA, ZA_HEIGHT_2M, ZA_HEIGHT_10M, &
                                    & ZA_LAKE_BOTTOM, ZA_LAKE_BOTTOM_HALF, ZA_MIX_LAYER, ZA_SEDIMENT_BOTTOM_TW_HALF, &
                                    & ZA_DEPTH_BELOW_SEA, ZA_DEPTH_BELOW_SEA_HALF, ZA_GENERIC_ICE, ZA_OCEAN_SEDIMENT, ZA_COUNT
  USE mo_restart_util,          ONLY: t_v_grid, t_restart_cdi_ids, set_vertical_grid, setGeneralRestartAttributes, &
                                    & setDynamicPatchRestartAttributes, setPhysicsRestartAttributes, create_restart_file_link, &
                                    & getRestartFilename
  USE mo_restart_var_data,      ONLY: getLevelPointers, has_valid_time_level
  USE mo_util_string,           ONLY: int2string, separator
  USE mo_restart_attributes,    ONLY: t_RestartAttributeList, RestartAttributeList_make
  USE mo_restart_descriptor,    ONLY: t_RestartDescriptor
  USE mo_restart_patch_description, ONLY: t_restart_patch_description
  USE mo_restart_var_data,      ONLY: t_RestartVarData, createRestartVarData
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

  USE mo_model_domain,          ONLY: t_patch, p_patch
  USE mo_mpi,                   ONLY: my_process_is_mpi_workroot, my_process_is_mpi_test
  USE mo_communication,         ONLY: t_comm_gather_pattern, exchange_data

#ifndef __NO_ICON__OCEAN
  USE mo_ocean_nml,                ONLY: lhamocc
  USE mo_sedmnt,                ONLY: ks, ksp, dzsed
  USE mo_math_utilities,        ONLY: set_zlev
#endif

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'

  PUBLIC :: create_restart_file
  PUBLIC :: t_SyncRestartDescriptor

  ! combine the DATA that describes a patch for restart purposes with the infos required for the asynchronous fetching of the DATA from the compute PEs
  TYPE t_PatchData
    TYPE(t_restart_patch_description) :: description
    TYPE(t_RestartVarData), POINTER :: varData(:)
  CONTAINS
    PROCEDURE :: writeFile => patchData_writeFile
  END TYPE t_PatchData

  TYPE, EXTENDS(t_RestartDescriptor) :: t_SyncRestartDescriptor
    TYPE(t_PatchData), ALLOCATABLE :: patchData(:)
  CONTAINS
    PROCEDURE :: construct => syncRestartDescriptor_construct   ! override
    PROCEDURE :: updatePatch => syncRestartDescriptor_updatePatch   ! override
    PROCEDURE :: writeRestart => syncRestartDescriptor_writeRestart ! override
    PROCEDURE :: destruct => syncRestartDescriptor_destruct ! override

    PROCEDURE, PRIVATE :: defineRestartAttributes => syncRestartDescriptor_defineRestartAttributes
  END TYPE t_SyncRestartDescriptor

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
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_sync_restart'

  !------------------------------------------------------------------------------------------------
CONTAINS

  SUBROUTINE syncRestartDescriptor_construct(me)
    CLASS(t_SyncRestartDescriptor), INTENT(INOUT) :: me
    INTEGER :: jg, error
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':restartDescriptor_construct'

    ! allocate patch data structure
    ALLOCATE(me%patchData(n_dom), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

    ! initialize the patch data structures
    DO jg = 1, n_dom
        ! construct the subobjects
        CALL me%patchData(jg)%description%init(p_patch(jg))
        me%patchData(jg)%varData => createRestartVarData(jg)
    END DO
  END SUBROUTINE syncRestartDescriptor_construct

  SUBROUTINE syncRestartDescriptor_updatePatch(me, patch, opt_pvct, opt_t_elapsed_phy, opt_lcall_phy, opt_sim_time, &
                                        &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth, opt_depth_lnd, &
                                        &opt_nlev_snow, opt_nice_class, opt_ndom)
    CLASS(t_SyncRestartDescriptor), INTENT(INOUT) :: me
    TYPE(t_patch), INTENT(IN) :: patch
    INTEGER, INTENT(IN), OPTIONAL :: opt_depth, opt_depth_lnd, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, &
                                   & opt_nlev_snow, opt_nice_class, opt_ndom
    REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time, opt_pvct(:), opt_t_elapsed_phy(:)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_lcall_phy(:)

    INTEGER :: jg
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartDescriptor_updatePatch"

    jg = patch%id
    IF(jg < 1 .OR. jg > SIZE(me%patchData)) CALL finish(routine, "assertion failed: patch id is out of range")
    IF(me%patchData(jg)%description%id /= jg) CALL finish(routine, "assertion failed: patch id doesn't match its array index")
    CALL me%patchData(jg)%description%update(patch, opt_pvct, opt_t_elapsed_phy, opt_lcall_phy, opt_sim_time, &
                                             &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth, opt_depth_lnd, &
                                             &opt_nlev_snow, opt_nice_class, opt_ndom)
  END SUBROUTINE syncRestartDescriptor_updatePatch

  SUBROUTINE syncRestartDescriptor_defineRestartAttributes(me, restartAttributes, datetime, jstep, opt_output_jfile)
    CLASS(t_SyncRestartDescriptor), TARGET, INTENT(IN) :: me
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    TYPE(t_datetime), INTENT(IN) :: datetime
    INTEGER, VALUE :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

    INTEGER :: i, jg, effectiveDomainCount
    CHARACTER(LEN = 2) :: jgString
    CHARACTER(LEN = :), ALLOCATABLE :: prefix
    TYPE(t_restart_patch_description), POINTER :: curDescription
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":syncRestartDescriptor_defineRestartAttributes"

    ! first the attributes that are independent of the domain
    effectiveDomainCount = 1
    IF(me%patchData(1)%description%l_opt_ndom) effectiveDomainCount = me%patchData(1)%description%opt_ndom
    CALL setGeneralRestartAttributes(restartAttributes, datetime, effectiveDomainCount, jstep, opt_output_jfile)

    ! now the stuff that depends on the domain
    DO jg = 1, n_dom
        curDescription => me%patchData(jg)%description
        jgString = TRIM(int2string(jg, "(i2.2)"))
        CALL setDynamicPatchRestartAttributes(restartAttributes, jg, nold(jg), nnow(jg), nnew(jg), nnow_rcf(jg), nnew_rcf(jg))

        !----------------
        ! additional restart-output for nonhydrostatic model
        IF(curDescription%l_opt_sim_time) CALL restartAttributes%setReal('sim_time_DOM'//jgString, curDescription%opt_sim_time )

        !-------------------------------------------------------------
        ! DR
        ! WORKAROUND FOR FIELDS WHICH NEED TO GO INTO THE RESTART FILE,
        ! BUT SO FAR CANNOT BE HANDELED CORRECTLY BY ADD_VAR OR
        ! SET_RESTART_ATTRIBUTE
        !-------------------------------------------------------------

        IF(curDescription%l_opt_ndyn_substeps) THEN
            CALL restartAttributes%setInteger('ndyn_substeps_DOM'//jgString, curDescription%opt_ndyn_substeps)
        END IF
        IF(curDescription%l_opt_jstep_adv_marchuk_order) THEN
            CALL restartAttributes%setInteger('jstep_adv_marchuk_order_DOM'//jgString, curDescription%opt_jstep_adv_marchuk_order)
        END IF

        IF(ALLOCATED(curDescription%opt_t_elapsed_phy) .AND. ALLOCATED(curDescription%opt_lcall_phy)) THEN
            CALL setPhysicsRestartAttributes(restartAttributes, jg, curDescription%opt_t_elapsed_phy(:), &
                                                                  & curDescription%opt_lcall_phy(:))
        END IF
    END DO
  END SUBROUTINE syncRestartDescriptor_defineRestartAttributes

  SUBROUTINE patchData_writeFile(me, restartAttributes, datetime, jstep, modelType, opt_output_jfile)
    CLASS(t_PatchData), INTENT(IN) :: me
    TYPE(t_RestartAttributeList) :: restartAttributes
    TYPE(t_datetime), INTENT(IN) :: datetime
    INTEGER, INTENT(IN) :: jstep                ! simulation step
    CHARACTER(LEN = *), INTENT(IN) :: modelType
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

    INTEGER :: inlev_soil, inlev_snow, i, nice_class, error
    INTEGER :: ndepth    ! depth of n
    REAL(wp), ALLOCATABLE :: zlevels_full(:), zlevels_half(:)
    CHARACTER(len=MAX_CHAR_LENGTH)  :: string
    TYPE(t_restart_cdi_ids), ALLOCATABLE :: cdiIds(:)

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":patchData_writeFile"

    IF(ALLOCATED(me%description%opt_pvct)) CALL set_restart_vct(me%description%opt_pvct)  ! Vertical coordinate (A's and B's)
    IF(me%description%l_opt_depth_lnd .AND. me%description%opt_depth_lnd > 0) THEN            ! geometrical depth for land module
        !This part is only called if me%description%opt_depth_lnd > 0
        inlev_soil = me%description%opt_depth_lnd
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
    ENDIF
    IF(me%description%l_opt_nlev_snow .AND. me%description%opt_nlev_snow > 0) THEN  ! number of snow levels (multi layer snow model)
        !This part is only called if me%description%opt_nlev_snow > 0
        inlev_snow = me%description%opt_nlev_snow
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
!DR end preliminary fix

    ndepth = 0
    IF(ALLOCATED(me%description%opt_ocean_zheight_cellMiddle)) THEN
      IF(.NOT. ALLOCATED(me%description%opt_ocean_Zheight_CellInterfaces) .OR. .NOT. me%description%l_opt_ocean_Zlevels) THEN
          CALL finish('patchData_writeFile','Ocean level parameteres not complete')
      END IF
      IF (.not. use_ocean_levels) THEN
        ! initialize the ocean levels
        IF (ALLOCATED(private_depth_half))  &
          &  DEALLOCATE(private_depth_half, private_depth_full)
        ALLOCATE(private_depth_half(me%description%opt_ocean_Zlevels+1), private_depth_full(me%description%opt_ocean_Zlevels))
        private_depth_half(:) = me%description%opt_ocean_Zheight_CellInterfaces(:)
        private_depth_full(:) = me%description%opt_ocean_Zheight_CellMiddle(:)
      END IF
      ndepth = me%description%opt_ocean_Zlevels
      use_ocean_levels = .TRUE.
     END IF

    nice_class = 1
    IF(me%description%l_opt_nice_class) nice_class = me%description%opt_nice_class

    CALL init_restart( me%description%n_patch_cells_g,             &! total # of cells
                     & me%description%n_patch_verts_g,             &! total # of vertices
                     & me%description%n_patch_edges_g,             &! total # of cells
                     & me%description%nlev,              &! total # of vertical layers
                     & ndepth,            &! total # of depths below sea
                     & inlev_soil,        &! total # of depths below land (TERRA or JSBACH)
                     & inlev_snow,        &! total # of vertical snow layers (TERRA)
                     & nice_class)         ! total # of ice classes (sea ice)
    ALLOCATE(cdiIds(nvar_lists), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    DO i = 1, SIZE(cdiIds, 1)
        CALL cdiIds(i)%init()
    END DO

    private_restart_time = iso8601(datetime)
    string = getRestartFilename(me%description%base_filename, me%description%id, datetime, modelType)

    CALL open_writing_restart_files(me%description%id, me%description%n_patch_cells_g, me%description%n_patch_verts_g, &
                                   &me%description%n_patch_edges_g, me%description%cell_type, TRIM(string), restartAttributes, &
                                   &cdiIds, datetime)

    CALL write_restart(me%description%id, cdiIds)

    IF(me%description%l_opt_ndom) THEN
        CALL close_writing_restart_files(me%description%id, cdiIds, me%description%opt_ndom)
    ELSE
        CALL close_writing_restart_files(me%description%id, cdiIds)
    END IF
    CALL finish_restart()

    DEALLOCATE(cdiIds)
  END SUBROUTINE patchData_writeFile

  SUBROUTINE syncRestartDescriptor_writeRestart(me, datetime, jstep, modelType, opt_output_jfile)
    CLASS(t_SyncRestartDescriptor), INTENT(INOUT) :: me
    TYPE(t_datetime), INTENT(IN) :: datetime
    INTEGER, INTENT(IN) :: jstep
    CHARACTER(LEN = *), INTENT(IN) :: modelType
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

    INTEGER :: jg
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes

    IF(ltimer) CALL timer_start(timer_write_restart_file)

    restartAttributes => RestartAttributeList_make()
    CALL me%defineRestartAttributes(restartAttributes, datetime, jstep, opt_output_jfile)

    DO jg = 1, n_dom
        CALL me%patchData(jg)%writeFile(restartAttributes, datetime, jstep, modelType, opt_output_jfile)
    END DO

    CALL restartAttributes%destruct()
    DEALLOCATE(restartAttributes)

    IF(ltimer) CALL timer_stop(timer_write_restart_file)
  END SUBROUTINE syncRestartDescriptor_writeRestart

  SUBROUTINE syncRestartDescriptor_destruct(me)
    CLASS(t_SyncRestartDescriptor), INTENT(INOUT) :: me

    DEALLOCATE(me%patchData)
  END SUBROUTINE syncRestartDescriptor_destruct

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
  SUBROUTINE open_writing_restart_files(domain, cellCount, vertCount, edgeCount, cellType, restart_filename, restartAttributes, &
                                       &cdiIds, datetime)
    INTEGER, VALUE :: domain, cellCount, vertCount, edgeCount, cellType  ! TODO[NH]: pass a t_restart_patch_description instead
    CHARACTER(LEN=*), INTENT(IN) :: restart_filename
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    TYPE(t_restart_cdi_ids), INTENT(INOUT) :: cdiIds(:)
    TYPE(t_datetime), INTENT(IN) :: datetime

    INTEGER :: status, i ,j, ia, ihg, ivg

    CHARACTER(len=64) :: attribute_name
    CHARACTER(len=256) :: text_attribute
    REAL(wp) :: real_attribute(1)
    INTEGER :: int_attribute(1), axisIds(ZA_COUNT)
    LOGICAL :: bool_attribute
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":open_writing_restart_files"

    IF (my_process_is_mpi_test()) RETURN

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
      IF (var_lists(i)%p%patch_id /= domain) CYCLE

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
                                        &cellCount, vertCount, edgeCount, cellType, vgrid_def(1:nv_grids), private_vct)
        ELSE
            CALL cdiIds(i)%openRestartAndCreateIds(TRIM(restart_filename), var_lists(i)%p%restart_type, restartAttributes, &
                                        &cellCount, vertCount, edgeCount, cellType, vgrid_def(1:nv_grids))
        END IF
        ! set the related fields IN the var_lists

        var_lists(i)%p%filename = TRIM(restart_filename)
        var_lists(i)%p%restart_opened = .TRUE.

        ! 6. add variables

        CALL addVarListToVlist(var_lists(i), cdiIds(i), domain)

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
        IF (var_lists(j)%p%patch_id /= domain) CYCLE

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

            CALL addVarListToVlist(var_lists(j), cdiIds(i), domain)

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

    TYPE (t_list_element), POINTER :: element
    TYPE (t_list_element), TARGET  :: start_with

    INTEGER :: time_level
    LOGICAL :: lskip_timelev, lskip_extra_timelevs

    element => start_with
    element%next_list_element => this_list%p%first_list_element

    DO
        element => element%next_list_element
        IF(.NOT.ASSOCIATED(element)) EXIT

        IF(has_valid_time_level(element%field%info, jg, nnew(jg), nnew_rcf(jg))) CALL cdiIds%defineVariable(element%field%info)
    END DO
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

    CALL open_writing_restart_files(jg, kcell, kvert, kedge, patch%geometry_info%cell_type, TRIM(string), restartAttributes, &
                                   &cdiIds, datetime)

    CALL write_restart(jg, cdiIds)

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

  SUBROUTINE write_restart(domain, cdiIds)
    INTEGER, VALUE :: domain
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
          IF (var_lists(j)%p%patch_id /= domain) CYCLE

          IF (cdiIds(j)%file == cdiIds(i)%file) THEN


            ! write variables

            CALL write_restart_var_list(var_lists(j), domain, cdiIds(i)%file)

          ENDIF
        ENDDO
      ENDIF
    ENDDO
!PR
  CALL message('','Finished Write netCDF2 restart for : '//TRIM(private_restart_time))

  END SUBROUTINE write_restart

  !------------------------------------------------------------------------------------------------

  ! write variables of a list for restart

  SUBROUTINE write_restart_var_list(this_list, domain, fileId)
    TYPE (t_var_list) ,INTENT(in) :: this_list
    INTEGER, VALUE :: domain, fileId

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
    LOGICAL :: lskip_timelev, lskip_extra_timelevs

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
      time_level = get_var_timelevel(element%field%info)

      lskip_timelev = .FALSE.
      lskip_extra_timelevs = .FALSE.
#ifndef __NO_ICON_ATMO__
      IF (iequations == INH_ATMOSPHERE .AND. .NOT. (l_limited_area .AND. domain == 1)) THEN
        lskip_extra_timelevs = .TRUE.
      ENDIF

      ! get information about timelevel to be skipped for current field
      IF (element%field%info%tlev_source == TLEV_NNOW ) THEN
        IF (time_level == nnew(domain))                    lskip_timelev = .TRUE.
        ! this is needed to skip the extra time levels allocated for nesting
        IF (lskip_extra_timelevs .AND. time_level > 2) lskip_timelev = .TRUE.
      ELSE IF (element%field%info%tlev_source == TLEV_NNOW_RCF) THEN
        IF (time_level == nnew_rcf(domain))  lskip_timelev = .TRUE.
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
        gather_pattern => p_patch(domain)%comm_pat_gather_c
      CASE (GRID_UNSTRUCTURED_VERT)
        private_n = private_nv
        gather_pattern => p_patch(domain)%comm_pat_gather_v
      CASE (GRID_UNSTRUCTURED_EDGE)
        private_n = private_ne
        gather_pattern => p_patch(domain)%comm_pat_gather_e
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

END MODULE mo_sync_restart
