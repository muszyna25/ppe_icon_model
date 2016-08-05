!>
!! A CLASS that bundles all a patch's metadata together which IS relevant for restart purposes.
!! ----------------------------------------------------------
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_restart_patch_description
    USE mo_cdi_constants, ONLY: ZA_SURFACE, ZA_HYBRID, ZA_HYBRID_HALF, ZA_HEIGHT_2M, ZA_HEIGHT_10M, ZA_TOA, ZA_LAKE_BOTTOM, &
                              & ZA_MIX_LAYER, ZA_LAKE_BOTTOM_HALF, ZA_SEDIMENT_BOTTOM_TW_HALF, ZA_GENERIC_ICE, ZA_DEPTH_RUNOFF_S, &
                              & ZA_DEPTH_RUNOFF_G, ZA_DEPTH_BELOW_LAND, ZA_DEPTH_BELOW_LAND_P1, ZA_SNOW, ZA_SNOW_HALF, &
                              & ZA_DEPTH_BELOW_SEA, ZA_DEPTH_BELOW_SEA_HALF, ZA_OCEAN_SEDIMENT, ZA_COUNT, GRID_UNSTRUCTURED_CELL, &
                              & GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE
    USE mo_cdi_ids, ONLY: t_CdiIds, set_vertical_grid, t_Vgrid
    USE mo_communication, ONLY: t_comm_gather_pattern
    USE mo_dynamics_config, ONLY: nold, nnow, nnew, nnew_rcf, nnow_rcf
    USE mo_exception, ONLY: finish
    USE mo_fortran_tools, ONLY: assign_if_present_allocatable
    USE mo_restart_attributes, ONLY: t_RestartAttributeList
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_io_units, ONLY: filename_max
    USE mo_kind, ONLY: wp
    USE mo_model_domain, ONLY: p_patch, t_patch
    USE mo_mpi, ONLY: p_pe_work, my_process_is_work
    USE mo_packed_message, ONLY: t_PackedMessage
    USE mo_restart_util, ONLY: setDynamicPatchRestartAttributes, setPhysicsRestartAttributes
    USE mo_util_string, ONLY: int2string

#ifndef __NO_ICON__OCEAN
    USE mo_ocean_nml, ONLY: lhamocc
    USE mo_sedmnt, ONLY: ks, dzsed
    USE mo_math_utilities, ONLY: set_zlev
#endif

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: t_restart_patch_description

    ! TYPE t_restart_patch_description contains all the DATA that describes a patch for restart purposes
    TYPE t_restart_patch_description
        ! vertical grid definitions
        TYPE(t_Vgrid) :: v_grid_defs(ZA_COUNT)
        INTEGER :: v_grid_count

        ! logical patch id
        INTEGER :: id

        ! current model domain activity flag
        LOGICAL :: l_dom_active

        ! number of full levels
        INTEGER :: nlev

        ! cell type
        INTEGER :: cell_type

        ! total # of cells, # of vertices per cell
        INTEGER :: n_patch_cells_g
        ! total # of cells, shape of control volume for edge
        INTEGER :: n_patch_edges_g
        ! total # of vertices, # of vertices per dual cell
        INTEGER :: n_patch_verts_g

        ! process id
        INTEGER :: restart_proc_id

        ! id of PE0 of working group (/= 0 in case of processor splitting)
        INTEGER :: work_pe0_id

        ! base file name contains already logical patch ident
        CHARACTER(LEN = filename_max) :: base_filename

        ! dynamic patch arguments (mandatory)
        INTEGER :: nold,nnow,nnew,nnew_rcf,nnow_rcf

        ! dynamic patch arguments (optionally)
        INTEGER, ALLOCATABLE :: opt_depth_lnd, opt_nlev_snow, opt_nice_class, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, &
                              & opt_ndom, opt_ocean_zlevels
        REAL(wp), ALLOCATABLE :: opt_sim_time

        REAL(wp), ALLOCATABLE :: opt_pvct(:)
        LOGICAL, ALLOCATABLE :: opt_lcall_phy(:)
        REAL(wp), ALLOCATABLE :: opt_t_elapsed_phy(:)
        REAL(wp), ALLOCATABLE :: opt_ocean_zheight_cellMiddle(:)
        REAL(wp), ALLOCATABLE :: opt_ocean_zheight_cellInterfaces(:)

        ! these are used for synchronous restart writing
        TYPE(t_comm_gather_pattern), POINTER :: cellGatherPattern, vertGatherPattern, edgeGatherPattern
    CONTAINS
        PROCEDURE :: init => restartPatchDescription_init
        PROCEDURE :: update => restartPatchDescription_update   ! called to set the DATA that may change from restart to restart
        PROCEDURE :: setTimeLevels => restartPatchDescription_setTimeLevels ! set the time level fields (nold, ...) to match the respective global variables
        PROCEDURE :: packer => restartPatchDescription_packer
        PROCEDURE :: defineVGrids => restartPatchDescription_defineVGrids
        PROCEDURE :: setRestartAttributes => restartPatchDescription_setRestartAttributes
        PROCEDURE :: getGatherPattern => restartPatchDescription_getGatherPattern
        PROCEDURE :: getGlobalGridSize => restartPatchDescription_getGlobalGridSize
    END TYPE t_restart_patch_description

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_patch_description"

CONTAINS

    SUBROUTINE restartPatchDescription_init(me, domain)
        CLASS(t_restart_patch_description), INTENT(INOUT) :: me
        INTEGER, VALUE :: domain

        ! DEFAULT initialization of all variables
        me%id = domain
        me%work_pe0_id = -1
        me%nlev = -1
        me%cell_type = -1
        me%base_filename = ''
        me%n_patch_cells_g = -1
        me%n_patch_verts_g = -1
        me%n_patch_edges_g = -1
        me%v_grid_count = 0
        me%l_dom_active = .FALSE.
        me%restart_proc_id = 0
        me%nold = nold(domain)
        me%nnow = nnow(domain)
        me%nnew = nnew(domain)
        me%nnew_rcf = nnew_rcf(domain)
        me%nnow_rcf = nnow_rcf(domain)
        me%cellGatherPattern => NULL()
        me%vertGatherPattern => NULL()
        me%edgeGatherPattern => NULL()

        ! patch dependent info, p_patch IS NOT available on restart PEs
        IF(my_process_is_work()) THEN
            me%work_pe0_id = p_patch(domain)%proc0
            me%nlev = p_patch(domain)%nlev
            me%cell_type = p_patch(domain)%geometry_info%cell_type
            me%base_filename = TRIM(p_patch(domain)%grid_filename)
            me%n_patch_cells_g = p_patch(domain)%n_patch_cells_g
            me%n_patch_verts_g = p_patch(domain)%n_patch_verts_g
            me%n_patch_edges_g = p_patch(domain)%n_patch_edges_g
            me%cellGatherPattern => p_patch(domain)%comm_pat_gather_c
            me%vertGatherPattern => p_patch(domain)%comm_pat_gather_v
            me%edgeGatherPattern => p_patch(domain)%comm_pat_gather_e
        END IF
    END SUBROUTINE restartPatchDescription_init

    SUBROUTINE restartPatchDescription_update(me, patch, opt_pvct, opt_t_elapsed_phy, opt_lcall_phy, opt_sim_time, &
                                             &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth_lnd, &
                                             &opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels, &
                                             &opt_ocean_zheight_cellMiddle, opt_ocean_zheight_cellInterfaces)
        CLASS(t_restart_patch_description), INTENT(INOUT) :: me
        TYPE(t_patch), INTENT(IN) :: patch
        INTEGER, INTENT(IN), OPTIONAL :: opt_depth_lnd, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, &
                                       & opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels
        REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time, opt_pvct(:), opt_t_elapsed_phy(:), opt_ocean_zheight_cellMiddle(:), &
                                        & opt_ocean_zheight_cellInterfaces(:)
        LOGICAL, INTENT(IN), OPTIONAL :: opt_lcall_phy(:)

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartPatchDescription_update"

        IF(me%id /= patch%id) CALL finish(routine, "assertion failed: wrong patch passed to update()")

        ! update activity flag - this needs to be done on all compute PEs, otherwise starting a nest during runtime would lead to incomplete restart files
        me%l_dom_active = patch%ldom_active

        ! otherwise, only the patch master process needs the dynamic restart arguments, these will be communicated to the other processes anyway
        IF(p_pe_work == 0 .OR. p_pe_work == me%work_pe0_id) THEN
            ! Patch-dependent attributes
            CALL assign_if_present_allocatable(me%opt_pvct, opt_pvct)
            CALL assign_if_present_allocatable(me%opt_t_elapsed_phy, opt_t_elapsed_phy)
            CALL assign_if_present_allocatable(me%opt_lcall_phy, opt_lcall_phy)
            CALL assign_if_present_allocatable(me%opt_ocean_zheight_cellMiddle, opt_ocean_zheight_cellMiddle)
            CALL assign_if_present_allocatable(me%opt_ocean_zheight_cellInterfaces, opt_ocean_zheight_cellInterfaces)
            CALL assign_if_present_allocatable(me%opt_ndyn_substeps, opt_ndyn_substeps)
            CALL assign_if_present_allocatable(me%opt_jstep_adv_marchuk_order, opt_jstep_adv_marchuk_order)
            CALL assign_if_present_allocatable(me%opt_depth_lnd, opt_depth_lnd)
            CALL assign_if_present_allocatable(me%opt_nlev_snow, opt_nlev_snow)
            CALL assign_if_present_allocatable(me%opt_nice_class, opt_nice_class)
            CALL assign_if_present_allocatable(me%opt_sim_time, opt_sim_time)
            CALL assign_if_present_allocatable(me%opt_ndom, opt_ndom)
            CALL assign_if_present_allocatable(me%opt_ocean_zlevels, opt_ocean_zlevels)

            ! consistency check for OPTIONAL ocean variables
            IF(ALLOCATED(me%opt_ocean_zheight_cellMiddle)) THEN
                IF(.NOT. ALLOCATED(me%opt_ocean_Zheight_CellInterfaces) .OR. .NOT. ALLOCATED(me%opt_ocean_Zlevels)) THEN
                    CALL finish(routine, 'Ocean level parameteres not complete')
                END IF
            END IF
        END IF
    END SUBROUTINE restartPatchDescription_update

    SUBROUTINE restartPatchDescription_setTimeLevels(me)
        CLASS(t_restart_patch_description), INTENT(INOUT) :: me

        me%nold = nold(me%id)
        me%nnow = nnow(me%id)
        me%nnow_rcf = nnow_rcf(me%id)
        me%nnew = nnew(me%id)
        me%nnew_rcf = nnew_rcf(me%id)
    END SUBROUTINE restartPatchDescription_setTimeLevels

    SUBROUTINE restartPatchDescription_packer(me, operation, message)
        INTEGER, VALUE :: operation
        CLASS(t_restart_patch_description), INTENT(INOUT) :: me
        TYPE(t_PackedMessage), INTENT(INOUT) :: message

        ! patch id AND activity flag
        CALL message%packer(operation, me%id)
        CALL message%packer(operation, me%l_dom_active)
        CALL message%packer(operation, me%nlev)
        CALL message%packer(operation, me%cell_type)
        CALL message%packer(operation, me%n_patch_cells_g)
        CALL message%packer(operation, me%n_patch_edges_g)
        CALL message%packer(operation, me%n_patch_verts_g)
        CALL message%packer(operation, me%restart_proc_id)
        CALL message%packer(operation, me%work_pe0_id)
        CALL message%packer(operation, me%base_filename)

        ! time levels
        CALL message%packer(operation, me%nold)
        CALL message%packer(operation, me%nnow)
        CALL message%packer(operation, me%nnow_rcf)
        CALL message%packer(operation, me%nnew)
        CALL message%packer(operation, me%nnew_rcf)

        ! optional parameter values
        CALL message%packerAllocatable(operation, me%opt_depth_lnd)
        CALL message%packerAllocatable(operation, me%opt_nlev_snow)
        CALL message%packerAllocatable(operation, me%opt_nice_class)
        CALL message%packerAllocatable(operation, me%opt_ndyn_substeps)
        CALL message%packerAllocatable(operation, me%opt_jstep_adv_marchuk_order)
        CALL message%packerAllocatable(operation, me%opt_sim_time)
        CALL message%packerAllocatable(operation, me%opt_ndom)
        CALL message%packerAllocatable(operation, me%opt_ocean_zlevels)

        ! optional parameter arrays
        CALL message%packer(operation, me%opt_pvct)
        CALL message%packer(operation, me%opt_lcall_phy)
        CALL message%packer(operation, me%opt_t_elapsed_phy)
    END SUBROUTINE restartPatchDescription_packer

    !  Set vertical grid definition.
    SUBROUTINE restartPatchDescription_defineVGrids(me)
        CLASS(t_restart_patch_description), TARGET, INTENT(INOUT) :: me

        INTEGER :: nlev_soil, nlev_snow, nlev_ocean, nice_class
        INTEGER :: error
        REAL(wp), ALLOCATABLE :: levels(:), levels_sp(:)
        CHARACTER(*), PARAMETER :: routine = modname//":restartPatchDescription_defineVGrids"

        ! DEFAULT values for the level counts
        nlev_soil = 0
        nlev_snow = 0
        nlev_ocean = 0
        nice_class = 1

        ! replace DEFAULT values by the overrides provided IN the me
        IF(ALLOCATED(me%opt_depth_lnd)) nlev_soil = me%opt_depth_lnd
        IF(ALLOCATED(me%opt_nlev_snow)) nlev_snow = me%opt_nlev_snow
        IF(ALLOCATED(me%opt_ocean_zlevels)) nlev_ocean = me%opt_ocean_zlevels
        IF(ALLOCATED(me%opt_nice_class)) nice_class = me%opt_nice_class

        ! set vertical grid definitions
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_SURFACE, 0._wp)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_HYBRID, me%nlev)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_HYBRID_HALF, me%nlev+1)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_HEIGHT_2M, 2._wp)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_HEIGHT_10M, 10._wp)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_TOA, 1._wp)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_LAKE_BOTTOM, 1._wp)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_MIX_LAYER, 1._wp)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_LAKE_BOTTOM_HALF, 1._wp)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_SEDIMENT_BOTTOM_TW_HALF, 0._wp)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_GENERIC_ICE, nice_class)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_DEPTH_RUNOFF_S, 1)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_DEPTH_RUNOFF_G, 1)
        IF(ALLOCATED(me%opt_depth_lnd)) THEN
            CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_DEPTH_BELOW_LAND, nlev_soil)
            CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_DEPTH_BELOW_LAND_P1, nlev_soil+1)
        END IF
        IF(ALLOCATED(me%opt_nlev_snow)) THEN
            CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_SNOW, nlev_snow)
            CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_SNOW_HALF, nlev_snow+1)
        END IF
        IF(ALLOCATED(me%opt_ocean_zlevels)) THEN
            CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_DEPTH_BELOW_SEA, nlev_ocean)
            CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_DEPTH_BELOW_SEA_HALF, nlev_ocean+1)
        END IF
#ifndef __NO_ICON__OCEAN
        IF(lhamocc) THEN
            ! HAMOCC sediment
            ALLOCATE(levels(ks), levels_sp(ks + 1), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

            CALL set_zlev(levels_sp, levels, ks, dzsed*1000._wp)
            CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_OCEAN_SEDIMENT, REAL(levels,wp))

            DEALLOCATE(levels, levels_sp)
        END IF
#endif
    END SUBROUTINE restartPatchDescription_defineVGrids

    SUBROUTINE restartPatchDescription_setRestartAttributes(me, restartAttributes)
        CLASS(t_restart_patch_description), INTENT(IN) :: me
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes

        CHARACTER(LEN = 2) :: domainString

        domainString = TRIM(int2string(me%id, '(i2.2)'))

        ! set time levels
        CALL setDynamicPatchRestartAttributes(restartAttributes, me%id, me%nold, me%nnow, me%nnew, me%nnow_rcf, me%nnew_rcf)

        ! additional restart-output for nonhydrostatic model
        IF(ALLOCATED(me%opt_sim_time)) CALL restartAttributes%setReal ('sim_time_DOM'//domainString, me%opt_sim_time)

        !-------------------------------------------------------------
        ! DR
        ! WORKAROUND FOR FIELDS WHICH NEED TO GO INTO THE RESTART FILE,
        ! BUT SO FAR CANNOT BE HANDELED CORRECTLY BY ADD_VAR OR
        ! SET_RESTART_ATTRIBUTE
        !-------------------------------------------------------------
        IF(ALLOCATED(me%opt_ndyn_substeps)) THEN
            CALL restartAttributes%setInteger('ndyn_substeps_DOM'//domainString, me%opt_ndyn_substeps)
        END IF
        IF(ALLOCATED(me%opt_jstep_adv_marchuk_order)) THEN
            CALL restartAttributes%setInteger('jstep_adv_marchuk_order_DOM'//domainString, me%opt_jstep_adv_marchuk_order)
        END IF

        IF (ALLOCATED(me%opt_t_elapsed_phy) .AND. ALLOCATED(me%opt_lcall_phy)) THEN
            CALL setPhysicsRestartAttributes(restartAttributes, me%id, me%opt_t_elapsed_phy, me%opt_lcall_phy)
        END IF
    END SUBROUTINE restartPatchDescription_setRestartAttributes

    FUNCTION restartPatchDescription_getGatherPattern(me, gridType) RESULT(RESULT)
        CLASS(t_restart_patch_description), INTENT(IN) :: me
        INTEGER, VALUE :: gridType
        TYPE(t_comm_gather_pattern), POINTER :: RESULT

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartPatchDescription_getGatherPattern"

        SELECT CASE(gridType)
            CASE(GRID_UNSTRUCTURED_CELL)
                RESULT => me%cellGatherPattern
            CASE(GRID_UNSTRUCTURED_VERT)
                RESULT => me%vertGatherPattern
            CASE(GRID_UNSTRUCTURED_EDGE)
                RESULT => me%edgeGatherPattern
            CASE DEFAULT
                CALL finish(routine, "illegal gridType argument")
        END SELECT
    END FUNCTION restartPatchDescription_getGatherPattern

    INTEGER FUNCTION restartPatchDescription_getGlobalGridSize(me, gridType) RESULT(RESULT)
        CLASS(t_restart_patch_description), INTENT(IN) :: me
        INTEGER, VALUE :: gridType

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartPatchDescription_getGlobalGridSize"

        SELECT CASE(gridType)
            CASE(GRID_UNSTRUCTURED_CELL)
                RESULT = me%n_patch_cells_g
            CASE(GRID_UNSTRUCTURED_VERT)
                RESULT = me%n_patch_verts_g
            CASE(GRID_UNSTRUCTURED_EDGE)
                RESULT = me%n_patch_edges_g
            CASE DEFAULT
                CALL finish(routine, "illegal gridType argument")
        END SELECT
    END FUNCTION restartPatchDescription_getGlobalGridSize

END MODULE mo_restart_patch_description
