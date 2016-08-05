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
                              & ZA_DEPTH_BELOW_SEA, ZA_DEPTH_BELOW_SEA_HALF, ZA_COUNT
    USE mo_dynamics_config, ONLY: nold, nnow, nnew, nnew_rcf, nnow_rcf
    USE mo_exception, ONLY: finish
    USE mo_fortran_tools, ONLY: assign_if_present, assign_if_present_allocatable
    USE mo_io_units, ONLY: filename_max
    USE mo_kind, ONLY: wp
    USE mo_model_domain, ONLY: t_patch
    USE mo_mpi, ONLY: p_pe_work
    USE mo_packed_message, ONLY: t_PackedMessage
    USE mo_util_restart, ONLY: t_v_grid, set_vertical_grid

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: t_restart_patch_description

    ! TYPE t_restart_patch_description contains all the DATA that describes a patch for restart purposes
    TYPE t_restart_patch_description
        ! vertical grid definitions
        TYPE(t_v_grid) :: v_grid_defs(ZA_COUNT)
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
        LOGICAL :: l_opt_depth
        INTEGER :: opt_depth
        LOGICAL :: l_opt_depth_lnd
        INTEGER :: opt_depth_lnd
        LOGICAL :: l_opt_nlev_snow
        INTEGER :: opt_nlev_snow
        LOGICAL :: l_opt_nice_class
        INTEGER :: opt_nice_class
        LOGICAL :: l_opt_ndyn_substeps
        INTEGER :: opt_ndyn_substeps
        LOGICAL :: l_opt_jstep_adv_marchuk_order
        INTEGER :: opt_jstep_adv_marchuk_order
        LOGICAL :: l_opt_sim_time
        REAL(wp) :: opt_sim_time
        LOGICAL :: l_opt_ndom
        INTEGER :: opt_ndom

        REAL(wp), ALLOCATABLE :: opt_pvct(:)
        LOGICAL, ALLOCATABLE :: opt_lcall_phy(:)
        REAL(wp), ALLOCATABLE :: opt_t_elapsed_phy(:)
    CONTAINS
        PROCEDURE :: init => restartPatchDescription_init
        PROCEDURE :: update => restartPatchDescription_update   ! called to set the DATA that may change from restart to restart
        PROCEDURE :: setTimeLevels => restartPatchDescription_setTimeLevels ! set the time level fields (nold, ...) to match the respective global variables
        PROCEDURE :: packer => restartPatchDescription_packer
        PROCEDURE :: defineVGrids => restartPatchDescription_defineVGrids
    END TYPE t_restart_patch_description

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_patch_description"

CONTAINS

    SUBROUTINE restartPatchDescription_init(me, p_patch)
        CLASS(t_restart_patch_description), INTENT(INOUT) :: me
        TYPE(t_patch), INTENT(IN) :: p_patch

        me%id = p_patch%id
        me%work_pe0_id = p_patch%proc0
        me%nlev = p_patch%nlev
        me%cell_type = p_patch%geometry_info%cell_type
        me%base_filename = TRIM(p_patch%grid_filename)
        me%n_patch_cells_g = p_patch%n_patch_cells_g
        me%n_patch_verts_g = p_patch%n_patch_verts_g
        me%n_patch_edges_g = p_patch%n_patch_edges_g
        me%v_grid_count = 0
        me%l_dom_active = .FALSE.
        me%restart_proc_id = 0
        me%nold = nold(p_patch%id)
        me%nnow = nnow(p_patch%id)
        me%nnew = nnew(p_patch%id)
        me%nnew_rcf = nnew_rcf(p_patch%id)
        me%nnow_rcf = nnow_rcf(p_patch%id)
        me%l_opt_depth = .FALSE.
        me%opt_depth = 0
        me%l_opt_depth_lnd = .FALSE.
        me%opt_depth_lnd = 0
        me%l_opt_nlev_snow = .FALSE.
        me%opt_nlev_snow = 0
        me%l_opt_nice_class = .FALSE.
        me%opt_nice_class = 0
        me%l_opt_ndyn_substeps = .FALSE.
        me%opt_ndyn_substeps = 0
        me%l_opt_jstep_adv_marchuk_order = .FALSE.
        me%opt_jstep_adv_marchuk_order = 0
        me%l_opt_ndom = .FALSE.
        me%opt_ndom = 0
        me%l_opt_sim_time = .FALSE.
        me%opt_sim_time = 0.
    END SUBROUTINE restartPatchDescription_init

    SUBROUTINE restartPatchDescription_update(me, patch, opt_pvct, opt_t_elapsed_phy, opt_lcall_phy, opt_sim_time, &
                                             &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth, opt_depth_lnd, &
                                             &opt_nlev_snow, opt_nice_class, opt_ndom)
        CLASS(t_restart_patch_description), INTENT(INOUT) :: me
        TYPE(t_patch), INTENT(IN) :: patch
        INTEGER, INTENT(IN), OPTIONAL :: opt_depth, opt_depth_lnd, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, &
                                       & opt_nlev_snow, opt_nice_class, opt_ndom
        REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time, opt_pvct(:), opt_t_elapsed_phy(:)
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

            CALL assign_if_present(me%opt_ndyn_substeps, opt_ndyn_substeps)
            me%l_opt_ndyn_substeps = PRESENT(opt_ndyn_substeps)

            CALL assign_if_present(me%opt_jstep_adv_marchuk_order, opt_jstep_adv_marchuk_order)
            me%l_opt_jstep_adv_marchuk_order = PRESENT(opt_jstep_adv_marchuk_order)

            CALL assign_if_present(me%opt_depth, opt_depth)
            me%l_opt_depth = PRESENT(opt_depth)

            CALL assign_if_present(me%opt_depth_lnd, opt_depth_lnd)
            me%l_opt_depth_lnd = PRESENT(opt_depth_lnd)

            CALL assign_if_present(me%opt_nlev_snow, opt_nlev_snow)
            me%l_opt_nlev_snow = PRESENT(opt_nlev_snow)

            CALL assign_if_present(me%opt_nice_class, opt_nice_class)
            me%l_opt_nice_class = PRESENT(opt_nice_class)

            CALL assign_if_present(me%opt_sim_time, opt_sim_time)
            me%l_opt_sim_time = PRESENT(opt_sim_time)

            CALL assign_if_present(me%opt_ndom, opt_ndom)
            me%l_opt_ndom = PRESENT(opt_ndom)
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
        CALL message%execute(operation, me%id)
        CALL message%execute(operation, me%l_dom_active)
        CALL message%execute(operation, me%nlev)
        CALL message%execute(operation, me%cell_type)
        CALL message%execute(operation, me%n_patch_cells_g)
        CALL message%execute(operation, me%n_patch_edges_g)
        CALL message%execute(operation, me%n_patch_verts_g)
        CALL message%execute(operation, me%restart_proc_id)
        CALL message%execute(operation, me%work_pe0_id)
        CALL message%execute(operation, me%base_filename)

        ! time levels
        CALL message%execute(operation, me%nold)
        CALL message%execute(operation, me%nnow)
        CALL message%execute(operation, me%nnow_rcf)
        CALL message%execute(operation, me%nnew)
        CALL message%execute(operation, me%nnew_rcf)

        ! optional parameter values
        CALL message%execute(operation, me%l_opt_depth)
        CALL message%execute(operation, me%opt_depth)
        CALL message%execute(operation, me%l_opt_depth_lnd)
        CALL message%execute(operation, me%opt_depth_lnd)
        CALL message%execute(operation, me%l_opt_nlev_snow)
        CALL message%execute(operation, me%opt_nlev_snow)
        CALL message%execute(operation, me%l_opt_nice_class)
        CALL message%execute(operation, me%opt_nice_class)
        CALL message%execute(operation, me%l_opt_ndyn_substeps)
        CALL message%execute(operation, me%opt_ndyn_substeps)
        CALL message%execute(operation, me%l_opt_jstep_adv_marchuk_order)
        CALL message%execute(operation, me%opt_jstep_adv_marchuk_order)
        CALL message%execute(operation, me%l_opt_sim_time)
        CALL message%execute(operation, me%opt_sim_time)
        CALL message%execute(operation, me%l_opt_ndom)
        CALL message%execute(operation, me%opt_ndom)

        ! optional parameter arrays
        CALL message%execute(operation, me%opt_pvct)
        CALL message%execute(operation, me%opt_lcall_phy)
        CALL message%execute(operation, me%opt_t_elapsed_phy)
    END SUBROUTINE restartPatchDescription_packer

    !  Set vertical grid definition.
    SUBROUTINE restartPatchDescription_defineVGrids(me)
        CLASS(t_restart_patch_description), TARGET, INTENT(INOUT) :: me

        INTEGER :: nlev_soil, nlev_snow, nlev_ocean, nice_class, ierrstat

        ! DEFAULT values for the level counts
        nlev_soil = 0
        nlev_snow = 0
        nlev_ocean = 0
        nice_class = 1

        ! replace DEFAULT values by the overrides provided IN the me
        IF(me%l_opt_depth_lnd) nlev_soil = me%opt_depth_lnd
        IF(me%l_opt_nlev_snow) nlev_snow = me%opt_nlev_snow
        IF(me%l_opt_depth) nlev_ocean = me%opt_depth
        IF(me%l_opt_nice_class) nice_class = me%opt_nice_class

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
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_GENERIC_ICE, 1._wp)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_DEPTH_RUNOFF_S, 1)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_DEPTH_RUNOFF_G, 1)
        IF(me%l_opt_depth_lnd) CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_DEPTH_BELOW_LAND, nlev_soil)
        IF(me%l_opt_depth_lnd) CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_DEPTH_BELOW_LAND_P1, nlev_soil+1)
        IF(me%l_opt_nlev_snow) CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_SNOW, nlev_snow)
        IF(me%l_opt_nlev_snow) CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_SNOW_HALF, nlev_snow+1)
        IF(me%l_opt_depth) CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_DEPTH_BELOW_SEA, nlev_ocean)
        IF(me%l_opt_depth) CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_DEPTH_BELOW_SEA_HALF, nlev_ocean+1)
    END SUBROUTINE restartPatchDescription_defineVGrids

END MODULE mo_restart_patch_description
