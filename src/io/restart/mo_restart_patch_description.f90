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
    USE mo_zaxis_type, ONLY: ZA_SURFACE, ZA_REFERENCE, ZA_REFERENCE_HALF, ZA_HEIGHT_2M, ZA_HEIGHT_10M, ZA_TOA, ZA_LAKE_BOTTOM, &
                              & ZA_MIX_LAYER, ZA_LAKE_BOTTOM_HALF, ZA_SEDIMENT_BOTTOM_TW_HALF, ZA_GENERIC_ICE, ZA_DEPTH_RUNOFF_S, &
                              & ZA_DEPTH_RUNOFF_G, ZA_DEPTH_BELOW_LAND, ZA_DEPTH_BELOW_LAND_P1, ZA_SNOW, ZA_SNOW_HALF, &
                              & ZA_DEPTH_BELOW_SEA, ZA_DEPTH_BELOW_SEA_HALF, ZA_OCEAN_SEDIMENT, ZA_TROPOPAUSE, zaxisTypeList
    USE mo_cdi_ids, ONLY: set_vertical_grid, t_Vgrid
    USE mo_communication, ONLY: t_comm_gather_pattern
    USE mo_dynamics_config, ONLY: nold, nnow, nnew, nnew_rcf, nnow_rcf
    USE mo_exception, ONLY: finish
    USE mo_fortran_tools, ONLY: assign_if_present_allocatable, assign_if_present
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE
    USE mo_io_units, ONLY: filename_max
    USE mo_kind, ONLY: wp
    USE mo_model_domain, ONLY: p_patch, t_patch
    USE mo_mpi, ONLY: p_pe_work, my_process_is_work, process_mpi_all_comm, p_pe_work, process_mpi_root_id, &
                    & my_process_is_mpi_workroot
    USE mo_packed_message, ONLY: t_PackedMessage, kPackOp, kUnpackOp
    USE mo_upatmo_flowevent_utils, ONLY: t_upatmoRestartAttributes, upatmoRestartAttributesAssign, upatmoRestartAttributesPack

#ifndef __NO_ICON_OCEAN__
    USE mo_ocean_nml, ONLY: lhamocc
    USE mo_hamocc_nml, ONLY: ks, dzsed
    USE mo_math_utilities, ONLY: set_zlev
#endif

#ifndef __NO_JSBACH__
    USE mo_echam_phy_config,  ONLY: echam_phy_config
    USE mo_jsb_vertical_axes, ONLY: set_vertical_grids_jsbach
#endif

    IMPLICIT NONE

    PRIVATE


    TYPE optional_integer
      INTEGER :: v !< value if present
      LOGICAL :: present = .FALSE.
    END TYPE optional_integer
    INTERFACE assign_if_present
      MODULE PROCEDURE assign_if_present_optional_i
    END INTERFACE assign_if_present

    PUBLIC :: t_restart_patch_description

    ! TYPE t_restart_patch_description contains all the DATA that
    ! describes a patch for restart purposes
    !
    ! @todo There does not seem to be a destructor routine for this
    !       object?!
    !
    TYPE t_restart_patch_description
        ! vertical grid definitions
        TYPE(t_Vgrid), ALLOCATABLE :: v_grid_defs(:)
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

        ! id of PE0 of working group (/= 0 in case of processor splitting)
        INTEGER :: work_pe0_id

        ! base file name contains already logical patch ident
        CHARACTER(LEN = filename_max) :: base_filename

        ! dynamic patch arguments (mandatory)
        INTEGER :: nold,nnow,nnew,nnew_rcf,nnow_rcf

        INTEGER :: opt_ndom = 1
        !> only considered if >0
        INTEGER :: opt_depth_lnd = -1
        ! dynamic patch arguments (optionally)
        TYPE(optional_integer) :: opt_nlev_snow, opt_nice_class, opt_ndyn_substeps, &
             opt_jstep_adv_marchuk_order, opt_ocean_zlevels

        REAL(wp), ALLOCATABLE :: opt_pvct(:)
        REAL(wp), ALLOCATABLE :: opt_t_elapsed_phy(:)
        REAL(wp), ALLOCATABLE :: opt_ocean_zheight_cellMiddle(:)
        REAL(wp), ALLOCATABLE :: opt_ocean_zheight_cellInterfaces(:)

        TYPE(t_upatmoRestartAttributes) :: opt_upatmo_restart_atts

        ! these are used for synchronous restart writing
        TYPE(t_comm_gather_pattern), POINTER :: cellGatherPattern, vertGatherPattern, edgeGatherPattern
    CONTAINS
        PROCEDURE :: init => restartPatchDescription_init

        ! called to set the DATA that may change from restart to restart
        PROCEDURE :: update => restartPatchDescription_update

        ! set the time level fields (nold, ...) to match the
        ! respective global variables
        PROCEDURE :: setTimeLevels => restartPatchDescription_setTimeLevels

        PROCEDURE :: packer => restartPatchDescription_packer
        PROCEDURE :: updateOnMaster => restartPatchDescription_updateOnMaster
        PROCEDURE :: updateVGrids => restartPatchDescription_updateVGrids
        PROCEDURE :: getGatherPattern => restartPatchDescription_getGatherPattern
        PROCEDURE :: getGlobalGridSize => restartPatchDescription_getGlobalGridSize
    END TYPE t_restart_patch_description

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_patch_description"

CONTAINS

  SUBROUTINE restartPatchDescription_init(me, jg)
    CLASS(t_restart_patch_description), INTENT(INOUT) :: me
    INTEGER, INTENT(in) :: jg

    ! DEFAULT initialization of all variables
    me%id = jg
    me%v_grid_count = 0
    me%l_dom_active = .FALSE.
    me%nold = nold(jg)
    me%nnow = nnow(jg)
    me%nnew = nnew(jg)
    me%nnew_rcf = nnew_rcf(jg)
    me%nnow_rcf = nnow_rcf(jg)
    ALLOCATE(me%v_grid_defs(zaxisTypeList%za_count()))

    ! patch dependent info, p_patch IS NOT available on restart PEs
    IF(my_process_is_work()) THEN
      me%work_pe0_id = p_patch(jg)%proc0
      me%nlev = p_patch(jg)%nlev
      me%cell_type = p_patch(jg)%geometry_info%cell_type
      me%base_filename = TRIM(p_patch(jg)%grid_filename)
      me%n_patch_cells_g = p_patch(jg)%n_patch_cells_g
      me%n_patch_verts_g = p_patch(jg)%n_patch_verts_g
      me%n_patch_edges_g = p_patch(jg)%n_patch_edges_g
      me%cellGatherPattern => p_patch(jg)%comm_pat_gather_c
      me%vertGatherPattern => p_patch(jg)%comm_pat_gather_v
      me%edgeGatherPattern => p_patch(jg)%comm_pat_gather_e
    ELSE
      me%work_pe0_id = -1
      me%nlev = -1
      me%cell_type = -1
      me%base_filename = ''
      me%n_patch_cells_g = -1
      me%n_patch_verts_g = -1
      me%n_patch_edges_g = -1
      NULLIFY(me%cellGatherPattern, me%vertGatherPattern, me%edgeGatherPattern)
    END IF
  END SUBROUTINE restartPatchDescription_init

    SUBROUTINE restartPatchDescription_update(me, patch, opt_pvct, opt_t_elapsed_phy, &
                                             &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth_lnd, &
                                             &opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels, &
                                             &opt_ocean_zheight_cellMiddle, opt_ocean_zheight_cellInterfaces, &
                                             &opt_upatmo_restart_atts)
        CLASS(t_restart_patch_description), INTENT(INOUT) :: me
        TYPE(t_patch), INTENT(IN) :: patch
        INTEGER, INTENT(IN), OPTIONAL :: opt_depth_lnd, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, &
                                       & opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels
        REAL(wp), INTENT(IN), OPTIONAL :: opt_pvct(:), opt_t_elapsed_phy(:), opt_ocean_zheight_cellMiddle(:), &
             & opt_ocean_zheight_cellInterfaces(:)
        TYPE(t_upatmoRestartAttributes), INTENT(IN), OPTIONAL :: opt_upatmo_restart_atts

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartPatchDescription_update"

        IF(me%id /= patch%id) CALL finish(routine, "assertion failed: wrong patch passed to update()")

        ! update activity flag - this needs to be done on all compute
        ! PEs, otherwise starting a nest during runtime would lead to
        ! incomplete restart files
        me%l_dom_active = patch%ldom_active

        ! otherwise, only the patch master process needs the dynamic
        ! restart arguments, these will be communicated to the other
        ! processes anyway
        IF(p_pe_work == 0 .OR. p_pe_work == me%work_pe0_id) THEN
            ! Patch-dependent attributes
            CALL assign_if_present_allocatable(me%opt_pvct, opt_pvct)
            CALL assign_if_present_allocatable(me%opt_t_elapsed_phy, opt_t_elapsed_phy)
            CALL assign_if_present_allocatable(me%opt_ocean_zheight_cellMiddle, opt_ocean_zheight_cellMiddle)
            CALL assign_if_present_allocatable(me%opt_ocean_zheight_cellInterfaces, opt_ocean_zheight_cellInterfaces)
            CALL assign_if_present(me%opt_ndyn_substeps, opt_ndyn_substeps)
            CALL assign_if_present(me%opt_jstep_adv_marchuk_order, opt_jstep_adv_marchuk_order)
            CALL assign_if_present(me%opt_depth_lnd, opt_depth_lnd)
            CALL assign_if_present(me%opt_nlev_snow, opt_nlev_snow)
            CALL assign_if_present(me%opt_nice_class, opt_nice_class)
            CALL assign_if_present(me%opt_ndom, opt_ndom)
            CALL assign_if_present(me%opt_ocean_zlevels, opt_ocean_zlevels)
            CALL upatmoRestartAttributesAssign(me%id, me%opt_upatmo_restart_atts, opt_upatmo_restart_atts)

            ! consistency check for OPTIONAL ocean variables
            IF(ALLOCATED(me%opt_ocean_zheight_cellMiddle)) THEN
                IF(.NOT. ALLOCATED(me%opt_ocean_Zheight_CellInterfaces) .OR. .NOT. me%opt_ocean_Zlevels%present) THEN
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

    SUBROUTINE restartPatchDescription_packer(me, operation, packedMessage)
        INTEGER, INTENT(in) :: operation
        CLASS(t_restart_patch_description), INTENT(INOUT) :: me
        CLASS(t_PackedMessage), INTENT(INOUT) :: packedMessage

        ! patch id AND activity flag
        CALL packedMessage%packer(operation, me%id)
        CALL packedMessage%packer(operation, me%l_dom_active)
        CALL packedMessage%packer(operation, me%nlev)
        CALL packedMessage%packer(operation, me%cell_type)
        CALL packedMessage%packer(operation, me%n_patch_cells_g)
        CALL packedMessage%packer(operation, me%n_patch_edges_g)
        CALL packedMessage%packer(operation, me%n_patch_verts_g)
        CALL packedMessage%packer(operation, me%work_pe0_id)
        CALL packedMessage%packer(operation, me%base_filename)

        ! time levels
        CALL packedMessage%packer(operation, me%nold)
        CALL packedMessage%packer(operation, me%nnow)
        CALL packedMessage%packer(operation, me%nnow_rcf)
        CALL packedMessage%packer(operation, me%nnew)
        CALL packedMessage%packer(operation, me%nnew_rcf)

        ! optional parameter values
        CALL packedMessage%packer(operation, me%opt_depth_lnd)
        CALL pack_optional_integer(packedMessage, operation, me%opt_nlev_snow)
        CALL pack_optional_integer(packedMessage, operation, me%opt_nice_class)
        CALL pack_optional_integer(packedMessage, operation, me%opt_ndyn_substeps)
        CALL pack_optional_integer(packedMessage, operation, me%opt_jstep_adv_marchuk_order)
        CALL packedMessage%packer(operation, me%opt_ndom)
        CALL pack_optional_integer(packedMessage, operation, me%opt_ocean_zlevels)

        ! optional parameter arrays
        CALL packedMessage%packer(operation, me%opt_pvct)
        CALL packedMessage%packer(operation, me%opt_t_elapsed_phy)

        CALL upatmoRestartAttributesPack(me%id, me%opt_upatmo_restart_atts, packedMessage, operation)
      END SUBROUTINE restartPatchDescription_packer

  SUBROUTINE assign_if_present_optional_i(o, opt_arg)
    TYPE(optional_integer), INTENT(out) :: o
    INTEGER, OPTIONAL, INTENT(in) :: opt_arg
    o%present = PRESENT(opt_arg)
    IF (o%present) THEN
      o%v = opt_arg
    ELSE
      o%v = -HUGE(o%v)
    END IF
  END SUBROUTINE assign_if_present_optional_i

  SUBROUTINE pack_optional_integer(msg, op, oi)
    CLASS(t_PackedMessage), INTENT(INOUT) :: msg
    INTEGER, INTENT(in) :: op
    TYPE(optional_integer) :: oi
    CALL msg%packer(op, oi%present)
    IF (oi%present) CALL msg%packer(op, oi%v)
  END SUBROUTINE pack_optional_integer

    ! This ensures that the work master has complete up-to-date
    ! knowledge of the patch description.  To this END, this calls
    ! setTimeLevels() on the patchDescription, AND THEN sends the
    ! description from the subset master to the work master.
    SUBROUTINE restartPatchDescription_updateOnMaster(me)
        CLASS(t_restart_patch_description), INTENT(INOUT) :: me
        TYPE(t_PackedMessage) :: packedMessage

        ! First ensure that the patch description IS up to date.
        CALL me%setTimeLevels() ! copy the global variables (nold, ...) to the patch description

        ! Then communicate it to the work master.
        IF(me%work_pe0_id == process_mpi_root_id) RETURN   ! nothing to communicate IF PE0 IS already the subset master

        IF(my_process_is_mpi_workroot()) THEN
            ! receive the package for this patch
            CALL packedMessage%recv(me%work_pe0_id, 0, process_mpi_all_comm)
            CALL me%packer(kUnpackOp, packedMessage)
        ELSE IF (p_pe_work == me%work_pe0_id) THEN
            ! send the time dependent DATA to process 0
            CALL me%packer(kPackOp, packedMessage)
            CALL packedMessage%send(process_mpi_root_id, 0, process_mpi_all_comm)
        END IF
    END SUBROUTINE restartPatchDescription_updateOnMaster

    !  Set vertical grid definition.
    SUBROUTINE restartPatchDescription_updateVGrids(me)
        CLASS(t_restart_patch_description), TARGET, INTENT(INOUT) :: me
        INTEGER :: nice_class
#ifndef __NO_ICON_OCEAN__
        INTEGER :: error
        REAL(wp), ALLOCATABLE :: levels(:), levels_sp(:)
#endif
        CHARACTER(*), PARAMETER :: routine = modname//":restartPatchDescription_updateVGrids"

        ! DEFAULT values for the level counts
        nice_class = 1

        ! Reset counter for vertical axes (needed for writing
        ! checkpoint files at multiple times)
        me%v_grid_count = 0

        ! replace DEFAULT values by the overrides provided IN the me
        IF (me%opt_nice_class%present) nice_class = me%opt_nice_class%v

        ! set vertical grid definitions
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_SURFACE, 0._wp)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_REFERENCE, me%nlev)
        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_REFERENCE_HALF, me%nlev+1)
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
        IF (me%opt_depth_lnd > 0) THEN
          CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, &
               ZA_DEPTH_BELOW_LAND, me%opt_depth_lnd)
          CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, &
               ZA_DEPTH_BELOW_LAND_P1, me%opt_depth_lnd+1)
        END IF
        IF(me%opt_nlev_snow%present) THEN
          CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, &
               ZA_SNOW, me%opt_nlev_snow%v)
          CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, &
               ZA_SNOW_HALF, me%opt_nlev_snow%v+1)
        END IF
        IF (me%opt_ocean_zlevels%present) THEN
          CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, &
               ZA_DEPTH_BELOW_SEA, me%opt_ocean_zlevels%v)
          CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, &
               ZA_DEPTH_BELOW_SEA_HALF, me%opt_ocean_zlevels%v+1)
        END IF

        CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_TROPOPAUSE, 1)

#ifndef __NO_JSBACH__
        IF (ANY(echam_phy_config(:)%ljsb)) CALL set_vertical_grids_jsbach(me%v_grid_defs, me%v_grid_count)
#endif

#ifndef __NO_ICON_OCEAN__
        IF(lhamocc) THEN
            ! HAMOCC sediment
            ALLOCATE(levels(ks), levels_sp(ks + 1), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

            CALL set_zlev(levels_sp, levels, ks, dzsed*1000._wp)
            CALL set_vertical_grid(me%v_grid_defs, me%v_grid_count, ZA_OCEAN_SEDIMENT, REAL(levels,wp))

            DEALLOCATE(levels, levels_sp)
        END IF
#endif

    END SUBROUTINE restartPatchDescription_updateVGrids

    FUNCTION restartPatchDescription_getGatherPattern(me, gridType) RESULT(resultVar)
        CLASS(t_restart_patch_description), INTENT(IN) :: me
        INTEGER, INTENT(IN) :: gridType
        TYPE(t_comm_gather_pattern), POINTER :: resultVar

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartPatchDescription_getGatherPattern"

        SELECT CASE(gridType)
            CASE(GRID_UNSTRUCTURED_CELL)
                resultVar => me%cellGatherPattern
            CASE(GRID_UNSTRUCTURED_VERT)
                resultVar => me%vertGatherPattern
            CASE(GRID_UNSTRUCTURED_EDGE)
                resultVar => me%edgeGatherPattern
            CASE DEFAULT
                CALL finish(routine, "illegal gridType argument")
        END SELECT
    END FUNCTION restartPatchDescription_getGatherPattern

    INTEGER FUNCTION restartPatchDescription_getGlobalGridSize(me, gridType) RESULT(resultVar)
        CLASS(t_restart_patch_description), INTENT(IN) :: me
        INTEGER, INTENT(in) :: gridType

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartPatchDescription_getGlobalGridSize"

        SELECT CASE(gridType)
            CASE(GRID_UNSTRUCTURED_CELL)
                resultVar = me%n_patch_cells_g
            CASE(GRID_UNSTRUCTURED_VERT)
                resultVar = me%n_patch_verts_g
            CASE(GRID_UNSTRUCTURED_EDGE)
                resultVar = me%n_patch_edges_g
            CASE DEFAULT
                CALL finish(routine, "illegal gridType argument")
        END SELECT
    END FUNCTION restartPatchDescription_getGlobalGridSize

END MODULE mo_restart_patch_description
