!>
!! Contains common helper routines for(a)synchronous restart
!! ----------------------------------------------------------
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_util_restart

    USE mo_exception, ONLY: finish, message, message_text
    USE mo_cdi, ONLY: CDI_UNDEFID, GRID_UNSTRUCTURED, gridCreate, gridDefNvertex, gridDefXname, gridDefXlongname, gridDefXunits, &
                    & gridDefYname, gridDefYlongname, gridDefYunits, zaxisCreate, zaxisDefLevels, streamOpenWrite, &
                    & vlistCreate, taxisCreate, vlistDefTaxis, TAXIS_ABSOLUTE, vlistDefVar, vlistDefVarDatatype, vlistDefVarName, &
                    & vlistDefVarLongname, vlistDefVarUnits, vlistDefVarMissval, TIME_VARIABLE, DATATYPE_FLT64, taxisDefVdate, &
                    & taxisDefVtime, cdiEncodeDate, cdiEncodeTime, streamDefTimestep
    USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_COUNT, GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT, &
                              & ZA_COUNT, ZA_DEPTH_BELOW_LAND, ZA_DEPTH_BELOW_LAND_P1, ZA_DEPTH_BELOW_SEA, &
                              & ZA_DEPTH_BELOW_SEA_HALF, ZA_DEPTH_RUNOFF_G, ZA_DEPTH_RUNOFF_S, ZA_GENERIC_ICE, ZA_HEIGHT_10M, &
                              & ZA_HEIGHT_2M, ZA_HYBRID, ZA_HYBRID_HALF, ZA_LAKE_BOTTOM, ZA_LAKE_BOTTOM_HALF, ZA_MIX_LAYER, &
                              & ZA_SEDIMENT_BOTTOM_TW_HALF, ZA_SNOW, ZA_SNOW_HALF, ZA_SURFACE, ZA_TOA, cdi_zaxis_types

    USE mo_cf_convention, ONLY: cf_global_info
    USE mo_datetime, ONLY: t_datetime, iso8601, iso8601extended
    USE mo_dynamics_config, ONLY: nold, nnow, nnew, nnew_rcf, nnow_rcf
    USE mo_exception, ONLY: get_filename_noext
    USE mo_fortran_tools, ONLY: assign_if_present, assign_if_present_allocatable
    USE mo_impl_constants, ONLY: SUCCESS, MAX_CHAR_LENGTH
    USE mo_io_restart_attributes, ONLY: t_RestartAttributeList
    USE mo_io_restart_namelist, ONLY: RestartNamelist_writeToFile
    USE mo_io_units, ONLY: filename_max
    USE mo_kind, ONLY: wp
    USE mo_model_domain, ONLY: t_patch
    USE mo_mpi, ONLY: p_pe_work
    USE mo_packed_message, ONLY: t_PackedMessage
    USE mo_run_config, ONLY: restart_filename
    USE mo_util_cdi, ONLY: cdiGetStringError
    USE mo_util_file, ONLY: util_symlink, util_islink, util_unlink
    USE mo_util_string, ONLY: int2string, associate_keyword, with_keywords, t_keyword_list
    USE mo_var_metadata_types, ONLY: t_var_metadata

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: t_v_grid
    PUBLIC :: t_restart_cdi_ids
    PUBLIC :: t_var_data
    PUBLIC :: t_restart_patch_description

    PUBLIC :: getRestartFilename
    PUBLIC :: setGeneralRestartAttributes
    PUBLIC :: setDynamicPatchRestartAttributes
    PUBLIC :: setPhysicsRestartAttributes
    PUBLIC :: set_vertical_grid
    PUBLIC :: create_restart_file_link

    ! maximumm number of verticale axes
    INTEGER, PARAMETER :: MAX_VERTICAL_AXES = 19

    ! TYPE t_v_grid contains the data of a vertical grid definition.
    TYPE t_v_grid
        INTEGER :: type
        REAL(wp), ALLOCATABLE :: levels(:)
    END type t_v_grid

    ! TYPE t_restart_cdi_ids IS just a simple container for all the different CDI IDs connected to a single restart file.
    TYPE t_restart_cdi_ids
        INTEGER :: file, vlist, taxis, hgrids(GRID_UNSTRUCTURED_COUNT), vgrids(ZA_COUNT)
    CONTAINS
        PROCEDURE :: init => restartCdiIds_init
        PROCEDURE :: openRestartAndCreateIds => restartCdiIds_openRestartAndCreateIds
        PROCEDURE :: finalizeVlist => restartCdiIds_finalizeVlist
        PROCEDURE :: defineVariable => restartCdiIds_defineVariable
        PROCEDURE :: closeAndDestroyIds => restartCdiIds_closeAndDestroyIds
    END TYPE t_restart_cdi_ids

    ! This takes a buffer for grid definitions, to which one entry IS appended by incrementing the count of used elements that's passed as well.
    INTERFACE set_vertical_grid
        MODULE PROCEDURE set_vertical_grid_array
        MODULE PROCEDURE set_vertical_grid_counted
        MODULE PROCEDURE set_vertical_grid_single
    END INTERFACE

    ! TYPE t_var_data (restart variable)
    TYPE t_var_data
        REAL(wp), POINTER :: r_ptr(:,:,:,:,:)
        TYPE(t_var_metadata) :: info
    END TYPE t_var_data

    ! TYPE t_restart_patch_description contains all the DATA that describes a patch for restart purposes
    TYPE t_restart_patch_description
        ! vertical grid definitions
        TYPE(t_v_grid), POINTER :: v_grid_defs(:)
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

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_util_restart"

CONTAINS

    FUNCTION getRestartFilename(baseName, domain, datetime, modelTypeName) RESULT(RESULT)
        CHARACTER(LEN = *), INTENT(IN) :: baseName, modelTypeName
        INTEGER, VALUE :: domain
        TYPE(t_datetime), INTENT(IN) :: datetime
        CHARACTER(LEN = :), ALLOCATABLE :: RESULT

        CHARACTER(LEN=32) :: datetimeString
        TYPE(t_keyword_list), POINTER :: keywords => NULL()

        datetimeString = iso8601(datetime)

        ! build the keyword list
        CALL associate_keyword("<gridfile>", TRIM(get_filename_noext(baseName)), keywords)
        CALL associate_keyword("<idom>", TRIM(int2string(domain, "(i2.2)")), keywords)
        CALL associate_keyword("<rsttime>", TRIM(datetimeString), keywords)
        CALL associate_keyword("<mtype>", TRIM(modelTypeName), keywords)

        ! replace keywords in file name
        RESULT = TRIM(with_keywords(keywords, TRIM(restart_filename)))
    END FUNCTION getRestartFilename

    SUBROUTINE setGeneralRestartAttributes(restartAttributes, datetime, n_dom, jstep, opt_output_jfile)
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
        TYPE(t_datetime), INTENT(IN) :: datetime
        INTEGER, VALUE :: n_dom, jstep
        INTEGER, OPTIONAL, INTENT(IN) :: opt_output_jfile(:)

        INTEGER :: i

        ! set CF-Convention required restart attributes
        CALL restartAttributes%setText('title',       TRIM(cf_global_info%title))
        CALL restartAttributes%setText('institution', TRIM(cf_global_info%institution))
        CALL restartAttributes%setText('source',      TRIM(cf_global_info%source))
        CALL restartAttributes%setText('history',     TRIM(cf_global_info%history))
        CALL restartAttributes%setText('references',  TRIM(cf_global_info%references))
        CALL restartAttributes%setText('comment',     TRIM(cf_global_info%comment))

        CALL restartAttributes%setReal( 'current_caltime', datetime%caltime )
        CALL restartAttributes%setInteger( 'current_calday' , INT(datetime%calday) )   !FIXME: Either it IS a bug that calday IS a 64bit INTEGER, OR it IS a bug that ONLY 32 bit of it are stored IN the restart file. Either way this needs to be fixed.
        CALL restartAttributes%setReal( 'current_daysec' , datetime%daysec )
        CALL restartAttributes%setText('tc_startdate', iso8601extended(datetime))   ! in preparation for move to mtime

        ! no. of domains AND simulation step
        CALL restartAttributes%setInteger( 'n_dom', n_dom)
        CALL restartAttributes%setInteger( 'jstep', jstep )

        IF(PRESENT(opt_output_jfile)) THEN
            DO i = 1, SIZE(opt_output_jfile)
                CALL restartAttributes%setInteger('output_jfile_'//TRIM(int2string(i, '(i2.2)')), opt_output_jfile(i) )
            END DO
        END IF
    END SUBROUTINE setGeneralRestartAttributes

    SUBROUTINE setDynamicPatchRestartAttributes(restartAttributes, jg, nold, nnow, nnew, nnow_rcf, nnew_rcf)
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
        INTEGER, VALUE :: jg, nold, nnow, nnew, nnow_rcf, nnew_rcf

        CHARACTER(LEN = 2) :: jgString

        jgString = TRIM(int2string(jg, "(i2.2)"))

        CALL restartAttributes%setInteger('nold_DOM'//jgString, nold)
        CALL restartAttributes%setInteger('nnow_DOM'//jgString, nnow)
        CALL restartAttributes%setInteger('nnew_DOM'//jgString, nnew)
        CALL restartAttributes%setInteger('nnow_rcf_DOM'//jgString, nnow_rcf)
        CALL restartAttributes%setInteger('nnew_rcf_DOM'//jgString, nnew_rcf)
    END SUBROUTINE setDynamicPatchRestartAttributes

    SUBROUTINE setPhysicsRestartAttributes(restartAttributes, jg, t_elapsed_phy, lcall_phy)
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
        INTEGER, VALUE :: jg
        REAL(wp), INTENT(IN) :: t_elapsed_phy(:)
        LOGICAL, INTENT(IN) :: lcall_phy(:)

        INTEGER :: i
        CHARACTER(LEN = :), ALLOCATABLE :: prefix

        prefix = 't_elapsed_phy_DOM'//TRIM(int2string(jg, "(i2.2)"))//'_PHY'
        DO i = 1, SIZE(t_elapsed_phy)
            CALL restartAttributes%setReal(prefix//TRIM(int2string(i, '(i2.2)')), t_elapsed_phy(i) )
        END DO

        prefix = 'lcall_phy_DOM'//TRIM(int2string(jg, "(i2.2)"))//'_PHY'
        DO i = 1, SIZE(lcall_phy)
            CALL restartAttributes%setLogical(prefix//TRIM(int2string(i, '(i2.2)')), lcall_phy(i) )
        END DO
    END SUBROUTINE setPhysicsRestartAttributes

    ! Creates a horizontal grid definition from the given parameters, returns the new CDI gridId.
    INTEGER FUNCTION create_cdi_hgrid_def(iCnt, iNVert, cNameX, cLNameX, cUnitsX, &
        &                                               cNameY, cLNameY, cUnitsY) RESULT(RESULT)

        INTEGER, VALUE :: iCnt, iNVert
        CHARACTER(LEN = *), INTENT(IN) :: cNameX, cLNameX, cUnitsX
        CHARACTER(LEN = *), INTENT(IN) :: cNameY, cLNameY, cUnitsY
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":create_cdi_hgrid_def"

        RESULT = gridCreate(GRID_UNSTRUCTURED, iCnt)
        IF(RESULT == CDI_UNDEFID) CALL finish(routine, "error creating CDI grid")
        CALL gridDefNvertex(RESULT, iNVert)

        CALL gridDefXname(RESULT, TRIM(cNameX))
        CALL gridDefXlongname(RESULT, TRIM(cLNameX))
        CALL gridDefXunits(RESULT, TRIM(cUnitsX))

        CALL gridDefYname(RESULT, TRIM(cNameY))
        CALL gridDefYlongname(RESULT, TRIM(cLNameY))
        CALL gridDefYunits(RESULT, TRIM(cUnitsY))
    END FUNCTION create_cdi_hgrid_def

    FUNCTION createHgrids(cellCount, vertexCount, edgeCount, cellType) RESULT(RESULT)
        INTEGER, VALUE :: cellCount, vertexCount, edgeCount, cellType
        INTEGER :: RESULT(GRID_UNSTRUCTURED_COUNT)

        RESULT(GRID_UNSTRUCTURED_CELL) = create_cdi_hgrid_def(cellCount, cellType, &
                                                             &'clon', 'center longitude', 'radian', &
                                                             &'clat', 'center latitude', 'radian')

        RESULT(GRID_UNSTRUCTURED_VERT) = create_cdi_hgrid_def(vertexCount, 9 - cellType, &
                                                             &'vlon', 'vertex longitude', 'radian', &
                                                             &'vlat', 'vertex latitude', 'radian')

        RESULT(GRID_UNSTRUCTURED_EDGE) = create_cdi_hgrid_def(edgeCount, 4, &
                                                             &'elon', 'edge midpoint longitude', 'radian', &
                                                             &'elat', 'edge midpoint latitude', 'radian')

    END FUNCTION createHgrids

    ! If no opt_levelValues are given, this defaults to numbering the levels from 1 to levelCount.
    INTEGER FUNCTION defineVAxis(cdiAxisType, levelValues) RESULT(RESULT)
        INTEGER, VALUE :: cdiAxisType
        REAL(KIND = wp), INTENT(IN) :: levelValues(:)

        INTEGER :: levelCount

        levelCount = SIZE(levelValues, 1)
        RESULT = zaxisCreate(cdiAxisType, levelCount)
        CALL zaxisDefLevels(RESULT, levelValues)
    END FUNCTION defineVAxis

    SUBROUTINE createVgrids(axisIds, gridDescriptions, opt_vct)
        INTEGER, INTENT(INOUT) :: axisIds(ZA_COUNT)
        TYPE(t_v_grid), INTENT(IN) :: gridDescriptions(:)
        REAL(wp), OPTIONAL, INTENT(IN) :: opt_vct(:)

        INTEGER :: i, nlevp1, gridId
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":createVgrids"

        axisIds(:) = CDI_UNDEFID
        DO i = 1, SIZE(gridDescriptions, 1)
            IF(cdi_zaxis_types(gridDescriptions(i)%TYPE) == CDI_UNDEFID) CALL finish(routine,"no CDI zaxis TYPE defined for vgrid")

            gridId = defineVAxis(cdi_zaxis_types(gridDescriptions(i)%TYPE), gridDescriptions(i)%levels)
            IF(gridId == CDI_UNDEFID) CALL finish(routine, "defineVAxis() returned an error")
            axisIds(gridDescriptions(i)%TYPE) = gridId

            SELECT CASE (gridDescriptions(i)%type)
                CASE (ZA_HYBRID)
                    IF(.NOT.PRESENT(opt_vct)) CYCLE
                    nlevp1 = SIZE(gridDescriptions(i)%levels, 1) + 1
                    CALL zaxisDefVct(gridId, 2*nlevp1, opt_vct(1:2*nlevp1))

                CASE (ZA_HYBRID_HALF)
                    IF(.NOT.PRESENT(opt_vct)) CYCLE
                    nlevp1 = SIZE(gridDescriptions(i)%levels, 1)
                    CALL zaxisDefVct(gridId, 2*nlevp1, opt_vct(1:2*nlevp1))

                CASE (ZA_LAKE_BOTTOM, ZA_MIX_LAYER)
                    CALL zaxisDefLbounds(gridId, [1._wp]) !necessary for GRIB2
                    CALL zaxisDefUbounds(gridId, [0._wp]) !necessary for GRIB2
                    CALL zaxisDefUnits  (gridId, "m")

                CASE (ZA_LAKE_BOTTOM_HALF, ZA_SEDIMENT_BOTTOM_TW_HALF)
                    CALL zaxisDefUnits(gridId, "m")
            END SELECT
        ENDDO
    END SUBROUTINE createVgrids

    SUBROUTINE restartCdiIds_init(me)
        CLASS(t_restart_cdi_ids), INTENT(INOUT) :: me

        me%file = CDI_UNDEFID
        me%vlist = CDI_UNDEFID
        me%taxis = CDI_UNDEFID
        me%hgrids(:) = CDI_UNDEFID
        me%vgrids(:) = CDI_UNDEFID
    END SUBROUTINE restartCdiIds_init

    SUBROUTINE restartCdiIds_openRestartAndCreateIds(me, filename, restartType, restartAttributes, cellCount, vertCount, &
                                                    &edgeCount, cellType, vgridDefs, opt_vct)
        CLASS(t_restart_cdi_ids), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: filename
        INTEGER, VALUE :: restartType, cellCount, vertCount, edgeCount, cellType
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
        TYPE(t_v_grid), INTENT(IN) :: vgridDefs(:)
        REAL(wp), INTENT(IN), OPTIONAL :: opt_vct(:)

        CHARACTER(LEN = MAX_CHAR_LENGTH) :: cdiErrorText
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartCdiIds_openRestartAndCreateIds"

        ! open the file
        me%file = streamOpenWrite(filename, restartType)

        IF(me%file < 0) THEN
            CALL cdiGetStringError(me%file, cdiErrorText)
            WRITE(message_text,'(a)') TRIM(cdiErrorText)
            CALL message('',message_text)
            CALL finish(routine, 'open failed on '//filename)
        END IF

        ! create the CDI IDs we need

        ! 1. vlist
        me%vlist = vlistCreate()

        ! 2. global attributes
        CALL RestartNamelist_writeToFile(me%vlist)
        CALL restartAttributes%writeToFile(me%vlist)

        ! 3. horizontal grids
        me%hgrids = createHgrids(cellCount, vertCount, edgeCount, cellType)

        ! 4. vertical grids
        CALL createVgrids(me%vgrids, vgridDefs, opt_vct)

        ! 5. time axis (always absolute time for restart files)
        me%taxis = taxisCreate(TAXIS_ABSOLUTE)
        CALL vlistDefTaxis(me%vlist, me%taxis)
    END SUBROUTINE restartCdiIds_openRestartAndCreateIds

    ! CDI reqires the vlist of a stream to be set before the timestep can be defined, which IS why we combine these two operations into one SUBROUTINE.
    SUBROUTINE restartCdiIds_finalizeVlist(me, datetime)
        CLASS(t_restart_cdi_ids), INTENT(INOUT) :: me
        TYPE(t_datetime), INTENT(IN) :: datetime

        INTEGER :: trash

        ! define the vlist, so that we are allowed to define the timestep
        CALL streamDefVlist(me%file, me%vlist)

        ! define the timestep
        CALL taxisDefVdate(me%taxis, cdiEncodeDate(datetime%year, datetime%month, datetime%day))
        CALL taxisDefVtime(me%taxis, cdiEncodeTime(datetime%hour, datetime%minute, NINT(datetime%second)))
        trash = streamDefTimestep(me%file, 0)
    END SUBROUTINE restartCdiIds_finalizeVlist

    ! Encapsulates the CDI calls to define a variable.
    ! lIsInteger AND lIsLogical reflect the TYPE of the variable, IF neither IS set, the variable IS assumed to be of TYPE REAL.
    SUBROUTINE restartCdiIds_defineVariable(me, info, lIsInteger, lIsLogical)
        CLASS(t_restart_cdi_ids), INTENT(IN) :: me
        TYPE(t_var_metadata), INTENT(INOUT) :: info
        LOGICAL, VALUE :: lIsInteger, lIslogical

        INTEGER :: varId, gridId, zaxisId
        REAL(wp) :: casted_missval
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartCdiIds_defineVariable"

        IF(lIsInteger.AND.lIsLogical) THEN
            CALL finish(routine, "assertion failed: attempt to define a variable both as INTEGER and LOGICAL")
        END IF

        ! get the horizontal grid ID
        gridId = info%cdiGridID
        SELECT CASE (info%hgrid)
            CASE(GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE)
                gridId = me%hgrids(info%hgrid)
        END SELECT
        IF(gridId == CDI_UNDEFID) CALL finish(routine, 'Grid type not defined for field '//TRIM(info%name))

        ! get the vertical axis ID
        zaxisId = info%cdiZaxisID
        if(zaxisId < 0) zaxisId = me%vgrids(info%vgrid)
        IF(zaxisId == CDI_UNDEFID) CALL finish(routine, 'Z axis not defined for field '//TRIM(info%name))

        ! define the variable with the required info
        varId = vlistDefVar(me%vlist, gridId, zaxisId, TIME_VARIABLE)
        IF(varID == CDI_UNDEFID) CALL finish(routine, 'error WHILE defining CDI variable "'//TRIM(info%name)//'"')
        info%cdiVarID = varId
        CALL vlistDefVarDatatype(me%vlist, varId, DATATYPE_FLT64)
        CALL vlistDefVarName(me%vlist, varId, TRIM(info%name))

        ! then add the three optional fields
        IF(info%cf%long_name /= '') CALL vlistDefVarLongname(me%vlist, varId, TRIM(info%cf%long_name))
        IF(info%cf%units /= '') CALL vlistDefVarUnits(me%vlist, varId, TRIM(info%cf%units))
        IF(info%lmiss) THEN
            casted_missval = info%missval%rval
            IF(lIsInteger) casted_missval = REAL(info%missval%ival, wp)
            IF(lIsLogical) THEN
                casted_missval = 0.0_wp
                IF(info%missval%lval) casted_missval = 1.0_wp
            ENDIF
            CALL vlistDefVarMissval(me%vlist, varId, casted_missval)
        ENDIF
    END SUBROUTINE restartCdiIds_defineVariable

    SUBROUTINE restartCdiIds_closeAndDestroyIds(me)
        CLASS(t_restart_cdi_ids), INTENT(INOUT) :: me

        INTEGER :: i

        ! close/destroy all open CDI IDs
        IF(me%file /= CDI_UNDEFID) CALL streamClose(me%file)
        IF(me%vlist /= CDI_UNDEFID) CALL vlistDestroy(me%vlist)
        IF(me%taxis /= CDI_UNDEFID) CALL taxisDestroy(me%taxis)
        DO i = 1, SIZE(me%hgrids, 1)
            IF(me%hgrids(i) /= CDI_UNDEFID) CALL gridDestroy(me%hgrids(i))
        END DO
        DO i = 1, SIZE(me%vgrids, 1)
            IF(me%vgrids(i) /= CDI_UNDEFID) CALL zaxisDestroy(me%vgrids(i))
        END DO

        ! reset the IDs
        CALL me%init()
    END SUBROUTINE restartCdiIds_closeAndDestroyIds

    SUBROUTINE set_vertical_grid_array(gridDefinitions, gridCount, type, levels)
        TYPE(t_v_grid), INTENT(INOUT) :: gridDefinitions(:)
        INTEGER, INTENT(INOUT) :: gridCount
        INTEGER, VALUE :: type
        REAL(wp), INTENT(in) :: levels(:)

        INTEGER :: levelCount, error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":set_vertical_grid_array"

        gridCount = gridCount+1
        IF(gridCount > SIZE(gridDefinitions, 1)) CALL finish(routine, "insufficient space in the gridDefinitions array, please &
                                                                     &increase the size of the array that is passed")

        gridDefinitions(gridCount)%type = type

        levelCount = SIZE(levels, 1)
        IF(ALLOCATED(gridDefinitions(gridCount)%levels)) DEALLOCATE(gridDefinitions(gridCount)%levels)
        ALLOCATE(gridDefinitions(gridCount)%levels(levelCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        gridDefinitions(gridCount)%levels(:) = levels(:)
    END SUBROUTINE set_vertical_grid_array

    SUBROUTINE set_vertical_grid_counted(gridDefinitions, gridCount, TYPE, levelCount)
        TYPE(t_v_grid), INTENT(INOUT) :: gridDefinitions(:)
        INTEGER, INTENT(INOUT) :: gridCount
        INTEGER, VALUE :: TYPE, levelCount

        INTEGER :: i

        CALL set_vertical_grid(gridDefinitions, gridCount, TYPE, [(REAL(i, wp), i = 1, levelCount)])
    END SUBROUTINE set_vertical_grid_counted

    SUBROUTINE set_vertical_grid_single(gridDefinitions, gridCount, TYPE, levelValue)
        TYPE(t_v_grid), INTENT(INOUT) :: gridDefinitions(:)
        INTEGER, INTENT(INOUT) :: gridCount
        INTEGER, VALUE :: TYPE
        REAL(wp), VALUE :: levelValue

        CALL set_vertical_grid(gridDefinitions, gridCount, TYPE, [levelValue])
    END SUBROUTINE set_vertical_grid_single

    SUBROUTINE create_restart_file_link(filename, modelType, proc_id, jg, opt_ndom)
        CHARACTER(LEN = *), INTENT(IN) :: filename, modelType
        INTEGER, VALUE :: proc_id, jg
        INTEGER, INTENT(IN), OPTIONAL :: opt_ndom

        INTEGER :: iret, ndom
        CHARACTER(LEN = 12) :: procIdString
        CHARACTER(LEN = 64) :: linkname
        CHARACTER(LEN=*), PARAMETER :: routine = modname//':create_restart_file_link'

        ! we need to add a process dependent part to the link NAME IF there are several restart processes
        procIdString = ''
        IF(proc_id /= 0) procIdString = TRIM(int2string(proc_id))

        ! IN CASE we have ONLY a single domain / no domain information, USE "_DOM01" IN the link NAME
        ndom = 1
        CALL assign_if_present(ndom, opt_ndom)
        IF(ndom == 1) jg = 1

        ! build link name
        linkname = 'restart'//TRIM(procIdString)//'_'//modelType//"_DOM"//TRIM(int2string(jg, "(i2.2)"))//'.nc'

        ! delete old symbolic link, if exists
        ! FIXME[NH]: handle the CASE that we have a file at that location which IS NOT a symlink
        IF(util_islink(TRIM(linkname))) THEN
            iret = util_unlink(TRIM(linkname))
            IF(iret /= SUCCESS) WRITE(0, *) routine//': cannot unlink "'//TRIM(linkname)//'"'
        ENDIF

        ! create a new symbolic link
        iret = util_symlink(filename,TRIM(linkname))
        IF(iret /= SUCCESS) WRITE(0, *) routine//': cannot create symbolic link "'//TRIM(linkname)//'" for "'//filename//'"'
    END SUBROUTINE create_restart_file_link

    SUBROUTINE restartPatchDescription_init(description, p_patch)
        CLASS(t_restart_patch_description), INTENT(INOUT) :: description
        TYPE(t_patch), INTENT(IN) :: p_patch

        description%id = p_patch%id
        description%work_pe0_id = p_patch%proc0
        description%nlev = p_patch%nlev
        description%cell_type = p_patch%geometry_info%cell_type
        description%base_filename = TRIM(p_patch%grid_filename)
        description%n_patch_cells_g = p_patch%n_patch_cells_g
        description%n_patch_verts_g = p_patch%n_patch_verts_g
        description%n_patch_edges_g = p_patch%n_patch_edges_g
        description%v_grid_defs => NULL()
        description%v_grid_count = 0
        description%l_dom_active = .FALSE.
        description%restart_proc_id = 0
        description%nold = nold(p_patch%id)
        description%nnow = nnow(p_patch%id)
        description%nnew = nnew(p_patch%id)
        description%nnew_rcf = nnew_rcf(p_patch%id)
        description%nnow_rcf = nnow_rcf(p_patch%id)
        description%l_opt_depth = .FALSE.
        description%opt_depth = 0
        description%l_opt_depth_lnd = .FALSE.
        description%opt_depth_lnd = 0
        description%l_opt_nlev_snow = .FALSE.
        description%opt_nlev_snow = 0
        description%l_opt_nice_class = .FALSE.
        description%opt_nice_class = 0
        description%l_opt_ndyn_substeps = .FALSE.
        description%opt_ndyn_substeps = 0
        description%l_opt_jstep_adv_marchuk_order = .FALSE.
        description%opt_jstep_adv_marchuk_order = 0
        description%l_opt_ndom = .FALSE.
        description%opt_ndom = 0
        description%l_opt_sim_time = .FALSE.
        description%opt_sim_time = 0.
    END SUBROUTINE restartPatchDescription_init

    SUBROUTINE restartPatchDescription_update(description, patch, opt_pvct, opt_t_elapsed_phy, opt_lcall_phy, opt_sim_time, &
                                             &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth, opt_depth_lnd, &
                                             &opt_nlev_snow, opt_nice_class, opt_ndom)
        CLASS(t_restart_patch_description), INTENT(INOUT) :: description
        TYPE(t_patch), INTENT(IN) :: patch
        INTEGER, INTENT(IN), OPTIONAL :: opt_depth, opt_depth_lnd, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, &
                                       & opt_nlev_snow, opt_nice_class, opt_ndom
        REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time, opt_pvct(:), opt_t_elapsed_phy(:)
        LOGICAL, INTENT(IN), OPTIONAL :: opt_lcall_phy(:)

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartPatchDescription_update"

        IF(description%id /= patch%id) CALL finish(routine, "assertion failed: wrong patch passed to update()")

        ! update activity flag - this needs to be done on all compute PEs, otherwise starting a nest during runtime would lead to incomplete restart files
        description%l_dom_active = patch%ldom_active

        ! otherwise, only the patch master process needs the dynamic restart arguments, these will be communicated to the other processes anyway
        IF(p_pe_work == 0 .OR. p_pe_work == description%work_pe0_id) THEN
            ! Patch-dependent attributes
            CALL assign_if_present_allocatable(description%opt_pvct, opt_pvct)
            CALL assign_if_present_allocatable(description%opt_t_elapsed_phy, opt_t_elapsed_phy)
            CALL assign_if_present_allocatable(description%opt_lcall_phy, opt_lcall_phy)

            CALL assign_if_present(description%opt_ndyn_substeps, opt_ndyn_substeps)
            description%l_opt_ndyn_substeps = PRESENT(opt_ndyn_substeps)

            CALL assign_if_present(description%opt_jstep_adv_marchuk_order, opt_jstep_adv_marchuk_order)
            description%l_opt_jstep_adv_marchuk_order = PRESENT(opt_jstep_adv_marchuk_order)

            CALL assign_if_present(description%opt_depth, opt_depth)
            description%l_opt_depth = PRESENT(opt_depth)

            CALL assign_if_present(description%opt_depth_lnd, opt_depth_lnd)
            description%l_opt_depth_lnd = PRESENT(opt_depth_lnd)

            CALL assign_if_present(description%opt_nlev_snow, opt_nlev_snow)
            description%l_opt_nlev_snow = PRESENT(opt_nlev_snow)

            CALL assign_if_present(description%opt_nice_class, opt_nice_class)
            description%l_opt_nice_class = PRESENT(opt_nice_class)

            CALL assign_if_present(description%opt_sim_time, opt_sim_time)
            description%l_opt_sim_time = PRESENT(opt_sim_time)

            CALL assign_if_present(description%opt_ndom, opt_ndom)
            description%l_opt_ndom = PRESENT(opt_ndom)
        END IF
    END SUBROUTINE restartPatchDescription_update

    SUBROUTINE restartPatchDescription_setTimeLevels(description)
        CLASS(t_restart_patch_description), INTENT(INOUT) :: description

        description%nold = nold(description%id)
        description%nnow = nnow(description%id)
        description%nnow_rcf = nnow_rcf(description%id)
        description%nnew = nnew(description%id)
        description%nnew_rcf = nnew_rcf(description%id)
    END SUBROUTINE restartPatchDescription_setTimeLevels

    SUBROUTINE restartPatchDescription_packer(description, operation, message)
        INTEGER, VALUE :: operation
        CLASS(t_restart_patch_description), INTENT(INOUT) :: description
        TYPE(t_PackedMessage), INTENT(INOUT) :: message

        ! patch id AND activity flag
        CALL message%execute(operation, description%id)
        CALL message%execute(operation, description%l_dom_active)
        CALL message%execute(operation, description%nlev)
        CALL message%execute(operation, description%cell_type)
        CALL message%execute(operation, description%n_patch_cells_g)
        CALL message%execute(operation, description%n_patch_edges_g)
        CALL message%execute(operation, description%n_patch_verts_g)
        CALL message%execute(operation, description%restart_proc_id)
        CALL message%execute(operation, description%work_pe0_id)
        CALL message%execute(operation, description%base_filename)

        ! time levels
        CALL message%execute(operation, description%nold)
        CALL message%execute(operation, description%nnow)
        CALL message%execute(operation, description%nnow_rcf)
        CALL message%execute(operation, description%nnew)
        CALL message%execute(operation, description%nnew_rcf)

        ! optional parameter values
        CALL message%execute(operation, description%l_opt_depth)
        CALL message%execute(operation, description%opt_depth)
        CALL message%execute(operation, description%l_opt_depth_lnd)
        CALL message%execute(operation, description%opt_depth_lnd)
        CALL message%execute(operation, description%l_opt_nlev_snow)
        CALL message%execute(operation, description%opt_nlev_snow)
        CALL message%execute(operation, description%l_opt_nice_class)
        CALL message%execute(operation, description%opt_nice_class)
        CALL message%execute(operation, description%l_opt_ndyn_substeps)
        CALL message%execute(operation, description%opt_ndyn_substeps)
        CALL message%execute(operation, description%l_opt_jstep_adv_marchuk_order)
        CALL message%execute(operation, description%opt_jstep_adv_marchuk_order)
        CALL message%execute(operation, description%l_opt_sim_time)
        CALL message%execute(operation, description%opt_sim_time)
        CALL message%execute(operation, description%l_opt_ndom)
        CALL message%execute(operation, description%opt_ndom)

        ! optional parameter arrays
        CALL message%execute(operation, description%opt_pvct)
        CALL message%execute(operation, description%opt_lcall_phy)
        CALL message%execute(operation, description%opt_t_elapsed_phy)
    END SUBROUTINE restartPatchDescription_packer

    !  Set vertical grid definition.
    SUBROUTINE restartPatchDescription_defineVGrids(p_desc)
        CLASS(t_restart_patch_description), TARGET, INTENT(INOUT) :: p_desc

        INTEGER :: nlev_soil, nlev_snow, nlev_ocean, nice_class, ierrstat
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartPatchDescription_defineVGrids"

        ! DEFAULT values for the level counts
        nlev_soil = 0
        nlev_snow = 0
        nlev_ocean = 0
        nice_class = 1

        ! replace DEFAULT values by the overrides provided IN the p_desc
        IF(p_desc%l_opt_depth_lnd) nlev_soil = p_desc%opt_depth_lnd
        IF(p_desc%l_opt_nlev_snow) nlev_snow = p_desc%opt_nlev_snow
        IF(p_desc%l_opt_depth) nlev_ocean = p_desc%opt_depth
        IF(p_desc%l_opt_nice_class) nice_class = p_desc%opt_nice_class

        ! set vertical grid definitions
        ALLOCATE(p_desc%v_grid_defs(MAX_VERTICAL_AXES), STAT=ierrstat)
        p_desc%v_grid_count = 0
        IF(ierrstat /= SUCCESS) CALL finish(routine, "memory allocation failed")

        CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_SURFACE, 0._wp)
        CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_HYBRID, p_desc%nlev)
        CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_HYBRID_HALF, p_desc%nlev+1)
        CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_HEIGHT_2M, 2._wp)
        CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_HEIGHT_10M, 10._wp)
        CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_TOA, 1._wp)
        CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_LAKE_BOTTOM, 1._wp)
        CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_MIX_LAYER, 1._wp)
        CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_LAKE_BOTTOM_HALF, 1._wp)
        CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_SEDIMENT_BOTTOM_TW_HALF, 0._wp)
        CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_GENERIC_ICE, 1._wp)
        CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_DEPTH_RUNOFF_S, 1)
        CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_DEPTH_RUNOFF_G, 1)
        IF(p_desc%l_opt_depth_lnd) CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_DEPTH_BELOW_LAND, nlev_soil)
        IF(p_desc%l_opt_depth_lnd) CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_DEPTH_BELOW_LAND_P1, &
                                                         &nlev_soil+1)
        IF(p_desc%l_opt_nlev_snow) CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_SNOW, nlev_snow)
        IF(p_desc%l_opt_nlev_snow) CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_SNOW_HALF, nlev_snow+1)
        IF(p_desc%l_opt_depth) CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_DEPTH_BELOW_SEA, nlev_ocean)
        IF(p_desc%l_opt_depth) CALL set_vertical_grid(p_desc%v_grid_defs, p_desc%v_grid_count, ZA_DEPTH_BELOW_SEA_HALF, &
                                                     &nlev_ocean+1)
    END SUBROUTINE restartPatchDescription_defineVGrids

END MODULE mo_util_restart
