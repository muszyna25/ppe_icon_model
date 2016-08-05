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
    USE mo_cdi, ONLY: CDI_UNDEFID, GRID_UNSTRUCTURED, gridCreate, gridDefNvertex, gridDefXname, gridDefXlongname, gridDefXunits, &
                    & gridDefYname, gridDefYlongname, gridDefYunits, zaxisCreate, zaxisDefLevels, streamOpenWrite, &
                    & vlistCreate, taxisCreate, vlistDefTaxis, TAXIS_ABSOLUTE, vlistDefVar, vlistDefVarDatatype, vlistDefVarName, &
                    & vlistDefVarLongname, vlistDefVarUnits, vlistDefVarMissval, TIME_VARIABLE, DATATYPE_FLT64, taxisDefVdate, &
                    & taxisDefVtime, cdiEncodeDate, cdiEncodeTime, streamDefTimestep
    USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_COUNT, GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT, &
                              & ZA_COUNT, ZA_HYBRID, ZA_HYBRID_HALF, ZA_LAKE_BOTTOM, ZA_LAKE_BOTTOM_HALF, ZA_MIX_LAYER, &
                              & ZA_SEDIMENT_BOTTOM_TW_HALF, cdi_zaxis_types
    USE mo_cf_convention, ONLY: cf_global_info
    USE mo_datetime, ONLY: t_datetime, iso8601, iso8601extended
    USE mo_exception, ONLY: get_filename_noext, finish, message, message_text
    USE mo_fortran_tools, ONLY: assign_if_present, assign_if_present_allocatable, t_ptr_2d
    USE mo_impl_constants, ONLY: SUCCESS, MAX_CHAR_LENGTH
    USE mo_io_restart_attributes, ONLY: t_RestartAttributeList
    USE mo_io_restart_namelist, ONLY: RestartNamelist_writeToFile
    USE mo_kind, ONLY: wp, i8
    USE mo_packed_message, ONLY: t_PackedMessage, kPackOp, kUnpackOp
    USE mo_parallel_config, ONLY: nproma
    USE mo_run_config, ONLY: restart_filename
    USE mo_util_cdi, ONLY: cdiGetStringError
    USE mo_util_file, ONLY: util_symlink, util_islink, util_unlink
    USE mo_util_string, ONLY: int2string, real2string, associate_keyword, with_keywords, t_keyword_list
    USE mo_var_metadata_types, ONLY: t_var_metadata

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: t_v_grid
    PUBLIC :: t_restart_cdi_ids
    PUBLIC :: t_var_data
    PUBLIC :: t_restart_args

    PUBLIC :: getRestartFilename
    PUBLIC :: setGeneralRestartAttributes
    PUBLIC :: setDynamicPatchRestartAttributes
    PUBLIC :: setPhysicsRestartAttributes
    PUBLIC :: set_vertical_grid
    PUBLIC :: create_restart_file_link
    PUBLIC :: getLevelPointers

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

    ! patch independent restart arguments
    TYPE t_restart_args
        TYPE(t_datetime) :: datetime
        INTEGER :: jstep
        INTEGER, ALLOCATABLE :: output_jfile(:)
    CONTAINS
        PROCEDURE :: construct => restartArgs_construct
        PROCEDURE :: packer => restartArgs_packer   ! unpacking IS considered construction
        PROCEDURE :: print => restartArgs_print
        PROCEDURE :: destruct => restartArgs_destruct
    END TYPE t_restart_args

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

    SUBROUTINE restartArgs_construct(me, datetime, jstep, opt_output_jfile)
        CLASS(t_restart_args), INTENT(INOUT) :: me
        TYPE(t_datetime), INTENT(IN) :: datetime
        INTEGER, VALUE :: jstep
        INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

        me%datetime = datetime
        me%jstep = jstep
        CALL assign_if_present_allocatable(me%output_jfile, opt_output_jfile)
    END SUBROUTINE restartArgs_construct

    SUBROUTINE restartArgs_packer(me, operation, message)
        CLASS(t_restart_args), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        TYPE(t_PackedMessage), INTENT(INOUT) :: message

        INTEGER :: calday

        CALL message%packer(operation, me%datetime%year)
        CALL message%packer(operation, me%datetime%month)
        CALL message%packer(operation, me%datetime%day)
        CALL message%packer(operation, me%datetime%hour)
        CALL message%packer(operation, me%datetime%minute)
        CALL message%packer(operation, me%datetime%second)
        CALL message%packer(operation, me%datetime%caltime)
        IF(operation == kPackOp) calday = INT(me%datetime%calday)   !The IF IS needed to avoid overflow when me%datetime%calday IS uninitialized.
        CALL message%packer(operation, calday)
        IF(operation == kUnpackOp) me%datetime%calday = INT(calday,i8)
        CALL message%packer(operation, me%datetime%daysec)
        CALL message%packer(operation, me%jstep)
        CALL message%packer(operation, me%output_jfile)
    END SUBROUTINE restartArgs_packer

    SUBROUTINE restartArgs_print(me, prefix)
        CLASS(t_restart_args), INTENT(IN) :: me
        CHARACTER(LEN = *), INTENT(IN) :: prefix

        PRINT*, prefix//'current_calday='//TRIM(int2string(INT(me%datetime%calday)))
        PRINT*, prefix//'current_caltime='//TRIM(real2string(me%datetime%caltime))
        PRINT*, prefix//'current_daysec='//TRIM(real2string(me%datetime%daysec))
    END SUBROUTINE restartArgs_print

    SUBROUTINE restartArgs_destruct(me)
        CLASS(t_restart_args), INTENT(INOUT) :: me

        IF(ALLOCATED(me%output_jfile)) DEALLOCATE(me%output_jfile)
    END SUBROUTINE restartArgs_destruct

    SUBROUTINE getLevelPointers(varMetadata, varData, levelPointers)
        TYPE(t_var_metadata), INTENT(IN) :: varMetadata
        REAL(wp), TARGET :: varData(:,:,:,:,:)
        TYPE(t_ptr_2d), ALLOCATABLE, INTENT(INOUT) :: levelPointers(:)

        INTEGER :: nindex, nlevs, var_ref_pos, error, jk
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":getLevelPointers"

        ! Check if first dimension of array is nproma.
        ! Otherwise we got an array which is not suitable for this output scheme.
        IF (varMetadata%used_dimensions(1) /= nproma) CALL finish(routine,'1st dim is not nproma: '//TRIM(varMetadata%name))

        ! get data index
        nindex = 1
        IF(varMetadata%lcontained) nindex = varMetadata%ncontained

        ! get number of data levels
        nlevs = 1
        IF(varMetadata%ndims /= 2) nlevs = varMetadata%used_dimensions(2)

        var_ref_pos = varMetadata%ndims + 1
        IF(varMetadata%lcontained) var_ref_pos = varMetadata%var_ref_pos

        ! ALLOCATE the POINTER array
        IF(ALLOCATED(levelPointers)) THEN
            ! can we reuse the old allocation?
            IF(SIZE(levelPointers) /= nlevs) DEALLOCATE(levelPointers)
        END IF
        IF(.NOT.ALLOCATED(levelPointers)) THEN
            ! no suitable old allocation that can be reused
            ALLOCATE(levelPointers(nlevs), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        END IF

        ! get data pointers
        SELECT CASE (varMetadata%ndims)
            CASE (2)
                ! make a 3D copy of the array
                SELECT CASE(var_ref_pos)
                    CASE (1)
                        levelPointers(1)%p => varData(nindex,:,:,1,1)
                    CASE (2)
                        levelPointers(1)%p => varData(:,nindex,:,1,1)
                    CASE (3)
                        levelPointers(1)%p => varData(:,:,nindex,1,1)
                    CASE default
                        CALL finish(routine, "internal error!")
                END SELECT
            CASE (3)
                ! copy the pointer
                DO jk = 1, nlevs
                    SELECT CASE(var_ref_pos)
                        CASE (1)
                            levelPointers(jk)%p => varData(nindex,:,jk,:,1)
                        CASE (2)
                            levelPointers(jk)%p => varData(:,nindex,jk,:,1)
                        CASE (3)
                            levelPointers(jk)%p => varData(:,jk,nindex,:,1)
                        CASE (4)
                            levelPointers(jk)%p => varData(:,jk,:,nindex,1)
                        CASE default
                            CALL finish(routine, "internal error!")
                    END SELECT
                END DO
            CASE DEFAULT
              CALL finish(routine, "'"//TRIM(varMetadata%NAME)//"': "//&
                                 & TRIM(int2string(varMetadata%ndims))//"d arrays not handled yet")
        END SELECT
    END SUBROUTINE getLevelPointers

END MODULE mo_util_restart
