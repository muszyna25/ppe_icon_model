!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This MODULE provides a TYPE that binds all the relevant CDI IDs connected to a single file together.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_cdi_ids
    USE mo_cdi, ONLY: zaxisCreate, zaxisDefLevels, streamOpenWrite, vlistCreate, taxisCreate, vlistDefTaxis, TAXIS_ABSOLUTE, &
                    & vlistDefVar, vlistDefVarDatatype, vlistDefVarName, vlistDefVarLongname, vlistDefVarUnits, &
                    & vlistDefVarMissval, TIME_VARIABLE, DATATYPE_FLT64, taxisDefVdate, taxisDefVtime, cdiEncodeDate, &
                    & cdiEncodeTime, streamDefTimestep, FILETYPE_NC2, FILETYPE_NC4, CDI_UNDEFID, gridCreate, gridDefNvertex, &
                    & gridDefXname, gridDefXlongname, gridDefXunits, gridDefYname, gridDefYlongname, gridDefYunits, &
                    & GRID_UNSTRUCTURED
    USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT
    USE mo_cdi_constants, ONLY: ZA_COUNT, ZA_HYBRID, ZA_HYBRID_HALF, ZA_LAKE_BOTTOM, ZA_LAKE_BOTTOM_HALF, ZA_MIX_LAYER, &
                              & ZA_SEDIMENT_BOTTOM_TW_HALF, cdi_zaxis_types
    USE mtime, ONLY: datetime
    USE mo_exception, ONLY: finish, message
    USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH, SUCCESS
    USE mo_kind, ONLY: wp
    USE mo_util_cdi, ONLY: cdiGetStringError
    USE mo_var_metadata_types, ONLY: t_var_metadata
    IMPLICIT NONE

    PRIVATE

    PUBLIC :: t_Vgrid
    PUBLIC :: t_CdiIds
    PUBLIC :: set_vertical_grid

    ! TYPE t_Vgrid contains the data of a vertical grid definition.
    TYPE t_Vgrid
        INTEGER :: type
        REAL(wp), ALLOCATABLE :: levels(:)
    END type t_Vgrid

    ! TYPE t_CdiIds IS just a simple container for all the different CDI IDs connected to a single restart file.
    TYPE t_CdiIds
        INTEGER :: file, vlist, taxis, hgrids(3), vgrids(ZA_COUNT)
    CONTAINS
        PROCEDURE :: init => restartCdiIds_init
        PROCEDURE :: openRestartAndCreateIds => restartCdiIds_openRestartAndCreateIds
        PROCEDURE :: finalizeVlist => restartCdiIds_finalizeVlist
        PROCEDURE :: defineVariable => restartCdiIds_defineVariable
        PROCEDURE :: closeAndDestroyIds => restartCdiIds_closeAndDestroyIds
    END TYPE t_CdiIds

    ! This takes a buffer for grid definitions, to which one entry IS appended by incrementing the count of used elements that's passed as well.
    INTERFACE set_vertical_grid
        MODULE PROCEDURE set_vertical_grid_array
        MODULE PROCEDURE set_vertical_grid_counted
        MODULE PROCEDURE set_vertical_grid_single
    END INTERFACE

    CHARACTER(*), PARAMETER :: modname = "mo_cdi_ids"

CONTAINS

    ! If no opt_levelValues are given, this defaults to numbering the levels from 1 to levelCount.
    INTEGER FUNCTION defineVAxis(cdiAxisType, levelValues) RESULT(resultVar)
        INTEGER, VALUE :: cdiAxisType
        REAL(KIND = wp), INTENT(IN) :: levelValues(:)

        INTEGER :: levelCount

        levelCount = SIZE(levelValues, 1)
        resultVar = zaxisCreate(cdiAxisType, levelCount)
        CALL zaxisDefLevels(resultVar, levelValues)
    END FUNCTION defineVAxis

    SUBROUTINE createVgrids(axisIds, gridDescriptions, opt_vct)
        INTEGER, INTENT(INOUT) :: axisIds(ZA_COUNT)
        TYPE(t_Vgrid), INTENT(IN) :: gridDescriptions(:)
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
        CLASS(t_CdiIds), INTENT(INOUT) :: me

        me%file = CDI_UNDEFID
        me%vlist = CDI_UNDEFID
        me%taxis = CDI_UNDEFID
        me%hgrids(:) = CDI_UNDEFID
        me%vgrids(:) = CDI_UNDEFID
    END SUBROUTINE restartCdiIds_init

    ! Creates a horizontal grid definition from the given parameters, returns the new CDI gridId.
    INTEGER FUNCTION create_cdi_hgrid_def(iCnt, iNVert, cNameX, cLNameX, cUnitsX, &
        &                                               cNameY, cLNameY, cUnitsY) RESULT(resultVar)

        INTEGER, VALUE :: iCnt, iNVert
        CHARACTER(LEN = *), INTENT(IN) :: cNameX, cLNameX, cUnitsX
        CHARACTER(LEN = *), INTENT(IN) :: cNameY, cLNameY, cUnitsY
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":create_cdi_hgrid_def"

        resultVar = gridCreate(GRID_UNSTRUCTURED, iCnt)
        IF(resultVar == CDI_UNDEFID) CALL finish(routine, "error creating CDI grid")
        CALL gridDefNvertex(resultVar, iNVert)

        CALL gridDefXname(resultVar, TRIM(cNameX))
        CALL gridDefXlongname(resultVar, TRIM(cLNameX))
        CALL gridDefXunits(resultVar, TRIM(cUnitsX))

        CALL gridDefYname(resultVar, TRIM(cNameY))
        CALL gridDefYlongname(resultVar, TRIM(cLNameY))
        CALL gridDefYunits(resultVar, TRIM(cUnitsY))
    END FUNCTION create_cdi_hgrid_def

    FUNCTION createHgrids(cellCount, vertexCount, edgeCount, cellType) RESULT(resultVar)
        INTEGER, VALUE :: cellCount, vertexCount, edgeCount, cellType
        INTEGER :: resultVar(3)

        resultVar(GRID_UNSTRUCTURED_CELL) = create_cdi_hgrid_def(cellCount, cellType, &
                                                             &'clon', 'center longitude', 'radian', &
                                                             &'clat', 'center latitude', 'radian')

        resultVar(GRID_UNSTRUCTURED_VERT) = create_cdi_hgrid_def(vertexCount, 9 - cellType, &
                                                             &'vlon', 'vertex longitude', 'radian', &
                                                             &'vlat', 'vertex latitude', 'radian')

        resultVar(GRID_UNSTRUCTURED_EDGE) = create_cdi_hgrid_def(edgeCount, 4, &
                                                             &'elon', 'edge midpoint longitude', 'radian', &
                                                             &'elat', 'edge midpoint latitude', 'radian')

    END FUNCTION createHgrids

    SUBROUTINE restartCdiIds_openRestartAndCreateIds(me, filename, restartType, cellCount, vertCount, edgeCount, cellType, &
                                                    &vgridDefs, opt_vct)
        CLASS(t_CdiIds), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: filename
        INTEGER, VALUE :: restartType, cellCount, vertCount, edgeCount, cellType
        TYPE(t_Vgrid), INTENT(IN) :: vgridDefs(:)
        REAL(wp), INTENT(IN), OPTIONAL :: opt_vct(:)

        CHARACTER(LEN = MAX_CHAR_LENGTH) :: cdiErrorText
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartCdiIds_openRestartAndCreateIds"

        ! check whether the given filetype IS OK
        SELECT CASE(restartType)
            CASE(FILETYPE_NC2, FILETYPE_NC4)
                ! These are ok, both formats allow for files larger than 2GiB.
            CASE DEFAULT
                CALL finish(routine, "unsupported restart file type")
        END SELECT

        ! open the file
        me%file = streamOpenWrite(filename, restartType)

        IF(me%file < 0) THEN
            CALL cdiGetStringError(me%file, cdiErrorText)
            CALL message('', TRIM(cdiErrorText))
            CALL finish(routine, 'open failed on '//filename)
        END IF

        ! create the CDI IDs we need

        ! 1. vlist
        me%vlist = vlistCreate()

        ! 3. horizontal grids
        me%hgrids = createHgrids(cellCount, vertCount, edgeCount, cellType)

        ! 4. vertical grids
        CALL createVgrids(me%vgrids, vgridDefs, opt_vct)

        ! 5. time axis (always absolute time for restart files)
        me%taxis = taxisCreate(TAXIS_ABSOLUTE)
        CALL vlistDefTaxis(me%vlist, me%taxis)
    END SUBROUTINE restartCdiIds_openRestartAndCreateIds

    ! CDI reqires the vlist of a stream to be set before the timestep can be defined, which IS why we combine these two operations into one SUBROUTINE.
    SUBROUTINE restartCdiIds_finalizeVlist(me, this_datetime)
        CLASS(t_CdiIds), INTENT(INOUT) :: me
        TYPE(datetime), POINTER, INTENT(IN) :: this_datetime

        INTEGER :: trash

        ! define the vlist, so that we are allowed to define the timestep
        CALL streamDefVlist(me%file, me%vlist)

        ! define the timestep
        CALL taxisDefVdate(me%taxis, &
             &             cdiEncodeDate(INT(this_datetime%date%year), this_datetime%date%month, this_datetime%date%day))
        CALL taxisDefVtime(me%taxis, &
             &             cdiEncodeTime(this_datetime%time%hour, this_datetime%time%minute, this_datetime%time%second))
        trash = streamDefTimestep(me%file, 0)
    END SUBROUTINE restartCdiIds_finalizeVlist

    ! Encapsulates the CDI calls to define a variable.
    SUBROUTINE restartCdiIds_defineVariable(me, info)
        CLASS(t_CdiIds), INTENT(IN) :: me
        TYPE(t_var_metadata), INTENT(INOUT) :: info

        INTEGER :: varId, gridId, zaxisId
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartCdiIds_defineVariable"

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
        IF(info%lmiss) CALL vlistDefVarMissval(me%vlist, varId, info%missval%rval)
    END SUBROUTINE restartCdiIds_defineVariable

    SUBROUTINE restartCdiIds_closeAndDestroyIds(me)
        CLASS(t_CdiIds), INTENT(INOUT) :: me

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
        TYPE(t_Vgrid), INTENT(INOUT) :: gridDefinitions(:)
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
        TYPE(t_Vgrid), INTENT(INOUT) :: gridDefinitions(:)
        INTEGER, INTENT(INOUT) :: gridCount
        INTEGER, VALUE :: TYPE, levelCount

        INTEGER :: i

        CALL set_vertical_grid(gridDefinitions, gridCount, TYPE, [(REAL(i, wp), i = 1, levelCount)])
    END SUBROUTINE set_vertical_grid_counted

    SUBROUTINE set_vertical_grid_single(gridDefinitions, gridCount, TYPE, levelValue)
        TYPE(t_Vgrid), INTENT(INOUT) :: gridDefinitions(:)
        INTEGER, INTENT(INOUT) :: gridCount
        INTEGER, VALUE :: TYPE
        REAL(wp), VALUE :: levelValue

        CALL set_vertical_grid(gridDefinitions, gridCount, TYPE, [levelValue])
    END SUBROUTINE set_vertical_grid_single

END MODULE mo_cdi_ids
