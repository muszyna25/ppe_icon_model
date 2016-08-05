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
                    & vlistDefVarLongname, vlistDefVarUnits, vlistDefVarMissval, TIME_VARIABLE, DATATYPE_FLT64
    USE mo_cdi_constants, ONLY: ZA_HYBRID, ZA_HYBRID_HALF, ZA_LAKE_BOTTOM, ZA_MIX_LAYER, ZA_LAKE_BOTTOM_HALF, &
                              & ZA_SEDIMENT_BOTTOM_TW_HALF, ZA_COUNT, cdi_zaxis_types, GRID_UNSTRUCTURED_CELL, &
                              & GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_COUNT
    USE mo_impl_constants, ONLY: SUCCESS, MAX_CHAR_LENGTH
    USE mo_io_restart_attributes, ONLY: t_RestartAttributeList
    USE mo_io_restart_namelist, ONLY: RestartNamelist_writeToFile
    USE mo_kind, ONLY: wp
    USE mo_util_cdi, ONLY: cdiGetStringError
    USE mo_var_metadata_types, ONLY: t_var_metadata

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: t_v_grid
    PUBLIC :: t_restart_cdi_ids

    PUBLIC :: openRestartAndCreateIds
    PUBLIC :: defineVariable
    PUBLIC :: closeAndDestroyIds
    PUBLIC :: set_vertical_grid

    ! TYPE t_v_grid contains the data of a vertical grid definition.
    TYPE t_v_grid
        INTEGER :: type
        REAL(wp), ALLOCATABLE :: levels(:)
    END type t_v_grid

    ! TYPE t_restart_cdi_ids IS just a simple container for all the different CDI IDs connected to a single restart file.
    TYPE t_restart_cdi_ids
        INTEGER :: file, vlist, taxis, hgrids(GRID_UNSTRUCTURED_COUNT), vgrids(ZA_COUNT)
    END TYPE t_restart_cdi_ids

    ! This takes a buffer for grid definitions, to which one entry IS appended by incrementing the count of used elements that's passed as well.
    INTERFACE set_vertical_grid
        MODULE PROCEDURE set_vertical_grid_array
        MODULE PROCEDURE set_vertical_grid_counted
        MODULE PROCEDURE set_vertical_grid_single
    END INTERFACE

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_util_restart"

CONTAINS

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
            IF(cdi_zaxis_types(gridDescriptions(i)%TYPE) == CDI_UNDEFID) CALL finish(routine, "no CDI zaxis TYPE defined for vgrid")

            gridId = defineVAxis(cdi_zaxis_types(gridDescriptions(i)%TYPE), gridDescriptions(i)%levels)
            IF(gridId == CDI_UNDEFID) CALL finish(routine, "defineVAxis() returned an error")
            axisIds(gridDescriptions(i)%TYPE) = gridId

            SELECT CASE (gridDescriptions(i)%type)
                CASE (ZA_HYBRID)
                    IF (.NOT.PRESENT(opt_vct)) CYCLE
                    nlevp1 = SIZE(gridDescriptions(i)%levels, 1) + 1
                    CALL zaxisDefVct(gridId, 2*nlevp1, opt_vct(1:2*nlevp1))

                CASE (ZA_HYBRID_HALF)
                    IF (.NOT.PRESENT(opt_vct)) CYCLE
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

    SUBROUTINE openRestartAndCreateIds(cdiIds, filename, restartType, restartAttributes, cellCount, vertCount, edgeCount, &
                                      &cellType, vgridDefs, opt_vct)
        TYPE(t_restart_cdi_ids), INTENT(INOUT) :: cdiIds
        CHARACTER(LEN = *), INTENT(IN) :: filename
        INTEGER, VALUE :: restartType, cellCount, vertCount, edgeCount, cellType
        TYPE(t_RestartAttributeList), POINTER, INTENT(INOUT) :: restartAttributes
        TYPE(t_v_grid), INTENT(IN) :: vgridDefs(:)
        REAL(wp), INTENT(IN), OPTIONAL :: opt_vct(:)

        CHARACTER(LEN = MAX_CHAR_LENGTH) :: cdiErrorText
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":openRestartAndCreateIds"

        ! open the file
        cdiIds%file = streamOpenWrite(filename, restartType)

        IF(cdiIds%file < 0) THEN
            CALL cdiGetStringError(cdiIds%file, cdiErrorText)
            WRITE(message_text,'(a)') TRIM(cdiErrorText)
            CALL message('',message_text)
            CALL finish(routine, 'open failed on '//filename)
        END IF

        ! create the CDI IDs we need

        ! 1. vlist
        cdiIds%vlist = vlistCreate()

        ! 2. global attributes
        CALL RestartNamelist_writeToFile(cdiIds%vlist)
        CALL restartAttributes%writeToFile(cdiIds%vlist)

        ! 3. horizontal grids
        cdiIds%hgrids = createHgrids(cellCount, vertCount, edgeCount, cellType)

        ! 4. vertical grids
        CALL createVgrids(cdiIds%vgrids, vgridDefs, opt_vct)

        ! 5. time axis (always absolute time for restart files)
        cdiIds%taxis = taxisCreate(TAXIS_ABSOLUTE)
        CALL vlistDefTaxis(cdiIds%vlist, cdiIds%taxis)
    END SUBROUTINE openRestartAndCreateIds

    ! Encapsulates the CDI calls to define a variable.
    ! lIsInteger AND lIsLogical reflect the TYPE of the variable, IF neither IS set, the variable IS assumed to be of TYPE REAL.
    SUBROUTINE defineVariable(cdiIds, info, lIsInteger, lIsLogical)
        TYPE(t_restart_cdi_ids), INTENT(IN) :: cdiIds
        TYPE(t_var_metadata), INTENT(INOUT) :: info
        LOGICAL, VALUE :: lIsInteger, lIslogical

        INTEGER :: varId, gridId, zaxisId
        REAL(wp) :: casted_missval
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":defineVariable"

        IF(lIsInteger.AND.lIsLogical) THEN
            CALL finish(routine, "assertion failed: attempt to define a variable both as INTEGER and LOGICAL")
        END IF

        ! get the horizontal grid ID
        gridId = info%cdiGridID
        SELECT CASE (info%hgrid)
            CASE(GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE)
                gridId = cdiIds%hgrids(info%hgrid)
        END SELECT
        IF (gridId == CDI_UNDEFID) CALL finish(routine, 'Grid type not defined for field '//TRIM(info%name))

        ! get the vertical axis ID
        zaxisId = info%cdiZaxisID
        if(zaxisId < 0) zaxisId = cdiIds%vgrids(info%vgrid)
        IF (zaxisId == CDI_UNDEFID) CALL finish(routine, 'Z axis not defined for field '//TRIM(info%name))

        ! define the variable with the required info
        varId = vlistDefVar(cdiIds%vlist, gridId, zaxisId, TIME_VARIABLE)
        IF (varID == CDI_UNDEFID) CALL finish(routine, 'error WHILE defining CDI variable "'//TRIM(info%name)//'"')
        info%cdiVarID = varId
        CALL vlistDefVarDatatype(cdiIds%vlist, varId, DATATYPE_FLT64)
        CALL vlistDefVarName(cdiIds%vlist, varId, TRIM(info%name))

        ! then add the three optional fields
        IF(info%cf%long_name /= '') CALL vlistDefVarLongname(cdiIds%vlist, varId, TRIM(info%cf%long_name))
        IF(info%cf%units /= '') CALL vlistDefVarUnits(cdiIds%vlist, varId, TRIM(info%cf%units))
        IF(info%lmiss) THEN
            casted_missval = info%missval%rval
            IF(lIsInteger) casted_missval = REAL(info%missval%ival, wp)
            IF(lIsLogical) THEN
                casted_missval = 0.0_wp
                IF(info%missval%lval) casted_missval = 1.0_wp
            ENDIF
            CALL vlistDefVarMissval(cdiIds%vlist, varId, casted_missval)
        ENDIF
    END SUBROUTINE defineVariable

    SUBROUTINE closeAndDestroyIds(cdiIds)
        TYPE(t_restart_cdi_ids), INTENT(INOUT) :: cdiIds

        INTEGER :: i

        ! close/destroy all open CDI IDs
        IF(cdiIds%file /= CDI_UNDEFID) CALL streamClose(cdiIds%file)
        IF(cdiIds%vlist /= CDI_UNDEFID) CALL vlistDestroy(cdiIds%vlist)
        IF(cdiIds%taxis /= CDI_UNDEFID) CALL taxisDestroy(cdiIds%taxis)
        DO i = 1, SIZE(cdiIds%hgrids, 1)
            IF(cdiIds%hgrids(i) /= CDI_UNDEFID) CALL gridDestroy(cdiIds%hgrids(i))
        END DO
        DO i = 1, SIZE(cdiIds%vgrids, 1)
            IF(cdiIds%vgrids(i) /= CDI_UNDEFID) CALL zaxisDestroy(cdiIds%vgrids(i))
        END DO

        ! reset the IDs
        cdiIds%file = CDI_UNDEFID
        cdiIds%vlist = CDI_UNDEFID
        cdiIds%taxis = CDI_UNDEFID
        cdiIds%hgrids(:) = CDI_UNDEFID
        cdiIds%vgrids(:) = CDI_UNDEFID
    END SUBROUTINE closeAndDestroyIds

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

END MODULE mo_util_restart
