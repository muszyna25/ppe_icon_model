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

    USE mo_exception, ONLY: finish
    USE mo_cdi, ONLY: CDI_UNDEFID, GRID_UNSTRUCTURED, gridCreate, gridDefNvertex, gridDefXname, gridDefXlongname, gridDefXunits, &
                    & gridDefYname, gridDefYlongname, gridDefYunits, zaxisCreate, zaxisDefLevels
    USE mo_kind, ONLY: wp

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: createHgrids
    PUBLIC :: defineVAxis
    PUBLIC :: defineSingleLevelAxis

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

    SUBROUTINE createHgrids(cellCount, outCellGridId, &
                           &vertexCount, outVertexGridId, &
                           &edgeCount, outEdgeGridId, cellType)
        INTEGER, VALUE :: cellCount, vertexCount, edgeCount, cellType
        INTEGER, INTENT(OUT) :: outCellGridId, outVertexGridId, outEdgeGridId

        outCellGridId = create_cdi_hgrid_def(cellCount, cellType, &
                                            &'clon', 'center longitude', 'radian', &
                                            &'clat', 'center latitude', 'radian')

        outVertexGridId = create_cdi_hgrid_def(vertexCount, 9 - cellType, &
                                              &'vlon', 'vertex longitude', 'radian', &
                                              &'vlat', 'vertex latitude', 'radian')

        outEdgeGridId = create_cdi_hgrid_def(edgeCount, 4, &
                                            &'elon', 'edge midpoint longitude', 'radian',  &
                                            &'elat', 'edge midpoint latitude', 'radian')

    END SUBROUTINE createHgrids

    ! If no opt_levelValues are given, this defaults to numbering the levels from 1 to levelCount.
    INTEGER FUNCTION defineVAxis(cdiAxisType, levelCount, opt_levelValues) RESULT(RESULT)
        INTEGER, VALUE :: cdiAxisType, levelCount
        REAL(KIND = wp), OPTIONAL, INTENT(IN) :: opt_levelValues(:)

        INTEGER :: i

        RESULT = zaxisCreate(cdiAxisType, levelCount)
        IF(PRESENT(opt_levelValues)) THEN
            CALL zaxisDefLevels(RESULT, opt_levelValues)
        ELSE
            CALL zaxisDefLevels(RESULT, [(REAL(i, wp), i = 1, levelCount)])
        END IF
    END FUNCTION defineVAxis

    INTEGER FUNCTION defineSingleLevelAxis(cdiAxisType, levelValue) RESULT(RESULT)
        INTEGER, VALUE :: cdiAxisType
        REAL(KIND = wp), VALUE :: levelValue

        RESULT = defineVAxis(cdiAxisType, 1, [levelValue])
    END FUNCTION defineSingleLevelAxis

END MODULE mo_util_restart
