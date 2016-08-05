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
                    & gridDefYname, gridDefYlongname, gridDefYunits

    public :: create_cdi_hgrid_def

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

END MODULE mo_util_restart
