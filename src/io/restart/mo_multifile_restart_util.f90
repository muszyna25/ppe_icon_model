!> Utilities for multifile restart writing.
!! These are mostly convenience functions for accessing the restart configuration.
!!
!! Initial implementation: Nathanael Hübbe
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_multifile_restart_util
    USE mo_exception, ONLY: finish
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_std_c_lib, ONLY: strerror
    USE mo_util_file, ONLY: createSymlink
    USE mo_util_string, ONLY: int2string

    IMPLICIT NONE

    PUBLIC :: kVarName_globalCellIndex
    PUBLIC :: kVarName_globalEdgeIndex
    PUBLIC :: kVarName_globalVertIndex

    PUBLIC :: multifileRestartLinkName
    PUBLIC :: createMultifileRestartLink
    PUBLIC :: multifileAttributesPath
    PUBLIC :: multifileMetadataPath
    PUBLIC :: multifilePayloadPath

    PRIVATE

    CHARACTER(*), PARAMETER :: kVarName_globalCellIndex = "global_cell_indices"
    CHARACTER(*), PARAMETER :: kVarName_globalEdgeIndex = "global_edge_indices"
    CHARACTER(*), PARAMETER :: kVarName_globalVertIndex = "global_vert_indices"

    CHARACTER(*), PARAMETER :: modname = "mo_multifile_restart_util"

CONTAINS

    FUNCTION multifileRestartLinkName(modelType) RESULT(resultVar)
        CHARACTER(:), ALLOCATABLE :: resultVar
        CHARACTER(*), INTENT(IN) :: modelType

        resultVar = "multifile_restart_"//modelType//".mfr"
    END FUNCTION multifileRestartLinkName

    SUBROUTINE createMultifileRestartLink(filename, modelType)
        CHARACTER(*), INTENT(IN) :: filename, modelType

        INTEGER :: error
        CHARACTER(*), PARAMETER :: routine = modname//":createMultifileRestartLink"

        error = createSymlink(filename, multifileRestartLinkName(modelType))
        IF(error /= SUCCESS) CALL finish(routine, "error creating symlink to restart file: '"//strerror(error)//"'")
    END SUBROUTINE createMultifileRestartLink

    !XXX: this IS NOT the ONLY place where this path IS defined, it IS also generated/recognized IN c_restart_util.c
    FUNCTION multifileAttributesPath(multifilePath) RESULT(resultVar)
        CHARACTER(*), INTENT(IN) :: multifilePath
        CHARACTER(*), PARAMETER :: suffix = "/attributes.nc"
        CHARACTER(LEN(multifilePath) + LEN(suffix)) :: resultVar

        resultVar = multifilePath//suffix
    END FUNCTION multifileAttributesPath

    !XXX: this IS NOT the ONLY place where this path IS defined, it IS also generated/recognized IN c_restart_util.c
    FUNCTION multifileMetadataPath(multifilePath, domain) RESULT(resultVar)
        CHARACTER(:), ALLOCATABLE :: resultVar
        CHARACTER(*), INTENT(IN) :: multifilePath
        INTEGER, VALUE :: domain

        resultVar = multifilePath//"/patch"//TRIM(int2string(domain))//"_metadata"
    END FUNCTION multifileMetadataPath

    !XXX: this IS NOT the ONLY place where this path IS defined, it IS also generated/recognized IN c_restart_util.c
    FUNCTION multifilePayloadPath(multifilePath, domain, procId) RESULT(resultVar)
        CHARACTER(:), ALLOCATABLE :: resultVar
        CHARACTER(*), INTENT(IN) :: multifilePath
        INTEGER, VALUE :: domain, procId

        resultVar = multifilePath//"/patch"//TRIM(int2string(domain))//"_"//TRIM(int2string(procId))//".nc"
    END FUNCTION multifilePayloadPath

END MODULE mo_multifile_restart_util
