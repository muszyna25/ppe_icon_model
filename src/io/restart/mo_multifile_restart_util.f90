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
  USE mo_kind,           ONLY: dp, sp
  USE mo_exception,      ONLY: finish
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_std_c_lib,      ONLY: strerror
  USE mo_util_file,      ONLY: createSymlink
  USE mo_util_string,    ONLY: int2string
  USE mo_restart_util,   ONLY: alloc_string

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

  SUBROUTINE multifileRestartLinkName(modelType, resultVar)
    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: resultVar
    CHARACTER(*), INTENT(IN) :: modelType
    INTEGER :: fn_len

    fn_len = LEN_TRIM("multifile_restart_"//modelType//".mfr")
    CALL alloc_string(fn_len, resultVar)
    resultVar = "multifile_restart_"//modelType//".mfr"
  END SUBROUTINE multifileRestartLinkName

  SUBROUTINE createMultifileRestartLink(filename, modelType)
    CHARACTER(*), INTENT(IN) :: filename, modelType
    CHARACTER(:), ALLOCATABLE :: linkname
    INTEGER :: error
    CHARACTER(*), PARAMETER :: routine = modname//":createMultifileRestartLink"

    CALL multifileRestartLinkName(modelType, linkname)
    error = createSymlink(filename, linkname)
    IF(error /= SUCCESS) CALL finish(routine, "error creating symlink to restart file: '"//strerror(error)//"'")
  END SUBROUTINE createMultifileRestartLink

  !XXX: this IS NOT the ONLY place where this path IS defined, it IS also generated/recognized IN c_restart_util.c
  SUBROUTINE multifileAttributesPath(multifilePath, resultVar)
    CHARACTER(*), INTENT(IN) :: multifilePath
    CHARACTER(*), PARAMETER :: suffix = "/attributes.nc"
    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: resultVar
    INTEGER :: fn_len

    fn_len = LEN(multifilePath) + LEN(suffix)
    CALL alloc_string(fn_len, resultVar)
    resultVar = multifilePath//suffix
  END SUBROUTINE multifileAttributesPath

  !XXX: this IS NOT the ONLY place where this path IS defined, it IS also generated/recognized IN c_restart_util.c
  SUBROUTINE multifileMetadataPath(multifilePath, domain, resultVar)
    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: resultVar
    CHARACTER(*), INTENT(IN) :: multifilePath
    INTEGER, VALUE :: domain
    INTEGER :: fn_len

    fn_len=LEN_TRIM(multifilePath//"/patch"//TRIM(int2string(domain))//"_metadata")
    CALL alloc_string(fn_len, resultVar)
    resultVar = multifilePath//"/patch"//TRIM(int2string(domain))//"_metadata"
  END SUBROUTINE multifileMetadataPath

  !XXX: this IS NOT the ONLY place where this path IS defined, it IS also generated/recognized IN c_restart_util.c
  SUBROUTINE multifilePayloadPath(multifilePath, domain, procId, resultVar)
    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: resultVar
    CHARACTER(*), INTENT(IN) :: multifilePath
    INTEGER, VALUE :: domain, procId
    INTEGER :: fn_len

    fn_len = LEN_TRIM(multifilePath//"/patch"//TRIM(int2string(domain))//"_"//TRIM(int2string(procId))//".nc")
    CALL alloc_string(fn_len, resultVar)
    resultVar = multifilePath//"/patch"//TRIM(int2string(domain))//"_"//TRIM(int2string(procId))//".nc"
  END SUBROUTINE multifilePayloadPath

END MODULE mo_multifile_restart_util
