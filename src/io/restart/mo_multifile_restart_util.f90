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
  USE mo_io_config,      ONLY: restartWritingParameters, ALL_WORKERS_INVOLVED
  USE mo_mpi,            ONLY: num_work_procs, p_comm_work_restart, p_comm_rank
  USE mo_fortran_tools,  ONLY: t_ptr_2d, t_ptr_2d_sp, t_ptr_2d_int

  IMPLICIT NONE

  PUBLIC :: vNames_glbIdx
  PUBLIC :: multifileRestartLinkName
  PUBLIC :: createMultifileRestartLink
  PUBLIC :: multifileAttributesPath
  PUBLIC :: multifilePayloadPath
  PUBLIC :: isAsync, restartProcCount
  PUBLIC :: rBuddy, rGroup, commonBuf_t, dataPtrs_t
  PUBLIC :: iAmRestartMaster, iAmRestartWriter

  PRIVATE

  TYPE :: commonBuf_t
    REAL(dp), POINTER :: d(:)
    REAL(sp), POINTER :: s(:)
    INTEGER,  POINTER :: i(:)
  END TYPE commonBuf_t

  TYPE :: dataPtrs_t
    TYPE(t_ptr_2d),     ALLOCATABLE :: d(:)
    TYPE(t_ptr_2d_sp),  ALLOCATABLE :: s(:)
    TYPE(t_ptr_2d_int), ALLOCATABLE :: i(:)
  CONTAINS
    PROCEDURE :: free => dataPtrs_t_free
  END TYPE dataPtrs_t

  CHARACTER(*), PARAMETER :: vNames_glbIdx(3) = (/"global_cell_indices", &
    &                                             "global_edge_indices", &
    &                                             "global_vert_indices"/)
  CHARACTER(*), PARAMETER :: modname = "mo_multifile_restart_util"

CONTAINS

  SUBROUTINE dataPtrs_t_free(this)
    CLASS(dataPtrs_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%d)) DEALLOCATE(this%d)
    IF (ALLOCATED(this%s)) DEALLOCATE(this%s)
    IF (ALLOCATED(this%i)) DEALLOCATE(this%i)
  END SUBROUTINE dataPtrs_t_free

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
  SUBROUTINE multifilePayloadPath(multifilePath, domain, procId, resultVar)
    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: resultVar
    CHARACTER(*), INTENT(IN) :: multifilePath
    INTEGER, VALUE :: domain, procId
    INTEGER :: fn_len

    fn_len = LEN_TRIM(multifilePath//"/patch"//TRIM(int2string(domain))//"_"//TRIM(int2string(procId))//".nc")
    CALL alloc_string(fn_len, resultVar)
    resultVar = multifilePath//"/patch"//TRIM(int2string(domain))//"_"//TRIM(int2string(procId))//".nc"
  END SUBROUTINE multifilePayloadPath

  INTEGER FUNCTION restartProcCount() RESULT(resultVar)
    CALL restartWritingParameters(opt_restartProcCount = resultVar)
    IF (resultVar == ALL_WORKERS_INVOLVED)  resultVar = num_work_procs
  END FUNCTION restartProcCount

  LOGICAL FUNCTION isAsync() RESULT(resultVar)
    CALL restartWritingParameters(opt_lDedicatedProcMode = resultVar)
  END FUNCTION isAsync

  INTEGER FUNCTION rBuddy(pe_in, nr_in, nw_in) RESULT(res)
    INTEGER, OPTIONAL,INTENT(IN) :: pe_in, nr_in, nw_in
    INTEGER :: pe, nr, nw

    nw = num_work_procs
    nr = restartProcCount()
    pe = p_comm_rank(p_comm_work_restart)
    IF (PRESENT(nr_in)) nr = nr_in
    IF (PRESENT(pe_in)) pe = pe_in
    IF (PRESENT(nw_in)) nw = nw_in
    res = rGroup(pe_in=pe, nr_in=nr, nw_in=nw)
    IF (isAsync() .AND. .NOT.PRESENT(nr_in)) THEN
      res = nw + res
    ELSE
      res = CEILING(REAL(res * nw)/REAL(nr))
    END IF
  END FUNCTION rBuddy

  INTEGER FUNCTION rGroup(pe_in, nr_in, nw_in) RESULT(res)
    INTEGER, OPTIONAL,INTENT(IN) :: pe_in, nr_in, nw_in
    INTEGER :: pe, nr, nw

    nw = num_work_procs
    nr = restartProcCount()
    pe = p_comm_rank(p_comm_work_restart)
    IF (PRESENT(nr_in)) nr = nr_in
    IF (PRESENT(pe_in)) pe = pe_in
    IF (PRESENT(nw_in)) nw = nw_in
    res = MERGE(pe - nw, &
      &         FLOOR(REAL(nr * pe)/REAL(nw)), &
      &         isAsync() .AND. pe .GE. nw)
  END FUNCTION rGroup

  LOGICAL FUNCTION iAmRestartMaster() RESULT(resultVar)

    resultVar = rGroup() == rGroup(pe_in=0) .AND. &
      &         iAmRestartWriter()
  END FUNCTION iAmRestartMaster

  LOGICAL FUNCTION iAmRestartWriter() RESULT(resultVar)
    INTEGER :: the_pe

    the_pe = p_comm_rank(p_comm_work_restart)
    resultVar = rBuddy() == the_pe
  END FUNCTION iAmRestartWriter

END MODULE mo_multifile_restart_util
