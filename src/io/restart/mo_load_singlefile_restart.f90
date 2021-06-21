!> Module for reading single file restart files
!!
!! Note: The asynchronous implementation of the restart output can be
!!       found in the module "mo_async_restart"
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
#include "omp_definitions.inc"
MODULE mo_load_singlefile_restart
  USE mo_cdi,                ONLY: vlistInqTaxis, taxisInqVdate, taxisInqVtime, streamOpenRead, &
    &                              streamInqVlist, vlistNvars, vlistInqVarName, streamClose
  USE mo_exception,          ONLY: message, warning, finish
  USE mo_impl_constants,     ONLY: SINGLE_T, REAL_T, INT_T
  USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE
  USE mo_kind,               ONLY: dp, sp
  USE mo_model_domain,       ONLY: t_patch
  USE mo_mpi,                ONLY: my_process_is_mpi_workroot, p_bcast, p_comm_work
  USE mo_restart_nml_and_att,ONLY: getAttributesForRestarting
  USE mo_key_value_store,    ONLY: t_key_value_store
  USE mo_restart_util,       ONLY: restartSymlinkName
  USE mo_restart_var_data,   ONLY: get_var_3d_ptr
  USE mo_var,                ONLY: t_var_ptr
  USE mo_timer, ONLY: timer_start, timer_stop, timer_load_restart_io, timers_level
  USE mo_util_string,        ONLY: int2string
  USE mo_read_netcdf_distributed, ONLY: t_distrib_read_data, distrib_nf_open, &
    & distrib_read, distrib_nf_close, idx_lvl_blk
  USE mo_fortran_tools,      ONLY: t_ptr_3d, t_ptr_3d_int, t_ptr_3d_sp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: singlefileCheckRestartFiles, singlefileReadPatch

  CHARACTER(*), PARAMETER :: modname = "mo_load_singlefile_restart"

CONTAINS

  ! This IS just an existence check for the different patch files, a
  ! warning IS generated each domain for which no corresponding
  ! restart file IS found.
  !
  ! XXX: The code that I found READ the restart attributes from all
  ! the files, *appending* them to the list of attributes. Since all
  ! restart files contain the same attributes, this ONLY served to
  ! duplicate them. The current code just ignores the attributes IN
  ! all but the first file, just like the original code ignored the
  ! namelists IN all but the first file.  However, it might be a good
  ! idea to add some consistency checking on the attributes/namelists
  ! of the other files.
  SUBROUTINE singlefileCheckRestartFiles(modelType)
    CHARACTER(*), INTENT(IN) :: modelType
    INTEGER :: jg, n_dom
    LOGICAL :: lexists
    CHARACTER(:), ALLOCATABLE :: filename
    TYPE(t_key_value_store), POINTER :: restartAttributes
    CHARACTER(*), PARAMETER :: routine = modname//":singlefileCheckRestartFiles"

    CALL getAttributesForRestarting(restartAttributes)
    CALL restartAttributes%get('n_dom', n_dom)
    ! check whether we have all the restart files we expect
    DO jg = 1, n_dom
        CALL restartSymlinkName(modelType, jg, filename, n_dom)
        INQUIRE(file = filename, exist = lexists)
        IF (lexists) THEN
            CALL message(routine, "found restart file at '"//filename//"'")
        ELSE
            CALL warning(routine, 'domain '//TRIM(int2string(jg))//' is not active at restart time')
        END IF
    END DO
  END SUBROUTINE singlefileCheckRestartFiles

  SUBROUTINE singlefileReadPatch(varData, modelType, p_patch, opt_ndom)
    TYPE(t_var_ptr), INTENT(in) :: varData(:)
    CHARACTER(*), INTENT(IN) :: modelType
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
    INTEGER, INTENT(IN) :: opt_ndom
    TYPE(t_ptr_3d) :: r(1)
    TYPE(t_ptr_3d_sp) :: s(1)
    TYPE(t_ptr_3d_int) :: i(1)
    CHARACTER(len=80) :: vname
    CHARACTER(:), ALLOCATABLE :: restart_filename
    INTEGER :: fID, dfID, vlID, vID, tID, nV, nread, vIdx, dt
    TYPE(t_distrib_read_data), POINTER :: dio
    LOGICAL :: is_root
    CHARACTER(*), PARAMETER :: routine = modname // ":singlefileReadPatch"

    CALL restartSymlinkName(modelType, p_patch%id, restart_filename, opt_ndom)
    is_root = my_process_is_mpi_workroot()
    IF(is_root) THEN
      WRITE(0, "(a)") "streamOpenRead " // restart_filename
      fID  = streamOpenRead(restart_filename)
      IF(fID < 0) CALL finish(routine, "could not open file '" // restart_filename // "'")
      vlID = streamInqVlist(fID)
      tID = vlistInqTaxis(vlID)
      WRITE(0, "(a)") routine // ': Read restart for: ' // &
        & TRIM(int2string(taxisInqVdate(tID))) // "T" // &
        & TRIM(int2string(taxisInqVtime(tID))) // "Z" // ' from ' // restart_filename
      nV = vlistNvars(vlID)
    END IF
    dfID = distrib_nf_open(restart_filename) 
    CALL p_bcast(nV, 0, comm=p_comm_work)
    nread = 0
    DO vID = 0, nV - 1
      IF(is_root) CALL vlistInqVarName(vlID, vID, vname)
      CALL p_bcast(vname, 0, comm=p_comm_work)
      !lookup the varData corresponding to this variable AND get the corresponding DATA pointers
      DO vIdx = 1, SIZE(varData)
        IF(vname == varData(vIdx)%p%info%NAME) EXIT
      END DO
      IF (vIdx > SIZE(varData)) THEN
        CALL warning(routine, "variable '" // TRIM(vname) // "' in restart file, but not in list of restart variables")
        CYCLE
      END IF
      SELECT CASE(varData(vIdx)%p%info%hgrid)
      CASE (GRID_UNSTRUCTURED_CELL)
        dio => p_patch%cells%dist_io_data
      CASE (GRID_UNSTRUCTURED_VERT)
        dio => p_patch%verts%dist_io_data
      CASE (GRID_UNSTRUCTURED_EDGE)
        dio => p_patch%edges%dist_io_data
      CASE default
        CALL finish(routine, 'unknown grid type')
      END SELECT
      dt = varData(vIdx)%p%info%data_type
      IF (timers_level >= 7) CALL timer_start(timer_load_restart_io)
      SELECT CASE(dt)
      CASE(REAL_T, INT_T)
        IF (dt .EQ. INT_T) THEN
          CALL get_var_3d_ptr(varData(vIdx)%p, i(1)%p)
          ALLOCATE(r(1)%p(SIZE(i(1)%p, 1), SIZE(i(1)%p, 2), SIZE(i(1)%p, 3)))
        ELSE
          CALL get_var_3d_ptr(varData(vIdx)%p, r(1)%p)
        END IF
        IF (SIZE(r(1)%p, 3) .GT. 0) r(1)%p(:,:,SIZE(r(1)%p, 3)) = 0._dp
        CALL distrib_read(dfID, vname, r, (/dio/), &
          & edim=(/SIZE(r(1)%p, 2)/), dimo=idx_lvl_blk)
        IF (SIZE(r(1)%p) .GT. 0 .AND. dt .EQ. INT_T) THEN
!ICON_OMP PARALLEL WORKSHARE
          i(1)%p(:,:,:) = INT(r(1)%p(:,:,:))
!ICON_OMP END PARALLEL WORKSHARE
        END IF
        IF (dt .EQ. INT_T) DEALLOCATE(r(1)%p)
      CASE(SINGLE_T)
        CALL get_var_3d_ptr(varData(vIdx)%p, s(1)%p)
        IF (SIZE(s(1)%p, 3) .GT. 0) s(1)%p(:,:,SIZE(s(1)%p, 3)) = 0._sp
        CALL distrib_read(dfID, vname, s, (/dio/), &
          & edim=(/SIZE(s(1)%p, 2)/), dimo=idx_lvl_blk)
      CASE DEFAULT
        CALL finish(routine, "Internal error! Variable " // TRIM(vname))
      END SELECT
      IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
      nread = nread + 1
#ifdef DEBUG
      IF (is_root) WRITE (0,*) ' ... read ' // TRIM(vname)
#endif
    END DO
    IF (is_root) WRITE (0,'(a,i5,a)') ' ... read ', nread, ' variables'
    IF (is_root)  CALL streamClose(fID)
    CALL distrib_nf_close(dfID)
  END SUBROUTINE singlefileReadPatch

END MODULE mo_load_singlefile_restart
