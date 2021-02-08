!>
!! Contains routines for asynchronous restart Output
!! --------------------------------------------------------
!!
!! Note: The synchronous implementation of the restart output can be
!!       found in the module "mo_sync_restart". See module header for
!!       more details on generated files.
!!
!! @par Revision History
!! Initial implementation by Joerg Benkenstein (2013-01-15)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_async_restart_patch_data
#ifndef NOMPI
  USE ISO_C_BINDING,                ONLY: C_PTR, C_F_POINTER, C_LOC
  USE mo_async_restart_comm_data,   ONLY: t_AsyncRestartCommData
  USE mo_exception,                 ONLY: finish
  USE mo_kind,                      ONLY: wp, dp, sp
  USE mo_impl_constants,            ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
  USE mo_cdi_constants,             ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE
  USE mo_parallel_config,           ONLY: config_restart_chunk_size => restart_chunk_size
  USE mo_run_config,                ONLY: msg_level
  USE mo_restart_var_data,          ONLY: has_valid_time_level
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_write_restart_io, timers_level
  USE mo_mpi,                       ONLY: num_work_procs, my_process_is_restart, p_mpi_wtime,  &
    &                                     p_real_sp_byte, p_real_dp_byte
#ifdef __SUNPRO_F95
  INCLUDE "mpif.h"
#else
  USE mpi,                          ONLY: MPI_ADDRESS_KIND
#endif
  USE mo_cdi,                       ONLY: streamWriteVarSlice, streamWriteVarSliceF
  USE mo_restart_patch_data,        ONLY: t_RestartPatchData
  USE mo_var_list_register_utils,   ONLY: vlr_select_restart_vars

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_AsyncPatchData

  TYPE, EXTENDS(t_RestartPatchData) :: t_AsyncPatchData
    TYPE(t_AsyncRestartCommData) :: commData
  CONTAINS
    PROCEDURE :: construct => asyncPatchData_construct
    PROCEDURE :: writeData => asyncPatchData_writeData
    PROCEDURE :: destruct  => asyncPatchData_destruct
  END TYPE t_AsyncPatchData

  CHARACTER(*), PARAMETER :: modname = 'mo_async_restart_patch_data'

CONTAINS

  ! collective across restart AND worker PEs
  SUBROUTINE asyncPatchData_construct(me, modelType, jg)
    CLASS(t_AsyncPatchData), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: modelType
    INTEGER, INTENT(IN) :: jg

    CALL me%description%init(jg)
    CALL vlr_select_restart_vars(me%varData, jg, modelType, me%restartType)
    CALL me%commData%construct(jg, me%varData)
    CALL me%transferToRestart()
  END SUBROUTINE asyncPatchData_construct

  SUBROUTINE asyncPatchData_destruct(me)
    CLASS(t_AsyncPatchData), INTENT(INOUT) :: me

    CALL me%commData%destruct()
    IF(ALLOCATED(me%varData)) DEALLOCATE(me%varData)
  END SUBROUTINE asyncPatchData_destruct

  !------------------------------------------------------------------------------------------------
  !
  ! Write restart variable list for a restart PE.
  !
  SUBROUTINE asyncPatchData_writeData(me, file_handle)
    CLASS(t_AsyncPatchData), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: file_handle
    INTEGER(MPI_ADDRESS_KIND) :: ioff(0:num_work_procs-1), bytesGet, bytesWrite
    REAL(dp), ALLOCATABLE, TARGET :: buffer_dp(:,:)
    REAL(sp), POINTER :: buffer_sp(:,:)
    INTEGER :: ichunk, nchunks, chunk_start, chunk_end, restart_chunk_size, &
      & max_nlevs, iv, nval, ierrstat, nlevs, ilev, pointCount
    REAL(dp) :: t_get, t_write
    LOGICAL :: flag_dp
    CHARACTER(*), PARAMETER :: routine = modname//':asyncPatchData_writeData'
    TYPE(c_ptr) :: cptr

    IF (.NOT. my_process_is_restart()) CALL finish(routine, 'Must be called on a restart PE!')
    t_get   = 0.d0
    t_write = 0.d0
    bytesGet = 0_mpi_address_kind
    bytesWrite = 0_mpi_address_kind
    nval = me%commData%maxLevelSize
    max_nlevs = 0
    VAR_NLEV_LOOP : DO iv = 1, SIZE(me%varData)
      IF(me%varData(iv)%p%info%ndims == 2) THEN
        nlevs = 1
      ELSE
        nlevs = me%varData(iv)%p%info%used_dimensions(2)
      ENDIF
      max_nlevs = MAX(max_nlevs, nlevs)
    END DO VAR_NLEV_LOOP
    restart_chunk_size = MERGE(MIN(config_restart_chunk_size, max_nlevs), max_nlevs, &
      &                        config_restart_chunk_size > 0)
    ALLOCATE(buffer_dp(nval,restart_chunk_size), STAT=ierrstat)
! HB: use C magic to reuse dp-buffer for sp stuff, and thus save some memory
    cptr = C_LOC(buffer_dp(1,1))
    CALL C_F_POINTER(cptr, buffer_sp, [nval,restart_chunk_size])
    IF (ierrstat /= SUCCESS) CALL finish (routine, "memory allocation failure")
    ioff(:) = 0
    ! go over the all restart variables in the associated array
    VAR_LOOP : DO iv = 1, SIZE(me%varData)
      IF (.NOT. has_valid_time_level(me%varData(iv)%p%info, me%description%id, &
        &                            me%description%nnew, me%description%nnew_rcf)) CYCLE
      IF(me%varData(iv)%p%info%ndims == 2) THEN
        nlevs = 1
      ELSE
        nlevs = me%varData(iv)%p%info%used_dimensions(2)
      ENDIF

      SELECT CASE (me%varData(iv)%p%info%hgrid)
      CASE (GRID_UNSTRUCTURED_CELL)
        pointCount = me%description%n_patch_cells_g
      CASE (GRID_UNSTRUCTURED_VERT)
        pointCount = me%description%n_patch_verts_g
      CASE (GRID_UNSTRUCTURED_EDGE)
        pointCount = me%description%n_patch_edges_g
      CASE DEFAULT
        CALL finish(routine, "Internal error: unexpected hgrid for variable "//TRIM(me%varData(iv)%p%info%name))
      END SELECT
      ! check if this is single or double precision:
      SELECT CASE(me%varData(iv)%p%info%data_type)
      CASE(REAL_T, INT_T)
        ! INTEGER fields: we write them as REAL-valued arrays
        flag_dp = .TRUE.
      CASE(SINGLE_T)
        flag_dp = .FALSE.
      CASE DEFAULT
        CALL finish(routine, "Internal error! Variable "//TRIM(me%varData(iv)%p%info%name))
      END SELECT
      ! no. of chunks of levels (each of size "restart_chunk_size"):
      nchunks = (nlevs-1)/restart_chunk_size + 1
      ! loop over all chunks (of levels)
      LEVELS : DO ichunk=1,nchunks
        chunk_start = (ichunk-1)*restart_chunk_size + 1
        chunk_end = MIN(chunk_start+restart_chunk_size-1, nlevs)
        IF (flag_dp) THEN
          CALL me%commData%collectData(me%varData(iv)%p%info%hgrid, chunk_end - chunk_start + 1, &
            &                          buffer_dp, ioff, t_get, bytesGet)
        ELSE
          CALL me%commData%collectData(me%varData(iv)%p%info%hgrid, chunk_end - chunk_start + 1, &
            &                          buffer_sp, ioff, t_get, bytesGet)
        END IF
        ! write field content into a file
        t_write = t_write - p_mpi_wtime()
        IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)
        IF (flag_dp) THEN
          DO ilev=chunk_start, chunk_end
            CALL streamWriteVarSlice(file_handle, me%varData(iv)%p%info%cdiVarID, ilev - 1, buffer_dp(:, ilev - chunk_start + 1), 0)
            bytesWrite = bytesWrite + pointCount*p_real_dp_byte
          END DO
        ELSE
          DO ilev=chunk_start, chunk_end
            CALL streamWriteVarSliceF(file_handle, me%varData(iv)%p%info%cdiVarID, ilev - 1, &
                 buffer_sp(:, ilev - chunk_start + 1), 0)
            bytesWrite = bytesWrite + pointCount*p_real_sp_byte
          END DO
        END IF
        IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
        t_write = t_write + p_mpi_wtime()
      ENDDO LEVELS

    ENDDO VAR_LOOP
    NULLIFY(buffer_sp)
    DEALLOCATE(buffer_dp, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed!')
    IF (msg_level >= 7) &
      & WRITE (0,'(10(a,f10.3))') ' Restart: Got ', REAL(bytesGet, dp)*1.d-6, ' MB, time get: ', t_get, ' s [', &
           & REAL(bytesGet, dp)*1.d-6/MAX(1.e-6_wp, t_get), ' MB/s], time write: ', t_write, ' s [', &
           & REAL(bytesWrite, dp)*1.d-6/MAX(1.e-6_wp,t_write), ' MB/s]'
  END SUBROUTINE asyncPatchData_writeData
#endif
END MODULE mo_async_restart_patch_data
