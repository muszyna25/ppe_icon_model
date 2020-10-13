!>
!! This class implements the asynchronous communication of the payload data to the restart processes.
!!
!!
!! 2018-08: Major revision / revamp / refactoring : Harald Braun (Atos SE)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

#include <handle_mpi_error.inc>
#include <icon_contiguous_defines.h>
#include <omp_definitions.inc>
MODULE mo_async_restart_comm_data
! There is no point in pretending this module is usable if NOMPI is defined.
#ifndef NOMPI

  USE ISO_C_BINDING,           ONLY: C_PTR, C_F_POINTER
  USE mo_cdi_constants,        ONLY: GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_CELL
  USE mo_reorder_info,         ONLY: t_reorder_info, mask2reorder_info, &
    &                                ri_cpy_blk2part, ri_cpy_part2whole, &
    &                                transfer_reorder_info, release_reorder_info
  USE mo_exception,            ONLY: finish
  USE mo_kind,                 ONLY: dp, i8, sp
  USE mo_model_domain,         ONLY: p_patch
  USE mo_mpi,                  ONLY: p_real_dp, p_real_dp_byte, num_work_procs, &
    &                                p_mpi_wtime, my_process_is_work, &
    &                                my_process_is_restart, p_pe_work, &
    &                                p_comm_work, p_comm_work_2_restart, &
    &                                p_comm_work_restart
  HANDLE_MPI_ERROR_USE
  USE mo_var,                  ONLY: t_var_ptr
  USE mo_parallel_config,      ONLY: config_restart_chunk_size => restart_chunk_size
  USE mo_restart_util,   ONLY: restartbcastroot
  USE mo_timer,                ONLY: timer_start, timer_stop, timer_write_restart_communication, timers_level
  USE mpi, ONLY: MPI_ADDRESS_KIND, MPI_INFO_NULL, MPI_WIN_NULL, MPI_LOCK_EXCLUSIVE, MPI_MODE_NOCHECK
#ifdef NO_MPI_RGET
  USE mpi,                     ONLY: MPI_LOCK_SHARED
#else
  USE mpi,                     ONLY: MPI_REQUEST_NULL, MPI_STATUS_IGNORE
#endif
#ifdef DEBUG
  USE mo_mpi,                  ONLY: p_pe
#endif

  IMPLICIT NONE

  PUBLIC :: t_AsyncRestartCommData

  PRIVATE

#if ICON_MPI_VERSION < 3 || ICON_MPI_VERSION == 3 && ICON_MPI_SUBVERSION == 0
  INTEGER, PARAMETER        :: max_inflight = 1
#else
  INTEGER, PARAMETER        :: max_inflight = 32
#endif

  TYPE inbuffer_t
    INTEGER :: srcPE, itype, nlev
    REAL(sp), CONTIGUOUS_POINTER :: dest_sp(:,:)
    REAL(dp), CONTIGUOUS_POINTER :: dest_dp(:,:), inbuf(:)
    TYPE(t_reorder_info), POINTER :: ri
  END TYPE inbuffer_t

  TYPE t_AsyncRestartCommData
    PRIVATE
    REAL(dp), CONTIGUOUS_POINTER :: win_buf(:)
    INTEGER :: in_use = -1, win 
    INTEGER, PUBLIC :: maxLevelSize
#ifdef NO_MPI_RGET
    INTEGER :: next_activate
#else
    INTEGER :: req_pool(max_inflight)
#endif
    TYPE(inbuffer_t) :: inbuffers(max_inflight) 
    TYPE(t_reorder_info), PUBLIC :: cells, edges, verts
  CONTAINS
    PROCEDURE, PUBLIC :: construct    => asyncRestartCommData_construct
    PROCEDURE, PRIVATE :: getPacker    => asyncRestartCommData_getPacker
    PROCEDURE :: postData_dp  => asyncRestartCommData_postData_dp
    PROCEDURE :: postData_sp  => asyncRestartCommData_postData_sp
    PROCEDURE :: postData_int => asyncRestartCommData_postData_int
    GENERIC, PUBLIC :: postData => postData_dp, postData_sp, postData_int
    PROCEDURE :: collectData_dp  => asyncRestartCommData_collectData_dp
    PROCEDURE :: collectData_sp  => asyncRestartCommData_collectData_sp
    GENERIC, PUBLIC :: collectData => collectData_dp, collectData_sp
    PROCEDURE, PUBLIC :: sync => asyncRestartCommData_sync
    PROCEDURE, PUBLIC :: destruct => asyncRestartCommData_destruct
    PROCEDURE, PRIVATE :: iterate => asyncRestartCommData_inbuffer_iterate
  END TYPE t_AsyncRestartCommData

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_async_restart_comm_data"

CONTAINS

  SUBROUTINE asyncRestartCommData_inbuffer_iterate(this, ri, nlev, elap, bytes, &
                                                   off, dest_sp, dest_dp)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: this
    TYPE(t_reorder_info), TARGET,  INTENT(IN)    :: ri
    INTEGER,                       INTENT(IN)    :: nlev
    REAL(dp),                      INTENT(INOUT) :: elap
    INTEGER(MPI_ADDRESS_KIND),     INTENT(INOUT) :: off(0:), bytes
    REAL(sp), CONTIGUOUS_POINTER,  INTENT(INOUT) :: dest_sp(:,:)
    REAL(dp), CONTIGUOUS_POINTER,  INTENT(INOUT) :: dest_dp(:,:)
    INTEGER :: srcPE

    DO srcPE = 0, num_work_procs-1
      IF(ri%pe_own(srcPE) == 0) CYCLE
      CALL acquire()
    ENDDO
    CALL waitall()
  CONTAINS

  SUBROUTINE acquire()
    CHARACTER(*), PARAMETER :: routine = modname//'::acquire'
    INTEGER :: rx_size, my_slot, ierror
    REAL(dp) :: time_mpi, time_start
    TYPE(inbuffer_t), POINTER :: inbuf_selected

    NULLIFY(inbuf_selected)
    IF (this%in_use < 0) CALL finish(routine, 'cannot use uninitialized buffer')
    time_mpi = 0._dp
    IF (this%in_use .LT. max_inflight) THEN
      IF (srcPE .GE. 0) THEN
        this%in_use = this%in_use  + 1
        my_slot = this%in_use
      ENDIF
    ELSE
      time_start = p_mpi_wtime()
#ifndef NO_MPI_RGET
      CALL MPI_Waitany(max_inflight, this%req_pool, my_slot, &
       &               MPI_STATUS_IGNORE, ierror)
      HANDLE_MPI_ERROR(ierror, 'MPI_Waitany')
#else
      my_slot=this%next_activate
      this%next_activate=MOD(my_slot,max_inflight)+1
      CALL MPI_Win_unlock(this%inbuffers(my_slot)%srcPE, this%win, ierror)
      HANDLE_MPI_ERROR(ierror, 'MPI_Win_unlock')
#endif
      time_mpi = p_mpi_wtime() - time_start
      CALL apply(this%inbuffers(my_slot))
    ENDIF
    inbuf_selected => this%inbuffers(my_slot)
    rx_size = ri%pe_own(srcPE) * nlev
    inbuf_selected%srcPE = srcPE
    inbuf_selected%nlev   = nlev
    time_start = p_mpi_wtime()
#ifndef NO_MPI_RGET
    CALL MPI_Rget(inbuf_selected%inbuf(1), rx_size, p_real_dp, &
     &            srcPE, off(srcPE), rx_size, p_real_dp, &
     &            this%win, this%req_pool(my_slot), ierror)
    HANDLE_MPI_ERROR(ierror, 'MPI_Rget')
#else
    CALL MPI_Win_lock(MPI_LOCK_SHARED, srcPE, MPI_MODE_NOCHECK, this%win, &
     &                ierror)
    HANDLE_MPI_ERROR(ierror, 'MPI_Win_lock')
    CALL MPI_Get(inbuf_selected%inbuf(1), rx_size, p_real_dp, &
     &           srcPE, off(srcPE), rx_size, p_real_dp, &
     &           this%win, ierror)
    HANDLE_MPI_ERROR(ierror, 'MPI_Get')
#endif
    off(srcPE) = off(srcPE) + INT(rx_size, MPI_ADDRESS_KIND)
    time_mpi = time_mpi + (p_mpi_wtime() - time_start)
    IF (ASSOCIATED(dest_sp)) THEN
      inbuf_selected%itype = 1
      inbuf_selected%dest_sp => dest_sp
      NULLIFY(inbuf_selected%dest_dp)
    ELSE IF (ASSOCIATED(dest_dp)) THEN
      inbuf_selected%itype = 2
      inbuf_selected%dest_dp => dest_dp
      NULLIFY(inbuf_selected%dest_sp)
    ELSE
      CALL finish(routine, ': no destination pointer')
    ENDIF
    inbuf_selected%ri => ri
    elap = elap + time_mpi
    bytes = bytes + INT(rx_size,i8) * 8_i8
  END SUBROUTINE acquire

  SUBROUTINE waitall()
    INTEGER :: remaining, err, acquired
#ifndef NO_MPI_RGET
    INTEGER :: inbuffer_map(this%in_use), i_slot, my_slot
#endif
    CHARACTER(*), PARAMETER :: routine = modname//'::inbuffer_waitall'
    REAL(dp) :: time_start, time_mpi

#ifndef NO_MPI_RGET
    DO i_slot = 1, this%in_use
      inbuffer_map(i_slot) = i_slot
    ENDDO
#endif
    time_mpi = 0._dp
    DO remaining = this%in_use, 1, -1
      time_start = p_mpi_wtime()
#ifndef NO_MPI_RGET
      CALL MPI_Waitany(remaining, this%req_pool, my_slot, &
       &               MPI_STATUS_IGNORE, err)
      HANDLE_MPI_ERROR(err, 'MPI_Waitany')
      acquired=inbuffer_map(my_slot)
#else
      acquired=MOD(this%next_activate-remaining-1+max_inflight,max_inflight)+1
      CALL MPI_Win_unlock(this%inbuffers(acquired)%srcPE, this%win, err)
      HANDLE_MPI_ERROR(err, 'MPI_Win_unlock')
#endif
      time_mpi = time_mpi + (p_mpi_wtime() - time_start)
      CALL apply(this%inbuffers(acquired))
#ifndef NO_MPI_RGET
      inbuffer_map(my_slot) = inbuffer_map(remaining)
      inbuffer_map(remaining) = acquired
      this%req_pool(my_slot) = this%req_pool(remaining)
      this%req_pool(remaining) = MPI_REQUEST_NULL
#endif
    ENDDO
    elap = elap + time_mpi
    this%in_use = 0
  END SUBROUTINE waitall

  SUBROUTINE apply(inbuf)
    TYPE(inbuffer_t), INTENT(inout) :: inbuf
    INTEGER :: n_own
    CHARACTER(*), PARAMETER :: routine = modname//'::apply'
    REAL(dp), CONTIGUOUS_POINTER :: buffer_2d(:,:)
    n_own = inbuf%ri%pe_own(inbuf%srcPE)
    buffer_2d(1:n_own, 1:inbuf%nlev) => inbuf%inbuf
!ICON_OMP PARALLEL
    SELECT CASE(inbuf%itype)
    CASE(1)
      CALL ri_cpy_part2whole(inbuf%ri, inbuf%srcPE, buffer_2d, inbuf%dest_sp)
    CASE(2)
      CALL ri_cpy_part2whole(inbuf%ri, inbuf%srcPE, buffer_2d, inbuf%dest_dp)
    END SELECT
!ICON_OMP END PARALLEL
    CALL inbuffer_reset(inbuf)
  END SUBROUTINE apply

  END SUBROUTINE asyncRestartCommData_inbuffer_iterate

  SUBROUTINE inbuffer_reset(this)
    TYPE(inbuffer_t), INTENT(INOUT) :: this

    this%srcPE = -1
    this%itype = -1
    this%nlev  = 0
    NULLIFY(this%dest_sp, this%dest_dp, this%ri)
  END SUBROUTINE inbuffer_reset

  ! collective across restart AND worker PEs
  ! returns the memory window offset for the next t_AsyncRestartCommData object
  SUBROUTINE asyncRestartCommData_construct(me, jg, var_data)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
    INTEGER, INTENT(in) :: jg
    TYPE(t_var_ptr), INTENT(IN) :: var_data(:)
    INTEGER :: i, nlevs, memWindowSize, max_nlevs, restart_chunk_size
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":asyncRestartCommData_construct"
    TYPE(t_reorder_info), POINTER :: ri
    LOGICAL :: is_restart

    is_restart = my_process_is_restart()
    CALL initReorderinfo()
    memWindowSize = 0
    max_nlevs = 0
    DO i = 1, SIZE(var_data)
      IF(var_data(i)%p%info%ndims == 2) THEN
        nlevs = 1
      ELSE
        nlevs = var_data(i)%p%info%used_dimensions(2)
      ENDIF
      max_nlevs = MAX(max_nlevs, nlevs)
      ri => me%getPacker(var_data(i)%p%info%hgrid, routine)
      memWindowSize = memWindowSize + nlevs * ri%n_own
    END DO
    CALL openMpiWindow()
    IF (is_restart) THEN
      restart_chunk_size = MERGE(MIN(config_restart_chunk_size, max_nlevs), &
        & max_nlevs, config_restart_chunk_size .GT. 0)
      IF (me%in_use >= 0) CALL finish(routine, 'buffer already initialized')
      me%in_use = 0
#ifdef NO_MPI_RGET
      me%next_activate = 1
#endif
      DO i = 1, max_inflight
#ifndef NO_MPI_RGET
        me%req_pool(i) = MPI_REQUEST_NULL
#endif
        CALL inbuffer_reset(me%inbuffers(i))
        CALL alloc_dp_for_rma(restart_chunk_size * MAX( &
          & MAXVAL(me%cells%pe_own), MAXVAL(me%verts%pe_own), &
          & MAXVAL(me%edges%pe_own)), me%inbuffers(i)%inbuf)
      ENDDO
    END IF
  CONTAINS

  SUBROUTINE initReorderinfo()
    INTEGER :: bcast_root

    IF (my_process_is_work()) THEN
      CALL mask2reorder_info(me%cells, &
        &    RESHAPE(p_patch(jg)%cells%decomp_info%owner_mask, &
        &            (/ p_patch(jg)%n_patch_cells /) ), &
        &    p_patch(jg)%n_patch_cells_g, &
        &    p_patch(jg)%cells%decomp_info%glb_index, p_comm_work)
      CALL mask2reorder_info(me%edges, &
        &    RESHAPE(p_patch(jg)%edges%decomp_info%owner_mask, &
        &            (/ p_patch(jg)%n_patch_edges /) ), &
        &    p_patch(jg)%n_patch_edges_g, &
        &    p_patch(jg)%edges%decomp_info%glb_index, p_comm_work)
      CALL mask2reorder_info(me%verts, &
        &    RESHAPE(p_patch(jg)%verts%decomp_info%owner_mask, &
        &            (/ p_patch(jg)%n_patch_verts /) ), &
        &    p_patch(jg)%n_patch_verts_g, &
        &    p_patch(jg)%verts%decomp_info%glb_index, p_comm_work)
    END IF
    bcast_root = restartbcastroot()
    CALL transfer_reorder_info(me%cells, is_restart, &
      &                        bcast_root, p_comm_work_2_restart)
    CALL transfer_reorder_info(me%edges, is_restart, &
      &                        bcast_root, p_comm_work_2_restart)
    CALL transfer_reorder_info(me%verts, is_restart, &
      &                        bcast_root, p_comm_work_2_restart)
    me%maxLevelSize = MAX(me%cells%n_glb, me%edges%n_glb, me%verts%n_glb)
  END SUBROUTINE initReorderinfo

  SUBROUTINE openMpiWindow()
    INTEGER :: rma_cache_hint, ierr
    INTEGER(MPI_ADDRESS_KIND) :: mem_bytes
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':openMpiWindow'

    me%win = MPI_WIN_NULL
    mem_bytes = INT(memWindowSize,  MPI_ADDRESS_KIND) * &
      &         INT(p_real_dp_byte, MPI_ADDRESS_KIND)
    CALL alloc_dp_for_rma(memWindowSize, me%win_buf)
    rma_cache_hint = MPI_INFO_NULL
#ifdef __xlC__
    ! IBM specific RMA hint, that we don't want window caching
    CALL MPI_Info_create(rma_cache_hint, ierr);
    HANDLE_MPI_ERROR(ierr, 'MPI_Info_create returned')
    CALL MPI_Info_set(rma_cache_hint, "IBM_win_cache", "0", ierr)
    HANDLE_MPI_ERROR(ierr, 'MPI_Info_set')
#endif
    CALL MPI_Win_create(me%win_buf, mem_bytes, INT(p_real_dp_byte), &
      &                 rma_cache_hint, p_comm_work_restart, me%win, ierr)
    HANDLE_MPI_ERROR(ierr, 'MPI_Win_create')
#ifdef __xlC__
    CALL MPI_Info_free(rma_cache_hint, ierr);
    HANDLE_MPI_ERROR(ierr, 'MPI_Info_free')
#endif
  END SUBROUTINE openMpiWindow

  SUBROUTINE alloc_dp_for_rma(n_dp, mem_ptr_dp)
    INTEGER, INTENT(IN) :: n_dp
    REAL(dp), CONTIGUOUS_POINTER, INTENT(OUT) :: mem_ptr_dp(:)
    INTEGER(MPI_ADDRESS_KIND) :: mem_bytes
    TYPE(C_PTR) :: c_mem_ptr
    INTEGER :: ierror
    CHARACTER(len=*), PARAMETER :: routine = modname//'::alloc_dp_for_rma'
    REAL(dp), TARGET, SAVE :: dummy(0)

    IF (n_dp > 0) THEN
      mem_bytes = INT(n_dp, MPI_ADDRESS_KIND) * INT(p_real_dp_byte, MPI_ADDRESS_KIND)
      CALL mpi_alloc_mem(mem_bytes, MPI_INFO_NULL, c_mem_ptr, ierror)
      HANDLE_MPI_ERROR(ierror, 'MPI_Alloc_mem')
      CALL C_F_POINTER(c_mem_ptr, mem_ptr_dp, (/ n_dp /))
    ELSE
      mem_ptr_dp => dummy
    END IF
  END SUBROUTINE alloc_dp_for_rma

  END SUBROUTINE asyncRestartCommData_construct

  ! Returns the pointer of the reorder data for the given field.
  FUNCTION asyncRestartCommData_getPacker(me, hgridType, routine) RESULT(ri_ptr)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(IN) :: me
    INTEGER, INTENT(IN) :: hgridType
    CHARACTER(LEN=*), INTENT(IN) :: routine
    TYPE(t_reorder_info), POINTER :: ri_ptr

    SELECT CASE(hgridType)
      CASE(GRID_UNSTRUCTURED_CELL)
        ri_ptr => me%cells
      CASE(GRID_UNSTRUCTURED_EDGE)
        ri_ptr => me%edges
      CASE(GRID_UNSTRUCTURED_VERT)
        ri_ptr => me%verts
      CASE default
        CALL finish(routine, "assertion failed: unexpected grid type")
    END SELECT
  END FUNCTION asyncRestartCommData_getPacker

  SUBROUTINE asyncRestartCommData_postData_dp(me, hgridType, ptr_3d, offset)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: hgridType
    REAL(dp), INTENT(IN) :: ptr_3d(:,:,:)
    INTEGER, INTENT(INOUT) :: offset
    TYPE(t_reorder_info), POINTER :: ri
    INTEGER :: ofs
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":asyncRestartCommData_postData_dp"

    ri => me%getPacker(hgridType, routine)
    ofs = offset
!ICON_OMP PARALLEL FIRSTPRIVATE(ofs)
    CALL ri_cpy_blk2part(ri, ptr_3d, me%win_buf, ofs)
!ICON_OMP MASTER
    offset = ofs
!ICON_OMP END MASTER
!ICON_OMP END PARALLEL
  END SUBROUTINE asyncRestartCommData_postData_dp

  SUBROUTINE asyncRestartCommData_postData_sp(me, hgridType, ptr_3d, offset)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: hgridType
    REAL(sp), INTENT(IN) :: ptr_3d(:,:,:)
    INTEGER, INTENT(INOUT) :: offset
    TYPE(t_reorder_info), POINTER :: ri
    INTEGER :: ofs
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":asyncRestartCommData_postData_sp"

    ri => me%getPacker(hgridType, routine)
    ofs = offset
!ICON_OMP PARALLEL FIRSTPRIVATE(ofs)
    CALL ri_cpy_blk2part(ri, ptr_3d, me%win_buf, ofs)
!ICON_OMP MASTER
    offset = ofs
!ICON_OMP END MASTER
!ICON_OMP END PARALLEL
  END SUBROUTINE asyncRestartCommData_postData_sp

  SUBROUTINE asyncRestartCommData_postData_int(me, hgridType, ptr_3d, offset)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: hgridType
    INTEGER, INTENT(in) :: ptr_3d(:,:,:)
    INTEGER, INTENT(INOUT) :: offset
    TYPE(t_reorder_info), POINTER :: ri
    INTEGER :: ofs
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":asyncRestartCommData_postData_int"

    ri => me%getPacker(hgridType, routine)
    ofs = offset
!ICON_OMP PARALLEL FIRSTPRIVATE(ofs)
    CALL ri_cpy_blk2part(ri, ptr_3d, me%win_buf, ofs)
!ICON_OMP MASTER
    offset = ofs
!ICON_OMP END MASTER
!ICON_OMP END PARALLEL
  END SUBROUTINE asyncRestartCommData_postData_int

  SUBROUTINE asyncRestartCommData_collectData_dp(me, hgridType, nlev, dest, off, &
    &                                            elapsedTime, bytes)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: hgridType, nlev
    REAL(dp), CONTIGUOUS_TARGET, INTENT(INOUT) :: dest(:,:)
    INTEGER(MPI_ADDRESS_KIND), INTENT(INOUT) :: off(0:), bytes
    REAL(dp), INTENT(INOUT) :: elapsedTime
    TYPE(t_reorder_info), POINTER :: ri
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":asyncRestartCommData_collectData_dp"
    REAL(KIND=sp), CONTIGUOUS_POINTER :: dest_sp(:,:) => NULL()
    REAL(KIND=dp), CONTIGUOUS_POINTER :: dest_dp(:,:)

    dest_dp => dest
    IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
    ri => me%getPacker(hgridType, routine)
    CALL me%iterate(ri, nlev, elapsedTime, bytes, off, dest_sp, dest_dp)
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
  END SUBROUTINE asyncRestartCommData_collectData_dp

  SUBROUTINE asyncRestartCommData_collectData_sp(me, hgridType, nlev, dest, off, &
    &                                            elapsedTime, bytes)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: hgridType, nlev
    REAL(sp), CONTIGUOUS_TARGET, INTENT(INOUT) :: dest(:,:)
    INTEGER(MPI_ADDRESS_KIND), INTENT(INOUT) :: off(0:), bytes
    REAL(dp), INTENT(INOUT) :: elapsedTime
    TYPE(t_reorder_info), POINTER :: ri
    CHARACTER(*), PARAMETER :: routine = modname//":asyncRestartCommData_collectData_sp"
    REAL(KIND=dp), CONTIGUOUS_POINTER :: dest_dp(:,:) => NULL()
    REAL(KIND=sp), CONTIGUOUS_POINTER :: dest_sp(:,:)

    dest_sp => dest
    IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
    ri => me%getPacker(hgridType, routine)
    CALL me%iterate(ri, nlev, elapsedTime, bytes, off, dest_sp, dest_dp)
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
  END SUBROUTINE asyncRestartCommData_collectData_sp

  SUBROUTINE asyncRestartCommData_destruct(this)
    CLASS(t_AsyncRestartCommData), INTENT(INOUT) :: this
    INTEGER :: i, ierr
    CHARACTER(*), PARAMETER :: routine = modname//":asyncRestartCommData_destruct"

    CALL MPI_Win_free(this%win, ierr)
    HANDLE_MPI_ERROR(ierr, 'MPI_Win_free')
    IF (SIZE(this%win_buf) > 0) THEN
      CALL MPI_Free_mem(this%win_buf, ierr)
      HANDLE_MPI_ERROR(ierr, 'MPI_Free_mem')
    END IF
    NULLIFY(this%win_buf)
    this%win = MPI_WIN_NULL
    DO i = 1, max_inflight
#ifndef NO_MPI_RGET
      this%req_pool(i) = MPI_REQUEST_NULL
#endif
      CALL inbuffer_reset(this%inbuffers(i))
      IF (my_process_is_restart()) THEN
        IF (this%in_use .NE. 0) &
          & CALL finish(routine, 'dirty buffer')
        IF (SIZE(this%inbuffers(i)%inbuf) > 0) THEN
          CALL MPI_Free_mem(this%inbuffers(i)%inbuf, ierr)
          HANDLE_MPI_ERROR(ierr, 'MPI_Free_mem')
        END IF
        NULLIFY(this%inbuffers(i)%inbuf)
      ENDIF
    END DO
    this%in_use = -1
    CALL release_reorder_info(this%cells)
    CALL release_reorder_info(this%verts)
    CALL release_reorder_info(this%edges)
  END SUBROUTINE asyncRestartCommData_destruct

  SUBROUTINE asyncRestartCommData_sync(this, after_update)
    CLASS(t_AsyncRestartCommData), INTENT(INOUT) :: this
    LOGICAL, INTENT(IN) :: after_update
    INTEGER :: ierror
    CHARACTER(*), PARAMETER :: routine = modname//":asyncRestartCommData_sync"

    IF (.NOT. my_process_is_restart()) THEN
      IF (after_update) THEN
        CALL mpi_win_unlock(p_pe_work, this%win, ierror)
        HANDLE_MPI_ERROR(ierror, 'MPI_Win_unlock')
      ELSE
        CALL mpi_win_lock(MPI_LOCK_EXCLUSIVE, p_pe_work, MPI_MODE_NOCHECK, &
             this%win, ierror)
        HANDLE_MPI_ERROR(ierror, 'mpi_win_lock')
      END IF
#ifndef NO_MPI_RGET
    ELSE
      IF (after_update) THEN
        CALL MPI_Win_unlock_all(this%win, ierror)
        HANDLE_MPI_ERROR(ierror, 'MPI_Win_unlock_all')
      ELSE
        CALL MPI_Win_lock_all(0, this%win, ierror)
        HANDLE_MPI_ERROR(ierror, 'MPI_Win_start')
      END IF
#endif
    END IF
  END SUBROUTINE asyncRestartCommData_sync
#endif
END MODULE mo_async_restart_comm_data
