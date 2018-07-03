!>
!! This class implements the asynchronous communication of the payload data to the restart processes.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_async_restart_comm_data
! There is no point in pretending this module is usable if NOMPI is defined.
#ifndef NOMPI

  USE ISO_C_BINDING,           ONLY: C_PTR, C_INTPTR_T, C_F_POINTER
  USE mo_async_restart_packer, ONLY: t_AsyncRestartPacker
  USE mo_cdi_constants,        ONLY: GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_CELL
  USE mo_decomposition_tools,  ONLY: t_grid_domain_decomp_info
  USE mo_exception,            ONLY: finish, message, em_warn
  USE mo_fortran_tools,        ONLY: t_ptr_2d, t_ptr_2d_sp, t_ptr_2d_int, ensureSize, no_copy
  USE mo_kind,                 ONLY: dp, i8, sp
  USE mo_model_domain,         ONLY: p_patch
  USE mo_mpi,                  ONLY: p_real_dp, p_comm_work_restart, p_pe_work, num_work_procs, &
    &                                p_mpi_wtime, my_process_is_work
  USE mo_restart_var_data,     ONLY: t_RestartVarData
  USE mo_timer,                ONLY: timer_start, timer_stop, timer_write_restart_communication, timers_level
  USE mo_util_string,          ONLY: int2string
  USE mpi,                     ONLY: MPI_ADDRESS_KIND, MPI_INFO_NULL, MPI_LOCK_SHARED,       &
    &                                MPI_MODE_NOCHECK, MPI_MODE_NOSTORE, MPI_MODE_NOPRECEDE, &
    &                                MPI_WIN_NULL, MPI_LOCK_EXCLUSIVE, MPI_SUCCESS,          &
    &                                MPI_REQUEST_NULL, MPI_MODE_NOSUCCEED, MPI_MODE_NOPUT,   &
    &                                MPI_STATUS_IGNORE, MPI_UNDEFINED, MPI_INTEGER, MPI_ERRORS_RETURN
#ifdef DEBUG
  USE mo_mpi,                  ONLY: p_pe
#endif

  IMPLICIT NONE

  PUBLIC :: t_AsyncRestartCommData

  PRIVATE

  INTEGER, PARAMETER        :: max_inflight = 32
  INTEGER, PARAMETER        :: reset_hard   = 1

  TYPE inbuffer_t
    PRIVATE
    INTEGER                             :: srcProc, itype, levelCount
    REAL(sp),                   POINTER :: dest_sp(:,:)
    REAL(dp),                   POINTER :: dest_dp(:,:)
    TYPE(t_AsyncRestartPacker), POINTER :: reorderData
    REAL(dp),                   POINTER :: inbuffer(:)
  CONTAINS
    PROCEDURE :: reset => asyncRestartCommData_inbuffer_reset
  END TYPE inbuffer_t


  TYPE inbufhandler_t
    PRIVATE
    INTEGER                   :: in_use, domain, win
    INTEGER                   :: req_pool(max_inflight)
    TYPE(inbuffer_t)          :: inbuffers(max_inflight)
    INTEGER                   :: server_group, client_group
    LOGICAL                   :: handshaked, hello, is_init
  CONTAINS
    PROCEDURE         :: prepare   => asyncRestartCommData_inbuffer_prepare
    PROCEDURE         :: clear     => asyncRestartCommData_inbuffer_clear
    PROCEDURE         :: aquire    => asyncRestartCommData_inbuffer_aquire
    PROCEDURE         :: waitall   => asyncRestartCommData_inbuffer_waitall
    PROCEDURE         :: apply     => asyncRestartCommData_inbuffer_apply
    PROCEDURE         :: iterate   => asyncRestartCommData_inbuffer_iterate
    PROCEDURE         :: handshake => asyncRestartCommData_inbuffer_handshake
    PROCEDURE         :: sync      => asyncRestartCommData_inbuffer_sync
  END TYPE inbufhandler_t


  ! this combines all the DATA that's relevant for transfering the
  ! payload DATA from the worker PEs to the restart PEs
  TYPE t_AsyncRestartCommData
    PRIVATE
    ! DATA for remote memory access
    INTEGER           :: win
    REAL(dp), POINTER :: win_buf(:)
    TYPE(inbufhandler_t), ALLOCATABLE :: inbufhandler
    ! reorder data
    TYPE(t_AsyncRestartPacker), PUBLIC :: cells
    TYPE(t_AsyncRestartPacker), PUBLIC :: edges
    TYPE(t_AsyncRestartPacker), PUBLIC :: verts
  CONTAINS
    PROCEDURE :: construct    => asyncRestartCommData_construct
    ! called to get the required buffer SIZE on the restart processes:
    PROCEDURE :: maxLevelSize => asyncRestartCommData_maxLevelSize
    ! RETURN the relevant t_AsyncRestartPacker object
    PROCEDURE :: getPacker    => asyncRestartCommData_getPacker
    ! called by the compute processes to write their data to their memory window:
    PROCEDURE :: postData_dp  => asyncRestartCommData_postData_dp  
    PROCEDURE :: postData_sp  => asyncRestartCommData_postData_sp  
    PROCEDURE :: postData_int => asyncRestartCommData_postData_int
    GENERIC, PUBLIC :: postData => postData_dp, postData_sp, postData_int
    ! called by the restart processes to fetch the DATA from the
    ! compute processes:
    PROCEDURE :: collectData_dp  => asyncRestartCommData_collectData_dp 
    PROCEDURE :: collectData_sp  => asyncRestartCommData_collectData_sp 
    GENERIC, PUBLIC :: collectData => collectData_dp, collectData_sp
    PROCEDURE, PUBLIC :: sync => asyncRestartCommData_sync
    PROCEDURE :: destruct     => asyncRestartCommData_destruct
  END TYPE t_AsyncRestartCommData

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_async_restart_comm_data"

CONTAINS

  SUBROUTINE asyncRestartCommData_inbuffer_iterate(this,                       &
                                                   reorderData_in, levelCount, &
                                                   elapsedTime, bytesFetched,  &
                                                   offsets,                    &
                                                   dest_sp_in, dest_dp_in)
    CLASS(inbufhandler_t), TARGET,       INTENT(INOUT) :: this
    TYPE(t_AsyncRestartPacker), POINTER, INTENT(IN)    :: reorderData_in
    INTEGER,                             INTENT(IN)    :: levelCount
    REAL(KIND=dp),                       INTENT(INOUT) :: elapsedTime
    INTEGER(KIND=i8),                    INTENT(INOUT) :: bytesFetched
    INTEGER(KIND=MPI_ADDRESS_KIND),      INTENT(INOUT) :: offsets(0:num_work_procs-1)
    REAL(KIND=sp), DIMENSION(:,:), POINTER, INTENT(INOUT) :: dest_sp_in
    REAL(KIND=dp), DIMENSION(:,:), POINTER, INTENT(INOUT) :: dest_dp_in
    CHARACTER(LEN=*), PARAMETER                        :: procedure_name = &
                                modname//'::inbufferhandler_iterate'
    INTEGER                                            :: srcProc
    REAL(KIND=sp), DIMENSION(:,:), POINTER                :: dest_sp
    REAL(KIND=dp), DIMENSION(:,:), POINTER                :: dest_dp

    IF(ASSOCIATED(dest_sp_in)) THEN
      dest_sp => dest_sp_in
      NULLIFY(dest_dp)
    ENDIF
    IF(ASSOCIATED(dest_dp_in)) THEN
      dest_dp => dest_dp_in
      NULLIFY(dest_sp)
    ENDIF
    DO srcProc = 0, num_work_procs-1
      IF(reorderData_in%pe_own(srcProc) == 0) CYCLE
      CALL this%aquire(srcProc, levelCount, elapsedTime, bytesFetched, &
       &               offsets, reorderData_in, dest_sp, dest_dp)
    ENDDO
    CALL this%waitall(elapsedTime)
  END SUBROUTINE asyncRestartCommData_inbuffer_iterate

  SUBROUTINE asyncRestartCommData_inbuffer_reset(this, hard)
    CLASS(inbuffer_t),        INTENT(INOUT)          :: this
    INTEGER,                  INTENT(IN),   OPTIONAL :: hard

    this%srcProc = -1
    this%itype        = -1
    this%levelCount   = 0
    NULLIFY(this%dest_sp)
    NULLIFY(this%dest_dp)
    NULLIFY(this%reorderData)
    IF (PRESENT(hard)) THEN
      IF (hard .EQ. reset_hard) THEN
        IF (ASSOCIATED(this%inbuffer)) THEN
          IF (ALLOCATED(this%inbuffer)) THEN
            DEALLOCATE(this%inbuffer)
            NULLIFY(this%inbuffer)
          ELSE
            NULLIFY(this%inbuffer)
          ENDIF
        ELSE
          NULLIFY(this%inbuffer)
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE asyncRestartCommData_inbuffer_reset

  SUBROUTINE asyncRestartCommData_inbuffer_waitall(this, elapsedTime)
    CLASS(inbufhandler_t), TARGET,  INTENT(INOUT) :: this
    REAL(dp),                       INTENT(INOUT) :: elapsedTime
    INTEGER                                       :: i_slot, my_slot, remaining, &
     &                                               err, inbuffer_map(this%in_use)
    CHARACTER(LEN=*), PARAMETER                   :: procedure_name = &
     &                          modname//'::inbufferhandler_waitall'
    REAL(KIND=dp)                                 :: time_start, time_mpi

    DO i_slot = 1, this%in_use
      inbuffer_map(i_slot) = i_slot
    ENDDO
    time_mpi = 0._dp
    DO remaining = this%in_use, 1, -1
      time_start = p_mpi_wtime()
      CALL MPI_Waitany(remaining, this%req_pool, my_slot, &
       &               MPI_STATUS_IGNORE, err)
      CALL ar_checkmpi(err, procedure_name, ': MPI_Waitany failed')
      time_mpi = time_mpi + (p_mpi_wtime() - time_start)
      CALL this%apply(inbuffer_map(my_slot))
      CALL this%inbuffers(inbuffer_map(my_slot))%reset()
      i_slot = inbuffer_map(my_slot)
      inbuffer_map(my_slot) = inbuffer_map(remaining) 
      inbuffer_map(remaining) = i_slot
      this%req_pool(my_slot) = this%req_pool(remaining)
      this%req_pool(remaining) = MPI_REQUEST_NULL
    ENDDO
    elapsedTime = elapsedTime + time_mpi
    this%in_use = 0
  END SUBROUTINE asyncRestartCommData_inbuffer_waitall

  SUBROUTINE asyncRestartCommData_inbuffer_apply(this, my_slot)
    CLASS(inbufhandler_t), TARGET,  INTENT(INOUT) :: this
    TYPE(inbuffer_t), POINTER                     :: my_inbuf
    INTEGER,                        INTENT(IN)    :: my_slot
    INTEGER                                       :: reorderOffset, &
                                                     ilevel, reorderOffsetNext
    CHARACTER(LEN=*), PARAMETER                   :: procedure_name = &
     &                          modname//'::inbufferhandler_apply'

    reorderOffset = 1
    my_inbuf => this%inbuffers(my_slot)
    DO ilevel = 1, my_inbuf%levelCount
      reorderOffsetNext = reorderOffset - 1 + &
        &                  my_inbuf%reorderData%pe_own(my_inbuf%srcProc)
      SELECT CASE(my_inbuf%itype)
      CASE(1)
        CALL my_inbuf%reorderData%unpackLevelFromPe( &
          &       ilevel, my_inbuf%srcProc,          &
          &       my_inbuf%inbuffer(reorderOffset:reorderOffsetNext), &
          &       my_inbuf%dest_sp)
      CASE(2)
        CALL my_inbuf%reorderData%unpackLevelFromPe( &
          &       ilevel, my_inbuf%srcProc,          &
          &       my_inbuf%inbuffer(reorderOffset:reorderOffsetNext), &
          &       my_inbuf%dest_dp)
      reorderOffset = reorderOffsetNext + 1
      END SELECT
    ENDDO
  END SUBROUTINE asyncRestartCommData_inbuffer_apply

  SUBROUTINE ar_checkmpi(mpierr, caller, msg)
    INTEGER,                        INTENT(IN) :: mpierr
    CHARACTER(LEN=*),               INTENT(IN) :: caller, msg

    IF (MPI_SUCCESS .NE. mpierr) THEN
      CALL finish(caller, msg)
    ENDIF
  END SUBROUTINE ar_checkmpi

  SUBROUTINE asyncRestartCommData_inbuffer_aquire(this, srcProc_in, levelCount_in, elapsedTime, bytesFetched, &
   &                                              offsets, reorderData_in, dest_sp, dest_dp)
    CLASS(inbufhandler_t), TARGET,          INTENT(INOUT)   :: this
    REAL(KIND=sp), DIMENSION(:,:), POINTER, INTENT(INOUT)   :: dest_sp
    REAL(KIND=dp), DIMENSION(:,:), POINTER, INTENT(INOUT)   :: dest_dp
    INTEGER,                                INTENT(IN)      :: srcProc_in, levelCount_in
    REAL(KIND=dp),                          INTENT(INOUT)   :: elapsedTime
    INTEGER(KIND=i8),                       INTENT(INOUT)   :: bytesFetched
    INTEGER(KIND=MPI_ADDRESS_KIND),         INTENT(INOUT)   :: offsets(0:num_work_procs-1)
    TYPE(t_AsyncRestartPacker),    POINTER                  :: reorderData_in
    CHARACTER(LEN=*), PARAMETER                             :: procedure_name = &
     &                          modname//'::inbufferhandler_aquire'
    INTEGER                                                 :: rx_size, my_slot, mpierr
    REAL(KIND=dp) :: time_mpi,time_start
    TYPE(inbuffer_t),             POINTER                   :: inbuf_selected
    NULLIFY(inbuf_selected)
    IF (.NOT.this%is_init) THEN
      CALL ar_checkmpi(-99, procedure_name, ': cannot use uninitialized buffer')
    ENDIF
    time_mpi = 0._dp
    IF (this%in_use .LT. max_inflight) THEN
      IF (srcProc_in .GE. 0) THEN
        this%in_use = this%in_use  + 1
        my_slot = this%in_use
      ENDIF
    ELSE
      time_start = p_mpi_wtime()
      CALL MPI_Waitany(max_inflight, this%req_pool, my_slot, &
       &               MPI_STATUS_IGNORE, mpierr)
      CALL ar_checkmpi(mpierr, procedure_name, ': MPI_Waitany failed')
      time_mpi = p_mpi_wtime() - time_start
      CALL this%apply(my_slot)
      CALL this%inbuffers(my_slot)%reset()
    ENDIF
    inbuf_selected => this%inbuffers(my_slot)
    rx_size = reorderData_in%pe_own(srcProc_in) * levelCount_in
    CALL ensureSize(inbuf_selected%inbuffer, rx_size, no_copy)
    inbuf_selected%srcProc = srcProc_in
    inbuf_selected%levelCount   = levelCount_in
    time_start = p_mpi_wtime()
    CALL MPI_Rget(inbuf_selected%inbuffer(1),      rx_size, p_real_dp, &
     &            srcProc_in, offsets(srcProc_in), rx_size, p_real_dp, &
     &            this%win, this%req_pool(my_slot), mpierr)
    CALL ar_checkmpi(mpierr, procedure_name, ': MPI_Rget failed')
    offsets(srcProc_in) = offsets(srcProc_in) + INT(rx_size, MPI_ADDRESS_KIND)
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
      CALL ar_checkmpi(-99, procedure_name, ': no destiation pointer')
    ENDIF
    IF(ASSOCIATED(reorderData_in)) THEN
      inbuf_selected%reorderData => reorderData_in
    ELSE
      CALL ar_checkmpi(-99, procedure_name, ': no reorderData pointer')
    ENDIF
    elapsedTime = elapsedTime + time_mpi
    bytesFetched = bytesFetched + INT(rx_size,i8) * 8_i8
  END SUBROUTINE asyncRestartCommData_inbuffer_aquire

  SUBROUTINE asyncRestartCommData_inbuffer_prepare(this, win_in, domain)
    CLASS(inbufhandler_t), TARGET, INTENT(INOUT) :: this
    INTEGER,               TARGET, INTENT(IN)    :: win_in
    INTEGER,                       INTENT(IN)    :: domain
    CHARACTER(LEN=*), PARAMETER                  :: procedure_name = &
     &                          modname//'::inbufferhandler_prepare'
    INTEGER                                      :: i
    IF (this%is_init) THEN
      CALL ar_checkmpi(-99, procedure_name, ': buffer already initialized')
    ENDIF
    this%win        = win_in
    this%in_use     = 0
    this%handshaked = .false.
    this%domain     = domain
    DO i = 1, max_inflight
      this%req_pool(i) = MPI_REQUEST_NULL
      CALL this%inbuffers(i)%reset()
      NULLIFY(this%inbuffers(i)%inbuffer)
    ENDDO
    this%is_init = .true.
    this%hello   = .true.
  END SUBROUTINE asyncRestartCommData_inbuffer_prepare

  SUBROUTINE asyncRestartCommData_inbuffer_clear(this)
    CLASS(inbufhandler_t), TARGET, INTENT(INOUT) :: this
    CHARACTER(LEN=*), PARAMETER :: procedure_name = &
     &                          modname//'::inbufferhandler_clear'
    INTEGER                                      :: i
    IF (this%in_use .NE. 0) THEN
      CALL ar_checkmpi(-99, procedure_name, ' dirty buffer')
    ENDIF
    this%win = MPI_WIN_NULL
    DO i = 1, max_inflight
      this%req_pool(i) = MPI_REQUEST_NULL
      CALL this%inbuffers(i)%reset(reset_hard)
    ENDDO
    this%is_init = .false.
  END SUBROUTINE asyncRestartCommData_inbuffer_clear

  ! Opens an MPI memory window for the given amount of double
  ! values, returning both the ALLOCATED buffer AND the MPI window
  ! handle.
  SUBROUTINE openMpiWindow(doubleCount, communicator, mem_ptr_dp, the_win)
    INTEGER(KIND=MPI_ADDRESS_KIND), VALUE :: doubleCount   ! the requested memory window SIZE IN doubles
    INTEGER, VALUE :: communicator
    REAL(dp), POINTER, INTENT(OUT) :: mem_ptr_dp(:) ! this returns a fortran POINTER to the memory buffer
    INTEGER, INTENT(OUT) :: the_win ! this returns the handle to the MPI window
#ifdef __xlC__
    INTEGER :: rma_cache_hint
#endif
    INTEGER :: nbytes_real, mpi_error
    INTEGER(KIND=MPI_ADDRESS_KIND) :: mem_bytes
    TYPE(C_PTR) :: c_mem_ptr
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':openMpiWindow'

#ifdef DEBUG
    WRITE(nerr, FORMAT_VALS3)routine, ' is called for p_pe=', p_pe
#endif
    the_win = MPI_WIN_NULL
    ! doubleCount is calculated as number of variables above, get number of bytes
    ! get the amount of bytes per REAL*8 variable (as used in MPI communication)
    CALL MPI_Type_extent(p_real_dp, nbytes_real, mpi_error)
    CALL ar_checkmpi(mpi_error, routine, 'MPI_Type_extent returned error '//TRIM(int2string(mpi_error)))
    ! for the restart PEs the amount of memory needed is 0 - allocate at least 1 word there:
    mem_bytes = MAX(doubleCount, 1_i8)*INT(nbytes_real, i8)
    ! allocate amount of memory needed with MPI_Alloc_mem
    ! 
    ! Depending on wether the Fortran 2003 C interoperability features
    ! are available, one needs to use non-standard language extensions
    ! for calls from Fortran, namely Cray Pointers, since
    ! MPI_Alloc_mem wants a C pointer argument.
    !
    ! see, for example: http://www.lrz.de/services/software/parallel/mpi/onesided/
    ! TYPE(C_PTR) and INTEGER(KIND=MPI_ADDRESS_KIND) do NOT necessarily have the same size!!!
    ! So check if at least C_INTPTR_T and MPI_ADDRESS_KIND are the same, else we may get
    ! into deep, deep troubles!
    ! There is still a slight probability that TYPE(C_PTR) does not have the size indicated
    ! by C_INTPTR_T since the standard only requires C_INTPTR_T is big enough to hold pointers
    ! (so it may be bigger than a pointer), but I hope no vendor screws up its ISO_C_BINDING
    ! in such a way!!!
    ! If C_INTPTR_T<=0, this type is not defined and we can't do this check, of course.
    IF(C_INTPTR_T > 0 .AND. C_INTPTR_T /= MPI_ADDRESS_KIND) THEN
        CALL finish(routine, 'C_INTPTR_T /= MPI_ADDRESS_KIND, too dangerous to proceed!')
    END IF
    CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, c_mem_ptr, mpi_error)
    CALL ar_checkmpi(mpi_error, routine, 'MPI_Alloc_mem returned error '//TRIM(int2string(mpi_error)))
    NULLIFY(mem_ptr_dp)
    CALL C_F_POINTER(c_mem_ptr, mem_ptr_dp, [doubleCount] )
#ifdef __xlC__
    ! IBM specific RMA hint, that we don't want window caching
    rma_cache_hint = MPI_INFO_NULL
    CALL MPI_Info_create(rma_cache_hint, mpi_error);
    CALL ar_checkmpi(mpi_error, routine, 'MPI_Info_create returned error '//TRIM(int2string(mpi_error)))
    CALL MPI_Info_set(rma_cache_hint, "IBM_win_cache", "0", mpi_error)
    CALL ar_checkmpi(mpi_error, routine, 'MPI_Info_set returned error '//TRIM(int2string(mpi_error)))
#endif
    ! create memory window for communication
    mem_ptr_dp(:) = 0._dp
    CALL MPI_Win_create(mem_ptr_dp, mem_bytes, nbytes_real, MPI_INFO_NULL, communicator, the_win, mpi_error)
    CALL ar_checkmpi(mpi_error, routine, 'MPI_Win_create returned error '//TRIM(int2string(mpi_error)))
#ifdef __xlC__
    CALL MPI_Info_free(rma_cache_hint, mpi_error);
    ar_checkmpi(mpi_error, CALL finish(routine, 'MPI_Info_free returned error '//TRIM(int2string(mpi_error)))
#endif
  END SUBROUTINE openMpiWindow

  SUBROUTINE closeMpiWindow(windowPtr, mpiWindow)
    REAL(dp), POINTER, INTENT(INOUT) :: windowPtr(:)
    INTEGER, INTENT(INOUT) :: mpiWindow
    INTEGER :: mpi_error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":closeMpiWindow"

    CALL MPI_Win_free(mpiWindow, mpi_error)
    CALL warnMpiError('MPI_Win_free')
    CALL MPI_Free_mem(windowPtr, mpi_error)
    CALL warnMpiError('MPI_Free_mem')
    mpiWindow = MPI_WIN_NULL
    NULLIFY(windowPtr)
  CONTAINS
    SUBROUTINE warnMpiError(mpiCall)
      CHARACTER(LEN = *), INTENT(IN) :: mpiCall
      IF(mpi_error /= MPI_SUCCESS) THEN
        CALL message(routine, mpiCall//"() returned error "//TRIM(int2string(mpi_error)), &
          &          level = em_warn, all_print = .TRUE.)
      END IF
    END SUBROUTINE warnMpiError
  END SUBROUTINE closeMpiWindow

  ! collective across restart AND worker PEs
  ! returns the memory window offset for the next t_AsyncRestartCommData object
  SUBROUTINE asyncRestartCommData_construct(me, jg, var_data)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
    INTEGER, VALUE :: jg
    TYPE(t_RestartVarData), POINTER, INTENT(IN) :: var_data(:)
    INTEGER :: iv, nlevs
    TYPE(t_grid_domain_decomp_info) :: dummyInfo
    INTEGER(KIND = MPI_ADDRESS_KIND) :: memWindowSize
    TYPE(t_AsyncRestartPacker), POINTER :: reorderData
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":asyncRestartCommData_construct"

    ! initialize the reorder info
    IF(my_process_is_work()) THEN
      CALL me%cells%construct(p_patch(jg)%n_patch_cells_g, p_patch(jg)%n_patch_cells, &
        &                     p_patch(jg)%cells%decomp_info)
      CALL me%edges%construct(p_patch(jg)%n_patch_edges_g, p_patch(jg)%n_patch_edges, &
        &                     p_patch(jg)%edges%decomp_info)
      CALL me%verts%construct(p_patch(jg)%n_patch_verts_g, p_patch(jg)%n_patch_verts, &
        &                     p_patch(jg)%verts%decomp_info)
    ELSE
      !must not access p_patch on the restart processes,
      !nevertheless the restart processes need to take part in
      !the collective calls
      CALL me%cells%construct(0, 0, dummyInfo)
      CALL me%edges%construct(0, 0, dummyInfo)
      CALL me%verts%construct(0, 0, dummyInfo)
    END IF
    ! initialize the memory window offsets
    memWindowSize = 0
    IF(ASSOCIATED(var_data)) THEN
      ! compute the SIZE of the memory window
      DO iv = 1, SIZE(var_data)
        nlevs = 1
        IF(var_data(iv)%info%ndims /= 2) nlevs = var_data(iv)%info%used_dimensions(2)
        reorderData => me%getPacker(var_data(iv)%info%hgrid, routine)
        memWindowSize = memWindowSize + INT(nlevs*reorderData%n_own, i8)
      END DO
    END IF
    ! actually open the memory window
    CALL openMpiWindow(memWindowSize, p_comm_work_restart, me%win_buf, me%win)
    ALLOCATE(me%inbufhandler)
    me%inbufhandler%is_init = .false.
    CALL me%inbufhandler%prepare(me%win,jg)
    CALL me%inbufhandler%handshake(.NOT.my_process_is_work())
  END SUBROUTINE asyncRestartCommData_construct

  INTEGER FUNCTION asyncRestartCommData_maxLevelSize(me) RESULT(resultVar)
    CLASS(t_AsyncRestartCommData), INTENT(IN) :: me

    resultVar = MAX(me%cells%n_glb, me%edges%n_glb, me%verts%n_glb)
  END FUNCTION asyncRestartCommData_maxLevelSize

  ! Returns the pointer of the reorder data for the given field.
  FUNCTION asyncRestartCommData_getPacker(me, hgridType, routine) RESULT(resultVar)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(IN) :: me
    INTEGER, VALUE :: hgridType
    CHARACTER(LEN=*), INTENT(IN) :: routine
    TYPE(t_AsyncRestartPacker), POINTER :: resultVar

    NULLIFY(resultVar)
    SELECT CASE(hgridType)
      CASE(GRID_UNSTRUCTURED_CELL)
        resultVar => me%cells
      CASE(GRID_UNSTRUCTURED_EDGE)
        resultVar => me%edges
      CASE(GRID_UNSTRUCTURED_VERT)
        resultVar => me%verts
      CASE default
        CALL finish(routine, "assertion failed: unexpected grid type")
    END SELECT
  END FUNCTION asyncRestartCommData_getPacker

  SUBROUTINE asyncRestartCommData_postData_dp(me, hgridType, dataPointers, offset)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
    INTEGER, VALUE :: hgridType
    TYPE(t_ptr_2d), ALLOCATABLE, INTENT(IN) :: dataPointers(:)
    INTEGER(i8), INTENT(INOUT) :: offset
    TYPE(t_AsyncRestartPacker), POINTER :: reorderData
    INTEGER :: i
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":asyncRestartCommData_postData_dp"

    reorderData => me%getPacker(hgridType, routine)
    DO i = 1, SIZE(dataPointers)
      CALL reorderData%packLevel(dataPointers(i)%p, me%win_buf, offset)
    END DO
  END SUBROUTINE asyncRestartCommData_postData_dp

  SUBROUTINE asyncRestartCommData_postData_sp(me, hgridType, dataPointers, offset)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
    INTEGER, VALUE :: hgridType
    TYPE(t_ptr_2d_sp), ALLOCATABLE, INTENT(IN) :: dataPointers(:)
    INTEGER(i8), INTENT(INOUT) :: offset
    TYPE(t_AsyncRestartPacker), POINTER :: reorderData
    INTEGER :: i
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":asyncRestartCommData_postData_sp"

    reorderData => me%getPacker(hgridType, routine)
    DO i = 1, SIZE(dataPointers)
        CALL reorderData%packLevel(dataPointers(i)%p, me%win_buf, offset)
    END DO
  END SUBROUTINE asyncRestartCommData_postData_sp

  SUBROUTINE asyncRestartCommData_postData_int(me, hgridType, dataPointers, offset)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
    INTEGER, VALUE :: hgridType
    TYPE(t_ptr_2d_int), ALLOCATABLE, INTENT(IN) :: dataPointers(:)
    INTEGER(i8), INTENT(INOUT) :: offset
    TYPE(t_AsyncRestartPacker), POINTER :: reorderData
    INTEGER :: i
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":asyncRestartCommData_postData_int"

    reorderData => me%getPacker(hgridType, routine)
    DO i = 1, SIZE(dataPointers)
        CALL reorderData%packLevel(dataPointers(i)%p, me%win_buf, offset)
    END DO
  END SUBROUTINE asyncRestartCommData_postData_int

  SUBROUTINE asyncRestartCommData_collectData_dp(me, hgridType, levelCount, dest, offsets, &
    &                                            elapsedTime, bytesFetched)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
    INTEGER, VALUE :: hgridType, levelCount
    REAL(dp), TARGET, INTENT(INOUT) :: dest(:,:)
    INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(INOUT) :: offsets(0:num_work_procs-1)
    REAL(dp), INTENT(INOUT) :: elapsedTime
    INTEGER(i8), INTENT(INOUT) :: bytesFetched
    TYPE(t_AsyncRestartPacker), POINTER :: reorderData
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":asyncRestartCommData_collectData_dp"
    REAL(KIND=sp), DIMENSION(:,:), POINTER                :: dest_sp
    REAL(KIND=dp), DIMENSION(:,:), POINTER                :: dest_dp

    dest_dp => dest
    NULLIFY(dest_sp)
    IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
    reorderData => me%getPacker(hgridType, routine)
    CALL me%inbufhandler%iterate(reorderData, levelCount, elapsedTime, bytesFetched, offsets, dest_sp, dest_dp)
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
  END SUBROUTINE asyncRestartCommData_collectData_dp

  SUBROUTINE asyncRestartCommData_collectData_sp(me, hgridType, levelCount, dest, offsets, &
    &                                            elapsedTime, bytesFetched)
    CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
    INTEGER, VALUE :: hgridType, levelCount
    REAL(sp), TARGET, INTENT(INOUT) :: dest(:,:)
    INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(INOUT) :: offsets(0:num_work_procs-1)
    REAL(dp), INTENT(INOUT) :: elapsedTime
    INTEGER(i8), INTENT(INOUT) :: bytesFetched
    TYPE(t_AsyncRestartPacker), POINTER :: reorderData
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":asyncRestartCommData_collectData_sp"
    REAL(KIND=sp), DIMENSION(:,:), POINTER                :: dest_sp
    REAL(KIND=dp), DIMENSION(:,:), POINTER                :: dest_dp

    dest_sp => dest
    NULLIFY(dest_dp)
    IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
    reorderData => me%getPacker(hgridType, routine)
    CALL me%inbufhandler%iterate(reorderData, levelCount, elapsedTime, bytesFetched, offsets, dest_sp, dest_dp)
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
  END SUBROUTINE asyncRestartCommData_collectData_sp

  SUBROUTINE asyncRestartCommData_destruct(this)
    CLASS(t_AsyncRestartCommData), INTENT(INOUT) :: this
    CALL this%inbufhandler%sync(.NOT.my_process_is_work(), .NOT.my_process_is_work(), .true.)
    CALL closeMpiWindow(this%win_buf, this%win)
    IF (ALLOCATED(this%inbufhandler)) THEN
      CALL this%inbufhandler%clear()
      DEALLOCATE(this%inbufhandler)
    ENDIF
    CALL this%cells%destruct()
    CALL this%verts%destruct()
    CALL this%edges%destruct()
  END SUBROUTINE asyncRestartCommData_destruct

  SUBROUTINE asyncRestartCommData_sync(this, i_are_baboon, after_update, goodbye)
    CLASS(t_AsyncRestartCommData), INTENT(INOUT) :: this
    LOGICAL,                       INTENT(IN)    :: i_are_baboon, after_update
    LOGICAL,      OPTIONAL,        INTENT(IN)    :: goodbye
    LOGICAL                                      :: goodbye_l

    IF (PRESENT(goodbye)) THEN
      goodbye_l = goodbye
    ELSE
      goodbye_l = .false.
    ENDIF
    CALL this%inbufhandler%sync(i_are_baboon, after_update, goodbye_l)
  END SUBROUTINE asyncRestartCommData_sync

  SUBROUTINE asyncRestartCommData_inbuffer_sync(this, i_are_baboon, after_update, goodbye)
    CLASS(inbufhandler_t),         INTENT(INOUT) :: this
    LOGICAL,                       INTENT(IN)    :: i_are_baboon, after_update, goodbye
    INTEGER                                      :: err, assert
    LOGICAL                                      :: hello_l
    CHARACTER(LEN=*), PARAMETER                  :: procedure_name = &
     &                          modname//'::inbufferhandler_sync'
    CHARACTER(LEN=64)                            :: occasion
    WRITE(occasion,'(a,i2)') 'DOM', this%domain
    hello_l = .false.
    IF (goodbye) THEN
      occasion = TRIM(occasion)//'<last>'
    ENDIF
    IF (this%hello) THEN
      hello_l = .true.
      occasion = TRIM(occasion)//'<first>'
      this%hello = .false.
    ELSE
      hello_l = .false.
    ENDIF
    assert = MPI_MODE_NOPUT
    IF (i_are_baboon .AND. .NOT.goodbye) THEN
      occasion = TRIM(occasion)//'<server>'
      IF (after_update .AND. .NOT. hello_l) THEN
        occasion = TRIM(occasion)//'<ready_for_unexpose>'
        CALL MPI_Win_complete(this%win, err)
        CALL ar_checkmpi(err, procedure_name, ': MPI_Win_complete('//TRIM(occasion)//') failed')
      ELSE
        occasion = TRIM(occasion)//'<ready_for_expose>'
        CALL MPI_Win_start(this%client_group, assert, this%win, err)
        CALL ar_checkmpi(err, procedure_name, ': MPI_Win_start('//TRIM(occasion)//') failed')
        this%hello = .false.
      ENDIF
    ELSE IF (.NOT. i_are_baboon) THEN
      occasion = TRIM(occasion)//'<client>'
      IF (after_update) THEN
        occasion = TRIM(occasion)//'<window_expose>'
        IF (.NOT.goodbye) THEN
          CALL MPI_Win_post(this%server_group, assert, this%win, err)
          CALL ar_checkmpi(err, procedure_name, ': MPI_Win_post('//TRIM(occasion)//') failed')
          this%hello = .false.
        ENDIF
      ELSE
        occasion = TRIM(occasion)//'<window_unexpose>'
        IF (.NOT.hello_l) THEN
          CALL MPI_Win_wait(this%win, err)
          CALL ar_checkmpi(err, procedure_name, ': MPI_Win_wait('//TRIM(occasion)//') failed')
        ENDIF
      ENDIF
    ENDIF
    IF (goodbye) THEN
      CALL MPI_Group_free(this%client_group, err)
      CALL ar_checkmpi(err, procedure_name, ': MPI_Group_free failed')
      CALL MPI_Group_free(this%server_group, err)
      CALL ar_checkmpi(err, procedure_name, ': MPI_Group_free failed')
      this%handshaked = .false.
    ENDIF
  END SUBROUTINE asyncRestartCommData_inbuffer_sync

  SUBROUTINE asyncRestartCommData_inbuffer_handshake(this, i_are_baboon)
    CLASS(inbufhandler_t), INTENT(INOUT) :: this
    LOGICAL,               INTENT(IN)    :: i_are_baboon
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: flags
    INTEGER                              :: my_flags(2), err, wgrp, &
      &                                     wgrp_size, wgrp_rank
    CHARACTER(LEN=*), PARAMETER                  :: procedure_name =      &
     &                          modname//'::handshake'
    IF (.NOT.this%is_init) THEN
      CALL ar_checkmpi(-99, procedure_name, ": window+buffers not initialized!")
    ENDIF
    CALL MPI_Win_get_group(this%win, wgrp, err)
    CALL ar_checkmpi(err, procedure_name, ": MPI_Win_get_group failed")
    CALL MPI_Group_size(wgrp, wgrp_size, err)
    CALL ar_checkmpi(err, procedure_name, ": MPI_Comm_size failed")
    CALL MPI_Group_rank(wgrp, wgrp_rank, err)
    CALL ar_checkmpi(err, procedure_name, ": MPI_Comm_rank failed")
    my_flags(1) = MERGE(wgrp_rank, MPI_UNDEFINED, i_are_baboon)
    my_flags(2) = MERGE(wgrp_rank, MPI_UNDEFINED, my_process_is_work())
    ALLOCATE(flags(wgrp_size*2))
    CALL MPI_Allgather(my_flags,  2, MPI_INTEGER, flags,  2, MPI_INTEGER, &
     &                 p_comm_work_restart, err)
    CALL ar_checkmpi(err, procedure_name, ": MPI_Allgather failed")
    IF (0 .NE. COUNT((flags(1::2).NE.MPI_UNDEFINED)  .AND. &
                     (flags(2::2).NE.MPI_UNDEFINED)) ) THEN
      CALL ar_checkmpi(-99, procedure_name, ": roles not consistent")
    ENDIF
    this%server_group = create_tgrp(1)
    this%client_group = create_tgrp(2)
    CALL MPI_Group_free(wgrp, err)
    CALL ar_checkmpi(err, procedure_name, ": MPI_Group_free failed")
    this%handshaked = .true.
  CONTAINS

    INTEGER FUNCTION create_tgrp(switch) RESULT(tgrp)
      INTEGER, INTENT(IN)  :: switch
      INTEGER              :: i, ii, tgrp_size
      INTEGER, ALLOCATABLE :: ranks_list(:)

      tgrp_size = COUNT(flags(switch::2).NE.MPI_UNDEFINED)
      ALLOCATE(ranks_list(tgrp_size))
      i = 0
      DO ii = 1, wgrp_size
        IF ( flags(switch + (ii - 1) * 2) .NE. MPI_UNDEFINED) THEN
          i = i + 1
          ranks_list(i) = flags(switch + (ii - 1) * 2)
        ENDIF
      ENDDO
      CALL MPI_Group_incl(wgrp, tgrp_size, ranks_list, tgrp, err)
      CALL ar_checkmpi(err, procedure_name, ": MPI_Group_incl failed")
      DEALLOCATE(ranks_list)
    END FUNCTION create_tgrp
  END SUBROUTINE asyncRestartCommData_inbuffer_handshake

#endif
END MODULE mo_async_restart_comm_data
