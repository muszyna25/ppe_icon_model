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
    USE ISO_C_BINDING, ONLY: C_PTR, C_INTPTR_T, C_F_POINTER
    USE mo_async_restart_packer, ONLY: t_AsyncRestartPacker
    USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_CELL
    USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info
    USE mo_exception, ONLY: finish, message, em_warn
    USE mo_fortran_tools, ONLY: t_ptr_2d, ensureSize
    USE mo_kind, ONLY: dp, i8
    USE mo_model_domain, ONLY: p_patch
    USE mo_mpi, ONLY: p_real_dp, p_comm_work_restart, p_pe_work, num_work_procs, p_mpi_wtime, my_process_is_work
    USE mo_restart_var_data, ONLY: t_RestartVarData
    USE mo_util_string, ONLY: int2string
#ifndef NOMPI
    USE mpi, ONLY: MPI_ADDRESS_KIND, MPI_INFO_NULL, MPI_LOCK_SHARED, MPI_MODE_NOCHECK, MPI_WIN_NULL, MPI_LOCK_EXCLUSIVE, &
                 & MPI_SUCCESS
#endif

    IMPLICIT NONE

    PUBLIC :: t_AsyncRestartCommData

    PRIVATE

    ! this combines all the DATA that's relevant for transfering the payload DATA from the worker PEs to the restart PEs
    TYPE t_AsyncRestartCommData
        PRIVATE

        ! DATA for remote memory access
        INTEGER :: mpiWindow
        REAL(dp), POINTER :: windowPtr(:)

        ! reorder data
        TYPE(t_AsyncRestartPacker), PUBLIC :: cells
        TYPE(t_AsyncRestartPacker), PUBLIC :: edges
        TYPE(t_AsyncRestartPacker), PUBLIC :: verts
    CONTAINS
! There is no point in pretending this is a usable class if NOMPI is defined.
#ifndef NOMPI
        PROCEDURE :: construct => asyncRestartCommData_construct
        PROCEDURE :: maxLevelSize => asyncRestartCommData_maxLevelSize  ! called to get the required buffer SIZE on the restart processes
        PROCEDURE :: getPacker => asyncRestartCommData_getPacker    ! RETURN the relevant t_AsyncRestartPacker object
        PROCEDURE :: postData => asyncRestartCommData_postData  ! called by the compute processes to WRITE their DATA to their memory window
        PROCEDURE :: collectData => asyncRestartCommData_collectData    ! called by the restart processes to fetch the DATA from the compute processes
        PROCEDURE :: destruct => asyncRestartCommData_destruct
#endif
    END TYPE t_AsyncRestartCommData

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_async_restart_comm_data"

CONTAINS

#ifndef NOMPI

    ! Opens an MPI memory window for the given amount of double values, returning both the ALLOCATED buffer AND the MPI window handle.
    SUBROUTINE openMpiWindow(doubleCount, communicator, mem_ptr_dp, mpi_win)
        INTEGER(KIND=MPI_ADDRESS_KIND), VALUE :: doubleCount   ! the requested memory window SIZE IN doubles
        INTEGER, VALUE :: communicator
        REAL(dp), POINTER, INTENT(OUT) :: mem_ptr_dp(:) ! this returns a fortran POINTER to the memory buffer
        INTEGER, INTENT(OUT) :: mpi_win ! this returns the handle to the MPI window

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

        mpi_win = MPI_WIN_NULL

        ! doubleCount is calculated as number of variables above, get number of bytes
        ! get the amount of bytes per REAL*8 variable (as used in MPI communication)
        CALL MPI_Type_extent(p_real_dp, nbytes_real, mpi_error)
        IF(mpi_error /= MPI_SUCCESS) CALL finish(routine, 'MPI_Type_extent returned error '//TRIM(int2string(mpi_error)))

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
        IF(mpi_error /= MPI_SUCCESS) CALL finish(routine, 'MPI_Alloc_mem returned error '//TRIM(int2string(mpi_error)))

        NULLIFY(mem_ptr_dp)
        CALL C_F_POINTER(c_mem_ptr, mem_ptr_dp, [doubleCount] )

#ifdef __xlC__
        ! IBM specific RMA hint, that we don't want window caching
        rma_cache_hint = MPI_INFO_NULL
        CALL MPI_Info_create(rma_cache_hint, mpi_error);
        IF(mpi_error /= MPI_SUCCESS) CALL finish(routine, 'MPI_Info_create returned error '//TRIM(int2string(mpi_error)))
        CALL MPI_Info_set(rma_cache_hint, "IBM_win_cache", "0", mpi_error)
        IF(mpi_error /= MPI_SUCCESS) CALL finish(routine, 'MPI_Info_set returned error '//TRIM(int2string(mpi_error)))
#endif

        ! create memory window for communication
        mem_ptr_dp(:) = 0._dp
        CALL MPI_Win_create(mem_ptr_dp, mem_bytes, nbytes_real, MPI_INFO_NULL, communicator, mpi_win, mpi_error)
        IF(mpi_error /= MPI_SUCCESS) CALL finish(routine, 'MPI_Win_create returned error '//TRIM(int2string(mpi_error)))

#ifdef __xlC__
        CALL MPI_Info_free(rma_cache_hint, mpi_error);
        IF(mpi_error /= MPI_SUCCESS) CALL finish(routine, 'MPI_Info_free returned error '//TRIM(int2string(mpi_error)))
#endif

    END SUBROUTINE openMpiWindow

    SUBROUTINE closeMpiWindow(windowPtr, mpiWindow)
        REAL(dp), POINTER, INTENT(INOUT) :: windowPtr(:)
        INTEGER, INTENT(INOUT) :: mpiWindow

        INTEGER :: mpi_error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":closeMpiWindow"

#ifdef NOMPI
        CALL finish(routine, "assertion failed: "//routine//"() must NOT be called without MPI")
#else
        ! release RMA window
        CALL MPI_Win_fence(0, mpiWindow, mpi_error)
        CALL warnMpiError('MPI_Win_fence')
        CALL MPI_Win_free(mpiWindow, mpi_error)
        CALL warnMpiError('MPI_Win_free')

        ! release RMA memory
        CALL MPI_Free_mem(windowPtr, mpi_error)
        CALL warnMpiError('MPI_Free_mem')

        ! safety
        mpiWindow = MPI_WIN_NULL
        windowPtr => NULL()

    CONTAINS
        SUBROUTINE warnMpiError(mpiCall)
            CHARACTER(LEN = *), INTENT(IN) :: mpiCall
            IF(mpi_error /= MPI_SUCCESS) CALL message(routine, mpiCall//"() returned error "//TRIM(int2string(mpi_error)), &
                                                     &level = em_warn, all_print = .TRUE.)
        END SUBROUTINE warnMpiError
#endif
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
            CALL me%cells%construct(p_patch(jg)%n_patch_cells_g, p_patch(jg)%n_patch_cells, p_patch(jg)%cells%decomp_info)
            CALL me%edges%construct(p_patch(jg)%n_patch_edges_g, p_patch(jg)%n_patch_edges, p_patch(jg)%edges%decomp_info)
            CALL me%verts%construct(p_patch(jg)%n_patch_verts_g, p_patch(jg)%n_patch_verts, p_patch(jg)%verts%decomp_info)
        ELSE
            !must not access p_patch on the restart processes, nevertheless the restart processes need to take part in the collective calls
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

            ! actually open the memory window
            CALL openMpiWindow(memWindowSize, p_comm_work_restart, me%windowPtr, me%mpiWindow)
        END IF
    END SUBROUTINE asyncRestartCommData_construct

    INTEGER FUNCTION asyncRestartCommData_maxLevelSize(me) RESULT(RESULT)
        CLASS(t_AsyncRestartCommData), INTENT(IN) :: me

        RESULT = MAX(me%cells%n_glb, me%edges%n_glb, me%verts%n_glb)
    END FUNCTION asyncRestartCommData_maxLevelSize

    ! Returns the pointer of the reorder data for the given field.
    FUNCTION asyncRestartCommData_getPacker(me, hgridType, routine) RESULT(RESULT)
        CLASS(t_AsyncRestartCommData), TARGET, INTENT(IN) :: me
        INTEGER, VALUE :: hgridType
        CHARACTER(LEN=*), INTENT(IN) :: routine

        TYPE(t_AsyncRestartPacker), POINTER :: RESULT

        RESULT => NULL()

        SELECT CASE(hgridType)
            CASE(GRID_UNSTRUCTURED_CELL)
                RESULT => me%cells
            CASE(GRID_UNSTRUCTURED_EDGE)
                RESULT => me%edges
            CASE(GRID_UNSTRUCTURED_VERT)
                RESULT => me%verts
            CASE default
                CALL finish(routine, "assertion failed: unexpected grid type")
        END SELECT

    END FUNCTION asyncRestartCommData_getPacker

    SUBROUTINE asyncRestartCommData_postData(me, hgridType, dataPointers, offset)
        CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
        INTEGER, VALUE :: hgridType
        TYPE(t_ptr_2d), ALLOCATABLE, INTENT(IN) :: dataPointers(:)
        !TODO[NH]: Replace this by an index within the t_AsyncRestartCommData structure.
        INTEGER(i8), INTENT(INOUT) :: offset

        TYPE(t_AsyncRestartPacker), POINTER :: reorderData
        INTEGER :: mpi_error, i
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":asyncRestartCommData_postData"

        ! get pointer to reorder data
        reorderData => me%getPacker(hgridType, routine)

        ! lock own window before writing to it
        CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, p_pe_work, MPI_MODE_NOCHECK, me%mpiWindow, mpi_error)
        IF(mpi_error /= MPI_SUCCESS) CALL finish(routine, 'MPI_Win_lock returned error '//TRIM(int2string(mpi_error)))

        ! just copy the OWN data points to the memory window
        DO i = 1, SIZE(dataPointers)
            CALL reorderData%packLevel(dataPointers(i)%p, me%windowPtr, offset)
        END DO

        ! unlock RMA window
        CALL MPI_Win_unlock(p_pe_work, me%mpiWindow, mpi_error)
        IF(mpi_error /= MPI_SUCCESS) CALL finish(routine, 'MPI_Win_unlock returned error '//TRIM(int2string(mpi_error)))
    END SUBROUTINE asyncRestartCommData_postData

    SUBROUTINE asyncRestartCommData_collectData(me, hgridType, levelCount, dest, offsets, elapsedTime, bytesFetched)
        CLASS(t_AsyncRestartCommData), TARGET, INTENT(INOUT) :: me
        INTEGER, VALUE :: hgridType, levelCount
        REAL(dp), INTENT(INOUT) :: dest(:,:)
        !TODO[NH]: Replace this by an index within the t_AsyncRestartCommData structure.
        INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(INOUT) :: offsets(0:num_work_procs-1)   ! This remembers the current offsets within all the different processes RMA windows. Pass a zero initialized array to the first collectData() CALL, THEN keep passing the array without external modifications.
        REAL(dp), INTENT(INOUT) :: elapsedTime
        INTEGER(i8), INTENT(INOUT) :: bytesFetched

        INTEGER :: sourceProc, doubleCount, reorderOffset, curLevel, mpi_error
        TYPE(t_AsyncRestartPacker), POINTER :: reorderData
        REAL(dp), POINTER, SAVE :: buffer(:) => NULL()
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":asyncRestartCommData_collectData"

        reorderData => me%getPacker(hgridType, routine)

        ! retrieve part of variable from every worker PE using MPI_Get
        DO sourceProc = 0, num_work_procs-1
            IF(reorderData%pe_own(sourceProc) == 0) CYCLE
            CALL ensureSize(buffer, reorderData%pe_own(sourceProc))

            ! number of words to transfer
            doubleCount = reorderData%pe_own(sourceProc) * levelCount
            elapsedTime = elapsedTime -  p_mpi_wtime()
            CALL MPI_Win_lock(MPI_LOCK_SHARED, sourceProc, MPI_MODE_NOCHECK, me%mpiWindow, mpi_error)
            IF(mpi_error /= MPI_SUCCESS) CALL finish(routine, 'MPI_Win_lock returned error '//TRIM(int2string(mpi_error)))

            CALL MPI_Get(buffer(1), doubleCount, p_real_dp, sourceProc, offsets(sourceProc), &
                                  & doubleCount, p_real_dp, me%mpiWindow, mpi_error)
            IF(mpi_error /= MPI_SUCCESS) CALL finish(routine, 'MPI_Get returned error '//TRIM(int2string(mpi_error)))

            CALL MPI_Win_unlock(sourceProc, me%mpiWindow, mpi_error)
            IF(mpi_error /= MPI_SUCCESS) CALL finish(routine, 'MPI_Win_unlock returned error '//TRIM(int2string(mpi_error)))

            elapsedTime  = elapsedTime + p_mpi_wtime()
            bytesFetched = bytesFetched + INT(doubleCount, i8)*8

            ! update the offset in the RMA window on compute PEs
            offsets(sourceProc) = offsets(sourceProc) + INT(doubleCount, i8)

            ! separate the levels received from PE "sourceProc":
            reorderOffset = 1
            DO curLevel = 1, levelCount
                CALL reorderData%unpackLevelFromPe(curLevel, sourceProc, buffer(reorderOffset:), dest)
                reorderOffset = reorderOffset + reorderData%pe_own(sourceProc)
            END DO
        END DO
    END SUBROUTINE asyncRestartCommData_collectData

    SUBROUTINE asyncRestartCommData_destruct(me)
        CLASS(t_AsyncRestartCommData), INTENT(INOUT) :: me

        CALL closeMpiWindow(me%windowPtr, me%mpiWindow)

        CALL me%cells%destruct()
        CALL me%verts%destruct()
        CALL me%edges%destruct()
    END SUBROUTINE asyncRestartCommData_destruct

#endif

END MODULE mo_async_restart_comm_data
