!>
!! A CLASS that handles the serialization of DATA into the async restart RMA windows, AND its subsequent deserialization into a global array.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_async_restart_packer
    USE mo_communication, ONLY: idx_no, blk_no
    USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info
    USE mo_exception, ONLY: finish
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_kind, ONLY: wp, dp, i8
    USE mo_mpi, ONLY: my_process_is_work, my_process_is_restart, p_n_work, p_int, p_comm_work, p_comm_work_2_restart, p_pe_work, &
                    & num_work_procs, p_bcast
#ifndef NOMPI
    USE mpi, ONLY: MPI_SUCCESS, MPI_ROOT, MPI_PROC_NULL
#endif

    IMPLICIT NONE

    PUBLIC :: t_asyncRestartPacker, restartBcastRoot

    PRIVATE

    !------------------------------------------------------------------------------------------------
    ! TYPE t_asyncRestartPacker describes how local cells/edges/verts
    ! have to be reordered to get the global array.
    ! Below, "points" refers to either cells, edges or verts.
    TYPE t_asyncRestartPacker
        PRIVATE
        INTEGER, PUBLIC :: n_glb  ! Global number of points per logical patch
        INTEGER, PUBLIC :: n_own  ! Number of own points (without halo, only belonging to logical patch). Only set on compute PEs, set to 0 on restart PEs.
        INTEGER, PUBLIC, ALLOCATABLE :: pe_own(:)   ! n_own, gathered for all compute PEs (set on all PEs)

        INTEGER, ALLOCATABLE :: own_idx(:), own_blk(:)  ! idx and blk for own points, only set on compute PEs

        ! Index how to reorder the contributions of all compute PEs into the global array (set on all PEs)
        INTEGER, ALLOCATABLE :: inverse_reorder_index(:)    ! given a gathered array index, this gives the index within the global array
    CONTAINS
! There is no point in pretending this is a usable class when NOMPI is defined.
#ifndef NOMPI
        PROCEDURE :: construct => asyncRestartPacker_construct    ! Collective across work AND restart. This IS a two step process utilizing the two methods below.

        PROCEDURE :: packLevel => asyncRestartPacker_packLevel    ! pack the contents of a single level/variable into our memory window
        PROCEDURE :: unpackLevelFromPe => asyncRestartPacker_unpackLevelFromPe  ! sort the DATA of a single level from a single PE into a global array.

        PROCEDURE :: destruct => asyncRestartPacker_destruct

        PROCEDURE, PRIVATE :: constructCompute => asyncRestartPacker_constructCompute  ! The part of the constructor that runs ONLY on the compute processes.
        PROCEDURE, PRIVATE :: transferToRestart => asyncRestartPacker_transferToRestart    ! Constructs the reorder DATA on the restart processes.
#endif
    END TYPE t_asyncRestartPacker

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_async_restart_packer"

CONTAINS

    ! Broadcast root for intercommunicator broadcasts from compute PEs to restart PEs using p_comm_work_2_restart.
    INTEGER FUNCTION restartBcastRoot() RESULT(RESULT)
#ifdef NOMPI
        RESULT = 0
#else
        IF(my_process_is_restart()) THEN
            ! root is proc 0 on the compute PEs
            RESULT = 0
        ELSE
            ! Special root setting for intercommunicators:
            ! The PE really sending must use MPI_ROOT, the others MPI_PROC_NULL.
            IF(p_pe_work == 0) THEN
                RESULT = MPI_ROOT
            ELSE
                RESULT = MPI_PROC_NULL
            END IF
        END IF
#endif
    END FUNCTION restartBcastRoot

#ifndef NOMPI
    ! the part of the constructor that runs on the compute processes
    SUBROUTINE asyncRestartPacker_constructCompute(me, n_points_g, n_points, decomp_info)
        CLASS(t_asyncRestartPacker), INTENT(INOUT) :: me ! Result: reorder data
        INTEGER, INTENT(IN) :: n_points_g      ! Global number of cells/edges/verts in logical patch
        INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
        TYPE(t_grid_domain_decomp_info), INTENT(IN) :: decomp_info

        INTEGER :: i, n, error
        LOGICAL, ALLOCATABLE :: owner_mask_1d(:) ! non-blocked owner mask for (logical) patch
        INTEGER, ALLOCATABLE :: glbidx_own(:), glbidx_glb(:), reorder_index_log_dom(:), pe_off(:), reorder_index(:)

        CHARACTER(LEN=*), PARAMETER :: routine = modname//':asyncRestartPacker_constructCompute'

#ifdef DEBUG
        WRITE(nerr, FORMAT_VALS3)routine, ' is called for p_pe=', p_pe
#endif

        ! just for safety
        IF(my_process_is_restart()) CALL finish(routine, "must not be called by restart processes")

        ! set the non-blocked patch owner mask
        ALLOCATE(owner_mask_1d(n_points), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "allocation of owner_mask_1d failed")
        DO i = 1, n_points
            owner_mask_1d(i) = decomp_info%owner_mask(idx_no(i), blk_no(i))
        ENDDO

        ! get number of owned cells/edges/verts (without halos)
        me%n_own = COUNT(owner_mask_1d(:))

        ! set index arrays to own cells/edges/verts
        ALLOCATE(me%own_idx(me%n_own), STAT=error)
        IF(error /= SUCCESS) CALL finish(routine, "allocation of me%own_idx failed")
        ALLOCATE(me%own_blk(me%n_own), STAT=error)
        IF(error /= SUCCESS) CALL finish(routine, "allocation of me%own_blk failed")

        ! global index of my own points
        ALLOCATE(glbidx_own(me%n_own), STAT=error)
        IF(error /= SUCCESS) CALL finish(routine, "allocation of glbidx_own failed")

        n = 0
        DO i = 1, n_points
            IF(owner_mask_1d(i)) THEN
                n = n+1
                me%own_idx(n) = idx_no(i)
                me%own_blk(n) = blk_no(i)
                glbidx_own(n)  = decomp_info%glb_index(i)
            ENDIF
        ENDDO
        DEALLOCATE(owner_mask_1d)

        ! gather the number of own points for every PE into me%pe_own
        ALLOCATE(me%pe_own(0:p_n_work-1), STAT=error)
        IF(error /= SUCCESS) CALL finish(routine, "allocation of me%pe_own failed")
        ALLOCATE(pe_off(0:p_n_work-1), STAT=error)
        IF(error /= SUCCESS) CALL finish(routine, "allocation of pe_off failed")

        CALL MPI_Allgather(me%n_own,  1, p_int, &
                          &me%pe_own, 1, p_int, p_comm_work, error)
        IF(error /= MPI_SUCCESS) CALL finish(routine, "MPI_Allgather returned an error")

        ! get offset within result array
        pe_off(0) = 0
        DO i = 1, p_n_work-1
            pe_off(i) = pe_off(i-1) + me%pe_own(i-1)
        ENDDO

        ! get global number of points for current patch
        me%n_glb = SUM(me%pe_own(:))

        ! Get the global index numbers of the data when it is gathered on PE 0
        ! exactly in the same order as it is retrieved later during restart.
        ALLOCATE(glbidx_glb(me%n_glb), STAT=error)
        IF(error /= SUCCESS) CALL finish(routine, "allocation of glbidx_glb failed")

        CALL MPI_Allgatherv(glbidx_own, me%n_own, p_int, &
                           &glbidx_glb, me%pe_own, pe_off, p_int, p_comm_work, error)
        IF(error /= MPI_SUCCESS) CALL finish(routine, "MPI_Allgatherv returned an error")
        DEALLOCATE(glbidx_own)
        DEALLOCATE(pe_off)

        ! spans the complete logical domain
        ALLOCATE(reorder_index_log_dom(n_points_g), STAT=error)
        IF(error /= SUCCESS) CALL finish(routine, "allocation of reorder_index_log_dom failed")
        reorder_index_log_dom(:) = 0

        DO i = 1, me%n_glb
            ! reorder_index_log_dom stores where a global point in logical domain comes from.
            ! It is nonzero only at the patch locations
            reorder_index_log_dom(glbidx_glb(i)) = i
        ENDDO
        DEALLOCATE(glbidx_glb)

        ! remove the zero entries from reorder_index_log_dom to form the reorder_index
        ALLOCATE(reorder_index(me%n_glb), STAT=error)
        IF(error /= SUCCESS) CALL finish(routine, "allocation of reorder_index failed")

        n = 0
        DO i = 1, n_points_g
            IF(reorder_index_log_dom(i)>0) THEN
                n = n+1
                reorder_index(n) = reorder_index_log_dom(i)
            ENDIF
        ENDDO
        DEALLOCATE(reorder_index_log_dom)
        IF(n /= me%n_glb) CALL finish(routine, "assertion failed: reorder_index is not surjective")

        ! compute the inverse of the reorder_index
        ALLOCATE(me%inverse_reorder_index(me%n_glb), STAT=error)
        IF(error /= SUCCESS) CALL finish(routine, "allocation of me%inverse_reorder_index failed")

        me%inverse_reorder_index(:) = 0
        DO i = 1, me%n_glb
            IF(me%inverse_reorder_index(reorder_index(i)) /= 0) THEN
                CALL finish(routine, "assertion failed: reorder_index is not injective")
            END IF
            me%inverse_reorder_index(reorder_index(i)) = i
        END DO
        DEALLOCATE(reorder_index)
    END SUBROUTINE asyncRestartPacker_constructCompute

    ! the part of the constructor that transfers the DATA from the compute processes to the restart processes
    SUBROUTINE asyncRestartPacker_transferToRestart(me)
        CLASS(t_asyncRestartPacker), INTENT(INOUT) :: me
        INTEGER :: error

        CHARACTER(LEN=*), PARAMETER :: routine = modname//':asyncRestartPacker_transferToRestart'

        ! transfer the global number of points, this is not yet known on restart PEs
        CALL p_bcast(me%n_glb,  restartBcastRoot(), p_comm_work_2_restart)

        IF(my_process_is_restart()) THEN
            me%n_own = 0  ! on restart PEs: n_own = 0, own_idx and own_blk are not allocated

            ! pe_own must be allocated for num_work_procs, not for p_n_work
            ALLOCATE(me%pe_own(0:num_work_procs-1), STAT=error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

            ALLOCATE(me%inverse_reorder_index(me%n_glb), STAT=error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        ENDIF

        CALL p_bcast(me%pe_own, restartBcastRoot(), p_comm_work_2_restart)
        CALL p_bcast(me%inverse_reorder_index, restartBcastRoot(), p_comm_work_2_restart)
    END SUBROUTINE asyncRestartPacker_transferToRestart

    SUBROUTINE asyncRestartPacker_construct(me, n_points_g, n_points, decomp_info)
        CLASS(t_asyncRestartPacker), INTENT(INOUT) :: me
        INTEGER, INTENT(IN) :: n_points_g      ! Global number of cells/edges/verts in logical patch
        INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
        TYPE(t_grid_domain_decomp_info), INTENT(IN) :: decomp_info

        IF(my_process_is_work()) THEN
            CALL me%constructCompute(n_points_g, n_points, decomp_info)
        END IF
        CALL me%transferToRestart()
    END SUBROUTINE asyncRestartPacker_construct

    SUBROUTINE asyncRestartPacker_packLevel(me, source, dest, offset)
        CLASS(t_asyncRestartPacker), INTENT(IN) :: me
        REAL(wp), INTENT(IN) :: source(:,:)
        REAL(dp), INTENT(INOUT) :: dest(:)
        !TODO[NH]: Replace this by an INTENT(IN) level count
        INTEGER(i8), INTENT(INOUT) :: offset

        INTEGER :: i

        DO i = 1, me%n_own
            offset = offset + 1
            dest(offset) = REAL(source(me%own_idx(i), me%own_blk(i)), dp)
        END DO
    END SUBROUTINE asyncRestartPacker_packLevel

    SUBROUTINE asyncRestartPacker_unpackLevelFromPe(me, level, pe, dataIn, globalArray)
        CLASS(t_asyncRestartPacker), INTENT(IN) :: me
        INTEGER, VALUE :: level, pe
        REAL(dp), INTENT(IN) :: dataIn(:)
        REAL(dp), INTENT(INOUT) :: globalArray(:,:)

        INTEGER :: offset, i

        offset = SUM(me%pe_own(0:pe - 1))
        DO i = 1, me%pe_own(pe)
            globalArray(me%inverse_reorder_index(offset + i), level) = dataIn(i)
        END DO
    END SUBROUTINE asyncRestartPacker_unpackLevelFromPe

    SUBROUTINE asyncRestartPacker_destruct(me)
        CLASS(t_asyncRestartPacker), INTENT(INOUT) :: me

        IF(ALLOCATED(me%own_idx)) DEALLOCATE(me%own_idx)
        IF(ALLOCATED(me%own_blk)) DEALLOCATE(me%own_blk)
        DEALLOCATE(me%pe_own)
        DEALLOCATE(me%inverse_reorder_index)
    END SUBROUTINE asyncRestartPacker_destruct

#endif

END MODULE mo_async_restart_packer
