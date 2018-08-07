!> Module for reading multifile restart files
!!
!! Note: The single file implementation of the restart input can be
!!       found in the module "mo_load_singlefile_restart"
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! ----------------------------------------------------------------------------------------------------
!!
!! All work processes are used to read the data from the
!! multifile. Each processor handles all domains of its horizontal
!! "chunk", since the domain loop is located outside of this routine,
!! in the calling "src/drivers" routine.
!!
!!
!! Implementation note [NH]: Unfortunately, t_comm_pattern uses ONLY
!!     p_comm_work AND does NOT allow for ANY other communicators.
!!     Extending t_comm_pattern to allow this seems prohibitive
!!     considering the small amount of time that's left to my work on
!!     this.  As such, we cannot make USE of the configured restart
!!     processes, as they are NOT IN the same p_comm_work IF they are
!!     dedicated restart processes.  Instead, all work processes are
!!     used to READ the DATA from the multifile.

MODULE mo_load_multifile_restart
    USE ISO_C_BINDING,             ONLY: C_CHAR, C_INT
    USE mo_broker_communication,   ONLY: t_BrokerCommunicationPattern
    USE mo_c_restart_util,         ONLY: checkMultifileDir
    USE mo_impl_constants,         ONLY: SUCCESS, VARNAME_LEN, SINGLE_T, REAL_T, INT_T
    USE mo_cdi_constants,          ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT
    USE mo_cdi,                    ONLY: streamOpenRead, streamInqVlist, vlistNvars, vlistCopyVarName, &
      &                                  streamClose, streamReadVar,                                   &
      &                                  streamReadVarSlice, streamReadVarSliceF, CDI_GLOBAL,          &
      &                                  cdiInqAttInt
    USE mo_communication_orig,     ONLY: t_comm_pattern_orig, exchange_data_noblk
    USE mo_decomposition_tools,    ONLY: t_glb2loc_index_lookup, init_glb2loc_index_lookup, set_inner_glb_index, &
      &                                  deallocate_glb2loc_index_lookup
    USE mo_dynamics_config,        ONLY: nnew, nnew_rcf
    USE mo_exception,              ONLY: finish
    USE mo_fortran_tools,          ONLY: t_ptr_2d_sp, t_ptr_2d, t_ptr_2d_int
    USE mo_kind,                   ONLY: sp, dp, i8
    USE mo_model_domain,           ONLY: t_patch
    USE mo_mpi,                    ONLY: p_barrier, p_comm_work, p_comm_size, p_comm_rank, my_process_is_work, &
      &                                  p_allreduce, p_sum_op, p_mpi_wtime, my_process_is_stdio, p_bcast,     &
      &                                  my_process_is_mpi_workroot
    USE mo_multifile_restart_util, ONLY: multifilePayloadPath, kVarName_globalCellIndex,      &
      &                                  kVarName_globalEdgeIndex, kVarName_globalVertIndex,  &
      &                                  t_ptr_1d_generic
    USE mo_parallel_config,        ONLY: nproma, idx_no, blk_no
    USE mo_restart_attributes,     ONLY: t_RestartAttributeList, getAttributesForRestarting
    USE mo_restart_var_data,       ONLY: t_restartVarData, getLevelPointers, has_valid_time_level
    USE mo_timer,                  ONLY: timer_start, timer_stop, timer_load_restart_io, &
      &                                  timer_load_restart_comm_setup, timer_load_restart_communication, &
      &                                  timer_load_restart_get_var_id, timers_level
    USE mo_util_cdi,               ONLY: get_cdi_varID
    USE mo_util_string,            ONLY: charArray_equal, real2string, int2string

    IMPLICIT NONE

    PUBLIC :: multifileCheckRestartFiles, multifileReadPatch

    TYPE t_PayloadFile
        INTEGER :: streamId, vlistId, varCount
        INTEGER :: cellCount, edgeCount, vertCount
        INTEGER :: cellIndexVar, edgeIndexVar, vertIndexVar
    CONTAINS
      ! open the file AND determine the cell/edge/vert counts
      PROCEDURE :: construct => payloadFile_construct

      ! READ one of the three variables holding the global
      ! cell/edge/vert indices of the points within this file
      PROCEDURE :: readIndexVar => payloadFile_readIndexVar

      ! READ a named variable into the supplied buffer
      PROCEDURE :: readVarLevel => payloadFile_readVarLevel

      ! close the file
      PROCEDURE :: destruct => payloadFile_destruct
    END TYPE t_PayloadFile

    ! This IS both a buffer AND a tool for communication:
    !
    ! First, the buffer IS filled locally, THEN the DATA IS
    ! redistributed among the processes AND written to the provided
    ! destination.
    TYPE t_ReadBuffer
      ! one based, must be equal to SIZE(readBuffer_1d_X)+1 when
      ! redistribute IS called, AND will be reset to 1 IN such a
      ! CALL.:
      INTEGER :: curWindowOffset

      ! should be pointers to tell Fortran, that these _will_ be
      ! aliased
      REAL(sp), POINTER :: readBuffer_1d_s(:)
      REAL(dp), POINTER :: readBuffer_1d_d(:)
      INTEGER,  POINTER :: readBuffer_1d_int(:)

      TYPE(t_comm_pattern_orig), POINTER :: commPattern
    CONTAINS
        PROCEDURE :: construct => readBuffer_construct

        ! This IS used to fill the buffer: It returns a POINTER to a
        ! window within the buffer for the caller to WRITE.  The
        ! caller IS expected to CALL nextReadWindow exactly as many
        ! times as necessary to completely fill the buffer before
        ! calling redistribute(), after which the CYCLE may start over
        ! again.
        PROCEDURE :: nextWindow => readBuffer_nextWindow

        ! This kicks of the communication that redistributes the DATA
        ! to the processes that actually require it.  It also has the
        ! effect of reseting the window so that nextWindow() can be
        ! called again.  This CALL also initializes the halo points of
        ! the given array.
        PROCEDURE :: redistribute_sp => readBuffer_redistribute_sp
        PROCEDURE :: redistribute_dp => readBuffer_redistribute_dp
        PROCEDURE :: redistribute_int => readBuffer_redistribute_int
        GENERIC :: redistribute => redistribute_sp, redistribute_dp, redistribute_int

        PROCEDURE :: destruct => readBuffer_destruct
    END TYPE t_ReadBuffer

    ! A MultifilePatchReader IS constructed for a single patch, AND IS
    ! responsible to READ all the ASSOCIATED patch<domain>_<N>.nc
    ! files, redistributing their contents to the correct processes.
    TYPE t_MultifilePatchReader
        INTEGER :: domain
        TYPE(t_ReadBuffer) :: readBufferCells, readBufferEdges, readBufferVerts
        TYPE(t_PayloadFile), ALLOCATABLE :: files(:)
    CONTAINS
        PROCEDURE :: construct => multifilePatchReader_construct
        PROCEDURE :: readData => multifilePatchReader_readData
        PROCEDURE :: destruct => multifilePatchReader_destruct
    END TYPE t_MultifilePatchReader

    CHARACTER(*), PARAMETER :: modname = "mo_load_multifile_restart"

CONTAINS

    SUBROUTINE payloadFile_construct(me, multifilePath, domain, partId)
        CLASS(t_PayloadFile), INTENT(INOUT) :: me
        CHARACTER(*), INTENT(IN) :: multifilePath
        !partId IS the process ID of the restart writer that that
        !wrote the payload file that IS to be READ
        INTEGER, VALUE :: domain, partId

        INTEGER :: varId, gridId, dummy, icount(1)
        LOGICAL :: haveCellCount, haveEdgeCount, haveVertCount
        CHARACTER(KIND = C_CHAR), POINTER :: nameCharArray(:)
        CHARACTER(*), PARAMETER :: routine = modname//":payloadFile_construct"

        IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)

        ! open the file AND cache the most important CDI IDs
        me%streamId = streamOpenRead(multifilePayloadPath(multifilePath, domain, partId))
        IF (me%streamId < 0) THEN
          WRITE (0,*) "failed to open ", multifilePayloadPath(multifilePath, domain, partId), " for reading"
        END IF
        me%vlistId = streamInqVlist(me%streamId)
        me%varCount = vlistNvars(me%vlistId)

        ! search for the global_*_indices variables AND determine their sizes
        haveCellCount = .FALSE.
        haveEdgeCount = .FALSE.
        haveVertCount = .FALSE.

        DO varId = 0, me%varCount - 1

            nameCharArray => vlistCopyVarName(me%vlistId, varId)

            IF(charArray_equal(nameCharArray, kVarName_globalCellIndex)) THEN
                IF(haveCellCount) THEN
                    CALL finish(routine, "corrupted restart file detected: double definition of "//kVarName_globalCellIndex)
                END IF

                ! Important Note: The cell count contained in the
                ! payload file may be 0!
                !
                ! The problem is as follows: If the restart file has
                ! been created in a run with processor splitting
                ! enabled, then it may be that a payload file is
                ! virtually "empty" (0 cells). However, in NetCDF it
                ! is not possible to define a zero dimension without
                ! declaring it to be an "unlimited dimension". This,
                ! again, is interpreted by the CDI as the time
                ! variable!
                !
                ! In order to overcome this problem without too much
                ! hassle, we define size-0 cell counts as size-1 grids
                ! with an additional integer attribute "cellCount" ...

                ! gridId = cdiInqVarGrid(me%vlistId, varId)
                ! me%cellCount = gridInqSize(gridId)
                
                dummy =  cdiInqAttInt(me%vlistId, CDI_GLOBAL, "cellCount", 1, icount);
                me%cellCount = icount(1)
                me%cellIndexVar = varId
                haveCellCount = .TRUE.
            ELSE IF(charArray_equal(nameCharArray, kVarName_globalEdgeIndex)) THEN
                IF(haveEdgeCount) THEN
                    CALL finish(routine, "corrupted restart file detected: double definition of "//kVarName_globalEdgeIndex)
                END IF

                ! gridId = cdiInqVarGrid(me%vlistId, varId) 
                ! me%edgeCount = gridInqSize(gridId)
                dummy =  cdiInqAttInt(me%vlistId, CDI_GLOBAL, "edgeCount", 1, icount);
                me%edgeCount = icount(1)

                me%edgeIndexVar = varId
                haveEdgeCount = .TRUE.
            ELSE IF(charArray_equal(nameCharArray, kVarName_globalVertIndex)) THEN
                IF(haveVertCount) THEN
                    CALL finish(routine, "corrupted restart file detected: double definition of "//kVarName_globalVertIndex)
                END IF

                ! gridId = cdiInqVarGrid(me%vlistId, varId)
                ! me%vertCount = gridInqSize(gridId)
                dummy =  cdiInqAttInt(me%vlistId, CDI_GLOBAL, "vertCount", 1, icount);
                me%vertCount = icount(1)

                me%vertIndexVar = varId
                haveVertCount = .TRUE.
            END IF
            DEALLOCATE(nameCharArray)

          END DO

        IF(timers_level >= 7) CALL timer_stop(timer_load_restart_io)

        ! sanity check
        IF(.NOT.haveCellCount) THEN
          CALL finish(routine, "corrupted restart file detected: "//kVarName_globalCellIndex//" IS missing")
        END IF
        IF(.NOT.haveEdgeCount) THEN
          CALL finish(routine, "corrupted restart file detected: "//kVarName_globalEdgeIndex//" IS missing")
        END IF
        IF(.NOT.haveVertCount) THEN
          CALL finish(routine, "corrupted restart file detected: "//kVarName_globalVertIndex//" IS missing")
        END IF
    END SUBROUTINE payloadFile_construct

    ! gridType IS either GRID_UNSTRUCTURED_CELL,
    ! GRID_UNSTRUCTURED_EDGE, OR GRID_UNSTRUCTURED_VERT the SIZE of
    ! globalIndices must agree with the respective me%XXXXCount
    ! variable
    SUBROUTINE payloadFile_readIndexVar(me, gridType, globalIndices)
        CLASS(t_PayloadFile), INTENT(IN) :: me
        INTEGER, VALUE :: gridType
        INTEGER, INTENT(INOUT) :: globalIndices(:)

        INTEGER :: pointCount, varId, error, i
        INTEGER(C_INT) :: trash
        REAL(dp), ALLOCATABLE :: buffer(:)
        CHARACTER(*), PARAMETER :: routine = modname//":payloadFile_readIndexVar"

        !get the parameters for the operation
        SELECT CASE(gridType)
            CASE(GRID_UNSTRUCTURED_CELL)
                pointCount = me%cellCount
                varId = me%cellIndexVar
            CASE(GRID_UNSTRUCTURED_EDGE)
                pointCount = me%edgeCount
                varId = me%edgeIndexVar
            CASE(GRID_UNSTRUCTURED_VERT)
                pointCount = me%vertCount
                varId = me%vertIndexVar
            CASE DEFAULT
                CALL finish(routine, "assertion failed: illegal value of gridType argument")
        END SELECT

        !check the SIZE of the RESULT buffer
        IF(SIZE(globalIndices) /= pointCount) THEN
          CALL finish(routine, "assertion failed: globalIndices(:) has the wrong size")
        END IF

        IF (pointCount > 0) THEN

          !ALLOCATE a buffer to READ the REAL values from the file
          ALLOCATE(buffer(pointCount), STAT = error)
          IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
          
          !READ the indices into the REAL buffer
          IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
          CALL streamReadVar(me%streamId, varId, buffer, trash)
          IF(timers_level >= 7) CALL timer_stop(timer_load_restart_io)
          
          !convert to INTEGER
          DO i = 1, pointCount
            globalIndices(i) = INT(buffer(i))
          END DO

          !cleanup
          DEALLOCATE(buffer)

        END IF
    END SUBROUTINE payloadFile_readIndexVar

    ! The offset IS 1 based AND will be updated to point to the next
    ! region within the buffer that should be filled.
    SUBROUTINE payloadFile_readVarLevel(me, varData, level, buffer)
        CLASS(t_PayloadFile), INTENT(INOUT) :: me
        TYPE(t_restartVarData), INTENT(IN) :: varData
        INTEGER, VALUE :: level
        TYPE(t_ReadBuffer), INTENT(INOUT) :: buffer

        INTEGER :: varId, pointCount, trash
        TYPE(t_ptr_1d_generic) :: curWindow
        CHARACTER(*), PARAMETER :: routine = modname//":payloadFile_readVarLevel"
        REAL(dp), ALLOCATABLE :: buffer_dp(:)

        ! get the correct point count for this variable
        SELECT CASE(varData%info%hgrid)
            CASE(GRID_UNSTRUCTURED_CELL)
                pointCount = me%cellCount
            CASE(GRID_UNSTRUCTURED_EDGE)
                pointCount = me%edgeCount
            CASE(GRID_UNSTRUCTURED_VERT)
                pointCount = me%vertCount
            CASE DEFAULT
                CALL finish(routine, "assertion failed: illegal value of vardata%info%hgrid")
        END SELECT

        IF (pointCount > 0) THEN

          ! get the cdi var ID (unfortunately, CDI itself does NOT provide NAME based lookup
          IF(timers_level >= 7) CALL timer_start(timer_load_restart_get_var_id)
          varId = get_cdi_varID(me%streamId, TRIM(varData%info%NAME))
          IF(timers_level >= 7) CALL timer_stop(timer_load_restart_get_var_id)

          ! get the buffer window
          curWindow = buffer%nextWindow(pointCount)

          ! actually READ the DATA
          IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)

          SELECT CASE(varData%info%data_type)
          CASE(REAL_T)
            IF(.NOT.ASSOCIATED(curWindow%dp)) CALL finish(routine, "assertion failed: buffer not associated")
            IF(pointCount /= SIZE(curWindow%dp)) THEN
              CALL finish(routine, "assertion failed: supplied buffer has wrong size")
            END IF
            CALL streamReadVarSlice(me%streamId, varId, level, curWindow%dp, trash)
            !
          CASE(SINGLE_T)
            IF(.NOT.ASSOCIATED(curWindow%sp)) CALL finish(routine, "assertion failed: buffer not associated")
            IF(pointCount /= SIZE(curWindow%sp)) THEN
              CALL finish(routine, "assertion failed: supplied buffer has wrong size")
            END IF
            CALL streamReadVarSliceF(me%streamId, varId, level, curWindow%sp, trash)
            !
          CASE(INT_T)
            IF(.NOT.ASSOCIATED(curWindow%int)) CALL finish(routine, "assertion failed: buffer not associated")
            IF(pointCount /= SIZE(curWindow%int)) THEN
              CALL finish(routine, "assertion failed: supplied buffer has wrong size")
            END IF
            ALLOCATE(buffer_dp(pointCount))
            CALL streamReadVarSlice(me%streamId, varId, level, buffer_dp, trash)
            curWindow%int(:) = INT(buffer_dp(:))
            DEALLOCATE(buffer_dp)
            !
          CASE DEFAULT
            CALL finish(routine, "Internal error! Variable "//TRIM(varData%info%name))
          END SELECT
          IF(timers_level >= 7) CALL timer_stop(timer_load_restart_io)

        END IF
    END SUBROUTINE payloadFile_readVarLevel

    SUBROUTINE payloadFile_destruct(me)
        CLASS(t_PayloadFile), INTENT(INOUT) :: me

        IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
        CALL streamClose(me%streamId)
        IF(timers_level >= 7) CALL timer_stop(timer_load_restart_io)
    END SUBROUTINE payloadFile_destruct

    SUBROUTINE readBuffer_construct(me, name, providedGlobalIndices, requiredGlobalIndices)
        CLASS(t_ReadBuffer), INTENT(INOUT) :: me
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER, INTENT(IN) :: providedGlobalIndices(:), requiredGlobalIndices(:)

        INTEGER :: pointCount, error
        CHARACTER(*), PARAMETER :: routine = modname//":readBuffer_construct"

        ! Build the communication pattern.
        ! Note: These communication patterns also set the halo points.
        me%commPattern => makeRedistributionPattern(name, providedGlobalIndices, requiredGlobalIndices)

        ! Allocate the READ buffers.
        pointCount = SIZE(providedGlobalIndices)
        ALLOCATE(me%readBuffer_1d_d(pointCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(me%readBuffer_1d_s(pointCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(me%readBuffer_1d_int(pointCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        ! Init the window offset
        me%curWindowOffset = 1
    END SUBROUTINE readBuffer_construct

    FUNCTION readBuffer_nextWindow(me, pointCount) RESULT(resultVar)
        CLASS(t_ReadBuffer), INTENT(INOUT) :: me
        INTEGER, VALUE :: pointCount
        TYPE(t_ptr_1d_generic) :: resultVar

        INTEGER :: windowEnd
        CHARACTER(*), PARAMETER :: routine = modname//":readBuffer_nextWindow"

        ! sanity check
        windowEnd = me%curWindowOffset + pointCount - 1
        IF(windowEnd > SIZE(me%readBuffer_1d_d)) THEN
          WRITE (0,*) "windowEnd = ", windowEnd, "; SIZE(me%readBuffer_1d_d) = ", SIZE(me%readBuffer_1d_d)
          CALL finish(routine, "assertion failed: attempt to overfill buffer")
        END IF

        ! set the pointers to grant access to the requested window
        resultVar%sp  => me%readBuffer_1d_s(me%curWindowOffset:windowEnd)
        resultVar%dp  => me%readBuffer_1d_d(me%curWindowOffset:windowEnd)
        resultVar%int => me%readBuffer_1d_int(me%curWindowOffset:windowEnd)

        ! remember that we have filled this part of the buffer
        me%curWindowOffset = windowEnd + 1
    END FUNCTION readBuffer_nextWindow

    SUBROUTINE readBuffer_redistribute_sp(me, destination)
        CLASS(t_ReadBuffer), INTENT(INOUT) :: me
        REAL(sp), INTENT(INOUT) :: destination(:,:)

        INTEGER :: i, pointCount
        CHARACTER(*), PARAMETER :: routine = modname//":readBuffer_redistribute_sp"

        ! sanity check
        pointCount = SIZE(me%readBuffer_1d_s)
        IF(me%curWindowOffset - 1 /= pointCount) THEN
            CALL finish(routine, "assertion failed: redistribute called on a non-full buffer")
        END IF

        ! redistribute the DATA
        IF(timers_level >= 7) CALL timer_start(timer_load_restart_communication)
        CALL exchange_data_noblk(me%commPattern, destination, me%readBuffer_1d_s)
        IF(timers_level >= 7) CALL timer_stop(timer_load_restart_communication)

        ! reset the window
        me%curWindowOffset = 1
    END SUBROUTINE readBuffer_redistribute_sp

    SUBROUTINE readBuffer_redistribute_dp(me, destination)
        CLASS(t_ReadBuffer), INTENT(INOUT) :: me
        REAL(dp), INTENT(INOUT) :: destination(:,:)

        INTEGER :: i, pointCount
        CHARACTER(*), PARAMETER :: routine = modname//":readBuffer_redistribute_dp"

        ! sanity check
        pointCount = SIZE(me%readBuffer_1d_d)
        IF(me%curWindowOffset - 1 /= pointCount) THEN
            CALL finish(routine, "assertion failed: redistribute called on a non-full buffer")
        END IF

        ! redistribute the DATA
        IF(timers_level >= 7) CALL timer_start(timer_load_restart_communication)
        CALL exchange_data_noblk(me%commPattern, destination, me%readBuffer_1d_d)
        IF(timers_level >= 7) CALL timer_stop(timer_load_restart_communication)

        ! reset the window
        me%curWindowOffset = 1
    END SUBROUTINE readBuffer_redistribute_dp

    SUBROUTINE readBuffer_redistribute_int(me, destination)
        CLASS(t_ReadBuffer), INTENT(INOUT) :: me
        INTEGER, INTENT(INOUT) :: destination(:,:)

        INTEGER :: i, pointCount
        CHARACTER(*), PARAMETER :: routine = modname//":readBuffer_redistribute_int"

        ! sanity check
        pointCount = SIZE(me%readBuffer_1d_int)
        IF(me%curWindowOffset - 1 /= pointCount) THEN
            CALL finish(routine, "assertion failed: redistribute called on a non-full buffer")
        END IF

        ! redistribute the DATA
        IF(timers_level >= 7) CALL timer_start(timer_load_restart_communication)
        CALL exchange_data_noblk(me%commPattern, destination, me%readBuffer_1d_int)
        IF(timers_level >= 7) CALL timer_stop(timer_load_restart_communication)

        ! reset the window
        me%curWindowOffset = 1
    END SUBROUTINE readBuffer_redistribute_int

    SUBROUTINE readBuffer_destruct(me)
        CLASS(t_ReadBuffer), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":readBuffer_destruct"

        IF(me%curWindowOffset /= 1) CALL finish(routine, "assertion failed: buffer NOT empty on destruction")

        DEALLOCATE(me%readBuffer_1d_d, me%readBuffer_1d_s, me%readBuffer_1d_int)
        CALL me%commPattern%delete()
        DEALLOCATE(me%commPattern)
    END SUBROUTINE readBuffer_destruct

    SUBROUTINE multifileCheckRestartFiles(filename)
        CHARACTER(*), INTENT(IN) :: filename

        TYPE(t_RestartAttributeList), POINTER :: restartAttributes
        INTEGER :: n_dom, multifile_file_count
        CHARACTER(*), PARAMETER :: routine = modname//":multifileCheckRestartFiles"

        !get the expected domain AND file counts from the restart attributes
        restartAttributes => getAttributesForRestarting()

        ! If not all domains are active, then the multifile restart
        ! directory contains less payload files than there are
        ! domains:
        n_dom = restartAttributes%getInteger('multifile_n_dom_active')

        multifile_file_count = restartAttributes%getInteger('multifile_file_count')

        ! Must NOT start the timer here because the timers have NOT been initialized at this point.
!       IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)

        !check the multifile
        IF(checkMultifileDir(filename, n_dom, multifile_file_count) /= SUCCESS) THEN
            CALL finish(routine, "'"//filename//"' is not a valid restart multifile")
        END IF
!       IF(timers_level >= 7) CALL timer_stop(timer_load_restart_io)
    END SUBROUTINE multifileCheckRestartFiles

    ! Opens all the payload files that should be opened by this
    ! process, creating t_PayloadFile objects for them, AND
    ! calculating the total numbers of cells/edges/vertices contained
    ! IN these files.
    !
    ! This FUNCTION performs the distribution of the payload files
    ! among the restart processes, subsequent code just has to iterate
    ! over the array returned by this FUNCTION.
    SUBROUTINE openPayloadFiles(multifilePath, multifile_file_count, files, domain, out_totalCellCount, &
      &                      out_totalEdgeCount, out_totalVertCount)
        TYPE(t_PayloadFile), ALLOCATABLE :: files(:)
        CHARACTER(*), INTENT(IN) :: multifilePath
        INTEGER, VALUE :: multifile_file_count, domain
        INTEGER, INTENT(OUT) :: out_totalCellCount, out_totalEdgeCount, out_totalVertCount

        INTEGER :: procCount, myRank, myFileCount, error, curFileIndex, curFileId, proc_file_ratio
        CHARACTER(*), PARAMETER :: routine = modname//":openPayloadFiles"

        IF(.NOT.my_process_is_work()) CALL finish(routine, "assertion failed: routine must ONLY be called on work processes")

        ! Determine how many files this process needs to READ.
        !
        ! Sorry for this cryptic formula; it's intended effect IS to
        ! compute exactly how many files each process gets when
        ! distributing them round-robin style.  The asserts IN/after
        ! the loop at the END of this FUNCTION ensure that this
        ! calculation agrees with the loop control.
        procCount = p_comm_size(p_comm_work)
        myRank = p_comm_rank(p_comm_work)
        proc_file_ratio = procCount/multifile_file_count
        IF (proc_file_ratio > 2) THEN
          IF (MOD(myRank,proc_file_ratio) == 0 .AND. myRank/proc_file_ratio < multifile_file_count) THEN
            myFileCount = 1
          ELSE
            myFileCount = 0
          ENDIF
        ELSE
          myFileCount = (multifile_file_count - myRank + procCount - 1)/procCount
        ENDIF
        ! Allocate space to open all our files IN parallel.
        ALLOCATE(files(myFileCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        ! Open our restart files.
        !
        ! This loop assigns the restart files to the restart processes
        ! IN round-robin fashion, AND its loop control must agree with
        ! myFileCount calculated above.  curFileId AND
        ! restartWorkProcId() are zero based, so we need to stop at
        ! multifile_file_count - 1.
        curFileIndex = 0    ! necessary for the CASE that myFileCount == 0
        out_totalCellCount = 0
        out_totalEdgeCount = 0
        out_totalVertCount = 0
        IF (proc_file_ratio <= 2) THEN
          DO curFileId = myRank, multifile_file_count - 1, procCount
            curFileIndex = curFileId/procCount + 1
            IF(curFileIndex > myFileCount) CALL finish(routine, "assertion failed: myFileCount IS too small")

            CALL files(curFileIndex)%construct(multifilePath, domain, curFileId)
            out_totalCellCount = out_totalCellCount + files(curFileIndex)%cellCount
            out_totalEdgeCount = out_totalEdgeCount + files(curFileIndex)%edgeCount
            out_totalVertCount = out_totalVertCount + files(curFileIndex)%vertCount
          END DO
        ELSE IF (myFileCount == 1) THEN
          curFileId = myRank/proc_file_ratio
          curFileIndex = 1

          CALL files(curFileIndex)%construct(multifilePath, domain, curFileId)
          out_totalCellCount = out_totalCellCount + files(curFileIndex)%cellCount
          out_totalEdgeCount = out_totalEdgeCount + files(curFileIndex)%edgeCount
          out_totalVertCount = out_totalVertCount + files(curFileIndex)%vertCount
        ENDIF
        IF(curFileIndex /= myFileCount) CALL finish(routine, "assertion failed: myFileCount IS too large")
    END SUBROUTINE openPayloadFiles

    ! globalSize: Total count of different global indices.
    ! providedGlobalIndices(N): Global indices of the points provided
    ! by this PE IN the order IN which they are provided. All entries
    ! must be IN the range [1, globalSize].  requiredGlobalIndices(M):
    ! Global indices of the points wanted by this PE IN the order IN
    ! which they are to be placed IN the resulting array. All entries
    ! must be IN the range [1, globalSize].  owners(M): Rank of the PE
    ! providing each of the required points.
    FUNCTION makeCommPattern(globalSize, providedGlobalIndices, requiredGlobalIndices, owners) RESULT(resultVar)
      TYPE(t_comm_pattern_orig), POINTER :: resultVar
        INTEGER, VALUE :: globalSize
        INTEGER, INTENT(IN) :: providedGlobalIndices(:), requiredGlobalIndices(:), owners(:)

        INTEGER :: i, nreq, nown
        INTEGER :: local_index(SIZE(providedGlobalIndices)), &
             local_owner(SIZE(providedGlobalIndices))
        TYPE(t_glb2loc_index_lookup) :: lookupTable
        CHARACTER(*), PARAMETER :: routine = modname//":makeCommPattern"

        nreq = SIZE(requiredGlobalIndices)
        IF(nreq /= SIZE(owners)) THEN
          CALL finish(routine, "assertion failed: requiredGlobalIndices(:) &
               &and owners(:) must have the same number of entries")
        END IF

        CALL init_glb2loc_index_lookup(lookupTable, globalSize)
        nown = SIZE(providedGlobalIndices)
        DO i = 1, nown
          local_index(i) = i
        END DO
        CALL set_inner_glb_index(lookupTable, providedGlobalIndices, local_index)
        local_owner = p_comm_rank(p_comm_work)

        ALLOCATE(resultvar)
        CALL resultVar%setup(nreq, owners, requiredGlobalIndices, &
             lookupTable, nown, local_owner, providedGlobalIndices)
        CALL deallocate_glb2loc_index_lookup(lookupTable)
    END FUNCTION makeCommPattern

    ! This asserts that the given pattern redistributes a given array
    ! exactly as described by the global indices arrays.
    SUBROUTINE checkRedistributionPattern(pattern, providedGlobalIndices, requiredGlobalIndices)
        TYPE(t_comm_pattern_orig), INTENT(IN) :: pattern
        INTEGER, INTENT(IN) :: providedGlobalIndices(:), requiredGlobalIndices(:)

        INTEGER :: inputSize, outputSize, i, error
        INTEGER, ALLOCATABLE :: input(:,:), output(:,:)
        CHARACTER(*), PARAMETER :: routine = modname//":checkRedistributionPattern"

        ! ALLOCATE some memory
        inputSize = SIZE(providedGlobalIndices)
        outputSize = SIZE(requiredGlobalIndices)
        ALLOCATE(input(nproma, blk_no(inputSize)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(output(nproma, blk_no(outputSize)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        ! perform an exchange according to the given communication pattern
        output(:,:) = -1
        CALL exchange_data_noblk(pattern, output, providedGlobalIndices)

        ! assert that the RESULT IS as expected
        DO i = 1, outputSize
            IF(output(idx_no(i), blk_no(i)) /= REAL(requiredGlobalIndices(i))) THEN
                CALL finish(routine, "assertion failed: redistribution pattern does not work as expected")
            END IF
        END DO

        ! cleanup
        DEALLOCATE(input, output)
    END SUBROUTINE checkRedistributionPattern

    ! Create a t_comm_pattern that IS able to redistribute the DATA as
    ! it IS READ from the multifile to the domain decomposition of the
    ! work processes.  The tricky part about this IS, that
    ! setup_comm_pattern() requires an owner(:) array as input, which
    ! we have to create first.
    !
    ! The resulting t_comm_pattern will also initialize the halo points.
    FUNCTION makeRedistributionPattern(name, providedGlobalIndices, requiredGlobalIndices) RESULT(resultVar)
        TYPE(t_comm_pattern_orig), POINTER :: resultVar
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER, INTENT(IN) :: providedGlobalIndices(:), requiredGlobalIndices(:)

        INTEGER :: globalSize, procCount, myRank, error, i
        INTEGER, ALLOCATABLE :: brokerBounds(:), providerBuffer(:), brokerBuffer(:), providerOfRequiredPoint(:)
        TYPE(t_BrokerCommunicationPattern) :: providerToBroker, consumerToBroker
        CHARACTER(*), PARAMETER :: routine = modname//":makeRedistributionPattern"

        IF(timers_level >= 7) CALL timer_start(timer_load_restart_comm_setup)

        ! compute the global SIZE of the field
        globalSize = p_allreduce(SIZE(providedGlobalIndices), p_sum_op(), p_comm_work)
        procCount = p_comm_size(p_comm_work)
        myRank = p_comm_rank(p_comm_work)

        ! consistency checks
        IF(1 > MINVAL(providedGlobalIndices)) THEN
          WRITE (0,*) "MINVAL(providedGlobalIndices) = ", MINVAL(providedGlobalIndices)
          CALL finish(routine//"-"//TRIM(name), "assertion failed: a global index is out of RANGE (<1)")
        END IF
        IF(globalSize < MAXVAL(providedGlobalIndices)) THEN
          WRITE (0,*) "MAXVAL(providedGlobalIndices) = ", MAXVAL(providedGlobalIndices), " > ", globalSize
          CALL finish(routine//"-"//TRIM(name), "assertion failed: a global index is out of range (>globalSize)")
        END IF

        ! define simple broker decomposition that IS known perfectly
        ! on all PEs (consecutive number ranges of similar SIZE)
        ! process rank IS responsible for the global index range
        ! [brokerBounds(rank) + 1, brokerBounds(rank + 1)]
        ALLOCATE(brokerBounds(0:procCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine//"-"//TRIM(name), "memory allocation failure")
        DO i = 0, procCount
            brokerBounds(i) = INT(INT(i, i8)*INT(globalSize, i8)/INT(procCount, i8))
        END DO

        ! create communication patterns to interact with that broker decomposition
        CALL providerToBroker%construct(providedGlobalIndices, brokerBounds)
        CALL consumerToBroker%construct(requiredGlobalIndices, brokerBounds)

        ! communicate the provider processes to the consumers
        ALLOCATE(providerBuffer(SIZE(providedGlobalIndices)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine//"-"//TRIM(name), "memory allocation failure")
        providerBuffer(:) = myRank
        ALLOCATE(brokerBuffer(brokerBounds(myRank + 1) - brokerBounds(myRank)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine//"-"//TRIM(name), "memory allocation failure")
        ALLOCATE(providerOfRequiredPoint(SIZE(requiredGlobalIndices)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine//"-"//TRIM(name), "memory allocation failure")
        CALL providerToBroker%communicateToBroker(providerBuffer, brokerBuffer)
        CALL consumerToBroker%communicateFromBroker(brokerBuffer, providerOfRequiredPoint)

        resultVar => makeCommPattern(globalSize, providedGlobalIndices, requiredGlobalIndices, &
          &                         providerOfRequiredPoint)

        ! cleanup
        CALL providerToBroker%destruct()
        CALL consumerToBroker%destruct()
        DEALLOCATE(brokerBounds, providerBuffer, brokerBuffer, providerOfRequiredPoint)

        ! sanity check (yes, this _is_ defensive)
        CALL checkRedistributionPattern(resultVar, providedGlobalIndices, requiredGlobalIndices)

        IF(timers_level >= 7) CALL timer_stop(timer_load_restart_comm_setup)
    END FUNCTION makeRedistributionPattern

    SUBROUTINE multifilePatchReader_construct(me, p_patch, multifilePath)
        CLASS(t_MultifilePatchReader), TARGET, INTENT(INOUT) :: me
        TYPE(t_patch), INTENT(in) :: p_patch
        CHARACTER(*), INTENT(IN) :: multifilePath

        INTEGER :: multifile_file_count, error, curFileIndex
        INTEGER :: totalCellCount, totalEdgeCount, totalVertCount

        ! these contain the global indices of all cells/edges/verts
        ! contained IN all of the files READ by this process
        INTEGER, ALLOCATABLE, TARGET :: globalCellIndices(:), globalEdgeIndices(:), globalVertIndices(:)

        INTEGER, POINTER :: indexWindow(:)
        INTEGER :: curCellOffset, curEdgeOffset, curVertOffset
        TYPE(t_PayloadFile), POINTER :: curFile
        TYPE(t_RestartAttributeList), POINTER :: restartAttributes
        CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchReader_construct"

        restartAttributes => getAttributesForRestarting()
        multifile_file_count = restartAttributes%getInteger('multifile_file_count')

        IF (my_process_is_mpi_workroot()) THEN
          WRITE (0,*) "reading from ", multifile_file_count, " files/patch."
        END IF

        me%domain = p_patch%id
        call openPayloadFiles(multifilePath, multifile_file_count, me%files, p_patch%id, &
                                   &totalCellCount, totalEdgeCount, totalVertCount)

        ! Allocate the global index arrays.
        ALLOCATE(globalCellIndices(totalCellCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(globalEdgeIndices(totalEdgeCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(globalVertIndices(totalVertCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        ! Read the global indices.

        curCellOffset = 1
        curEdgeOffset = 1
        curVertOffset = 1
        DO curFileIndex = 1, SIZE(me%files)
            curFile => me%files(curFileIndex)

            ! If the restart file has been created in a run with
            ! processor splitting enabled, then it may be that a
            ! payload file is virtually "empty" (0 cells).

            IF (curFile%cellCount > 0) THEN
              indexWindow => globalCellIndices(curCellOffset : curCellOffset + curFile%cellCount - 1)
              CALL curFile%readIndexVar(GRID_UNSTRUCTURED_CELL, indexWindow)
              curCellOffset = curCellOffset + curFile%cellCount
            END IF

            IF (curFile%edgeCount > 0) THEN
              indexWindow => globalEdgeIndices(curEdgeOffset : curEdgeOffset + curFile%EdgeCount - 1)
              CALL curFile%readIndexVar(GRID_UNSTRUCTURED_EDGE, indexWindow)
              curEdgeOffset = curEdgeOffset + curFile%EdgeCount
            END IF

            IF (curFile%vertCount > 0) THEN
              indexWindow => globalVertIndices(curVertOffset : curVertOffset + curFile%VertCount - 1)
              CALL curFile%readIndexVar(GRID_UNSTRUCTURED_VERT, indexWindow)
              curVertOffset = curVertOffset + curFile%VertCount
            END IF

        END DO

        ! Init the READ buffers / DATA distributors.
        CALL me%readBufferCells%construct("cells", globalCellIndices, p_patch%cells%decomp_info%glb_index)
        CALL me%readBufferEdges%construct("edges", globalEdgeIndices, p_patch%edges%decomp_info%glb_index)
        CALL me%readBufferVerts%construct("verts", globalVertIndices, p_patch%verts%decomp_info%glb_index)

        ! Cleanup.
        DEALLOCATE(globalCellIndices, globalEdgeIndices, globalVertIndices)
    END SUBROUTINE multifilePatchReader_construct

    SUBROUTINE multifilePatchReader_readData(me, varData)
        CLASS(t_MultifilePatchReader), TARGET, INTENT(INOUT) :: me
        TYPE(t_restartVarData), INTENT(INOUT) :: varData(:)

        INTEGER                         :: varIndex, level, levelCount, fileIndex
        TYPE(t_ReadBuffer), POINTER     :: readBuffer
        TYPE(t_ptr_2d_sp), ALLOCATABLE  :: levelPointers_s(:)
        TYPE(t_ptr_2d), ALLOCATABLE     :: levelPointers_d(:)
        TYPE(t_ptr_2d_int), ALLOCATABLE :: levelPointers_int(:)
        CHARACTER(LEN = VARNAME_LEN)    :: curVarname
        CHARACTER(*), PARAMETER         :: routine = modname//":multifilePatchReader_readData"

        DO varIndex = 1, SIZE(varData)
            ! Skip variables that are NOT to be expected IN the restart file.
            IF(.NOT.has_valid_time_level(varData(varIndex)%info, me%domain, nnew(me%domain), &
              &                          nnew_rcf(me%domain))) CYCLE

            ! Check that all processes have a consistent order of variables IN varData(:).
            curVarname = varData(varIndex)%info%NAME

            CALL p_bcast(curVarname, 0, p_comm_work)
            IF(curVarname /= varData(varIndex)%info%NAME) THEN
                CALL finish(routine, "assertion failed: inconsistent order of varData(:) array")
            END IF

            ! Get the correct READ buffer object.
            SELECT CASE(varData(varIndex)%info%hgrid)
                CASE(GRID_UNSTRUCTURED_CELL)
                    readBuffer => me%readBufferCells
                CASE(GRID_UNSTRUCTURED_EDGE)
                    readBuffer => me%readBufferEdges
                CASE(GRID_UNSTRUCTURED_VERT)
                    readBuffer => me%readBufferVerts
                CASE DEFAULT
                    CALL finish(routine, "assertion failed: illegal VALUE of varData(varIndex)%info%hgrid")
            END SELECT

            ! Get the pointers to the final destinations of the data.
            SELECT CASE(varData(varIndex)%info%data_type)
            CASE(REAL_T)
              !
              CALL getLevelPointers(varData(varIndex)%info, varData(varIndex)%r_ptr, levelPointers_d)
              levelCount = SIZE(levelPointers_d)
              !
            CASE(SINGLE_T)
              CALL getLevelPointers(varData(varIndex)%info, varData(varIndex)%s_ptr, levelPointers_s)
              levelCount = SIZE(levelPointers_s)
              !
            CASE(INT_T)
              CALL getLevelPointers(varData(varIndex)%info, varData(varIndex)%i_ptr, levelPointers_int)
              levelCount = SIZE(levelPointers_int)
              !
            CASE DEFAULT
              CALL finish(routine, "Internal error! Variable "//TRIM(varData(varIndex)%info%name))
            END SELECT

            DO level = 1, levelCount
                ! Read the DATA of one level into the READ buffer.
                DO fileIndex = 1, SIZE(me%files)
                    CALL me%files(fileIndex)%readVarLevel(varData(varIndex), level - 1, readBuffer)
                END DO

                ! Redistribute the DATA, writing it to its FINAL destination. This also sets the halo points.
                SELECT CASE(varData(varIndex)%info%data_type)
                CASE(REAL_T)
                  CALL readBuffer%redistribute(levelPointers_d(level)%p)
                  !
                CASE(SINGLE_T)
                  CALL readBuffer%redistribute(levelPointers_s(level)%p)
                  !
                CASE(INT_T)
                  CALL readBuffer%redistribute(levelPointers_int(level)%p)
                  !
                CASE DEFAULT
                  CALL finish(routine, "Internal error! Variable "//TRIM(varData(varIndex)%info%name))
                END SELECT
            END DO
        END DO
    END SUBROUTINE multifilePatchReader_readData

    SUBROUTINE multifilePatchReader_destruct(me)
        CLASS(t_MultifilePatchReader), INTENT(INOUT) :: me

        INTEGER :: i

        CALL me%readBufferCells%destruct()
        CALL me%readBufferEdges%destruct()
        CALL me%readBufferVerts%destruct()
        DO i = 1, SIZE(me%files)
            CALL me%files(i)%destruct()
        END DO

        DEALLOCATE(me%files)
    END SUBROUTINE multifilePatchReader_destruct

    SUBROUTINE multifileReadPatch(varData, p_patch, multifilePath)
        TYPE(t_restartVarData), INTENT(INOUT) :: varData(:)
        TYPE(t_patch), INTENT(in) :: p_patch
        CHARACTER(*), INTENT(IN) :: multifilePath

        TYPE(t_MultifilePatchReader) :: reader

        CALL reader%construct(p_patch, multifilePath)
        CALL reader%readData(varData)
        CALL reader%destruct()
      END SUBROUTINE multifileReadPatch



END MODULE mo_load_multifile_restart
