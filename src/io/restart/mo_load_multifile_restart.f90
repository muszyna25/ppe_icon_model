!> Module for reading multifile restart files
!!
!! 2018-08: Major revision / revamp / refactoring : Harald Braun (Atos SE)
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
    &                                  streamClose, streamReadVar, cdiInqAttInt,  CDI_GLOBAL,      &
    &                                  streamReadVarSlice, streamReadVarSliceF
! TODO : get rid of mo_communication_orig -> use mo_communication instead
  USE mo_communication_orig,     ONLY: t_comm_pattern_orig, exchange_data_noblk
  USE mo_decomposition_tools,    ONLY: t_glb2loc_index_lookup, init_glb2loc_index_lookup, set_inner_glb_index, &
    &                                  deallocate_glb2loc_index_lookup
  USE mo_dynamics_config,        ONLY: nnew, nnew_rcf
  USE mo_exception,              ONLY: finish, warning
  USE mo_kind,                   ONLY: sp, dp, i8
  USE mo_model_domain,           ONLY: t_patch
  USE mo_mpi,                    ONLY: p_barrier, p_comm_work, p_comm_size, p_comm_rank, my_process_is_work, &
    &                                  p_allreduce, p_sum_op, p_mpi_wtime, my_process_is_stdio, p_bcast,     &
    &                                  my_process_is_mpi_workroot
  USE mo_multifile_restart_util, ONLY: multifilePayloadPath, commonBuf_t, dataPtrs_t, rBuddy, rGroup, &
    &                                  vNames_glbIdx
  USE mo_parallel_config,        ONLY: nproma, idx_no, blk_no
  USE mo_restart_attributes,     ONLY: t_RestartAttributeList, getAttributesForRestarting, ocean_initFromRestart_OVERRIDE
  USE mo_restart_var_data,       ONLY: t_restartVarData, getLevelPointers, has_valid_time_level
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_load_restart_io, timers_level, &
    &                                  timer_load_restart_comm_setup, timer_load_restart_communication, &
    &                                  timer_load_restart_get_var_id
  USE mo_util_cdi,               ONLY: get_cdi_varID, test_cdi_varID
  USE mo_util_string,            ONLY: charArray_equal, real2string

  IMPLICIT NONE

  PUBLIC :: multifileCheckRestartFiles, multifileReadPatch

  TYPE t_PayloadFile
    INTEGER :: streamId, vlistId, varCount
    INTEGER :: iCnts(3), iVarIds(3)
  CONTAINS
    PROCEDURE :: construct    => payloadFile_construct
    PROCEDURE :: readIndexVar => payloadFile_readIndexVar
    PROCEDURE :: readVarLevel => payloadFile_readVarLevel
    PROCEDURE :: destruct     => payloadFile_destruct
  END TYPE t_PayloadFile

  ! This IS both a buffer AND a tool for communication:
  ! First, the buffer IS filled locally, THEN the DATA IS
  ! redistributed among the processes AND written to the provided
  ! destination.
  TYPE t_ReadBuffer
    INTEGER, PRIVATE :: curWindowOffset
    TYPE(commonBuf_t) :: readBuf_1d
    TYPE(t_comm_pattern_orig), POINTER, PRIVATE :: commPattern 
  CONTAINS
    PROCEDURE :: construct => readBuffer_construct
    PROCEDURE, PRIVATE :: redistribute => readBuffer_redistribute
    PROCEDURE :: destruct => readBuffer_destruct
  END TYPE t_ReadBuffer

  ! A MultifilePatchReader IS constructed for a single patch, AND IS
  ! responsible to READ all the ASSOCIATED patch<domain>_<N>.nc
  ! files, redistributing their contents to the correct processes.
  TYPE t_MultifilePatchReader
    TYPE(t_ReadBuffer), PRIVATE :: readBufCells, readBufEdges, readBufVerts
    TYPE(t_PayloadFile), PRIVATE, ALLOCATABLE :: files(:)
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
    INTEGER, VALUE :: domain, partId
    INTEGER :: varId, dummy, icount(1), iG
    LOGICAL :: haveCount(3)
    CHARACTER(KIND = C_CHAR), POINTER :: nameCharArray(:)
    CHARACTER(*), PARAMETER :: routine = modname//":payloadFile_construct"
    CHARACTER(*), PARAMETER :: attNames(3) = (/"cellCount", "edgeCount", "vertCount"/)
    CHARACTER(:), ALLOCATABLE :: pathname

    IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
    ! open the file AND cache the most important CDI IDs
    CALL multifilePayloadPath(multifilePath, domain, partId, pathname)
    me%streamId = streamOpenRead(pathname)
    IF (me%streamId < 0) WRITE (0,*) "failed to open ", pathname, " for reading"
    me%vlistId = streamInqVlist(me%streamId)
    me%varCount = vlistNvars(me%vlistId)
    haveCount(:) = .false.
    DO varId = 0, me%varCount - 1
      IF(.NOT.ANY(.NOT.haveCount(:))) CYCLE
      nameCharArray => vlistCopyVarName(me%vlistId, varId)
      DO iG = 1, 3
        IF(charArray_equal(nameCharArray, vNames_glbIdx(iG))) THEN
          IF(haveCount(iG)) &
            & CALL finish(routine, &
              & "corrupted restart file: double def "//vNames_glbIdx(iG))
          dummy = cdiInqAttInt(me%vlistId, CDI_GLOBAL, attNames(iG), 1, icount);
          me%iCnts(iG) = icount(1)
          me%iVarIds(iG) = varId
          haveCount(iG) = .true.
        END IF
      END DO
      DEALLOCATE(nameCharArray)
    END DO
    IF(timers_level >= 7) CALL timer_stop(timer_load_restart_io)
    IF(ANY(.NOT.haveCount(:))) &
      & CALL finish(routine, "corrupted restart file: an index variable is missing")
  END SUBROUTINE payloadFile_construct

  SUBROUTINE payloadFile_readIndexVar(me, pCnt, vId, glbIdces)
    CLASS(t_PayloadFile), INTENT(IN) :: me
    INTEGER, INTENT(IN) :: pCnt, vId
    INTEGER, INTENT(INOUT), CONTIGUOUS :: glbIdces(:)
    INTEGER :: error, i
    INTEGER(C_INT) :: trash
    REAL(dp), ALLOCATABLE :: buffer(:)
    CHARACTER(*), PARAMETER :: routine = modname//":payloadFile_readIndexVar"

    IF (SIZE(glbIdces) /= pCnt) &
      & CALL finish(routine, "globalIndices(:) has the wrong size")
    IF (pCnt .LE. 0) RETURN !nothing to do
    ALLOCATE(buffer(pCnt), STAT = error)
    IF (error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    IF (timers_level >= 7) CALL timer_start(timer_load_restart_io)
    CALL streamReadVar(me%streamId, vId, buffer, trash)
    IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
!$OMP PARALLEL DO SCHEDULE(STATIC)
    DO i = 1, pCnt
      glbIdces(i) = INT(buffer(i))
    END DO
    DEALLOCATE(buffer)
  END SUBROUTINE payloadFile_readIndexVar

  ! The offset IS 1 based AND will be updated to point to the next
  ! region within the buffer that should be filled.
  SUBROUTINE payloadFile_readVarLevel(me, varData, varId, lev, buf, pCt)
    CLASS(t_PayloadFile), INTENT(INOUT) :: me
    TYPE(t_restartVarData), INTENT(IN) :: varData
    INTEGER, VALUE :: lev, varId, pCt
    TYPE(t_ReadBuffer), INTENT(INOUT) :: buf
    INTEGER :: i, trash, windowEnd
    REAL(KIND=dp), POINTER :: pdp(:), buffer_dp(:)
    REAL(KIND=sp), POINTER :: psp(:)
    INTEGER,       POINTER :: pint(:)
    CHARACTER(*), PARAMETER :: routine = modname//":payloadFile_readVarLevel"

    IF (pCt .LE. 0) RETURN ! nothing to do
    windowEnd = buf%curWindowOffset + pCt - 1
    IF(windowEnd > SIZE(buf%readBuf_1d%d)) THEN
      WRITE (0,*) "windowEnd = ", windowEnd, &
        &         "; SIZE(buf%readBuf_1d%d) = ", SIZE(buf%readBuf_1d%d)
      CALL finish(routine, "assertion failed: attempt to overfill buffer")
    END IF
    IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
    SELECT CASE(varData%info%data_type)
    CASE(REAL_T)
      pdp  => buf%readBuf_1d%d(buf%curWindowOffset:windowEnd)
      CALL streamReadVarSlice(me%streamId, varId, lev, pdp, trash)
    CASE(SINGLE_T)
      psp  => buf%readBuf_1d%s(buf%curWindowOffset:windowEnd)
      CALL streamReadVarSliceF(me%streamId, varId, lev, psp, trash)
    CASE(INT_T)
      pint => buf%readBuf_1d%i(buf%curWindowOffset:windowEnd)
      ALLOCATE(buffer_dp(pCt))
      CALL streamReadVarSlice(me%streamId, varId, lev, buffer_dp, trash)
!$OMP PARALLEL DO SCHEDULE(STATIC)
      DO i = 1, pCt
        pint(i) = INT(buffer_dp(i))
      END DO
      DEALLOCATE(buffer_dp)
    CASE DEFAULT
      CALL finish(routine, "Internal error! Variable "//TRIM(varData%info%name))
    END SELECT
    IF(timers_level >= 7) CALL timer_stop(timer_load_restart_io)
    buf%curWindowOffset = windowEnd + 1
  END SUBROUTINE payloadFile_readVarLevel

  SUBROUTINE payloadFile_destruct(me)
    CLASS(t_PayloadFile), INTENT(INOUT) :: me

    IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
    CALL streamClose(me%streamId)
    IF(timers_level >= 7) CALL timer_stop(timer_load_restart_io)
  END SUBROUTINE payloadFile_destruct

  SUBROUTINE readBuffer_construct(me, provGlbIdces, reqdGlbIdces)
    CLASS(t_ReadBuffer), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: provGlbIdces(:), reqdGlbIdces(:)
    INTEGER :: pCnt, ierr
    CHARACTER(*), PARAMETER :: routine = modname//":readBuffer_construct"

    ! Note: These communication patterns also set the halo points.
    me%commPattern => makeRedistributionPattern(provGlbIdces, reqdGlbIdces)
    pCnt = SIZE(provGlbIdces)
    ALLOCATE(me%readBuf_1d%d(pCnt), me%readBuf_1d%s(pCnt), me%readBuf_1d%i(pCnt), STAT = ierr)
    IF(ierr /= SUCCESS) CALL finish(routine, "memory allocation failure")
    me%curWindowOffset = 1
  END SUBROUTINE readBuffer_construct

  SUBROUTINE readBuffer_redistribute(me, dest, lev)
    CLASS(t_ReadBuffer), INTENT(INOUT) :: me
    TYPE(dataPtrs_t), INTENT(INOUT) :: dest
    INTEGER, INTENT(IN) :: lev

    IF(timers_level >= 7) CALL timer_start(timer_load_restart_communication)
    IF(ALLOCATED(dest%s)) &
      CALL exchange_data_noblk(me%commPattern, dest%s(lev)%p, me%readBuf_1d%s)
    IF(ALLOCATED(dest%d)) &
      CALL exchange_data_noblk(me%commPattern, dest%d(lev)%p, me%readBuf_1d%d)
    IF(ALLOCATED(dest%i)) &
      CALL exchange_data_noblk(me%commPattern, dest%i(lev)%p, me%readBuf_1d%i)
    IF(timers_level >= 7) CALL timer_stop(timer_load_restart_communication)
    me%curWindowOffset = 1
  END SUBROUTINE readBuffer_redistribute

  SUBROUTINE readBuffer_destruct(me)
    CLASS(t_ReadBuffer), INTENT(INOUT) :: me
    CHARACTER(*), PARAMETER :: routine = modname//":readBuffer_destruct"

    IF(me%curWindowOffset /= 1) CALL finish(routine, "buffer NOT empty on destruction")
    CALL me%commPattern%delete()
    DEALLOCATE(me%commPattern)
  END SUBROUTINE readBuffer_destruct

  SUBROUTINE multifileCheckRestartFiles(filename)
    CHARACTER(*), INTENT(IN) :: filename
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    INTEGER :: n_dom, multifile_file_count
    CHARACTER(*), PARAMETER :: routine = modname//":multifileCheckRestartFiles"

    restartAttributes => getAttributesForRestarting()
    n_dom = restartAttributes%getInteger('multifile_n_dom_active')
    multifile_file_count = restartAttributes%getInteger('multifile_file_count')
    IF(checkMultifileDir(filename, n_dom, multifile_file_count) /= SUCCESS) THEN
      CALL finish(routine, "'"//filename//"' is not a valid restart multifile")
    END IF
  END SUBROUTINE multifileCheckRestartFiles

  ! Opens all the payload files that should be opened by this
  ! process, creating t_PayloadFile objects for them, AND
  ! calculating the total numbers of cells/edges/vertices contained
  ! IN these files.
  !
  ! This FUNCTION performs the distribution of the payload files
  ! among the restart processes, subsequent code just has to iterate
  ! over the array returned by this FUNCTION.
  SUBROUTINE openPayloadFiles(mfPath, mfCnt, dom, tCnt, f)
    TYPE(t_PayloadFile), ALLOCATABLE, TARGET, INTENT(INOUT) :: f(:)
    CHARACTER(*), INTENT(IN) :: mfPath
    INTEGER, INTENT(IN) :: mfCnt, dom
    INTEGER, INTENT(OUT) :: tCnt(3)
    INTEGER :: pCnt, myR, myFCnt, myFF, i, cFId
    CHARACTER(*), PARAMETER :: routine = modname//":openPayloadFiles"

    IF(.NOT.my_process_is_work()) CALL finish(routine, "this is not a work process")
    pCnt = p_comm_size(p_comm_work)
    myR = p_comm_rank(p_comm_work)
    IF (pCnt .GE. mfCnt) THEN
      myFCnt = MERGE(1, 0, rBuddy(nr_in=mfCnt) .EQ. myR)
      myFF = MERGE(rGroup(nr_in=mfCnt), -1, myFCnt .GT. 0)
    ELSE
      myFCnt = COUNT( (/ (rGroup(nr_in=pCnt, nw_in=mfCnt, pe_in=i) .EQ. myR, &
        &                    i = 0, mfCnt-1) /) )
      myFF = -1
      DO i = 0, mfCnt-1
        IF(myFF .GE. 0) CYCLE
        myFF = MERGE(i, -1, rGroup(nr_in=pCnt, nw_in=mfCnt, pe_in=i).EQ. myR)
      END DO
    END IF
#ifdef DEBUG
    WRITE(0, "(3(a,i4),a)") "restart: rank ", myR, " will start reading from file ", &
      &                      myFF, ", ", myFCnt, " files in total"
#endif
    IF (ALLOCATED(f)) DEALLOCATE(f)
    ALLOCATE(f(myFCnt), STAT = i)
    IF(i /= SUCCESS) CALL finish(routine, "memory allocation failure")
    tCnt(:) = 0
    DO i = 1, myFCnt
      cFId = myFF - 1 + i
      CALL f(i)%construct(mfPath, dom, cFId)
      tCnt(:) = tCnt(:) + f(i)%iCnts(:)
    END DO
  END SUBROUTINE openPayloadFiles

  ! This asserts that the given pattern redistributes a given array
  ! exactly as described by the global indices arrays.
  SUBROUTINE checkRedistributionPattern(pattern, provGlbIdces, reqdGlbIdces)
    TYPE(t_comm_pattern_orig), INTENT(IN) :: pattern
    INTEGER, INTENT(IN) :: provGlbIdces(:), reqdGlbIdces(:)
    INTEGER :: iSize, oSize, i, error
    INTEGER, ALLOCATABLE :: input(:,:), output(:,:)
    CHARACTER(*), PARAMETER :: routine = modname//":checkRedistributionPattern"

    ! ALLOCATE some memory
    iSize = SIZE(provGlbIdces)
    oSize = SIZE(reqdGlbIdces)
    ALLOCATE(input(nproma, blk_no(iSize)), output(nproma, blk_no(oSize)), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
!$OMP PARALLEL DO SCHEDULE(STATIC)
    DO i = 1, blk_no(oSize)
      output(:,i) = -1
    END DO
    CALL exchange_data_noblk(pattern, output, provGlbIdces)
    ! assert that the RESULT IS as expected
!$OMP PARALLEL DO SCHEDULE(STATIC)
    DO i = 1, oSize
      IF(output(idx_no(i), blk_no(i)) /= REAL(reqdGlbIdces(i))) THEN
!$OMP CRITICAL
        CALL finish(routine, &
          & "assertion failed: redistribution pattern does not work as expected")
!$OMP END CRITICAL
      END IF
    END DO
    DEALLOCATE(input, output)
  END SUBROUTINE checkRedistributionPattern

  ! Create a t_comm_pattern that IS able to redistribute the DATA as
  ! it IS READ from the multifile to the domain decomposition of the
  ! work processes.  The tricky part about this IS, that
  ! setup_comm_pattern() requires an owner(:) array as input, which
  ! we have to create first.
  !
  ! The resulting t_comm_pattern will also initialize the halo points.
  FUNCTION makeRedistributionPattern(provGlbIdces, reqdGlbIdces) RESULT(resVar)
    TYPE(t_comm_pattern_orig), POINTER :: resVar
    INTEGER, INTENT(IN) :: provGlbIdces(:), reqdGlbIdces(:)
    INTEGER :: globalSize, procCount, myRank, error, i
    INTEGER, ALLOCATABLE :: brokBnds(:), provBuf(:), brokBuf(:), &
      &  provOfReqdPt(:), local_owner(:)
    TYPE(t_BrokerCommunicationPattern) :: providerToBroker, consumerToBroker
    TYPE(t_glb2loc_index_lookup) :: lookupTable
    CHARACTER(*), PARAMETER :: routine = modname//":makeRedistributionPattern"

    IF(timers_level >= 7) CALL timer_start(timer_load_restart_comm_setup)
    ! compute the global SIZE of the field
    globalSize = p_allreduce(SIZE(provGlbIdces), p_sum_op(), p_comm_work)
    procCount = p_comm_size(p_comm_work)
    myRank = p_comm_rank(p_comm_work)
#ifdef DEBUG 
    ! consistency checks
    IF(1 > MINVAL(provGlbIdces)) &
      & CALL finish(routine, "global index is out of RANGE (<1)")
    IF(globalSize < MAXVAL(provGlbIdces)) &
      &  CALL finish(routine, "global index is out of range (>globalSize)")
#endif
    ! define simple broker decomposition that IS known perfectly
    ! on all PEs (consecutive number ranges of similar SIZE)
    ! process rank IS responsible for the global index range
    ! [brokerBounds(rank) + 1, brokerBounds(rank + 1)]
    ALLOCATE(brokBnds(0:procCount), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    DO i = 0, procCount
      brokBnds(i) = INT(INT(i, i8)*INT(globalSize, i8)/INT(procCount, i8))
    END DO
    ! create communication patterns to interact with that broker decomposition
    CALL providerToBroker%construct(provGlbIdces, brokBnds)
    CALL consumerToBroker%construct(reqdGlbIdces, brokBnds)
    ! communicate the provider processes to the consumers
    ALLOCATE(brokBuf(brokBnds(myRank + 1) - brokBnds(myRank)), &
      & provOfReqdPt(SIZE(reqdGlbIdces)), provBuf(SIZE(provGlbIdces)), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    provBuf(:) = myRank
    CALL providerToBroker%communicateToBroker(provBuf, brokBuf)
    CALL consumerToBroker%communicateFromBroker(brokBuf, provOfReqdPt)

    IF (SIZE(reqdGlbIdces) /= SIZE(provOfReqdPt)) &
      & CALL finish(routine, "reqdGlbIdces(:) and owners(:) " // &
          & "must be of same size")
    CALL init_glb2loc_index_lookup(lookupTable, globalSize)
    CALL set_inner_glb_index(lookupTable, provGlbIdces, &
      & [(i, i = 1, SIZE(provGlbIdces))])
    ALLOCATE(local_owner(SIZE(provGlbIdces)))
    local_owner(:) = p_comm_rank(p_comm_work)
    ALLOCATE(resVar)
    CALL resVar%setup(SIZE(reqdGlbIdces), provOfReqdPt, reqdGlbIdces, &
      & lookupTable, SIZE(provGlbIdces), local_owner, provGlbIdces)
    CALL deallocate_glb2loc_index_lookup(lookupTable)
    CALL providerToBroker%destruct()
    CALL consumerToBroker%destruct()
    DEALLOCATE(brokBnds, provBuf, brokBuf, provOfReqdPt, local_owner)
#ifdef DEBUG
    ! sanity check (yes, this _is_ defensive)
    CALL checkRedistributionPattern(resVar, provGlbIdces, reqdGlbIdces)
#endif
    IF(timers_level >= 7) CALL timer_stop(timer_load_restart_comm_setup)
  END FUNCTION makeRedistributionPattern

  SUBROUTINE multifilePatchReader_construct(me, p_patch, multifilePath)
    CLASS(t_MultifilePatchReader), TARGET, INTENT(INOUT) :: me
    TYPE(t_patch), INTENT(in) :: p_patch
    CHARACTER(*), INTENT(IN) :: multifilePath
    INTEGER :: mfileCnt, ierr, cFId, tCnt(3), cOff(3), iG
    TYPE i_ptr_arr_t
      INTEGER, POINTER :: p(:)
    END TYPE i_ptr_arr_t
    TYPE(i_ptr_arr_t) :: glbIdx(3)
    INTEGER, POINTER :: idxPtr(:)
    TYPE(t_PayloadFile), POINTER :: cFile
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchReader_construct"

    restartAttributes => getAttributesForRestarting()
    mfileCnt = restartAttributes%getInteger('multifile_file_count')
    IF (my_process_is_mpi_workroot()) &
      & WRITE(0, *) "reading from ", mfileCnt, " files/patch."
    CALL openPayloadFiles(multifilePath, mfileCnt, p_patch%id, tCnt, me%files)
    ALLOCATE(glbIdx(1)%p(tCnt(1)), glbIdx(2)%p(tCnt(2)), glbIdx(3)%p(tCnt(3)), STAT = ierr)
    IF(ierr /= SUCCESS) CALL finish(routine, "memory allocation failure")
    cOff(:) = 1
    DO cFId = 1, SIZE(me%files)
      cFile => me%files(cFId)
      DO iG = 1, 3
        IF (cFile%iCnts(iG) > 0) THEN
          idxPtr => glbIdx(iG)%p(cOff(iG) : cOff(iG) + cFile%iCnts(iG) - 1)
          CALL cFile%readIndexVar(cFile%iCnts(iG), cFile%iVarIds(iG), idxPtr)
          cOff(iG) = cOff(iG) + cFile%iCnts(iG)
        END IF
      END DO
    END DO
    CALL me%readBufCells%construct(glbIdx(1)%p, p_patch%cells%decomp_info%glb_index)
    CALL me%readBufEdges%construct(glbIdx(2)%p, p_patch%edges%decomp_info%glb_index)
    CALL me%readBufVerts%construct(glbIdx(3)%p, p_patch%verts%decomp_info%glb_index)
    DEALLOCATE(glbIdx(1)%p, glbIdx(2)%p, glbIdx(3)%p)
  END SUBROUTINE multifilePatchReader_construct

  SUBROUTINE multifilePatchReader_readData(me, vDat, dom)
    CLASS(t_MultifilePatchReader), TARGET, INTENT(INOUT) :: me
    TYPE(t_restartVarData), INTENT(INOUT) :: vDat(:)
    INTEGER, INTENT(IN) :: dom
    INTEGER :: vId, lId, lCt, fId, varID, svDat, pCt(SIZE(me%files))
    TYPE(t_ReadBuffer), POINTER :: rBuf
    TYPE(dataPtrs_t) :: lPtrs
#ifdef DEBUG
    CHARACTER(LEN=VARNAME_LEN) :: cVname
#endif
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchReader_readData"

    svDat = SIZE(vDat)
    DO vId = 1, SIZE(vDat)
      IF(.NOT.has_valid_time_level(vDat(vId)%info, dom, nnew(dom), nnew_rcf(dom))) CYCLE
      ! Check that all processes have a consistent order of variables IN varData(:).
      IF(timers_level >= 7) CALL timer_start(timer_load_restart_get_var_id)
      IF(SIZE(me%files) .GT. 0) THEN
! fatal hack from coding hell to make init_fromRestart=.true. (ocean) work
        IF(ocean_initFromRestart_OVERRIDE) THEN 
          varId = test_cdi_varID(me%files(1)%streamId, TRIM(vDat(vId)%info%NAME))
          IF(varId .eq. -1) THEN
            IF(timers_level >= 7) CALL timer_stop(timer_load_restart_get_var_id)
            CALL warning(routine, "variable '" // TRIM(vDat(vId)%info%NAME) // &
              & "' from restart file not found in the list of restart variables")
            CALL warning(routine, &
              & "that MAY be intended if initialize_fromRestart=.true.")
            CYCLE
          END IF
        ELSE
          varId = test_cdi_varID(me%files(1)%streamId, TRIM(vDat(vId)%info%NAME))
        END IF
      END IF
      IF(timers_level >= 7) CALL timer_stop(timer_load_restart_get_var_id)
#ifdef DEBUG
      cVname = vDat(vId)%info%NAME
      CALL p_bcast(cVname, 0, p_comm_work)
      IF(cVname /= vDat(vId)%info%NAME) &
        & CALL finish(routine, "inconsistent order of varData(:) array")
      varId0 = varId
      CALL p_bcast(varId0, 0, p_comm_work)
      IF(varId0 /= varId .AND. SIZE(me%files) .GT. 0) &
        & CALL finish(routine, "inconsistent order of varData(:) array")
#endif
      SELECT CASE(vDat(vId)%info%hgrid)
      CASE(GRID_UNSTRUCTURED_CELL)
        rBuf => me%readBufCells
        pCt = me%files(:)%iCnts(1)
      CASE(GRID_UNSTRUCTURED_EDGE)
        rBuf => me%readBufEdges
        pCt = me%files(:)%iCnts(2)
      CASE(GRID_UNSTRUCTURED_VERT)
        rBuf => me%readBufVerts
        pCt = me%files(:)%iCnts(3)
      CASE DEFAULT
        CALL finish(routine, "illegal varData(varIndex)%info%hgrid")
      END SELECT
      CALL getLevelPointers(vDat(vId)%info, vDat(vId), lPtrs, lCnt=lCt)
      DO lId = 1, lCt
        DO fId = 1, SIZE(me%files)
          CALL me%files(fId)%readVarLevel(vDat(vId), varId, &
            &           lId - 1, rBuf, pCt(fId))
        END DO
        CALL rBuf%redistribute(lPtrs, lId)
      END DO
    END DO
    CALL lPtrs%free()
  END SUBROUTINE multifilePatchReader_readData

  SUBROUTINE multifilePatchReader_destruct(me)
    CLASS(t_MultifilePatchReader), INTENT(INOUT) :: me
    INTEGER :: i

    CALL me%readBufCells%destruct()
    CALL me%readBufEdges%destruct()
    CALL me%readBufVerts%destruct()
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
    CALL reader%readData(varData, p_patch%id)
    CALL reader%destruct()
  END SUBROUTINE multifileReadPatch

END MODULE mo_load_multifile_restart
