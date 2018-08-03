!> Module for writing restart as a multifile (synchronously)
!!
!! Note: Other implementations of the restart writing interface can be
!! found in mo_sync_restart and mo_async_restart.
!!
!! Initial implementation: Nathanael Huebbe
!! Refactoring (point-to-point comm. replaced by one-sided comm.): F. Prill
!! 2018-08: Major revision / revamp / refactoring : Harald Braun (Atos SE)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! --------------------------------------------------------------------------------
!!
!! Generated files:
!! ================
!!
!! Each restart is contained within a dedicated directory. The
!! directory itself is called the restart file.  The structure within
!! the directory is entirely internal to ICON and should be of no
!! concern to users.  By default, this directory is given the filename
!! extension "mfr" for "MultiFile Restart".
!!
!! The internal structure is as follows:
!!
!!   * restartFile.mfr/attributes.nc: 
!!     NetCDF file that contains the restart attributes and serialized
!!     namelists, this contains the domain count
!!
!!   * restartFile.mfr/patch<JG>_<N>.nc: 
!!     Data, stored as 2D or 3D arrays of doubles/singles. Each
!!     generating process streams its data to exactly one restart
!!     process, which writes it to exactly one file. Each file
!!     contains data from exactly one restart process, which received
!!     it from floor(num_work_procs/num_restart_procs) to
!!     ceil(num_work_procs/num_restart_procs) generating processes.
!!     Three integer arrays are added with the global indices for each
!!     cell/vert/edge.
!!
!!     Note: CDI does not allow writing of integer data, so those
!!     index arrays are written as double precision floats.
!!
!!     TODO: Add the NetCDF attribute 'icon-restart-uuid' to all files
!!     (same value everywhere), and sanity check their values when
!!     loading the restart.
!!
!! An alternative approach would be to include a tree of hashes with
!! its root in restartFile.mfr/attributes.nc, however this would add
!! considerable complexity to both the restart writing and reading, so
!! the simpler UUID approach should be used.
!!
!! A symbolic link to the latest successfully written restart file (=
!! the *.mfr directory) is added/updated.
!!
!!
!! Restart writing procedure:
!! ==========================
!!
!!  1. The current patch descriptions are transferred to the restart
!!     processes.  This is a two step process where the first step is
!!     to aggregate full up-to-date information on the work master
!!     (processor splitting!), the second step broadcasts the
!!     resulting descriptions to the restart processes.
!!
!!  2. The restart master process creates the restartFile.mfr/
!!     directory and ensures that it's empty.  If the directory
!!     contains only parts of an old multifile restart, we silently
!!     delete the contents, providing a warning if the multifile was
!!     incomplete.  If any unknown file is encountered in the process,
!!     we immediately fail with a hard error to avoid clobbering data.
!!
!!  3. The restart master process writes the restartFile.mfr/attributes.nc .
!!
!!  4. The restart processes write the
!!     restartFile.mfr/patch<JG>_metadata files.  This happens in
!!     round-robin fashion, restart process 1 writes
!!     restartFile.mfr/patch1_metadata, process 2 writes
!!     restartFile.mfr/patch2_metadata, and so on.  Since process 0 is
!!     generally not a part of this, this happens in parallel to the
!!     writing of the restartFile.mfr/attributes.nc file.
!!
!!  5. For each patch in order:
!!
!!      1. The restart processes write the global indices of the
!!         cells/verts/edges they handle to the the
!!         restartFile.mfr/patch<JG>_<N> files.
!!
!!      2. Data is transferred from the work processes to their
!!         respective restart process via direct one-sided MPI
!!         communication and written to the
!!         restartFile.mfr/patch<JG>_<N> files.
!!         See the notes on MPI communication below for details.
!!
!!
!! Restart reading procedure:
!! ==========================
!!
!! There may be a different count of reading work processes than there
!! were writing restart processes, so each reading work process may
!! need to handle the files produced by 0 to n writing restart
!! processes.
!!
!!  1. The restart master process reads and broadcasts the contents of
!!     restartFile.mfr/attributes.nc, the
!!     restartFile.mfr/patch<JG>_metadata files are ignored.  This
!!     step actually reuses the single-file restart loading code.
!!
!!  2. Each writing restart process is associated with exactly one
!!     reading work process in round robin fashion.  That is, if there
!!     were N writing restart processes and there are M reading work
!!     processes, then each reading work process will handle the data
!!     from floor(N/M) to ceil(N/M) writing restart processes.
!!
!!     Note that the reading is always done by the work processes only.
!!
!!  3. For each patch in order:
!!
!!      1. Each reading work process reads the integer arrays
!!         containing the global indices of the points in its files.
!!
!!      2. A t_comm_pattern is created that redistributes the data
!!         according to the global indices provided from the restart
!!         files and the global indices required by each work PE.
!!         This is implemented via the t_BrokerCommunicationPattern
!!         which uses an implicit naive domain decomposition to
!!         mediate the necessary information from the reading
!!         processes to the work processes that need the data.  This
!!         step shows reasonable performance on the Cray at the DWD.
!!
!!      3. The restart processes read their files and forward the
!!         contained data via the t_comm_pattern to their respective work
!!         processes.

!!         XXX/TODO: Measurements on the Cray at the DWD show that
!!         this is reasonably efficient when restarting with the same
!!         process count from a restart written in the joint-procs
!!         mode - in this case each process reads pretty much the same
!!         data as it needs. However, this step is quite expensive if
!!         significant amounts of data need to be transferred between
!!         the processes, making the advantage of the parallel reading
!!         disappear on the Cray. I am not sure how this could be
!!         solved, maybe YAXT performs better in this setting than
!!         t_comm_pattern, it might be worthwhile to try. But since
!!         I'm writing this during the last week of my employment at
!!         the DWD, I definitely lack the time to do it.
!!
!!
!! Notes on the MPI data transfer:
!! ===============================
!!
!! SENDER SIDE:
!!
!! All worker processes (i.e. the data senders) allocate local buffers
!! for every restart variable and every level. Thus, by a call of the
!! subroutine "multifilePatchData_exposeData", the restart data can be
!! completely copied from the prognostic fields to these buffers and -
!! in "dedicated proc mode" - the worker PEs can immediately proceed
!! with the next time step.
!!
!! Note: The send buffer type "t_CollectorSendBuffer" is identified
!! with an MPI memory window. When accessing the buffer locally for a
!! specific variable and level, we use objects of type
!! "t_MultifileRestartCollector" which point into this buffer.
!!
!! RECEIVER SIDE:
!!
!! On the writer PEs a receive buffer is allocated
!! ("t_CollectorRecvBuffer"). The size of this buffer corresponds to a
!! single vertical level and a horizontal chunk of the grid, see the
!! notes on the generated files above.  Such a buffer exists for
!! cell-based, edge-based and vertex-based variables each. When
!! calling "multifilePatchData_collectData", then a one-sided MPI_GET
!! is launched which fetches the data from the remote send buffers.
!!
!! The writer PEs then write the received level to the file via CDI
!! calls. In "dedicated proc mode" this does not halt the other PEs,
!! these may proceed with their work. In "joint proc mode", workers
!! and restart writers are the same PEs, which means that the process
!! is blocked until the data has been received and written to file.
!!
!!
!! Further notes:
!! ==============
!!
!! I have been warned against using the INQUIRE intrinsic for
!! anything, it is not reliable.  As such, I have opted to implement
!! the directory scanning code directly in C
!! (support/util_multifile_restart.c).
!!
!!

MODULE mo_multifile_restart
  USE mo_async_restart_packer,         ONLY: restartBcastRoot
  USE mo_c_restart_util,               ONLY: createEmptyMultifileDir
  USE mo_cdi,                          ONLY: streamOpenWrite, vlistCreate, streamDefVlist, streamClose, &
    &                                        vlistDestroy, FILETYPE_NC4, streamDefRecord,               &
    &                                        gridCreate, GRID_GENERIC, zaxisCreate, ZAXIS_GENERIC,      &
    &                                        vlistDefVar, TSTEP_CONSTANT, gridDestroy, zaxisDestroy
  USE mo_cdi_ids,                      ONLY: t_CdiIds
  USE mo_exception,                    ONLY: finish, message, message_text
  USE mo_grid_config,                  ONLY: n_dom
  USE mo_impl_constants,               ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_io_config,                    ONLY: restartWritingParameters
  USE mo_kind,                         ONLY: dp, i8
  USE mtime,                           ONLY: datetime
  USE mo_mpi,                          ONLY: p_bcast, my_process_is_work, my_process_is_restart,        &
    &                                        p_comm_work_2_restart, p_comm_work,                        &
    &                                        p_mpi_wtime, p_comm_work_restart, num_work_procs,          &
    &                                        my_process_is_mpi_workroot, p_reduce, p_sum_op
  USE mo_multifile_restart_patch_data, ONLY: t_MultifilePatchData, toMultifilePatchData
  USE mo_multifile_restart_util,       ONLY: createMultifileRestartLink, multifileAttributesPath,       &
    &                                        isAsync, rBuddy, rGroup,            &
    &                                        iAmRestartMaster, iAmRestartWriter, restartProcCount
  USE mo_packed_message,               ONLY: t_PackedMessage, kPackOp, kUnpackOp
  USE mo_restart_attributes,           ONLY: t_RestartAttributeList, RestartAttributeList_make
  USE mo_restart_descriptor,           ONLY: t_RestartDescriptor
  USE mo_restart_namelist,             ONLY: t_NamelistArchive, namelistArchive
  USE mo_restart_util,                 ONLY: t_restart_args, getRestartFilename
  USE mo_timer,                        ONLY: timer_start, timer_stop, timer_write_restart,               &
    &                                        timer_write_restart_communication, timer_write_restart_io,  &
    &                                        timers_level, timer_write_restart_setup,                    &
    &                                        timer_write_restart_wait
  USE mo_util_cdi,                     ONLY: cdiGetStringError
  USE mo_util_string,                  ONLY: int2string, real2string

  IMPLICIT NONE

  PUBLIC :: t_MultifileRestartDescriptor
  PUBLIC :: multifileRestart_mainLoop

  PRIVATE

  TYPE, EXTENDS(t_RestartDescriptor) :: t_MultifileRestartDescriptor
  CONTAINS
    PROCEDURE :: construct => multifileRestartDescriptor_construct
    PROCEDURE :: writeRestart => multifileRestartDescriptor_writeRestart
    PROCEDURE :: destruct => multifileRestartDescriptor_destruct
    PROCEDURE, PRIVATE :: updatePatchData => multifileRestartDescriptor_updatePatchData
    PROCEDURE, PRIVATE :: writeRestartInternal => multifileRestartDescriptor_writeRestartInternal
    PROCEDURE, PRIVATE :: writeAttributeFile => multifileRestartDescriptor_writeAttributeFile
  END TYPE t_MultifileRestartDescriptor

  !Constants for communicating the current operation to the restart processes.
  ENUM, BIND(C)
    ENUMERATOR :: kIllegalOp = 1, kWriteRestartOp, kShutdownOp
  END ENUM
  CHARACTER(*), PARAMETER :: modname = "mo_multifile_restart"

CONTAINS

  !Multifile restart initialization. This IS called at the very
  !beginning via mo_restart:createRestartDescriptor().
  !
  !This initializes the restart writers, assigning their roles IN
  !the joint procs case, AND creates the restart patch DATA
  !structures.
  SUBROUTINE multifileRestartDescriptor_construct(me, modelType)
    CLASS(t_MultifileRestartDescriptor), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: modelType
    CHARACTER(*), PARAMETER             :: routine = modname//":multifileRestartDescriptor_construct"
    LOGICAL                             :: lthis_pe_active
    INTEGER, ALLOCATABLE                :: srcRanks(:)
    TYPE(t_MultifilePatchData), POINTER :: patchData
    INTEGER                             :: error, jg, myProcId, myWrt, i, sRStrt, sREnd, jg0, &
      &                                    jfile, nStreams, nSrcRanks, ierr, ckSize

    IF(.NOT.my_process_is_work() .AND. .NOT.my_process_is_restart()) RETURN
    IF(timers_level >= 5) CALL timer_start(timer_write_restart)
    CALL me%restartDescriptor_construct(modelType)
    CALL restartWritingParameters(opt_nrestart_streams = nStreams)
    IF(isAsync()) CALL me%transferGlobalParameters()
    !set some local variables describing the communication that needs to be done
    CALL MPI_Comm_rank(p_comm_work_restart, myProcId, ierr)!restartWorkProcId()
    myWrt  = rBuddy() 
    nSrcRanks = COUNT((/(rBuddy(pe_in=i) == myWrt, i = 0, num_work_procs -1)/))
    ALLOCATE(srcRanks(nSrcRanks))
    srcRanks = PACK((/(i, i=0,num_work_procs-1)/), &
                     (/(rBuddy(pe_in=i)==myWrt, i=0,num_work_procs-1)/))
    IF (nStreams > nSrcRanks) &
      CALL finish(routine, "more horizontal multifile chunks worker ranks requested!")
    ! allocate patch data structure
    ALLOCATE(t_MultifilePatchData :: me%patchData(n_dom*nStreams), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    ckSize = (nSrcRanks + nStreams - 1)/nStreams
    DO jfile = 1, nStreams
      sRStrt = (jfile-1)*ckSize + 1
      sREnd  = MIN(nSrcRanks, jfile*ckSize)
      ! restart writer #writerRank: putting ranks
      ! sourceRank_start:sourceRank_end into a separate file.
      jg0 = (jfile-1)*n_dom
      DO jg = 1, n_dom
        patchData => toMultifilePatchData(me%patchData(jg0+jg))
        CALL patchData%construct(me%modelType, jg)
        ! initialize the patch DATA structures
        IF (.NOT.iAmRestartWriter()) THEN 
          lthis_pe_active = ANY(myProcId == srcRanks(sRStrt:sREnd))
          CALL patchData%createCollectors(myWrt, srcRanks(1:0), lthis_pe_active)
        ELSE
          CALL patchData%createCollectors(myWrt, srcRanks(sRStrt:sREnd), .TRUE.)
        END IF
      END DO
    END DO
    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE multifileRestartDescriptor_construct

  !This creates the restartArgs object, communicating it to the dedicated restart processes as appropriate.
  !On the restart side, this IS paired with a CALL to createRestartArgs_restart().
  FUNCTION createRestartArgs_compute(modelType, this_datetime, jstep, opt_output_jfile) &
    & RESULT(resultVar)
    TYPE(t_restart_args) :: resultVar
    CHARACTER(*), INTENT(IN) :: modelType
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
    INTEGER, INTENT(IN) :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)
    TYPE(t_PackedMessage) :: packedMessage
    CHARACTER(*), PARAMETER :: routine = modname//":createRestartArgs_compute"
    REAL(KIND=dp) :: time_start, time_end

    CALL resultVar%construct(this_datetime, jstep, modelType, opt_output_jfile)
    !In the CASE of dedicated proc mode, we need to inform the restart processes.
    IF(.NOT.isAsync()) RETURN
    IF(timers_level >= 7) CALL timer_start(timer_write_restart_wait)
    CALL packedMessage%construct()
    CALL packedMessage%pack(kWriteRestartOp)
    CALL resultVar%packer(kPackOp, packedMessage)
    IF(my_process_is_mpi_workroot()) time_start = p_mpi_wtime()
    CALL packedMessage%bcast(restartBcastRoot(), p_comm_work_2_restart)
    IF(my_process_is_mpi_workroot()) time_end = p_mpi_wtime()
    CALL packedMessage%destruct()
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
    IF(my_process_is_mpi_workroot()) THEN
      WRITE(message_text, "(a,e9.2,a)") 'compute procs waited ', time_end-time_start, &
                                        ' sec for dedicated restart procs to be ready...'
      CALL message(routine, message_text)             
    END IF
  END FUNCTION createRestartArgs_compute

  !This sends the stop message to the dedicated restart processes.
  !Like createRestartArgs_compute(), this also pairs with a CALL to createRestartArgs_restart().
  SUBROUTINE sendStopToRestart()
    TYPE(t_PackedMessage) :: packedMessage
    CHARACTER(*), PARAMETER :: routine = modname//":sendStopToRestart"
    REAL(KIND=dp) :: time_start, time_end

    IF(.NOT.isAsync()) RETURN
    IF(timers_level >= 7) CALL timer_start(timer_write_restart_wait)
    CALL packedMessage%construct()
    CALL packedMessage%pack(kShutdownOp)
    IF(my_process_is_mpi_workroot()) time_start = p_mpi_wtime()
    CALL packedMessage%bcast(restartBcastRoot(), p_comm_work_2_restart)
    IF(my_process_is_mpi_workroot()) time_end = p_mpi_wtime()
    CALL packedMessage%destruct()
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
    IF(my_process_is_mpi_workroot()) THEN
      WRITE(message_text, "(a,e9.2,a)") 'compute procs waited ', time_end-time_start, &
                                        ' sec for dedicated restart procs to finish...'
      CALL message(routine, message_text)
    END IF
  END SUBROUTINE sendStopToRestart

  ! This creates the restartArgs object by receiving it from the
  ! compute processes.
  !
  ! This pairs with either a CALL to createRestartArgs_compute() OR
  ! sendStopToRestart() on the work processes' side, which CALL it
  ! was paired with IS returned IN out_restartOp.
  !
  ! Important: In the CASE of a kShutdownOp, the RESULT IS NOT
  ! constructed, so out_restartOp must be analyzed before the
  ! FUNCTION RESULT IS used!
  FUNCTION createRestartArgs_restart(out_restartOp) RESULT(resultVar)
    TYPE(t_restart_args) :: resultVar
    INTEGER, INTENT(OUT) :: out_restartOp
    TYPE(t_PackedMessage) :: packedMessage
    CHARACTER(*), PARAMETER :: routine = modname//":createRestartArgs_restart"
    REAL(KIND=dp) :: time_start, time_end

    !receive the message from the work processes
    IF(timers_level >= 7) CALL timer_start(timer_write_restart_wait)
    CALL packedMessage%construct()
    IF(iAmRestartMaster()) time_start = p_mpi_wtime()
    CALL packedMessage%bcast(restartBcastRoot(), p_comm_work_2_restart)
    IF(iAmRestartMaster()) time_end = p_mpi_wtime()
    !unpack it
    out_restartOp = kIllegalOp
    CALL packedMessage%unpack(out_restartOp)
    IF(out_restartOp == kWriteRestartOp) CALL resultVar%packer(kUnpackOp, packedMessage)
    CALL packedMessage%destruct()
    IF(iAmRestartMaster()) THEN
      WRITE(message_text, "(a,e9.2,a)") 'dedicated restart procs waited ', &
                                       time_end-time_start, &
                                     ' sec for next event triggered by computes...'
      CALL message(routine, message_text)
    END IF
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
  END FUNCTION createRestartArgs_restart

  !The entry point to restart writing for the work processes.
  SUBROUTINE multifileRestartDescriptor_writeRestart(me, this_datetime, jstep, opt_output_jfile)
    CLASS(t_MultifileRestartDescriptor), INTENT(INOUT) :: me
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
    INTEGER, INTENT(IN) :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)
    TYPE(t_restart_args) :: restartArgs

    IF(timers_level >= 5) CALL timer_start(timer_write_restart)
    restartArgs = createRestartArgs_compute(me%modelType, this_datetime, jstep, opt_output_jfile)
    !This CALL IS entered by work AND restart processes IN lockstep:
    CALL me%writeRestartInternal(restartArgs)
    CALL restartArgs%destruct()
    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE multifileRestartDescriptor_writeRestart

  ! Ensure that all processes have full up-to-date knowledge of all patches.
  ! 
  ! This IS done by first letting the respective subset masters send
  ! their current state to the work master, AND THEN broadcasting
  ! that information to both work AND restart processes.
  !
  ! Collective across work AND restart PEs.
  SUBROUTINE multifileRestartDescriptor_updatePatchData(me)
    CLASS(t_MultifileRestartDescriptor), INTENT(INOUT) :: me
    INTEGER :: i
    TYPE(t_PackedMessage) :: packedMessage

    IF(timers_level >= 7) CALL timer_start(timer_write_restart_setup)
    !Make sure that the master has full up-to-date knowledge of all patches.
    IF (.NOT.isAsync() .OR. .NOT. my_process_is_restart()) THEN
      DO i = 1, SIZE(me%patchData)
        CALL me%patchData(i)%description%updateOnMaster()
      END DO
    END IF
    !Create a packedMessage with the full description on the master.
    CALL packedMessage%construct()
    IF(my_process_is_mpi_workroot()) THEN
       DO i = 1, SIZE(me%patchData)
         CALL me%patchData(i)%description%packer(kPackOp, packedMessage)
       END DO
    END IF
    !Distribute among workers.
    IF(my_process_is_work()) CALL packedMessage%bcast(0, p_comm_work)
    !Distribute to the dedicated restart processes (IF we have them).
    IF(isAsync()) CALL packedMessage%bcast(restartBcastRoot(), p_comm_work_2_restart)
    !Unpack everywhere.
    DO i = 1, SIZE(me%patchData)
       CALL me%patchData(i)%description%packer(kUnpackOp, packedMessage)
    END DO
    CALL packedMessage%destruct()
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_setup)
    !Update the vgrids.
    DO i = 1, SIZE(me%patchData)
      CALL me%patchData(i)%description%updateVGrids()
    END DO
  END SUBROUTINE multifileRestartDescriptor_updatePatchData

  ! This method IS entered by all work AND restart processes IN
  ! lockstep.
  !
  ! It first transfers all the DATA synchronously to the respective
  ! restart processes, THEN creates the directory AND metadata
  ! files, AND finally writes the payload DATA.
  SUBROUTINE multifileRestartDescriptor_writeRestartInternal(me, restartArgs)
    CLASS(t_MultifileRestartDescriptor), INTENT(INOUT) :: me
    TYPE(t_restart_args), INTENT(IN) :: restartArgs
    CHARACTER(*), PARAMETER               :: routine = ":multifileRestartDescriptor_writeRestartInternal"
    INTEGER                               :: jg, dummy, n_dom_active, nrestart_streams, findex
    INTEGER(i8)                           :: totBWritten, bWritten
    REAL(dp)                              :: gbWritten, dpTime
    CHARACTER(:), ALLOCATABLE             :: filename
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    TYPE(t_NamelistArchive), POINTER      :: namelists
    TYPE(t_CdiIds)                        :: cdiIds(SIZE(me%patchData))
    TYPE(t_MultifilePatchData), POINTER   :: patchData

    namelists => namelistArchive()
    dpTime = p_mpi_wtime()
    bWritten = 0_i8
    CALL restartWritingParameters(opt_nrestart_streams=nrestart_streams)
    !ensure that all processes have up-to-date patch DATA
    CALL me%updatePatchData()
    !create the multifile directory
    CALL getRestartFilename('multifile', 0, restartArgs, filename)
    IF(iAmRestartMaster()) THEN
        IF(createEmptyMultifileDir(filename) /= SUCCESS) CALL finish(routine, "error creating multifile-dir")
    END IF
    IF(timers_level >= 7) CALL timer_start(timer_write_restart_wait)
    !ensure that the other processes DO NOT continue before the
    !directory has been created:
    CALL p_bcast(dummy, rBuddy(0), p_comm_work_restart) 
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
    IF(iAmRestartMaster()) THEN
      restartAttributes => RestartAttributeList_make()
      CALL restartAttributes%setInteger('multifile_file_count', nrestart_streams*restartProcCount())
      n_dom_active = 0
      DO jg = 1, SIZE(me%patchData)
        patchData => toMultifilePatchData(me%patchData(jg))
        IF (patchData%description%l_dom_active) THEN
          n_dom_active = n_dom_active + 1
        END IF
      END DO
      CALL restartAttributes%setInteger('multifile_n_dom_active', n_dom_active)
      CALL me%defineRestartAttributes(restartAttributes, restartArgs)
      CALL me%writeAttributeFile(filename, restartAttributes, namelists)
      CALL restartAttributes%destruct()
      DEALLOCATE(restartAttributes)
    END IF
    ! COLLECT the payload data
    DO jg = 1, SIZE(me%patchData)
      patchData => toMultifilePatchData(me%patchData(jg))
      IF (patchData%description%l_dom_active) THEN
        CALL patchData%start_local_access()
      END IF
    END DO
    IF (iAmRestartWriter()) THEN
      DO jg = 1, SIZE(me%patchData)
        patchData => toMultifilePatchData(me%patchData(jg))
        IF (patchData%description%l_dom_active .AND. ASSOCIATED(patchData%varData)) THEN
          findex = nrestart_streams*rGroup() + (jg-1)/n_dom
          cdiIds(jg) = patchData%openPayloadFile(filename, findex, &
                                                 restartArgs%restart_datetime, bWritten)
        END IF
      END DO
    END IF
    IF (my_process_is_work()) THEN
      DO jg = 1, SIZE(me%patchData)
        patchData => toMultifilePatchData(me%patchData(jg))
        IF (patchData%description%l_dom_active) THEN
          CALL patchData%exposeData()
        END IF
      END DO
    END IF
    ! WRITE the payload files
    DO jg = 1, SIZE(me%patchData)
      patchData => toMultifilePatchData(me%patchData(jg))
      IF (patchData%description%l_dom_active) THEN
        CALL patchData%start_remote_access()
      END IF
    END DO
    IF(iAmRestartWriter()) THEN
      DO jg = 1, SIZE(me%patchData)
        patchData => toMultifilePatchData(me%patchData(jg))
        IF (patchData%description%l_dom_active) THEN
          CALL patchData%collectData(cdiIds(jg)%fHndl, bWritten)
          CALL cdiIds(jg)%closeAndDestroyIds()
        END IF
      END DO
    END IF
    IF(iAmRestartMaster()) CALL createMultifileRestartLink(filename, TRIM(restartArgs%modelType))
#ifndef NOMPI
    IF(isAsync().AND.my_process_is_work()) THEN
      !dedicated proc mode: work processes
      dpTime = p_mpi_wtime() - dpTime
      IF(my_process_is_mpi_workroot()) THEN
        WRITE(message_text, *) "restart: preparing checkpoint-data took "//TRIM(real2string(dpTime))//"s"
        CALL message(routine, message_text)
      END IF
    ELSE
      IF(timers_level >= 7) CALL timer_start(timer_write_restart_wait)
      IF(my_process_is_restart()) THEN
        !dedicated proc mode: restart processes
        totBWritten = p_reduce(bWritten, p_sum_op(), 0, p_comm_work)
      ELSE
        !joint proc mode: all processes
        totBWritten = p_reduce(bWritten, p_sum_op(), 0, p_comm_work_restart)
      END IF
      IF(timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
      !the previous CALL has synchronized us, so now it's safe
      !to stop the timer
      dpTime = p_mpi_wtime() - dpTime
      IF(iAmRestartMaster()) THEN
        gbWritten = REAL(totBWritten, dp)/1024.0_dp/1024.0_dp/1024.0_dp
        WRITE(message_text,*) "restart: finished: "//TRIM(real2string(gbWritten))//"GB of data in "//&
                   &TRIM(real2string(dpTime))//"s ("//TRIM(real2string(gbWritten/dpTime))//"GB/s)"
        CALL message(routine, message_text, all_print=.true.)
      END IF
    END IF
#endif
  END SUBROUTINE multifileRestartDescriptor_writeRestartInternal

  SUBROUTINE multifileRestartDescriptor_writeAttributeFile(me, filename, restartAttributes, namelists)
    CLASS(t_MultifileRestartDescriptor), INTENT(IN) :: me
    CHARACTER(*), INTENT(IN) :: filename
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    TYPE(t_NamelistArchive), INTENT(INOUT) :: namelists
    INTEGER :: metaFile, grid, zaxis, variable, vlist
    CHARACTER(:), ALLOCATABLE :: effectiveFilename
    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartDescriptor_writeAttributeFile"

    IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)
    CALL multifileAttributesPath(filename, effectiveFilename)
    metaFile = streamOpenWrite(effectiveFilename, FILETYPE_NC4)
    CALL checkCdiId(metaFile, "error opening file '"//effectiveFilename//"' for writing")
    vlist = vlistCreate()
    CALL checkCdiId(vlist, "error creating CDI vlist")
    ! to keep CDI happy ...
    grid = gridCreate(GRID_GENERIC, 1)
    CALL checkCdiId(grid, "error creating CDI grid")
    zaxis = zaxisCreate(ZAXIS_GENERIC, 1)
    CALL checkCdiId(zaxis, "error creating CDI zaxis")
    variable = vlistDefVar(vlist, grid, zaxis, TSTEP_CONSTANT)
    CALL checkCdiId(variable, "error creating CDI variable")
    CALL namelists%writeToCdiVlist(vlist)
    CALL restartAttributes%writeToCdiVlist(vlist)
    CALL streamDefVlist(metaFile, vlist)
    CALL streamDefRecord(metaFile, variable, 0)
    CALL streamClose(metaFile)
    CALL vlistDestroy(vlist)
    CALL gridDestroy(grid)
    CALL zaxisDestroy(zaxis)
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
  CONTAINS

   SUBROUTINE checkCdiId(cdiId, errorMessage)
     INTEGER, VALUE :: cdiId
     CHARACTER(*), INTENT(IN) :: errorMessage
     CHARACTER(MAX_CHAR_LENGTH) :: cdiErrorText

     IF(cdiId < 0) THEN
       CALL cdiGetStringError(cdiId, cdiErrorText)
       CALL finish(routine, errorMessage//": "//TRIM(cdiErrorText))
     END IF
   END SUBROUTINE checkCdiId
  END SUBROUTINE multifileRestartDescriptor_writeAttributeFile

  !Shuts down the dedicated restart writers IF they exist. And does some cleanup.
  SUBROUTINE multifileRestartDescriptor_destruct(me)
      CLASS(t_MultifileRestartDescriptor), INTENT(INOUT) :: me
      INTEGER :: i
      CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartDescriptor_destruct"

      IF(timers_level >= 5) CALL timer_start(timer_write_restart)
      IF(my_process_is_work()) CALL sendStopToRestart()
      DO i = 1, SIZE(me%patchData)
        CALL me%patchData(i)%destruct()
      END DO
      DEALLOCATE(me%patchData)
      IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE multifileRestartDescriptor_destruct

  ! This IS ONLY called on dedicated restart writers. It's an event
  ! loop that's driven by the messages from the compute processes.
  SUBROUTINE multifileRestart_mainLoop()
    TYPE(t_MultifileRestartDescriptor) :: restartDescriptor
    INTEGER :: restartOp
    TYPE(t_restart_args) :: restartArgs
    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestart_mainLoop"
    LOGICAL :: atService

    IF(.NOT.my_process_is_restart()) CALL finish(routine, "should not be here")
    CALL restartDescriptor%construct('')    !the model TYPE will be supplied by the work processes
    atService = .true.
    DO WHILE(atService)
      !get the next restart command
      restartOp = kIllegalOp
      restartArgs = createRestartArgs_restart(restartOp)
      SELECT CASE(restartOp)
        CASE(kWriteRestartOp)
          IF(timers_level >= 5) CALL timer_start(timer_write_restart)
          CALL restartDescriptor%writeRestartInternal(restartArgs)
          IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
        CASE(kShutdownOp)
          CALL restartDescriptor%destruct()
          atService = .false.
        CASE DEFAULT
          CALL finish(routine, "illegal restart operation received from work processes (the received &
                               &value is "//TRIM(int2string(restartOp))//")")
      END SELECT
    END DO
  END SUBROUTINE multifileRestart_mainLoop

END MODULE mo_multifile_restart
