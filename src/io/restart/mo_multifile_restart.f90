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
MODULE mo_multifile_restart
  USE mo_c_restart_util,               ONLY: createEmptyMultifileDir
  USE mo_cdi,                          ONLY: streamOpenWrite, vlistCreate, streamDefVlist, streamClose, &
    &                                        vlistDestroy, FILETYPE_NC4, streamDefRecord, CDI_UNDEFID,  &
    &                                        gridCreate, GRID_GENERIC, zaxisCreate, ZAXIS_GENERIC,      &
    &                                        vlistDefVar, TSTEP_CONSTANT, gridDestroy, zaxisDestroy
  USE mo_cdi_ids,                      ONLY: t_CdiIds
  USE mo_exception,                    ONLY: finish, message
  USE mo_grid_config,                  ONLY: n_dom
  USE mo_impl_constants,               ONLY: SUCCESS
  USE mo_io_config,                    ONLY: restartWritingParameters
  USE mo_kind,                         ONLY: dp, i8
  USE mtime,                           ONLY: datetime
  USE mo_mpi,                          ONLY: my_process_is_work, my_process_is_restart,                 &
    &                                        p_comm_work_2_restart, p_comm_work, p_comm_rank,           &
    &                                        p_mpi_wtime, p_comm_work_restart, num_work_procs,          &
    &                                        my_process_is_mpi_workroot, p_reduce, mpi_sum, p_barrier
  USE mo_multifile_restart_patch_data, ONLY: t_MultifilePatchData, toMultifilePatchData
  USE mo_multifile_restart_util,       ONLY: createMultifileRestartLink, isAsync, rBuddy, rGroup, &
    &                                        initUtil, iAmRestartMaster, iAmRestartWriter, restartProcCount
  USE mo_packed_message,               ONLY: t_PackedMessage, kPackOp, kUnpackOp
  USE mo_key_value_store,              ONLY: t_key_value_store
  USE mo_restart_descriptor,           ONLY: t_RestartDescriptor
  USE mo_restart_util,                 ONLY: t_restart_args, getRestartFilename, restartBcastRoot
  USE mo_timer,                        ONLY: timer_start, timer_stop, timer_write_restart,               &
    &                                        timer_write_restart_communication, timer_write_restart_io,  &
    &                                        timers_level, timer_write_restart_setup,                    &
    &                                        timer_write_restart_wait
  USE mo_util_string,                  ONLY: int2string, real2string
  USE mo_restart_nml_and_att,          ONLY: restartAttributeList_write_to_cdi

  IMPLICIT NONE

  PUBLIC :: t_MultifileRestartDescriptor
  PUBLIC :: multifileRestart_mainLoop

  PRIVATE

  TYPE, EXTENDS(t_RestartDescriptor) :: t_MultifileRestartDescriptor
  CONTAINS
    PROCEDURE :: construct => multifileRestartDescriptor_construct
    PROCEDURE :: writeRestart => multifileRestartDescriptor_writeRestart
    PROCEDURE :: destruct => multifileRestartDescriptor_destruct
    PROCEDURE, PRIVATE :: writeRestartInternal => multifileRestartDescriptor_writeRestartInternal
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
    CLASS(t_MultifileRestartDescriptor), INTENT(INOUT), TARGET :: me
    CHARACTER(*), INTENT(IN) :: modelType
    CHARACTER(*), PARAMETER             :: routine = modname//":multifileRestartDescriptor_construct"
    LOGICAL                             :: lthis_pe_active
    INTEGER, ALLOCATABLE                :: srcRanks(:)
    TYPE(t_MultifilePatchData), POINTER :: patchData(:)
    INTEGER                             :: error, jg, myProcId, myWrt, i, sRStrt, sREnd, jg0, &
      &                                    jfile, nStreams, nSrcRanks, ckSize

    IF(.NOT.my_process_is_work() .AND. .NOT.my_process_is_restart()) RETURN
    IF(timers_level >= 5) CALL timer_start(timer_write_restart)
    CALL initUtil()
    me%modelType = modelType
    CALL restartWritingParameters(opt_nrestart_streams = nStreams)
    IF(isAsync()) CALL me%transferGlobalParameters()
    !set some local variables describing the communication that needs to be done
    myProcId = p_comm_rank(p_comm_work_restart)
    myWrt  = rBuddy() 
    nSrcRanks = COUNT((/(rBuddy(pe_in=i) == myWrt, i = 0, num_work_procs -1)/))
    ALLOCATE(srcRanks(nSrcRanks))
    srcRanks = PACK((/(i, i=0,num_work_procs-1)/), &
                     (/(rBuddy(pe_in=i)==myWrt, i=0,num_work_procs-1)/))
    IF (nStreams > nSrcRanks) &
      CALL finish(routine, "more horizontal multifile chunks than worker ranks requested!")
    ! allocate patch data structure
    ALLOCATE(t_MultifilePatchData :: me%patchData(n_dom*nStreams), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    ckSize = (nSrcRanks + nStreams - 1)/nStreams
    patchData => toMultifilePatchData(me%patchData)
    DO jfile = 1, nStreams
      sRStrt = (jfile-1)*ckSize + 1
      sREnd  = MIN(nSrcRanks, jfile*ckSize)
      ! restart writer #writerRank: putting ranks
      ! sourceRank_start:sourceRank_end into a separate file.
      jg0 = (jfile-1)*n_dom
      DO jg = 1, n_dom
        CALL patchData(jg0+jg)%construct(me%modelType, jg)
        ! initialize the patch DATA structures
        IF (.NOT.iAmRestartWriter()) THEN
          lthis_pe_active = ANY(myProcId == srcRanks(sRStrt:sREnd))
          CALL patchData(jg0+jg)%createCollectors(myWrt, srcRanks(1:0), lthis_pe_active)
          ! to keep mo_timer happy
          IF (timers_level >= 7) CALL timer_start(timer_write_restart_communication)
          IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
        ELSE
          CALL patchData(jg0+jg)%createCollectors(myWrt, srcRanks(sRStrt:sREnd), .TRUE.)
        END IF
      END DO
    END DO
    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE multifileRestartDescriptor_construct

  !The entry point to restart writing for the work processes.
  SUBROUTINE multifileRestartDescriptor_writeRestart(me, this_datetime, jstep, opt_output_jfile)
    CLASS(t_MultifileRestartDescriptor), INTENT(INOUT), TARGET :: me
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
    INTEGER, INTENT(IN) :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)
    TYPE(t_restart_args) :: rArgs

    IF(timers_level >= 5) CALL timer_start(timer_write_restart)
    CALL createRestartArgs_compute()
    CALL me%writeRestartInternal(rArgs)
    CALL rArgs%destruct()
    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  CONTAINS

    SUBROUTINE createRestartArgs_compute()
      TYPE(t_PackedMessage) :: pmsg
      CHARACTER(*), PARAMETER :: routine = modname//":createRestartArgs_compute"
      REAL(KIND=dp) :: timing

      CALL rArgs%construct(this_datetime, jstep, me%modelType, opt_output_jfile)
#ifndef NOMPI
      IF(isAsync()) THEN
        IF (timers_level >= 7) CALL timer_start(timer_write_restart_wait)
        CALL pmsg%pack(kWriteRestartOp)
        CALL rArgs%packer(kPackOp, pmsg)
        IF (my_process_is_mpi_workroot()) timing = p_mpi_wtime()
        CALL pmsg%bcast(restartBcastRoot(), p_comm_work_2_restart)
        IF (timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
        IF (my_process_is_mpi_workroot()) &
          & CALL message(routine, 'compute procs waited ' // &
            & TRIM(real2string(p_mpi_wtime() - timing))//' sec for dedicated restart procs to be ready...')
      END IF
#endif
    END SUBROUTINE createRestartArgs_compute
  END SUBROUTINE multifileRestartDescriptor_writeRestart

  SUBROUTINE multifileRestartDescriptor_writeRestartInternal(me, restartArgs)
    CLASS(t_MultifileRestartDescriptor), INTENT(INOUT), TARGET :: me
    TYPE(t_restart_args), INTENT(IN) :: restartArgs
    CHARACTER(*), PARAMETER               :: routine = ":writeRestartInternal"
    INTEGER                               :: jg, n_rstreams, findex
    INTEGER(i8)                           :: totBWritten, bWritten
    REAL(dp)                              :: gbWritten, dpTime
    CHARACTER(:), ALLOCATABLE             :: filename
    TYPE(t_CdiIds)                        :: cdiIds(SIZE(me%patchData))
    TYPE(t_MultifilePatchData), POINTER   :: patchData(:)

    dpTime = p_mpi_wtime()
    bWritten = 0_i8
    CALL restartWritingParameters(opt_nrestart_streams=n_rstreams)
    CALL updatePatchData()
    patchData => toMultifilePatchData(me%patchData)
    DO jg = 1, SIZE(me%patchData)
      IF (patchData(jg)%description%l_dom_active) &
        & CALL patchData(jg)%start_local_access()
    END DO
    IF (my_process_is_work()) THEN
      DO jg = 1, SIZE(me%patchData)
        IF (patchData(jg)%description%l_dom_active) THEN
          CALL patchData(jg)%exposeData()
        END IF
      END DO
    END IF
    DO jg = 1, SIZE(me%patchData)
      IF (patchData(jg)%description%l_dom_active) THEN
        CALL patchData(jg)%start_remote_access()
      END IF
    END DO
    IF(iAmRestartMaster()) THEN
      CALL getRestartFilename('multifile', 0, restartArgs, filename)
      IF (createEmptyMultifileDir(filename) /= SUCCESS) &
        & CALL finish(routine, "error creating multifile-dir")
      CALL p_barrier(p_comm_work)
      CALL writeAttributeFile()
    ELSE IF(.NOT.isAsync().OR.iAmRestartWriter()) THEN
      CALL p_barrier(p_comm_work)
    ELSE 
      dpTime = p_mpi_wtime() - dpTime
      IF(my_process_is_mpi_workroot()) &
        CALL message(routine, "restart: preparing checkpoint-data took " &
          & //TRIM(real2string(dpTime))//"s")
    END IF
    IF (iAmRestartWriter()) THEN
      CALL getRestartFilename('multifile', 0, restartArgs, filename)
      DO jg = 1, SIZE(me%patchData)
        IF (patchData(jg)%description%l_dom_active .AND. SIZE(patchData(jg)%varData) > 0) THEN
          findex = n_rstreams*rGroup() + (jg-1)/n_dom
          cdiIds(jg) = patchData(jg)%openPayloadFile(filename, findex, &
            & restartArgs%restart_datetime, bWritten)
          CALL patchData(jg)%collectData(cdiIds(jg)%fHndl, bWritten)
          CALL cdiIds(jg)%closeAndDestroyIds()
        END IF
      END DO
    END IF
    IF(iAmRestartMaster()) CALL createMultifileRestartLink(filename, TRIM(restartArgs%modelType))
#ifndef NOMPI
    IF(.NOT.(isAsync().AND.my_process_is_work())) THEN
      IF(timers_level >= 7) CALL timer_start(timer_write_restart_wait)
      IF(my_process_is_restart()) THEN
        !dedicated proc mode: restart processes
        totBWritten = p_reduce(bWritten, mpi_sum, 0, p_comm_work)
      ELSE
        !joint proc mode: all processes
        totBWritten = p_reduce(bWritten, mpi_sum, 0, p_comm_work_restart)
      END IF
      IF(timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
      dpTime = p_mpi_wtime() - dpTime
      IF(iAmRestartMaster()) THEN
        gbWritten = REAL(totBWritten, dp)/1024.0_dp/1024.0_dp/1024.0_dp
        CALL message(routine, "restart: finished: " // TRIM(real2string(gbWritten)) // &
          & "GB of data in " // TRIM(real2string(dpTime)) // "s (" // &
          & TRIM(real2string(gbWritten/dpTime)) // "GB/s)", all_print=.true.)
      END IF
    END IF
#endif
  CONTAINS

    SUBROUTINE writeAttributeFile()
      INTEGER :: metaFile, grid, zaxis, var, vlist, n_dom_active
      CHARACTER(:), ALLOCATABLE :: mafname
      CHARACTER(*), PARAMETER :: routine = modname//":writeAttributeFile"
      TYPE(t_key_value_store), ALLOCATABLE :: rAttribs

      CALL me%defineRestartAttributes(rAttribs, restartArgs)
      CALL rAttribs%put('multifile_file_count', n_rstreams*restartProcCount())
      n_dom_active = COUNT(patchData(:)%description%l_dom_active)
      CALL rAttribs%put('multifile_n_dom_active', n_dom_active)
      IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)
      mafname = filename//"/attributes.nc"
      metaFile = streamOpenWrite(mafname, FILETYPE_NC4)
      IF (metaFile .EQ. CDI_UNDEFID) &
        & CALL finish(routine, "error opening file '" // mafname // "'")
      vlist = vlistCreate() ! ... to keep CDI happy ...
      grid = gridCreate(GRID_GENERIC, 1)
      zaxis = zaxisCreate(ZAXIS_GENERIC, 1)
      var = vlistDefVar(vlist, grid, zaxis, TSTEP_CONSTANT)
      CALL restartAttributeList_write_to_cdi(rAttribs, vlist)
      CALL streamDefVlist(metaFile, vlist)
      CALL streamDefRecord(metaFile, var, 0)
      CALL streamClose(metaFile)
      CALL vlistDestroy(vlist)
      CALL gridDestroy(grid)
      CALL zaxisDestroy(zaxis)
      IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
      CALL rAttribs%destruct()
      DEALLOCATE(rAttribs)
    END SUBROUTINE writeAttributeFile

    SUBROUTINE updatePatchData()
      INTEGER :: i
      TYPE(t_PackedMessage) :: packedMessage

      IF(timers_level >= 7) CALL timer_start(timer_write_restart_setup)
      IF (.NOT.isAsync() .OR. .NOT. my_process_is_restart()) THEN
        DO i = 1, SIZE(me%patchData)
          CALL me%patchData(i)%description%updateOnMaster()
        END DO
      END IF
      IF(my_process_is_mpi_workroot()) THEN
         DO i = 1, SIZE(me%patchData)
           CALL me%patchData(i)%description%packer(kPackOp, packedMessage)
         END DO
      END IF
      IF(my_process_is_work()) CALL packedMessage%bcast(0, p_comm_work)
#ifndef NOMPI
      IF(isAsync()) CALL packedMessage%bcast(restartBcastRoot(), p_comm_work_2_restart)
#endif
      DO i = 1, SIZE(me%patchData)
         CALL me%patchData(i)%description%packer(kUnpackOp, packedMessage)
      END DO
      IF(timers_level >= 7) CALL timer_stop(timer_write_restart_setup)
      DO i = 1, SIZE(me%patchData)
        CALL me%patchData(i)%description%updateVGrids()
      END DO
    END SUBROUTINE updatePatchData
  END SUBROUTINE multifileRestartDescriptor_writeRestartInternal

  !Shuts down the dedicated restart writers IF they exist. And does some cleanup.
  SUBROUTINE multifileRestartDescriptor_destruct(me)
    CLASS(t_MultifileRestartDescriptor), INTENT(INOUT) :: me
    INTEGER :: i

    IF(timers_level >= 5) CALL timer_start(timer_write_restart)
    IF(my_process_is_work()) CALL sendStopToRestart()
    DO i = 1, SIZE(me%patchData)
      CALL me%patchData(i)%destruct()
    END DO
    DEALLOCATE(me%patchData)
    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  CONTAINS

    SUBROUTINE sendStopToRestart()
      TYPE(t_PackedMessage) :: pmsg
      CHARACTER(*), PARAMETER :: routine = modname//":sendStopToRestart"
      REAL(KIND=dp) :: timing

#ifndef NOMPI
      IF(.NOT.isAsync()) RETURN
      IF(timers_level >= 7) CALL timer_start(timer_write_restart_wait)
      CALL pmsg%pack(kShutdownOp)
      IF(my_process_is_mpi_workroot()) timing = p_mpi_wtime()
      CALL pmsg%bcast(restartBcastRoot(), p_comm_work_2_restart)
      IF(timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
      IF(my_process_is_mpi_workroot()) THEN
        timing = p_mpi_wtime() - timing
        CALL message(routine, 'compute procs waited ' // TRIM(real2string(timing)) // &
                              ' sec for dedicated restart procs to finish...')
      END IF
#endif
    END SUBROUTINE sendStopToRestart
  END SUBROUTINE multifileRestartDescriptor_destruct

  ! This IS ONLY called on dedicated restart writers. It's an event
  ! loop that's driven by the messages from the compute processes.
  SUBROUTINE multifileRestart_mainLoop()
#ifndef NOMPI
    TYPE(t_MultifileRestartDescriptor) :: restartDescriptor
    INTEGER :: op_code
    TYPE(t_restart_args) :: rArgs
    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestart_mainLoop"
    LOGICAL :: atService

    IF(.NOT.my_process_is_restart()) CALL finish(routine, "should not be here")
    CALL restartDescriptor%construct('')    !the model TYPE will be supplied by the work processes
    atService = .true.
    DO WHILE(atService)
      CALL createRestartArgs_restart()
      SELECT CASE(op_code)
      CASE(kWriteRestartOp)
        IF(timers_level >= 5) CALL timer_start(timer_write_restart)
        CALL restartDescriptor%writeRestartInternal(rArgs)
        IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
      CASE(kShutdownOp)
        atService = .false.
      CASE DEFAULT
        CALL finish(routine, "illegal restart operation received from work processes (the received &
                             &value is "//TRIM(int2string(op_code))//")")
      END SELECT
    END DO
    CALL restartDescriptor%destruct()
  CONTAINS

    SUBROUTINE createRestartArgs_restart()
      TYPE(t_PackedMessage) :: pmsg
      REAL(dp) :: timing
  
      IF(timers_level >= 7) CALL timer_start(timer_write_restart_wait)
      IF(iAmRestartMaster()) timing = p_mpi_wtime()
      CALL pmsg%bcast(restartBcastRoot(), p_comm_work_2_restart)
      op_code = kIllegalOp
      CALL pmsg%unpack(op_code)
      IF (op_code == kWriteRestartOp) CALL rArgs%packer(kUnpackOp, pmsg)
      IF(iAmRestartMaster()) THEN
        timing = p_mpi_wtime() - timing
        CALL message(routine, 'dedicated restart procs waited ' // &
          & TRIM(real2string(timing)) // ' s for next event triggered by computes...')
      END IF
      IF(timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
    END SUBROUTINE createRestartArgs_restart
#endif
  END SUBROUTINE multifileRestart_mainLoop

END MODULE mo_multifile_restart
