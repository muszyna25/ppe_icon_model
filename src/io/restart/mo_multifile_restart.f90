!> Module for writing restart as a multifile (synchronously)
!!
!! Note: Other implementations of the restart writing interface can be
!! found in mo_sync_restart and mo_async_restart.
!!
!! Initial implementation: Nathanael Huebbe
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
!!   * restartFile.mfr/patch<JG>_metadata: 
!!     Full description for each patch. This especially includes the
!!     info in t_RestartPatchDescription.
!!
!!     TODO: While this is currently written, it appears pointless as
!!     I have not had any need to read it back in. Consequently, this
!!     should be removed as it only complicates the code. Note that
!!     removing this also allows removing mo_serialized_data entirely.
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
!!         respective restart process via direct point to point
!!         communication and written to the
!!         restartFile.mfr/patch<JG>_<N> files.
!!
!!         In joint-proc mode, the data is written immediately after
!!         it has been gathered on the respective restart writer.
!!         This is to keep the memory footprint down.  (All processes
!!         are synchronizing with the restart writers anyway, so there
!!         is no room for functional parallelism in this case.)
!!
!!         Things are different with dedicated restart processes: In
!!         this case, all the restart data is first transfered to
!!         memory buffers in the restart writers.  After that, the
!!         work processes carry on with their work immediately, while
!!         the restart writers perform the actual restart writing
!!         asynchronously.
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
!! Further notes:
!! ==============
!!
!! I have been warned against using the INQUIRE intrinsic for
!! anything, it is not reliable.  As such, I have opted to implement
!! the directory scanning code directly in C
!! (support/util_multifile_restart.c).
!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Dedicated proc mode: The current implementation creates a
!! bottle neck: in "openPayloadFile", a loop over all variables
!! and all levels collects all restart data. The collecting
!! process therefore runs into memory problems, before it is
!! able to write the data out. Thus it would be necessary to
!! put the variabel loop outside of the the collect&write
!! process. This, again, is not possible, since the collect*
!! implementation uses blocking sends and receives!
!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!

MODULE mo_multifile_restart
    USE ISO_C_BINDING,                   ONLY: C_INT
    USE mo_async_restart_packer,         ONLY: restartBcastRoot
    USE mo_c_restart_util,               ONLY: createEmptyMultifileDir
    USE mo_cdi,                          ONLY: streamOpenWrite, vlistCreate, streamDefVlist, streamClose, &
      &                                        vlistDestroy, FILETYPE_NC4, streamDefRecord,               &
      &                                        gridCreate, GRID_GENERIC, zaxisCreate, ZAXIS_GENERIC,      &
      &                                        vlistDefVar, TSTEP_CONSTANT, gridDestroy, zaxisDestroy
    USE mo_cdi_ids,                      ONLY: t_CdiIds
    USE mo_exception,                    ONLY: finish
    USE mo_grid_config,                  ONLY: n_dom
    USE mo_impl_constants,               ONLY: SUCCESS, MAX_CHAR_LENGTH
    USE mo_io_config,                    ONLY: restartWritingParameters
    USE mo_kind,                         ONLY: dp, i8
    USE mtime,                           ONLY: datetime
    USE mo_mpi,                          ONLY: p_bcast, my_process_is_work, my_process_is_restart,        &
      &                                        p_comm_work_2_restart, p_pe_work, p_comm_work, p_barrier,  &
      &                                        p_mpi_wtime, p_comm_work_restart,                          &
      &                                        my_process_is_mpi_workroot, p_reduce, p_sum_op
    USE mo_multifile_restart_patch_data, ONLY: t_MultifilePatchData, toMultifilePatchData
    USE mo_multifile_restart_util,       ONLY: createMultifileRestartLink, multifileAttributesPath,       &
      &                                        multifileMetadataPath
    USE mo_packed_message,               ONLY: t_PackedMessage, kPackOp, kUnpackOp
    USE mo_restart_attributes,           ONLY: t_RestartAttributeList, RestartAttributeList_make
    USE mo_restart_descriptor,           ONLY: t_RestartDescriptor
    USE mo_restart_namelist,             ONLY: t_NamelistArchive, namelistArchive
    USE mo_restart_util,                 ONLY: t_restart_args, workProcCount, dedicatedRestartProcCount,   &
      &                                        restartProcCount, isDedicatedProcMode, restartWriterId,     &
      &                                        my_process_is_restart_master, my_process_is_restart_writer, &
      &                                        restartWorkProcId, restartWorkProcId2Rank, getRestartFilename
    USE mo_serialized_data,              ONLY: t_SerializedData
    USE mo_timer,                        ONLY: timer_start, timer_stop, timer_write_restart,               &
      &                                        timer_write_restart_communication, timer_write_restart_io,  &
      &                                        timers_level
    USE mo_util_cdi,                     ONLY: cdiGetStringError
    USE mo_util_file,                    ONLY: putFile
    USE mo_util_string,                  ONLY: int2string, real2string

    IMPLICIT NONE

    PUBLIC :: t_MultifileRestartDescriptor
    PUBLIC :: multifileRestart_mainLoop

    PRIVATE

    TYPE, EXTENDS(t_RestartDescriptor) :: t_MultifileRestartDescriptor
        !this IS for measuring the throughput
        REAL(dp) :: startTime
        INTEGER(i8) :: bytesWritten
    CONTAINS
      ! override, collective across both work AND restart procs
        PROCEDURE :: construct => multifileRestartDescriptor_construct

        PROCEDURE :: writeRestart => multifileRestartDescriptor_writeRestart    ! override

        PROCEDURE, PRIVATE :: patchDataPacker => multifileRestartDescriptor_patchDataPacker

        ! collective across both work AND restart procs
        PROCEDURE, PRIVATE :: updatePatchData => multifileRestartDescriptor_updatePatchData

        ! collective across both work AND restart procs
        PROCEDURE, PRIVATE :: writeRestartInternal => multifileRestartDescriptor_writeRestartInternal

        ! ONLY executed by the restart master
        PROCEDURE, PRIVATE :: writeAttributeFile => multifileRestartDescriptor_writeAttributeFile

        ! non-collective, executed once by a single PE for each patch
        PROCEDURE, PRIVATE :: writeMetadataFile => multifileRestartDescriptor_writeMetadataFile

        ! collective across both work AND restart procs
        PROCEDURE, PRIVATE :: openPayloadFile => multifileRestartDescriptor_openPayloadFile

        ! no communication
        PROCEDURE, PRIVATE :: closePayloadFile => multifileRestartDescriptor_closePayloadFile

        PROCEDURE :: destruct => multifileRestartDescriptor_destruct    ! override
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

        LOGICAL :: lDedicatedProcMode
        INTEGER :: error, jg, myProcId, writerCount, writerProcId, writerRank, i, sourceRankCount
        INTEGER, ALLOCATABLE :: sourceRanks(:)
        TYPE(t_MultifilePatchData), POINTER :: patchData
        CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartDescriptor_construct"

        IF(.NOT.my_process_is_work() .AND. .NOT.my_process_is_restart()) RETURN

        IF(timers_level >= 5) CALL timer_start(timer_write_restart)

        CALL me%restartDescriptor_construct(modelType)

        CALL restartWritingParameters(opt_lDedicatedProcMode = lDedicatedProcMode)
        IF(lDedicatedProcMode) CALL me%transferGlobalParameters()

        !set some local variables describing the communication that needs to be done
        myProcId = restartWorkProcId()
        writerCount = restartProcCount()
        !TODO: It might be better to use a different scheme here to
        !      get adjacent processes into the same restart writing
        !      group.  It would be ideal if we could guarantee that
        !      all compute nodes housed the same amount of writer
        !      processes.  In either case, restartWriteId(),
        !      restartWorkProcId(), restartWorkProcId2Rank(),
        !      my_process_is_restart_master(), and
        !      my_process_is_restart_writer() in mo_restart_util would
        !      need to be adapted accordingly.
        writerProcId = MODULO(myProcId, writerCount)
        writerRank = restartWorkProcId2Rank(writerProcId)
        sourceRanks = [(restartWorkProcId2Rank(i), i = writerProcId, workProcCount() + dedicatedRestartProcCount() - 1, &
                                                     & writerCount)]
        sourceRankCount = SIZE(sourceRanks)
        IF(.NOT.my_process_is_restart_writer()) sourceRankCount = 0 !no source processes IF this IS a NOT a writer PE

        ! ALLOCATE patch DATA structure
        ALLOCATE(t_MultifilePatchData :: me%patchData(n_dom), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        ! initialize the patch DATA structures
        DO jg = 1, n_dom
            patchData => toMultifilePatchData(me%patchData(jg))
            CALL patchData%construct(me%modelType, jg)
            CALL patchData%createCollectors(writerRank, sourceRanks(1:sourceRankCount))
        END DO

        IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
    END SUBROUTINE multifileRestartDescriptor_construct

    !This creates the restartArgs object, communicating it to the dedicated restart processes as appropriate.
    !On the restart side, this IS paired with a CALL to createRestartArgs_restart().
    FUNCTION createRestartArgs_compute(modelType, this_datetime, jstep, opt_output_jfile) RESULT(resultVar)
        TYPE(t_restart_args) :: resultVar
        CHARACTER(*), INTENT(IN) :: modelType
        TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
        INTEGER, INTENT(IN) :: jstep
        INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

        LOGICAL :: lDedicatedProcMode
        TYPE(t_PackedMessage) :: packedMessage
        CHARACTER(*), PARAMETER :: routine = modname//":createRestartArgs_compute"

        CALL resultVar%construct(this_datetime, jstep, modelType, opt_output_jfile)

        !In the CASE of dedicated proc mode, we need to inform the restart processes.
        CALL restartWritingParameters(opt_lDedicatedProcMode = lDedicatedProcMode)
        IF(lDedicatedProcMode) THEN
            IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
            CALL packedMessage%construct()
            CALL packedMessage%pack(kWriteRestartOp)
            CALL resultVar%packer(kPackOp, packedMessage)
            CALL packedMessage%bcast(restartBcastRoot(), p_comm_work_2_restart)
            CALL packedMessage%destruct()
            IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
        END IF
    END FUNCTION createRestartArgs_compute

    !This sends the stop message to the dedicated restart processes.
    !Like createRestartArgs_compute(), this also pairs with a CALL to createRestartArgs_restart().
    SUBROUTINE sendStopToRestart()
        LOGICAL :: lDedicatedProcMode
        TYPE(t_PackedMessage) :: packedMessage
        CHARACTER(*), PARAMETER :: routine = modname//":sendStopToRestart"

        CALL restartWritingParameters(opt_lDedicatedProcMode = lDedicatedProcMode)
        IF(lDedicatedProcMode) THEN
            IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
            CALL packedMessage%construct()
            CALL packedMessage%pack(kShutdownOp)
            CALL packedMessage%bcast(restartBcastRoot(), p_comm_work_2_restart)
            CALL packedMessage%destruct()
            IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
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

        IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)

        !receive the message from the work processes
        CALL packedMessage%construct()
        CALL packedMessage%bcast(restartBcastRoot(), p_comm_work_2_restart)

        !unpack it
        out_restartOp = kIllegalOp
        CALL packedMessage%unpack(out_restartOp)
        IF(out_restartOp == kWriteRestartOp) CALL resultVar%packer(kUnpackOp, packedMessage)
        CALL packedMessage%destruct()

        IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
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

    !Serialization of the patchData array.
    SUBROUTINE multifileRestartDescriptor_patchDataPacker(me, operation, packedMessage)
        CLASS(t_MultifileRestartDescriptor), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        TYPE(t_PackedMessage), INTENT(INOUT) :: packedMessage

        INTEGER :: i

        DO i = 1, SIZE(me%patchData, 1)
            CALL me%patchData(i)%description%packer(operation, packedMessage)
        END DO
    END SUBROUTINE multifileRestartDescriptor_patchDataPacker

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
        LOGICAL :: lDedicatedProcMode

        IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)

        CALL restartWritingParameters(opt_lDedicatedProcMode = lDedicatedProcMode)

        !Make sure that the master has full up-to-date knowledge of all patches.
        IF (.NOT. lDedicatedProcMode .OR. .NOT. my_process_is_restart()) THEN
          DO i = 1, SIZE(me%patchData, 1)
            CALL me%patchData(i)%description%updateOnMaster()
          END DO
        END IF

        !Create a packedMessage with the full description on the master.
        CALL packedMessage%construct()
        IF(p_pe_work == 0) CALL me%patchDataPacker(kPackOp, packedMessage)

        !Distribute among workers.
        IF(my_process_is_work()) CALL packedMessage%bcast(0, p_comm_work)

        !Distribute to the dedicated restart processes (IF we have them).
        IF(lDedicatedProcMode) CALL packedMessage%bcast(restartBcastRoot(), p_comm_work_2_restart)

        !Unpack everywhere.
        CALL me%patchDataPacker(kUnpackOp, packedMessage)
        CALL packedMessage%destruct()

        IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)

        !Update the vgrids.
        DO i = 1, SIZE(me%patchData, 1)
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

        INTEGER :: jg, myFirstPatch, dummy
        INTEGER(i8) :: totalBytesWritten
        REAL(dp) :: gibibytesWritten, elapsedTime
        CHARACTER(:), ALLOCATABLE :: filename
        TYPE(t_RestartAttributeList), POINTER :: restartAttributes
        TYPE(t_NamelistArchive), POINTER :: namelists
        TYPE(t_CdiIds) :: cdiIds(SIZE(me%patchData))
        CHARACTER(*), PARAMETER :: routine = ":multifileRestartDescriptor_writeRestartInternal"

        namelists => namelistArchive()

        !start the timer
        IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
        CALL p_barrier(p_comm_work_restart)
        IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
        me%startTime = p_mpi_wtime()
        me%bytesWritten = 0_i8

        !ensure that all processes have up-to-date patch DATA
        CALL me%updatePatchData()

        !create the multifile directory
        filename = getRestartFilename('multifile', 0, restartArgs)
        IF(my_process_is_restart_master()) THEN
            IF(createEmptyMultifileDir(filename) /= SUCCESS) CALL finish(routine, "error creating restart multifile")
        END IF
        IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
        !ensure that the other processes DO NOT continue before the
        !directory has been created:
        CALL p_bcast(dummy, restartWorkProcId2Rank(0), p_comm_work_restart) 
        IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)

        !WRITE the attributes.nc file
        IF(my_process_is_restart_master()) THEN
            restartAttributes => RestartAttributeList_make()
            !this IS for checking the completeness of the multifile
            !when loading the restart:
            CALL restartAttributes%setInteger('multifile_file_count', restartProcCount())
            CALL me%defineRestartAttributes(restartAttributes, restartArgs)

            CALL me%writeAttributeFile(filename, restartAttributes, namelists)

            CALL restartAttributes%destruct()
            DEALLOCATE(restartAttributes)
        END IF

        !WRITE the patchN_metadata files
        IF(my_process_is_restart_writer()) THEN
            ! This does NOT USE the straight-forward
            ! `restartWriterId() + 1` because that would cause the
            ! restart master to also be the first process to WRITE a
            ! _metadata file.  The way, we DO it here, puts the
            ! restart master at the END of the line of processes to
            ! WRITE a _metadata file, allowing all the _metadata files
            ! to be written IN parallel to the attributes.nc file IN
            ! the normal CASE where we have more restart processes
            ! than patches.
            myFirstPatch = restartWriterId()    ! 0 ... restartProcCount-1
            IF(myFirstPatch == 0) myFirstPatch = restartProcCount() ! 1 ... restartProcCount
            DO jg = myFirstPatch, SIZE(me%patchData), restartProcCount()
                CALL me%writeMetadataFile(filename, jg)
            END DO
        END IF

        !WRITE the payload files
        !
        !IN the CASE of dedicated proc mode, the first loop just
        !streams the DATA to the restart PEs.  The actual writing
        !happens IN the second loop which the work PEs can leave IN a
        !hurry as its communication free, WHILE the restart writers
        !take significantly longer to complete the close operation.
        !Accordingly, it IS important, NOT to mix the contents of the
        !two loop into one - that would force the workers to wait for
        !the actual writing to complete.
        !
        !In the CASE of joint proc mode, the writing IS performed
        !interleaved with the collecting of DATA inside the first
        !loop.
        !
        DO jg = 1, SIZE(me%patchData)
            cdiIds(jg) = me%openPayloadFile(filename, jg, restartArgs%restart_datetime)
        END DO
        DO jg = 1, SIZE(me%patchData)
            CALL me%closePayloadFile(jg, cdiIds(jg))
        END DO

        IF(my_process_is_restart_master()) CALL createMultifileRestartLink(filename, TRIM(restartArgs%modelType))

        !stop the timer AND print the resulting speed
        !
        !While this code could be compiled AND run IN the NOMPI CASE,
        !it would print nonsense since p_mpi_wtime() always returns
        !zero IN this CASE.
#ifndef NOMPI
        IF(isDedicatedProcMode().AND.my_process_is_work()) THEN
          !IN dedicated proc mode we must avoid blocking the work
          !processes

            !dedicated proc mode: work processes
            IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
            CALL p_barrier(p_comm_work)
            IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
            !the previous CALL has synchronized us, so now it's safe
            !to stop the timer
            elapsedTime = p_mpi_wtime() - me%startTime
            IF(my_process_is_mpi_workroot()) THEN
                WRITE(0,*) "restart: finished streaming of DATA to restart processes, took "//TRIM(real2string(elapsedTime))//"s"
            END IF
        ELSE
            IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
            IF(my_process_is_restart()) THEN
              !dedicated proc mode: restart processes
              totalBytesWritten = p_reduce(me%bytesWritten, p_sum_op(), 0, p_comm_work)
            ELSE
              !joint proc mode: all processes
              totalBytesWritten = p_reduce(me%bytesWritten, p_sum_op(), 0, p_comm_work_restart)
            END IF
            IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
            !the previous CALL has synchronized us, so now it's safe
            !to stop the timer
            elapsedTime = p_mpi_wtime() - me%startTime
            IF(my_process_is_restart_master()) THEN
                gibibytesWritten = REAL(totalBytesWritten, dp)/1024.0_dp/1024.0_dp/1024.0_dp
                WRITE(0,*) "restart: finished writing, "//TRIM(real2string(gibibytesWritten))//"GiB of data in "//&
                           &TRIM(real2string(elapsedTime))//"s ("//TRIM(real2string(gibibytesWritten/elapsedTime))//"GiB/s)"
            END IF
        END IF
#endif
    END SUBROUTINE multifileRestartDescriptor_writeRestartInternal

    !This dumps the restart attributes AND the namelist archive into
    !the restartFile.mfr/attributes.nc file.  Should ONLY be called on
    !a single process.
    SUBROUTINE multifileRestartDescriptor_writeAttributeFile(me, filename, restartAttributes, namelists)
        CLASS(t_MultifileRestartDescriptor), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: filename
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
        TYPE(t_NamelistArchive), INTENT(INOUT) :: namelists

        INTEGER :: file, grid, zaxis, variable, vlist
        CHARACTER(:), ALLOCATABLE :: effectiveFilename
        CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartDescriptor_writeAttributeFile"

        IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)

        effectiveFilename = multifileAttributesPath(filename)

        !open the file as a CDI stream AND create a vlist to carry our attributes
        file = streamOpenWrite(effectiveFilename, FILETYPE_NC4)
        CALL checkCdiId(file, "error opening file '"//effectiveFilename//"' for writing")

        vlist = vlistCreate()
        CALL checkCdiId(vlist, "error creating CDI vlist")

        !XXX: Unfortunately, CDI does NOT WRITE attributes to the file
        !     unless we at least define a record, for which we need at
        !     least a variable, for which we need at least a grid AND
        !     a zaxis. So we need to create some dummy objects here.
        !     The alternative would have been to bypass the CDI layer
        !     AND USE NetCDF directly, however that would have
        !     resulted IN duplicating roughly 150 lines of code IN
        !     mo_restart_attributes AND mo_restart_namelists because
        !     we still need the CDI-based attribute writing for the
        !     legacy restart modules.  So, WHILE it feels really wrong
        !     to USE CDI AND THEN work around the problems that
        !     creates, I guess it's still the sensible thing to DO.
        grid = gridCreate(GRID_GENERIC, 1)
        CALL checkCdiId(grid, "error creating CDI grid")

        zaxis = zaxisCreate(ZAXIS_GENERIC, 1)
        CALL checkCdiId(zaxis, "error creating CDI zaxis")

        variable = vlistDefVar(vlist, grid, zaxis, TSTEP_CONSTANT)
        CALL checkCdiId(variable, "error creating CDI variable")

        !write the restart attributes AND the namelist archive to the vlist attributes
        CALL namelists%writeToCdiVlist(vlist)
        CALL restartAttributes%writeToCdiVlist(vlist)
        CALL streamDefVlist(file, vlist)
        !this IS the point where we actually convince CDI to output
        !the attributes
        CALL streamDefRecord(file, variable, 0)

        !cleanup
        CALL streamClose(file)
        CALL vlistDestroy(vlist)
        CALL gridDestroy(grid)
        CALL zaxisDestroy(zaxis)

        IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)

    CONTAINS

        !check whether the given ID IS valid, calling finish with both
        !the given error message AND the error message retrieved from
        !CDI IN CASE of error
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

    !This dumps the contents of a patchData entry into the
    !restartFile.mfr/patch<JG>_metadata file.
    SUBROUTINE multifileRestartDescriptor_writeMetadataFile(me, baseFilename, jg)
        CLASS(t_MultifileRestartDescriptor), INTENT(INOUT) :: me
        CHARACTER(*), INTENT(IN) :: baseFilename
        INTEGER, VALUE :: jg

        CHARACTER(:), ALLOCATABLE :: effectiveFilename
        TYPE(t_SerializedData) :: metadataPack
        CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartDescriptor_writeMetadataFile"

        IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)

        effectiveFilename = multifileMetadataPath(baseFilename, jg)
        CALL metadataPack%construct()
        CALL me%patchData(jg)%description%packer(kPackOp, metadataPack)
        IF(0 /= putFile(effectiveFilename, metadataPack%getData(), INT(o'640', C_INT))) THEN
            CALL finish(routine, "error while writing metadata file at '"//effectiveFilename//"'")
        END IF
        CALL metadataPack%destruct()

        IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
    END SUBROUTINE multifileRestartDescriptor_writeMetadataFile

    !Collective CALL to create a batch of
    !restartFile.mfr/patch<JG>_<N>.nc files.  In the joint procs CASE,
    !this also writes the payload DATA; with dedicated writer
    !processes, this ONLY collects the DATA into buffers on the writer
    !PEs.
    FUNCTION multifileRestartDescriptor_openPayloadFile(me, filename, jg, this_datetime) RESULT(cdiIds)
        CLASS(t_MultifileRestartDescriptor), INTENT(INOUT) :: me
        CHARACTER(*), INTENT(IN) :: filename
        INTEGER, VALUE :: jg
        TYPE(datetime), POINTER, INTENT(IN) :: this_datetime

        TYPE(t_MultifilePatchData), POINTER :: patchData
        TYPE(t_CdiIds) :: cdiIds
        CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartDescriptor_openPayloadFile"

        patchData => toMultifilePatchData(me%patchData(jg))
        IF(.NOT.ASSOCIATED(patchData%varData)) RETURN !no variables to WRITE for this patch -> no restart file


        !the writer procs open a CDI stream AND WRITE the global index
        !arrays to it
        IF(my_process_is_restart_writer()) THEN
            cdiIds = patchData%openPayloadFile(filename, this_datetime, me%bytesWritten)
        END IF

        IF(my_process_is_restart_writer() .AND. .NOT. isDedicatedProcMode()) THEN
            CALL patchData%collectData(cdiIds%file, me%bytesWritten)
        ELSE
            CALL patchData%collectData()
        END IF
    END FUNCTION multifileRestartDescriptor_openPayloadFile

    ! Communication free CALL to finish writing a batch of
    ! restartFile.mfr/patch<JG>_<N>.nc files.  If the process IS a
    ! dedicated restart writer, this actually writes the payload DATA
    ! from the IN-memory buffers to the file.  On the work processes,
    ! this IS basically a NOOP, so that they may resume their work
    ! immediately WHILE the dedicated restart writers perform the
    ! actual writing.
    SUBROUTINE multifileRestartDescriptor_closePayloadFile(me, jg, cdiIds)
        CLASS(t_MultifileRestartDescriptor), INTENT(INOUT) :: me
        INTEGER, VALUE :: jg
        TYPE(t_CdiIds) :: cdiIds

        TYPE(t_MultifilePatchData), POINTER :: patchData

        patchData => toMultifilePatchData(me%patchData(jg))

        IF(my_process_is_restart_writer()) THEN
            IF(isDedicatedProcMode()) CALL patchData%writeFromBuffer(cdiIds%file, me%bytesWritten)
            CALL cdiIds%closeAndDestroyIds()
        END IF
    END SUBROUTINE multifileRestartDescriptor_closePayloadFile

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

        IF(.NOT.my_process_is_restart()) CALL finish(routine, "assertion failed: routine entered by wrong process")

        CALL restartDescriptor%construct('')    !the model TYPE will be supplied by the work processes

        DO
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
                    RETURN

                CASE DEFAULT
                    CALL finish(routine, "assertion failed: illegal restart operation received from work processes (the received &
                                         &value is "//TRIM(int2string(restartOp))//")")

            END SELECT
        END DO
    END SUBROUTINE multifileRestart_mainLoop

END MODULE mo_multifile_restart
