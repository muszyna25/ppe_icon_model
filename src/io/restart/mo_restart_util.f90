!>
!! Contains common helper routines for(a)synchronous restart
!! ----------------------------------------------------------
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_restart_util
    USE mo_cf_convention,      ONLY: cf_global_info
    USE mo_exception,          ONLY: get_filename_noext, finish
    USE mo_fortran_tools,      ONLY: assign_if_present, assign_if_present_allocatable
    USE mo_impl_constants,     ONLY: SUCCESS
    USE mo_io_config,          ONLY: restartWritingParameters, kMultifileRestartModule, ALL_WORKERS_INVOLVED
    USE mo_kind,               ONLY: wp, i8
    USE mo_mpi,                ONLY: num_work_procs, my_process_is_restart, my_process_is_work, p_pe, &
      &                              p_restart_pe0, p_work_pe0, stop_mpi, p_n_work
    USE mo_packed_message,     ONLY: t_PackedMessage, kPackOp
    USE mo_restart_attributes, ONLY: t_RestartAttributeList
    USE mo_run_config,         ONLY: restart_filename
    USE mo_std_c_lib,          ONLY: strerror
    USE mo_timer,              ONLY: ltimer,timer_stop, timer_model_init, print_timer
    USE mo_util_file,          ONLY: createSymlink
    USE mo_util_string,        ONLY: int2string, associate_keyword, with_keywords, t_keyword_list
    USE mtime,                 ONLY: datetime, newDatetime, deallocateDatetime, datetimeToString, &
      &                              MAX_DATETIME_STR_LEN
    
    IMPLICIT NONE

    PRIVATE

    PUBLIC :: t_restart_args
    PUBLIC :: alloc_string

    PUBLIC :: workProcCount
    PUBLIC :: dedicatedRestartProcCount
    PUBLIC :: restartProcCount
    PUBLIC :: isDedicatedProcMode
    PUBLIC :: restartWriterId
    PUBLIC :: my_process_is_restart_master
    PUBLIC :: my_process_is_restart_writer
    PUBLIC :: restartWorkProcId
    PUBLIC :: restartWorkProcId2Rank

    PUBLIC :: becomeDedicatedRestartProc
    PUBLIC :: shutdownRestartProc

    PUBLIC :: getRestartFilename
    PUBLIC :: setGeneralRestartAttributes
    PUBLIC :: setDynamicPatchRestartAttributes
    PUBLIC :: setPhysicsRestartAttributes
    PUBLIC :: restartSymlinkName
    PUBLIC :: create_restart_file_link

    ! patch independent restart arguments
    TYPE t_restart_args
        TYPE(datetime), POINTER :: restart_datetime => NULL()
        INTEGER :: jstep
        CHARACTER(LEN = 32) :: modelType
        INTEGER, ALLOCATABLE :: output_jfile(:)
    CONTAINS
        PROCEDURE :: construct => restartArgs_construct
        PROCEDURE :: packer => restartArgs_packer   ! unpacking IS considered construction
        PROCEDURE :: destruct => restartArgs_destruct
    END TYPE t_restart_args

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_util"

CONTAINS

    SUBROUTINE alloc_string(str_len, str)
        INTEGER, INTENT(IN) :: str_len
        CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: str

        IF(ALLOCATED(str)) THEN
          IF(LEN(str) .NE. str_len) THEN
            DEALLOCATE(str)
            ALLOCATE(CHARACTER(LEN=str_len) :: str)
          ENDIF
        ELSE
          ALLOCATE(CHARACTER(LEN=str_len) :: str)
        ENDIF
    END SUBROUTINE alloc_string

    INTEGER FUNCTION workProcCount() RESULT(resultVar)
        resultVar = num_work_procs
    END FUNCTION workProcCount

    INTEGER FUNCTION dedicatedRestartProcCount() RESULT(resultVar)
        CALL restartWritingParameters(opt_dedicatedProcCount = resultVar)
    END FUNCTION dedicatedRestartProcCount

    INTEGER FUNCTION restartProcCount() RESULT(resultVar)
        CALL restartWritingParameters(opt_restartProcCount = resultVar)
        IF (resultVar == ALL_WORKERS_INVOLVED)  resultVar = p_n_work
    END FUNCTION restartProcCount

    LOGICAL FUNCTION isDedicatedProcMode() RESULT(resultVar)
        CALL restartWritingParameters(opt_lDedicatedProcMode = resultVar)
    END FUNCTION isDedicatedProcMode

    !Returns a VALUE IN the range [0, restartProcCount - 1] IF the
    !process IS a restart writer (either a work OR a dedicated restart
    !process according to lDedicatedProcMode), OR -1 on all other
    !processes.
    INTEGER FUNCTION restartWriterId() RESULT(resultVar)
        LOGICAL :: lDedicatedProcMode

        CALL restartWritingParameters(opt_lDedicatedProcMode = lDedicatedProcMode)
        resultVar = -1
        IF(lDedicatedProcMode) THEN
            IF(my_process_is_restart()) resultVar = p_pe - p_restart_pe0
        ELSE
            IF(my_process_is_work()) resultVar = p_pe - p_work_pe0
            IF(resultVar >= restartProcCount()) resultVar = -1
        END IF
    END FUNCTION restartWriterId

    !In dedicated procs mode, this returns TRUE on the first restart PE,
    !iN joint procs mode, this returns TRUE on the first work PE.
    LOGICAL FUNCTION my_process_is_restart_master() RESULT(resultVar)
        resultVar = restartWriterId() == 0
    END FUNCTION my_process_is_restart_master

    !In dedicated procs mode, this returns TRUE on the restart PEs,
    !IN joint procs mode, this returns TRUE on the first few work PEs (according to restartProcCount).
    LOGICAL FUNCTION my_process_is_restart_writer() RESULT(resultVar)
        resultVar = restartWriterId() >= 0
    END FUNCTION my_process_is_restart_writer

    !Returns a unique process Id for all restart AND work
    !processes. The processes are numbered as follows:
    !0 ... restartProcCount-1: restart processes
    !restartProcCount ... restartWorkProcCount-1: work processes that are NOT also restart processes
    !
    !TODO: Implement using another communicator AND
    !      MPI_Group_translate_ranks().  It would be good to DO this
    !      together with choosing a better restart proc placement, see
    !      comment IN
    !      mo_multifile_restart:multifileRestartDescriptor_construct().
    INTEGER FUNCTION restartWorkProcId() RESULT(resultVar)
        LOGICAL :: lDedicatedProcMode
        INTEGER :: restartProcCount

        CALL restartWritingParameters(opt_lDedicatedProcMode = lDedicatedProcMode, &
          &                           opt_restartProcCount = restartProcCount)
        IF (restartProcCount == ALL_WORKERS_INVOLVED)  restartProcCount = p_n_work

        IF(lDedicatedProcMode .AND. my_process_is_restart()) THEN
            resultVar = p_pe - p_restart_pe0
        ELSE IF(lDedicatedProcMode) THEN
            resultVar = restartProcCount + p_pe - p_work_pe0
        ELSE
            resultVar = p_pe - p_work_pe0
        END IF
    END FUNCTION restartWorkProcId

    !Unfortunately, the restart PEs have ranks higher than the work
    !PEs, which clashes with the way we defined the
    !restartWorkProcId() above.
    !Nevetheless, it IS good to have the restartWorkProcId() as
    !defined above, because there we know that we always have the
    !restart writers IN the same position.  So, this FUNCTION maps our
    !restartWorkProcId()s to ranks IN p_comm_work_restart.
    INTEGER FUNCTION restartWorkProcId2Rank(procId) RESULT(rank)
        INTEGER, VALUE :: procId

        LOGICAL :: lDedicatedProcMode
        INTEGER :: restartProcCount

        CALL restartWritingParameters(opt_lDedicatedProcMode = lDedicatedProcMode, &
          &                           opt_restartProcCount = restartProcCount)
        IF (restartProcCount == ALL_WORKERS_INVOLVED)  restartProcCount = p_n_work

        IF(lDedicatedProcMode .AND. procId < restartProcCount) THEN
            rank = num_work_procs + procId
        ELSE IF(lDedicatedProcMode) THEN
            rank = procId - restartProcCount
        ELSE
            rank = procId
        END IF
    END FUNCTION restartWorkProcId2Rank

    ! Performs all actions needed to decouple the calling process from the rest of the program.
    SUBROUTINE becomeDedicatedRestartProc()
        ! Dedicated restart processes are detached during the
        ! model_init phase, so the corresponding timer IS still
        ! running.  It must be stopped before the timer table IS
        ! printed, preferably before we start the restart timers, so
        ! better DO it right away.
        IF(ltimer) CALL timer_stop(timer_model_init)
    END SUBROUTINE becomeDedicatedRestartProc

    ! Does all the tidying up that's necessary before executing a STOP.
    ! Does NOT RETURN.
    SUBROUTINE shutdownRestartProc()
        IF(ltimer) CALL print_timer
        CALL stop_mpi

        STOP
    END SUBROUTINE shutdownRestartProc

    SUBROUTINE getRestartFilename(baseName, domain, restartArgs, resultVar)
        USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH
        CHARACTER(LEN = *), INTENT(IN) :: baseName
        INTEGER, VALUE :: domain
        TYPE(t_restart_args), INTENT(IN) :: restartArgs
        CHARACTER(LEN = :), ALLOCATABLE, INTENT(INOUT) :: resultVar

        CHARACTER(LEN=32) :: datetimeString
        INTEGER :: restartModule, fn_len
        TYPE(t_keyword_list), POINTER :: keywords
        TYPE(datetime), POINTER :: dt
        CHARACTER(len=MAX_CHAR_LENGTH) :: tempname

        dt => restartArgs%restart_datetime
        WRITE (datetimeString,'(i4.4,2(i2.2),a,3(i2.2),a)')    &
             & dt%date%year, dt%date%month, dt%date%day , 'T', &
             & dt%time%hour, dt%time%minute, dt%time%second, 'Z'
        NULLIFY(keywords)
        ! build the keyword list
        CALL associate_keyword("<gridfile>", TRIM(get_filename_noext(baseName)), keywords)
        CALL associate_keyword("<idom>", TRIM(int2string(domain, "(i2.2)")), keywords)
        CALL associate_keyword("<rsttime>", TRIM(datetimeString), keywords)
        CALL associate_keyword("<mtype>", TRIM(restartArgs%modelType), keywords)
        CALL restartWritingParameters(opt_restartModule = restartModule)
        IF(restartModule == kMultifileRestartModule) THEN
            CALL associate_keyword("<extension>", "mfr", keywords)
        ELSE
            CALL associate_keyword("<extension>", "nc", keywords)
        END IF

        ! replace keywords in file name
        WRITE(tempname, "(a)") TRIM(with_keywords(keywords, TRIM(restart_filename)))
        fn_len = LEN_TRIM(tempname)
        CALL alloc_string(fn_len, resultVar)
        WRITE(resultVar, '(a)') TRIM(tempname)
    END SUBROUTINE getRestartFilename

    SUBROUTINE setGeneralRestartAttributes(restartAttributes, this_datetime, n_dom, jstep, opt_output_jfile)
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
        TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
        INTEGER, VALUE :: n_dom, jstep
        INTEGER, OPTIONAL, INTENT(IN) :: opt_output_jfile(:)

        CHARACTER(len=MAX_DATETIME_STR_LEN) :: dstring
        INTEGER :: i
        
        ! set CF-Convention required restart attributes
        CALL restartAttributes%setText('title',       TRIM(cf_global_info%title))
        CALL restartAttributes%setText('institution', TRIM(cf_global_info%institution))
        CALL restartAttributes%setText('source',      TRIM(cf_global_info%source))
        CALL restartAttributes%setText('history',     TRIM(cf_global_info%history))
        CALL restartAttributes%setText('references',  TRIM(cf_global_info%references))
        CALL restartAttributes%setText('comment',     TRIM(cf_global_info%comment))
        CALL datetimeToString(this_datetime, dstring)
        CALL restartAttributes%setText('tc_startdate', TRIM(dstring))   ! in preparation for move to mtime

        ! no. of domains AND simulation step
        CALL restartAttributes%setInteger( 'n_dom', n_dom)
        CALL restartAttributes%setInteger( 'jstep', jstep )

        IF(PRESENT(opt_output_jfile)) THEN
            DO i = 1, SIZE(opt_output_jfile)
                CALL restartAttributes%setInteger('output_jfile_'//TRIM(int2string(i, '(i2.2)')), opt_output_jfile(i) )
            END DO
        END IF
    END SUBROUTINE setGeneralRestartAttributes

    SUBROUTINE setDynamicPatchRestartAttributes(restartAttributes, jg, nold, nnow, nnew, nnow_rcf, nnew_rcf)
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
        INTEGER, VALUE :: jg, nold, nnow, nnew, nnow_rcf, nnew_rcf

        CHARACTER(LEN = 2) :: jgString

        jgString = int2string(jg, "(i2.2)")

        CALL restartAttributes%setInteger('nold_DOM'//jgString, nold)
        CALL restartAttributes%setInteger('nnow_DOM'//jgString, nnow)
        CALL restartAttributes%setInteger('nnew_DOM'//jgString, nnew)
        CALL restartAttributes%setInteger('nnow_rcf_DOM'//jgString, nnow_rcf)
        CALL restartAttributes%setInteger('nnew_rcf_DOM'//jgString, nnew_rcf)
    END SUBROUTINE setDynamicPatchRestartAttributes


    SUBROUTINE setPhysicsRestartAttributes(restartAttributes, jg, opt_t_elapsed_phy)
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
        INTEGER, VALUE :: jg
        REAL(wp), OPTIONAL, INTENT(IN) :: opt_t_elapsed_phy(:)

        INTEGER :: i, fn_len
        CHARACTER(LEN = :), ALLOCATABLE :: prefix

        ! F2008: A null pointer or unallocated allocatable can be used to denote an absent optional argument
        IF (PRESENT(opt_t_elapsed_phy)) THEN
          fn_len = LEN_TRIM('t_elapsed_phy_DOM'//TRIM(int2string(jg, "(i2.2)"))//'_PHY')
          ALLOCATE(CHARACTER(LEN=fn_len) :: prefix)
          prefix = 't_elapsed_phy_DOM'//TRIM(int2string(jg, "(i2.2)"))//'_PHY'
          DO i = 1, SIZE(opt_t_elapsed_phy)
              CALL restartAttributes%setReal(prefix//TRIM(int2string(i, '(i2.2)')), opt_t_elapsed_phy(i) )
          END DO
        ENDIF

    END SUBROUTINE setPhysicsRestartAttributes

    SUBROUTINE restartSymlinkName(modelType, jg, resultVar, opt_ndom)
        CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: resultVar
        CHARACTER(*), INTENT(IN) :: modelType
        INTEGER, VALUE :: jg
        INTEGER, INTENT(IN), OPTIONAL :: opt_ndom

        INTEGER :: ndom, fn_len

        ! IN CASE we have ONLY a single domain / no domain information, USE "_DOM01" IN the link NAME
        ndom = 1
        CALL assign_if_present(ndom, opt_ndom)
        IF(ndom == 1) jg = 1

        ! build link name
        fn_len = LEN_TRIM('restart_'//modelType//"_DOM"//TRIM(int2string(jg, "(i2.2)"))//'.nc')
        CALL alloc_string(fn_len, resultVar)
        resultVar = 'restart_'//modelType//"_DOM"//TRIM(int2string(jg, "(i2.2)"))//'.nc'
    END SUBROUTINE restartSymlinkName

    SUBROUTINE create_restart_file_link(filename, modelType, jg, opt_ndom)
        CHARACTER(LEN = *), INTENT(IN) :: filename, modelType
        INTEGER, VALUE :: jg
        INTEGER, INTENT(IN), OPTIONAL :: opt_ndom

        INTEGER :: error
        CHARACTER(:), ALLOCATABLE :: linkname
        CHARACTER(LEN=*), PARAMETER :: routine = modname//':create_restart_file_link'

        CALL restartSymlinkName(modelType, jg, linkname, opt_ndom)
        error = createSymlink(filename, linkname)
        IF(error /= SUCCESS) CALL finish(routine, "error creating symlink at '"//linkname//"': "//strerror(error))
    END SUBROUTINE create_restart_file_link

    SUBROUTINE restartArgs_construct(me, this_datetime, jstep, modelType, opt_output_jfile)
        CLASS(t_restart_args), INTENT(INOUT) :: me
        TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
        INTEGER, VALUE :: jstep
        CHARACTER(LEN = *), INTENT(IN) :: modelType
        INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

        integer :: ierr

        me%restart_datetime => newDatetime(this_datetime%date%year, this_datetime%date%month,    &
          &                                this_datetime%date%day, this_datetime%time%hour,      &
          &                                this_datetime%time%minute, this_datetime%time%second, &
          &                                this_datetime%time%ms, ierr)
        me%jstep = jstep
        me%modelType = modelType
        CALL assign_if_present_allocatable(me%output_jfile, opt_output_jfile)
    END SUBROUTINE restartArgs_construct

    SUBROUTINE restartArgs_packer(me, operation, packedMessage)
        CLASS(t_restart_args), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        TYPE(t_PackedMessage), INTENT(INOUT) :: packedMessage
        CHARACTER(len=*), PARAMETER ::  routine = modname//':restartArgs_packer'
        
        IF (operation == kPackOp .AND. .NOT. ASSOCIATED(me%restart_datetime)) THEN
          CALL finish(routine, 'Assertion failed: cannot pack unconstructed object.')
        ENDIF
        IF (.NOT. ASSOCIATED(me%restart_datetime)) THEN
          me%restart_datetime => newDatetime(1878_i8,1,1,0,0,0,0)
        ENDIF
        
        CALL packedMessage%packer(operation, me%restart_datetime%date%year)
        CALL packedMessage%packer(operation, me%restart_datetime%date%month)
        CALL packedMessage%packer(operation, me%restart_datetime%date%day)
        CALL packedMessage%packer(operation, me%restart_datetime%time%hour)
        CALL packedMessage%packer(operation, me%restart_datetime%time%minute)
        CALL packedMessage%packer(operation, me%restart_datetime%time%second)
        CALL packedMessage%packer(operation, me%jstep)
        CALL packedMessage%packer(operation, me%modelType)
        CALL packedMessage%packer(operation, me%output_jfile)
    END SUBROUTINE restartArgs_packer

    SUBROUTINE restartArgs_destruct(me)
        CLASS(t_restart_args), INTENT(INOUT) :: me

        IF(ASSOCIATED(me%restart_datetime)) CALL deallocateDatetime(me%restart_datetime)
        IF(ALLOCATED(me%output_jfile)) DEALLOCATE(me%output_jfile)
    END SUBROUTINE restartArgs_destruct

END MODULE mo_restart_util
