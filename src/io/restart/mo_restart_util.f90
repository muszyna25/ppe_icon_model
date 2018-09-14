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
  USE mo_io_config,          ONLY: restartWritingParameters, kMultifileRestartModule
  USE mo_kind,               ONLY: wp, i8
  USE mo_mpi,                ONLY: stop_mpi
  USE mo_packed_message,     ONLY: t_PackedMessage, kPackOp
  USE mo_restart_attributes, ONLY: t_RestartAttributeList
  USE mo_run_config,         ONLY: restart_filename
  USE mo_std_c_lib,          ONLY: strerror
  USE mo_timer,              ONLY: ltimer, print_timer, timer_stop, timer_model_init
  USE mo_util_file,          ONLY: createSymlink
  USE mo_util_string,        ONLY: int2string, associate_keyword, with_keywords, t_keyword_list
  USE mtime,                 ONLY: datetime, newDatetime, deallocateDatetime, datetimeToString, &
    &                              MAX_DATETIME_STR_LEN
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_restart_args
  PUBLIC :: alloc_string

  PUBLIC :: getRestartFilename
  PUBLIC :: setGeneralRestartAttributes
  PUBLIC :: setDynamicPatchRestartAttributes
  PUBLIC :: setPhysicsRestartAttributes
  PUBLIC :: restartSymlinkName, create_restart_file_link
  PUBLIC :: becomeDedicatedRestartProc, shutdownRestartProc

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

  SUBROUTINE becomeDedicatedRestartProc()

    IF(ltimer) CALL timer_stop(timer_model_init)
  END SUBROUTINE becomeDedicatedRestartProc

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

    jgString = TRIM(int2string(jg, "(i2.2)"))
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

    ndom = 1
    CALL assign_if_present(ndom, opt_ndom)
    IF(ndom == 1) jg = 1
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
