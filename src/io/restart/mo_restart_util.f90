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
    USE mo_cdi, ONLY: CDI_UNDEFID
    USE mo_cf_convention, ONLY: cf_global_info
    USE mo_datetime, ONLY: t_datetime, iso8601, iso8601extended
    USE mo_exception, ONLY: get_filename_noext, finish, message
    USE mo_fortran_tools, ONLY: assign_if_present, assign_if_present_allocatable
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_restart_attributes, ONLY: t_RestartAttributeList
    USE mo_kind, ONLY: wp, i8
    USE mo_packed_message, ONLY: t_PackedMessage, kPackOp, kUnpackOp
    USE mo_run_config, ONLY: restart_filename
    USE mo_util_file, ONLY: util_symlink, util_islink, util_unlink
    USE mo_util_string, ONLY: int2string, real2string, associate_keyword, with_keywords, t_keyword_list

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: t_restart_args

    PUBLIC :: getRestartFilename
    PUBLIC :: setGeneralRestartAttributes
    PUBLIC :: setDynamicPatchRestartAttributes
    PUBLIC :: setPhysicsRestartAttributes
    PUBLIC :: create_restart_file_link

    ! patch independent restart arguments
    TYPE t_restart_args
        TYPE(t_datetime) :: datetime
        INTEGER :: jstep
        CHARACTER(LEN = 32) :: modelType
        INTEGER, ALLOCATABLE :: output_jfile(:)
    CONTAINS
        PROCEDURE :: construct => restartArgs_construct
        PROCEDURE :: packer => restartArgs_packer   ! unpacking IS considered construction
        PROCEDURE :: print => restartArgs_print
        PROCEDURE :: destruct => restartArgs_destruct
    END TYPE t_restart_args

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_util"

CONTAINS

    FUNCTION getRestartFilename(baseName, domain, datetime, modelTypeName) RESULT(resultVar)
        CHARACTER(LEN = *), INTENT(IN) :: baseName, modelTypeName
        INTEGER, VALUE :: domain
        TYPE(t_datetime), INTENT(IN) :: datetime
        CHARACTER(LEN = :), ALLOCATABLE :: resultVar

        CHARACTER(LEN=32) :: datetimeString
        TYPE(t_keyword_list), POINTER :: keywords => NULL()

        datetimeString = iso8601(datetime)

        ! build the keyword list
        CALL associate_keyword("<gridfile>", TRIM(get_filename_noext(baseName)), keywords)
        CALL associate_keyword("<idom>", TRIM(int2string(domain, "(i2.2)")), keywords)
        CALL associate_keyword("<rsttime>", TRIM(datetimeString), keywords)
        CALL associate_keyword("<mtype>", TRIM(modelTypeName), keywords)

        ! replace keywords in file name
        resultVar = TRIM(with_keywords(keywords, TRIM(restart_filename)))
    END FUNCTION getRestartFilename

    SUBROUTINE setGeneralRestartAttributes(restartAttributes, datetime, n_dom, jstep, opt_output_jfile)
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
        TYPE(t_datetime), INTENT(IN) :: datetime
        INTEGER, VALUE :: n_dom, jstep
        INTEGER, OPTIONAL, INTENT(IN) :: opt_output_jfile(:)

        INTEGER :: i

        ! set CF-Convention required restart attributes
        CALL restartAttributes%setText('title',       TRIM(cf_global_info%title))
        CALL restartAttributes%setText('institution', TRIM(cf_global_info%institution))
        CALL restartAttributes%setText('source',      TRIM(cf_global_info%source))
        CALL restartAttributes%setText('history',     TRIM(cf_global_info%history))
        CALL restartAttributes%setText('references',  TRIM(cf_global_info%references))
        CALL restartAttributes%setText('comment',     TRIM(cf_global_info%comment))

        CALL restartAttributes%setReal( 'current_caltime', datetime%caltime )
        CALL restartAttributes%setInteger( 'current_calday' , INT(datetime%calday) )   !FIXME: Either it IS a bug that calday IS a 64bit INTEGER, OR it IS a bug that ONLY 32 bit of it are stored IN the restart file. Either way this needs to be fixed.
        CALL restartAttributes%setReal( 'current_daysec' , datetime%daysec )
        CALL restartAttributes%setText('tc_startdate', iso8601extended(datetime))   ! in preparation for move to mtime

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

    SUBROUTINE setPhysicsRestartAttributes(restartAttributes, jg, t_elapsed_phy, lcall_phy)
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
        INTEGER, VALUE :: jg
        REAL(wp), INTENT(IN) :: t_elapsed_phy(:)
        LOGICAL, INTENT(IN) :: lcall_phy(:)

        INTEGER :: i
        CHARACTER(LEN = :), ALLOCATABLE :: prefix

        prefix = 't_elapsed_phy_DOM'//TRIM(int2string(jg, "(i2.2)"))//'_PHY'
        DO i = 1, SIZE(t_elapsed_phy)
            CALL restartAttributes%setReal(prefix//TRIM(int2string(i, '(i2.2)')), t_elapsed_phy(i) )
        END DO

        prefix = 'lcall_phy_DOM'//TRIM(int2string(jg, "(i2.2)"))//'_PHY'
        DO i = 1, SIZE(lcall_phy)
            CALL restartAttributes%setLogical(prefix//TRIM(int2string(i, '(i2.2)')), lcall_phy(i) )
        END DO
    END SUBROUTINE setPhysicsRestartAttributes

    SUBROUTINE create_restart_file_link(filename, modelType, proc_id, jg, opt_ndom)
        CHARACTER(LEN = *), INTENT(IN) :: filename, modelType
        INTEGER, VALUE :: proc_id, jg
        INTEGER, INTENT(IN), OPTIONAL :: opt_ndom

        INTEGER :: iret, ndom
        CHARACTER(LEN = 12) :: procIdString
        CHARACTER(LEN = 64) :: linkname
        CHARACTER(LEN=*), PARAMETER :: routine = modname//':create_restart_file_link'

        ! we need to add a process dependent part to the link NAME IF there are several restart processes
        procIdString = ''
        IF(proc_id /= 0) procIdString = TRIM(int2string(proc_id))

        ! IN CASE we have ONLY a single domain / no domain information, USE "_DOM01" IN the link NAME
        ndom = 1
        CALL assign_if_present(ndom, opt_ndom)
        IF(ndom == 1) jg = 1

        ! build link name
        linkname = 'restart'//TRIM(procIdString)//'_'//modelType//"_DOM"//TRIM(int2string(jg, "(i2.2)"))//'.nc'

        ! delete old symbolic link, if exists
        ! FIXME[NH]: handle the CASE that we have a file at that location which IS NOT a symlink
        IF(util_islink(TRIM(linkname))) THEN
            iret = util_unlink(TRIM(linkname))
            IF(iret /= SUCCESS) WRITE(0, *) routine//': cannot unlink "'//TRIM(linkname)//'"'
        ENDIF

        ! create a new symbolic link
        iret = util_symlink(filename,TRIM(linkname))
        IF(iret /= SUCCESS) WRITE(0, *) routine//': cannot create symbolic link "'//TRIM(linkname)//'" for "'//filename//'"'
    END SUBROUTINE create_restart_file_link

    SUBROUTINE restartArgs_construct(me, datetime, jstep, modelType, opt_output_jfile)
        CLASS(t_restart_args), INTENT(INOUT) :: me
        TYPE(t_datetime), INTENT(IN) :: datetime
        INTEGER, VALUE :: jstep
        CHARACTER(LEN = *), INTENT(IN) :: modelType
        INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

        me%datetime = datetime
        me%jstep = jstep
        me%modelType = modelType
        CALL assign_if_present_allocatable(me%output_jfile, opt_output_jfile)
    END SUBROUTINE restartArgs_construct

    SUBROUTINE restartArgs_packer(me, operation, message)
        CLASS(t_restart_args), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        TYPE(t_PackedMessage), INTENT(INOUT) :: message

        INTEGER :: calday

        CALL message%packer(operation, me%datetime%year)
        CALL message%packer(operation, me%datetime%month)
        CALL message%packer(operation, me%datetime%day)
        CALL message%packer(operation, me%datetime%hour)
        CALL message%packer(operation, me%datetime%minute)
        CALL message%packer(operation, me%datetime%second)
        CALL message%packer(operation, me%datetime%caltime)
        IF(operation == kPackOp) calday = INT(me%datetime%calday)   !The IF IS needed to avoid overflow when me%datetime%calday IS uninitialized.
        CALL message%packer(operation, calday)
        IF(operation == kUnpackOp) me%datetime%calday = INT(calday,i8)
        CALL message%packer(operation, me%datetime%daysec)
        CALL message%packer(operation, me%jstep)
        CALL message%packer(operation, me%modelType)
        CALL message%packer(operation, me%output_jfile)
    END SUBROUTINE restartArgs_packer

    SUBROUTINE restartArgs_print(me, prefix)
        CLASS(t_restart_args), INTENT(IN) :: me
        CHARACTER(LEN = *), INTENT(IN) :: prefix

        PRINT*, prefix//'current_calday='//TRIM(int2string(INT(me%datetime%calday)))
        PRINT*, prefix//'current_caltime='//TRIM(real2string(me%datetime%caltime))
        PRINT*, prefix//'current_daysec='//TRIM(real2string(me%datetime%daysec))
    END SUBROUTINE restartArgs_print

    SUBROUTINE restartArgs_destruct(me)
        CLASS(t_restart_args), INTENT(INOUT) :: me

        IF(ALLOCATED(me%output_jfile)) DEALLOCATE(me%output_jfile)
    END SUBROUTINE restartArgs_destruct

END MODULE mo_restart_util
