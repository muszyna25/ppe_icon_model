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
  USE mo_exception,          ONLY: get_filename_noext, finish
  USE mo_fortran_tools,      ONLY: assign_if_present_allocatable
  USE mo_impl_constants,     ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
  USE mo_io_config,          ONLY: restartWritingParameters, kMultifileRestartModule
  USE mo_kind,               ONLY: i8
  USE mo_packed_message,     ONLY: t_PackedMessage, kPackOp
  USE mo_run_config,         ONLY: restart_filename
  USE mo_std_c_lib,          ONLY: strerror
  USE mo_util_file,          ONLY: createSymlink
  USE mo_util_string,        ONLY: int2string, associate_keyword, with_keywords, t_keyword_list
  USE mtime,                 ONLY: datetime, newDatetime, deallocateDatetime
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_read_netcdf_distributed, ONLY: nf
#ifndef NOMPI
  USE mo_mpi, ONLY: p_pe_work, my_process_is_restart
  USE mpi, ONLY: MPI_PROC_NULL, MPI_ROOT
#endif

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'netcdf.inc'

  PUBLIC :: t_restart_args, t_rfids
  PUBLIC :: getRestartFilename, restartSymlinkName, create_restart_file_link
  PUBLIC :: restartBcastRoot
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

  TYPE t_rfids
    LOGICAL :: isinit = .FALSE.
    INTEGER :: ncid, nlayids, ftid, gids(3)
    INTEGER, ALLOCATABLE :: layids(:), nlays(:)
  CONTAINS
    PROCEDURE :: init => rfids_init
    PROCEDURE :: def_ncdfvar => rfids_def_ncdfvar
  END TYPE t_rfids

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_util"

CONTAINS

  ! Broadcast root for intercommunicator broadcasts from compute PEs to restart
  ! PEs using p_comm_work_2_restart.
  INTEGER FUNCTION restartBcastRoot() RESULT(resultVar)
    resultVar = 0
#ifndef NOMPI
    IF(.NOT.my_process_is_restart()) THEN
      ! Special root setting for intercommunicators:
      ! The PE really sending must use MPI_ROOT, the others MPI_PROC_NULL.
      IF(p_pe_work == 0) THEN
        resultVar = MPI_ROOT
      ELSE
        resultVar = MPI_PROC_NULL
      END IF
    END IF
#endif
  END FUNCTION restartBcastRoot

  SUBROUTINE getRestartFilename(baseName, jg, restartArgs, resultVar)
    CHARACTER(LEN = *), INTENT(IN) :: baseName
    INTEGER, INTENT(in) :: jg
    TYPE(t_restart_args), INTENT(IN) :: restartArgs
    CHARACTER(LEN = :), ALLOCATABLE, INTENT(INOUT) :: resultVar
    CHARACTER(LEN=32) :: datetimeString
    INTEGER :: restartModule
    TYPE(t_keyword_list), POINTER :: keywords
    TYPE(datetime), POINTER :: dt

    dt => restartArgs%restart_datetime
    WRITE (datetimeString,'(i4.4,2(i2.2),a,3(i2.2),a)')    &
       & dt%date%year, dt%date%month, dt%date%day , 'T', &
       & dt%time%hour, dt%time%minute, dt%time%second, 'Z'
    NULLIFY(keywords)
    ! build the keyword list
    CALL associate_keyword("<gridfile>", TRIM(get_filename_noext(baseName)), keywords)
    CALL associate_keyword("<idom>", TRIM(int2string(jg, "(i2.2)")), keywords)
    CALL associate_keyword("<rsttime>", TRIM(datetimeString), keywords)
    CALL associate_keyword("<mtype>", TRIM(restartArgs%modelType), keywords)
    CALL restartWritingParameters(opt_restartModule = restartModule)
    IF(restartModule == kMultifileRestartModule) THEN
      CALL associate_keyword("<extension>", "mfr", keywords)
    ELSE
      CALL associate_keyword("<extension>", "nc", keywords)
    END IF
    ! replace keywords in file name
    resultVar = TRIM(with_keywords(keywords, TRIM(restart_filename)))
  END SUBROUTINE getRestartFilename

  SUBROUTINE restartSymlinkName(modelType, jg, resultVar, opt_ndom)
    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: resultVar
    CHARACTER(*), INTENT(IN) :: modelType
    INTEGER, VALUE :: jg
    INTEGER, INTENT(IN) :: opt_ndom
    CHARACTER(2) :: sg

    IF (opt_ndom == 1) jg = 1
    WRITE (sg, '(i2.2)') jg
    resultVar = 'restart_'//modelType//"_DOM"//sg//'.nc'
  END SUBROUTINE restartSymlinkName

  SUBROUTINE create_restart_file_link(filename, modelType, jg, opt_ndom)
    CHARACTER(LEN = *), INTENT(IN) :: filename, modelType
    INTEGER, INTENT(in) :: jg
    INTEGER, INTENT(IN) :: opt_ndom
    INTEGER :: ierr
    CHARACTER(:), ALLOCATABLE :: linkname
    CHARACTER(*), PARAMETER :: routine = modname//':create_restart_file_link'

    CALL restartSymlinkName(modelType, jg, linkname, opt_ndom)
    ierr = createSymlink(filename, linkname)
    IF(ierr /= SUCCESS) CALL finish(routine, "error creating symlink at '"//linkname//"': "//strerror(ierr))
  END SUBROUTINE create_restart_file_link

  SUBROUTINE restartArgs_construct(me, this_datetime, jstep, modelType, opt_output_jfile)
    CLASS(t_restart_args), INTENT(INOUT) :: me
    TYPE(datetime), INTENT(IN) :: this_datetime
    INTEGER, INTENT(in) :: jstep
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
    INTEGER, INTENT(in) :: operation
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

  SUBROUTINE rfids_init(rfids, ncid_in, nelem)
    CLASS(t_rfids), INTENT(OUT) :: rfids
    INTEGER, INTENT(IN) :: ncid_in, nelem(3)

    rfids%ncid = ncid_in
    IF (rfids%isinit) DEALLOCATE(rfids%layids, rfids%nlays)
    rfids%nlayids = 0
    ALLOCATE(rfids%layids(32), rfids%nlays(32))
    CALL nf(nf_def_dim(rfids%ncid, "fake_time", NF_UNLIMITED, rfids%ftid))
    CALL nf(nf_def_dim(rfids%ncid, "cells", nelem(1), rfids%gids(1)))
    CALL nf(nf_def_dim(rfids%ncid, "verts", nelem(2), rfids%gids(2)))
    CALL nf(nf_def_dim(rfids%ncid, "edges", nelem(3), rfids%gids(3)))
    rfids%isinit = .TRUE.
  END SUBROUTINE rfids_init

  SUBROUTINE rfids_def_ncdfvar(rfids, ci, gid)
    CLASS(t_rfids), INTENT(INOUT) :: rfids
    INTEGER, INTENT(IN) :: gid
    TYPE(t_var_metadata), POINTER, INTENT(INOUT) :: ci
    INTEGER :: dtid, ivid, iivid, dids(3), nd
    INTEGER, ALLOCATABLE :: tmpa(:), tmpb(:)
    CHARACTER(LEN=10) :: layname

    IF (.NOT.rfids%isinit) CALL finish(modname//":rfids_def_ncdfvar", "not initialized")
    nd = ci%ndims
    dids(1) = rfids%gids(gid)
    SELECT CASE(ci%data_type)
    CASE(REAL_T)
      dtid = NF_DOUBLE
    CASE(SINGLE_T)
      dtid = NF_REAL
    CASE(INT_T)
      dtid = NF_INT
    END SELECT
    SELECT CASE(nd)
    CASE(2)
      dids(2) = rfids%ftid
    CASE(3)
      ivid = 0
      DO iivid = 1, rfids%nlayids
        IF (rfids%nlays(iivid) .EQ. ci%used_dimensions(2)) THEN
          ivid = iivid
          EXIT
        END IF
      END DO
      IF (ivid .EQ. 0) THEN
        IF (rfids%nlayids .EQ. SIZE(rfids%nlays)) THEN
          CALL MOVE_ALLOC(rfids%nlays, tmpa)
          CALL MOVE_ALLOC(rfids%layids, tmpb)
          ALLOCATE(rfids%nlays(rfids%nlayids+32), rfids%layids(rfids%nlayids+32))
          rfids%nlays(:rfids%nlayids) = tmpa
          rfids%layids(:rfids%nlayids) = tmpb
        END IF
        rfids%nlayids = rfids%nlayids + 1
        rfids%nlays(rfids%nlayids) = ci%used_dimensions(2)
        WRITE(layname, "(a,i0)") "layers_", rfids%nlays(rfids%nlayids)
        CALL nf(nf_def_dim(rfids%ncid, TRIM(layname), &
          &                rfids%nlays(rfids%nlayids), rfids%layids(rfids%nlayids)))
        ivid = rfids%nlayids
      END IF
      dids(2) = rfids%layids(ivid)
      dids(3) = rfids%ftid
    END SELECT
    CALL nf(nf_def_var(rfids%ncid, TRIM(ci%name), dtid, nd, dids(:nd), ci%cdiVarId))
  END SUBROUTINE rfids_def_ncdfvar

END MODULE mo_restart_util
