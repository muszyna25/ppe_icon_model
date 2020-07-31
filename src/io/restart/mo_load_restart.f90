!> Module for reading restart files
!!
!! Note: The asynchronous implementation of the restart output can be
!!       found in the module "mo_async_restart"
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_load_restart
  USE mo_cdi,                ONLY: streamOpenRead, cdiStringError, streamInqVlist, streamClose
  USE mo_exception,          ONLY: message, finish, warning
  USE mo_load_multifile_restart, ONLY: multifileReadPatch, multifileCheckRestartFiles
  USE mo_load_singlefile_restart, ONLY: singlefileReadPatch, singlefileCheckRestartFiles
  USE mo_model_domain,       ONLY: t_patch
  USE mo_mpi,                ONLY: p_comm_work, p_comm_rank, my_process_is_mpi_workroot, my_process_is_stdio
  USE mo_multifile_restart_util, ONLY: multifileRestartLinkName
  USE mo_restart_nml_and_att,ONLY: restartAttributeList_read, ocean_initFromRestart_OVERRIDE
  USE mo_restart_util,       ONLY: restartSymlinkName
  USE mo_restart_var_data,   ONLY: createRestartVarData
  USE mo_var_list_element,   ONLY: t_p_var_list_element
  USE mo_timer,              ONLY: timer_start, timer_stop, timer_load_restart, timer_load_restart_io, &
    & timer_load_restart_comm_setup, timer_load_restart_communication, &
    & timer_load_restart_get_var_id, timers_level
  USE mo_util_string,        ONLY: separator, toCharacter
  USE mo_var_list_register,  ONLY: vl_iter
  USE mo_master_control,     ONLY: get_my_process_name

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_restart_files, read_restart_header

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_load_restart"

CONTAINS

  !Checks the existence of either a single- OR multifile restart
  !symlink, throwing a warning IF both exist.  Returns the path to the
  !symlink that IS to be used, out_lIsMultifile IS TRUE IF that
  !symlink points to a multifile.
  !
  !The RESULT of the first CALL IS cached, so that this FUNCTION may
  !be called several times without negative consequences.
  SUBROUTINE findRestartFile(modelType, out_lIsMultifile, linkname)
    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: linkname
    CHARACTER(*), INTENT(IN) :: modelType
    LOGICAL, INTENT(OUT) :: out_lIsMultifile

    LOGICAL, SAVE :: cache_isValid = .FALSE., cache_isMultifile
    CHARACTER(:), ALLOCATABLE :: singlefileLinkName, multifileLinkName
    LOGICAL :: haveSinglefileLink, haveMultifileLink
    CHARACTER(*), PARAMETER :: routine = modname//":findRestartFile"

    IF(modelType == '') CALL finish(routine, "assertion failed: invalid modelType")

    IF(.NOT.cache_isValid) THEN
      !get the two possible paths
      CALL restartSymlinkName(modelType, 1, singlefileLinkName, -1)
      CALL multifileRestartLinkName(modelType, multifileLinkName)
      !check whether the respective files exist
      INQUIRE(file = singlefileLinkName, exist = haveSinglefileLink)
#if defined (__INTEL_COMPILER)
      INQUIRE(directory = multifileLinkName, exist = haveMultifileLink)
#else
      INQUIRE(file = multifileLinkName, exist = haveMultifileLink)
#endif
      !determine which path to USE
      IF(haveMultifileLink) THEN
        cache_isMultifile = .TRUE.
        IF (haveSinglefileLink) CALL warning(routine, &
          & "both kinds, single- and multi-file restart-links, present, choosing multi-file")
      ELSE
        IF(haveSinglefileLink) THEN
          cache_isMultifile = .FALSE.
        ELSE
          WRITE (0,*) "singlefileLinkName = ", singlefileLinkName
          WRITE (0,*) "multifileLinkName = ", multifileLinkName
          CALL finish(routine, "no functional symlink found")
        END IF
      END IF
      cache_isValid = .TRUE.
    END IF
    out_lIsMultifile = cache_isMultifile
    IF(cache_isMultifile) THEN
      CALL multifileRestartLinkName(modelType, linkname)
    ELSE
      CALL restartSymlinkName(modelType, 1, linkname, -1)
    END IF
  END SUBROUTINE findRestartFile

  ! Read all namelists used in the previous run
  ! and store them in a buffer. These values will overwrite the
  ! model default, and will later be overwritten if the user has
  ! specified something different for this integration.
  !
  ! Collective across p_comm_work.
  !
  ! This IS used for both single- AND multifile restart files.
  SUBROUTINE readRestartAttributeFile(attributeFile)
    CHARACTER(*), INTENT(IN) :: attributeFile
    INTEGER :: fileId, vlistId
    CHARACTER(:), POINTER :: cdiErrorText 
    CHARACTER(*), PARAMETER :: routine = modname//":readRestartAttributeFile"
    LOGICAL :: isReader

    isReader = 0 == p_comm_rank(p_comm_work)
    IF(isReader) THEN
      fileId  = streamOpenRead(attributeFile)
      ! check if the file could be opened
      IF(fileId < 0) THEN
        cdiErrorText => toCharacter(cdiStringError(fileId))
        CALL finish(routine, 'File '//attributeFile//' cannot be opened: '//cdiErrorText)
      END IF
      vlistId = streamInqVlist(fileId)
    END IF
    CALL restartAttributeList_read(vlistId, 0, p_comm_work)
    IF(isReader) CALL streamClose(fileId)
    CALL message(routine, "read namelists and attributes from restart file")
  END SUBROUTINE readRestartAttributeFile

  ! Reads attributes and namelists for all available domains from restart file.
  SUBROUTINE read_restart_header(modelType)
    CHARACTER(LEN=*), INTENT(IN) :: modelType
    CHARACTER(:), ALLOCATABLE :: filename, mfaname
    LOGICAL :: lIsMultifile
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":read_restart_header"

    ! Must NOT USE a timer here as the timers are NOT initialized yet,
    CALL findRestartFile(modelType, lIsMultifile, filename)
    IF(lIsMultifile) THEN
      mfaname = filename//"/attributes.nc"
      CALL readRestartAttributeFile(mfaname)
      IF(my_process_is_stdio()) CALL multifileCheckRestartFiles(filename)
    ELSE
      CALL readRestartAttributeFile(filename)
      IF(my_process_is_stdio()) CALL singlefileCheckRestartFiles(modelType)
    END IF
  END SUBROUTINE read_restart_header

  SUBROUTINE read_restart_files(p_patch, opt_ndom)
    TYPE(t_patch), INTENT(in) :: p_patch
    INTEGER, OPTIONAL, INTENT(in) :: opt_ndom
    CHARACTER(:), ALLOCATABLE :: restartPath
    LOGICAL :: lIsMultifileRestart, lMultifileTimersInitialized
    TYPE mtype_ll
      CHARACTER(:), ALLOCATABLE :: a
      TYPE(mtype_ll), POINTER :: next => NULL()
    END TYPE mtype_ll
    TYPE(mtype_ll), POINTER :: cur_mType, entry_mType
    INTEGER :: ndom_deopt
    TYPE(t_p_var_list_element), ALLOCATABLE :: varData(:)

    IF(timers_level >= 5) CALL timer_start(timer_load_restart)

    IF (ocean_initFromRestart_OVERRIDE) CALL read_restart_header(TRIM(get_my_process_name()))
    ! Make sure that all the subcounters are recognized as subcounters on all work processes.
    IF(timers_level >= 7) THEN
        CALL timer_start(timer_load_restart_io)
        CALL timer_stop(timer_load_restart_io)
        CALL timer_start(timer_load_restart_communication)
        CALL timer_stop(timer_load_restart_communication)
    END IF
    lMultifileTimersInitialized = .FALSE.

    IF (my_process_is_mpi_workroot()) &
         WRITE(0,'(a,i0)') "restart: reading restart data for patch ", p_patch%id
    ALLOCATE(entry_mType)
    CALL getModelTypes()
    cur_mType => entry_mType
    DO WHILE(ASSOCIATED(cur_mType%next))
        CALL createRestartVarData(varData, p_patch%id, cur_mType%next%a)
        !determine whether we have a multifile to READ
        CALL findRestartFile(cur_mType%next%a, lIsMultifileRestart, restartPath)
        IF(lIsMultifileRestart) THEN
            !multifile loading also uses these two timers
            IF(timers_level >= 7 .AND..NOT.lMultifileTimersInitialized) THEN
                CALL timer_start(timer_load_restart_comm_setup)
                CALL timer_stop(timer_load_restart_comm_setup)
                CALL timer_start(timer_load_restart_get_var_id)
                CALL timer_stop(timer_load_restart_get_var_id)
                lMultifileTimersInitialized = .TRUE.
            END IF

            !load the multifile
            CALL multifileReadPatch(varData, p_patch, restartPath)
        ELSE
            !can't USE restartPath here because the path IS patch dependent IN the single file CASE
            ndom_deopt = -1
            IF (PRESENT(opt_ndom)) ndom_deopt = opt_ndom
            CALL singlefileReadPatch(varData, cur_mType%next%a, p_patch, ndom_deopt)
        END IF

        DEALLOCATE(varData)
        entry_mType => cur_mType
        cur_mType => cur_mType%next
        DEALLOCATE(entry_mType)
    END DO

    IF(timers_level >= 5) CALL timer_stop(timer_load_restart)

    CALL message('','')
    CALL message('',separator)
    CALL message('','')
  CONTAINS

    SUBROUTINE getModelTypes()
      INTEGER :: lm
      LOGICAL :: skip

      DO WHILE(vl_iter%next())
        IF (vl_iter%cur%p%patch_id .NE. p_patch%id) CYCLE
        lm = LEN_TRIM(vl_iter%cur%p%model_type)
        cur_mType => entry_mType
        skip = .FALSE.
        DO WHILE(ASSOCIATED(cur_mType%next))
          skip = vl_iter%cur%p%model_type(1:lm) == cur_mType%next%a
          IF (skip) EXIT
          cur_mType => cur_mType%next
        END DO
        IF (skip) CYCLE
        ALLOCATE(cur_mType%next)
        cur_mType%next%a = vl_iter%cur%p%model_type(1:lm)
      END DO
    END SUBROUTINE getModelTypes
  END SUBROUTINE read_restart_files

END MODULE mo_load_restart
