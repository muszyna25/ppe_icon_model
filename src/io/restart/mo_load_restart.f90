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
    USE mo_cdi,                ONLY: streamOpenRead, streamInqVlist, streamClose, streamOpenRead, &
      &                              streamReadVarSlice, streamReadVarSliceF,                     &
      &                              vlistInqTaxis, vlistNvars,                                   &
      &                              vlistInqVarName, vlistInqVarGrid, vlistInqVarZaxis,          &
      &                              taxisInqVdate, taxisInqVtime, zaxisInqType, zaxisInqSize,    &
      &                              gridInqSize, ZAXIS_SURFACE, cdiStringError
    USE mo_fortran_tools,      ONLY: t_alloc_character
    USE mo_dictionary,         ONLY: t_dictionary, dict_size, dict_getKey, dict_set, dict_init,   &
      &                              dict_finalize
    USE mo_communication,      ONLY: t_ScatterPattern
    USE mo_exception,          ONLY: message, finish, warning
    USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, SINGLE_T, REAL_T, INT_T, SUCCESS
    USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT,              &
      &                              GRID_UNSTRUCTURED_EDGE
    USE mo_load_multifile_restart, ONLY: multifileReadPatch, multifileCheckRestartFiles
    USE mo_load_singlefile_restart, ONLY: singlefileReadPatch, singlefileCheckRestartFiles
    USE mo_kind,               ONLY: dp, sp
    USE mo_linked_list,        ONLY: t_list_element
    USE mo_model_domain,       ONLY: t_patch
    USE mo_mpi,                ONLY: p_comm_work, p_comm_rank, my_process_is_mpi_workroot, my_process_is_stdio
    USE mo_multifile_restart_util, ONLY: multifileRestartLinkName, multifileAttributesPath
    USE mo_restart_attributes, ONLY: t_RestartAttributeList, RestartAttributeList_make, &
      & setAttributesForRestarting, getAttributesForRestarting, ocean_initFromRestart_OVERRIDE
    USE mo_restart_namelist,   ONLY: t_NamelistArchive, namelistArchive
    USE mo_restart_util,       ONLY: restartSymlinkName
    USE mo_restart_var_data,   ONLY: t_restartVarData, createRestartVarData
    USE mo_timer,              ONLY: timer_start, timer_stop, timer_load_restart, timer_load_restart_io, &
                                   & timer_load_restart_comm_setup, timer_load_restart_communication, &
                                   & timer_load_restart_get_var_id, timers_level
    USE mo_util_string,        ONLY: int2string, separator, toCharacter
    USE mo_var_list,           ONLY: nvar_lists, var_lists

    IMPLICIT NONE

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
        CALL restartSymlinkName(modelType, 1, singlefileLinkName)
        CALL multifileRestartLinkName(modelType, multifileLinkName)

        !check whether the respective files exist
        INQUIRE(file = singlefileLinkName, exist = haveSinglefileLink)
#if defined (__INTEL_COMPILER)
        INQUIRE(directory = multifileLinkName, exist = haveMultifileLink)
#else
        INQUIRE(file = multifileLinkName, exist = haveMultifileLink)
#endif

        !determine which path to USE
        IF(haveSinglefileLink) THEN
            IF(haveMultifileLink) THEN
                CALL warning(routine, "both a singlefile and a multifile restart file are present, using the multifile restart")
                cache_isMultifile = .TRUE.
            ELSE
                cache_isMultifile = .FALSE.
            END IF
        ELSE
            IF(haveMultifileLink) THEN
                cache_isMultifile = .TRUE.
            ELSE
              WRITE (0,*) "singlefileLinkName = ", singlefileLinkName
              WRITE (0,*) "multifileLinkName = ", multifileLinkName
                CALL finish(routine, "fatal error: could not locate the restart symlink to use")
            END IF
        END IF

        cache_isValid = .TRUE.
    END IF

    out_lIsMultifile = cache_isMultifile
    IF(cache_isMultifile) THEN
        CALL multifileRestartLinkName(modelType, linkname)
    ELSE
        CALL restartSymlinkName(modelType, 1, linkname)
    END IF
  END SUBROUTINE findRestartFile

  ! Read all namelists used in the previous run
  ! and store them in a buffer. These values will overwrite the
  ! model default, and will later be overwritten if the user has
  ! specified something different for this integration.
  !
  ! Note: We read the namelists AND attributes only once and assume that these
  !       are identical for all domains (which IS guaranteed by the way the restart files are written).
  !
  ! Collective across p_comm_work.
  !
  ! This IS used for both single- AND multifile restart files.
  SUBROUTINE readRestartAttributeFile(attributeFile)
    CHARACTER(*), INTENT(IN) :: attributeFile

    INTEGER :: myRank, fileId, vlistId
    TYPE(t_NamelistArchive), POINTER :: namelists
    CHARACTER(:), POINTER :: cdiErrorText
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    CHARACTER(*), PARAMETER :: routine = modname//":readRestartAttributeFile"

    myRank = p_comm_rank(p_comm_work)

    IF(.NOT.ocean_initFromRestart_OVERRIDE) &
! no namelists needed in case of initialize_FromRestart=true
! only the other attributes are actually of interest...
      & namelists => namelistArchive()
    IF(myRank == 0) THEN
        fileId  = streamOpenRead(attributeFile)
        ! check if the file could be opened
        IF(fileId < 0) THEN
            cdiErrorText => toCharacter(cdiStringError(fileId))
            CALL finish(routine, 'File '//attributeFile//' cannot be opened: '//cdiErrorText)
        END IF
        vlistId = streamInqVlist(fileId)
        IF(.NOT.ocean_initFromRestart_OVERRIDE) &
          & CALL namelists%readFromFile(vlistId)
    END IF
    IF(.NOT.ocean_initFromRestart_OVERRIDE) &
      & CALL namelists%bcast(0, p_comm_work)
    restartAttributes => RestartAttributeList_make(vlistId, 0, p_comm_work)
    CALL setAttributesForRestarting(restartAttributes)
    IF(myRank == 0) THEN
        CALL streamClose(fileId)
        WRITE(0,*) "restart: read namelists and attributes from restart file"
    END IF
  END SUBROUTINE readRestartAttributeFile

  ! Reads attributes and namelists for all available domains from restart file.
  SUBROUTINE read_restart_header(modelType)
    CHARACTER(LEN=*), INTENT(IN) :: modelType
    CHARACTER(:), ALLOCATABLE :: filename, mfaname
    LOGICAL :: lIsMultifile
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":read_restart_header"

    ! Must NOT USE a timer here as the timers are NOT initialized yet,
    ! AND I dare NOT to make such a significant order change within construct_atmo_model()
    ! for fear of silently breaking the initialization of some other global variable.
!   IF(timers_level >= 5) CALL timer_start(timer_load_restart)

    CALL findRestartFile(modelType, lIsMultifile, filename)
    IF(lIsMultifile) THEN
        CALL multifileAttributesPath(filename, mfaname)
        CALL readRestartAttributeFile(mfaname)
        IF(my_process_is_stdio()) CALL multifileCheckRestartFiles(filename)
    ELSE
        CALL readRestartAttributeFile(filename)
        IF(my_process_is_stdio()) CALL singlefileCheckRestartFiles(modelType)
    END IF

    ! See above.
!   IF(timers_level >= 5) CALL timer_stop(timer_load_restart)
  END SUBROUTINE read_restart_header

  !returns an array with all the different model types that are PRESENT IN var_lists
  SUBROUTINE getModelTypes(domain, resultVar) 
    TYPE(t_alloc_character), ALLOCATABLE, INTENT(INOUT) :: resultVar(:)
    INTEGER, INTENT(IN) :: domain

    TYPE(t_dictionary) :: modelTypes    !this dictionary IS actually used as a set, so the values are empty strings
    INTEGER :: i, modelTypeCount, error
    CHARACTER(*), PARAMETER :: routine = modname//":getModelTypes"

    CALL dict_init(modelTypes, .TRUE.)

    !insert all the model_type strings as keys into the dictionary
    DO i = 1, nvar_lists
        IF(var_lists(i)%p%patch_id == domain) THEN
          CALL dict_set(modelTypes, var_lists(i)%p%model_type, "")
        END IF
    END DO

    !convert the RESULT into an array
    modelTypeCount = dict_size(modelTypes)
    ALLOCATE(resultVar(modelTypeCount), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    DO i = 1, modelTypeCount
        resultVar(i)%a = TRIM(dict_getKey(modelTypes, i))
    END DO

    CALL dict_finalize(modelTypes)
  END SUBROUTINE getModelTypes

  SUBROUTINE read_restart_files(p_patch, opt_ndom)
    TYPE(t_patch), INTENT(in) :: p_patch
    INTEGER, OPTIONAL, INTENT(in) :: opt_ndom
    CHARACTER(:), ALLOCATABLE :: restartPath
    LOGICAL :: lIsMultifileRestart, lMultifileTimersInitialized
    TYPE(t_alloc_character), ALLOCATABLE :: modelTypes(:)
    INTEGER :: i, integerTrash
    TYPE(t_restartVarData), POINTER :: varData(:)

    IF(timers_level >= 5) CALL timer_start(timer_load_restart)

    IF (ocean_initFromRestart_OVERRIDE) CALL read_restart_header("oce")
    ! Make sure that all the subcounters are recognized as subcounters on all work processes.
    IF(timers_level >= 7) THEN
        CALL timer_start(timer_load_restart_io)
        CALL timer_stop(timer_load_restart_io)
        CALL timer_start(timer_load_restart_communication)
        CALL timer_stop(timer_load_restart_communication)
    END IF
    lMultifileTimersInitialized = .FALSE.

    IF(my_process_is_mpi_workroot()) THEN
        WRITE(0,*) "restart: reading restart data for patch "//TRIM(int2string(p_patch%id))// &
                 & ", nvar_lists = "//TRIM(int2string(nvar_lists))
    END IF

    ! get the list of all the different model types that we need to
    ! consider:
    CALL getModelTypes(p_patch%id, modelTypes)

    DO i = 1, SIZE(modelTypes)
        varData => createRestartVarData(p_patch%id, modelTypes(i)%a, integerTrash)

        !determine whether we have a multifile to READ
        CALL findRestartFile(modelTypes(i)%a, lIsMultifileRestart, restartPath)
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
            CALL singlefileReadPatch(varData, modelTypes(i)%a, p_patch, opt_ndom)
        END IF

        DEALLOCATE(varData)
    END DO

    IF(timers_level >= 5) CALL timer_stop(timer_load_restart)

    CALL message('','')
    CALL message('',separator)
    CALL message('','')
  END SUBROUTINE read_restart_files

END MODULE mo_load_restart
