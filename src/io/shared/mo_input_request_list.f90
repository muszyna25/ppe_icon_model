!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

! A list of requests for input.
!
! The use case for which this class has been developed is this:
!  1. A list is created with `myList = t_InputRequestList()`, and the names (only the names) of the requested fields are added with `request()`.
!  2. Files are read with readFile(). This reads all DATA ASSOCIATED with the requested variable names into memory (already distributing it to the worker PEs to keep memory footprint down).
!  3. User code inspects AND retrieves the fetched DATA.

MODULE mo_input_request_list
    USE ISO_C_BINDING, ONLY: C_CHAR, C_SIGNED_CHAR, C_INT, C_LONG, C_DOUBLE, C_NULL_PTR, C_NULL_CHAR, C_ASSOCIATED

    USE mo_cdi, ONLY: t_CdiIterator, cdiIterator_new, cdiIterator_nextField, cdiIterator_delete, cdiIterator_inqVTime, &
                    & cdiIterator_inqLevelType, cdiIterator_inqLevel, cdiIterator_inqGridId, cdiIterator_inqVariableName, &
                    & gridInqType, gridInqSize, CDI_UNDEFID, ZAXIS_SURFACE, ZAXIS_GENERIC, ZAXIS_HYBRID, ZAXIS_HYBRID_HALF, &
                    & ZAXIS_PRESSURE, ZAXIS_HEIGHT, ZAXIS_DEPTH_BELOW_SEA, ZAXIS_DEPTH_BELOW_LAND, ZAXIS_ISENTROPIC, &
                    & ZAXIS_TRAJECTORY, ZAXIS_ALTITUDE, ZAXIS_SIGMA, ZAXIS_MEANSEA, ZAXIS_TOA, ZAXIS_SEA_BOTTOM, &
                    & ZAXIS_ATMOSPHERE, ZAXIS_CLOUD_BASE, ZAXIS_CLOUD_TOP, ZAXIS_ISOTHERM_ZERO, ZAXIS_SNOW, ZAXIS_LAKE_BOTTOM, &
                    & ZAXIS_SEDIMENT_BOTTOM, ZAXIS_SEDIMENT_BOTTOM_TA, ZAXIS_SEDIMENT_BOTTOM_TW, ZAXIS_MIX_LAYER, &
                    & ZAXIS_REFERENCE, cdiIterator_inqTile, CDI_NOERR, CDI_EINVAL, GRID_UNSTRUCTURED, t_CdiParam, &
                    & cdiIterator_inqParamParts, gridInqNumber, gridInqPosition, cdiGribIterator_inqLongValue, t_CdiGribIterator, &
                    & cdiGribIterator_clone, cdiGribIterator_delete, cdiIterator_inqRTime, CDI_UUID_SIZE
    USE mo_communication, ONLY: t_ScatterPattern
    USE mo_dictionary, ONLY: t_dictionary, dict_copy, dict_init, dict_get, dict_finalize
    USE mo_exception, ONLY: message, finish
    USE mo_grid_config, ONLY: n_dom
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_initicon_config, ONLY: timeshift, lconsistency_checks
    USE mo_initicon_utils, ONLY: initicon_inverse_post_op
    USE mo_input_container, ONLY: t_InputContainer, inputContainer_make
    USE mo_kind, ONLY: wp, dp
    USE mo_lnd_nwp_config, ONLY: t_tile, select_tile, get_tile_suffix, find_tile_id
    USE mo_model_domain, ONLY: t_patch
    USE mo_mpi, ONLY: my_process_is_mpi_workroot, get_my_mpi_work_id, p_bcast, process_mpi_root_id, p_comm_work, &
                    & my_process_is_stdio, p_pe, p_isEqual
    USE mo_run_config, ONLY: msg_level
    USE mo_time_config, ONLY: time_config
    USE mo_util_cdi, ONLY: trivial_tileId
    USE mo_util_string, ONLY: int2string, toCharArray, toCharacter, charArray_equal, charArray_toLower, charArray_dup, one_of
    USE mo_util_table, ONLY: t_table, initialize_table, add_table_column, set_table_entry, print_table, finalize_table
    USE mo_util_uuid, ONLY: t_uuid, uuid_string_length, uuid_unparse, OPERATOR(==)
    USE mtime, ONLY: datetime, timedelta, newDatetime, newTimedelta, timedeltaToString, deallocateDatetime, deallocateTimedelta, &
                   & max_timedelta_str_len, OPERATOR(-), OPERATOR(+), OPERATOR(==)

    IMPLICIT NONE

PUBLIC :: t_InputRequestList, InputRequestList_create

    TYPE :: t_InputRequestList
        PRIVATE
        TYPE(t_ListEntry), POINTER :: list(:)
        INTEGER :: variableCount

    CONTAINS
        PROCEDURE :: request => InputRequestList_request    !< Require that a variable be read.
        PROCEDURE :: requestMultiple => InputRequestList_requestMultiple    !< Require that a list of variables. Unlike request() this will request the trimmed strings (because it's impossible to pass an array of strings of different LEN).
        PROCEDURE :: readFile => InputRequestList_readFile  !< Scan a file for input data to satisfy the requests.

        PROCEDURE :: getLevels => InputRequestList_getLevels    !< Get the count AND height values (elevation/presure/whatever) of all the levels encountered IN the file.

        !> The `fetchXXX()` methods simply RETURN FALSE IF the DATA could NOT be fetched entirely.
        !> This IS an atomic operation: Either the entire output array IS overwritten OR it IS NOT touched at all.
        PROCEDURE :: fetch2d => InputRequestList_fetch2d
        PROCEDURE :: fetch3d => InputRequestList_fetch3d
        PROCEDURE :: fetchSurface => InputRequestList_fetchSurface  !No level given, fail IF there are several levels.
        PROCEDURE :: fetchTiled2d => InputRequestList_fetchTiled2d
        PROCEDURE :: fetchTiled3d => InputRequestList_fetchTiled3d
        PROCEDURE :: fetchTiledSurface => InputRequestList_fetchTiledSurface  !No level given, fail IF there are several levels.

        !> The `fetchRequiredXXX()` methods will CALL `finish()` IF there are holes IN the DATA that was READ.
        PROCEDURE :: fetchRequired2d => InputRequestList_fetchRequired2d
        PROCEDURE :: fetchRequired3d => InputRequestList_fetchRequired3d
        PROCEDURE :: fetchRequiredSurface => InputRequestList_fetchRequiredSurface
        PROCEDURE :: fetchRequiredTiled2d => InputRequestList_fetchRequiredTiled2d
        PROCEDURE :: fetchRequiredTiled3d => InputRequestList_fetchRequiredTiled3d
        PROCEDURE :: fetchRequiredTiledSurface => InputRequestList_fetchRequiredTiledSurface

        PROCEDURE :: printInventory => InputRequestList_printInventory
        PROCEDURE :: checkRuntypeAndUuids => InputRequestList_checkRuntypeAndUuids

        PROCEDURE :: destruct => InputRequestList_destruct  !< Destructor.

        PROCEDURE, PRIVATE :: checkRequests => InputRequestList_checkRequests   !< Check that all processes IN the communicator have the same view on which variables are needed.
        PROCEDURE, PRIVATE :: translateNames => InputRequestList_translateNames !< Recalculates the translatedNames of all list entries using the given dictionary.
        PROCEDURE, PRIVATE :: findIconName => InputRequestList_findIconName    !< Retrieve a t_ListEntry for the given ICON variable name if it exists already.
        PROCEDURE, PRIVATE :: findTranslatedName => InputRequestList_findTranslatedName    !< As findIconName, but uses the translatedVarName.
        PROCEDURE, PRIVATE :: sendStopMessage => InputRequestList_sendStopMessage
        PROCEDURE, PRIVATE :: sendFieldMetadata => InputRequestList_sendFieldMetadata
        PROCEDURE, PRIVATE :: recieveFieldMetadata => InputRequestList_recieveFieldMetadata
        PROCEDURE, PRIVATE :: isRecordValid => InputRequestList_isRecordValid
        PROCEDURE, PRIVATE :: nextField => InputRequestList_nextField
    END TYPE

PRIVATE

    ! These objects are created via findDomainData(, , opt_lcreate = .TRUE.), which will already instanciate an empty container.
    ! On the I/O PE, it IS the job of InputRequestList_isRecordValid() to immediately add a MetadataCache, so that ANY DomainData object returned by findDomainData() CONTAINS both a valid InputContainer AND a valid MetadataCache.
    TYPE :: t_DomainData
        INTEGER :: domain
        TYPE(t_DomainData), POINTER :: next

        CLASS(t_InputContainer), POINTER :: container
        TYPE(t_MetadataCache), POINTER :: metadata  !< Some metadata connected with the variable, which IS only used for consistency checking AND printing of the inventory table.
    END TYPE

    TYPE :: t_ListEntry
        CHARACTER(KIND = C_CHAR), POINTER :: iconVarName(:)  !< The NAME as it has been requested.
        CHARACTER(KIND = C_CHAR), POINTER :: translatedVarName(:)   !< The NAME as it IS matched against the stored NAME. This has dictionary translation, trimming, AND CASE canonization applied.
        TYPE(t_DomainData), POINTER :: domainData   !< A linked list of an InputContainer AND a MetadataCache for each domain. Only accessed via findDomainData().
    END TYPE

    TYPE :: t_MetadataCache
        CHARACTER(KIND = C_CHAR), POINTER :: rtime(:), vtime(:)
        INTEGER :: levelType, gridNumber, gridPosition, runClass, experimentId, generatingProcessType
        TYPE(t_CdiParam) :: param
        TYPE(t_uuid) :: gridUuid

    CONTAINS
        PROCEDURE :: equalTo => MetadataCache_equalTo
    END TYPE

    CHARACTER(*), PARAMETER :: modname = "mo_input_request_list"
    LOGICAL, PARAMETER :: debugModule = .FALSE.

CONTAINS

    !Can't use a type constructor interface t_InputRequestList() since the cray compiler looses the list pointer while returning + assigning the function result.
    FUNCTION InputRequestList_create() RESULT(result)
        TYPE(t_InputRequestList), POINTER :: result

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_create"
        INTEGER :: error, i

        ALLOCATE(result, STAT = error)
        if(error /= SUCCESS) CALL finish(routine, "error allocating memory")

        result%variableCount = 0
        ALLOCATE(result%list(8), STAT = error)
        if(error /= SUCCESS) CALL finish(routine, "error allocating memory")
        DO i = 1, SIZE(RESULT%list, 1)
            RESULT%list(i)%iconVarName => NULL()
            RESULT%list(i)%translatedVarName => NULL()
            RESULT%list(i)%domainData => NULL()
        END DO
    END FUNCTION InputRequestList_create

    FUNCTION findDomainData(listEntry, domain, opt_lcreate) RESULT(RESULT)
        TYPE(t_ListEntry), POINTER, INTENT(INOUT) :: listEntry
        INTEGER, VALUE :: domain
        LOGICAL, OPTIONAL, VALUE :: opt_lcreate
        TYPE(t_DomainData), POINTER :: RESULT

        CHARACTER(*), PARAMETER :: routine = modname//":findDomainData"
        INTEGER :: error

        IF(.NOT.ASSOCIATED(listEntry)) CALL finish(routine, "assertion failed, listEntry IS NOT ASSOCIATED")

        ! Try to find a preexisting DomainData object.
        RESULT => listEntry%domainData
        DO
            IF(.NOT.ASSOCIATED(RESULT)) EXIT
            IF(RESULT%domain == domain) RETURN
            RESULT => RESULT%next
        END DO

        ! Nothing preexisting found, should we create a new one?
        if(PRESENT(opt_lcreate)) THEN
            IF(opt_lcreate) THEN
                ALLOCATE(RESULT, STAT = error)
                IF(error /= SUCCESS) CALL finish(routine, "error allocating memory")
                RESULT%domain = domain
                RESULT%next => listEntry%domainData
                RESULT%container => InputContainer_make()
                RESULT%metadata => NULL()
                listEntry%domainData => RESULT
            END IF
        END IF
    END FUNCTION findDomainData

    SUBROUTINE InputRequestList_translateNames(me, opt_dict)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        TYPE(t_dictionary), OPTIONAL, INTENT(IN) :: opt_dict

        INTEGER :: i
        CHARACTER(:), POINTER :: tempName

        DO i = 1, me%variableCount
            IF(ASSOCIATED(me%list(i)%translatedVarName)) DEALLOCATE(me%list(i)%translatedVarName)
            IF(PRESENT(opt_dict)) THEN
                tempName => toCharacter(me%list(i)%iconVarName)
                me%list(i)%translatedVarName => toCharArray(TRIM(dict_get(opt_dict, tempName, tempName)))
                DEALLOCATE(tempName)
            ELSE
                me%list(i)%translatedVarName => charArray_dup(me%list(i)%iconVarName)
            END IF
            CALL charArray_toLower(me%list(i)%translatedVarName)
        END DO
    END SUBROUTINE InputRequestList_translateNames

    FUNCTION InputRequestList_findIconName(me, fieldName, opt_lDebug) RESULT(RESULT)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: fieldName
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug
        TYPE(t_ListEntry), POINTER :: result

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_findIconName"
        INTEGER :: i
        LOGICAL :: debugInfo
        CHARACTER(:), POINTER :: tempName

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        result => NULL()
        DO i = 1, me%variableCount
            tempName => toCharacter(me%list(i)%iconVarName)
            IF(fieldName == tempName) THEN
                IF(debugInfo) CALL message(routine, fieldName//" == "//tempName)
                result => me%list(i)
                RETURN
            ELSE
                IF(debugInfo) CALL message(routine, fieldName//" /= "//tempName)
            END IF
            DEALLOCATE(tempName)
        END DO
    END FUNCTION InputRequestList_findIconName

    FUNCTION InputRequestList_findTranslatedName(me, fieldName) RESULT(result)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(KIND = C_CHAR), INTENT(INOUT) :: fieldName(:)
        TYPE(t_ListEntry), POINTER :: result

        INTEGER :: i

        CALL charArray_toLower(fieldName)
        result => NULL()
        DO i = 1, me%variableCount
            IF(charArray_equal(fieldName, me%list(i)%translatedVarName)) THEN
                result => me%list(i)
                RETURN
            END IF
        END DO
    END FUNCTION InputRequestList_findTranslatedName

    !XXX: This also ensures that the requests have been given in the same order on all processes. Not technically necessary, but easier to implement and I guess, if the order is not the same, that's a hint that there is a bug somewhere else.
    SUBROUTINE InputRequestList_checkRequests(me)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_checkRequests"
        INTEGER :: i, j, error, concatenatedSize, curSize, accumulatedSize
        CHARACTER(KIND = C_CHAR), ALLOCATABLE :: concatenatedNames(:)

        !compute the concatenation of all requested variables
        concatenatedSize = 0
        DO i = 1, me%variableCount
            concatenatedSize = concatenatedSize + SIZE(me%list(i)%iconVarName) + 1
        END DO
        ALLOCATE(concatenatedNames(concatenatedSize), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        accumulatedSize = 0
        DO i = 1, me%variableCount
            curSize = SIZE(me%list(i)%iconVarName)
            DO j = 1, curSize
                concatenatedNames(accumulatedSize + j) = me%list(i)%iconVarName(j)
            END DO
            concatenatedNames(accumulatedSize + curSize + 1) = C_NULL_CHAR
            accumulatedSize = accumulatedSize + curSize + 1
        END DO

        !check that all processes have the same concatenatedNames string
        IF(.NOT. p_isEqual(concatenatedNames, p_comm_work)) THEN
            print*, "process ", p_pe, " has the variable list: ", concatenatedNames
            CALL finish(routine, "not all processes have the same requests in their t_InputRequestList")
        END IF
    END SUBROUTINE InputRequestList_checkRequests

    SUBROUTINE InputRequestList_request(me, fieldName)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: fieldName

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_request"
        INTEGER :: i, listSize, error, jg
        TYPE(t_ListEntry), POINTER :: tempList(:), newEntry

        ! don't add a name twice
        IF(ASSOCIATED(me%findIconName(fieldName))) RETURN

        IF(debugModule .AND. my_process_is_mpi_workroot()) print*, 'Adding request for variable "', fieldName, '"'

        ! ensure space for the new container
        listSize = SIZE(me%list, 1)
        IF(me%variableCount == listSize) THEN
            ALLOCATE(tempList(2*listSize), STAT = error)
            if(error /= SUCCESS) CALL finish(routine, "error allocating memory")
            DO i = 1, listSize
                tempList(i) = me%list(i)
            END DO
            DO i = listSize + 1, 2*listSize
                tempList(i)%iconVarName => NULL()
                tempList(i)%translatedVarName => NULL()
                tempList(i)%domainData => NULL()
            END DO
            DEALLOCATE(me%list)
            me%list => tempList
        END IF

        ! add the entry to our list
        me%variableCount = me%variableCount + 1
        newEntry => me%list(me%variableCount)

        newEntry%iconVarName => toCharArray(fieldName)
        newEntry%translatedVarName => NULL()
    END SUBROUTINE InputRequestList_request

    SUBROUTINE InputRequestList_requestMultiple(me, fieldNames)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: fieldNames(:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_requestMultiple"
        INTEGER :: i

        DO i = 1, SIZE(fieldNames, 1)
            CALL me%request(TRIM(fieldNames(i)))
        END DO
    END SUBROUTINE InputRequestList_requestMultiple

    LOGICAL FUNCTION InputRequestList_isRecordValid(me, iterator, p_patch, level, tileId, variableName, lIsFg) RESULT(result)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        TYPE(t_CdiIterator) :: iterator
        TYPE(t_patch), INTENT(IN) :: p_patch
        REAL(dp), INTENT(OUT) :: level
        INTEGER, INTENT(OUT) :: tileId
        CHARACTER(KIND = C_CHAR), DIMENSION(:), POINTER, INTENT(INOUT) :: variableName
        LOGICAL, VALUE :: lIsFg

        INTEGER(KIND = C_INT) :: error, gridId, gridType, gridSize, tileIndex, tileAttribute
        REAL(KIND = C_DOUBLE) :: levelValue
        TYPE(t_MetadataCache), POINTER :: metadata
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        TYPE(t_CdiGribIterator) :: gribIterator
        TYPE(datetime) :: iniTime, startTime
        TYPE(datetime), POINTER :: tempTime
        INTEGER(KIND = C_SIGNED_CHAR) :: gridUuid(CDI_UUID_SIZE)
        CHARACTER(:), POINTER :: vtimeString
        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_isRecordValid"

        result = .TRUE.
        variableName => cdiIterator_inqVariableName(iterator)

        metadata => MetadataCache_create()
        CALL cdiIterator_inqParamParts(iterator, metadata%param%discipline, metadata%param%category, metadata%param%number)

        !Check the time.
        metadata%vtime => cdiIterator_inqVTime(iterator)
        metadata%rtime => cdiIterator_inqRTime(iterator)
        IF (.NOT. ASSOCIATED(metadata%vtime) .OR. .NOT. ASSOCIATED(metadata%rtime)) THEN
          CALL finish(routine, "Internal error!")
        END IF

        IF(lconsistency_checks) THEN
            vtimeString => toCharacter(metadata%vtime)
            tempTime => newDatetime(vtimeString)
            DEALLOCATE(vtimeString)
            IF(lIsFg) THEN
                iniTime = time_config%tc_startdate
                ! add timeshift to INI-datetime to get true starting time
                startTime = iniTime + timeshift%mtime_shift
                IF(.NOT.(tempTime == startTime)) CALL fail("vtime of first-guess field does not match model start time")
            ELSE
                iniTime = time_config%tc_startdate
                IF(.NOT.(tempTime == iniTime)) CALL fail("vtime of analysis field does not match model initialization time")
            END IF
            CALL deallocateDatetime(tempTime)
        END IF

        !We only check the primary (top) level (selector = 0). Usually that's the only one, but GRIB does allow a secondary lower boundary level.
        metadata%levelType = cdiIterator_inqLevelType(iterator, 0)
        SELECT CASE(metadata%levelType)
            !the level types that translate to a single height VALUE
            CASE(ZAXIS_SURFACE, ZAXIS_PRESSURE, ZAXIS_HEIGHT, ZAXIS_DEPTH_BELOW_SEA, ZAXIS_DEPTH_BELOW_LAND, ZAXIS_ALTITUDE, &
                &ZAXIS_REFERENCE, ZAXIS_SNOW)

                levelValue = 0.0
                error = cdiIterator_inqLevel(iterator, 1, outValue1 = levelValue)
                level = REAL(levelValue, dp)
                IF(error /= 0) CALL fail("cdiIterator_inqLevel() failed")
                !TODO[NH]: check the zaxis UUID

            !the level types for special levels
            CASE(ZAXIS_TOA, ZAXIS_ATMOSPHERE, ZAXIS_CLOUD_BASE, ZAXIS_CLOUD_TOP, ZAXIS_ISOTHERM_ZERO, &
                &ZAXIS_MEANSEA, ZAXIS_SEA_BOTTOM, ZAXIS_LAKE_BOTTOM, ZAXIS_SEDIMENT_BOTTOM, ZAXIS_SEDIMENT_BOTTOM_TA, &
                &ZAXIS_SEDIMENT_BOTTOM_TW, ZAXIS_MIX_LAYER)

                level = REAL(-metadata%levelType, dp)

            !the known z-axis types that are NOT handled by this code
            CASE(ZAXIS_GENERIC)
                CALL fail("z-axis type ZAXIS_GENERIC is not implemented")
            CASE(ZAXIS_HYBRID)
                CALL fail("z-axis type ZAXIS_HYBRID is not implemented")
            CASE(ZAXIS_HYBRID_HALF)
                CALL fail("z-axis type ZAXIS_HYBRID_HALF is not implemented")
            CASE(ZAXIS_ISENTROPIC)
                CALL fail("z-axis type ZAXIS_ISENTROPIC is not implemented")
            CASE(ZAXIS_TRAJECTORY)
                CALL fail("z-axis type ZAXIS_TRAJECTORY is not implemented")
            CASE(ZAXIS_SIGMA)
                CALL fail("z-axis type ZAXIS_SIGMA is not implemented")

            !fallback to catch future expansions of the list of available z-axis types
            CASE DEFAULT
                CALL fail("unknown z-axis TYPE ("//int2string(metadata%levelType)//")")
        END SELECT

        !Check the grid.
        gridId = cdiIterator_inqGridId(iterator)
        IF(gridId /= CDI_UNDEFID) THEN
            !XXX: I believe, it's enough sanity checking if we check the type and the size of the grid,
            !     I don't want to go into checking the lon/lat for all its vertices here...
            !     A test for the correct grid SIZE IS IMPLICIT IN the t_InputContainer when it selects the scatter pattern to USE.
            gridType = gridInqType(gridId)
            IF(gridType /= GRID_UNSTRUCTURED) THEN
                CALL fail("support for this gridtype is not implemented (CDI grid type = "//TRIM(int2string(gridType))//")")
            ELSE
                CALL gridInqUuid(gridId, metadata%gridUuid%DATA)
                metadata%gridNumber = gridInqNumber(gridID)
                metadata%gridPosition = gridInqPosition(gridID)
            END IF
        ELSE
            CALL fail("couldn't inquire grid ID")
        END IF

        error = cdiIterator_inqTile(iterator, tileIndex, tileAttribute);
        SELECT CASE(error)
            CASE(CDI_NOERR)
                tileId = find_tile_id(tileIndex, tileAttribute)

            CASE(CDI_EINVAL)
                !There IS no tile information connected to this field, so we USE the trivial tileId.
                tileId = trivial_tileId

            CASE DEFAULT
                CALL finish(routine, "unexpected error while reading tile information")
        END SELECT

        !Fetch some additional metadata.
        metadata%runClass = -1
        metadata%experimentId = -1
        metadata%generatingProcessType = -1
        IF(RESULT) THEN
            gribIterator = cdiGribIterator_clone(iterator)
            IF(C_ASSOCIATED(gribIterator%ptr)) THEN
                metadata%runClass = INT(cdiGribIterator_inqLongValue(gribIterator, "backgroundProcess"))
                metadata%experimentId = INT(cdiGribIterator_inqLongValue(gribIterator, "localNumberOfExperiment"))
                metadata%generatingProcessType = INT(cdiGribIterator_inqLongValue(gribIterator, "typeOfGeneratingProcess"))
                CALL cdiGribIterator_delete(gribIterator)
            END IF
        END IF

        !Check whether the metadata of this record IS consistent with the metadata we've already seen for this variable.
        listEntry => me%findTranslatedName(variableName)
        IF(.NOT.ASSOCIATED(listEntry)) THEN
            RESULT = .FALSE.    !We are NOT interested IN this variable.
            CALL MetadataCache_delete(metadata)
            RETURN
        END IF
        IF(RESULT) domainData => findDomainData(listEntry, p_patch%id)
        IF(RESULT .AND. ASSOCIATED(domainData)) RESULT = metadata%equalTo(domainData%metadata)

        !Commit AND cleanup.
        IF(RESULT) THEN
            IF(ASSOCIATED(domainData)) THEN
                CALL MetadataCache_delete(metadata)
            ELSE
                domainData => findDomainData(listEntry, p_patch%id, opt_lcreate = .TRUE.)
                domainData%metadata => metadata  !We don't have a metadata cache yet, so we just remember this one.
            END IF
        ELSE
            !The record was not valid.
            DEALLOCATE(variableName)
            variableName => NULL()
            CALL MetadataCache_delete(metadata)
        END IF

    CONTAINS

        SUBROUTINE fail(message)
            CHARACTER(LEN = *), INTENT(IN) :: message

            IF(msg_level >= 1) print*, 'invalid record for variable "', variableName, '" encountered: '//message
            RESULT = .FALSE.
        END SUBROUTINE fail

    END FUNCTION InputRequestList_isRecordValid

    ! Format of the messages which are broadcasted by the following three FUNCTION:
    !
    ! REAL(dp) :: message(3)
    ! message(1) = length of variable NAME, zero IF this is a stop message
    ! message(2) = level
    ! message(3) = tileId
    !
    ! In the case that the variable NAME length is nonzero, this is followed by another message containing the NAME itself.
    ! Note: The broadcasts IN recieveFieldMetadata() are matched with the broadcasts within sendFieldMetadata() and sendStopMessage().
    SUBROUTINE InputRequestList_sendStopMessage(me)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_sendStopMessage"
        REAL(dp) :: message(3)

        message(1) = 0.0_dp
        message(2) = 0.0_dp
        message(3) = 0.0_dp
        CALL p_bcast(message, process_mpi_root_id, p_comm_work)
    END SUBROUTINE InputRequestList_sendStopMessage

    SUBROUTINE InputRequestList_sendFieldMetadata(me, level, tileId, variableName)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        REAL(dp), VALUE :: level
        INTEGER, VALUE :: tileId
        CHARACTER(KIND = C_CHAR), DIMENSION(:), POINTER, INTENT(IN) :: variableName

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_sendFieldMetadata"
        REAL(dp) :: message(3)
        CHARACTER(:), POINTER :: tempName
        INTEGER :: error

        message(1) = REAL(SIZE(variableName, 1), dp)
        message(2) = level
        message(3) = REAL(tileId, dp)
        CALL p_bcast(message, process_mpi_root_id, p_comm_work)

        tempName => toCharacter(variableName)
        IF(debugModule) print*, 'Reading field for variable "'//tempName//'"'
        CALL p_bcast(tempName, process_mpi_root_id, p_comm_work)
        DEALLOCATE(tempName)
    END SUBROUTINE InputRequestList_sendFieldMetadata

    LOGICAL FUNCTION InputRequestList_recieveFieldMetadata(me, level, tileId, variableName) RESULT(RESULT)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        REAL(dp), INTENT(INOUT) :: level
        INTEGER, INTENT(INOUT) :: tileId
        CHARACTER(KIND = C_CHAR), DIMENSION(:), POINTER, INTENT(INOUT) :: variableName

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_recieveFieldMetadata"
        REAL(dp) :: message(3)
        CHARACTER(:), ALLOCATABLE :: tempName
        INTEGER :: error

        message(1) = 0.0_dp
        message(2) = 0.0_dp
        message(3) = 0.0_dp
        CALL p_bcast(message, process_mpi_root_id, p_comm_work)
        level = message(2)
        tileId = INT(message(3))
        RESULT = message(1) /= 0.0_dp
        variableName => NULL()
        IF(RESULT) THEN
            ALLOCATE(CHARACTER(LEN = INT(message(1))) :: tempName, STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "error allocating memory")
            CALL p_bcast(tempName, process_mpi_root_id, p_comm_work)
            variableName => toCharArray(tempName)
            DEALLOCATE(tempName)
        END IF
    END FUNCTION InputRequestList_recieveFieldMetadata

    ! Find the next field that we are interested IN.
    ! This FUNCTION is collective: either all processes RETURN .TRUE. or all RETURN .FALSE. .
    ! When this FUNCTION returns FALSE, THEN there is no further field IN the file.
    !
    ! ignoredRecords IS NOT reset by this FUNCTION, it IS ONLY incremented
    LOGICAL FUNCTION InputRequestList_nextField(me, iterator, p_patch, level, tileId, variableName, ignoredRecords, lIsFg) &
    &RESULT(RESULT)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        TYPE(t_CdiIterator), VALUE :: iterator
        TYPE(t_patch), INTENT(IN) :: p_patch
        REAL(dp), INTENT(OUT) :: level
        INTEGER, INTENT(OUT) :: tileId
        CHARACTER(KIND = C_CHAR), DIMENSION(:), POINTER, INTENT(INOUT) :: variableName
        INTEGER, INTENT(INOUT) :: ignoredRecords
        LOGICAL, VALUE :: lIsFg

        RESULT = .FALSE.
        IF(my_process_is_mpi_workroot()) THEN
            ! Scan the file until we find a field that we are interested in.
            DO
                IF(cdiIterator_nextField(iterator) /= 0) THEN
                    CALL me%sendStopMessage()
                    RETURN
                ELSE
                    IF(me%isRecordValid(iterator, p_patch, level, tileId, variableName, lIsFg)) THEN
                        IF(ASSOCIATED(me%findTranslatedName(variableName))) THEN
                            CALL me%sendFieldMetadata(level, tileId, variableName)
                            RESULT = .TRUE.
                            RETURN
                        END IF
                    ELSE
                        ignoredRecords = ignoredRecords + 1
                    END IF
                END IF
            END DO
        ELSE
            RESULT = me%recieveFieldMetadata(level, tileId, variableName)
        END IF
    END FUNCTION InputRequestList_nextField

    SUBROUTINE InputRequestList_readFile(me, p_patch, path, lIsFg, opt_dict)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        TYPE(t_patch), INTENT(IN) :: p_patch
        CHARACTER(LEN = *, KIND = C_CHAR), INTENT(IN) :: path
        TYPE(t_dictionary), OPTIONAL, INTENT(IN) :: opt_dict
        LOGICAL, VALUE :: lIsFg

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_readFile"
        TYPE(t_CdiIterator) :: iterator
        REAL(dp) :: level
        CHARACTER(KIND = C_CHAR), DIMENSION(:), POINTER :: vtime, variableName
        INTEGER :: i, tileId, recordsRead, recordsIgnored
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData

        recordsRead = 0
        recordsIgnored = 0

        CALL me%checkRequests() !sanity checks
        CALL me%translateNames(opt_dict)

        iterator%ptr = C_NULL_PTR
        IF(my_process_is_mpi_workroot()) THEN
            iterator = cdiIterator_new(path)
            IF(.NOT. C_ASSOCIATED(iterator%ptr)) CALL finish(routine, "can't open file "//'"'//path//'" for reading')
        END IF
        DO WHILE(me%nextField(iterator, p_patch, level, tileId, variableName, recordsIgnored, lIsFg))
            recordsRead = recordsRead + 1
            ! We have now found the next field that we are interested IN.
            listEntry => me%findTranslatedName(variableName)
            IF(.NOT.ASSOCIATED(listEntry)) CALL finish(routine, "Assertion failed: Processes have different input request lists!")
            domainData => findDomainData(listEntry, p_patch%id, opt_lcreate = .TRUE.)
            CALL domainData%container%readField(level, tileId, p_patch%id, iterator)
            DEALLOCATE(variableName)
        END DO
        IF(my_process_is_mpi_workroot()) THEN
            IF(msg_level > 4) WRITE(0, *) routine//": READ "//TRIM(int2string(recordsRead))//" records from file '"//path//"', &
                                          &ignoring "//TRIM(int2string(recordsIgnored))//" records"
            CALL cdiIterator_delete(iterator)
        END IF
    END SUBROUTINE InputRequestList_readFile

    FUNCTION InputRequestList_getLevels(me, varName, domain, opt_lDebug) RESULT(RESULT)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, VALUE :: domain
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug
        REAL(dp), POINTER :: RESULT(:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_getLevels"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch level data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, domain)
        RESULT => NULL()
        IF(ASSOCIATED(domainData)) RESULT => domainData%container%getLevels()
    END FUNCTION InputRequestList_getLevels

    LOGICAL FUNCTION InputRequestList_fetch2d(me, varName, level, tile, jg, outData, opt_lDebug) RESULT(RESULT)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        REAL(dp), VALUE :: level
        INTEGER, VALUE :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetch2d"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        TYPE(t_tile) :: tileinfo
        LOGICAL :: debugInfo

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        RESULT = ASSOCIATED(domainData)
        IF(RESULT) RESULT = domainData%container%fetch2d(level, tile, outData, opt_lDebug)
        IF(RESULT) THEN
            tileinfo = select_tile(tile)
            CALL initicon_inverse_post_op( &
            &   TRIM(varName//TRIM(get_tile_suffix(tileinfo%GRIB2_tile%tileIndex, tileinfo%GRIB2_att%tileAttribute))), &
            &   optvar_out2D=outData)
        ELSE IF(debugInfo) THEN
            CALL message(routine, "InputContainer_fetch2d() returned an error")
        END IF
    END FUNCTION InputRequestList_fetch2d

    LOGICAL FUNCTION InputRequestList_fetch3d(me, varName, tile, jg, outData, optLevelDimension, opt_lDebug) RESULT(RESULT)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, VALUE :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: optLevelDimension
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetch3d"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        TYPE(t_tile) :: tileinfo
        LOGICAL :: debugInfo

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        RESULT = ASSOCIATED(domainData)
        IF(RESULT) RESULT = domainData%container%fetch3d(tile, outData, optLevelDimension, opt_lDebug)
        IF(RESULT .AND. varName /= 'smi' .AND. varName /= 'SMI') THEN   !SMI IS NOT IN the ICON variable lists, so we need to skip inverse postprocessing for it manually.
            tileinfo = select_tile(tile)
            CALL initicon_inverse_post_op( &
            &   TRIM(varName//TRIM(get_tile_suffix(tileinfo%GRIB2_tile%tileIndex, tileinfo%GRIB2_att%tileAttribute))), &
            &   optvar_out3D=outData)
        ELSE IF(debugInfo) THEN
            CALL message(routine, "InputContainer_fetch3d() returned an error")
        END IF
    END FUNCTION InputRequestList_fetch3d

    LOGICAL FUNCTION InputRequestList_fetchSurface(me, varName, tile, jg, outData, opt_lDebug) RESULT(RESULT)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, VALUE :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchSurface"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        REAL(dp), POINTER :: levels(:)
        TYPE(t_tile) :: tileinfo
        LOGICAL :: debugInfo

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        RESULT = ASSOCIATED(domainData)
        IF(RESULT) THEN
            levels => domainData%container%getLevels()
            SELECT CASE(SIZE(levels, 1))
                CASE(0)
                    RESULT = .FALSE.
                    IF(debugInfo) CALL message(routine, "no levels found")
                CASE(1)
                    RESULT = domainData%container%fetch2d(levels(1), tile, outData, opt_lDebug)
                    IF(debugInfo .AND. .NOT. RESULT) CALL message(routine, "InputContainer_fetch2d() returned an error")
                CASE DEFAULT
                    CALL finish(routine, "trying to read '"//varName//"' as a surface variable, but the file contains several &
                                         &levels of this variable")
            END SELECT
        END IF
        IF(RESULT) THEN
            tileinfo = select_tile(tile)
            CALL initicon_inverse_post_op( &
            &   TRIM(varName//TRIM(get_tile_suffix(tileinfo%GRIB2_tile%tileIndex, tileinfo%GRIB2_att%tileAttribute))), &
            &   optvar_out2D=outData)
        END IF
    END FUNCTION InputRequestList_fetchSurface

    LOGICAL FUNCTION InputRequestList_fetchTiled2d(me, varName, level, jg, outData, opt_lDebug) RESULT(RESULT)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        REAL(dp), VALUE :: level
        INTEGER, VALUE :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchTiled2d"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        INTEGER :: i
        TYPE(t_tile) :: tileinfo
        LOGICAL :: debugInfo

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        RESULT = ASSOCIATED(domainData)
        IF(RESULT) RESULT = domainData%container%fetchTiled2d(level, outData, opt_lDebug)
        IF(RESULT) THEN
            DO i = 1, SIZE(outData, 3)
                tileinfo = select_tile(i)
                CALL initicon_inverse_post_op(TRIM(varName//TRIM(get_tile_suffix(tileinfo%GRIB2_tile%tileIndex, &
                &   tileinfo%GRIB2_att%tileAttribute))), optvar_out2D=outData(:,:,i))
            END DO
        ELSE IF(debugInfo) THEN
            CALL message(routine, "InputContainer_fetchTiled2d() returned an error")
        END IF
    END FUNCTION InputRequestList_fetchTiled2d

    LOGICAL FUNCTION InputRequestList_fetchTiled3d(me, varName, jg, outData, optLevelDimension, opt_lDebug) RESULT(RESULT)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, VALUE :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: optLevelDimension
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchTiled3d"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        INTEGER :: i
        TYPE(t_tile) :: tileinfo
        LOGICAL :: debugInfo

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        RESULT = ASSOCIATED(domainData)
        IF(RESULT) RESULT = domainData%container%fetchTiled3d(outData, optLevelDimension, opt_lDebug)
        IF(RESULT .AND. varName /= 'smi' .AND. varName /= 'SMI') THEN   !SMI IS NOT IN the ICON variable lists, so we need to skip inverse postprocessing for it manually.
            DO i = 1, SIZE(outData, 4)
                tileinfo = select_tile(i)
                CALL initicon_inverse_post_op(TRIM(varName//TRIM(get_tile_suffix(tileinfo%GRIB2_tile%tileIndex, &
                &   tileinfo%GRIB2_att%tileAttribute))), optvar_out3D=outData(:,:,:,i))
            END DO
        ELSE IF(debugInfo) THEN
            CALL message(routine, "InputContainer_fetchTiled3d() returned an error")
        END IF
    END FUNCTION InputRequestList_fetchTiled3d

    LOGICAL FUNCTION InputRequestList_fetchTiledSurface(me, varName, jg, outData, opt_lDebug) RESULT(RESULT)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, VALUE :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchTiledSurface"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        REAL(dp), POINTER :: levels(:)
        INTEGER :: i
        TYPE(t_tile) :: tileinfo
        LOGICAL :: debugInfo

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        RESULT = ASSOCIATED(domainData)
        IF(RESULT) THEN
            levels => domainData%container%getLevels()
            SELECT CASE(SIZE(levels, 1))
                CASE(0)
                    RESULT = .FALSE.
                    IF(debugInfo) CALL message(routine, "no levels found")
                CASE(1)
                    RESULT = domainData%container%fetchTiled2d(levels(1), outData, opt_lDebug)
                    IF(debugInfo .AND. .NOT. RESULT) CALL message(routine, "InputContainer_fetch2d() returned an error")
                CASE DEFAULT
                    CALL finish(routine, "trying to read '"//varName//"' as a surface variable, but the file contains several &
                                         &levels of this variable")
            END SELECT
        END IF
        IF(RESULT) THEN
            DO i = 1, SIZE(outData, 3)
                tileinfo = select_tile(i)
                CALL initicon_inverse_post_op(TRIM(varName//TRIM(get_tile_suffix(tileinfo%GRIB2_tile%tileIndex, &
                &   tileinfo%GRIB2_att%tileAttribute))), optvar_out2D=outData(:,:,i))
            END DO
        END IF
    END FUNCTION InputRequestList_fetchTiledSurface

    SUBROUTINE InputRequestList_fetchRequired2d(me, varName, level, tile, jg, outData)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        REAL(dp), VALUE :: level
        INTEGER, VALUE :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequired2d"

        IF(.NOT. me%fetch2d(varName, level, tile, jg, outData)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequired2d

    SUBROUTINE InputRequestList_fetchRequired3d(me, varName, tile, jg, outData, optLevelDimension)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, VALUE :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: optLevelDimension

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequired3d"

        IF(.NOT. me%fetch3d(varName, tile, jg, outData, optLevelDimension)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequired3d

    SUBROUTINE InputRequestList_fetchRequiredSurface(me, varName, tile, jg, outData)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, VALUE :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequiredSurface"

        IF(.NOT. me%fetchSurface(varName, tile, jg, outData)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequiredSurface

    SUBROUTINE InputRequestList_fetchRequiredTiled2d(me, varName, level, jg, outData)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        REAL(dp), VALUE :: level
        INTEGER, VALUE :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequiredTiled2d"

        IF(.NOT. me%fetchTiled2d(varName, level, jg, outData)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequiredTiled2d

    SUBROUTINE InputRequestList_fetchRequiredTiled3d(me, varName, jg, outData, optLevelDimension)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, VALUE :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: optLevelDimension

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequiredTiled3d"

        IF(.NOT. me%fetchTiled3d(varName, jg, outData, optLevelDimension)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequiredTiled3d

    SUBROUTINE InputRequestList_fetchRequiredTiledSurface(me, varName, jg, outData)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, VALUE :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequiredTiledSurface"

        IF(.NOT. me%fetchTiledSurface(varName, jg, outData)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequiredTiledSurface

    SUBROUTINE InputRequestList_checkRuntypeAndUuids(me, incrementVariables, gridUuids, lIsFg, lHardCheckUuids)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: incrementVariables(:)
        TYPE(t_uuid), INTENT(IN) :: gridUuids(:)    !< gridUuids(n_dom)
        LOGICAL, VALUE :: lIsFg, lHardCheckUuids

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_checkRuntypeAndUuids"
        INTEGER :: i, jg, expectedRuntype
        TYPE(t_ListEntry), POINTER :: curVar
        TYPE(t_DomainData), POINTER :: curDomain
        CHARACTER(:), POINTER :: varnameString
        CHARACTER(LEN = uuid_string_length) :: expectedUuid, foundUuid

        IF(.NOT.my_process_is_mpi_workroot()) CALL finish(routine, "assertion failed")
        IF(SIZE(gridUuids, 1) /= n_dom) CALL finish(routine, "assertion failed")
        DO jg = 1, n_dom
            DO i = 1, me%variableCount
                curVar => me%list(i)
                curDomain => findDomainData(curVar, jg)
                IF(.NOT.ASSOCIATED(curDomain)) CYCLE
                varnameString => toCharacter(curVar%iconVarName)

                ! first check the TYPE of the generating process of the DATA
                IF(.NOT.lIsFg) THEN
                    IF(one_of(varnameString, incrementVariables) > 0) THEN
                        expectedRuntype = 201   ! analysis increment variables
                    ELSE
                        expectedRuntype = 0 ! analysis full variables
                    END IF
                    IF(expectedRuntype /= curDomain%metadata%generatingProcessType) THEN
                        CALL finish(routine, "detected wrong type of generating process on variable '"//varnameString//"', &
                                             &expected "//TRIM(int2string(expectedRuntype))//", &
                                             &found "//TRIM(int2string(curDomain%metadata%generatingProcessType)))
                    END IF
                END IF

                ! second check the UUID of the grid
                IF(.NOT.(gridUuids(jg) == curDomain%metadata%gridUuid)) THEN
                    CALL uuid_unparse(gridUuids(jg), expectedUuid)
                    CALL uuid_unparse(curDomain%metadata%gridUuid, foundUuid)
                    IF(lHardCheckUuids) THEN
                        CALL finish(routine, "detected wrong UUID of grid for variable '"//varnameString//"', &
                                             &expected "//expectedUuid//", &
                                             &found "//foundUuid)
                    ELSE
                        CALL message(routine, "warning: unexpected UUID of grid for variable '"//varnameString//"', &
                                             &expected "//expectedUuid//", &
                                             &found "//foundUuid)
                    END IF
                END IF
                DEALLOCATE(varnameString)
            END DO
        END DO
    END SUBROUTINE InputRequestList_checkRuntypeAndUuids

    SUBROUTINE InputRequestList_printInventory(me)
        CLASS(t_InputRequestList), INTENT(IN) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_printInventory"
        INTEGER :: i, jg, curRow, levelCount, tileCount
        LOGICAL :: lUntiledData
        TYPE(t_ListEntry), POINTER :: curVar
        TYPE(t_DomainData), POINTER :: curDomain
        CHARACTER(:), POINTER :: varnameString, rtimeString, vtimeString
        TYPE(t_table) :: table
        CHARACTER(*), PARAMETER :: domainCol = "jg", &
                                 & variableCol = "variable", &
                                 & tripleCol = "triple", &
                                 & vtimeCol = "validity time", &
                                 & levelTypeCol = "level type", &
                                 & levelCountCol = "level count", &
                                 & tileCountCol = "tileCount", &
                                 & untiledCol = "untiled data", &
                                 & runtypeCol = "runtype", &
                                 & vvmmCol = "vvmm", &
                                 & clasCol = "clas", &
                                 & expidCol = "expid", &
                                 & gridCol = "grid", &
                                 & rgridCol = "rgrid"
        CHARACTER(LEN = 3*3+2) :: parameterString
        TYPE(datetime), POINTER :: rtime, vtime
        TYPE(timedelta), POINTER :: forecastTime
        CHARACTER(len=max_timedelta_str_len) :: forecastTimeString

        CALL initialize_table(table)

        CALL add_table_column(table, domainCol)
        CALL add_table_column(table, variableCol)
        CALL add_table_column(table, tripleCol)
        CALL add_table_column(table, vtimeCol)
        CALL add_table_column(table, vvmmCol)
        CALL add_table_column(table, levelTypeCol)
        CALL add_table_column(table, levelCountCol)
        CALL add_table_column(table, tileCountCol)
        CALL add_table_column(table, untiledCol)
        CALL add_table_column(table, runtypeCol)
        CALL add_table_column(table, clasCol)
        CALL add_table_column(table, expidCol)
        CALL add_table_column(table, gridCol)
        CALL add_table_column(table, rgridCol)

        IF(.NOT.my_process_is_mpi_workroot()) CALL finish(routine, "assertion failed")
        curRow = 1  !we can have zero to n_dom rows for each variable, so we can't USE the loop counter for the rows
        DO jg = 1, n_dom
            DO i = 1, me%variableCount
                curVar => me%list(i)
                curDomain => findDomainData(curVar, jg)
                IF(.NOT.ASSOCIATED(curDomain)) CYCLE
                CALL curDomain%container%getCounts(levelCount, tileCount, lUntiledData)

                !domain, NAME, AND triple columns
                CALL set_table_entry(table, curRow, domainCol, TRIM(int2string(curDomain%domain)))
                varnameString => toCharacter(curVar%iconVarName)
                CALL set_table_entry(table, curRow, variableCol, varnameString)
                DEALLOCATE(varnameString)
                WRITE(parameterString, '(3(I3,:,"."))') curDomain%metadata%param%discipline, curDomain%metadata%param%category, &
                &                                    curDomain%metadata%param%number
                CALL set_table_entry(table, curRow, tripleCol, parameterString)


                !date AND forecast time columns
                rtimeString => toCharacter(curDomain%metadata%rtime)
                vtimeString => toCharacter(curDomain%metadata%vtime)
                CALL set_table_entry(table, curRow, vtimeCol, vtimeString)

                rtime => newDatetime(rtimeString)
                vtime => newDatetime(vtimeString)

                forecastTime => newTimedelta("PT00H")  ! this 'initialization' IS necessary, IN order to correctly deal with timedelta=0.
                forecastTime = vtime - rtime
                CALL timedeltaToString(forecastTime, forecastTimeString)
                CALL set_table_entry(table, curRow, vvmmCol, TRIM(forecastTimeString))

                CALL deallocateDatetime(rtime)
                CALL deallocateDatetime(vtime)
                CALL deallocateTimedelta(forecastTime)
                DEALLOCATE(rtimeString)
                DEALLOCATE(vtimeString)


                !the simpler columns
                CALL set_table_entry(table, curRow, levelTypeCol, TRIM(int2string(curDomain%metadata%levelType)))
                CALL set_table_entry(table, curRow, levelCountCol, TRIM(int2string(levelCount)))
                IF(tileCount /= 0) CALL set_table_entry(table, curRow, tileCountCol, TRIM(int2string(tileCount)))
                IF(lUntiledData) THEN
                    CALL set_table_entry(table, curRow, untiledCol, "yes")
                ELSE
                    CALL set_table_entry(table, curRow, untiledCol, "no")
                END IF
                CALL set_table_entry(table, curRow, clasCol, TRIM(int2string(curDomain%metadata%runClass)))
                CALL set_table_entry(table, curRow, expidCol, TRIM(int2string(curDomain%metadata%experimentId)))
                IF(curDomain%metadata%generatingProcessType /= -1) THEN
                    CALL set_table_entry(table, curRow, runtypeCol, TRIM(int2string(curDomain%metadata%generatingProcessType)))
                END IF
                IF(curDomain%metadata%gridNumber /= -1) THEN
                    CALL set_table_entry(table, curRow, gridCol, TRIM(int2string(curDomain%metadata%gridNumber)))
                END IF
                IF(curDomain%metadata%gridPosition /= -1) THEN
                    CALL set_table_entry(table, curRow, rgridCol, TRIM(int2string(curDomain%metadata%gridPosition)))
                END IF


                !next row
                curRow = curRow + 1
            END DO
        END DO

        CALL print_table(table, opt_delimiter = " | ")
        CALL finalize_table(table)

    END SUBROUTINE InputRequestList_printInventory


    SUBROUTINE InputRequestList_destruct(me)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        INTEGER :: i
        TYPE(t_ListEntry), POINTER :: currentEntry
        TYPE(t_DomainData), POINTER :: domainData, domainDataTemp

        DO i = 1, me%variableCount
            currentEntry => me%list(i)
            DEALLOCATE(currentEntry%iconVarName)
            IF(ASSOCIATED(currentEntry%translatedVarName)) DEALLOCATE(currentEntry%translatedVarName)
            domainData => currentEntry%domainData
            DO
                IF(.NOT.ASSOCIATED(domainData)) EXIT
                IF(ASSOCIATED(domainData%container)) THEN
                    CALL domainData%container%destruct()
                    DEALLOCATE(domainData%container)
                END IF
                IF(ASSOCIATED(domainData%metadata)) CALL MetadataCache_delete(domainData%metadata)
                domainDataTemp => domainData%next
                DEALLOCATE(domainData)
                domainData => domainDataTemp
            END DO
        END DO
        DEALLOCATE(me%list)
    END SUBROUTINE InputRequestList_destruct

    FUNCTION MetadataCache_create() RESULT(RESULT)
        TYPE(t_MetadataCache), POINTER :: RESULT

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":MetadataCache_create"
        INTEGER :: error

        ALLOCATE(RESULT, STAT = error)
        IF(error /= success) CALL finish(routine, "memory allocation error")
        RESULT%vtime => NULL()
        RESULT%rtime => NULL()
    END FUNCTION MetadataCache_create

    LOGICAL FUNCTION MetadataCache_equalTo(me, other) RESULT(RESULT)
        CLASS(t_MetadataCache), INTENT(IN) :: me, other

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":MetadataCache_create"

        INTEGER :: gridNumber, gridPosition, runClass, experimentId, generatingProcessType

        RESULT = .FALSE.

        !compare the time strings
        IF(.NOT.ASSOCIATED(me%rtime).OR..NOT.ASSOCIATED(other%rtime)) THEN
            CALL finish(routine, "internal error, please report this bug")
        END IF
        IF(SIZE(me%rtime) /= SIZE(other%rtime) .OR. ANY(me%rtime /= other%rtime)) THEN
            CALL message(routine, "inconsistent rtime detected")
            RETURN
        END IF

        IF(.NOT.ASSOCIATED(me%vtime).OR..NOT.ASSOCIATED(other%vtime)) THEN
            CALL finish(routine, "internal error, please report this bug")
        END IF
        IF(SIZE(me%vtime) /= SIZE(other%vtime) .OR. ANY(me%vtime /= other%vtime)) THEN
            CALL message(routine, "inconsistent vtime detected")
            RETURN
        END IF

        !compare the parameters
        IF(me%param%discipline /= other%param%discipline) THEN
            CALL message(routine, "inconsistent discipline detected")
            RETURN
        END IF
        IF(me%param%category /= other%param%category) THEN
            CALL message(routine, "inconsistent category detected")
            RETURN
        END IF
        IF(me%param%number /= other%param%number) THEN
            CALL message(routine, "inconsistent number detected")
            RETURN
        END IF

        !compare the other fields
        IF(me%levelType /= other%levelType) THEN
            CALL message(routine, "inconsistent level type detected")
            RETURN
        END IF

        IF(me%gridNumber /= other%gridNumber) THEN
            CALL message(routine, "inconsistent number of grids detected")
            RETURN
        END IF

        IF(me%gridPosition /= other%gridPosition) THEN
            CALL message(routine, "inconsistent grid index detected")
            RETURN
        END IF

        IF(me%runClass /= other%runClass) THEN
            CALL message(routine, "inconsistent run CLASS detected")
            RETURN
        END IF

        IF(me%experimentId /= other%experimentId) THEN
            CALL message(routine, "inconsistent experiment ID detected")
            RETURN
        END IF

        IF(me%generatingProcessType /= other%generatingProcessType) THEN
            CALL message(routine, "inconsistent type of generating process detected")
            RETURN
        END IF

        IF(.NOT.(me%gridUuid == other%gridUuid)) THEN
            CALL message(routine, "inconsistent UUID of grid detected")
            RETURN
        END IF

        RESULT = .TRUE.
    END FUNCTION MetadataCache_equalTo

    SUBROUTINE MetadataCache_delete(me)
        TYPE(t_MetadataCache), POINTER, INTENT(INOUT) :: me

        IF(ASSOCIATED(me%vtime)) DEALLOCATE(me%vtime)
        DEALLOCATE(me)
    END SUBROUTINE MetadataCache_delete

END MODULE mo_input_request_list
