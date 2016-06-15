!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_input_container
    USE ISO_C_BINDING, ONLY: C_INT32_T, C_INT64_T, C_DOUBLE, C_FLOAT, C_ASSOCIATED
    USE mo_cdi, ONLY: t_CdiIterator, gridInqSize, cdiIterator_inqGridId, cdiIterator_readField, cdiIterator_readFieldF, &
                    & cdiIterator_inqDatatype, DATATYPE_PACK23, DATATYPE_PACK32, DATATYPE_FLT64, DATATYPE_INT32
    USE mo_communication, ONLY: t_ScatterPattern
    USE mo_exception, ONLY: message, finish
    USE mo_fortran_tools, ONLY: assign_if_present, t_Destructible
    USE mo_hash_table, ONLY: t_HashTable, hashTable_make
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_kind, ONLY: wp, dp
    USE mo_mpi, ONLY: p_bcast, p_comm_work, my_process_is_stdio
    USE mo_parallel_config, ONLY: blk_no, nproma
    USE mo_scatter_pattern_base, ONLY: lookupScatterPattern
    USE mo_util_cdi, ONLY: trivial_tileId
    USE mo_util_string, ONLY: int2string, real2string

    IMPLICIT NONE

PUBLIC :: t_InputContainer, InputContainer_make

    TYPE :: t_ValueList
        PRIVATE
        REAL(dp), POINTER :: values(:)
        INTEGER :: valueCount
    CONTAINS
        PROCEDURE :: init => ValueList_init
        PROCEDURE :: addValue => ValueList_addValue
        PROCEDURE :: haveValue => ValueList_haveValue
        PROCEDURE :: ensureSpace => ValueList_ensureSpace
        PROCEDURE :: destruct => ValueList_destruct
    END TYPE

    ! A t_InputContainer is used for file driven I/O to cache input DATA IN a distributed fashion
    ! until we have enough information to ALLOCATE the arrays to hold this DATA.
    ! This is especially needed for the fields READ by request_dwdfg_atm_ii() AND fetch_dwdfg_atm_ii()
    ! for which we have no prior knowledge about the count of levels PRESENT IN the input file.
    !
    ! t_InputContainers are always created via a t_InputRequestList which distributes the DATA to the different t_InputContainer objects that it created.
    !
    ! TODO (optional): make it possible to use a different communicator than p_comm_work.
    TYPE :: t_InputContainer
        PRIVATE
        !the NAME of the variable is NOT a part of the t_InputContainer since it's ONLY needed by the t_InputRequestList
        TYPE(t_HashTable), POINTER :: fields    !this is a collection of all the different 2D levels that have been READ (tile, level)
        INTEGER :: fieldCount
        TYPE(t_ValueList) :: tiles, levels

    CONTAINS
        PROCEDURE :: getLevels => InputContainer_getLevels
        PROCEDURE :: getCounts => InputContainer_getCounts
        PROCEDURE :: destruct => InputContainer_destruct

        !the fetch operations are atomical, the output array IS NOT touched upon failure
        PROCEDURE :: fetch2D => InputContainer_fetch2D  !returns .TRUE. on success
        PROCEDURE :: fetch3D => InputContainer_fetch3D  !returns .TRUE. on success
        PROCEDURE :: fetchTiled2D => InputContainer_fetchTiled2D  !returns .TRUE. on success
        PROCEDURE :: fetchTiled3D => InputContainer_fetchTiled3D  !returns .TRUE. on success

        PROCEDURE :: readField => InputContainer_readField
        PROCEDURE, PRIVATE :: dataAvailable => InputContainer_dataAvailable
    END TYPE

PRIVATE

    TYPE, EXTENDS(t_Destructible) :: t_LevelKey
        PRIVATE
        REAL(dp) :: levelValue
        INTEGER :: tileId   !< the ICON-internal tile id

    CONTAINS
        PROCEDURE :: destruct => LevelKey_destruct  !< override
    END TYPE

    TYPE, EXTENDS(t_Destructible) :: t_LevelPointer
        PRIVATE
        REAL(wp), POINTER :: ptr(:,:)   !dimensions: (nproma, nblks_*)

    CONTAINS
        PROCEDURE :: destruct => LevelPointer_destruct  !< override
    END TYPE

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_input_container"
    REAL(dp), PARAMETER :: k_levelComparisonInterval = 1e-5

CONTAINS

    INTEGER(C_INT32_T) FUNCTION InputContainer_hashKey(key) RESULT(result)
        CLASS(t_Destructible), POINTER, INTENT(IN) :: key

        ! Just some large primes to produce some good pseudorandom bits in the hashes.
        INTEGER(C_INT64_T), PARAMETER :: prime1 = 729326603, prime2 = 941095657
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":InputContainer_hashKey"

        INTEGER(C_INT64_T) :: temp
        INTEGER(KIND = C_INT64_T), PARAMETER :: mask = 2_C_INT64_T**31 - 1_C_INT64_T    ! a bitmask for the 31 low order bits

        SELECT TYPE(key)
            TYPE IS(t_LevelKey)
                temp = prime1*INT(key%levelValue, C_INT64_T)
                temp = temp + prime2*INT(key%tileId, C_INT64_T)
                RESULT = INT(IAND(temp, mask), C_INT32_T)
            CLASS DEFAULT
                CALL finish(routine, "illegal argument type")
        END SELECT
    END FUNCTION InputContainer_hashKey

    LOGICAL FUNCTION InputContainer_equalKeysFunction(keyA, keyB) RESULT(result)
        CLASS(t_Destructible), POINTER, INTENT(IN) :: keyA, keyB

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":InputContainer_equalKeysFunction"

        SELECT TYPE(keyA)
            TYPE IS(t_LevelKey)
                SELECT TYPE(keyB)
                    TYPE IS(t_LevelKey)
                        result = keyA%levelValue > keyB%levelValue - k_levelComparisonInterval .and. &
                               & keyA%levelValue < keyB%levelValue + k_levelComparisonInterval
                        result = result .and. keyA%tileId == keyB%tileId
                    CLASS DEFAULT
                        CALL finish(routine, "illegal argument type")
                END SELECT
            CLASS DEFAULT
                CALL finish(routine, "illegal argument type")
        END SELECT
    END FUNCTION InputContainer_equalKeysFunction

    SUBROUTINE LevelKey_destruct(me)
        CLASS(t_LevelKey), INTENT(INOUT) :: me
    END SUBROUTINE LevelKey_destruct

    SUBROUTINE LevelPointer_destruct(me)
        CLASS(t_LevelPointer), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":LevelPointer_destruct"

        DEALLOCATE(me%ptr)
    END SUBROUTINE LevelPointer_destruct

    FUNCTION InputContainer_make() RESULT(RESULT)
        CLASS(t_InputContainer), POINTER :: RESULT

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_make"
        INTEGER :: error

        ALLOCATE(RESULT, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        RESULT%fields => hashTable_make(InputContainer_hashKey, InputContainer_equalKeysFunction)
        RESULT%fieldCount = 0
        CALL RESULT%tiles%init()
        CALL RESULT%levels%init()
    END FUNCTION InputContainer_make

    FUNCTION InputContainer_getLevels(me) RESULT(RESULT)
        CLASS(t_InputContainer), INTENT(IN) :: me
        REAL(dp), POINTER :: RESULT(:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_getLevels"
        INTEGER :: i, error

        ALLOCATE(RESULT(me%levels%valueCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        DO i = 1, me%levels%valueCount
            RESULT(i) = me%levels%values(i)
        END DO
    END FUNCTION InputContainer_getLevels

    SUBROUTINE InputContainer_getCounts(me, levelCount, tileCount, lHaveUntiledData)
        CLASS(t_InputContainer), INTENT(IN) :: me
        INTEGER, INTENT(OUT) :: levelCount, tileCount
        LOGICAL, INTENT(OUT) :: lHaveUntiledData

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_getCounts"

        levelCount = me%levels%valueCount
        tileCount = me%tiles%valueCount
        lHaveUntiledData = me%tiles%haveValue(REAL(trivial_tileId, dp))
        IF(lHaveUntiledData) tileCount = tileCount - 1
    END SUBROUTINE InputContainer_getCounts

    LOGICAL FUNCTION InputContainer_dataAvailable(me, opt_level, opt_tile, opt_jg) RESULT(RESULT)
        CLASS(t_InputContainer), INTENT(IN) :: me
        REAL(dp), OPTIONAL, INTENT(IN) :: opt_level
        INTEGER, OPTIONAL, INTENT(IN) :: opt_tile, opt_jg

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_dataAvailable"
        CLASS(t_Destructible), POINTER :: key, polymorphicKey, VALUE
        INTEGER :: levelCount, tileCount, levelIndex, tileIndex, error

        RESULT = .FALSE.
        IF(me%levels%valueCount <= 0 .OR. me%tiles%valueCount <= 0) RETURN

        levelCount = me%levels%valueCount
        IF(PRESENT(opt_level)) levelCount = 1
        tileCount = me%tiles%valueCount
        IF(PRESENT(opt_tile)) tileCount = 1

        ALLOCATE(t_LevelKey :: key, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        polymorphicKey => key   !Fortran needs this X|
        SELECT TYPE(key)
            TYPE IS(t_LevelKey)

                RESULT = .TRUE.
                outerLoop: DO levelIndex = 1, levelCount
                    key%levelValue = me%levels%values(levelIndex)
                    IF(PRESENT(opt_level)) key%levelValue = opt_level
                    DO tileIndex = 1, tileCount
                        key%tileId = INT(me%tiles%values(tileIndex))
                        IF(PRESENT(opt_tile)) key%tileId = opt_tile

                        VALUE => me%fields%getEntry(polymorphicKey)
                        IF(.NOT.ASSOCIATED(VALUE)) THEN
                            RESULT = .FALSE.
                            EXIT outerLoop
                        END IF
                    END DO
                END DO outerLoop

            CLASS DEFAULT
                CALL finish(routine, "assertion failed")
        END SELECT

        CALL key%destruct()
        DEALLOCATE(key)
    END FUNCTION InputContainer_dataAvailable

    LOGICAL FUNCTION InputContainer_fetch2D(me, level, tile, outData, opt_lDebug) RESULT(RESULT)
        CLASS(t_InputContainer), INTENT(IN) :: me
        REAL(dp), VALUE :: level
        INTEGER, VALUE :: tile
        REAL(wp), INTENT(INOUT) :: outData(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_fetch2D"
        INTEGER :: error
        CLASS(t_Destructible), POINTER :: key, VALUE
        LOGICAL :: debugInfo

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug
        IF(debugInfo) CALL message(routine, "level = "//TRIM(real2string(level))//", tile = "//TRIM(int2string(tile)))
        RESULT = me%dataAvailable(opt_level = level, opt_tile = tile)
        IF(.NOT.RESULT) RETURN

        ALLOCATE(t_LevelKey :: key, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        SELECT TYPE(key)
            TYPE IS(t_LevelKey)
                key%levelValue = level
                key%tileId = tile
            CLASS DEFAULT
                CALL finish(routine, "assertion failed")
        END SELECT
        VALUE => me%fields%getEntry(key)

        IF(.NOT.ASSOCIATED(VALUE)) CALL finish(routine, "internal error: dataAvailable() returned TRUE, but getEntry() failed")
        SELECT TYPE(VALUE)
            TYPE IS(t_LevelPointer)
                IF(SIZE(outData, 1) /= SIZE(VALUE%ptr, 1) .OR.  SIZE(outData, 2) /= SIZE(VALUE%ptr, 2)) THEN
                   CALL finish(routine, "dimensions of output array do not match the dimensions of the data read from file")
                END IF
                outData(:,:) = VALUE%ptr(:,:)
            CLASS DEFAULT
                CALL finish(routine, "assertion failed")
        END SELECT

        CALL key%destruct()
        DEALLOCATE(key)
    END FUNCTION InputContainer_fetch2D

    LOGICAL FUNCTION InputContainer_fetch3D(me, tile, outData, optLevelDimension, opt_lDebug) RESULT(RESULT)
        CLASS(t_InputContainer), INTENT(IN) :: me
        INTEGER, VALUE :: tile
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        INTEGER, INTENT(IN), OPTIONAL :: optLevelDimension
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_fetch3D"
        INTEGER :: levelDimension, i
        LOGICAL :: debugInfo

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        levelDimension = 2
        CALL assign_if_present(levelDimension, optLevelDimension)
        IF(SIZE(outData, levelDimension) /= me%levels%valueCount) THEN
            RESULT = .FALSE.
            IF(debugInfo) CALL message(routine, "wrong level count")
            RETURN
        END IF

        RESULT = me%dataAvailable(opt_tile = tile)
        IF(.NOT.RESULT) RETURN

        DO i = 1, me%levels%valueCount
            SELECT CASE(levelDimension)
                CASE(2)
                    RESULT = me%fetch2D(me%levels%values(i), tile, outData(:, i, :))
                CASE(3)
                    RESULT = me%fetch2D(me%levels%values(i), tile, outData(:, :, i))
                CASE DEFAULT
                    CALL finish(routine, "illegal value in optleveldimension")
            END SELECT
            IF(.NOT. RESULT) CALL finish(routine, "internal error: dataAvailable() returned TRUE, but a fetch2D() CALL failed")
        END DO
    END FUNCTION InputContainer_fetch3D

    LOGICAL FUNCTION InputContainer_fetchTiled2D(me, level, outData, opt_lDebug) RESULT(RESULT)
        CLASS(t_InputContainer), INTENT(IN) :: me
        REAL(dp), VALUE :: level
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_fetchTiled2D"
        INTEGER :: tileId, expectedTileCount
        LOGICAL :: debugInfo, useUntiledData

        useUntiledData = .FALSE.
        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug
        RESULT = me%dataAvailable(opt_level = level)
        IF(.NOT.RESULT) RETURN

        !check whether we have the right count of tiles, ignoring untiled fields of the same variable that may be PRESENT as well
        IF(me%tiles%haveValue(REAL(trivial_tileId, dp))) THEN
            RESULT = SIZE(outData, 3) == me%tiles%valueCount - 1
            IF(.NOT.RESULT.AND.SIZE(outData, 3) == 1) THEN
                CALL message(routine, "Warning: No tiled input data is present. Reading single tile from untiled data instead.")
                useUntiledData = .TRUE.
                RESULT = .TRUE.
            END IF
            IF(debugInfo .AND. .NOT. RESULT) THEN
                CALL message(routine, "Wrong tile count. &
                                      &Tiles from file = "//TRIM(int2string(me%tiles%valueCount - 1))//" &
                                      &(+ one untiled field), &
                                      &expected tiles = "//TRIM(int2string(SIZE(outData, 3))))
            END IF
        ELSE
            RESULT = SIZE(outData, 3) == me%tiles%valueCount
            IF(debugInfo .AND. .NOT. RESULT) THEN
                CALL message(routine, "Wrong tile count. &
                                      &Tiles from file = "//TRIM(int2string(me%tiles%valueCount))//" &
                                      &expected tiles = "//TRIM(int2string(SIZE(outData, 3))))
            END IF
        END IF
        IF(.NOT. RESULT) RETURN

        DO tileId = 1, SIZE(outData, 3)
            IF(useUntiledData) THEN
                RESULT = me%fetch2D(level, trivial_tileId, outData(:, :, tileId))
            ELSE
                RESULT = me%fetch2D(level, tileId, outData(:, :, tileId))
            END IF
            IF(.NOT. RESULT) CALL finish(routine, "internal error: dataAvailable() returned TRUE, but a fetch2D() CALL failed")
        END DO
    END FUNCTION InputContainer_fetchTiled2D

    LOGICAL FUNCTION InputContainer_fetchTiled3D(me, outData, optLevelDimension, opt_lDebug) RESULT(RESULT)
        CLASS(t_InputContainer), INTENT(IN) :: me
        REAL(wp), INTENT(INOUT) :: outData(:,:,:,:)
        INTEGER, INTENT(IN), OPTIONAL :: optLevelDimension
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_fetchTiled3D"
        INTEGER :: tileId
        LOGICAL :: debugInfo, useUntiledData

        useUntiledData = .FALSE.
        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        !XXX: This IS more strict than strictly necessary:
        !     It will check that untiled DATA that may also be PRESENT has the same amount of levels as the tiled DATA.
        !     I guess, this additional strictness shouldn't be a problem, but I may be wrong.
        RESULT = me%dataAvailable()
        IF(.NOT.RESULT) RETURN

        !check whether we have the right count of tiles, ignoring untiled fields of the same variable that may be PRESENT as well
        IF(me%tiles%haveValue(REAL(trivial_tileId, dp))) THEN
            RESULT = SIZE(outData, 4) == me%tiles%valueCount - 1
            IF(.NOT.RESULT.AND.SIZE(outData, 4) == 1) THEN
                CALL message(routine, "Warning: No tiled input data is present. Reading single tile from untiled data instead.")
                useUntiledData = .TRUE.
                RESULT = .TRUE.
            END IF
            IF(debugInfo .AND. .NOT. RESULT) THEN
                CALL message(routine, "Wrong tile count. &
                                      &Tiles from file = "//TRIM(int2string(me%tiles%valueCount - 1))//" &
                                      &(+ one untiled field), &
                                      &expected tiles = "//TRIM(int2string(SIZE(outData, 4))))
            END IF
        ELSE
            RESULT = SIZE(outData, 4) == me%tiles%valueCount
            IF(debugInfo .AND. .NOT. RESULT) THEN
                CALL message(routine, "Wrong tile count. &
                                      &Tiles from file = "//TRIM(int2string(me%tiles%valueCount))//" &
                                      &expected tiles = "//TRIM(int2string(SIZE(outData, 4))))
            END IF
        END IF
        IF(.NOT. RESULT) RETURN

        DO tileId = 1, SIZE(outData, 4)
            IF(useUntiledData) THEN
                RESULT = me%fetch3D(trivial_tileId, outData(:, :, :, tileId), optLevelDimension)
            ELSE
                RESULT = me%fetch3D(tileId, outData(:, :, :, tileId), optLevelDimension)
            END IF
            IF(.NOT. RESULT) CALL finish(routine, "internal error: dataAvailable() returned TRUE, but a fetch3D() CALL failed")
        END DO
    END FUNCTION InputContainer_fetchTiled3D

    SUBROUTINE InputContainer_destruct(me)
        CLASS(t_InputContainer) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_destruct"

        CALL me%fields%destruct()
        DEALLOCATE(me%fields)
        CALL me%tiles%destruct()
        CALL me%levels%destruct()
    END SUBROUTINE InputContainer_destruct

    SUBROUTINE InputContainer_readField(me, variableName, level, tile, jg, iterator)
        CLASS(t_InputContainer), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: variableName
        REAL(dp), VALUE :: level
        INTEGER, VALUE :: tile, jg
        TYPE(t_CdiIterator), VALUE :: iterator

        CHARACTER(*), PARAMETER :: routine = modname//":InputContainer_readField"
        CLASS(t_Destructible), POINTER :: key, VALUE
        REAL(C_DOUBLE), POINTER :: bufferD(:)
        REAL(C_FLOAT), POINTER :: bufferS(:)
        INTEGER :: gridSize, datatype, packedMessage(2), error    !packedMessage(1) == gridSize, packedMessage(2) == datatype
        CLASS(t_ScatterPattern), POINTER :: distribution

        !sanity check: fail IF this field has already been READ
        ALLOCATE(t_LevelKey :: key, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        SELECT TYPE(key)
            TYPE IS(t_LevelKey)
                key%levelValue = level
                key%tileId = tile
            CLASS DEFAULT
                CALL finish(routine, "assertion failed")
        END SELECT
        IF(ASSOCIATED(me%fields%getEntry(key))) THEN
            CALL finish(routine, "double definition of variable '"//variableName//"' in an input file")
        END IF

        !Inquire buffer SIZE information AND broadcast it.
        IF(C_ASSOCIATED(iterator%ptr)) THEN
            gridSize = gridInqSize(cdiIterator_inqGridId(iterator))
            datatype = cdiIterator_inqDatatype(iterator)
            packedMessage(1) = gridSize
            packedMessage(2) = datatype
        END IF
        CALL p_bcast(packedMessage, 0, p_comm_work)
        gridSize = packedMessage(1)
        datatype = packedMessage(2)

        !Get the corresponding ScatterPattern AND initialize our hash table entry.
        distribution => lookupScatterPattern(jg, gridSize)
        IF(.NOT. ASSOCIATED(distribution)) CALL finish(routine, "could not find scatter pattern to distribute input field")
        ALLOCATE(t_LevelPointer :: VALUE, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        SELECT TYPE(VALUE)
            TYPE IS(t_LevelPointer)

                !ALLOCATE the local buffer
                ALLOCATE(VALUE%ptr(nproma, blk_no(distribution%localSize())), STAT = error)
                IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
                VALUE%ptr(:,:) = 0.0_wp     !Avoid nondeterministic values IN the unused points for checksumming.
                SELECT CASE(datatype)
                    CASE(DATATYPE_PACK23:DATATYPE_PACK32, DATATYPE_FLT64, DATATYPE_INT32)
                        !ALLOCATE the global buffer
                        IF(C_ASSOCIATED(iterator%ptr)) THEN
                            ALLOCATE(bufferD(gridSize), STAT = error)
                        ELSE
                            ALLOCATE(bufferD(1), STAT = error)
                        END IF
                        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")

                        !READ the DATA
                        IF(C_ASSOCIATED(iterator%ptr)) CALL cdiIterator_readField(iterator, bufferD)
                        CALL distribution%distribute(bufferD, VALUE%ptr(:, :), .FALSE.)

                        DEALLOCATE(bufferD)

                    CASE DEFAULT
                        !ALLOCATE the global buffer
                        IF(C_ASSOCIATED(iterator%ptr)) THEN
                            ALLOCATE(bufferS(gridSize), STAT = error)
                        ELSE
                            ALLOCATE(bufferS(1), STAT = error)
                        END IF
                        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")

                        !READ the DATA
                        IF(C_ASSOCIATED(iterator%ptr)) CALL cdiIterator_readFieldF(iterator, bufferS)
                        CALL distribution%distribute(bufferS, VALUE%ptr(:, :), .FALSE.)

                        DEALLOCATE(bufferS)

                END SELECT

            CLASS DEFAULT
                CALL finish(routine, "assertion failed")
        END SELECT

        !store the DATA IN our hash table
        CALL me%fields%setEntry(key, VALUE) !this will DEALLOCATE both the VALUE AND the key eventually
        me%fieldCount = me%fieldCount + 1
        CALL me%tiles%addValue(REAL(tile, dp))
        CALL me%levels%addValue(level)
    END SUBROUTINE InputContainer_readField

    SUBROUTINE ValueList_init(me)
        CLASS(t_ValueList), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":ValueList_init"
        INTEGER :: i, error

        ALLOCATE(me%values(8), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        me%valueCount = 0
        DO i = 1, SIZE(me%values, 1)
            me%values(i) = 0.0
        END DO
    END SUBROUTINE ValueList_init

    SUBROUTINE ValueList_addValue(me, VALUE)
        CLASS(t_ValueList), INTENT(INOUT) :: me
        REAL(dp), VALUE :: VALUE

        CHARACTER(*), PARAMETER :: routine = modname//":ValueList_addValue"
        INTEGER :: i, j

        !check whether the VALUE is already IN the list
        DO i = 1, me%valueCount
            IF(me%values(i) == VALUE) RETURN    !nothing to DO, we already have this value
            IF(me%values(i) > VALUE) EXIT
        END DO

        !insert the VALUE at the point that i is pointing to
        CALL me%ensureSpace()
        DO j = me%valueCount, i, -1
            me%values(j+1) = me%values(j)
        END DO
        me%values(i) = VALUE
        me%valueCount = me%valueCount + 1
    END SUBROUTINE ValueList_addValue

    FUNCTION ValueList_haveValue(me, VALUE) RESULT(RESULT)
        CLASS(t_ValueList), INTENT(IN) :: me
        REAL(dp), VALUE :: VALUE
        LOGICAL :: RESULT

        CHARACTER(*), PARAMETER :: routine = modname//":ValueList_haveValue"
        INTEGER :: i

        RESULT = .TRUE.
        DO i = 1, me%valueCount
            IF(me%values(i) == VALUE) RETURN    !found it, RETURN success
        END DO
        RESULT = .FALSE.
    END FUNCTION ValueList_haveValue

    SUBROUTINE ValueList_ensureSpace(me)
        CLASS(t_ValueList), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":ValueList_ensureSpace"
        INTEGER :: i, error
        REAL(dp), POINTER :: temp(:)

        IF(me%valueCount == SIZE(me%values, 1)) THEN
            temp => me%values
            ALLOCATE(me%values(2*SIZE(me%values, 1)), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
            DO i = 1, me%valueCount
                me%values(i) = temp(i)
            END DO
            DO i = me%valueCount + 1, SIZE(me%values, 1)
                me%values(i) = 0.0
            END DO
            DEALLOCATE(temp)
        END IF
    END SUBROUTINE ValueList_ensureSpace

    SUBROUTINE ValueList_destruct(me)
        CLASS(t_ValueList) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":ValueList_destruct"

        DEALLOCATE(me%values)
    END SUBROUTINE ValueList_destruct

END MODULE mo_input_container
