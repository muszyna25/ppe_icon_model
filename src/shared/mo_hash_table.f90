!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A generic hash table
! A container for any kind of objects that are stored under any kind of keys, with O(1) complexity for insertion, lookup, and removal.
!
! Implementation note:
! The loops look more complex than they are:
! Their bodies are executed only once in most cases, sometimes zero times, and more than once in even less cases.
! This is due to the fact that we grow the hash table to be at least as large as we have entries.
! As such, the lists have at most one entry on average.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_hash_table
    USE ISO_C_BINDING, ONLY: C_INT32_T
    USE mo_exception, ONLY: finish
    USE mo_fortran_tools, ONLY: t_Destructible
    USE mo_impl_constants, ONLY: SUCCESS

    IMPLICIT NONE

    PUBLIC :: t_HashTable, hashTable_make
    PRIVATE

    ABSTRACT INTERFACE
        INTEGER(C_INT32_T) FUNCTION f_hashFunction(key)
            IMPORT C_INT32_T, t_Destructible
            CLASS(t_Destructible), POINTER, INTENT(IN) :: key
        END FUNCTION f_hashFunction

        LOGICAL FUNCTION f_equalKeysFunction(keyA, keyB)
            IMPORT t_Destructible
            CLASS(t_Destructible), POINTER, INTENT(IN) :: keyA, keyB
        END FUNCTION f_equalKeysFunction
    END INTERFACE

    TYPE :: t_HashTable
        PRIVATE
        PROCEDURE(f_hashFunction), NOPASS, POINTER :: getHash
        PROCEDURE(f_equalKeysFunction), NOPASS, POINTER :: equalKeys
        TYPE(t_HashEntryPtr), POINTER :: table(:)   ! grows in powers of two
        INTEGER(C_INT32_T) :: entryCount   ! the actual entry count (<= table size && >= non-NULL pointers in table due to the possible chaining of entries)
        INTEGER :: hashBits ! the current count of bits that are used to index the hash table
    CONTAINS
        PROCEDURE :: setEntry => hashTable_setEntry
        PROCEDURE :: removeEntry => hashTable_removeEntry
        PROCEDURE :: getEntry => hashTable_getEntry
        PROCEDURE :: destruct => hashTable_destruct

        PROCEDURE, PRIVATE :: findBin => hashTable_findBin
        PROCEDURE, PRIVATE :: growTable => hashTable_growTable
        PROCEDURE, PRIVATE :: removeFromList => hashTable_removeFromList
    END TYPE

    TYPE :: t_HashEntryPtr
        TYPE(t_HashEntry), POINTER :: ptr
    END TYPE
    TYPE :: t_HashEntry
        CLASS(t_Destructible), POINTER :: key, value
        TYPE(t_HashEntryPtr) :: next
        INTEGER(C_INT32_T) :: hash
    END TYPE

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_hash_table"

CONTAINS

    FUNCTION hashTable_make(hashFunction, compareFunction) RESULT(result)
        PROCEDURE(f_hashFunction) :: hashFunction
        PROCEDURE(f_equalKeysFunction) :: compareFunction
        TYPE(t_HashTable), POINTER :: result

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":hashTable_make"
        INTEGER :: error, i

        ALLOCATE(result, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        result%getHash => hashFunction
        result%equalKeys => compareFunction
        result%entryCount = 0
        result%hashBits = 5

        ALLOCATE(result%table(2**result%hashBits), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        DO i = 1, 2**result%hashBits
            result%table(i)%ptr => NULL()
        END DO
    END FUNCTION hashTable_make

    FUNCTION hashTable_findBin(me, hash) RESULT(result)
        CLASS(t_HashTable), INTENT(IN) :: me
        INTEGER(C_INT32_T), VALUE :: hash
        TYPE(t_HashEntryPtr), POINTER :: result

        INTEGER(C_INT32_T) :: reducedHash, i

        IF(hash < 0) hash = NOT(hash)   !fortran has no unsigned types
        reducedHash = 0
        DO i = 1, (31 + me%hashBits - 1)/me%hashBits
            reducedHash = IEOR(reducedHash, hash)
            hash = ISHFT(hash, -me%hashBits)
        END DO
        reducedHash = IAND(reducedHash, 2**me%hashBits - 1) + 1
        result => me%table(reducedHash)
    END FUNCTION hashTable_findBin

    SUBROUTINE hashTable_growTable(me)
        CLASS(t_HashTable), INTENT(INOUT) :: me

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":hashTable_growTable"
        TYPE(t_HashEntry), POINTER :: curEntry, nextEntry
        TYPE(t_HashEntryPtr), POINTER :: oldTable(:), bin
        INTEGER :: i, error

        oldTable => me%table

        ! Create a new empty table.
        me%hashBits = me%hashBits + 1
        ALLOCATE(me%table(2**me%hashBits), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        DO i = 1, 2**me%hashBits
            me%table(i)%ptr => NULL()
        END DO

        ! Move over the contents of the old table.
        DO i = 1, SIZE(oldTable, 1)
            curEntry => oldTable(i)%ptr
            DO WHILE(ASSOCIATED(curEntry))
                nextEntry => curEntry%next%ptr

                ! Insert curEntry into the new table.
                bin => me%findBin(curEntry%hash)
                curEntry%next%ptr => bin%ptr
                bin%ptr => curEntry

                curEntry => nextEntry
            END DO
        END DO

        ! Cleanup the old table.
        DEALLOCATE(oldTable)
    END SUBROUTINE hashTable_growTable

    SUBROUTINE hashTable_removeFromList(me, list, key, hash)
        CLASS(t_HashTable), INTENT(INOUT) :: me
        TYPE(t_HashEntryPtr), POINTER, INTENT(INOUT) :: list
        CLASS(t_Destructible), POINTER, INTENT(IN) :: key
        INTEGER(C_INT32_T), VALUE :: hash

        TYPE(t_HashEntry), POINTER :: curEntry
        TYPE(t_HashEntryPtr), POINTER :: iterator

        iterator => list
        DO WHILE(ASSOCIATED(iterator%ptr))
            curEntry => iterator%ptr
            IF(curEntry%hash == hash) THEN
                IF(me%equalKeys(curEntry%key, key)) THEN
                    ! remove from list
                    iterator%ptr => curEntry%next%ptr
                    me%entryCount = me%entryCount - 1

                    ! destroy the entry
                    CALL curEntry%key%destruct()
                    DEALLOCATE(curEntry%key)
                    CALL curEntry%VALUE%destruct()
                    DEALLOCATE(curEntry%value)
                    DEALLOCATE(curEntry)
                    CYCLE
                END IF
            END IF
            iterator => curEntry%next    ! point to next entry
        END DO
    END SUBROUTINE hashTable_removeFromList

    ! The hash table takes possession of both the key and the value and will DEALLOCATE() them eventually.
    SUBROUTINE hashTable_setEntry(me, key, value)
        CLASS(t_HashTable), INTENT(INOUT) :: me
        CLASS(t_Destructible), POINTER, INTENT(IN) :: key, value

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":hashTable_setEntry"
        TYPE(t_HashEntry), POINTER :: newEntry
        TYPE(t_HashEntryPtr), POINTER :: bin
        INTEGER :: error

        ! Prepare the new entry.
        ALLOCATE(newEntry, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        newEntry%key => key
        newEntry%value => value
        newEntry%hash = me%getHash(key)
        newEntry%next%ptr => NULL()

        ! If there is a prexeisting entry for this key, remove it.
        bin => me%findBin(newEntry%hash)
        CALL me%removeFromList(bin, key, newEntry%hash)

        ! Add the new entry to the list.
        newEntry%next%ptr => bin%ptr
        bin%ptr => newEntry
        me%entryCount = me%entryCount + 1

        ! Check whether we need to grow the table.
        IF(me%entryCount == SIZE(me%table, 1)) CALL me%growTable()
    END SUBROUTINE hashTable_setEntry

    SUBROUTINE hashTable_removeEntry(me, key)
        CLASS(t_HashTable), INTENT(INOUT) :: me
        CLASS(t_Destructible), POINTER, INTENT(IN) :: key

        INTEGER(C_INT32_T) :: hash
        TYPE(t_HashEntryPtr), POINTER :: bin

        hash = me%getHash(key)
        bin => me%findBin(hash)
        CALL me%removeFromList(bin, key, hash)
    END SUBROUTINE hashTable_removeEntry

    FUNCTION hashTable_getEntry(me, key) RESULT(result)
        CLASS(t_HashTable), INTENT(IN) :: me
        CLASS(t_Destructible), POINTER, INTENT(IN) :: key
        CLASS(t_Destructible), POINTER :: result

        INTEGER(C_INT32_T) :: hash
        TYPE(t_HashEntryPtr), POINTER :: bin
        TYPE(t_HashEntry), POINTER :: curEntry

        result => NULL()
        hash = me%getHash(key)
        bin => me%findBin(hash)
        curEntry => bin%ptr
        DO WHILE(ASSOCIATED(curEntry))
            IF(curEntry%hash == hash) THEN
                IF(me%equalKeys(curEntry%key, key)) THEN
                    result => curEntry%value
                    RETURN
                END IF
            END IF
            curEntry => curEntry%next%ptr
        END DO
    END FUNCTION hashTable_getEntry

    SUBROUTINE hashTable_destruct(me)
        CLASS(t_HashTable), INTENT(INOUT) :: me

        INTEGER :: i
        TYPE(t_HashEntry), POINTER :: curEntry, nextEntry

        DO i = 1, SIZE(me%table)
            curEntry => me%table(i)%ptr
            DO WHILE(ASSOCIATED(curEntry))
                nextEntry => curEntry%next%ptr
                CALL curEntry%key%destruct()
                DEALLOCATE(curEntry%key)
                CALL curEntry%VALUE%destruct()
                DEALLOCATE(curEntry%value)
                DEALLOCATE(curEntry)
                curEntry => nextEntry
            END DO
        END DO
        DEALLOCATE(me%table)
    END SUBROUTINE hashTable_destruct

END MODULE mo_hash_table
