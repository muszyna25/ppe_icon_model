!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A generic hash table
! A container for any kind of objects that are stored under any kind
! of keys, with O(1) complexity for insertion, lookup, and removal.
!
! Implementation note:
! The loops look more complex than they are:
! Their bodies are executed only once in most cases, sometimes zero
! times, and more than once in even less cases.  This is due to the
! fact that we grow the hash table to be at least as large as we have
! entries.
! As such, the lists have at most one entry on average.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE mo_hash_table
  USE ISO_C_BINDING,     ONLY: C_INT32_T, C_INT64_T
  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: finish
  USE mo_fortran_tools,  ONLY: t_Destructible
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_util_string,    ONLY: tolower
  
  IMPLICIT NONE
  
  PUBLIC :: t_HashTable_Base
  PUBLIC :: t_HashTable, hashTable_make, t_HashIterator, t_const_hashIterator
  PUBLIC :: stringVal_len
  PUBLIC :: t_scalarVal, t_stringVal, t_realVal, t_intVal, t_logVal
  PUBLIC :: storage_hashKey_DJB_cs, storage_hashKey_DJB_ci
  PUBLIC :: storage_equalKeysFunction_cs, storage_equalKeysFunction_ci
  PRIVATE
  
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_hash_table'
  
  ! --------------------------------------------------------------------------------
  ! auxiliary data types for basic types in hash table entries
  
  INTEGER, PARAMETER :: stringVal_len= 120  !< Maximum allowed length for string keys or values
  
  ! constants for hashing
  INTEGER(C_INT32_T), PARAMETER :: i32_t = 0
  INTEGER(C_INT64_T), PARAMETER :: p32 = INT(HUGE(i32_t), C_INT64_T) + 1 
  
  TYPE, EXTENDS(t_Destructible) :: t_scalarVal
  PRIVATE
  CONTAINS
    PROCEDURE :: destruct => scalarVal_destruct
  END TYPE t_scalarVal

  TYPE, EXTENDS(t_scalarVal) :: t_stringVal
    CHARACTER(LEN=stringVal_len) :: stringVal
  END TYPE t_stringVal

  TYPE, EXTENDS(t_scalarVal) :: t_realVal
    REAL(wp) :: realVal
  END TYPE t_realVal

  TYPE, EXTENDS(t_scalarVal) :: t_intVal
    INTEGER :: intVal
  END TYPE t_intVal

  TYPE, EXTENDS(t_scalarVal) :: t_logVal
    LOGICAL :: logVal
  END TYPE t_logVal


  ! --------------------------------------------------------------------------------

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

  ! Provide a base type to t_hashTable that throws a meaningful error if
  ! t_HashTable is used uninitialized
  TYPE :: t_HashTable_Base
    PRIVATE
    PROCEDURE(f_hashFunction), NOPASS, POINTER :: getHash
    PROCEDURE(f_equalKeysFunction), NOPASS, POINTER :: equalKeys
    TYPE(t_HashEntryPtr), POINTER :: table(:)   ! grows in powers of two

    ! the actual entry count (<= table size && >= non-NULL pointers in
    ! table due to the possible chaining of entries):
    INTEGER(C_INT32_T) :: entryCount

    INTEGER :: hashBits ! the current count of bits that are used to index the hash table
  CONTAINS
    PROCEDURE :: setEntry    => hashTable_setEntry_Base
    PROCEDURE :: removeEntry => hashTable_removeEntry_Base
    PROCEDURE :: getEntry    => hashTable_getEntry_Base
    PROCEDURE :: destruct    => hashTable_destruct_Base
    PROCEDURE :: nentries    => hashTable_nentries_Base
    
    PROCEDURE, PRIVATE :: findBin        => hashTable_findBin
    PROCEDURE, PRIVATE :: growTable      => hashTable_growTable
    PROCEDURE, PRIVATE :: removeFromList => hashTable_removeFromList
  END TYPE t_HashTable_Base

  TYPE, EXTENDS(t_HashTable_Base) :: t_HashTable
  CONTAINS
    PROCEDURE :: setEntry    => hashTable_setEntry
    PROCEDURE :: removeEntry => hashTable_removeEntry
    PROCEDURE :: getEntry    => hashTable_getEntry
    PROCEDURE :: destruct    => hashTable_destruct
    PROCEDURE :: nentries    => hashTable_nentries
  END TYPE t_HashTable


  ! Provides sequential access to all entries of a hash table
  !
  ! <<< WARNING >>>:
  !   iterators are invalidated by setEntry(), removeEntry(), and
  !   destruct() calls on the corresponding hash table
  !
  TYPE, ABSTRACT :: t_HashIteratorBase
    PRIVATE
    CLASS(t_HashTable_Base), POINTER :: table
    INTEGER :: curBin
    TYPE(t_HashEntry), POINTER :: curEntry
  CONTAINS
    PROCEDURE :: init => hashIteratorBase_init
  END TYPE t_HashIteratorBase

  ! Provides sequential access to all entries of a hash table
  !
  ! <<< WARNING >>>:
  !   iterators are invalidated by setEntry(), removeEntry(), and
  !   destruct() calls on the corresponding hash table
  !
  TYPE, EXTENDS(t_HashIteratorBase) :: t_HashIterator
  CONTAINS
    PROCEDURE :: nextEntry => hashIterator_nextEntry ! returns .TRUE. IF the operation was successful
  END TYPE t_HashIterator

  ! Provides sequential access to all entries of a hash table.
  !
  ! The "const" iterator is the read-only variant: When using a
  ! "t_HashIterator" for traversal of hash table entries, it is
  ! possible to change the contents of the hash table. There are
  ! situations, however, when the user needs to loop over the entries
  ! of a dictionary that is INTENT(IN).
  TYPE, EXTENDS(t_HashIteratorBase) :: t_const_hashIterator
  CONTAINS
    PROCEDURE :: nextEntry => const_hashIterator_nextEntry  ! returns .TRUE. IF the operation was successful
  END TYPE t_const_hashIterator


  TYPE :: t_HashEntryPtr
    TYPE(t_HashEntry), POINTER :: ptr
  END TYPE t_HashEntryPtr
  TYPE :: t_HashEntry
    CLASS(t_Destructible), POINTER :: key, val
    TYPE(t_HashEntryPtr) :: next
    INTEGER(C_INT32_T) :: hash
  END TYPE t_HashEntry
  
CONTAINS

  !--------------------------------------------------------------------------
  !> (Empty) destructor for hash table entries.
  !
  SUBROUTINE scalarVal_destruct(me)
    CLASS(t_scalarVal), INTENT(inout) :: me
  END SUBROUTINE scalarVal_destruct


  !--------------------------------------------------------------------------
  !> String hash key as proposed by Daniel J. Bernstein, case sensitive
  !
  INTEGER(C_INT32_T) FUNCTION storage_hashKey_DJB_cs(key) RESULT(res)
    CLASS(t_Destructible), POINTER, INTENT(in) :: key
    ! Local
    INTEGER :: i
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//":storage_hashKey_DJB"
    INTEGER(C_INT64_T) :: t 
    
    res = 5381
    
    SELECT TYPE(key)
    TYPE IS(t_stringVal)
      DO i=1,LEN(TRIM(key%stringVal))
        t = ISHFT(res,5)
        res = INT( MOD( (t + res) + ICHAR(key%stringVal(i:i)), p32 ), C_INT32_T)
      END DO
    CLASS DEFAULT
      CALL finish(routine, "Unknown type for key.")
    END SELECT
  END FUNCTION storage_hashKey_DJB_cs


  !--------------------------------------------------------------------------
  !> String hash key as proposed by Daniel J. Bernstein, case-insensitive
  !
  INTEGER(C_INT32_T) FUNCTION storage_hashKey_DJB_ci(key) RESULT(res)
    CLASS(t_Destructible), POINTER, INTENT(in) :: key
    ! Local
    INTEGER :: i, ic
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//":storage_hashKey_DJB"
    INTEGER(C_INT64_T) :: t 
  
    INTEGER, PARAMETER :: idel = ICHAR('a')-ICHAR('A')
    INTEGER, PARAMETER :: ia = ICHAR('A')
    INTEGER, PARAMETER :: iz = ICHAR('Z')

    res = 5381

    SELECT TYPE(key)
    TYPE IS(t_stringVal)
      DO i=1,LEN(TRIM(key%stringVal))
        t = ISHFT(res,5)
        ic = ICHAR(key%stringVal(i:i))
        IF (ic >= ia .AND. ic <= iz)   ic = ic + idel
        res = INT( MOD( (t + res) + ic, p32 ), C_INT32_T)
      END DO
    CLASS DEFAULT
      CALL finish(routine, "Unknown type for key.")
    END SELECT
  END FUNCTION storage_hashKey_DJB_ci


  !--------------------------------------------------------------------------
  !> Simple string comparison (case sensitive)
  LOGICAL FUNCTION storage_equalKeysFunction_cs(keyA, keyB) RESULT(res)
    CLASS(t_Destructible), POINTER, INTENT(in) :: keyA, keyB
    ! Local
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//":storage_equalKeysFunction"
  
    SELECT TYPE(keyA)
    TYPE IS(t_stringVal)
      SELECT TYPE(keyB)
      TYPE IS(t_stringVal)
        res =.FALSE.
        IF(TRIM(keyA%stringVal) == TRIM(keyB%stringVal)) res =.TRUE.
      CLASS DEFAULT
        CALL finish(routine, "Unknown type for keyB.")
      END SELECT
    CLASS DEFAULT
      CALL finish(routine, "Unknown type for keyA.")
    END SELECT
  END FUNCTION storage_equalKeysFunction_cs


  !--------------------------------------------------------------------------
  !> Simple string comparison (case insensitive)
  LOGICAL FUNCTION storage_equalKeysFunction_ci(keyA, keyB) RESULT(res)
    CLASS(t_Destructible), POINTER, INTENT(in) :: keyA, keyB
    ! Local
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//":storage_equalKeysFunction"

    SELECT TYPE(keyA)
    TYPE IS(t_stringVal)
      SELECT TYPE(keyB)
      TYPE IS(t_stringVal)
        res =.FALSE.
        IF(TRIM(tolower(keyA%stringVal)) == TRIM(tolower(keyB%stringVal))) res =.TRUE.
      CLASS DEFAULT
        CALL finish(routine, "Unknown type for keyB.")
      END SELECT
    CLASS DEFAULT
      CALL finish(routine, "Unknown type for keyA.")
    END SELECT
  END FUNCTION storage_equalKeysFunction_ci

  ! --------------------------------------------------------------------------------

  
  FUNCTION hashTable_make(hashFunction, compareFunction) RESULT(resultVar)
    PROCEDURE(f_hashFunction) :: hashFunction
    PROCEDURE(f_equalKeysFunction) :: compareFunction
    CLASS(t_HashTable_base), POINTER :: resultVar

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":hashTable_make"
    INTEGER :: error, i

    ALLOCATE(t_HashTable::resultVar, STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

    resultVar%getHash => hashFunction
    resultVar%equalKeys => compareFunction
    resultVar%entryCount = 0
    resultVar%hashBits = 5

    ALLOCATE(resultVar%table(2**resultVar%hashBits), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    DO i = 1, 2**resultVar%hashBits
      resultVar%table(i)%ptr => NULL()
    END DO
  END FUNCTION hashTable_make

  FUNCTION hashTable_findBin(me, hash) RESULT(resultVar)
    CLASS(t_HashTable_Base), INTENT(IN) :: me
    INTEGER(C_INT32_T), VALUE :: hash
    TYPE(t_HashEntryPtr), POINTER :: resultVar

    INTEGER(C_INT32_T) :: reducedHash, i

    IF(hash < 0) hash = NOT(hash)   !fortran has no unsigned types
    reducedHash = 0
    DO i = 1, (31 + me%hashBits - 1)/me%hashBits
      reducedHash = IEOR(reducedHash, hash)
      hash = ISHFT(hash, -me%hashBits)
    END DO
    reducedHash = IAND(reducedHash, 2**me%hashBits - 1) + 1
    resultVar => me%table(reducedHash)
  END FUNCTION hashTable_findBin

  SUBROUTINE hashTable_growTable(me)
    CLASS(t_HashTable_Base), INTENT(INOUT) :: me

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
    CLASS(t_HashTable_Base), INTENT(INOUT) :: me
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
          CALL curEntry%val%destruct()
          DEALLOCATE(curEntry%val)
          DEALLOCATE(curEntry)
          CYCLE
        END IF
      END IF
      iterator => curEntry%next    ! point to next entry
    END DO
  END SUBROUTINE hashTable_removeFromList

  ! The hash table takes possession of both the key and the val and will DEALLOCATE() them eventually.
  SUBROUTINE hashTable_setEntry(me, key, val)
    CLASS(t_HashTable), INTENT(INOUT) :: me
    CLASS(t_Destructible), POINTER, INTENT(IN) :: key, val

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":hashTable_setEntry"
    TYPE(t_HashEntry), POINTER :: newEntry
    TYPE(t_HashEntryPtr), POINTER :: bin
    INTEGER :: error

    ! Prepare the new entry.
    ALLOCATE(newEntry, STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    newEntry%key => key
    newEntry%val => val
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

  SUBROUTINE hashTable_setEntry_Base(me, key, val)
    CLASS(t_HashTable_Base), INTENT(INOUT)     :: me
    CLASS(t_Destructible), POINTER, INTENT(IN) :: key, val
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//":hashTable_setEntry_Base"

    CALL finish(routine, "setEntry for uninitialized hash tables not possible")
  END SUBROUTINE hashTable_setEntry_Base

  SUBROUTINE hashTable_removeEntry(me, key)
    CLASS(t_HashTable), INTENT(INOUT) :: me
    CLASS(t_Destructible), POINTER, INTENT(IN) :: key

    INTEGER(C_INT32_T) :: hash
    TYPE(t_HashEntryPtr), POINTER :: bin

    hash = me%getHash(key)
    bin => me%findBin(hash)
    CALL me%removeFromList(bin, key, hash)
  END SUBROUTINE hashTable_removeEntry

  SUBROUTINE hashTable_removeEntry_Base(me, key)
    CLASS(t_HashTable_Base), INTENT(INOUT)     :: me
    CLASS(t_Destructible), POINTER, INTENT(IN) :: key
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//":hashTable_removeEntry_Base"

    CALL finish(routine, "removeEntry for uninitialized hash tables not possible")
  END SUBROUTINE hashTable_removeEntry_Base

  FUNCTION hashTable_getEntry(me, key) RESULT(resultVar)
    CLASS(t_HashTable), INTENT(IN) :: me
    CLASS(t_Destructible), POINTER, INTENT(IN) :: key
    CLASS(t_Destructible), POINTER :: resultVar

    INTEGER(C_INT32_T) :: hash
    TYPE(t_HashEntryPtr), POINTER :: bin
    TYPE(t_HashEntry), POINTER :: curEntry

    resultVar => NULL()
    hash = me%getHash(key)
    bin => me%findBin(hash)
    curEntry => bin%ptr
    DO WHILE(ASSOCIATED(curEntry))
      IF(curEntry%hash == hash) THEN
        IF(me%equalKeys(curEntry%key, key)) THEN
          resultVar => curEntry%val
          RETURN
        END IF
      END IF
      curEntry => curEntry%next%ptr
    END DO
  END FUNCTION hashTable_getEntry

  FUNCTION hashTable_getEntry_Base(me, key) RESULT(resultVar)
    CLASS(t_HashTable_Base), INTENT(IN)        :: me
    CLASS(t_Destructible), POINTER, INTENT(IN) :: key
    CLASS(t_Destructible), POINTER :: resultVar
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//":hashTable_getEntry_Base"

    resultVar => NULL()
    CALL finish(routine, "getEntry for uninitialized hash tables not possible")
  END FUNCTION hashTable_getEntry_Base

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
        CALL curEntry%val%destruct()
        DEALLOCATE(curEntry%val)
        DEALLOCATE(curEntry)
        curEntry => nextEntry
      END DO
    END DO
    DEALLOCATE(me%table)
  END SUBROUTINE hashTable_destruct

  SUBROUTINE hashTable_destruct_Base(me)
    CLASS(t_HashTable_Base), INTENT(INOUT) :: me
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//":hashTable_destruct_Base"

    CALL finish(routine, "destruct for uninitialized hash tables not possible")
  END SUBROUTINE hashTable_destruct_Base

  INTEGER FUNCTION hashTable_nentries(me)
    CLASS(t_HashTable), INTENT(IN) :: me
    hashTable_nentries = me%entryCount
  END FUNCTION hashTable_nentries

  INTEGER FUNCTION hashTable_nentries_Base(me)
    CLASS(t_HashTable_Base), INTENT(IN) :: me
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//":hashTable_nentries_Base"
    hashTable_nentries_Base = -1
    CALL finish(routine, "nentries for uninitialized hash tables not possible")
  END FUNCTION hashTable_nentries_Base

  !> Initialize iterator for sequential access to all entries of a
  !  hash table.
  SUBROUTINE hashIteratorBase_init(me, table)
    CLASS(t_HashIteratorBase), INTENT(INOUT)     :: me
    CLASS(t_HashTable_Base), POINTER, INTENT(IN) :: table

    me%table => table
    me%curBin = 0   ! will be incremented in the first nextEntry() call
    me%curEntry => NULL()
  END SUBROUTINE hashIteratorBase_init


  LOGICAL FUNCTION hashIterator_nextEntry(me, key, val) RESULT(resultVar)
    CLASS(t_HashIterator), INTENT(INOUT) :: me
    CLASS(t_Destructible), POINTER, INTENT(INOUT) :: key, val

    key => NULL()
    val => NULL()
    resultVar = .FALSE.

    ! try to ADVANCE within the current bin
    IF(ASSOCIATED(me%curEntry)) me%curEntry => me%curEntry%next%ptr

    ! search for the next entry
    DO
      ! check whether we have found the next entry
      IF(ASSOCIATED(me%curEntry)) THEN
        key => me%curEntry%key
        val => me%curEntry%val
        resultVar = .TRUE.
        RETURN
      END IF

      ! check whether we can ADVANCE ANY further
      IF(me%curBin >= SIZE(me%table%table, 1)) RETURN

      ! ADVANCE to the next bin
      me%curBin = me%curBin + 1
      me%curEntry => me%table%table(me%curBin)%ptr
    END DO
  END FUNCTION hashIterator_nextEntry


  LOGICAL FUNCTION const_hashIterator_nextEntry(me, key, val) RESULT(resultVar)
    CLASS(t_const_hashIterator), INTENT(INOUT) :: me
    CLASS(t_Destructible), ALLOCATABLE, INTENT(OUT) :: key, val

    resultVar = .FALSE.

    ! try to advance within the current bin
    IF(ASSOCIATED(me%curEntry)) me%curEntry => me%curEntry%next%ptr

    ! search for the next entry
    DO
      ! check whether we have found the next entry
      IF(ASSOCIATED(me%curEntry)) THEN
        ALLOCATE(key, source=me%curEntry%key)
        ALLOCATE(val, source=me%curEntry%val)
        resultVar = .TRUE.
        RETURN
      END IF

      ! check whether we can ADVANCE ANY further
      IF(me%curBin >= SIZE(me%table%table, 1)) RETURN

      ! advance to the next bin
      me%curBin = me%curBin + 1
      me%curEntry => me%table%table(me%curBin)%ptr
    END DO
  END FUNCTION const_hashIterator_nextEntry

END MODULE mo_hash_table
