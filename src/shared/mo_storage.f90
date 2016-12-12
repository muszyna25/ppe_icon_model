!<
!! mo_storage
!! This module provides a key value storage
!!
!! Author: Daniel Rieger, KIT
!! Initial Release: 2016-12-12
!!
!! Modifications: 
!! YYYY-MM-DD: <name>,<institution>
!! - <description>
!!
!! Copyright and License
!! This code is subject to the KIT-Software-License-Agreement in
!! its most recent form.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!! This software is provided for non-commercial use only.
!!
!! See the LICENSE conditions:
!! http://icon-art.imk-tro.kit.edu/
!! Institute of Meteorology and Climate Research, KIT, Karlsruhe
!>
MODULE mo_storage

  USE ISO_C_BINDING,                    ONLY: C_INT32_T
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish
  USE mo_hash_table,                    ONLY: t_HashTable, hashTable_make
  USE mo_fortran_tools,                 ONLY: t_Destructible

  IMPLICIT NONE

  PRIVATE

  INTEGER,PARAMETER :: stringVal_len= 120  !< Maximum allowed length for string keys or values

  TYPE, EXTENDS(t_Destructible) :: t_stringVal
    PRIVATE
    CHARACTER(LEN=stringVal_len) :: stringVal
  CONTAINS
    PROCEDURE :: destruct => stringVal_destruct
  END TYPE t_stringVal

  TYPE, EXTENDS(t_Destructible) :: t_realVal
    PRIVATE
    REAL(wp) :: realVal
  CONTAINS
    PROCEDURE :: destruct => realVal_destruct
  END TYPE t_realVal

  TYPE, EXTENDS(t_Destructible) :: t_intVal
    PRIVATE
    INTEGER :: intVal
  CONTAINS
    PROCEDURE :: destruct => intVal_destruct
  END TYPE t_intVal

  TYPE t_storage
    TYPE(t_HashTable),PRIVATE :: &
      &  container
    CONTAINS
      PROCEDURE :: init => init_storage
      !PROCEDURE :: free => free_table
      PROCEDURE, PRIVATE :: put_real
      PROCEDURE, PRIVATE :: put_int
      PROCEDURE, PRIVATE :: put_string
      GENERIC,   PUBLIC  :: put => put_real, put_int, put_string
      PROCEDURE, PRIVATE :: get_real
      PROCEDURE, PRIVATE :: get_int
      PROCEDURE, PRIVATE :: get_string
      GENERIC,   PUBLIC  :: get => get_real, get_int, get_string
  END TYPE t_storage
  
  CHARACTER(len=*), PARAMETER :: modname = 'mo_storage'

  PUBLIC :: t_storage
  
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE init_storage(this_storage)
  CLASS(t_storage),INTENT(inout)    :: &
    &  this_storage

  this_storage%container = hashTable_make(storage_hashKey_DJB, storage_equalKeysFunction)

END SUBROUTINE init_storage
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE stringVal_destruct(me)
  CLASS(t_stringVal), INTENT(INOUT) :: me
END SUBROUTINE stringVal_destruct

SUBROUTINE realVal_destruct(me)
  CLASS(t_realVal), INTENT(INOUT) :: me
END SUBROUTINE realVal_destruct

SUBROUTINE intVal_destruct(me)
  CLASS(t_intVal), INTENT(INOUT) :: me
END SUBROUTINE intVal_destruct
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE put_real(this_storage, key, value)
  CLASS(t_storage),INTENT(INOUT) :: this_storage
  CHARACTER(LEN=*),INTENT(IN)    :: key
  REAL(wp)        ,INTENT(IN)    :: value
! Local
  CLASS(t_Destructible),POINTER :: p_key
  CLASS(t_Destructible),POINTER :: p_value

  ALLOCATE(t_stringVal :: p_key  )
  ALLOCATE(t_realVal   :: p_value)

  SELECT TYPE(p_key)
    TYPE IS(t_stringVal)
      p_key%stringVal = TRIM(key)
  END SELECT
  
  SELECT TYPE(p_value)
    TYPE IS(t_realVal)
      p_value%realVal = value
  END SELECT

  CALL this_storage%container%setEntry(p_key, p_value)

END SUBROUTINE put_real
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE put_int(this_storage, key, value)
  CLASS(t_storage),INTENT(INOUT) :: this_storage
  CHARACTER(LEN=*),INTENT(IN)    :: key
  INTEGER         ,INTENT(IN)    :: value
! Local
  CLASS(t_Destructible),POINTER :: p_key
  CLASS(t_Destructible),POINTER :: p_value

  ALLOCATE(t_stringVal :: p_key  )
  ALLOCATE(t_intVal    :: p_value)

  SELECT TYPE(p_key)
    TYPE IS(t_stringVal)
      p_key%stringVal = TRIM(key)
  END SELECT
  
  SELECT TYPE(p_value)
    TYPE IS(t_intVal)
      p_value%intVal = value
  END SELECT

  CALL this_storage%container%setEntry(p_key, p_value)

END SUBROUTINE put_int
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE put_string(this_storage, key, value)
  CLASS(t_storage),INTENT(INOUT) :: this_storage
  CHARACTER(LEN=*),INTENT(IN)    :: key
  CHARACTER(LEN=*),INTENT(IN)    :: value
! Local
  CLASS(t_Destructible),POINTER  :: p_key
  CLASS(t_Destructible),POINTER  :: p_value

  ALLOCATE(t_stringVal :: p_key  )
  ALLOCATE(t_stringVal :: p_value)

  SELECT TYPE(p_key)
    TYPE IS(t_stringVal)
      p_key%stringVal = TRIM(key)
  END SELECT
  
  SELECT TYPE(p_value)
    TYPE IS(t_stringVal)
      p_value%stringVal = TRIM(value)
  END SELECT

  CALL this_storage%container%setEntry(p_key, p_value)

END SUBROUTINE put_string
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_real(this_storage, key, value)
  CLASS(t_storage),INTENT(IN)    :: this_storage
  CHARACTER(LEN=*),INTENT(IN)    :: key
  REAL(wp)        ,INTENT(OUT)   :: value
! Local
  CHARACTER(LEN=*), PARAMETER    :: routine = modname//":get_real"
  CLASS(t_Destructible),POINTER  :: p_key
  CLASS(t_Destructible),POINTER  :: p_value

  ALLOCATE(t_stringVal :: p_key)

  SELECT TYPE(p_key)
    TYPE IS(t_stringVal)
      p_key%stringVal = TRIM(key)
  END SELECT

  p_value => this_storage%container%getEntry(p_key)

  SELECT TYPE (p_value)
    TYPE IS(t_realVal)
      value = p_value%realVal
    CLASS DEFAULT
      CALL finish(routine, "Wrong return type for "//TRIM(key)//".")
  END SELECT

END SUBROUTINE get_real
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_int(this_storage, key, value)
  CLASS(t_storage),INTENT(IN)    :: this_storage
  CHARACTER(LEN=*),INTENT(IN)    :: key
  INTEGER         ,INTENT(OUT)   :: value
! Local
  CHARACTER(LEN=*), PARAMETER    :: routine = modname//":get_int"
  CLASS(t_Destructible),POINTER  :: p_key
  CLASS(t_Destructible),POINTER  :: p_value

  ALLOCATE(t_stringVal :: p_key)

  SELECT TYPE(p_key)
    TYPE IS(t_stringVal)
      p_key%stringVal = TRIM(key)
  END SELECT

  p_value => this_storage%container%getEntry(p_key)

  SELECT TYPE (p_value)
    TYPE IS(t_intVal)
      value = p_value%intVal
    CLASS DEFAULT
      CALL finish(routine, "Wrong return type for "//TRIM(key)//".")
  END SELECT

END SUBROUTINE get_int
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_string(this_storage, key, value)
  CLASS(t_storage),INTENT(IN)    :: this_storage
  CHARACTER(LEN=*),INTENT(IN)    :: key
  CHARACTER(LEN=*),INTENT(OUT)   :: value
! Local
  CHARACTER(LEN=*), PARAMETER    :: routine = modname//":get_string"
  CLASS(t_Destructible),POINTER  :: p_key
  CLASS(t_Destructible),POINTER  :: p_value

  ALLOCATE(t_stringVal :: p_key)

  SELECT TYPE(p_key)
    TYPE IS(t_stringVal)
      p_key%stringVal = TRIM(key)
  END SELECT

  p_value => this_storage%container%getEntry(p_key)

  SELECT TYPE (p_value)
    TYPE IS(t_stringVal)
      value = TRIM(p_value%stringVal)
    CLASS DEFAULT
      CALL finish(routine, "Wrong return type for "//TRIM(key)//".")
  END SELECT

END SUBROUTINE get_string
!!
!!-------------------------------------------------------------------------
!!
INTEGER(C_INT32_T) FUNCTION storage_hashKey_DJB(key) RESULT(result)
! Create hash key as proposed by Daniel J. Bernstein
  CLASS(t_Destructible), POINTER, INTENT(IN) :: key
! Local
  integer :: i
  CHARACTER(LEN=*), PARAMETER    :: routine = modname//":storage_hashKey_DJB"

  result = 5381

  SELECT TYPE(key)
    TYPE IS(t_stringVal)
      do i=1,len(TRIM(key%stringVal))
        result = (ishft(result,5) + result) + ichar(key%stringVal(i:i))
      end do
    CLASS DEFAULT
      CALL finish(routine, "Unknown type for key.")
  END SELECT
END FUNCTION storage_hashKey_DJB
!!
!!-------------------------------------------------------------------------
!!
LOGICAL FUNCTION storage_equalKeysFunction(keyA, keyB) RESULT(result)
! Simple string comparison
  CLASS(t_Destructible), POINTER, INTENT(IN) :: keyA, keyB
! Local
  CHARACTER(LEN=*), PARAMETER    :: routine = modname//":storage_equalKeysFunction"

  SELECT TYPE(keyA)
    TYPE IS(t_stringVal)
      SELECT TYPE(keyB)
        TYPE IS(t_stringVal)
          result =.FALSE.
          IF(TRIM(keyA%stringVal) == TRIM(keyB%stringVal)) result =.TRUE.
        CLASS DEFAULT
          CALL finish(routine, "Unknown type for keyB.")
      END SELECT
    CLASS DEFAULT
      CALL finish(routine, "Unknown type for keyA.")
  END SELECT
END FUNCTION storage_equalKeysFunction
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_storage
