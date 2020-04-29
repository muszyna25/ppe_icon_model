!>
!! This module provides a simple key value storage based on string comparisons
!! using the generic hash tables provided by mo_hash_table. 
!! As method to calculate hash keys, the DJB algorithm is used.
!!
!!
!! @par Revision History
!! Initial release by Daniel Rieger, KIT (2016-12-13)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_storage

  USE ISO_C_BINDING,                    ONLY: C_INT32_T, C_INT64_T
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish, message, message_text
  USE mo_util_string,                   ONLY: tolower
  USE mo_hash_table,                    ONLY: t_HashTable, hashTable_make, t_DestructibleContainer
  USE mo_fortran_tools,                 ONLY: t_Destructible
  USE mo_impl_constants,                ONLY: SUCCESS

  IMPLICIT NONE

  PRIVATE

  INTEGER,PARAMETER :: stringVal_len= 120  !< Maximum allowed length for string keys or values

  INTEGER(C_INT32_T), PARAMETER :: i32_t = 0
  INTEGER(C_INT64_T), PARAMETER :: p32 = INT(HUGE(i32_t), C_INT64_T) + 1 


  TYPE, EXTENDS(t_Destructible) :: t_scalarVal
    PRIVATE
  CONTAINS
    PROCEDURE :: destruct => scalarVal_destruct
  END TYPE t_scalarVal

  TYPE, EXTENDS(t_scalarVal) :: t_stringVal
    PRIVATE
    CHARACTER(LEN=stringVal_len) :: stringVal
  END TYPE t_stringVal

  TYPE, EXTENDS(t_scalarVal) :: t_realVal
    PRIVATE
    REAL(wp) :: realVal
  END TYPE t_realVal

  TYPE, EXTENDS(t_scalarVal) :: t_intVal
    PRIVATE
    INTEGER :: intVal
  END TYPE t_intVal

  TYPE, EXTENDS(t_scalarVal) :: t_logVal
    PRIVATE
    LOGICAL :: logVal
  END TYPE t_logVal

  TYPE t_storage
    TYPE(t_HashTable),PRIVATE :: &
      &  container
    CONTAINS
      PROCEDURE, PUBLIC  :: init     => init_storage
      PROCEDURE, PUBLIC  :: destruct => destruct_storage
      PROCEDURE, PRIVATE :: put_real
      PROCEDURE, PRIVATE :: put_int
      PROCEDURE, PRIVATE :: put_string
      PROCEDURE, PRIVATE :: put_logical
      GENERIC,   PUBLIC  :: put      => put_real, put_int, put_string, put_logical
      PROCEDURE, PRIVATE :: get_real
      PROCEDURE, PRIVATE :: get_int
      PROCEDURE, PRIVATE :: get_string
      PROCEDURE, PRIVATE :: get_logical
      GENERIC,   PUBLIC  :: get      => get_real, get_int, get_string, get_logical
      PROCEDURE, PUBLIC  :: dump     => dump_storage
  END TYPE t_storage
  
  CHARACTER(len=*), PARAMETER :: modname = 'mo_storage'

  PUBLIC :: t_storage
  
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE destruct_storage(this_storage)
  CLASS(t_storage),INTENT(inout)    :: &
    &  this_storage

  CALL this_storage%container%destruct

END SUBROUTINE destruct_storage
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE init_storage(this_storage,lcase_sensitivity)
  CLASS(t_storage),INTENT(inout)    :: &
    &  this_storage
  LOGICAL,INTENT(in),OPTIONAL       :: &
    &  lcase_sensitivity
  !
  TYPE(t_HashTable), POINTER :: tmp_container

  IF (PRESENT(lcase_sensitivity))THEN
    IF (lcase_sensitivity) THEN
      tmp_container => hashTable_make(storage_hashKey_DJB_cs, storage_equalKeysFunction_cs)
    ELSE
      tmp_container => hashTable_make(storage_hashKey_DJB_ci, storage_equalKeysFunction_ci)
    ENDIF
  ELSE
    ! Not present: case insensitive by default
    tmp_container => hashTable_make(storage_hashKey_DJB_ci, storage_equalKeysFunction_ci)
  ENDIF
 this_storage%container = tmp_container
 DEALLOCATE(tmp_container)
END SUBROUTINE init_storage
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE dump_storage(this_storage,opt_label_in)
  CLASS(t_storage),INTENT(inout)           :: &
    &  this_storage
  CHARACTER(LEN = *), INTENT(IN),OPTIONAL  :: &
    &  opt_label_in
  CHARACTER(LEN = *), PARAMETER            :: routine = modname//":dump_storage"
  CHARACTER(LEN = 3)                       :: i_string
  CHARACTER(LEN = 128)                     :: opt_label
  INTEGER :: i
  CLASS(t_DestructibleContainer),POINTER :: keys(:),vals(:)

  opt_label = ""
  IF(PRESENT(opt_label_in)) THEN
    opt_label = TRIM(opt_label_in)
  ENDIF

  CALL this_storage%container%dumpPrepare(keys,vals)

  IF(.NOT.(ASSOCIATED(keys).AND.ASSOCIATED(vals))) THEN
    CALL message(routine,"Storage structure is empty")
  ELSE
    WRITE(message_text,"(A)") "START STORAGE DUMP ("//TRIM(opt_label)//")"
    CALL message(routine,message_text)
    DO i=1,SIZE(keys)
      WRITE(message_text,"(A,I2.2,A)") ">> ",i,": key="
      SELECT TYPE(p_key => keys(i)%dest)
        TYPE IS(t_stringVal)
          WRITE(message_text,"(A,A,A20,A)") TRIM(message_text),"'",p_key%stringVal,"'"
        CLASS DEFAULT
          WRITE(message_text,"(A,A20)") TRIM(message_text),"<invalid>"
      END SELECT
      WRITE(message_text,"(A)") TRIM(message_text)//" ; value="
      SELECT TYPE(p_val => vals(i)%dest)
        TYPE IS(t_stringVal)
          WRITE(message_text,"(A,A)") TRIM(message_text),"'"//TRIM(p_val%stringVal)//"'"
        TYPE IS(t_intVal)
          WRITE(message_text,"(A,I0)") TRIM(message_text),p_val%intVal
        TYPE IS(t_realVal)
          WRITE(message_text,"(A,ES13.4)") TRIM(message_text),p_val%realVal
        TYPE IS(t_logVal)
          WRITE(message_text,"(A,L5)") TRIM(message_text),p_val%logVal
        CLASS DEFAULT
          WRITE(message_text,"(A,A)") TRIM(message_text),"<invalid>"
      END SELECT
      CALL message(routine,message_text)
    END DO
    CALL message(routine,"END OF STORAGE DUMP ("//TRIM(opt_label)//")")
  END IF
END SUBROUTINE dump_storage
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE scalarVal_destruct(me)
  CLASS(t_scalarVal), INTENT(inout) :: me
END SUBROUTINE scalarVal_destruct
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE put_real(this_storage, key, value)
  CLASS(t_storage),INTENT(inout) :: this_storage
  CHARACTER(LEN=*),INTENT(in)    :: key
  REAL(wp)        ,INTENT(in)    :: value
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
  CLASS(t_storage),INTENT(inout) :: this_storage
  CHARACTER(LEN=*),INTENT(in)    :: key
  INTEGER         ,INTENT(in)    :: value
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
  CLASS(t_storage),INTENT(inout) :: this_storage
  CHARACTER(LEN=*),INTENT(in)    :: key
  CHARACTER(LEN=*),INTENT(in)    :: value
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
SUBROUTINE put_logical(this_storage, key, value)
  CLASS(t_storage),INTENT(inout) :: this_storage
  CHARACTER(LEN=*),INTENT(in)    :: key
  LOGICAL         ,INTENT(in)    :: value
! Local
  CLASS(t_Destructible),POINTER :: p_key
  CLASS(t_Destructible),POINTER :: p_value

  ALLOCATE(t_stringVal :: p_key  )
  ALLOCATE(t_logVal    :: p_value)

  SELECT TYPE(p_key)
    TYPE IS(t_stringVal)
      p_key%stringVal = TRIM(key)
  END SELECT
  
  SELECT TYPE(p_value)
    TYPE IS(t_logVal)
      p_value%logVal = value
  END SELECT

  CALL this_storage%container%setEntry(p_key, p_value)

END SUBROUTINE put_logical
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_real(this_storage, key, value, ierror)
! Implemented as subroutine rather than function as the interface "get" is not distinguishable otherwise.
  CLASS(t_storage),INTENT(in)    :: this_storage
  CHARACTER(LEN=*),INTENT(in)    :: key
  REAL(wp)        ,INTENT(out)   :: value
  INTEGER,OPTIONAL,INTENT(out)   :: ierror
! Local
  CHARACTER(LEN=*), PARAMETER    :: routine = modname//":get_real"
  CLASS(t_Destructible),POINTER  :: p_key
  CLASS(t_Destructible),POINTER  :: p_value

  ALLOCATE(t_stringVal :: p_key)
  IF (PRESENT(ierror)) ierror = SUCCESS

  SELECT TYPE(p_key)
    TYPE IS(t_stringVal)
      p_key%stringVal = TRIM(key)
  END SELECT

  p_value => this_storage%container%getEntry(p_key)

  IF (ASSOCIATED(p_value)) THEN
    SELECT TYPE (p_value)
      TYPE IS(t_realVal)
        value = p_value%realVal
      CLASS DEFAULT
        CALL finish(routine, "Wrong return type for "//TRIM(key)//".")
    END SELECT
  ELSE
    IF (PRESENT(ierror)) ierror = 1
  ENDIF

  DEALLOCATE(p_key)

END SUBROUTINE get_real
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_int(this_storage, key, value, ierror)
! Implemented as subroutine rather than function as the interface "get" is not distinguishable otherwise.
  CLASS(t_storage),INTENT(in)    :: this_storage
  CHARACTER(LEN=*),INTENT(in)    :: key
  INTEGER         ,INTENT(out)   :: value
  INTEGER,OPTIONAL,INTENT(out)   :: ierror
! Local
  CHARACTER(LEN=*), PARAMETER    :: routine = modname//":get_int"
  CLASS(t_Destructible),POINTER  :: p_key
  CLASS(t_Destructible),POINTER  :: p_value

  ALLOCATE(t_stringVal :: p_key)
  IF (PRESENT(ierror)) ierror = SUCCESS

  SELECT TYPE(p_key)
    TYPE IS(t_stringVal)
      p_key%stringVal = TRIM(key)
  END SELECT

  p_value => this_storage%container%getEntry(p_key)

  IF (ASSOCIATED(p_value)) THEN
    SELECT TYPE (p_value)
      TYPE IS(t_intVal)
        value = p_value%intVal
      CLASS DEFAULT
        CALL finish(routine, "Wrong return type for "//TRIM(key)//".")
    END SELECT
  ELSE
    IF (PRESENT(ierror)) ierror = 1
  ENDIF

  DEALLOCATE(p_key)

END SUBROUTINE get_int
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_string(this_storage, key, value, ierror)
! Implemented as subroutine rather than function as the interface "get" is not distinguishable otherwise.
  CLASS(t_storage),INTENT(in)    :: this_storage
  CHARACTER(LEN=*),INTENT(in)    :: key
  CHARACTER(LEN=*),INTENT(out)   :: value
  INTEGER,OPTIONAL,INTENT(out)   :: ierror
! Local
  CHARACTER(LEN=*), PARAMETER    :: routine = modname//":get_string"
  CLASS(t_Destructible),POINTER  :: p_key
  CLASS(t_Destructible),POINTER  :: p_value

  ALLOCATE(t_stringVal :: p_key)
  IF (PRESENT(ierror)) ierror = SUCCESS

  SELECT TYPE(p_key)
    TYPE IS(t_stringVal)
      p_key%stringVal = TRIM(key)
  END SELECT

  p_value => this_storage%container%getEntry(p_key)

  IF (ASSOCIATED(p_value)) THEN
    SELECT TYPE (p_value)
      TYPE IS(t_stringVal)
        value = TRIM(p_value%stringVal)
      CLASS DEFAULT
        CALL finish(routine, "Wrong return type for "//TRIM(key)//".")
    END SELECT
  ELSE
    IF (PRESENT(ierror)) ierror = 1
  ENDIF

  DEALLOCATE(p_key)

END SUBROUTINE get_string
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_logical(this_storage, key, value, ierror)
! Implemented as subroutine rather than function as the interface "get" is not distinguishable otherwise.
  CLASS(t_storage),INTENT(in)    :: this_storage
  CHARACTER(LEN=*),INTENT(in)    :: key
  LOGICAL         ,INTENT(out)   :: value
  INTEGER,OPTIONAL,INTENT(out)   :: ierror
! Local
  CHARACTER(LEN=*), PARAMETER    :: routine = modname//":get_logical"
  CLASS(t_Destructible),POINTER  :: p_key
  CLASS(t_Destructible),POINTER  :: p_value

  ALLOCATE(t_stringVal :: p_key)
  IF (PRESENT(ierror)) ierror = SUCCESS

  SELECT TYPE(p_key)
    TYPE IS(t_stringVal)
      p_key%stringVal = TRIM(key)
  END SELECT

  p_value => this_storage%container%getEntry(p_key)

  IF (ASSOCIATED(p_value)) THEN
    SELECT TYPE (p_value)
      TYPE IS(t_logVal)
        value = p_value%logVal
      CLASS DEFAULT
        CALL finish(routine, "Wrong return type for "//TRIM(key)//".")
    END SELECT
  ELSE
    IF (PRESENT(ierror)) ierror = 1
  ENDIF

  DEALLOCATE(p_key)

END SUBROUTINE get_logical
!!
!!-------------------------------------------------------------------------
!!
INTEGER(C_INT32_T) FUNCTION storage_hashKey_DJB_cs(key) RESULT(result)
! Create hash key as proposed by Daniel J. Bernstein
  CLASS(t_Destructible), POINTER, INTENT(in) :: key
! Local
  integer :: i
  CHARACTER(LEN=*), PARAMETER    :: routine = modname//":storage_hashKey_DJB"
  INTEGER(C_INT64_T) :: t 

  result = 5381

  SELECT TYPE(key)
    TYPE IS(t_stringVal)
      do i=1,len(TRIM(key%stringVal))
        t = ISHFT(RESULT,5)
        RESULT = MOD( (t + RESULT) + ICHAR(key%stringVal(i:i)), p32 ) 
      end do
    CLASS DEFAULT
      CALL finish(routine, "Unknown type for key.")
  END SELECT
END FUNCTION storage_hashKey_DJB_cs
!!
!!-------------------------------------------------------------------------
!!
INTEGER(C_INT32_T) FUNCTION storage_hashKey_DJB_ci(key) RESULT(result)
! Create hash key as proposed by Daniel J. Bernstein
  CLASS(t_Destructible), POINTER, INTENT(in) :: key
! Local
  integer :: i
  CHARACTER(LEN=*), PARAMETER    :: routine = modname//":storage_hashKey_DJB"
  INTEGER(C_INT64_T) :: t 

  result = 5381

  SELECT TYPE(key)
    TYPE IS(t_stringVal)
      do i=1,len(TRIM(key%stringVal))
        t = ISHFT(RESULT,5)
        RESULT = MOD( (t + RESULT) + ICHAR(tolower(key%stringVal(i:i))), p32 ) 
      end do
    CLASS DEFAULT
      CALL finish(routine, "Unknown type for key.")
  END SELECT
END FUNCTION storage_hashKey_DJB_ci
!!
!!-------------------------------------------------------------------------
!!
LOGICAL FUNCTION storage_equalKeysFunction_cs(keyA, keyB) RESULT(result)
! Simple string comparison (case sensitive)
  CLASS(t_Destructible), POINTER, INTENT(in) :: keyA, keyB
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
END FUNCTION storage_equalKeysFunction_cs
!!
!!-------------------------------------------------------------------------
!!
LOGICAL FUNCTION storage_equalKeysFunction_ci(keyA, keyB) RESULT(result)
! Simple string comparison (case insensitive)
  CLASS(t_Destructible), POINTER, INTENT(in) :: keyA, keyB
! Local
  CHARACTER(LEN=*), PARAMETER    :: routine = modname//":storage_equalKeysFunction"

  SELECT TYPE(keyA)
    TYPE IS(t_stringVal)
      SELECT TYPE(keyB)
        TYPE IS(t_stringVal)
          result =.FALSE.
          IF(TRIM(tolower(keyA%stringVal)) == TRIM(tolower(keyB%stringVal))) result =.TRUE.
        CLASS DEFAULT
          CALL finish(routine, "Unknown type for keyB.")
      END SELECT
    CLASS DEFAULT
      CALL finish(routine, "Unknown type for keyA.")
  END SELECT
END FUNCTION storage_equalKeysFunction_ci
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_storage
