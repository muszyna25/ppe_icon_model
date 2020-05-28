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
MODULE mo_key_value_store
  USE ISO_C_BINDING, ONLY: C_DOUBLE, C_INT, C_INT32_T, C_INT64_T
  USE mo_cdi, ONLY: DATATYPE_FLT64, DATATYPE_INT32, DATATYPE_TXT, CDI_GLOBAL, &
    & cdiInqNatts, cdiInqAtt, cdiInqAttFlt, cdiInqAttInt, cdiInqAttTxt, &
    & cdiDefAttTxt, cdiDefAttFlt, cdiDefAttInt, CDI_UNDEFID, CDI_MAX_NAME, &
    & DATATYPE_LONG
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_hash_table,            ONLY: t_HashTable, hashTable_make, t_HashIterator
  USE mo_impl_constants,        ONLY: SUCCESS
  USE mo_kind,                  ONLY: wp
  USE mo_mpi,                   ONLY: p_get_bcast_role
  USE mo_packed_message,        ONLY: t_PackedMessage
  USE mo_util_string,           ONLY: tolower

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_key_value_store

#ifdef __PGI
  TYPE t_char_workaround
    CHARACTER(:), ALLOCATABLE :: c
  END TYPE t_char_workaround
#endif

  TYPE :: t_key_value_store
    PRIVATE
    LOGICAL, PUBLIC :: is_init = .false.
    LOGICAL :: lcase_sensitive = .true.
    TYPE(t_HashTable), POINTER :: table => NULL()
  CONTAINS
    PROCEDURE, PRIVATE :: put_c => key_value_store_put_c
    PROCEDURE, PRIVATE :: put_r => key_value_store_put_r
    PROCEDURE, PRIVATE :: put_i => key_value_store_put_i
    PROCEDURE, PRIVATE :: put_l => key_value_store_put_l
    GENERIC, PUBLIC :: put => put_c, put_r, put_i, put_l
    PROCEDURE, PRIVATE :: get_c => key_value_store_get_c
    PROCEDURE, PRIVATE :: get_r => key_value_store_get_r
    PROCEDURE, PRIVATE :: get_i => key_value_store_get_i
    PROCEDURE, PRIVATE :: get_l => key_value_store_get_l
    GENERIC, PUBLIC :: get => get_c, get_r, get_i, get_l
    PROCEDURE, PRIVATE :: input_cdi => key_value_store_input_cdi
    PROCEDURE, PRIVATE :: input_pmsg => key_value_store_input_pmsg
    GENERIC, PUBLIC :: input => input_cdi, input_pmsg
    PROCEDURE, PUBLIC :: output => key_value_store_output
    PROCEDURE, PUBLIC :: init => key_value_store_init
    PROCEDURE, PUBLIC :: destruct => key_value_store_destruct
    PROCEDURE, PUBLIC :: bcast => key_value_store_bcast
    PROCEDURE, PUBLIC :: getEntryCount => key_value_store_getEntryCount
  END TYPE t_key_value_store

  CHARACTER(*), PARAMETER :: modname = "mo_key_value_store"
  INTEGER(C_INT64_T), PARAMETER :: p32 = INT(HUGE(0_C_INT32_T), C_INT64_T) + 1_C_INT64_T

CONTAINS

  FUNCTION sel_char(key, routine, err_msg) RESULT(ptr)
    CLASS(*), POINTER, INTENT(in) :: key
    CHARACTER(*), INTENT(IN) :: routine, err_msg
    CHARACTER(:), POINTER :: ptr

    SELECT TYPE(key)
#ifdef __PGI
    TYPE IS(t_char_workaround)
      ptr => key%c
#endif
    TYPE IS(CHARACTER(*))
      ptr => key
    CLASS DEFAULT
      CALL finish(routine, err_msg)
    END SELECT
  END FUNCTION sel_char

  INTEGER(C_INT32_T) FUNCTION text_hash_cs(key) RESULT(hash)
    CLASS(*), POINTER, INTENT(in) :: key
    INTEGER :: i
    CHARACTER(*), PARAMETER :: routine = modname//":text_hash_cs"
    INTEGER(C_INT64_T) :: t
    CHARACTER(:), POINTER :: key_p

    hash = 5381_C_INT32_T
    key_p => sel_char(key, routine, "Unknown type for key.")
    DO i = 1,LEN_TRIM(key_p)
      t = ISHFT(INT(hash, C_INT64_T), 5) + INT(hash, C_INT64_T)
      hash = INT(MOD(t + ICHAR(key_p(i:i), C_INT64_T), p32 ), C_INT32_T)
    END DO
  END FUNCTION text_hash_cs

  INTEGER(C_INT32_T) FUNCTION text_hash_ci(key) RESULT(hash)
    CLASS(*), POINTER, INTENT(in) :: key
    INTEGER :: i
    CHARACTER(*), PARAMETER :: routine = modname//":text_hash_ci"
    INTEGER(C_INT64_T) :: t
    CHARACTER(:), POINTER :: key_p

    hash = 5381_C_INT32_T
    key_p => sel_char(key, routine, "Unknown type for key.")
    DO i = 1,LEN_TRIM(key_p)
      t = ISHFT(INT(hash, C_INT64_T), 5) + INT(hash, C_INT64_T)
      hash = INT(MOD(t + ICHAR(tolower(key_p(i:i)), C_INT64_T), p32 ), C_INT32_T)
    END DO
  END FUNCTION text_hash_ci

  LOGICAL FUNCTION text_isEqual_cs(keyA, keyB) RESULT(is_equal)
    CLASS(*), POINTER, INTENT(in) :: keyA, keyB
    CHARACTER(*), PARAMETER :: routine = modname//":text_isEqual_cs"
    CHARACTER(:), POINTER :: keyA_p, keyB_p

    keyA_p => sel_char(keyA, routine, "Unknown type for keyA.")
    keyB_p => sel_char(keyB, routine, "Unknown type for keyB.")
    is_equal = keyA_p == keyB_p
  END FUNCTION text_isEqual_cs

  LOGICAL FUNCTION text_isEqual_ci(keyA, keyB) RESULT(is_equal)
    CLASS(*), POINTER, INTENT(in) :: keyA, keyB
    CHARACTER(*), PARAMETER :: routine = modname//":text_isEqual_ci"
    CHARACTER(:), POINTER :: keyA_p, keyB_p

    keyA_p => sel_char(keyA, routine, "Unknown type for keyA.")
    keyB_p => sel_char(keyB, routine, "Unknown type for keyB.")
    is_equal = tolower(keyA_p) == tolower(keyB_p)
  END FUNCTION text_isEqual_ci

  SUBROUTINE key_value_store_init(me, cs)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    LOGICAL, INTENT(IN) :: cs
    
    IF (me%is_init) THEN
      IF (cs .NEQV. me%lcase_sensitive .OR. me%table%getEntryCount() .GT. 0) &
        & CALL key_value_store_destruct(me)
    END IF
    IF (.NOT. me%is_init) THEN
      me%lcase_sensitive = cs
      IF (cs) THEN
        me%table => hashTable_make(text_hash_cs, text_isEqual_cs)
      ELSE
        me%table => hashTable_make(text_hash_ci, text_isEqual_ci)
      END IF
      me%is_init = .true.
    END IF
  END SUBROUTINE key_value_store_init

  SUBROUTINE key_value_store_bcast(me, root, comm)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: root, comm
    LOGICAL :: isRoot, isDest
    TYPE(t_PackedMessage) :: pmsg

    CALL p_get_bcast_role(root, comm, isRoot, isDest)
    IF (isRoot) CALL key_value_store_output(me, pmsg=pmsg)
    CALL pmsg%bcast(root, comm)
    IF (isDest) CALL key_value_store_input_pmsg(me, pmsg)
  END SUBROUTINE key_value_store_bcast

#ifdef __PGI
#define ALLOC_ASSIGN(__a) \
    TYPE(t_char_workaround), POINTER :: key_p; __a , POINTER :: val_p; ALLOCATE(key_p, val_p); \
    ALLOCATE(CHARACTER(LEN=LEN_TRIM(key)) :: key_p%c); WRITE(key_p%c, "(a)") TRIM(key); keyObj => key_p; \
    val_p = val; valObj => val_p
#define ALLOC_ASSIGN_C \
    TYPE(t_char_workaround), POINTER :: key_p, val_p; ALLOCATE(key_p, val_p);\
    ALLOCATE(CHARACTER(LEN=LEN_TRIM(key)) :: key_p%c); WRITE(key_p%c, "(a)") TRIM(key); keyObj => key_p; \
    ALLOCATE(CHARACTER(LEN=LEN_TRIM(val)) :: val_p%c); WRITE(val_p%c, "(a)") TRIM(val); valObj => val_p
#else
#define ALLOC_ASSIGN(__a) \
    ALLOCATE(keyObj , SOURCE=TRIM(key)); ALLOCATE(valObj , SOURCE=val)
#define ALLOC_ASSIGN_C ALLOC_ASSIGN(xxxx)
#endif

  SUBROUTINE key_value_store_put_c(me, key, val)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key, val
    CLASS(*), POINTER :: keyObj, valObj

    ALLOC_ASSIGN_C
    CALL me%table%setEntry(keyObj, valObj)
  END SUBROUTINE key_value_store_put_c

  SUBROUTINE key_value_store_put_r(me, key, val)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key
    REAL(wp), INTENT(IN) :: val
    CLASS(*), POINTER :: keyObj, valObj

    ALLOC_ASSIGN(REAL(wp))
    CALL me%table%setEntry(keyObj, valObj)
  END SUBROUTINE key_value_store_put_r

  SUBROUTINE key_value_store_put_i(me, key, val)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key
    INTEGER, INTENT(IN) :: val
    CLASS(*), POINTER :: keyObj, valObj

    ALLOC_ASSIGN(INTEGER)
    CALL me%table%setEntry(keyObj, valObj)
  END SUBROUTINE key_value_store_put_i

  SUBROUTINE key_value_store_put_l(me, key, val)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key
    LOGICAL, INTENT(IN) :: val
    CLASS(*), POINTER :: keyObj, valObj

    ALLOC_ASSIGN(LOGICAL)
    CALL me%table%setEntry(keyObj, valObj)
  END SUBROUTINE key_value_store_put_l

  SUBROUTINE key_value_store_get_internal(me, trimmed_key, opt_err, opt_c, opt_r, opt_i, opt_l)
    CLASS(t_key_value_store), INTENT(IN) :: me
    CHARACTER(*), INTENT(IN), TARGET :: trimmed_key
    CHARACTER(:), INTENT(OUT), ALLOCATABLE, OPTIONAL :: opt_c
    REAL(wp), INTENT(OUT), OPTIONAL :: opt_r
    INTEGER, INTENT(OUT), OPTIONAL :: opt_i, opt_err
    LOGICAL, INTENT(OUT), OPTIONAL :: opt_l
    CLASS(*), POINTER :: keyObj, valObj
    CHARACTER(*), PARAMETER :: routine = modname//"get"

    keyObj => trimmed_key
    valObj => me%table%getEntry(keyObj)
    IF(ASSOCIATED(valObj)) THEN
      IF (PRESENT(opt_err)) opt_err = 0
      SELECT TYPE(valObj)
#ifdef __PGI
      TYPE IS(t_char_workaround)
        IF(.NOT.PRESENT(opt_c)) CALL finish(routine, "type mismatch '"//trimmed_key//"', is CHARACTER")
        opt_c = valObj%c
#else
      TYPE IS(CHARACTER(*))
        IF(.NOT.PRESENT(opt_c)) CALL finish(routine, "type mismatch '"//trimmed_key//"', is CHARACTER")
        opt_c = valObj
#endif
      TYPE IS(REAL(wp))
        IF(.NOT.PRESENT(opt_r)) CALL finish(routine, "type mismatch '"//trimmed_key//"', is REAL(wp)")
        opt_r = valObj
      TYPE IS(INTEGER)
        IF(.NOT.PRESENT(opt_i)) CALL finish(routine, "type mismatch '"//trimmed_key//"', is INTEGER")
        opt_i = INT(valObj, C_INT)
      TYPE IS(LOGICAL)
        IF(.NOT.PRESENT(opt_l)) CALL finish(routine, "type mismatch '"//trimmed_key//"', is LOGICAL")
        opt_l = valObj
      CLASS DEFAULT
        CALL finish(routine, "assertion failed")
      END SELECT
    ELSE
      IF (PRESENT(opt_err)) THEN
        opt_err = 1
      ELSE
        CALL key_value_store_output(me)
        CALL finish(routine, "key '"//trimmed_key//"' not found")
      END IF
    END IF
  END SUBROUTINE key_value_store_get_internal

  SUBROUTINE key_value_store_get_c(me, key, res, opt_err)
    CLASS(t_key_value_store), INTENT(IN) :: me
    CHARACTER(:), ALLOCATABLE, INTENT(OUT) :: res
    CHARACTER(*), INTENT(IN) :: key
    INTEGER, INTENT(OUT), OPTIONAL :: opt_err

    CALL key_value_store_get_internal(me, TRIM(key), opt_c = res, opt_err = opt_err)
  END SUBROUTINE key_value_store_get_c

  SUBROUTINE key_value_store_get_r(me, key, res, opt_err)
    CLASS(t_key_value_store), INTENT(IN) :: me
    REAL(wp), INTENT(OUT) :: res
    CHARACTER(*), INTENT(IN) :: key
    INTEGER, INTENT(OUT), OPTIONAL :: opt_err

    CALL key_value_store_get_internal(me, TRIM(key), opt_r = res, opt_err = opt_err)
  END SUBROUTINE key_value_store_get_r

  SUBROUTINE key_value_store_get_i(me, key, res, opt_err)
    CLASS(t_key_value_store), INTENT(IN) :: me
    CHARACTER(*), INTENT(IN) :: key
    INTEGER, INTENT(OUT) :: res
    INTEGER, INTENT(OUT), OPTIONAL :: opt_err

    CALL key_value_store_get_internal(me, TRIM(key), opt_i = res, opt_err = opt_err)
  END SUBROUTINE key_value_store_get_i

  SUBROUTINE key_value_store_get_l(me, key, res, opt_err)
    CLASS(t_key_value_store), INTENT(IN) :: me
    CHARACTER(*), INTENT(IN) :: key
    LOGICAL, INTENT(OUT) :: res
    INTEGER, INTENT(OUT), OPTIONAL :: opt_err

    CALL key_value_store_get_internal(me, TRIM(key), opt_l = res, opt_err = opt_err)
  END SUBROUTINE key_value_store_get_l

  SUBROUTINE key_value_store_input_cdi(me, vlistID)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: vlistId
    INTEGER :: natt, i, ierr, attType, attLen, namLen
    CHARACTER(CDI_MAX_NAME) :: attName
    REAL(C_DOUBLE) :: oneDouble(1)
    INTEGER(C_INT) :: oneInt(1)
    CHARACTER(:), ALLOCATABLE :: attTxt
    CHARACTER(*), PARAMETER :: routine = "input_cdi"

    CALL key_value_store_init(me, .true.)
    IF (vlistID == CDI_UNDEFID) CALL finish(routine, "invalid vlistId")
    ALLOCATE(CHARACTER(LEN=4096) :: attTxt, STAT=ierr)
    IF(ierr /= SUCCESS) CALL finish(routine, "memory allocation error")
    ierr = cdiinqNatts(vlistID, CDI_GLOBAL, natt)
    IF (ierr /= SUCCESS) CALL finish(routine, "getting attribute count failed")
    DO i = 0, natt - 1
      ierr = cdiInqAtt(vlistId, CDI_GLOBAL, i, attName, attType, attLen)
      IF (ierr /= SUCCESS) CALL finish(routine, "failed reading attribute name")
      namLen = LEN_TRIM(attName)
      SELECT CASE(attType)
      CASE(DATATYPE_FLT64)
        ierr = cdiInqAttFlt(vlistID, CDI_GLOBAL, attName(1:namLen), 1, oneDouble)
        IF (ierr /= SUCCESS) CALL finish(routine, "failed reading attribute '"//attName(1:namLen)//"'")
        CALL key_value_store_put_r(me, attName(1:namLen), REAL(oneDouble(1), wp))
      CASE(DATATYPE_INT32)
        ierr = cdiInqAttInt(vlistID, CDI_GLOBAL, attName(1:namLen), 1, oneInt)
        IF (ierr /= SUCCESS) CALL finish(routine, "failed reading attribute '"//attName(1:namLen)//"'")
        IF ('bool_' == attName(1:MIN(5,namLen))) THEN
          CALL key_value_store_put_l(me, attName(6:namLen), INT(oneInt(1)) .NE. 0)
        ELSE
          CALL key_value_store_put_i(me, attName(1:namLen), INT(oneInt(1)))
        END IF
      CASE(DATATYPE_TXT)
        IF (attLen .GT. LEN(attTxt)) DEALLOCATE(attTxt)
        IF (.NOT.ALLOCATED(attTxt)) THEN
          ALLOCATE(CHARACTER(LEN=attLen) :: attTxt, STAT=ierr)
          IF (ierr /= SUCCESS) CALL finish(routine, "memory allocation error")
        END IF
        ierr = cdiInqAttTxt(vlistID, CDI_GLOBAL, attName(1:namLen), attLen, attTxt(1:attLen))
        IF (ierr /= SUCCESS) CALL finish(routine, "failed reading attribute '"//attName(1:namLen)//"'")
        CALL key_value_store_put_c(me, attName(1:namLen), attTxt(1:attLen))
      END SELECT
    ENDDO
  END SUBROUTINE key_value_store_input_cdi

  SUBROUTINE key_value_store_input_pmsg(me, pmsg)
    CLASS(t_key_value_store), INTENT(INOUT) :: me
    TYPE(t_PackedMessage), INTENT(INOUT) :: pmsg
    INTEGER :: natt, i, attType, attLen, namLen, alnLen, alvLen, aInteger
    REAL(wp) :: aDouble
    LOGICAL :: aLogical, lcase_sensitive
    CHARACTER(:), ALLOCATABLE :: attTxt, attName

    alvLen = -1
    alnLen = -1
    CALL pmsg%unpack(lcase_sensitive)
    CALL key_value_store_init(me, lcase_sensitive)
    CALL pmsg%unpack(natt)
    DO i = 1, natt
      CALL pmsg%unpack(namLen)
      IF (namLen .NE. alnLen) THEN
        alnLen = namLen
        IF (ALLOCATED(attName)) DEALLOCATE(attName)
        ALLOCATE(CHARACTER(alnLen) :: attName)
      END IF
      CALL pmsg%unpack(attName)
      CALL pmsg%unpack(attType)
      SELECT CASE(attType)
      CASE(DATATYPE_FLT64)
        CALL pmsg%unpack(aDouble)
        CALL key_value_store_put_r(me, attName, aDouble)
      CASE(DATATYPE_INT32)
        CALL pmsg%unpack(aInteger)
        CALL key_value_store_put_i(me, attName, aInteger)
      CASE(DATATYPE_LONG)
        CALL pmsg%unpack(aLogical)
        CALL key_value_store_put_l(me, attName, aLogical)
      CASE(DATATYPE_TXT)
        CALL pmsg%unpack(attLen)
        IF (attLen .NE. alvLen) THEN
          alvLen = attLen
          IF (ALLOCATED(attTxt)) DEALLOCATE(attTxt)
          ALLOCATE(CHARACTER(alvLen) :: attTxt)
        END IF
        CALL pmsg%unpack(attTxt)
        CALL key_value_store_put_c(me, attName, attTxt)
      END SELECT
    ENDDO
  END SUBROUTINE key_value_store_input_pmsg

  SUBROUTINE key_value_store_output(me, vlistID, pmsg, archive, copy, inverse, label)
    CLASS(t_key_value_store), INTENT(IN) :: me
    INTEGER, INTENT(IN), OPTIONAL :: vlistId
    TYPE(t_PackedMessage), INTENT(INOUT), OPTIONAL :: pmsg
    TYPE(t_key_value_store), INTENT(INOUT), OPTIONAL :: archive, copy, inverse
    CHARACTER(*), INTENT(IN), OPTIONAL :: label
    CHARACTER(*), PARAMETER :: routine = modname//"output"
    TYPE(t_HashIterator) :: iterator
    CLASS(*), POINTER :: curKey, curVal
    CHARACTER(:), POINTER :: ccKey, ccVal
    LOGICAL :: selector(5)

    selector = [PRESENT(vlistID), PRESENT(pmsg), PRESENT(archive), &
      & PRESENT(copy), PRESENT(inverse)]
    IF (COUNT(selector(:)) .GT. 1) CALL finish(routine, "too many arguments")
    IF (selector(2)) THEN
      IF (.NOT.me%is_init) RETURN
      CALL pmsg%pack(me%lcase_sensitive)
      CALL pmsg%pack(me%table%getEntryCount())
    ELSE IF(selector(3)) THEN
      CALL key_value_store_init(archive, me%lcase_sensitive)
    ELSE IF(selector(4)) THEN
      CALL key_value_store_init(copy, me%lcase_sensitive)
    ELSE IF(selector(5)) THEN
      CALL key_value_store_init(inverse, me%lcase_sensitive)
    END IF
    CALL iterator%init(me%table)
    DO WHILE (iterator%nextEntry(curKey, curVal))
      ccKey => sel_char(curKey, routine, "key: invalid type")
      IF (selector(1)) THEN
        CALL output_cdi()
      ELSE IF(selector(2)) THEN
        CALL output_pmsg()
      ELSE IF(selector(3)) THEN
        IF (ccKey(1:4) /= 'nml_') CYCLE
        ccVal => sel_char(curVal, routine, "nml-archive copy: only character valued pairs!")
        CALL key_value_store_put_c(archive, ccKey, ccVal)
      ELSE IF(selector(4)) THEN
        ccVal => sel_char(curVal, routine, "copy: only character valued pairs!")
        CALL key_value_store_put_c(copy, ccKey, ccVal)
      ELSE IF(selector(5)) THEN
        ccVal => sel_char(curVal, routine, "inversion: only character valued pairs!")
        CALL key_value_store_put_c(inverse, ccVal, ccKey)
      ELSE
        CALL output_stdio()
      END IF
    END DO
  CONTAINS

    SUBROUTINE output_cdi()
      INTEGER :: ierr
  
      SELECT TYPE(curVal)
#ifdef __PGI
      TYPE IS(t_char_workaround)
        ierr = cdidefAttTxt(vlistId, CDI_GLOBAL, ccKey, LEN(curVal%c), curVal%c)
#else
      TYPE IS(CHARACTER(*))
        ierr = cdidefAttTxt(vlistId, CDI_GLOBAL, ccKey, LEN(curVal), curVal)
#endif
      TYPE IS(REAL(wp))
        ierr = cdidefAttFlt(vlistId, CDI_GLOBAL, ccKey, DATATYPE_FLT64, 1, [REAL(curVal, C_DOUBLE)])
      TYPE IS(INTEGER)
        ierr = cdidefAttInt(vlistId, CDI_GLOBAL, ccKey, DATATYPE_INT32, 1, [INT(curVal, C_INT)])
      TYPE IS(LOGICAL)
        ierr = cdidefAttInt(vlistId, CDI_GLOBAL, 'bool_'//ccKey, DATATYPE_INT32, 1, [INT(MERGE(1,0,curVal), C_INT)])
      CLASS DEFAULT
        CALL finish(routine, "val: invalid type")
      END SELECT
      IF(ierr /= SUCCESS) CALL finish(routine, "error writing attribute")
    END SUBROUTINE output_cdi
  
    SUBROUTINE output_pmsg()
      INTEGER :: clen
  
      clen = LEN_TRIM(ccKey)
      CALL pmsg%pack(clen)
      CALL pmsg%pack(ccKey(1:clen))
      SELECT TYPE(curVal)
#ifdef __PGI
      TYPE IS(t_char_workaround)
        ccVal => curVal%c
#else
      TYPE IS(CHARACTER(*))
        ccVal => curVal
#endif
        CALL pmsg%pack(DATATYPE_TXT)
        clen = LEN_TRIM(ccVal)
        CALL pmsg%pack(clen)
        CALL pmsg%pack(ccVal(1:clen))
      TYPE IS(REAL(wp))
        CALL pmsg%pack(DATATYPE_FLT64)
        CALL pmsg%pack(curVal)
      TYPE IS(INTEGER)
        CALL pmsg%pack(DATATYPE_INT32)
        CALL pmsg%pack(curVal)
      TYPE IS(LOGICAL)
        CALL pmsg%pack(DATATYPE_LONG)
        CALL pmsg%pack(curVal)
      CLASS DEFAULT
        CALL finish(routine, "val: invalid type")
      END SELECT
    END SUBROUTINE output_pmsg

    SUBROUTINE output_stdio()
  
      message_text = ''
      SELECT TYPE(curVal)
#ifdef __PGI
      TYPE IS(t_char_workaround)
        WRITE(message_text, "(5a)") "key = >", ccKey, "< val = >", curVal%c, "<"
#else
      TYPE IS(CHARACTER(*))
        WRITE(message_text, "(5a)") "key = >", ccKey, "< val = >", curVal, "<"
#endif
      TYPE IS(REAL(wp))
        WRITE(message_text, "(3a,e12.5,a)") "key = >", ccKey, "< val = >", curVal, "<"
      TYPE IS(INTEGER)
        WRITE(message_text, "(3a,i6,a)") "key = >", ccKey, "< val = >", curVal, "<"
      TYPE IS(LOGICAL)
        WRITE(message_text, "(3a,l1,a)") "key = >", ccKey, "< val = >", curVal, "<"
      CLASS DEFAULT
        CALL finish(routine, "val: invalid type")
      END SELECT
      IF (PRESENT(label)) THEN
        CALL message(routine, message_text)
      ELSE
        PRINT "(a)", TRIM(message_text)
      END IF
    END SUBROUTINE output_stdio
  END SUBROUTINE key_value_store_output

  INTEGER FUNCTION key_value_store_getEntryCount(me) RESULT(n)
    CLASS(t_key_value_store), INTENT(IN) :: me

    n = me%table%getEntryCount()
  END FUNCTION key_value_store_getEntryCount

  SUBROUTINE key_value_store_destruct(me)
    CLASS(t_key_value_store), INTENT(INOUT) :: me

    IF (ASSOCIATED(me%table)) THEN
      CALL me%table%destruct()
      DEALLOCATE(me%table)
    END IF
    me%lcase_sensitive = .true.
    me%is_init = .false.
  END SUBROUTINE key_value_store_destruct

END MODULE mo_key_value_store
