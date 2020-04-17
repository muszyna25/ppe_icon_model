!> Module for associating pairs of strings,
!! e.g. for translating between variable names.
!!
!! @author F. Prill, DWD;
!! essential parts from mo_storage: Daniel Rieger, KIT (2016-12-13)
!!
!! @note When loading key-value pairs from an external text file, the
!!       current implementation is restricted to strings that do not
!!       contain spaces.
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2013-01-21)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! -----------------------------------------------------------------------------------
MODULE mo_dictionary

#ifdef __ICON__
  USE mo_exception,      ONLY: finish, message_text
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_io_units,       ONLY: find_next_free_unit
#else
  USE mo_utilities,      ONLY: finish, message_text, SUCCESS, &
    &                          find_next_free_unit
#endif
  USE mo_hash_table,     ONLY: t_HashTable_Base, hashTable_make, t_stringVal,                 &
    &                          stringVal_len, storage_hashKey_DJB_cs, storage_hashKey_DJB_ci, &
    &                          storage_equalKeysFunction_cs, storage_equalKeysFunction_ci,    &
    &                          t_const_hashIterator
  USE mo_fortran_tools,  ONLY: t_Destructible

  IMPLICIT NONE

  !--------------------------------------------------------------------------
  ! definition of constants:

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_dictionary'

  INTEGER, PARAMETER :: DICT_MAX_STRLEN = stringVal_len            !< maximum string length


  !--------------------------------------------------------------------------
  ! type definition:

  ! dictionary data type
  TYPE t_dictionary
    !> Flag. If .TRUE. key strings are handled case-sensitively
    LOGICAL :: lcase_sensitive
    !> key-value hashtable (and its inverse):
    CLASS(t_HashTable_Base), POINTER :: hashtable, hashtable_inv

  CONTAINS
    PROCEDURE :: init       => dict_init
    PROCEDURE :: finalize   => dict_finalize
    PROCEDURE :: set        => dict_set
    PROCEDURE :: get        => dict_get
    PROCEDURE :: get_size   => dict_get_size
    PROCEDURE :: loadfile   => dict_loadfile
    PROCEDURE :: copy       => dict_copy
    PROCEDURE :: to_array   => dict_to_array
    PROCEDURE :: from_array => dict_from_array
  END TYPE t_dictionary


  ! public interface definition
  PRIVATE
  ! data
  PUBLIC :: t_dictionary
  PUBLIC :: DICT_MAX_STRLEN


CONTAINS

  !--------------------------------------------------------------------------
  !> Initializes dictionary data structure.
  !
  SUBROUTINE dict_init(dict, lcase_sensitive)
    CLASS(t_dictionary), INTENT(INOUT) :: dict            !< dictionary data structure
    LOGICAL,             INTENT(IN)    :: lcase_sensitive !< Flag. If .TRUE. keys are handled case-sensitively
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::dict_init'

    dict%lcase_sensitive = lcase_sensitive
    IF (lcase_sensitive) THEN
      dict%hashtable     => hashTable_make(storage_hashKey_DJB_cs, storage_equalKeysFunction_cs)
      dict%hashtable_inv => hashTable_make(storage_hashKey_DJB_cs, storage_equalKeysFunction_cs)
    ELSE
      dict%hashtable     => hashTable_make(storage_hashKey_DJB_ci, storage_equalKeysFunction_ci)
      dict%hashtable_inv => hashTable_make(storage_hashKey_DJB_ci, storage_equalKeysFunction_ci)
    ENDIF
  END SUBROUTINE dict_init


  !--------------------------------------------------------------------------
  !> Destroys dictionary data structure and frees allocated memory.
  !
  SUBROUTINE dict_finalize(dict)
    CLASS(t_dictionary), INTENT(INOUT) :: dict !< dictionary data structure
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::dict_finalize'

    IF (ASSOCIATED(dict%hashtable)) THEN
      CALL dict%hashtable%destruct()
      DEALLOCATE(dict%hashtable)
    END IF
    IF (ASSOCIATED(dict%hashtable_inv)) THEN
      CALL dict%hashtable_inv%destruct()
      DEALLOCATE(dict%hashtable_inv)
    END IF
  END SUBROUTINE dict_finalize


  !--------------------------------------------------------------------------
  !> Deep-copy of dictionary data structure.
  !
  SUBROUTINE dict_copy(src_dict, dst_dict)
    CLASS(t_dictionary), INTENT(IN)    :: src_dict !< source dictionary data structure
    TYPE(t_dictionary),  INTENT(INOUT) :: dst_dict !< destination dictionary data structure
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::dict_copy'

    ! "const" iterator for traversing the hash table entries: Note
    ! that the dictionary is INTENT(IN) and we want to make clear that
    ! we do not change the contents of the hash table.
    TYPE(t_const_hashIterator) :: it

    CLASS(t_Destructible), ALLOCATABLE  :: p_key, p_value

    CALL dst_dict%init(src_dict%lcase_sensitive)
    CALL it%init(src_dict%hashtable)
    DO
      IF (.NOT. it%nextEntry(p_key, p_value))  EXIT

      SELECT TYPE (p_key)
      TYPE IS(t_stringVal)
        IF (ALLOCATED(p_value)) THEN
          SELECT TYPE (p_value)
          TYPE IS(t_stringVal)
            CALL dst_dict%set(p_key%stringVal, p_value%stringVal)
          END SELECT
        ENDIF
      CLASS DEFAULT
        CALL finish(routine, "Wrong entry type.")
      END SELECT

      IF (ALLOCATED(p_key))    DEALLOCATE(p_key)
      IF (ALLOCATED(p_value))  DEALLOCATE(p_value)
    END DO
  END SUBROUTINE dict_copy


  !--------------------------------------------------------------------------
  !> Copy dictionary data to an array.
  !
  SUBROUTINE dict_to_array(dict, array)
    !> source dictionary data structure
    CLASS(t_dictionary), INTENT(IN) :: dict
    !> dictionary data, dims: (2,nmax_entries)
    CHARACTER(LEN=DICT_MAX_STRLEN), ALLOCATABLE, INTENT(INOUT) :: array(:,:)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::dict_to_array'
    INTEGER :: i, ierror

    ! "const" iterator for traversing the hash table entries: Note
    ! that the dictionary is INTENT(IN) and we want to make clear that
    ! we do not change the contents of the hash table.
    TYPE(t_const_hashIterator) :: it

    CLASS(t_Destructible), ALLOCATABLE  :: p_key, p_value

    ALLOCATE(array(2,dict%hashtable%nentries()), STAT=ierror)
    IF (ierror /= SUCCESS) CALL finish(routine, "memory allocation failure")

    i = 0
    CALL it%init(dict%hashtable)
    DO
      IF (.NOT. it%nextEntry(p_key, p_value))  EXIT

      SELECT TYPE (p_key)
      TYPE IS(t_stringVal)
        IF (ALLOCATED(p_value)) THEN
          SELECT TYPE (p_value)
          TYPE IS(t_stringVal)
            i = i + 1
            array(1,i) = TRIM(p_key%stringVal)
            array(2,i) = TRIM(p_value%stringVal)
          END SELECT
        ENDIF
      END SELECT

      IF (ALLOCATED(p_key))    DEALLOCATE(p_key)
      IF (ALLOCATED(p_value))  DEALLOCATE(p_value)
    END DO
    ! consistency check
    IF (i /= SIZE(array,2)) THEN
      CALL finish(routine, "Internal error!")
    END IF
  END SUBROUTINE dict_to_array


  !--------------------------------------------------------------------------
  !> Initialize a dictionary with data from an array.
  !
  SUBROUTINE dict_from_array(dict, array, lcase_sensitive)
    !> target dictionary data structure
    CLASS(t_dictionary),            INTENT(INOUT) :: dict
    !> dictionary data, dims: (2,nmax_entries)
    CHARACTER(LEN=DICT_MAX_STRLEN), INTENT(IN) :: array(:,:)
    !> Flag. If .TRUE. keys are handled case-sensitively
    LOGICAL,                        INTENT(IN)    :: lcase_sensitive

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::dict_from_array'
    INTEGER :: i

    CALL dict%init(lcase_sensitive)
    DO i=1,SIZE(array,2)
      CALL dict%set(array(1,i), array(2,i))
    END DO
  END SUBROUTINE dict_from_array


  !--------------------------------------------------------------------------
  !> Auxiliary function: Add string to hash table.
  SUBROUTINE put_string(hashtable, key, value)
    CLASS(t_HashTable_Base),INTENT(inout),POINTER :: hashtable
    CHARACTER(LEN=*),INTENT(in)           :: key, value
    ! Local
    CLASS(t_Destructible), POINTER  :: p_key, p_value

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
    CALL hashtable%setEntry(p_key, p_value)
  END SUBROUTINE put_string


  !--------------------------------------------------------------------------
  !> Insert a new key into the given dictionary. Overwrites any
  !  existing key-value pair.
  !
  SUBROUTINE dict_set(dict, key, val)
    CLASS(t_dictionary), INTENT(INOUT) :: dict !< dictionary data structure
    CHARACTER(LEN=*),   INTENT(IN)     :: key  !< new search key
    CHARACTER(LEN=*),   INTENT(IN)     :: val  !< new value associated with key
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::dict_set'

    IF ((LEN(TRIM(key)) > DICT_MAX_STRLEN)  .OR.  &
      & (LEN(TRIM(val)) > DICT_MAX_STRLEN)) THEN
      CALL finish(routine, "Key/value pair does not fit into dictionary.")
    END IF
    CALL put_string(dict%hashtable, key, val)
    CALL put_string(dict%hashtable_inv, val, key)

    ! consistency check
    IF (dict%hashtable%nentries() /= dict%hashtable_inv%nentries()) THEN
      CALL finish(routine, "Dictionary is not invertible one-to-one!")
    END IF
  END SUBROUTINE dict_set


  !--------------------------------------------------------------------------
  !> Auxiliary function: Get an entry from string hash table.
  !
  ! Note: Implemented as subroutine rather than function as the
  !       interface "get" is not distinguishable otherwise.
  !
  SUBROUTINE get_string(hashtable, key, value, ierror)
    CLASS(t_HashTable_Base),INTENT(in) :: hashtable
    CHARACTER(LEN=*),INTENT(in)        :: key
    CHARACTER(LEN=*),INTENT(out)       :: value
    INTEGER,OPTIONAL,INTENT(out)       :: ierror
    ! Local
    CHARACTER(LEN=*), PARAMETER        :: routine = modname//":get_string"
    CLASS(t_Destructible),POINTER      :: p_key, p_value

    ALLOCATE(t_stringVal :: p_key)
    IF (PRESENT(ierror)) ierror = SUCCESS

    SELECT TYPE(p_key)
    TYPE IS(t_stringVal)
      p_key%stringVal = TRIM(key)
    END SELECT

    p_value => hashtable%getEntry(p_key)

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


  !--------------------------------------------------------------------------
  !> Retrieve a value from the given dictionary which is assciated
  !  with @p key.
  !
  !  @param[in] linverse   optional flag: If .TRUE., the dictionary is queried in inverse order.
  !
  !  @return If the key is not found, this function returns the
  !          "default" value. If this optional argument has not been
  !          provided, this function throws an error message.
  !
  FUNCTION dict_get(dict, key, default, linverse)
    CHARACTER(LEN=DICT_MAX_STRLEN)           :: dict_get !< return value.
    CLASS(t_dictionary),INTENT(IN)           :: dict     !< dictionary data structure
    CHARACTER(LEN=*),   INTENT(IN)           :: key      !< new search key
    CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: default  !< default value
    LOGICAL,            INTENT(IN), OPTIONAL :: linverse !< Flag. If .TRUE., the dictionary is queried in inverse order.
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::dict_get'
    INTEGER :: ierror
    LOGICAL :: flag_linverse

    flag_linverse = .FALSE.
    IF (PRESENT(linverse))  flag_linverse = linverse

    ierror = 0
    IF (flag_linverse) THEN
      CALL get_string(dict%hashtable_inv, key, dict_get, ierror)
    ELSE
      CALL get_string(dict%hashtable, key, dict_get, ierror)
    END IF

    IF (ierror /= 0) THEN
      IF (PRESENT(default)) THEN
        dict_get = default
      ELSE
        WRITE(message_text,'(a,a,a)') 'Requested dictionary key ', TRIM(key), ' not found!'
        CALL finish(routine, message_text)
      END IF
    END IF
  END FUNCTION dict_get


  !--------------------------------------------------------------------------
  !> Get the number of entries within the dictionary.
  !
  INTEGER FUNCTION dict_get_size(dict) RESULT(resultVar)
    CLASS(t_dictionary), INTENT(IN) :: dict
    resultVar = dict%hashtable%nentries()
  END FUNCTION dict_get_size
  

  !--------------------------------------------------------------------------
  !> Load the contents of a text file into the given dictionary data
  !  structure.
  !
  !  @note The file format is assumed to be as follows: Each line of
  !        the text file contains a key-value pair, where both strings
  !        are separated by one or more spaces. Comment lines
  !        (beginning with "#") are ignored.
  !
  SUBROUTINE dict_loadfile(dict, filename, linverse)
    CLASS(t_dictionary), INTENT(INOUT) :: dict      !< dictionary data structure
    CHARACTER(LEN=*),   INTENT(IN)    :: filename  !< text file name
    LOGICAL, OPTIONAL,  INTENT(IN)    :: linverse  !< Flag. If .TRUE., dictionary columns are read in inverse order.
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::dict_loadfile'
    INTEGER                        :: iunit, ist
    CHARACTER(LEN=256)             :: line
    CHARACTER(LEN=DICT_MAX_STRLEN) :: key, val
    LOGICAL                        :: valid, lread_inverse

    lread_inverse = .FALSE.
    IF (PRESENT(linverse)) THEN
      lread_inverse = linverse
    END IF

    iunit = find_next_free_unit(10,99)
    OPEN (unit=iunit,file=filename,access='SEQUENTIAL', &
      &  form='FORMATTED', action='READ', status='OLD', IOSTAT=ist)
    IF(ist/=0) CALL finish (routine, 'Open of dictionary file '//TRIM(filename)//' failed')
    DO
      line = ' '
      READ(iunit,'(a)',iostat=ist) line
      IF(ist < 0) EXIT ! No more lines
      IF(ist > 0) CALL finish(routine, 'Read error in dictionary file '//TRIM(filename))
      CALL parse_line(line,key,val,valid)
      IF(valid) THEN
        IF (.NOT. lread_inverse) THEN
          CALL  dict%set(TRIM(ADJUSTL(key)), TRIM(ADJUSTL(val)))
        ELSE
          CALL  dict%set(TRIM(ADJUSTL(val)), TRIM(ADJUSTL(key)))
        END IF
      END IF
    ENDDO
    ! close file
    CLOSE(unit=iunit)
  END SUBROUTINE dict_loadfile


  !--------------------------------------------------------------------------
  !> Utility function: Parses a line from map_file and returns key
  !  and value part (first two nonblank words separated by blanks or
  !  tabs).
  !
  !  Initial implementation: R. Johanni, 2011
  !
  SUBROUTINE parse_line(line, key, val, valid)
    CHARACTER(LEN=*), INTENT(IN)  :: line
    CHARACTER(LEN=*), INTENT(OUT) :: key, val
    LOGICAL,          INTENT(OUT) :: valid
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::parse_line'
    INTEGER :: ipos1, ipos2

    valid = .FALSE.
    key = ' '
    val = ' '
    ! Search first nonblank character
    DO ipos1 = 1, LEN(line)
      IF(.NOT.isblank(line(ipos1:ipos1))) EXIT
    ENDDO
    IF(ipos1 > LEN(line)) RETURN        ! completely empty line
    IF(line(ipos1:ipos1) == '#') RETURN ! comment line
    ! Search end of key
    DO ipos2 = ipos1, LEN(line)
      IF(isblank(line(ipos2:ipos2))) EXIT
    ENDDO
    IF ((ipos2-ipos1) > DICT_MAX_STRLEN) THEN
      CALL finish(routine, "Key does not fit into dictionary.")
    END IF
    key = line(ipos1:ipos2-1)
    ! Search next nonblank character
    DO ipos1 = ipos2, LEN(line)
      IF(.NOT.isblank(line(ipos1:ipos1))) EXIT
    ENDDO
    IF(ipos1 > LEN(line)) THEN ! line contains no value part
      CALL finish(routine, "Illegal line in name_map.")
      RETURN
    ENDIF
    ! Search end of val
    DO ipos2 = ipos1, LEN(line)
      IF(isblank(line(ipos2:ipos2))) EXIT
    ENDDO
    IF ((ipos2-ipos1) > DICT_MAX_STRLEN) THEN
      CALL finish(routine, "Value does not fit into dictionary.")
    END IF
    val = line(ipos1:ipos2-1)
    valid = .TRUE.

    CONTAINS

    ! Fortran equivalent to C isblank function
    LOGICAL FUNCTION isblank(c)
      CHARACTER(LEN=1) :: c
      IF(c==' ' .OR. ICHAR(c)==9 .OR. ICHAR(c)==0) THEN
        ! Blank or Tab
        isblank = .TRUE.
      ELSE
        isblank = .FALSE.
      ENDIF
    END FUNCTION isblank
  END SUBROUTINE parse_line

END MODULE mo_dictionary
