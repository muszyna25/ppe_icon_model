!>
!! Module for associating pairs of strings,
!! e.g. for translating between variable names.
!!
!! @author F. Prill, DWD
!!
!! @note This implementation is aimed at *small* dictionaries,
!!       e.g. translation between variable names. The algorithms
!!       contained in this module are therefore simple and will
!!       not scale with large data bases.
!!
!! @note When loading key-value pairs from an external text file, the
!!       current implementation is restricted to string that do not
!!       contain spaces.
!!
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2013-01-21)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!! -----------------------------------------------------------------------------------
MODULE mo_dictionary

  USE mo_exception,      ONLY: finish
  USE mo_util_string,    ONLY: tolower    
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_io_units,       ONLY: find_next_free_unit

  IMPLICIT NONE

  !--------------------------------------------------------------------------
  ! definition of constants:

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_grid_icon')  

  INTEGER, PARAMETER :: DICT_MAX_STRLEN = 100    !< maximum string length
  INTEGER, PARAMETER :: NINITIAL        = 100    !< initial dictionary size
  INTEGER, PARAMETER :: DICT_UNDEFID    =  -1    !< return value: undefined/not found


  !--------------------------------------------------------------------------
  ! type definition:

  TYPE t_dictionary
    INTEGER :: nentries         !< current no. of entries
    INTEGER :: nmax_entries     !< current array size
    LOGICAL :: lcase_sensitive  !< Flag. If .TRUE. key strings are handled case-sensitively
  
    CHARACTER(LEN=DICT_MAX_STRLEN), ALLOCATABLE :: array(:,:) !< dictionary data, dims: (2,nmax_entries)
  END TYPE t_dictionary

  ! public interface definition
  PRIVATE
  ! functions and subroutines
  PUBLIC :: dict_init
  PUBLIC :: dict_finalize
  PUBLIC :: dict_set
  PUBLIC :: dict_get
  PUBLIC :: dict_loadfile
  ! data
  PUBLIC :: t_dictionary
  PUBLIC :: DICT_MAX_STRLEN
 

CONTAINS

  !--------------------------------------------------------------------------
  !> Initializes dictionary data structure, pre-allocates memory.
  !  If the @p new_size argument has not been provided, the dictionary
  !  is pre-allocated with length NINITIAL.
  !
  SUBROUTINE dict_init(dict, lcase_sensitive, new_size)
    TYPE(t_dictionary), INTENT(INOUT) :: dict            !< dictionary data structure
    LOGICAL,            INTENT(IN)    :: lcase_sensitive !< Flag. If .TRUE. key strings are handled case-sensitively
    INTEGER,            INTENT(IN), OPTIONAL :: new_size !< New dictionary size
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: &
      &  routine    = TRIM(TRIM(modname)//'::dict_init')
    INTEGER :: new_nmax_entries
    
    new_nmax_entries = NINITIAL
    IF (PRESENT(new_size)) new_nmax_entries = new_size
    dict%nentries        = 0
    dict%lcase_sensitive = lcase_sensitive
    CALL dict_resize(dict, new_nmax_entries)
  END SUBROUTINE dict_init


  !--------------------------------------------------------------------------
  !> Destroys dictionary data structure and frees allocated memory.
  !
  SUBROUTINE dict_finalize(dict)
    TYPE(t_dictionary), INTENT(INOUT) :: dict !< dictionary data structure
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: &
      &  routine    = TRIM(TRIM(modname)//'::dict_finalize')
    INTEGER :: ierrstat

    IF (ALLOCATED(dict%array)) THEN
      DEALLOCATE(dict%array, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
    dict%nentries     = 0
    dict%nmax_entries = 0
  END SUBROUTINE dict_finalize


  !--------------------------------------------------------------------------
  !> Copies dictionary data structure to a new size. If the @p
  !  new_size argument has not been provided, the dictionary is copied
  !  to a data structure of twice the current size.
  !
  SUBROUTINE dict_resize(dict, new_size)
    TYPE(t_dictionary), INTENT(INOUT) :: dict !< dictionary data structure
    INTEGER,            INTENT(IN), OPTIONAL :: new_size !< New dictionary size.
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: &
      &  routine    = TRIM(TRIM(modname)//'::dict_resize')
    CHARACTER(LEN=DICT_MAX_STRLEN), ALLOCATABLE :: tmp(:,:)
    INTEGER :: new_nmax_entries, ierrstat
    
    new_nmax_entries = 2*dict%nmax_entries
    IF (PRESENT(new_size)) new_nmax_entries = new_size
    if (new_nmax_entries < dict%nentries) &
      &  CALL finish(routine, "Dictionary cannot be resized.")

    ! allocate temporary storage
    IF (ALLOCATED(dict%array)) THEN
      ALLOCATE(tmp(2,dict%nentries), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      ! copy existing data and resize
      tmp(:,1:dict%nentries) = dict%array(:,1:dict%nentries)
      DEALLOCATE(dict%array, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
    ALLOCATE(dict%array(2,new_nmax_entries), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    IF (ALLOCATED(tmp)) THEN
      dict%array(:,1:dict%nentries) = tmp(:,1:dict%nentries)
    END IF
    dict%nmax_entries = new_nmax_entries
    ! clean up
    IF (ALLOCATED(tmp)) THEN
      DEALLOCATE(tmp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
  END SUBROUTINE dict_resize


  !--------------------------------------------------------------------------
  !> Insert a new key into the given dictionary. Overwrites any
  !  existing key-value pair.
  !
  SUBROUTINE dict_set(dict, key, val)
    TYPE(t_dictionary), INTENT(INOUT) :: dict !< dictionary data structure
    CHARACTER(LEN=*),   INTENT(IN)    :: key  !< new search key
    CHARACTER(LEN=*),   INTENT(IN)    :: val  !< new value associated with key
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: &
      &  routine    = TRIM(TRIM(modname)//'::dict_set')
    INTEGER :: idx

    IF ((LEN(TRIM(key)) > DICT_MAX_STRLEN)  .OR.  &
      & (LEN(TRIM(val)) > DICT_MAX_STRLEN)) THEN
      CALL finish(routine, "Key/value pair does not fit into dictionary.")
    END IF
    idx = dict_find(dict, key)
    IF (idx == DICT_UNDEFID) THEN
      ! resize dictionary, if necessary
      idx = dict%nentries + 1
      IF (idx > dict%nmax_entries) CALL dict_resize(dict)
      dict%nentries = dict%nentries + 1
    END IF
    ! insert new key-value pair
    dict%array(1,idx) = TRIM(key)
    dict%array(2,idx) = TRIM(val)
  END SUBROUTINE dict_set


  !--------------------------------------------------------------------------
  !> Retrieve a value from the given dictionary which is assciated
  !  with @p key.
  !
  !  @return If the key is not found, this function returns the
  !          "default" value. If this optional argument has not been
  !          provided, this function throws an error message.
  !
  FUNCTION dict_get(dict, key, default)
    CHARACTER(LEN=DICT_MAX_STRLEN) dict_get !< return value.
    TYPE(t_dictionary), INTENT(IN)    :: dict !< dictionary data structure
    CHARACTER(LEN=*),   INTENT(IN)    :: key  !< new search key
    CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: default !< default value
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: &
      &  routine    = TRIM(TRIM(modname)//'::dict_get')
    INTEGER :: idx

    idx = dict_find(dict, key)
    IF (idx == DICT_UNDEFID) THEN
      IF (PRESENT(default)) THEN
        dict_get = TRIM(default)
      ELSE
        CALL finish(routine, "Requested dictionary key not found!")
      END IF
    ELSE
      dict_get = TRIM(dict%array(2,idx))
    END IF
  END FUNCTION dict_get


  !--------------------------------------------------------------------------
  !> Load the contents of a text file into the given dictionary data
  !  structure.
  ! 
  !  @note The file format is assumed to be as follows: Each line of
  !        the text file contains a key-value pair, where both strings
  !        are separated by one or more spaces. Comment lines
  !        (beginning with "#") are ignored.
  !
  SUBROUTINE dict_loadfile(dict, filename)
    TYPE(t_dictionary), INTENT(INOUT) :: dict      !< dictionary data structure
    CHARACTER(LEN=*),   INTENT(IN)    :: filename  !< text file name
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: &
      &  routine    = TRIM(TRIM(modname)//'::dict_loadfile')
    INTEGER                        :: iunit, ist
    CHARACTER(LEN=256)             :: line
    CHARACTER(LEN=DICT_MAX_STRLEN) :: key, val
    LOGICAL                        :: valid

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
      IF(valid) CALL  dict_set(dict, TRIM(key), TRIM(val))
    ENDDO
    ! close file
    CLOSE(unit=iunit)
  END SUBROUTINE dict_loadfile


  !--------------------------------------------------------------------------
  !> Find the first key-value pair in the given dictionary which is
  !  assciated with @p key.
  ! 
  !  @return Array position index or, if the key has not been found,
  !          DICT_UNDEFID.
  !
  FUNCTION dict_find(dict, key)
    INTEGER :: dict_find
    TYPE(t_dictionary), INTENT(IN) :: dict      !< dictionary data structure
    CHARACTER(LEN=*),   INTENT(IN) :: key
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: &
      &  routine    = TRIM(TRIM(modname)//'::dict_find')
    INTEGER :: idx

    ! simply loop over all dictionary entries
    dict_find = DICT_UNDEFID
    LOOP_IDX : DO idx=1,dict%nentries
      IF (.NOT. dict%lcase_sensitive) THEN
        IF (TRIM(tolower(TRIM(key))) == TRIM(tolower(TRIM(dict%array(1,idx))))) THEN
          dict_find = idx
          EXIT LOOP_IDX
        END IF
      ELSE
        IF (TRIM(key) ==TRIM(dict%array(1,idx))) THEN
          dict_find = idx
          EXIT LOOP_IDX
        END IF
      END IF
    END DO LOOP_IDX
  END FUNCTION dict_find


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
    CHARACTER(LEN=*), PARAMETER :: &
      &  routine    = TRIM(TRIM(modname)//'::parse_line')
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
      IF(c==' ' .OR. ICHAR(c)==9) THEN
        ! Blank or Tab
        isblank = .TRUE.
      ELSE
        isblank = .FALSE.
      ENDIF
    END FUNCTION isblank
  END SUBROUTINE parse_line

END MODULE mo_dictionary
