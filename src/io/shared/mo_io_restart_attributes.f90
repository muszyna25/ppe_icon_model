!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_io_restart_attributes
  !
  USE ISO_C_BINDING,            ONLY: C_DOUBLE, C_INT
  USE mo_kind,                  ONLY: wp, i8
  USE mo_exception,             ONLY: finish
  USE mo_mpi,                   ONLY: p_bcast
  USE mo_cdi,                   ONLY: DATATYPE_FLT64, DATATYPE_INT32, DATATYPE_TXT, CDI_GLOBAL, vlistInqNatts, vlistInqAtt, &
      &                               vlistInqAttFlt, vlistInqAttInt, vlistInqAttTxt
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: set_restart_attribute
  PUBLIC :: get_restart_attribute
  PUBLIC :: read_and_bcast_attributes
  PUBLIC :: delete_attributes
  !
  PUBLIC :: restart_attributes_count_text
  PUBLIC :: restart_attributes_count_real
  PUBLIC :: restart_attributes_count_int
  PUBLIC :: restart_attributes_count_bool
  !
  !------------------------------------------------------------------------------------------------
  !
  TYPE t_att_text
    CHARACTER(len= 64) :: name
    CHARACTER(len=256) :: val
  END TYPE t_att_text
  !
  TYPE t_att_real
    CHARACTER(len=64) :: name
    REAL(wp)          :: val
  END TYPE t_att_real
  TYPE t_att_int
    CHARACTER(len=64) :: name
    INTEGER           :: val
  END TYPE t_att_int
  TYPE t_att_bool
    CHARACTER(len=64) :: name
    LOGICAL           :: val
    INTEGER           :: store
  END TYPE t_att_bool
  !
  INTEGER, PARAMETER :: nmax_atts = 1024
  INTEGER, SAVE :: natts_text = 0
  INTEGER, SAVE :: natts_real = 0 
  INTEGER, SAVE :: natts_int  = 0 
  INTEGER, SAVE :: natts_bool = 0 
  TYPE(t_att_text), ALLOCATABLE :: restart_attributes_text(:)
  TYPE(t_att_real), ALLOCATABLE :: restart_attributes_real(:)
  TYPE(t_att_int),  ALLOCATABLE :: restart_attributes_int(:)
  TYPE(t_att_bool), ALLOCATABLE :: restart_attributes_bool(:)
  !
  !------------------------------------------------------------------------------------------------
  !
  INTERFACE set_restart_attribute
    MODULE PROCEDURE set_restart_attribute_text
    MODULE PROCEDURE set_restart_attribute_real
    MODULE PROCEDURE set_restart_attribute_int
    MODULE PROCEDURE set_restart_attribute_int8
    MODULE PROCEDURE set_restart_attribute_bool
  END INTERFACE set_restart_attribute
  !
  INTERFACE get_restart_attribute
    MODULE PROCEDURE get_restart_att_text_by_name
    MODULE PROCEDURE get_restart_att_text_by_index
    MODULE PROCEDURE get_restart_att_real_by_name
    MODULE PROCEDURE get_restart_att_real_by_index
    MODULE PROCEDURE get_restart_att_int_by_name
    MODULE PROCEDURE get_restart_att_int8_by_name
    MODULE PROCEDURE get_restart_att_int_by_index
    MODULE PROCEDURE get_restart_att_bool_by_name
    MODULE PROCEDURE get_restart_att_bool_by_index
  END INTERFACE get_restart_attribute
  !
CONTAINS
  !
  !------------------------------------------------------------------------------------------------
  !
  FUNCTION restart_attributes_count_text() RESULT(natts)
    INTEGER :: natts
    natts = natts_text
  END FUNCTION restart_attributes_count_text
  !
  FUNCTION restart_attributes_count_real() RESULT(natts)
    INTEGER :: natts
    natts = natts_real
  END FUNCTION restart_attributes_count_real
  !
  FUNCTION restart_attributes_count_int() RESULT(natts)
    INTEGER :: natts
    natts = natts_int
  END FUNCTION restart_attributes_count_int
  !
  FUNCTION restart_attributes_count_bool() RESULT(natts)
    INTEGER :: natts
    natts = natts_bool
  END FUNCTION restart_attributes_count_bool
  !
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE get_restart_att_text_by_name(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in)  :: attribute_name
    CHARACTER(len=*), INTENT(out) :: attribute_value
    INTEGER :: i
    DO i = 1, natts_text
      IF (TRIM(attribute_name) == TRIM(restart_attributes_text(i)%name)) THEN
        attribute_value = TRIM(restart_attributes_text(i)%val)
        RETURN
      ENDIF
    ENDDO
    CALL finish('','Attribute '//TRIM(attribute_name)//' not found.')
  END SUBROUTINE get_restart_att_text_by_name
  !
  SUBROUTINE get_restart_att_text_by_index(idx, attribute_name, attribute_value)
    INTEGER,          INTENT(in)  :: idx
    CHARACTER(len=*), INTENT(out) :: attribute_name
    CHARACTER(len=*), INTENT(out) :: attribute_value
    attribute_name  = TRIM(restart_attributes_text(idx)%name)
    attribute_value = TRIM(restart_attributes_text(idx)%val)
  END SUBROUTINE get_restart_att_text_by_index
  !
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE get_restart_att_real_by_name(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in)  :: attribute_name
    REAL(wp),         INTENT(out) :: attribute_value
    INTEGER :: i
    DO i = 1, natts_real
      IF (TRIM(attribute_name) == TRIM(restart_attributes_real(i)%name)) THEN
        attribute_value = restart_attributes_real(i)%val
        RETURN
      ENDIF
    ENDDO
    CALL finish('','Attribute '//TRIM(attribute_name)//' not found.')
  END SUBROUTINE get_restart_att_real_by_name
  !
  SUBROUTINE get_restart_att_real_by_index(idx, attribute_name, attribute_value)
    INTEGER,          INTENT(in)  :: idx
    CHARACTER(len=*), INTENT(out) :: attribute_name
    REAL(wp),         INTENT(out) :: attribute_value
    attribute_name  = TRIM(restart_attributes_real(idx)%name)
    attribute_value = restart_attributes_real(idx)%val
  END SUBROUTINE get_restart_att_real_by_index
  !
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE get_restart_att_int_by_name(attribute_name, attribute_value, &
    &                                    opt_default)
    CHARACTER(len=*), INTENT(in)  :: attribute_name
    INTEGER,          INTENT(out) :: attribute_value
    INTEGER,          INTENT(IN), OPTIONAL :: opt_default
    INTEGER :: i
    DO i = 1, natts_int
      IF (TRIM(attribute_name) == TRIM(restart_attributes_int(i)%name)) THEN
        attribute_value = restart_attributes_int(i)%val
        RETURN
      ENDIF
    ENDDO
    IF (PRESENT(opt_default)) THEN
      attribute_value = opt_default
    ELSE
      CALL finish('','Attribute '//TRIM(attribute_name)//' not found.')
    END IF
  END SUBROUTINE get_restart_att_int_by_name
  !
  SUBROUTINE get_restart_att_int8_by_name(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in)  :: attribute_name
    INTEGER(i8),      INTENT(out) :: attribute_value
    INTEGER :: i
    DO i = 1, natts_int
      IF (TRIM(attribute_name) == TRIM(restart_attributes_int(i)%name)) THEN
        attribute_value = INT(restart_attributes_int(i)%val, i8)
        RETURN
      ENDIF
    ENDDO
    CALL finish('','Attribute '//TRIM(attribute_name)//' not found.')
  END SUBROUTINE get_restart_att_int8_by_name
  !
  SUBROUTINE get_restart_att_int_by_index(idx, attribute_name, attribute_value)
    INTEGER,          INTENT(in)  :: idx
    CHARACTER(len=*), INTENT(out) :: attribute_name
    INTEGER,          INTENT(out) :: attribute_value
    attribute_name  = TRIM(restart_attributes_int(idx)%name)
    attribute_value = restart_attributes_int(idx)%val
  END SUBROUTINE get_restart_att_int_by_index
  !
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE get_restart_att_bool_by_name(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in)  :: attribute_name
    LOGICAL,          INTENT(out) :: attribute_value
    INTEGER :: i
    DO i = 1, natts_bool
      IF (TRIM(attribute_name) == TRIM(restart_attributes_bool(i)%name)) THEN
        attribute_value = restart_attributes_bool(i)%val
        RETURN
      ENDIF
    ENDDO
    CALL finish('','Attribute '//TRIM(attribute_name)//' not found.')
  END SUBROUTINE get_restart_att_bool_by_name
  !
  SUBROUTINE get_restart_att_bool_by_index(idx, attribute_name, attribute_value)
    INTEGER,          INTENT(in)  :: idx
    CHARACTER(len=*), INTENT(out) :: attribute_name
    LOGICAL,          INTENT(out) :: attribute_value
    attribute_name  = TRIM(restart_attributes_bool(idx)%name)
    attribute_value = restart_attributes_bool(idx)%val
  END SUBROUTINE get_restart_att_bool_by_index
  !
  !------------------------------------------------------------------------------------------------
  !
  !
  !> Searches in a given attribute name list "list" for the name which
  !> is to be set. If the name already exists (case sensitive search)
  !> then return the index "idx" s.t. the current value is
  !> overwritten, otherwise expand the list of attributes and return
  !> the index of its tail.
  !
  SUBROUTINE find_or_expand_list(name, list, nlist, idx)
    CHARACTER(LEN=*), INTENT(IN)    :: name
    CHARACTER(LEN=*), INTENT(IN)    :: list(:)
    INTEGER,          INTENT(INOUT) :: nlist
    INTEGER,          INTENT(INOUT) :: idx
    ! local variables
    INTEGER :: i


    idx = -1
    loop : DO i=1,nlist 
      IF (TRIM(name) == TRIM(list(i))) THEN
        idx=i
      END IF
    END DO loop

    IF (idx == -1) THEN
      nlist = nlist+1
      idx   = nlist
      IF (nlist > nmax_atts) THEN
        CALL finish('find_or_expand_list','Too many restart attributes for restart file')
      END IF
    END IF
  END SUBROUTINE find_or_expand_list
  !
  SUBROUTINE set_restart_attribute_text(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in) :: attribute_name
    CHARACTER(len=*), INTENT(in) :: attribute_value
    INTEGER :: idx
    IF (.NOT. ALLOCATED(restart_attributes_text)) THEN
      ALLOCATE(restart_attributes_text(nmax_atts))
    ENDIF
    CALL find_or_expand_list(attribute_name, restart_attributes_text(:)%name, natts_text, idx)
    restart_attributes_text(idx)%name = TRIM(attribute_name)
    restart_attributes_text(idx)%val  = attribute_value
  END SUBROUTINE set_restart_attribute_text
  !
  SUBROUTINE set_restart_attribute_real(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in) :: attribute_name
    REAL(wp),         INTENT(in) :: attribute_value
    INTEGER :: idx
    IF (.NOT. ALLOCATED(restart_attributes_real)) THEN
      ALLOCATE(restart_attributes_real(nmax_atts))
    ENDIF
    CALL find_or_expand_list(attribute_name, restart_attributes_real(:)%name, natts_real, idx)
    restart_attributes_real(idx)%name = TRIM(attribute_name)
    restart_attributes_real(idx)%val  = attribute_value
  END SUBROUTINE set_restart_attribute_real
  !
  SUBROUTINE set_restart_attribute_int(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in) :: attribute_name
    INTEGER,          INTENT(in) :: attribute_value
    INTEGER :: idx
    IF (.NOT. ALLOCATED(restart_attributes_int)) THEN
      ALLOCATE(restart_attributes_int(nmax_atts))
    ENDIF
    CALL find_or_expand_list(attribute_name, restart_attributes_int(:)%name, natts_int, idx)
    restart_attributes_int(idx)%name = TRIM(attribute_name)
    restart_attributes_int(idx)%val = attribute_value
  END SUBROUTINE set_restart_attribute_int
  !
  SUBROUTINE set_restart_attribute_int8(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in) :: attribute_name
    INTEGER(i8),      INTENT(in) :: attribute_value
    INTEGER :: idx
    IF (.NOT. ALLOCATED(restart_attributes_int)) THEN
      ALLOCATE(restart_attributes_int(nmax_atts))
    ENDIF
    CALL find_or_expand_list(attribute_name, restart_attributes_int(:)%name, natts_int, idx)
    restart_attributes_int(idx)%name = TRIM(attribute_name)
    restart_attributes_int(idx)%val  = INT(attribute_value)
  END SUBROUTINE set_restart_attribute_int8
  !
  SUBROUTINE set_restart_attribute_bool(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in) :: attribute_name
    LOGICAL,          INTENT(in) :: attribute_value
    INTEGER :: idx
    IF (.NOT. ALLOCATED(restart_attributes_bool)) THEN
      ALLOCATE(restart_attributes_bool(nmax_atts))
    ENDIF
    CALL find_or_expand_list(attribute_name, restart_attributes_bool(:)%name, natts_bool, idx)

    restart_attributes_bool(idx)%name = 'bool_'//TRIM(attribute_name)
    restart_attributes_bool(idx)%val  = attribute_value
    ! for storing follows the C convention: false = 0, true = 1
    IF (attribute_value) THEN
      restart_attributes_bool(idx)%store = 1        
    ELSE
      restart_attributes_bool(idx)%store = 0        
    ENDIF
  END SUBROUTINE set_restart_attribute_bool
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE read_and_bcast_attributes(vlistID, lread_pe, root_pe, comm)
    INTEGER, INTENT(IN) :: vlistID      !< CDI vlist ID
    LOGICAL, INTENT(IN) :: lread_pe     !< .TRUE., if current PE has opened the file for reading
    INTEGER, INTENT(IN) :: root_pe      !< rank of broadcast root PE
    INTEGER, INTENT(IN) :: comm         !< MPI communicator
    !
    CHARACTER(len=256) :: att_name
    INTEGER :: natts, att_type, att_len, status, text_len, i
    REAL(KIND = C_DOUBLE) :: oneDouble(1)
    INTEGER(KIND = C_INT) :: oneInt(1)
    !
    IF (.NOT. ALLOCATED(restart_attributes_text)) THEN
      ALLOCATE(restart_attributes_text(nmax_atts))
      natts_text = 0
    ENDIF
    !
    IF (.NOT. ALLOCATED(restart_attributes_real)) THEN
      ALLOCATE(restart_attributes_real(nmax_atts))
      natts_real = 0 
    ENDIF
    !
    IF (.NOT. ALLOCATED(restart_attributes_int)) THEN
      ALLOCATE(restart_attributes_int(nmax_atts))
      natts_int  = 0 
    ENDIF
    !
    IF (.NOT. ALLOCATED(restart_attributes_bool)) THEN
      ALLOCATE(restart_attributes_bool(nmax_atts))
      natts_bool = 0 
    ENDIF
    !
    IF (lread_pe) THEN
      status = vlistInqNatts(vlistID, CDI_GLOBAL, natts)
    END IF
    CALL p_bcast(natts, root_pe, comm)
    !
    DO i = 0, natts-1
      IF (lread_pe) THEN
        status = vlistInqAtt(vlistID, CDI_GLOBAL, i, att_name, att_type, att_len)
      END IF
      CALL p_bcast(att_name, root_pe, comm)
      IF ( att_name(1:4) == 'nml_') CYCLE ! skip this, it is a namelist 

      CALL p_bcast(att_type, root_pe, comm)
      SELECT CASE(att_type)
      CASE(DATATYPE_FLT64)
        natts_real = natts_real+1
        IF (natts_real > SIZE(restart_attributes_real))  CALL finish("", "Too many restart arguments!")
        restart_attributes_real(natts_real)%name = TRIM(att_name)
        IF (lread_pe) THEN
          status = vlistInqAttFlt(vlistID, CDI_GLOBAL, TRIM(att_name), 1, oneDouble)
          restart_attributes_real(natts_real)%val = oneDouble(1)
        END IF
        CALL p_bcast(restart_attributes_real(natts_real)%val, root_pe, comm)
      CASE(DATATYPE_INT32)
        IF (att_name(1:5) == 'bool_') THEN
          natts_bool = natts_bool+1
          IF (natts_bool > SIZE(restart_attributes_bool))  CALL finish("", "Too many restart arguments!")
          restart_attributes_bool(natts_bool)%name = TRIM(att_name(6:))
          IF (lread_pe) THEN
            status = vlistInqAttInt(vlistID, CDI_GLOBAL, TRIM(att_name), 1, oneInt)
            restart_attributes_bool(natts_bool)%store = oneInt(1)
          END IF
          CALL p_bcast(restart_attributes_bool(natts_bool)%store, root_pe, comm)
          restart_attributes_bool(natts_bool)%val = (restart_attributes_bool(natts_bool)%store == 1)
        ELSE
          natts_int = natts_int+1
          IF (natts_int > SIZE(restart_attributes_int))  CALL finish("", "Too many restart arguments!")
          restart_attributes_int(natts_int)%name = TRIM(att_name)
          IF (lread_pe) THEN
            status = vlistInqAttInt(vlistID, CDI_GLOBAL, TRIM(att_name), 1, oneInt)
            restart_attributes_int(natts_int)%val = oneInt(1)
          END IF
          CALL p_bcast(restart_attributes_int(natts_int)%val, root_pe, comm)
        ENDIF
      CASE(DATATYPE_TXT)
        natts_text = natts_text+1
        IF (natts_text > SIZE(restart_attributes_text))  CALL finish("", "Too many restart arguments!")
        restart_attributes_text(natts_text)%name = TRIM(att_name)
        restart_attributes_text(natts_text)%val  = ''
        IF (lread_pe) THEN
          text_len = att_len
          status = vlistInqAttTxt(vlistID, CDI_GLOBAL, &
            &                   TRIM(att_name), text_len, restart_attributes_text(natts_text)%val)
        END IF
        CALL p_bcast(restart_attributes_text(natts_text)%val, root_pe, comm)
      END SELECT
    ENDDO
    !
  END SUBROUTINE read_and_bcast_attributes
  !
  SUBROUTINE delete_attributes
    IF (ALLOCATED(restart_attributes_text)) DEALLOCATE (restart_attributes_text)
    natts_text = 0
    IF (ALLOCATED(restart_attributes_real)) DEALLOCATE (restart_attributes_real)
    natts_real = 0 
    IF (ALLOCATED(restart_attributes_int)) DEALLOCATE (restart_attributes_int)
    natts_int  = 0 
    IF (ALLOCATED(restart_attributes_bool)) DEALLOCATE (restart_attributes_bool)
    natts_bool = 0 
  END SUBROUTINE delete_attributes
  !
END MODULE mo_io_restart_attributes
