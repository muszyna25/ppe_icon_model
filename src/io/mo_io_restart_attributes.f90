MODULE mo_io_restart_attributes
  !
  USE mo_kind,                  ONLY: wp, i8
  USE mo_exception,             ONLY: finish
  USE mo_cdi_constants
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: set_restart_attribute
  PUBLIC :: get_restart_attribute
  PUBLIC :: read_restart_attributes
  PUBLIC :: read_attributes
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
    CHARACTER(len=64) :: name
    CHARACTER(len=64) :: val
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
  INTEGER, PARAMETER :: nmax_atts = 256
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
  SUBROUTINE get_restart_att_int_by_name(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in)  :: attribute_name
    INTEGER,          INTENT(out) :: attribute_value
    INTEGER :: i
    DO i = 1, natts_int
      IF (TRIM(attribute_name) == TRIM(restart_attributes_int(i)%name)) THEN
        attribute_value = restart_attributes_int(i)%val
        RETURN
      ENDIF
    ENDDO
    CALL finish('','Attribute '//TRIM(attribute_name)//' not found.')
  END SUBROUTINE get_restart_att_int_by_name
  !
  SUBROUTINE get_restart_att_int8_by_name(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in)  :: attribute_name
    INTEGER(i8),      INTENT(out) :: attribute_value
    INTEGER :: i
    DO i = 1, natts_int
      IF (TRIM(attribute_name) == TRIM(restart_attributes_int(i)%name)) THEN
        attribute_value = restart_attributes_int(i)%val
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
  SUBROUTINE set_restart_attribute_text(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in) :: attribute_name
    CHARACTER(len=*), INTENT(in) :: attribute_value
    IF (.NOT. ALLOCATED(restart_attributes_text)) THEN
      ALLOCATE(restart_attributes_text(nmax_atts))
    ENDIF
    natts_text = natts_text+1
    IF (natts_text > nmax_atts) THEN
      CALL finish('set_restart_attribute_text','too many restart attributes for restart file')
    ELSE
      restart_attributes_text(natts_text)%name = attribute_name
      restart_attributes_text(natts_text)%val = attribute_value
    ENDIF    
  END SUBROUTINE set_restart_attribute_text
  !
  SUBROUTINE set_restart_attribute_real(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in) :: attribute_name
    REAL(wp),         INTENT(in) :: attribute_value
    IF (.NOT. ALLOCATED(restart_attributes_real)) THEN
      ALLOCATE(restart_attributes_real(nmax_atts))
    ENDIF
    natts_real = natts_real+1
    IF (natts_real > nmax_atts) THEN
      CALL finish('set_restart_attribute_real','too many restart attributes for restart file')
    ELSE
      restart_attributes_real(natts_real)%name = attribute_name
      restart_attributes_real(natts_real)%val  = attribute_value
    ENDIF    
  END SUBROUTINE set_restart_attribute_real
  !
  SUBROUTINE set_restart_attribute_int(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in) :: attribute_name
    INTEGER,          INTENT(in) :: attribute_value
    IF (.NOT. ALLOCATED(restart_attributes_int)) THEN
      ALLOCATE(restart_attributes_int(nmax_atts))
    ENDIF
    natts_int = natts_int+1
    IF (natts_int > nmax_atts) THEN
      CALL finish('set_restart_attribute_int','too many restart attributes for restart file')
    ELSE
      restart_attributes_int(natts_int)%name = attribute_name
      restart_attributes_int(natts_int)%val = attribute_value
    ENDIF    
  END SUBROUTINE set_restart_attribute_int
  !
  SUBROUTINE set_restart_attribute_int8(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in) :: attribute_name
    INTEGER(i8),      INTENT(in) :: attribute_value
    IF (.NOT. ALLOCATED(restart_attributes_int)) THEN
      ALLOCATE(restart_attributes_int(nmax_atts))
    ENDIF
    natts_int = natts_int+1
    IF (natts_int > nmax_atts) THEN
      CALL finish('set_restart_attribute_int','too many restart attributes for restart file')
    ELSE
      restart_attributes_int(natts_int)%name = attribute_name
      restart_attributes_int(natts_int)%val = attribute_value
    ENDIF    
  END SUBROUTINE set_restart_attribute_int8
  !
  SUBROUTINE set_restart_attribute_bool(attribute_name, attribute_value)
    CHARACTER(len=*), INTENT(in) :: attribute_name
    LOGICAL,          INTENT(in) :: attribute_value
    IF (.NOT. ALLOCATED(restart_attributes_bool)) THEN
      ALLOCATE(restart_attributes_bool(nmax_atts))
    ENDIF
    natts_bool = natts_bool+1
    IF (natts_bool > nmax_atts) THEN
      CALL finish('set_restart_attribute_bools','too many restart attributes for restart file')
    ELSE
      restart_attributes_bool(natts_bool)%name = 'bool_'//TRIM(attribute_name)
      restart_attributes_bool(natts_bool)%val = attribute_value
      ! for storing follows the C convention: false = 0, true = 1
      IF (attribute_value) THEN
        restart_attributes_bool(natts_bool)%store = 1        
      ELSE
        restart_attributes_bool(natts_bool)%store = 0        
      ENDIF
    ENDIF    
  END SUBROUTINE set_restart_attribute_bool
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE read_restart_attributes(filename)
    CHARACTER(len=*), INTENT(in) :: filename
    INTEGER :: fileID, vlistID
    !
    fileID = streamOpenRead(TRIM(filename))
    vlistID = streamInqVlist(fileID)
    !
    CALL read_attributes(vlistID)
    !
    CALL streamClose(fileID)
    !
  END SUBROUTINE read_restart_attributes
  !
  SUBROUTINE read_attributes(vlistID)
    INTEGER, INTENT(in) :: vlistID
    !
    CHARACTER(len=64) :: att_name
    INTEGER :: natts, att_type, att_len
    INTEGER :: status, text_len, mlen
    !
    INTEGER :: i
    !
    status = vlistInqNatts(vlistID, CDI_GLOBAL, natts)
    !
    IF (.NOT. ALLOCATED(restart_attributes_text)) THEN
      ALLOCATE(restart_attributes_text(nmax_atts))
    ENDIF
    natts_text = 0
    !
    IF (.NOT. ALLOCATED(restart_attributes_real)) THEN
      ALLOCATE(restart_attributes_real(nmax_atts))
    ENDIF
    natts_real = 0 
    !
    IF (.NOT. ALLOCATED(restart_attributes_int)) THEN
      ALLOCATE(restart_attributes_int(nmax_atts))
    ENDIF
    natts_int  = 0 
    !
    IF (.NOT. ALLOCATED(restart_attributes_bool)) THEN
      ALLOCATE(restart_attributes_bool(nmax_atts))
    ENDIF
    natts_bool = 0 
    !
    DO i = 0, natts-1
      status = vlistInqAtt(vlistID, CDI_GLOBAL, i, att_name, att_type, att_len)
      IF ( att_name(1:4) == 'nml_') CYCLE ! skip this, it is a namelist 
      SELECT CASE(att_type)
      CASE(DATATYPE_FLT64)
        natts_real = natts_real+1
        restart_attributes_real(natts_real)%name = TRIM(att_name)
        mlen = 1
        status =  vlistInqAttFlt(vlistID, CDI_GLOBAL, &              
             &                   TRIM(att_name), mlen, restart_attributes_real(natts_real)%val)
      CASE(DATATYPE_INT32)
        IF (att_name(1:5) == 'bool_') THEN
          natts_bool = natts_bool+1
          restart_attributes_bool(natts_bool)%name = TRIM(att_name(6:))
          mlen = 1
          status =  vlistInqAttInt(vlistID, CDI_GLOBAL, &              
               &                   TRIM(att_name), mlen, restart_attributes_bool(natts_bool)%store)
          restart_attributes_bool(natts_bool)%val = &
               (restart_attributes_bool(natts_bool)%store == 1)
        ELSE
          natts_int = natts_int+1
          restart_attributes_int(natts_int)%name = TRIM(att_name)
          mlen = 1
          status =  vlistInqAttInt(vlistID, CDI_GLOBAL, &              
               &                   TRIM(att_name), mlen, restart_attributes_int(natts_int)%val)
        ENDIF
      CASE(DATATYPE_TXT)
        natts_text = natts_text+1
        restart_attributes_text(natts_text)%name = TRIM(att_name)
        text_len = att_len
        restart_attributes_text(natts_text)%val = ''
        status =  vlistInqAttTxt(vlistID, CDI_GLOBAL, &
             &                   TRIM(att_name), text_len, restart_attributes_text(natts_text)%val)
      END SELECT
    ENDDO
    !
  END SUBROUTINE read_attributes
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
