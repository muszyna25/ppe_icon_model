#if defined (__GNUC__) || defined (__INTEL_COMPILER) || defined (__SUNPRO_F95)
! do not use F2003
#else
#define HAVE_F2003
#endif
MODULE mo_io_restart_namelist
  !
  USE mo_util_file, ONLY: util_tmpnam, util_filesize, util_unlink
  USE mo_io_units,  ONLY: find_next_free_unit, filename_max
  USE mo_exception, ONLY: message, finish
  USE mo_cdi_constants
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: open_tmpfile
  PUBLIC :: store_and_close_namelist
  PUBLIC :: open_and_restore_namelist
  PUBLIC :: close_tmpfile
  PUBLIC :: read_restart_namelists
  PUBLIC :: nmls
  PUBLIC :: restart_namelist
#ifdef __SUNPRO_F95
  PUBLIC :: t_att_namelist  ! The Sun compiler does not allow a public 
                            ! variable declared as a private type
#endif
  !
  INTEGER, PARAMETER :: nmllen_max = 4096
  !
  TYPE t_att_namelist
    CHARACTER(len=64) :: name
#ifdef HAVE_F2003
    CHARACTER(len=:), ALLOCATABLE :: text
#else
    CHARACTER(len=nmllen_max) :: text
#endif
  END TYPE t_att_namelist
  !
  INTEGER, PARAMETER :: nmax_nmls = 64
  INTEGER, SAVE :: nmls = 0
  TYPE(t_att_namelist), ALLOCATABLE :: restart_namelist(:)
  !
CONTAINS
  !
  SUBROUTINE set_restart_namelist(namelist_name, namelist_text)
    CHARACTER(len=*), INTENT(in) :: namelist_name
    CHARACTER(len=*), INTENT(in) :: namelist_text
#ifdef HAVE_F2003
    INTEGER :: text_len
#endif
    IF (.NOT. ALLOCATED(restart_namelist)) THEN
      ALLOCATE(restart_namelist(nmax_nmls))
    ENDIF
    nmls = nmls+1
    IF (nmls > nmax_nmls) THEN
      CALL finish('set_restart_namelist', &
           &      'too many restart attributes for restart file')
    ELSE
      restart_namelist(nmls)%name = 'nml_'//namelist_name
#ifdef HAVE_F2003
      text_len = LEN_TRIM(namelist_text)
      IF (ALLOCATED(restart_namelist(nmls)%text)) DEALLOCATE(restart_namelist(nmls)%text) 
      ALLOCATE(CHARACTER(len=text_len) :: restart_namelist(nmls)%text)
#endif
      restart_namelist(nmls)%text = namelist_text
    ENDIF    
  END SUBROUTINE set_restart_namelist
  !
  SUBROUTINE get_restart_namelist(namelist_name, namelist_text)
    CHARACTER(len=*),              INTENT(in)  :: namelist_name
#ifdef HAVE_F2003
    INTEGER :: text_len
    CHARACTER(len=:), ALLOCATABLE, INTENT(out) :: namelist_text
#else
    CHARACTER(len=*),              INTENT(out) :: namelist_text
#endif
    !
    INTEGER :: i
    !
    DO i = 1, nmls
      IF (TRIM(namelist_name) == TRIM(restart_namelist(i)%name)) THEN
#ifdef HAVE_F2003
        text_len = LEN(restart_namelist(i)%text)
        IF (ALLOCATED(namelist_text)) DEALLOCATE(namelist_text)
        ALLOCATE(CHARACTER(len=text_len) :: namelist_text)
#endif
        namelist_text = restart_namelist(i)%text
        RETURN
      ENDIF
    ENDDO
    !
    CALL finish('','namelist '//TRIM(namelist_name)//' not available in restart file.')
    !
  END SUBROUTINE get_restart_namelist
  !
  FUNCTION open_tmpfile() RESULT(funit)
    INTEGER :: funit
    !
    INTEGER :: flen
    CHARACTER(len=filename_max) :: filename
    !
    flen = util_tmpnam(filename, filename_max)
    funit = find_next_free_unit(10,100)
    OPEN(UNIT=funit, FILE=filename(1:flen),   &
         ACTION='write', ACCESS='sequential', &
         DELIM='apostrophe')
    !
  END FUNCTION open_tmpfile
  !
  SUBROUTINE store_and_close_namelist(funit, name)
    INTEGER,          INTENT(in) :: funit
    CHARACTER(len=*), INTENT(in) :: name 
    !
    CHARACTER(len=filename_max) :: filename
    INTEGER :: nmllen
#ifdef HAVE_F2003
    CHARACTER(len=:), ALLOCATABLE :: nmlbuf
#else
    CHARACTER(len=nmllen_max) :: nmlbuf
#endif
    INTEGER :: iret
    !
    INQUIRE(funit,NAME=filename)
    !
    CLOSE(funit)
    !
    nmllen = util_filesize(filename)
    IF (nmllen == 0) THEN
      CALL message('','namelist '//TRIM(name)//' is empty, saving in restart file fails.')
    ENDIF
#ifdef HAVE_F2003
    ALLOCATE(CHARACTER(len=nmllen) :: nmlbuf)
#else
    IF (nmllen > nmllen_max) THEN
      CALL message('', &
           'The problem could be solved by increasing nmllen_max in mo_io_restart_namelist.f90.')
      CALL finish('','namelist '//TRIM(name)//' is too long, saving in restart file fails.')
    ENDIF
#endif
    !
    OPEN(UNIT=funit, FILE=TRIM(filename), ACTION='read', ACCESS='stream')
    READ(funit) nmlbuf(1:nmllen)
    CLOSE(funit)
    !
    CALL set_restart_namelist(TRIM(name), nmlbuf)
    !
#ifdef HAVE_F2003
    DEALLOCATE(nmlbuf)
#endif
    !
    iret = util_unlink(TRIM(filename))
    !
  END SUBROUTINE store_and_close_namelist
  !
  FUNCTION open_and_restore_namelist(name) RESULT(funit)
    INTEGER :: funit
    CHARACTER(len=*), INTENT(in) :: name 
    !
    INTEGER :: flen
    CHARACTER(len=filename_max) :: filename
    !
#ifdef HAVE_F2003
    CHARACTER(len=:), ALLOCATABLE :: nmlbuf
#else
    CHARACTER(len=nmllen_max) :: nmlbuf
#endif
    !
    CALL get_restart_namelist(name, nmlbuf)
    !
    funit = find_next_free_unit(10,100)
    flen = util_tmpnam(filename, filename_max)
    OPEN(UNIT=funit, FILE=filename(1:flen), ACTION='write', ACCESS='stream')
    WRITE(funit) TRIM(nmlbuf)
    CLOSE(funit)
    !
#ifdef HAVE_F2003
    DEALLOCATE(nmlbuf)
#endif
    !
    OPEN(UNIT=funit, FILE=filename(1:flen), ACTION='read', ACCESS='sequential')
    !    
  END FUNCTION open_and_restore_namelist
  !
  SUBROUTINE close_tmpfile(funit)
    INTEGER, INTENT(in) :: funit
    CHARACTER(len=filename_max) :: filename
    INTEGER :: iret
    !
    INQUIRE(funit,NAME=filename)
    !
    CLOSE(funit)
    !
    iret = util_unlink(TRIM(filename))
    !
  END SUBROUTINE close_tmpfile
  !
  SUBROUTINE read_restart_namelists(filename)
    CHARACTER(len=*), INTENT(in) :: filename
    INTEGER :: fileID, vlistID
    CHARACTER(len=64) :: att_name
    INTEGER :: natts, att_type, att_len
    INTEGER :: nmllen
#ifdef HAVE_F2003
    CHARACTER(len=:), ALLOCATABLE :: nmlbuf
#else
    CHARACTER(len=nmllen_max) :: nmlbuf
#endif
    INTEGER :: i, status
    !
    IF (.NOT. ALLOCATED(restart_namelist)) THEN
      ALLOCATE(restart_namelist(nmax_nmls))
    ENDIF
    nmls = 0
    !
    fileID  = streamOpenRead(TRIM(filename))
    vlistID = streamInqVlist(fileID)
    !
    status = vlistInqNatts(vlistID, CDI_GLOBAL, natts)
    !
    DO i = 0, natts-1
      status = vlistInqAtt(vlistID, CDI_GLOBAL, i, att_name, att_type, att_len)
      IF ( att_name(1:4) /= 'nml_') CYCLE ! skip this, it is not a namelist 
      nmllen = att_len
#ifdef HAVE_F2003
      ALLOCATE(CHARACTER(len=nmllen) :: nmlbuf)
#else
      IF (nmllen > nmllen_max) THEN
        CALL message('', 'The problem could be solved '//                     &
                    &'by increasing nmllen_max in mo_io_restart_namelist.f90.')
        CALL finish('','namelist '//TRIM(att_name)//                          &
                   &' is too long, reload from restart file fails.')
      ENDIF
#endif
      status = vlistInqAttTxt(vlistID, CDI_GLOBAL, TRIM(att_name), nmllen, nmlbuf)
      !
      nmls = nmls+1
      IF (nmls > nmax_nmls) THEN
        CALL finish('set_restart_namelist', &
             &      'too many restart attributes for restart file')
      ENDIF
      restart_namelist(nmls)%name = TRIM(att_name(5:))
#ifdef HAVE_F2003
      IF (ALLOCATED(restart_namelist(nmls)%text)) DEALLOCATE(restart_namelist(nmls)%text) 
      ALLOCATE(CHARACTER(len=nmllen) :: restart_namelist(nmls)%text)
#endif
      restart_namelist(nmls)%text = nmlbuf
#ifdef HAVE_F2003
      DEALLOCATE(nmlbuf)
#endif
    ENDDO
    CALL streamClose(fileID)
    !
  END SUBROUTINE read_restart_namelists

END MODULE mo_io_restart_namelist
