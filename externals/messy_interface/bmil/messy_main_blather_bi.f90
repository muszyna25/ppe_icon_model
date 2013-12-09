! ************************************************************************
MODULE messy_main_blather_bi
! ************************************************************************
  
  USE messy_main_blather

  IMPLICIT NONE
  
  PUBLIC :: start_message_bi    ! standard messages for start and end ...
  PUBLIC :: end_message_bi      ! .. of submodel-specific MESSy-routines
  PUBLIC :: info_bi
  PUBLIC :: warning_bi
  PUBLIC :: error_bi
  PUBLIC :: messy_blather_endfile_bi

CONTAINS

  ! --------------------------------------------------------------------
  SUBROUTINE start_message_bi(modstr, str, substr)

    ! ECHAM5
    USE messy_main_mpi_bi,    ONLY: p_parallel_io

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN) :: modstr, str, substr

    CALL start_message(modstr, str, substr, p_parallel_io)
    
  END SUBROUTINE start_message_bi
  ! --------------------------------------------------------------------
  
  ! --------------------------------------------------------------------
  SUBROUTINE end_message_bi(modstr, str, substr)
    
    ! ECHAM5
    USE messy_main_mpi_bi,    ONLY: p_parallel_io

    IMPLICIT NONE
    
    ! I/O
    CHARACTER(LEN=*), INTENT(IN) :: modstr, str, substr

    CALL end_message(modstr, str, substr, p_parallel_io)
    
  END SUBROUTINE end_message_bi
  ! --------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE info_bi(string, substr)

    USE messy_main_mpi_bi,    ONLY: p_parallel_io

    IMPLICIT NONE

    CHARACTER(LEN=*),           INTENT(IN) :: string
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: substr

    CALL info(string, substr, l_print=p_parallel_io)

  END SUBROUTINE info_bi
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
  SUBROUTINE warning_bi(string, substr)

    USE messy_main_mpi_bi,    ONLY: p_parallel_io

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: string, substr

    CALL warning(string, substr, l_print=p_parallel_io)

  END SUBROUTINE warning_bi
!-------------------------------------------------------------------------
  
!-------------------------------------------------------------------------
  SUBROUTINE error_bi(string, substr)

    USE messy_main_mpi_bi, ONLY: p_abort, p_parallel, p_pe
    USE messy_main_tools,  ONLY: find_next_free_unit

    IMPLICIT NONE
    INTRINSIC :: TRIM

    CHARACTER(LEN=*), INTENT(IN) :: string, substr
    INTEGER          :: iou
    LOGICAL          :: lex
    CHARACTER(LEN=7) :: efile ! e.g. "END0001"

    WRITE (*,'(/,78("*"),/)')
    WRITE (*,'("ERROR: p_pe   = ",I4)') p_pe
    WRITE (*,'("ERROR: substr = ",A)')  TRIM(substr)
    WRITE (*,'("ERROR: ",A)')           TRIM(string)
    WRITE (*,'(/,78("*"),/)')

    ! mz_pj_20030707+
    ! write file 'END' for breaking rerun-chain in runscript xmessy
    iou = find_next_free_unit(100,300)
    WRITE(efile,'(a3,i4.4)') 'END',p_pe
    INQUIRE(FILE=efile, EXIST=lex) ! check if file exists
    IF (lex) THEN
      OPEN(iou, FILE=efile, STATUS='OLD', POSITION='APPEND') ! existing file
    ELSE
      OPEN(iou, FILE=efile, STATUS='UNKNOWN', POSITION='APPEND') ! new file
    ENDIF
    WRITE (iou,'("ERROR: p_pe   = ",I4)') p_pe
    WRITE (iou,'("ERROR: substr = ",A)')  TRIM(substr)
    WRITE (iou,'("ERROR: ",A)')           TRIM(string)
    CLOSE(iou)
    ! mz_pj_20030707-

#ifndef COSMO
   IF (p_parallel) THEN 
      CALL p_abort
    ELSE
      STOP
    ENDIF
#else
    IF (p_parallel) THEN 
       CALL p_abort(string,substr)
    ELSE
       STOP
    ENDIF
#endif

  END SUBROUTINE error_bi
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE messy_blather_endfile_bi(string, psubstr)

    USE messy_main_mpi_bi,    ONLY: p_pe
    USE messy_main_tools,  ONLY: find_next_free_unit
    USE messy_main_constants_mem, ONLY: STRLEN_ULONG ! um_ak_20101110

    IMPLICIT NONE
    
    INTRINSIC TRIM ! um_ak_20101110

    CHARACTER(LEN=*),           INTENT(IN) :: string
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: psubstr ! um_ak_20101110 renamed

    INTEGER          :: iou
    LOGICAL          :: lex
    CHARACTER(LEN=7) :: efile ! e.g. "END0001"
    ! um_ak_20101110+
    CHARACTER(LEN=STRLEN_ULONG) :: substr = ''
    
    IF (PRESENT(psubstr)) THEN
       substr = TRIM(psubstr)
    ELSE
       substr = ''
    ENDIF       
    ! um_ak_20101110-

    CALL info_bi(string//" (endfile_bi) ", substr)

    ! write file 'END' for breaking rerun-chain in runscript xmessy
    iou = find_next_free_unit(100,300)
    WRITE(efile,'(a3,i4.4)') 'END',p_pe
    INQUIRE(FILE=efile, EXIST=lex) ! check if file exists
    IF (lex) THEN
      OPEN(iou, FILE=efile, STATUS='OLD', POSITION='APPEND') ! existing file
    ELSE
      OPEN(iou, FILE=efile, STATUS='UNKNOWN', POSITION='APPEND') ! new file
    ENDIF
    WRITE (iou,'("END: p_pe   = ",I4)') p_pe
    WRITE (iou,'("END: substr = ",A)')  TRIM(substr)
    WRITE (iou,'("END: ",A)')           TRIM(string)
    CLOSE(iou)

  END SUBROUTINE messy_blather_endfile_bi
!-------------------------------------------------------------------------

! ************************************************************************
END MODULE messy_main_blather_bi
! ************************************************************************
