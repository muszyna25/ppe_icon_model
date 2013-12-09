! **********************************************************************
MODULE messy_main_channel_attributes
! **********************************************************************

  ! MESSY DATA TRANSFER AND EXPORT INTERFACE (MEMORY MANAGEMENT)
  !
  ! Author: Patrick Joeckel, MPICH, May 2005

  USE messy_main_constants_mem, ONLY: DP  &
       , STRLEN_MEDIUM, STRLEN_ULONG

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  ! ATTRIBUTE TYPES
  INTEGER, PARAMETER, PUBLIC :: TYPE_UNKNOWN   = -1
  INTEGER, PARAMETER, PUBLIC :: TYPE_INTEGER   = 0
  INTEGER, PARAMETER, PUBLIC :: TYPE_STRING    = 1
  INTEGER, PARAMETER, PUBLIC :: TYPE_REAL_DP   = 2

  ! ATTRIBUTE FLAGS
  INTEGER, PARAMETER, PUBLIC :: AF_NONE = 0
!  INTEGER, PARAMETER, PUBLIC :: AF_RST_INF = 1
!  INTEGER, PARAMETER, PUBLIC :: AF_RST_ERR = 2
!  INTEGER, PARAMETER, PUBLIC :: AF_RST_INP = 3

  TYPE t_attribute
     CHARACTER(LEN=STRLEN_MEDIUM) :: name    = ''
     INTEGER                      :: type    = TYPE_UNKNOWN
     INTEGER                      :: iflag   = AF_NONE
     !
     INTEGER                      :: i = 0
     CHARACTER(LEN=STRLEN_ULONG)  :: c = ''
     REAL(DP)                     :: r = 0.0_DP
  END TYPE t_attribute

  TYPE t_attribute_list
!     PRIVATE
     TYPE(t_attribute)               :: this
     TYPE(t_attribute_list), POINTER :: next => NULL()
  END TYPE t_attribute_list

  ! TYPES
  PUBLIC :: t_attribute            ! TYPE: dimension
  PUBLIC :: t_attribute_list       ! TYPE: concat. list of attributes

  ! INTERFACES
  PUBLIC :: add_attribute          ! add attribute to list  
  INTERFACE write_attribute
     MODULE PROCEDURE write_attribute_one
     MODULE PROCEDURE write_attribute_all
  END INTERFACE
  PUBLIC :: write_attribute        ! print attribute (list) to standard output
  PUBLIC :: return_attribute       ! get attribute from list
  PUBLIC :: delete_attribute       ! remove attribute from list
  PUBLIC :: copy_attribute_list    ! copy complete attribute list
  PUBLIC :: clean_attribute_list   ! delete attribute list
  !PRIVATE :: loc_attribute        ! locate pointer to attribute
  
CONTAINS

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE add_attribute(status, list, name, i, c, r, loverwrite, iflag)
    
    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, LEN_TRIM, ASSOCIATED, PRESENT, TRIM

    ! I/O
    INTEGER,                   INTENT(OUT)          :: status
    TYPE(t_attribute_list),    POINTER              :: list
    CHARACTER(LEN=*),          INTENT(IN)           :: name
    INTEGER,                   INTENT(IN), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(IN), OPTIONAL :: c
    REAL(DP),                  INTENT(IN), OPTIONAL :: r
    LOGICAL,                   INTENT(IN), OPTIONAL :: loverwrite
    INTEGER,                   INTENT(IN), OPTIONAL :: iflag

    ! LOCAL
    INTEGER                         :: nopt  ! number of opt. parameters
    TYPE(t_attribute_list), POINTER :: ai  => NULL()
    TYPE(t_attribute_list), POINTER :: ae  => NULL()
    TYPE(t_attribute),      POINTER :: att => NULL()
    INTEGER                         :: zstat
    LOGICAL                         :: zloverwrite

    ! CHECKS
    nopt = 0
    IF (PRESENT(i)) nopt = nopt + 1
    IF (PRESENT(c)) nopt = nopt + 1
    IF (PRESENT(r)) nopt = nopt + 1
    IF (nopt /=1 ) THEN
       status = 801   ! ATTRIBUTE TYPE IS AMBIGUOUS
       RETURN
    END IF
    !
    IF (PRESENT(c)) THEN
       IF (LEN_TRIM(ADJUSTL(c)) > STRLEN_ULONG) THEN
          status = 803   ! CHARACTER ATTRIBUTE TOO LONG
          RETURN
       END IF
    END IF

    IF (PRESENT(loverwrite)) THEN
       zloverwrite = loverwrite
    ELSE
       zloverwrite = .FALSE.  ! DEFAULT
    END IF

    CALL loc_attribute(zstat, list, name, att)
    IF (zstat /= 805) THEN  ! ATTRIBUTE DOES NOT EXIST (IS OK HERE !)
       IF (zstat == 0) THEN ! ATTRIBUTE EXISTS ALREADY
          IF (zloverwrite) THEN
             CALL delete_attribute(status, list, name)
             IF (status /= 0) RETURN
          ELSE
             status = 804  ! ATTRIBUTE EXISTS ALREADY
             RETURN
          END IF
       ELSE
          status = zstat   ! ERROR
          RETURN
       END IF
    END IF

    ! GOTO END OF LIST
    ai => list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       ae => ai
       ai => ai%next
    END DO

    ! ADD NEW
    ALLOCATE(ai)
    NULLIFY(ai%next)
    IF (.NOT. ASSOCIATED(list)) THEN
       list => ai                      ! SET POINTER TO FIRST ELEMENT
    ELSE
       ae%next => ai                   ! SET NEXT POINTER OF LAST ELEMENT
       !                               ! TO NEW ELEMENT
    END IF
       
    ! SET VALUES
    ai%this%name = TRIM(ADJUSTL(name))
    !
    IF (PRESENT(iflag)) THEN
       ai%this%iflag = iflag
    END IF
    !
    IF (PRESENT(c)) THEN
       ai%this%type = TYPE_STRING
       ai%this%c    = TRIM(ADJUSTL(c))
    END IF
    !
    IF (PRESENT(i)) THEN
       ai%this%type = TYPE_INTEGER
       ai%this%i    = i
    END IF
    !
    IF (PRESENT(r)) THEN
       ai%this%type = TYPE_REAL_DP
       ai%this%r    = r
    END IF

    status = 0

  END SUBROUTINE add_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE return_attribute(status, list, name, i, c, r, iflag)

    IMPLICIT NONE

    INTRINSIC :: PRESENT

    ! I/O
    INTEGER,                   INTENT(OUT)           :: status
    TYPE(t_attribute_list),    POINTER               :: list
    CHARACTER(LEN=*),          INTENT(IN)            :: name
    INTEGER,                   INTENT(OUT), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(OUT), OPTIONAL :: c
    REAL(DP),                  INTENT(OUT), OPTIONAL :: r
    INTEGER,                   INTENT(OUT), OPTIONAL :: iflag

    ! LOCAL
    TYPE(t_attribute), POINTER :: att

    CALL loc_attribute(status, list, name, att)
    IF (status /= 0) RETURN

    IF (PRESENT(i)) i = att%i
    IF (PRESENT(c)) c = att%c
    IF (PRESENT(r)) r = att%r

    IF (PRESENT(iflag)) iflag = att%iflag

  END SUBROUTINE return_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE delete_attribute(status, list, name)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, TRIM

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_attribute_list), POINTER     :: list
    CHARACTER(LEN=*),       INTENT(IN)  :: name

    ! LOCAL
    TYPE(t_attribute_list), POINTER :: ai  => NULL()
    TYPE(t_attribute_list), POINTER :: ae  => NULL()
    TYPE(t_attribute),      POINTER :: att => NULL()

    CALL loc_attribute(status, list, name, att)
    IF (status /= 0) RETURN
    
    ! GOTO POSITION
    ai => list
    NULLIFY(ae)
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       IF (TRIM(ADJUSTL(name)) == TRIM(ai%this%name)) THEN
          EXIT
       END IF
       ae => ai
       ai => ai%next
    END DO    

    IF (ASSOCIATED(ae)) THEN
       ae%next => ai%next
    ELSE
       list => list%next
    END IF

    ! DELETE ELEMENT
    DEALLOCATE(ai)

    status = 0

  END SUBROUTINE delete_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE copy_attribute_list(status, list1, list2)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_attribute_list), POINTER     :: list1    ! source
    TYPE(t_attribute_list), POINTER     :: list2    ! destination

    ! LOCAL
    TYPE(t_attribute_list), POINTER :: ai  => NULL()
    TYPE(t_attribute),      POINTER :: att => NULL()

    ! CHECKS
    IF (ASSOCIATED(list2)) THEN
       status = 807 ! ATTRIBUTE LIST IS ALREADY ASSOCIATED
       RETURN
    END IF

    ai => list1
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT

       att => ai%this

       SELECT CASE(ai%this%type)
       CASE(TYPE_STRING)
          CALL add_attribute(status, list2, TRIM(att%name), c=TRIM(att%c) &
               , iflag=att%iflag)
       CASE(TYPE_INTEGER)
          CALL add_attribute(status, list2, TRIM(att%name), i=att%i &
               , iflag=att%iflag)
       CASE(TYPE_REAL_DP)
          CALL add_attribute(status, list2, TRIM(att%name), r=att%r &
               , iflag=att%iflag)
       CASE DEFAULT
          status = 806  ! UNKNOWN ATTRIBUTE TYPE
       END SELECT

       IF (status /= 0) RETURN

       ai => ai%next
    END DO    

    status = 0

  END SUBROUTINE copy_attribute_list
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE clean_attribute_list(status, list)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_attribute_list), POINTER     :: list

    ! LOCAL
    TYPE(t_attribute_list), POINTER :: ai => NULL()
    TYPE(t_attribute_list), POINTER :: ae => NULL()

    ai => list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       ae => ai
       ai => ai%next
       DEALLOCATE(ae)
       NULLIFY(ae)
    END DO    

    NULLIFY(list)
    
    status = 0
    
  END SUBROUTINE clean_attribute_list
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_attribute_all(status, list)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_attribute_list), POINTER     :: list

    ! LOCAL
    TYPE(t_attribute_list), POINTER :: ai => NULL()

    ! EMPTY LIST ?
    IF (.NOT. ASSOCIATED(list)) THEN

       WRITE(*,*) '   |    |->  *** ATTRIBUTE LIST IS EMPTY ***'

    ELSE

       ai => list
       DO
          IF (.NOT. ASSOCIATED(ai)) EXIT
          !
          CALL write_attribute_one(status, ai%this)
          IF (status /= 0) RETURN
          !
          ai => ai%next
       END DO
       
    END IF
    
    status = 0

  END SUBROUTINE write_attribute_all
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_attribute_one(status, att)

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,           INTENT(OUT) :: status
    TYPE(t_attribute), INTENT(IN)  :: att

    SELECT CASE(att%type)
    CASE(TYPE_STRING)
       WRITE(*,*) '   |    |    |-> ',att%name,': ', TRIM(att%c), &
            ' (',att%iflag,')'
    CASE(TYPE_INTEGER)
       WRITE(*,*) '   |    |    |-> ',att%name,': ', att%i, &
            ' (',att%iflag,')'
    CASE(TYPE_REAL_DP)
       WRITE(*,*) '   |    |    |-> ',att%name,': ', att%r, &
            ' (',att%iflag,')'
    CASE DEFAULT
       status = 806  ! UNKNOWN ATTRIBUTE TYPE
       RETURN
    END SELECT

    status = 0

  END SUBROUTINE write_attribute_one
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE loc_attribute(status, list, name, att)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, LEN_TRIM, TRIM, ASSOCIATED

    ! I/O
    INTEGER,                   INTENT(OUT)          :: status
    TYPE(t_attribute_list),    POINTER              :: list
    CHARACTER(LEN=*),          INTENT(IN)           :: name
    TYPE(t_attribute),         POINTER              :: att

    ! LOCAL
    TYPE(t_attribute_list), POINTER :: ai  => NULL()
    TYPE(t_attribute_list), POINTER :: ae  => NULL()
    LOGICAL                         :: lex

    ! INIT
    lex = .FALSE.

    ! CHECKS
    IF (LEN_TRIM(ADJUSTL(name)) > STRLEN_MEDIUM) THEN
       status = 802   ! ATTRIBUTE NAME TOO LONG
       RETURN
    END IF
    !
    IF (.NOT. ASSOCIATED(list)) THEN
       status = 805  ! ATTRIBUTE DOES NOT EXIST
       RETURN
    END IF

    ! CHECK, IF IT EXISTS
    ai => list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       IF (TRIM(ADJUSTL(name)) == TRIM(ai%this%name)) THEN
          lex = .TRUE.
          EXIT
       END IF
       ae => ai
       ai => ai%next
    END DO    

    IF (lex) THEN
       att => ai%this
    ELSE
       status = 805   ! ATTRIBUTE DOES NOT EXIST
       RETURN
    END IF

    status = 0

  END SUBROUTINE loc_attribute
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_main_channel_attributes
! **********************************************************************
