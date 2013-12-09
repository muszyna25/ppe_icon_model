! **********************************************************************
MODULE messy_main_channel_dimvar
! **********************************************************************

  ! MESSY DATA TRANSFER AND EXPORT INTERFACE (MEMORY MANAGEMENT)
  !
  ! Author: Patrick Joeckel, MPICH, May 2005

  USE messy_main_constants_mem,      ONLY: DP, STRLEN_MEDIUM
  USE messy_main_channel_attributes, ONLY: t_attribute_list

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  TYPE t_dimvar
     CHARACTER(LEN=STRLEN_MEDIUM)    :: name  = ''    ! NAME OF VARIABLE
     REAL(DP), DIMENSION(:), POINTER :: val => NULL() ! VALUES
     TYPE(t_attribute_list), POINTER :: att => NULL() ! ATTRIBUTE LIST
  END TYPE t_dimvar

  TYPE t_dimvar_list
     TYPE(t_dimvar)               :: this
     TYPE(t_dimvar_list), POINTER :: next => NULL()
  END TYPE t_dimvar_list

  ! TYPES
  PUBLIC :: t_dimvar            ! TYPE: dimension
  PUBLIC :: t_dimvar_list       ! TYPE: concat. list of dimvars

  ! INTERFACES
  PUBLIC :: add_dimvar          ! add dimvar to list
  PUBLIC :: add_dimvar_att      ! add dimvar attribute
  INTERFACE write_dimvar
     MODULE PROCEDURE write_dimvar_one
     MODULE PROCEDURE write_dimvar_all
  END INTERFACE
  PUBLIC :: write_dimvar        ! print dimvar (list) to standard output
  PUBLIC :: get_dimvar          ! locate pointer to dimvar (by name)
  PUBLIC :: delete_dimvar       ! remove dimvar from list
  PUBLIC :: clean_dimvar_list   ! delete dimvar list
  !PRIVATE :: loc_dimvar        ! locate pointer to attribute
  

CONTAINS

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE add_dimvar(status, list, name, val)
    
    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, TRIM, SIZE, ASSOCIATED

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_dimvar_list),    POINTER     :: list
    CHARACTER(LEN=*),       INTENT(IN)  :: name
    REAL(DP), DIMENSION(:), INTENT(IN)  :: val

    ! LOCAL
    TYPE(t_dimvar_list), POINTER :: ai     => NULL()
    TYPE(t_dimvar_list), POINTER :: ae     => NULL()
    TYPE(t_dimvar),      POINTER :: dimvar => NULL()
    INTEGER                      :: zstat

    CALL loc_dimvar(zstat, list, name, dimvar)
    IF (zstat /= 953) THEN  ! DIMENSION VARIABLE DOES NOT EXIST (IS OK HERE !)
       IF (zstat == 0) THEN ! DIMENSION VARIABLE EXISTS ALREADY
          status = 952      ! DIMENSION VARIABLE EXISTS ALREADY
       ELSE
          status = zstat    ! ERROR
       END IF
       RETURN
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
    ALLOCATE(ai%this%val(SIZE(val)), STAT=zstat)
    IF (zstat /= 0) THEN
       status = 1000 ! MEMORY ALLOCATION FAILED
       RETURN
    END IF
    ai%this%val(:) = val(:)

    status = 0

  END SUBROUTINE add_dimvar
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE add_dimvar_att(status, list, name, aname, i, c, r &
       , loverwrite, iflag)

    USE messy_main_channel_attributes, ONLY: add_attribute

    IMPLICIT NONE

    ! I/O
    INTEGER,                   INTENT(OUT)          :: status
    TYPE(t_dimvar_list),       POINTER              :: list
    CHARACTER(LEN=*),          INTENT(IN)           :: name
    CHARACTER(LEN=*),          INTENT(IN)           :: aname
    INTEGER,                   INTENT(IN), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(IN), OPTIONAL :: c
    REAL(DP),                  INTENT(IN), OPTIONAL :: r
    LOGICAL,                   INTENT(IN), OPTIONAL :: loverwrite
    INTEGER,                   INTENT(IN), OPTIONAL :: iflag

    ! LOCAL
    TYPE(t_dimvar), POINTER :: dimvar

    CALL loc_dimvar(status, list, name, dimvar)
    IF (status /= 0) RETURN

    CALL add_attribute(status, dimvar%att, aname, i, c, r, loverwrite, iflag)
    IF (status /= 0) RETURN

  END SUBROUTINE add_dimvar_att
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_dimvar(status, list, name, dimvar)

    IMPLICIT NONE

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_dimvar_list),    POINTER     :: list
    CHARACTER(LEN=*),       INTENT(IN)  :: name
    TYPE(t_dimvar),         POINTER     :: dimvar

    CALL loc_dimvar(status, list, name, dimvar)
    IF (status /= 0) RETURN

  END SUBROUTINE get_dimvar
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE delete_dimvar(status, list, name)

    USE messy_main_channel_attributes, ONLY: clean_attribute_list

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, TRIM, ASSOCIATED

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_dimvar_list),    POINTER     :: list
    CHARACTER(LEN=*),       INTENT(IN)  :: name

    ! LOCAL
    TYPE(t_dimvar_list), POINTER :: ai     => NULL()
    TYPE(t_dimvar_list), POINTER :: ae     => NULL()
    TYPE(t_dimvar),      POINTER :: dimvar => NULL()

    CALL loc_dimvar(status, list, name, dimvar)
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
    ! REMOVE CONTENTS
    CALL clean_attribute_list(status, ai%this%att)
    IF (status /= 0) RETURN
    IF (ASSOCIATED(ai%this%val)) DEALLOCATE(ai%this%val)

    ! DELETE ELEMENT
    DEALLOCATE(ai)
    
    status = 0
    
  END SUBROUTINE delete_dimvar
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE clean_dimvar_list(status, list)

    USE messy_main_channel_attributes, ONLY: clean_attribute_list

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED 

    ! I/O
    INTEGER,             INTENT(OUT) :: status
    TYPE(t_dimvar_list), POINTER     :: list

    ! LOCAL
    TYPE(t_dimvar_list), POINTER :: ai => NULL()
    TYPE(t_dimvar_list), POINTER :: ae => NULL()

    ai => list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       ae => ai
       ai => ai%next
       ! REMOVE CONTENTS
       CALL clean_attribute_list(status, ae%this%att)
       IF (status /= 0) RETURN
       IF (ASSOCIATED(ae%this%val)) DEALLOCATE(ae%this%val)      
       ! DELETE ELEMENT
       DEALLOCATE(ae)
       NULLIFY(ae)
    END DO

    NULLIFY(list)
    
    status = 0
    
  END SUBROUTINE clean_dimvar_list
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_dimvar_all(status, list)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED 

    ! I/O
    INTEGER,             INTENT(OUT) :: status
    TYPE(t_dimvar_list), POINTER     :: list

    ! LOCAL
    TYPE(t_dimvar_list), POINTER     :: ai => NULL()

    ! EMPTY LIST ?
    IF (.NOT. ASSOCIATED(list)) THEN

       WRITE(*,*) '   |    |->  *** DIMENSION VARIABLE LIST IS EMPTY ***'

    ELSE

       ai => list
       DO
          IF (.NOT. ASSOCIATED(ai)) EXIT
          !
          CALL write_dimvar_one(status, ai%this)
          IF (status /= 0) RETURN
          !
          ai => ai%next
       END DO
       
    END IF
    
    status = 0

  END SUBROUTINE write_dimvar_all
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_dimvar_one(status, dimvar)

    USE messy_main_channel_attributes, ONLY: write_attribute

    IMPLICIT NONE

    INTRINSIC :: MIN, MOD, SIZE, TRIM

    ! I/O
    INTEGER,        INTENT(OUT) :: status
    TYPE(t_dimvar), INTENT(IN)  :: dimvar

    ! LOCAL
    INTEGER, PARAMETER          :: vpl = 2  ! values per line
    INTEGER                     :: n1, ne
    INTEGER                     :: i        ! line counter
    INTEGER                     :: len

    WRITE(*,*) '   |    |-> ',TRIM(dimvar%name)

    WRITE(*,*) '   |    |   VALUES    : '
    len = SIZE(dimvar%val)
    n1 = 1
    ne = MIN(vpl, len)
    DO i=1, (len / vpl) + MIN(MOD(len, vpl), 1)
       WRITE(*,*) '   |    |    |-> ',dimvar%val(n1:ne)
       n1 = ne + 1
       ne = MIN((i+1)*vpl, len)
    END DO

    WRITE(*,*) '   |    |   ATTRIBUTES: '
    CALL write_attribute(status, dimvar%att)
 
  END SUBROUTINE write_dimvar_one
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE loc_dimvar(status, list, name, dimvar)

    IMPLICIT NONE
    
    INTRINSIC :: ADJUSTL, ASSOCIATED, LEN_TRIM, TRIM

    ! I/O
    INTEGER,                INTENT(OUT)          :: status
    TYPE(t_dimvar_list),    POINTER              :: list
    CHARACTER(LEN=*),       INTENT(IN)           :: name
    TYPE(t_dimvar),         POINTER              :: dimvar

    ! LOCAL
    TYPE(t_dimvar_list), POINTER :: ai  => NULL()
    TYPE(t_dimvar_list), POINTER :: ae  => NULL()
    LOGICAL                      :: lex

    ! INIT
    lex = .FALSE.

    ! CHECKS
    IF (LEN_TRIM(ADJUSTL(name)) > STRLEN_MEDIUM) THEN
       status = 951  ! DIMENSION VARIABLE NAME TOO LONG
       RETURN
    END IF
    !
    IF (.NOT. ASSOCIATED(list)) THEN
       status = 953  ! DIMENSION VARIABLE DOES NOT EXIST
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
       dimvar => ai%this
    ELSE
       status = 953  ! DIMENSION VARIABLE DOES NOT EXIST
       RETURN
    END IF

    status = 0

  END SUBROUTINE loc_dimvar
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_main_channel_dimvar
! **********************************************************************
