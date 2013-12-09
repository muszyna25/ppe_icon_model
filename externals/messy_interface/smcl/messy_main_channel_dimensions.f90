! **********************************************************************
MODULE messy_main_channel_dimensions
! **********************************************************************

  ! MESSY DATA TRANSFER AND EXPORT INTERFACE (MEMORY MANAGEMENT)
  !
  ! Author: Patrick Joeckel, MPICH, May 2005

  USE messy_main_constants_mem,     ONLY: DP, STRLEN_MEDIUM
  USE messy_main_channel_dimvar,    ONLY: t_dimvar_list

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  INTEGER, PARAMETER, PUBLIC :: DIMID_UNDEF = -1

  TYPE t_dimension
     CHARACTER(LEN=STRLEN_MEDIUM) :: name  = ''          ! NAME OF DIMENSION
     INTEGER                      :: id    = DIMID_UNDEF ! ID OF DIMENSION
     INTEGER                      :: len   = 0           ! LENGTH OF DIMENSION
     LOGICAL                      :: ltime = .FALSE.     ! FLAG FOR TIME DIM.
     !
     ! COORDINATE VARIABLE(S)
     TYPE(t_dimvar_list), POINTER :: var   => NULL()
  END TYPE t_dimension

  TYPE t_dimension_ptr
     TYPE(t_dimension), POINTER :: ptr => NULL()
  END TYPE t_dimension_ptr

  TYPE t_dimension_list
     PRIVATE
     TYPE(t_dimension)               :: this
     TYPE(t_dimension_list), POINTER :: next => NULL()
  END TYPE t_dimension_list

  ! ====================================================================
  ! CONCAT. LIST OF DIMENSIONS
  TYPE(t_dimension_list), POINTER, SAVE         :: GDIMLIST => NULL()
  ! NUMBER OF DIMENSIONS
  INTEGER,                         SAVE, PUBLIC :: NDIM = 0
  ! ====================================================================
!!$! um_ak_20090729+
!!$  ! DEFINE SPECIAL DIMENSIONS FOR AXIS RECOGNITION
!!$  INTEGER, SAVE, PUBLIC :: DIMID_LON     ! LONGITUDE ID
!!$  INTEGER, SAVE, PUBLIC :: DIMID_LAT     ! LATITUDE  ID
!!$  INTEGER, SAVE, PUBLIC :: DIMID_LEV     ! LEVEL     ID
!!$  INTEGER, SAVE, PUBLIC :: DIMID_ILEV    ! INTERFACE LEVEL ID
!!$! um_ak_20090729-
  ! ====================================================================

  ! TYPES
  PUBLIC :: t_dimension                  ! TYPE: dimension
  PUBLIC :: t_dimension_ptr              ! TYPE: dimension_ptr

  ! INTERFACES
  PUBLIC :: new_dimension                ! add dimension to (internal) list
  INTERFACE add_dimension_variable
     MODULE PROCEDURE add_dimvar_by_name
     MODULE PROCEDURE add_dimvar_by_id
  END INTERFACE
  PUBLIC :: add_dimension_variable
  PUBLIC :: update_dimension_variable    ! update values of dimension variable
                                         ! with ltime - flag
  INTERFACE add_dimension_variable_att
     MODULE PROCEDURE add_dimvar_att_by_name
     MODULE PROCEDURE add_dimvar_att_by_id
  END INTERFACE
  PUBLIC :: add_dimension_variable_att   ! add attribute to dimension variable
  INTERFACE get_dimension
     MODULE PROCEDURE get_dimension_by_name
     MODULE PROCEDURE get_dimension_by_id
  END INTERFACE
  PUBLIC :: get_dimension
! mz_ak_10091015+
  INTERFACE get_dimension_info
     MODULE PROCEDURE get_dimension_info_by_id   ! op_pj_20100121
     MODULE PROCEDURE get_dimension_info_by_name
  END INTERFACE
  PUBLIC :: get_dimension_info
! mz_ak_10091015+
  INTERFACE write_dimension
     MODULE PROCEDURE write_dimension_one
     MODULE PROCEDURE write_dimension_all
  END INTERFACE
  PUBLIC :: write_dimension              ! write dimension (list) to
  PUBLIC :: clean_dimensions             ! delete dimension list
  INTERFACE loc_dimension
     MODULE PROCEDURE loc_dimension_by_name
     MODULE PROCEDURE loc_dimension_by_id
  END INTERFACE
  !PRIVATE :: loc_dimension

CONTAINS

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE new_dimension(status, id, name, len, ltime)
    
    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, PRESENT, TRIM

    ! I/O
    INTEGER,                INTENT(OUT)          :: status
    INTEGER,                INTENT(OUT)          :: id
    CHARACTER(LEN=*),       INTENT(IN)           :: name
    INTEGER,                INTENT(IN)           :: len
    LOGICAL,                INTENT(IN), OPTIONAL :: ltime

    ! LOCAL
    TYPE(t_dimension_list), POINTER :: ai  => NULL()
    TYPE(t_dimension_list), POINTER :: ae  => NULL()
    TYPE(t_dimension),      POINTER :: dim => NULL()
    INTEGER                         :: zstat

    ! INIT
    id = 0

    CALL loc_dimension(zstat, GDIMLIST, name, dim)
    IF (zstat /= 905) THEN  ! DIMENSION DOES NOT EXIST (IS OK HERE !)
       IF (zstat == 0) THEN ! DIMENSION EXISTS ALREADY
          status = 904      ! DIMENSION EXISTS ALREADY
       ELSE
          status = zstat    ! ERROR
       END IF
       RETURN
    END IF

    ! GOTO END OF LIST
    ai => GDIMLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       ae => ai
       ai => ai%next
    END DO

    ! ADD NEW
    ALLOCATE(ai)
    NULLIFY(ai%next)
    IF (.NOT. ASSOCIATED(GDIMLIST)) THEN
       GDIMLIST => ai                  ! SET POINTER TO FIRST ELEMENT
    ELSE
       ae%next => ai                   ! SET NEXT POINTER OF LAST ELEMENT
       !                               ! TO NEW ELEMENT
    END IF

    ! SET VALUES
    ai%this%name = TRIM(ADJUSTL(name))
    ai%this%len  = len
    !
    IF (PRESENT(ltime)) THEN
       ai%this%ltime = ltime
    ELSE
       ai%this%ltime = .FALSE.  ! DEFAULT
    END IF
    !
    IF (ai%this%ltime) THEN
       ai%this%len  = 1         ! INTERNAL LENGHT FOR TIME ALWAS 1 !!!
    ELSE
       ai%this%len  = len
    END IF

    ! COUNT AND SET ID
    NDIM       = NDIM + 1
    ai%this%id = NDIM
    id         = NDIM

    status = 0

  END SUBROUTINE new_dimension
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE add_dimvar_by_name(status, dname, vname, val)

    USE messy_main_channel_dimvar, ONLY: add_dimvar

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    CHARACTER(LEN=*),       INTENT(IN)  :: dname   ! DIMENSION NAME
    CHARACTER(LEN=*),       INTENT(IN)  :: vname   ! VARIABLE NAME
    REAL(DP), DIMENSION(:), INTENT(IN)  :: val     ! VALUES

    ! LOCAL
    TYPE(t_dimension),      POINTER     :: dim

    CALL loc_dimension_by_name(status, GDIMLIST, dname, dim)
    IF (status /= 0) RETURN

    IF (SIZE(val) /= dim%len) THEN
       status = 954 ! VECTOR LENGTH NOT CONFORM WITH DIMENSION
       RETURN
    ENDIF

    CALL add_dimvar(status, dim%var, vname, val)

    status = 0
    
  END SUBROUTINE add_dimvar_by_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE add_dimvar_by_id(status, id, vname, val)

    USE messy_main_channel_dimvar, ONLY: add_dimvar

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    INTEGER,                INTENT(IN)  :: id      ! DIMENSION ID
    CHARACTER(LEN=*),       INTENT(IN)  :: vname   ! VARIABLE NAME
    REAL(DP), DIMENSION(:), INTENT(IN)  :: val     ! VALUES

    ! LOCAL
    TYPE(t_dimension),      POINTER     :: dim

    CALL loc_dimension_by_id(status, GDIMLIST, id, dim)
    IF (status /= 0) RETURN

    IF (SIZE(val) /= dim%len) THEN
       status = 954 ! VECTOR LENGTH NOT CONFORM WITH DIMENSION
       RETURN
    ENDIF

    CALL add_dimvar(status, dim%var, vname, val)

    status = 0
    
  END SUBROUTINE add_dimvar_by_id
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE update_dimension_variable(status, dname, vname, val)

    USE messy_main_channel_dimvar, ONLY: get_dimvar, t_dimvar

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    CHARACTER(LEN=*),       INTENT(IN)  :: dname   ! DIMENSION NAME
    CHARACTER(LEN=*),       INTENT(IN)  :: vname   ! VARIABLE NAME
    REAL(DP), DIMENSION(:), INTENT(IN)  :: val     ! VALUES

    ! LOCAL
    TYPE(t_dimension),      POINTER     :: dim
    TYPE(t_dimvar),         POINTER     :: dimvar

    CALL loc_dimension_by_name(status, GDIMLIST, dname, dim)
    IF (status /= 0) RETURN

    IF (.NOT. dim%ltime) THEN
       status = 955 ! UPADATE ONLY FOR TIME DIMENSION VARIABLES
       RETURN
    END IF

    IF (SIZE(val) /= dim%len) THEN    ! == 1
       status = 954 ! VECTOR LENGTH NOT CONFORM WITH DIMENSION
       RETURN
    ENDIF

    CALL get_dimvar(status, dim%var, vname, dimvar)
    IF (status /= 0) RETURN

    dimvar%val(:) = val(:)

    status = 0
    
  END SUBROUTINE update_dimension_variable
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE add_dimvar_att_by_name(status, dname, vname, aname, i, c, r &
       , loverwrite, iflag)

    USE messy_main_channel_dimvar,     ONLY: add_dimvar_att

    IMPLICIT NONE

    ! I/O
    INTEGER,                   INTENT(OUT)          :: status
    CHARACTER(LEN=*),          INTENT(IN)           :: dname
    CHARACTER(LEN=*),          INTENT(IN)           :: vname
    CHARACTER(LEN=*),          INTENT(IN)           :: aname
    INTEGER,                   INTENT(IN), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(IN), OPTIONAL :: c
    REAL(DP),                  INTENT(IN), OPTIONAL :: r
    LOGICAL,                   INTENT(IN), OPTIONAL :: loverwrite
    INTEGER,                   INTENT(IN), OPTIONAL :: iflag

    ! LOCAL
    TYPE(t_dimension), POINTER :: dim

    CALL loc_dimension_by_name(status, GDIMLIST, dname, dim)
    IF (status /= 0) RETURN

    CALL add_dimvar_att(status, dim%var, vname, aname, i, c, r &
         , loverwrite, iflag)
    IF (status /= 0) RETURN

  END SUBROUTINE add_dimvar_att_by_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE add_dimvar_att_by_id(status, id, vname, aname, i, c, r &
       , loverwrite, iflag)

    USE messy_main_channel_dimvar,     ONLY: add_dimvar_att

    IMPLICIT NONE

    ! I/O
    INTEGER,                   INTENT(OUT)          :: status
    INTEGER,                   INTENT(IN)           :: id
    CHARACTER(LEN=*),          INTENT(IN)           :: vname
    CHARACTER(LEN=*),          INTENT(IN)           :: aname
    INTEGER,                   INTENT(IN), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(IN), OPTIONAL :: c
    REAL(DP),                  INTENT(IN), OPTIONAL :: r
    LOGICAL,                   INTENT(IN), OPTIONAL :: loverwrite
    INTEGER,                   INTENT(IN), OPTIONAL :: iflag

    ! LOCAL
    TYPE(t_dimension), POINTER :: dim

    CALL loc_dimension_by_id(status, GDIMLIST, id, dim)
    IF (status /= 0) RETURN

    CALL add_dimvar_att(status, dim%var, vname, aname, i, c, r &
         , loverwrite, iflag)
    IF (status /= 0) RETURN

  END SUBROUTINE add_dimvar_att_by_id
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_dimension_by_name(status, name, dim)

    ! I/O
    INTEGER,                INTENT(OUT)          :: status
    CHARACTER(LEN=*),       INTENT(IN)           :: name
    TYPE(t_dimension),      POINTER              :: dim

    CALL loc_dimension_by_name(status, GDIMLIST, name, dim)

  END SUBROUTINE get_dimension_by_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_dimension_by_id(status, id, dim)

    ! I/O
    INTEGER,                INTENT(OUT)          :: status
    INTEGER,                INTENT(IN)           :: id
    TYPE(t_dimension),      POINTER              :: dim

    CALL loc_dimension_by_id(status, GDIMLIST, id, dim)

  END SUBROUTINE get_dimension_by_id
  ! -------------------------------------------------------------------
! um_ak_20091015+
  ! -------------------------------------------------------------------
  SUBROUTINE get_dimension_info_by_name(status, name, id, len)

    ! I/O
    INTEGER,                INTENT(OUT)           :: status
    CHARACTER(LEN=*),       INTENT(IN)            :: name
    INTEGER,                INTENT(OUT), OPTIONAL :: id
    INTEGER,                INTENT(OUT), OPTIONAL :: len
    ! LOCAL
    TYPE(t_dimension),      POINTER              :: dim

    CALL loc_dimension_by_name(status, GDIMLIST, name, dim)
    IF (status/=0) RETURN

    IF (PRESENT(id))  id  = dim%id
    IF (PRESENT(len)) len = dim%len

  END SUBROUTINE get_dimension_info_by_name
  ! -------------------------------------------------------------------
! um_ak_20091015-

! op_pj_20100121+
  ! -------------------------------------------------------------------
  SUBROUTINE get_dimension_info_by_id(status, dimid, name, len)

    ! I/O
    INTEGER,                INTENT(OUT)           :: status
    INTEGER,                INTENT(IN)            :: dimid
    CHARACTER(LEN=*),       INTENT(OUT), OPTIONAL :: name
    INTEGER,                INTENT(OUT), OPTIONAL :: len
    ! LOCAL
    TYPE(t_dimension),      POINTER              :: dim

    CALL loc_dimension_by_id(status, GDIMLIST, dimid, dim)
    IF (status/=0) RETURN

    IF (PRESENT(name))  name = dim%name
    IF (PRESENT(len))   len  = dim%len

  END SUBROUTINE get_dimension_info_by_id
  ! -------------------------------------------------------------------
! op_pj_20100121-

  ! -------------------------------------------------------------------
  SUBROUTINE clean_dimensions(status)

    USE messy_main_channel_dimvar, ONLY: clean_dimvar_list

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                INTENT(OUT) :: status

    ! LOCAL
    TYPE(t_dimension_list), POINTER :: ai => NULL()
    TYPE(t_dimension_list), POINTER :: ae => NULL()
    
    ai => GDIMLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       !
       ae => ai
       ai => ai%next
       !
       CALL clean_dimvar_list(status, ae%this%var)
       IF (status /= 0) RETURN
       !
       DEALLOCATE(ae)
       NULLIFY(ae)
       !
       ! COUNT
       NDIM = NDIM - 1
    END DO    

    NULLIFY(GDIMLIST)

    IF (NDIM /= 0) THEN
       status = 1002
       RETURN
    END IF

    status = 0

  END SUBROUTINE clean_dimensions
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_dimension_all(status)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                INTENT(OUT) :: status

    ! LOCAL
    TYPE(t_dimension_list), POINTER :: ai => NULL()

    WRITE(*,*) '=== DIMENSIONS: ===================================='

    ! EMPTY LIST ?
    IF (.NOT. ASSOCIATED(GDIMLIST)) THEN

       WRITE(*,*) '   *** DIMENSION LIST IS EMPTY ***'

    ELSE

       ai => GDIMLIST
       DO
          IF (.NOT. ASSOCIATED(ai)) EXIT
          !
          CALL write_dimension_one(status, ai%this)
          IF (status /= 0) RETURN
          !
          ai => ai%next
       END DO
       
    END IF
    
    WRITE(*,*) '===================================================='

    status = 0

  END SUBROUTINE write_dimension_all
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_dimension_one(status, dim)

    USE messy_main_channel_dimvar, ONLY: write_dimvar

    IMPLICIT NONE

    ! I/O
    INTEGER,           INTENT(OUT) :: status
    TYPE(t_dimension), INTENT(IN)  :: dim


    WRITE(*,*) '   |-> NAME      : ',dim%name
    WRITE(*,*) '   |   ID        : ',dim%id
    WRITE(*,*) '   |   LENGTH    : ',dim%len
    WRITE(*,*) '   |   IS_TIME ? : ',dim%ltime
    WRITE(*,*) '   |   VARIABLES : '
    CALL write_dimvar(status, dim%var)

  END SUBROUTINE write_dimension_one
  ! -------------------------------------------------------------------

  ! *******************************************************************
  ! PRIVATE ROUTINES
  ! *******************************************************************

  ! -------------------------------------------------------------------
  SUBROUTINE loc_dimension_by_name(status, list, name, dim)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, LEN_TRIM, TRIM 

    ! I/O
    INTEGER,                INTENT(OUT)          :: status
    TYPE(t_dimension_list), POINTER              :: list
    CHARACTER(LEN=*),       INTENT(IN)           :: name
    TYPE(t_dimension),      POINTER              :: dim

    ! LOCAL
    TYPE(t_dimension_list), POINTER :: ai  => NULL()
    TYPE(t_dimension_list), POINTER :: ae  => NULL()
    LOGICAL                         :: lex

    ! INIT
    lex = .FALSE.
    NULLIFY(dim)

    ! CHECKS
    IF (LEN_TRIM(ADJUSTL(name)) > STRLEN_MEDIUM) THEN
       status = 902  ! DIMENSION NAME TOO LONG
       RETURN
    END IF
    !
    IF (.NOT. ASSOCIATED(list)) THEN
       status = 905  ! DIMENSION DOES NOT EXIST
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
       dim => ai%this
    ELSE
       status = 905  ! DIMENSION DOES NOT EXIST
       RETURN
    END IF

    status = 0

  END SUBROUTINE loc_dimension_by_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE loc_dimension_by_id(status, list, id, dim)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                INTENT(OUT)          :: status
    TYPE(t_dimension_list), POINTER              :: list
    INTEGER,                INTENT(IN)           :: id
    TYPE(t_dimension),      POINTER              :: dim

    ! LOCAL
    TYPE(t_dimension_list), POINTER :: ai  => NULL()
    TYPE(t_dimension_list), POINTER :: ae  => NULL()
    LOGICAL                         :: lex

    ! INIT
    lex = .FALSE.
    NULLIFY(dim)

    ! CHECKS
    IF (id <=0) THEN
       status = 910 ! DIMENSION ID INVALID
       RETURN
    END IF
    !
    IF (.NOT. ASSOCIATED(list)) THEN
       status = 905  ! DIMENSION DOES NOT EXIST
       RETURN
    END IF

    ! CHECK, IF IT EXISTS
    ai => list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       IF (id == ai%this%id) THEN
          lex = .TRUE.
          EXIT
       END IF
       ae => ai
       ai => ai%next
    END DO    

    IF (lex) THEN
       dim => ai%this
    ELSE
       status = 905  ! DOES NOT EXIST
       RETURN
    END IF

    status = 0

  END SUBROUTINE loc_dimension_by_id
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_main_channel_dimensions
! **********************************************************************
