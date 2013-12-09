! **********************************************************************
MODULE messy_main_tracer
! **********************************************************************

  ! THIS MODULE IS THE MAIN SMCL MODULE OF THE GENERIC MESSy-SUBMODEL
  ! 'TRACER'

  USE messy_main_constants_mem, ONLY: DP &
       , STRLEN_MEDIUM, STRLEN_LONG, STRLEN_VLONG

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE

  CHARACTER(len=*), PARAMETER, PUBLIC :: modstr = 'tracer'
  CHARACTER(len=*), PARAMETER, PUBLIC :: modver = '2.3'

  PUBLIC :: DP

  ! GLOBAL NAMELIST SWITCHES (CTRL)
  LOGICAL, SAVE, PUBLIC :: l_family     = .FALSE.
  LOGICAL, SAVE, PUBLIC :: l_pdef       = .FALSE.

  ! SPECIAL STRING LENGTHs
  INTEGER, PARAMETER, PUBLIC  :: STRLEN_TRSET  = 10
  INTEGER, PARAMETER, PUBLIC  :: STRLEN_FNAME  = 2*STRLEN_MEDIUM + 1         

  ! GENERAL SWITCHES
  INTEGER, PARAMETER, PUBLIC :: UNDEFINED = -1
  INTEGER, PARAMETER, PUBLIC :: OFF       =  0
  INTEGER, PARAMETER, PUBLIC :: ON        =  1
  ! ... AEROSOL MODEL METHOD
  INTEGER, PARAMETER, PUBLIC :: MODAL     = 2
  INTEGER, PARAMETER, PUBLIC :: BIN       = 3

  ! TYPE
  INTEGER, PARAMETER, PUBLIC :: SINGLE  = 0
  INTEGER, PARAMETER, PUBLIC :: FAMILY  = 1
  INTEGER, PARAMETER, PUBLIC :: ISOTOPE = 2
  INTEGER, PARAMETER         :: MAX_TYPE = 2

  ! 'HOSTING' MEDIUM
  INTEGER, PARAMETER, PUBLIC :: AIR        = 1
  INTEGER, PARAMETER, PUBLIC :: AEROSOL    = 2
  INTEGER, PARAMETER, PUBLIC :: CLOUD      = 3
  INTEGER, PARAMETER, PUBLIC :: OCEAN      = 4
  INTEGER, PARAMETER, PUBLIC :: LAKE       = 5
  INTEGER, PARAMETER, PUBLIC :: RIVER      = 6
  INTEGER, PARAMETER, PUBLIC :: LANDICE    = 7
  INTEGER, PARAMETER, PUBLIC :: SEAICE     = 8
  INTEGER, PARAMETER, PUBLIC :: VEGETATION = 9
  INTEGER, PARAMETER         :: MAX_MEDIUM = 9

  ! QUANTITY
  INTEGER, PARAMETER, PUBLIC :: AMOUNTFRACTION = 1 ! = MOLAR MIXING RATIO
  INTEGER, PARAMETER, PUBLIC :: NUMBERDENSITY  = 2
  INTEGER, PARAMETER, PUBLIC :: CONCENTRATION  = 3
  INTEGER, PARAMETER         :: MAX_QUANTITY   = 3

  ! =================================================================

  ! INTEGER CONTAINERS
  INTEGER, PARAMETER, PUBLIC :: I_ADVECT     =  1  ! ADVECTION
  INTEGER, PARAMETER, PUBLIC :: I_CONVECT    =  2  ! CONVECTION
  INTEGER, PARAMETER, PUBLIC :: I_VDIFF      =  3  ! VERTICAL DIFFUSION
  INTEGER, PARAMETER, PUBLIC :: I_WETDEP     =  4  ! WET DEPOSITION
  INTEGER, PARAMETER, PUBLIC :: I_DRYDEP     =  5  ! DRY DEPOSITION
  INTEGER, PARAMETER, PUBLIC :: I_SEDI       =  6  ! SEDIMENTATION
  INTEGER, PARAMETER, PUBLIC :: I_SCAV       =  7  ! SCAVENGING
  INTEGER, PARAMETER, PUBLIC :: I_MIX        =  8  ! TURBULENT MIXING
  INTEGER, PARAMETER, PUBLIC :: I_FORCE_COL  =  9  ! FORCING IN COLUMN MODE
  INTEGER, PARAMETER, PUBLIC :: I_INTEGRATE  = 10  ! TIME INTEGRATION
  INTEGER, PARAMETER, PUBLIC :: I_TIMEFILTER = 11  ! TIME FILTER
  INTEGER, PARAMETER, PUBLIC :: I_FORCE_INIT = 12  ! FORCE INIT AFTER RESTART
  INTEGER, PARAMETER, PUBLIC :: I_AEROSOL_METHOD = 13 ! MODAL OR BIN
  INTEGER, PARAMETER, PUBLIC :: I_AEROSOL_MODE   = 14 ! MODE OR BIN NUMBER
  INTEGER, PARAMETER, PUBLIC :: I_AEROSOL_SOL    = 15 ! SOLUBLE ON/OFF
  INTEGER, PARAMETER, PUBLIC :: I_HDIFF          = 16 ! HORIZONTAL DIFFUSION
  INTEGER, PARAMETER, PUBLIC :: I_RELAX          = 17 ! BOUNDARY DATA AVAILABLE
                                                     ! i.e. relaxation possible
  INTEGER, PARAMETER, PUBLIC :: I_MMD_INIT       = 18
! op_pj_20100319+
  INTEGER, PARAMETER, PUBLIC :: I_TAG_REG_IDT    = 19 ! id of associated regular
  !                                                   ! species
  INTEGER, PARAMETER, PUBLIC :: I_TAG_SPECIFIC   = 20 ! flag for special
  !                                                   ! treatment
! op_pj_20100319-
  INTEGER, PARAMETER, PUBLIC :: MAX_CASK_I       = 20
  !
  INTEGER, PARAMETER :: NAMES_CASK_I_STRLEN = 14
  CHARACTER(LEN=NAMES_CASK_I_STRLEN), DIMENSION(MAX_CASK_I) &
       , PARAMETER, PUBLIC :: &
       NAMES_CASK_I = (/ &
       'advect        ', 'convect       ', 'vdiff         ' , &
       'wetdep        ', 'drydep        ', 'sedi          ' , &
       'scav          ', 'mix           ', 'force_col     ' , &
       'integrate     ', 'timefilter    ', 'force_init    ' , &
       'aerosol_method', 'aerosol_mode  ', 'aerosol_sol   ' , &
       'hori_diff     ', 'relaxation    ', 'mmd_init      ' , &
       'tag_reg_idt   ', 'tag_specific  ' /)
  ! NOTES: - CASK_I(I_TAG_REG_IDT) will default to TRACER-ID
  !          (see subroutine new_tracer)
  INTEGER, DIMENSION(MAX_CASK_I), PARAMETER, PUBLIC :: DEFAULT_CASK_I = &
       (/ ON, ON, ON, OFF, OFF, OFF, OFF, ON, OFF, ON, ON, OFF &
       , MODAL, 0, ON, OFF, ON, OFF, 0, 0/)

  ! =================================================================

  ! STRING CONTAINERS
  INTEGER, PARAMETER, PUBLIC :: S_AEROSOL_MODEL = 1
  INTEGER, PARAMETER, PUBLIC :: MAX_CASK_S      = 1
  !
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(MAX_CASK_S), PARAMETER, PUBLIC :: &
       DEFAULT_CASK_S = (/ '' /)
  !
  INTEGER, PARAMETER :: NAMES_CASK_S_STRLEN = 14
  CHARACTER(LEN=NAMES_CASK_S_STRLEN), DIMENSION(MAX_CASK_S) &
       , PARAMETER, PUBLIC :: &
       NAMES_CASK_S = (/ &
       'aerosol_model ' &
       /)  

  ! =================================================================

  ! REAL CONTAINERS
  INTEGER, PARAMETER, PUBLIC :: R_MOLARMASS       = 1
  INTEGER, PARAMETER, PUBLIC :: R_HENRY           = 2
  INTEGER, PARAMETER, PUBLIC :: R_DRYREAC_SF      = 3
  INTEGER, PARAMETER, PUBLIC :: R_VINI            = 4
  INTEGER, PARAMETER, PUBLIC :: R_AEROSOL_DENSITY = 5
  INTEGER, PARAMETER, PUBLIC :: MAX_CASK_R        = 5
  !
  REAL(DP), DIMENSION(MAX_CASK_R), PARAMETER, PUBLIC :: &
       DEFAULT_CASK_R = (/ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /)
  !
  INTEGER, PARAMETER :: NAMES_CASK_R_STRLEN = 15
  CHARACTER(LEN=NAMES_CASK_R_STRLEN), DIMENSION(MAX_CASK_R) &
       , PARAMETER, PUBLIC :: &
       NAMES_CASK_R = (/ &
       'molarmass      ', 'henry          ', 'dryreac_sf     ', &
       'vini           ', 'aerosol_density' &
       /)  

  ! SPECIAL ERROR NUMBERS
  INTEGER, PARAMETER, PUBLIC :: TR_EXIST  = 202  ! TRACER EXISTS
  INTEGER, PARAMETER, PUBLIC :: TR_NEXIST = 205  ! TRACER DOES NOT EXIST

  ! PUBLIC STRUCTURES
  PUBLIC :: t_ident         ! TRACER IDENTIFICATION
  PUBLIC :: t_meta          ! ADDITIONAL META INFORMATION
  !
  PUBLIC :: t_trinfo        ! META-STRUCT WITH ALL INFORMATION
  PUBLIC :: t_trinfo_tp     ! TRACER PROPERTIES
  PUBLIC :: t_trinfo_list   ! LIST OF META-STRUCT WITH ALL INFORMATION

  ! PUBLIC SUBROUTINES
  PUBLIC :: new_tracer_set   ! BML, BMIL     ! DEFINE NEW TRACER SET
  PUBLIC :: copy_tracer_set  ! BML, BMIL     ! COPY COMPLETE TRACER SET
  PUBLIC :: setup_tracer_set ! BML, BMIL     ! ALLOCATE MEMORY FOR TRACER SET
  PUBLIC :: get_tracer_set   ! BML, BMIL     ! SET REFERENCES TO TRACER SETS
  PUBLIC :: print_tracer_set ! BML, BMIL     ! PRINT TRACER SET SUMMARY
  PUBLIC :: print_tracer_set_val ! BML, BMIL ! PRINT TRACER VALUE RANGE
  PUBLIC :: clean_tracer_set ! BML, BMIL     ! REMOVE TRACER SET FROM MEMORY
  !
  PUBLIC :: new_tracer       ! SMIL          ! DEFINE NEW TRACER IN SET
  PUBLIC :: set_tracer       ! SMIL          ! DEFINE TRACER PROPERTIES
! mz_pj_20070817+ will become obsolete
  PUBLIC :: new_tracer_old   ! SMIL          ! DEFINE NEW TRACER IN SET
! mz_pj_20070817-
! mz_pj_20090507+
  INTERFACE get_tracer
     MODULE PROCEDURE get_tracer_by_name
     MODULE PROCEDURE get_tracer_by_id
     ! op_pj_20100325+
     MODULE PROCEDURE get_tracer_i
     MODULE PROCEDURE get_tracer_r
     MODULE PROCEDURE get_tracer_s
     ! op_pj_20100325-
  END INTERFACE
! mz_pj_20090507-
  PUBLIC :: get_tracer       ! SMIL          ! GET TRACER INFORMATION FROM SET
  PUBLIC :: get_tracer_list  ! SMIL          ! GET TRACERS WITH SAME BASENAME
  PUBLIC :: tracer_iniflag                   ! SET INITIALISATION FLAG
  !
  PUBLIC :: tracer_error_str                 ! RETURN STATUS INFORMATION
  PUBLIC :: param2string                     ! PARAMETER TO STRING CONVERSION
  PUBLIC :: full2base_sub                    ! fullname -> basename + subname
  PUBLIC :: main_tracer_read_nml_ctrl

  INTERFACE get_tracer_set
     MODULE PROCEDURE get_tracer_set_by_name
     MODULE PROCEDURE get_tracer_set_by_id
  END INTERFACE

  INTERFACE set_tracer
     MODULE PROCEDURE set_tracer_i
     MODULE PROCEDURE set_tracer_r
     MODULE PROCEDURE set_tracer_s
  END INTERFACE

  ! PRIVATE SUBROUTINES
  !PRIVATE :: get_tracer_set_id  ! get tracer set id

  ! STRUCTURE DEFINITIONS
  TYPE t_ident    ! IDENTIFICATION
     CHARACTER(LEN=STRLEN_MEDIUM)     :: basename    = '' ! name of tracer
     CHARACTER(LEN=STRLEN_MEDIUM)     :: subname     = '' ! OPTIONAL subname
     CHARACTER(LEN=STRLEN_FNAME)      :: fullname    = '' ! name_subname
     CHARACTER(LEN=STRLEN_LONG)       :: longname    = '' ! 
     CHARACTER(LEN=STRLEN_MEDIUM)     :: unit        = '' !
     CHARACTER(LEN=STRLEN_MEDIUM)     :: submodel    = '' ! requesting submodel
     INTEGER                          :: idx              ! tracer index in set
     INTEGER                          :: medium   = AIR   ! hosting medium
     INTEGER                          :: quantity = AMOUNTFRACTION
     INTEGER                          :: type     = SINGLE ! 
  END TYPE t_ident

  TYPE t_meta
     INTEGER,                      DIMENSION(MAX_CASK_I)  :: cask_i &
          = DEFAULT_CASK_I
     REAL(dp),                     DIMENSION(MAX_CASK_R)  :: cask_r &
          = DEFAULT_CASK_R
     CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(MAX_CASK_S)  :: cask_s &
          = DEFAULT_CASK_S
     LOGICAL :: linit = .FALSE.
  END TYPE t_meta

  TYPE t_trinfo
     TYPE(t_ident) :: ident       ! IDENTIFICATION
     TYPE(t_meta)  :: meta        ! ADDITIONAL META-INFORMATION
  END TYPE t_trinfo

  TYPE t_trinfo_list
     TYPE(t_trinfo)               :: info
     TYPE(t_trinfo_list), POINTER :: next
  END TYPE t_trinfo_list

  TYPE t_trinfo_tp
     TYPE(t_trinfo), POINTER :: tp
  END TYPE t_trinfo_tp

  ! ################## INTERNAL MEMORY AND POINTER MANAGEMENT ############
  TYPE T_TRACERSET
     ! NAME OF THE TRACER-SET FOR IDENTIFICATION
     CHARACTER(LEN=STRLEN_TRSET)              :: name
     !
     ! SWITCH TO ENABLE/DISABLE THIS TRACER SET, E.G., IF ITS
     ! EXISTENCE IS DEPENDENT ON A SPECIFIC SUBMODEL
     LOGICAL                                  :: l_enable = .TRUE.
     !
     ! SWITCH FOR TRIGGERING THIS SET TO BE INITIALIZED BY THE
     ! main_tracer_init_tracer (BMIL) SUBROUTINE
     LOGICAL                                  :: l_init = .TRUE.
     !
     ! NUMBER OF TRACERS IN SET
     INTEGER                                  :: ntrac  = 0
     !
     ! CONCATENATED LIST OF TRACER META INFORAMTION STRUCTURES
     TYPE(t_trinfo_list),             POINTER :: tilist => NULL()
     !
     ! ENUMERATED LIST OF TRACER META INFORAMTION STRUCTURES
     ! NOTE: THIS IS ONLY AVAILABLE AFTER CALL setup_tracer_set
     TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti     => NULL()
     !
     ! MEMORY FOR TRACER DATA
     REAL(DP), DIMENSION(:,:,:,:,:,:),  POINTER :: mem    => NULL()
     ! NUMBER OF (TIME-) LEVELS
     INTEGER                                    :: nt     = 0
     ! 3 STANDARD (TIME-) LEVELS OF MEMORY FOR TRACER DATA
     LOGICAL                                  :: l_tfstd = .FALSE.
     REAL(DP), DIMENSION(:,:,:,:,:),  POINTER :: xt      => NULL()
     REAL(DP), DIMENSION(:,:,:,:,:),  POINTER :: xtte    => NULL()
     REAL(DP), DIMENSION(:,:,:,:,:),  POINTER :: xtm1    => NULL()
     ! 'EXTENDED' MEMORY
     INTEGER                                  :: nex     = 0
     REAL(DP), DIMENSION(:,:,:,:,:),  POINTER :: xmem    => NULL()
  END TYPE T_TRACERSET
  PUBLIC :: T_TRACERSET
  !
  INTEGER, PARAMETER, PUBLIC :: NMAXSETID = 10  ! MAX. NUMBER OF TRACER SETS
  INTEGER,      SAVE, PUBLIC :: NSETID    = 0   ! ACT. NUMBER OF TRACER SETS
  TYPE(T_TRACERSET), DIMENSION(NMAXSETID), PUBLIC, SAVE :: TRSET
  ! ################## INTERNAL MEMORY AND POINTER MANAGEMENT ############

CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE new_tracer_set(status, setname, l_enable)

    IMPLICIT NONE
    INTRINSIC :: LEN, TRIM

    ! I/O
    INTEGER,               INTENT(OUT) :: status
    CHARACTER(LEN=*),      INTENT(IN)  :: setname
    LOGICAL,               INTENT(IN)  :: l_enable

    ! LOCAL
    INTEGER :: zid

    status = 100 ! ERROR: setname too long
    IF (LEN(setname) > STRLEN_TRSET) RETURN

    status = 101 ! ERROR: setname not unique
    DO zid = 1, NSETID
       IF (TRIM(trset(zid)%name) == TRIM(setname)) RETURN
    END DO

    status = 102 ! ERROR: no free tracer set available
    NSETID = NSETID + 1
    IF (NSETID > NMAXSETID) RETURN

    ! RESULT
    trset(NSETID)%name     = TRIM(setname)
    trset(NSETID)%l_enable = l_enable

    status = 0

  END SUBROUTINE new_tracer_set
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE copy_tracer_set(status, oldset, newset)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, LEN, TRIM

    ! I/O
    INTEGER,             INTENT(OUT) :: status
    CHARACTER(LEN=*),    INTENT(IN)  :: oldset
    CHARACTER(LEN=*),    INTENT(IN)  :: newset

    ! LOCAL
    INTEGER                      :: id_old, id_new
    TYPE(t_trinfo_list), POINTER :: ti_old => NULL()
    TYPE(t_trinfo_list), POINTER :: ti_new => NULL()
    TYPE(t_trinfo_list), POINTER :: te_new => NULL()

    status = 100 ! ERROR: setname too long
    IF (LEN(oldset) > STRLEN_TRSET) RETURN
    IF (LEN(newset) > STRLEN_TRSET) RETURN

    status = 101 ! ERROR: setname not unique
    DO id_new = 1, NSETID
       IF (TRIM(trset(id_new)%name) == TRIM(newset)) RETURN
    END DO

    status = 102 ! ERROR: no free tracer set available
    NSETID = NSETID + 1
    IF (NSETID > NMAXSETID) RETURN
    id_new = NSETID

    ! CHECK, IF OLD SET IS AVAILABLE; AND SET ID
    CALL get_tracer_set_id(status, oldset, id_old)
    IF (status /= 0) RETURN

    ! RESULT: NEW SET
    trset(id_new)%name = TRIM(newset)
    trset(id_new)%l_enable = trset(id_old)%l_enable

    ! COPY ALL TRACERS; LOOP OVER TRACERS IN OLD SET
    ti_old => trset(id_old)%tilist
    tracer_loop: DO 
       IF (.NOT. ASSOCIATED(ti_old)) EXIT

       ! ADD NEW TRACER TO NEW LIST
       ALLOCATE(ti_new)
       NULLIFY(ti_new%next)
       IF (trset(id_new)%ntrac == 0) THEN
          trset(id_new)%tilist => ti_new   ! SET POINTER TO FIRST ELEMENT
       ELSE
          te_new%next => ti_new            ! SET NEXT POINTER OF LAST ELEMENT
          !                                ! TO NEW ELEMENT
       END IF
       te_new => ti_new

       ! ADJUST TRACER INDEX
       trset(id_new)%ntrac = trset(id_new)%ntrac + 1

       ! SET TRACER INFORMATION
       ti_new%info = ti_old%info

       ti_old => ti_old%next
    END DO tracer_loop

    status = 0

  END SUBROUTINE copy_tracer_set
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE new_tracer(status, setname, basename, submodel             &
       , idx, subname, longname, unit, medium, quantity, type           &
       , cask_i, cask_r, cask_s)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, INDEX, LEN, PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    !
    CHARACTER(LEN=*), INTENT(IN)            :: basename
    CHARACTER(LEN=*), INTENT(IN)            :: submodel
    !
    INTEGER,          INTENT(OUT), OPTIONAL :: idx
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: subname
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: longname
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: unit
    INTEGER,          INTENT(IN),  OPTIONAL :: medium
    INTEGER,          INTENT(IN),  OPTIONAL :: quantity
    INTEGER,          INTENT(IN),  OPTIONAL :: type
    !
    INTEGER,  DIMENSION(MAX_CASK_I),   INTENT(IN),  OPTIONAL :: cask_i
    REAL(DP), DIMENSION(MAX_CASK_R),   INTENT(IN),  OPTIONAL :: cask_r
    CHARACTER(LEN=STRLEN_MEDIUM), &
         DIMENSION(MAX_CASK_S), INTENT(IN), OPTIONAL :: cask_s

    ! LOCAL
    INTEGER                        :: zid     = 0
    TYPE(t_trinfo_list), POINTER   :: ti      => NULL()
    TYPE(t_trinfo_list), POINTER   :: te      => NULL()
    CHARACTER(LEN=STRLEN_FNAME)    :: fullname

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! TRACER
    status = 200 ! ERROR: tracer basename too long
    IF (LEN(basename) > STRLEN_MEDIUM) RETURN

    status = 208 ! ERROR: basename MUST NOT contain '_'
    IF (INDEX(basename, '_') /= 0) RETURN

    IF (PRESENT(subname)) THEN
       status = 201 ! ERROR: tracer subname too long
       IF (LEN(subname) > STRLEN_MEDIUM) RETURN
       IF (TRIM(subname) /= '') THEN
          fullname = TRIM(basename)//'_'//TRIM(subname)
       ELSE
          fullname = TRIM(basename)
       END IF
    ELSE
       fullname = TRIM(basename)
    END IF

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       IF (PRESENT(idx)) idx = 0
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (TRIM(ti%info%ident%fullname) == TRIM(fullname)) THEN
          status = TR_EXIST  ! ERROR: tracer exists already
          RETURN
       END IF
       te => ti
       ti => ti%next
    END DO

    ! ADD NEW TRACER TO LIST
    ALLOCATE(ti)
    NULLIFY(ti%next)
    IF (trset(zid)%ntrac == 0) THEN
       trset(zid)%tilist => ti         ! SET POINTER TO FIRST ELEMENT
    ELSE
       te%next => ti                   ! SET NEXT POINTER OF LAST ELEMENT
       !                               ! TO NEW ELEMENT
    END IF

    ! ADJUST TRACER INDEX
    trset(zid)%ntrac = trset(zid)%ntrac + 1
    IF (PRESENT(idx)) idx = trset(zid)%ntrac

    ! SET TRACER INFORMATION
    ! - IDENTIFICATION
    ! -- MANDATORY
    ti%info%ident%basename                   = TRIM(basename)    ! BASENAME
    !
    IF (LEN(submodel) > STRLEN_MEDIUM) THEN
       status = 204 ! ERROR: unit too long
       RETURN
    END IF
    ti%info%ident%submodel                   = TRIM(submodel)    ! SUBMODEL
    !
    ti%info%ident%fullname                   = TRIM(fullname)    ! FULLNAME
    ti%info%ident%idx                        = trset(zid)%ntrac  ! SET INDEX

    ! -- OPTIONAL
    !
    IF (PRESENT(subname)) THEN           ! SUBNAME
       IF (LEN(subname) > STRLEN_MEDIUM) THEN
          status = 201 ! ERROR: subname too long
          RETURN
       END IF
       ti%info%ident%subname    = TRIM(subname)
    END IF
    !
    !                                    ! LONGNAME
    IF (PRESENT(longname)) THEN
       IF (LEN(longname) > STRLEN_LONG) THEN
          status = 207 ! ERROR: longname too long
          RETURN
       END IF
       ti%info%ident%longname  = TRIM(longname)
    END IF
    !
    IF (PRESENT(unit)) THEN              ! UNIT
       IF (LEN(unit) > STRLEN_MEDIUM) THEN
          status = 203 ! ERROR: unit too long
          RETURN
       END IF
       ti%info%ident%unit          = TRIM(unit)
    END IF
    !
    IF (PRESENT(medium)) THEN            ! MEDIUM
       IF ((medium < 0) .OR. (medium > MAX_MEDIUM) ) THEN
          status = 300  ! ERROR: unknown medium
          RETURN
       END IF
       ti%info%ident%medium   = medium
    END IF
    !
    IF (PRESENT(quantity)) THEN          ! QUANTITY
       IF ((quantity < 0) .OR. (quantity > MAX_QUANTITY) ) THEN
          status = 301  ! ERROR: unknown quantity
          RETURN
       END IF
       ti%info%ident%quantity = quantity
    END IF
    ! 
    IF (PRESENT(type)) THEN              ! TYPE
       IF ((type < 0) .OR. (type > MAX_TYPE  ) ) THEN
          status = 302  ! ERROR: unknown type
          RETURN
       END IF
       ti%info%ident%type = type
    END IF

    ! ASSIGN OPTIONAL PARAMETERS
    ! DEFAULT
    ti%info%meta%cask_i(:) = default_cask_i(:)
    ti%info%meta%cask_r(:) = default_cask_r(:)
    ti%info%meta%cask_s(:) = default_cask_s(:)
    ! SPECIFIC
    IF (PRESENT(cask_i)) ti%info%meta%cask_i(:) = cask_i(:)
    IF (PRESENT(cask_r)) ti%info%meta%cask_r(:) = cask_r(:)
    IF (PRESENT(cask_s)) ti%info%meta%cask_s(:) = cask_s(:)
    !
    ! special "dynamic default":
    IF (ti%info%meta%cask_i(I_TAG_REG_IDT) == 0) &
         ti%info%meta%cask_i(I_TAG_REG_IDT) = ti%info%ident%idx

    status = 0

  END SUBROUTINE new_tracer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_tracer_s(status, setname, idx, flag, s)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    INTEGER,          INTENT(IN)            :: flag
    CHARACTER(LEN=*), INTENT(IN)            :: s

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    IF ( (flag < 1) .OR. (flag > MAX_CASK_S) ) THEN
       status = 702 ! STRING FLAG INDEX OUT OF RANGE
       RETURN
    END IF
    ti%info%meta%cask_s(flag) = s

    status = 0

  END SUBROUTINE set_tracer_s
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_tracer_i(status, setname, idx, flag, i)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    INTEGER,          INTENT(IN)            :: flag
    INTEGER,          INTENT(IN)            :: i

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    IF ( (flag < 1) .OR. (flag > MAX_CASK_I) ) THEN
       status = 700 ! INTEGER FLAG INDEX OUT OF RANGE
       RETURN
    END IF
    ti%info%meta%cask_i(flag) = i

    status = 0

  END SUBROUTINE set_tracer_i
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_tracer_r(status, setname, idx, flag, r)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    INTEGER,          INTENT(IN)            :: flag
    REAL(DP),         INTENT(IN)            :: r

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    IF ( (flag < 1) .OR. (flag > MAX_CASK_R) ) THEN
       status = 701 ! REAL FLAG INDEX OUT OF RANGE
       RETURN
    END IF
    ti%info%meta%cask_r(flag) = r

    status = 0

  END SUBROUTINE set_tracer_r
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE new_tracer_old(status, setname, basename, submodel         &
       , idx, subname, longname, unit, medium, quantity, type           &
       , nadvect, nconvect, nvdiff, nwetdep, ndrydep, nsedi, nscav      &
       , nmix                                                           &
       , nintegrate, ntimefilter                                        &
       , vini, lforce_init                                              &
       , nforce_col                                                     &
       , molarmass, henry, dryreac_sf                                   &
       , m_aerosol_method, m_aerosol_modelname, m_aerosol_density       &
       , m_aerosol_mode, m_aerosol_nsol                                 &
       )

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, INDEX, LEN, PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    !
    CHARACTER(LEN=*), INTENT(IN)            :: basename
    CHARACTER(LEN=*), INTENT(IN)            :: submodel
    !
    INTEGER,          INTENT(OUT), OPTIONAL :: idx
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: subname
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: longname
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: unit
    INTEGER,          INTENT(IN),  OPTIONAL :: medium
    INTEGER,          INTENT(IN),  OPTIONAL :: quantity
    INTEGER,          INTENT(IN),  OPTIONAL :: type
    !
    INTEGER,          INTENT(IN),  OPTIONAL :: nadvect 
    INTEGER,          INTENT(IN),  OPTIONAL :: nconvect
    INTEGER,          INTENT(IN),  OPTIONAL :: nvdiff
    INTEGER,          INTENT(IN),  OPTIONAL :: nwetdep
    INTEGER,          INTENT(IN),  OPTIONAL :: ndrydep
    INTEGER,          INTENT(IN),  OPTIONAL :: nsedi
    INTEGER,          INTENT(IN),  OPTIONAL :: nscav
    INTEGER,          INTENT(IN),  OPTIONAL :: nmix
    ! SPECIAL FOR COLUMN MODE
    INTEGER,          INTENT(IN),  OPTIONAL :: nforce_col
    ! NUMERICAL
    INTEGER,          INTENT(IN),  OPTIONAL :: nintegrate
    INTEGER,          INTENT(IN),  OPTIONAL :: ntimefilter
    ! INFRASTRUCTURE
    REAL(DP),         INTENT(IN),  OPTIONAL :: vini
    LOGICAL,          INTENT(IN),  OPTIONAL :: lforce_init
    ! SINGLE
    REAL(DP),         INTENT(IN),  OPTIONAL :: molarmass
    REAL(DP),         INTENT(IN),  OPTIONAL :: henry
    REAL(DP),         INTENT(IN),  OPTIONAL :: dryreac_sf
    ! MEDIUM AEROSOL
    INTEGER,          INTENT(IN),  OPTIONAL :: m_aerosol_method
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: m_aerosol_modelname
    REAL(DP),         INTENT(IN),  OPTIONAL :: m_aerosol_density
    INTEGER,          INTENT(IN),  OPTIONAL :: m_aerosol_mode
    INTEGER,          INTENT(IN),  OPTIONAL :: m_aerosol_nsol

    ! LOCAL
    INTEGER                        :: zid     = 0
    TYPE(t_trinfo_list), POINTER   :: ti      => NULL()
    TYPE(t_trinfo_list), POINTER   :: te      => NULL()
    CHARACTER(LEN=STRLEN_FNAME)    :: fullname

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! TRACER
    status = 200 ! ERROR: tracer basename too long
    IF (LEN(basename) > STRLEN_MEDIUM) RETURN

    status = 208 ! ERROR: basename MUST NOT contain '_'
    IF (INDEX(basename, '_') /= 0) RETURN

    IF (PRESENT(subname)) THEN
       status = 201 ! ERROR: tracer subname too long
       IF (LEN(subname) > STRLEN_MEDIUM) RETURN
       IF (TRIM(subname) /= '') THEN
          fullname = TRIM(basename)//'_'//TRIM(subname)
       ELSE
          fullname = TRIM(basename)
       END IF
    ELSE
       fullname = TRIM(basename)
    END IF

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       IF (PRESENT(idx)) idx = 0
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (TRIM(ti%info%ident%fullname) == TRIM(fullname)) THEN
          status = TR_EXIST  ! ERROR: tracer exists already
          RETURN
       END IF
       te => ti
       ti => ti%next
    END DO

    ! ADD NEW TRACER TO LIST
    ALLOCATE(ti)
    NULLIFY(ti%next)
    IF (trset(zid)%ntrac == 0) THEN
       trset(zid)%tilist => ti         ! SET POINTER TO FIRST ELEMENT
    ELSE
       te%next => ti                   ! SET NEXT POINTER OF LAST ELEMENT
       !                               ! TO NEW ELEMENT
    END IF

    ! ADJUST TRACER INDEX
    trset(zid)%ntrac = trset(zid)%ntrac + 1
    IF (PRESENT(idx)) idx = trset(zid)%ntrac

    ! SET TRACER INFORMATION
    ! - IDENTIFICATION
    ! -- MANDATORY
    ti%info%ident%basename                   = TRIM(basename)    ! BASENAME
    !
    IF (LEN(submodel) > STRLEN_MEDIUM) THEN
       status = 204 ! ERROR: unit too long
       RETURN
    END IF
    ti%info%ident%submodel                   = TRIM(submodel)    ! SUBMODEL
    !
    ti%info%ident%fullname                   = TRIM(fullname)    ! FULLNAME
    ti%info%ident%idx                        = trset(zid)%ntrac  ! SET INDEX

    ! -- OPTIONAL
    !
    IF (PRESENT(subname)) THEN           ! SUBNAME
       IF (LEN(subname) > STRLEN_MEDIUM) THEN
          status = 201 ! ERROR: subname too long
          RETURN
       END IF
       ti%info%ident%subname    = TRIM(subname)
    END IF
    !
    !                                    ! LONGNAME
    IF (PRESENT(longname)) THEN
       IF (LEN(longname) > STRLEN_LONG) THEN
          status = 207 ! ERROR: longname too long
          RETURN
       END IF
       ti%info%ident%longname  = TRIM(longname)
    END IF
    !
    IF (PRESENT(unit)) THEN              ! UNIT
       IF (LEN(unit) > STRLEN_MEDIUM) THEN
          status = 203 ! ERROR: unit too long
          RETURN
       END IF
       ti%info%ident%unit          = TRIM(unit)
    END IF
    !
    IF (PRESENT(medium)) THEN            ! MEDIUM
       IF ((medium < 0) .OR. (medium > MAX_MEDIUM) ) THEN
          status = 300  ! ERROR: unknown medium
          RETURN
       END IF
       ti%info%ident%medium   = medium
    END IF
    !
    IF (PRESENT(quantity)) THEN          ! QUANTITY
       IF ((quantity < 0) .OR. (quantity > MAX_QUANTITY) ) THEN
          status = 301  ! ERROR: unknown quantity
          RETURN
       END IF
       ti%info%ident%quantity = quantity
    END IF
    ! 
    IF (PRESENT(type)) THEN              ! TYPE
       IF ((type < 0) .OR. (type > MAX_TYPE  ) ) THEN
          status = 302  ! ERROR: unknown type
          RETURN
       END IF
       ti%info%ident%type = type
    END IF

    ! ASSIGN OPTIONAL PARAMETERS
    ! DEFAULT
    ti%info%meta%cask_i(:) = default_cask_i(:)
    ti%info%meta%cask_r(:) = default_cask_r(:)
    ti%info%meta%cask_s(:) = default_cask_s(:)
    ! PHYSICAL PROCESSES
    IF (PRESENT(nadvect))     ti%info%meta%cask_i(I_advect)     = nadvect 
    IF (PRESENT(nconvect))    ti%info%meta%cask_i(I_convect)    = nconvect
    IF (PRESENT(nvdiff))      ti%info%meta%cask_i(I_vdiff)      = nvdiff
    IF (PRESENT(nwetdep))     ti%info%meta%cask_i(I_wetdep)     = nwetdep
    IF (PRESENT(ndrydep))     ti%info%meta%cask_i(I_drydep)     = ndrydep
    IF (PRESENT(nsedi))       ti%info%meta%cask_i(I_sedi)       = nsedi
    IF (PRESENT(nscav))       ti%info%meta%cask_i(I_scav)       = nscav
    IF (PRESENT(nmix))        ti%info%meta%cask_i(I_mix)        = nmix
    ! SPECIAL FOR COLUMN MODE
    IF (PRESENT(nforce_col))  ti%info%meta%cask_i(I_force_col)  = nforce_col
    ! NUMERICAL
    IF (PRESENT(nintegrate))  ti%info%meta%cask_i(I_integrate)  = nintegrate
    IF (PRESENT(ntimefilter)) ti%info%meta%cask_i(I_timefilter) = ntimefilter
    ! INFRASTRUCTURE
    IF (PRESENT(vini))        ti%info%meta%cask_r(R_vini)        = vini
    IF (PRESENT(lforce_init)) THEN
       IF (lforce_init) THEN
          ti%info%meta%cask_i(I_force_init) = ON
       ELSE
          ti%info%meta%cask_i(I_force_init) = OFF
       END IF
    END IF
    ! SINGLE
    IF (PRESENT(molarmass))   ti%info%meta%cask_r(R_molarmass)  = molarmass
    IF (PRESENT(henry))       ti%info%meta%cask_r(R_henry)      = henry
    IF (PRESENT(dryreac_sf))  ti%info%meta%cask_r(R_dryreac_sf) = dryreac_sf
    ! MEDIUM AEROSOL
    IF (PRESENT(m_aerosol_method))    ti%info%meta%cask_I(I_aerosol_method) = &
         m_aerosol_method
    IF (PRESENT(m_aerosol_modelname)) ti%info%meta%cask_S(S_aerosol_model) = &
         m_aerosol_modelname
    IF (PRESENT(m_aerosol_density))  ti%info%meta%cask_R(R_aerosol_density) = &
         m_aerosol_density
    IF (PRESENT(m_aerosol_mode))      ti%info%meta%cask_I(I_aerosol_mode) = &
         m_aerosol_mode
    IF (PRESENT(m_aerosol_nsol))      ti%info%meta%cask_I(I_aerosol_sol) = &
         m_aerosol_nsol

    ! special "dynamic default":
    IF (ti%info%meta%cask_i(I_TAG_REG_IDT) == 0) &
         ti%info%meta%cask_i(I_TAG_REG_IDT) = ti%info%ident%idx

    !
    status = 0

  END SUBROUTINE new_tracer_old
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_set_by_name(status, setname, trlist, ti, ntrac &
       , xt, xtte, xtm1, xmem, l_tfstd, l_init, l_enable)

    IMPLICIT NONE
    INTRINSIC :: LEN, PRESENT, TRIM

    ! I/O
    INTEGER,                         INTENT(OUT)           :: status
    CHARACTER(LEN=*),                INTENT(IN)            :: setname
    TYPE(t_trinfo_list),             POINTER,     OPTIONAL :: trlist
    TYPE(t_trinfo_tp), DIMENSION(:), POINTER,     OPTIONAL :: ti
    INTEGER,                         INTENT(OUT), OPTIONAL :: ntrac
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xt
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xtte
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xtm1
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xmem
    LOGICAL,                         INTENT(OUT), OPTIONAL :: l_tfstd
    LOGICAL,                         INTENT(OUT), OPTIONAL :: l_init
    LOGICAL,                         INTENT(OUT), OPTIONAL :: l_enable

    ! LOCAL
    INTEGER :: zid = 0

    status = 100 ! ERROR: setname too long
    IF (LEN(setname) > STRLEN_TRSET) RETURN

    status = 103 ! ERROR: set not available
    DO zid = 1, NSETID
       IF (TRIM(trset(zid)%name) == TRIM(setname)) THEN
          status = 0
          EXIT
       END IF
    END DO
    IF (status /= 0) RETURN

    IF (PRESENT(ntrac))  ntrac  =  trset(zid)%ntrac
    IF (PRESENT(trlist)) trlist => trset(zid)%tilist
    IF (PRESENT(ti))     ti     => trset(zid)%ti

    IF (PRESENT(xt))     xt     => trset(zid)%xt
    IF (PRESENT(xtte))   xtte   => trset(zid)%xtte
    IF (PRESENT(xtm1))   xtm1   => trset(zid)%xtm1
    IF (PRESENT(xmem))   xmem   => trset(zid)%xmem

    IF (PRESENT(l_tfstd))  l_tfstd  = trset(zid)%l_tfstd
    IF (PRESENT(l_init))   l_init   = trset(zid)%l_init
    IF (PRESENT(l_enable)) l_enable = trset(zid)%l_enable

    status = 0

  END SUBROUTINE get_tracer_set_by_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_set_by_id(status, setid, setname, trlist, ti, ntrac &
       , xt, xtte, xtm1, xmem, l_tfstd, l_init, l_enable)

    IMPLICIT NONE
    INTRINSIC :: PRESENT

    ! I/O
    INTEGER,                         INTENT(OUT)           :: status
    INTEGER,                         INTENT(IN)            :: setid
    CHARACTER(LEN=STRLEN_TRSET),     INTENT(OUT), OPTIONAL :: setname
    TYPE(t_trinfo_list),             POINTER,     OPTIONAL :: trlist
    TYPE(t_trinfo_tp), DIMENSION(:), POINTER,     OPTIONAL :: ti
    INTEGER,                         INTENT(OUT), OPTIONAL :: ntrac
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xt
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xtte
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xtm1
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xmem
    LOGICAL,                         INTENT(OUT), OPTIONAL :: l_tfstd
    LOGICAL,                         INTENT(OUT), OPTIONAL :: l_init
    LOGICAL,                         INTENT(OUT), OPTIONAL :: l_enable

    IF (setid > NSETID) THEN
       status = 104 ! TRACER-SET ID NOT EXISTENT
       RETURN
    ENDIF

    IF (PRESENT(setname)) setname   =  trset(setid)%name
    IF (PRESENT(trlist))  trlist    => trset(setid)%tilist
    IF (PRESENT(ti))      ti        => trset(setid)%ti
    IF (PRESENT(ntrac))   ntrac     =  trset(setid)%ntrac

    IF (PRESENT(xt))      xt        => trset(setid)%xt
    IF (PRESENT(xtte))    xtte      => trset(setid)%xtte
    IF (PRESENT(xtm1))    xtm1      => trset(setid)%xtm1
    IF (PRESENT(xmem))    xmem      => trset(setid)%xmem

    IF (PRESENT(l_tfstd))  l_tfstd  =  trset(setid)%l_tfstd
    IF (PRESENT(l_init))   l_init   =  trset(setid)%l_init
    IF (PRESENT(l_enable)) l_enable =  trset(setid)%l_enable

    status = 0

  END SUBROUTINE get_tracer_set_by_id
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_set_id(status, setname, id)

    IMPLICIT NONE
    INTRINSIC :: LEN, TRIM

    ! I/O
    INTEGER,                        INTENT(OUT)  :: status
    CHARACTER(LEN=*),               INTENT(IN)   :: setname
    INTEGER,                        INTENT(OUT)  :: id

    status = 100 ! ERROR: setname too long
    IF (LEN(setname) > STRLEN_TRSET) RETURN

    status = 103 ! ERROR: set not available
    DO id = 1, NSETID
       IF (TRIM(trset(id)%name) == TRIM(setname)) THEN
          status = 0
          EXIT
       END IF
    END DO
    IF (status /= 0) RETURN

    status = 0

  END SUBROUTINE get_tracer_set_id
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE setup_tracer_set(status, setname, dim, nt, l_tfstd, l_init &
       , nbounds) ! mz_ab_20100610
    
    IMPLICIT NONE
    INTRINSIC :: NULL

    ! I/O
    INTEGER,               INTENT(OUT) :: status
    CHARACTER(LEN=*),      INTENT(IN)  :: setname ! name of tracer set
    INTEGER, DIMENSION(3), INTENT(IN)  :: dim     ! 3 free dimensions
    INTEGER,               INTENT(IN)  :: nt      ! number of (time-) levels
    LOGICAL,               INTENT(IN)  :: l_tfstd ! level 2 and 3 specific
    LOGICAL,               INTENT(IN)  :: l_init  ! allow init
    ! mz_ab_20100509+
    INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL :: nbounds    
    ! mz_ab_20100509-

    ! LOCAL
    INTEGER :: zid, zstat
    TYPE(t_trinfo_list),   POINTER     :: til => NULL()
    INTEGER                            :: jt
    INTEGER, DIMENSION(3)              :: zdim ! mz_ab_20100509

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! SETUP ENUMERATED LIST OF TRACER META INFORMATION STRUCTURES
    ALLOCATE(trset(zid)%ti(trset(zid)%ntrac))
    til => trset(zid)%tilist
    DO jt=1, trset(zid)%ntrac
       IF (jt /= til%info%ident%idx) THEN
          status = 999 
          RETURN
       END IF
       trset(zid)%ti(jt)%tp => til%info
       til => til%next
    END DO

    ! CHECK LEVELS
    IF (nt < 1) THEN
       status = 105 ! WRONG NUMBER OF (TIME-) LEVELS
       RETURN
    END IF
    trset(zid)%nt = nt

    status = 1000 ! ERROR: MEMORY ALLOCATION FAILED
    ! mz_ab_20100528+
    IF (PRESENT(nbounds)) THEN
       zdim(:)=dim(:)+2*nbounds(:) 
    ELSE 
       zdim(:) = dim(:)
    END IF
    !!$ALLOCATE(trset(zid)%mem( dim(1),dim(2),trset(zid)%ntrac,dim(3),1,nt), &
    !!$     STAT=zstat)
    ALLOCATE(trset(zid)%mem( zdim(1),zdim(2),trset(zid)%ntrac,zdim(3),1,nt), &
         STAT=zstat)
    ! mz_ab_20100528-

    IF (zstat /= 0) RETURN
    trset(zid)%mem(:,:,:,:,:,:) = 0.0_DP

    trset(zid)%xt   => trset(zid)%mem(:,:,:,:,:,1)

    IF ((l_tfstd) .AND. (nt >= 3)) THEN
       trset(zid)%l_tfstd = .TRUE.
       trset(zid)%xtte => trset(zid)%mem(:,:,:,:,:,2)
       trset(zid)%xtm1 => trset(zid)%mem(:,:,:,:,:,3)
       IF (nt > 3) THEN
          trset(zid)%xmem => trset(zid)%mem(:,:,:,:,1,4:)
          trset(zid)%nex = 4
       ELSE
          trset(zid)%xmem => NULL() 
          trset(zid)%nex = 0
       END IF
    ELSE
       trset(zid)%l_tfstd = .FALSE.
       trset(zid)%xtte => NULL()
       trset(zid)%xtm1 => NULL()
       IF (nt > 1) THEN
          trset(zid)%xmem => trset(zid)%mem(:,:,:,:,1,2:)
          trset(zid)%nex = 2
       ELSE
          trset(zid)%xmem => NULL() 
          trset(zid)%nex = 0
       END IF
    END IF

    trset(zid)%l_init = l_init

    status = 0

  END SUBROUTINE setup_tracer_set
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE clean_tracer_set(status, setname)
    
    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,               INTENT(OUT) :: status
    CHARACTER(LEN=*),      INTENT(IN)  :: setname
    
    ! LOCAL
    INTEGER                      :: zid, zstat
    TYPE(t_trinfo_list), POINTER :: ti => NULL()
    TYPE(t_trinfo_list), POINTER :: te => NULL()
    INTEGER                      :: jt
    
    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    status = 1001 ! ERROR: MEMORY DEALLOCATION FAILED
    IF (ASSOCIATED(trset(zid)%mem)) THEN
      DEALLOCATE(trset(zid)%mem,   STAT=zstat)
      IF (zstat /= 0) RETURN
    END IF
    
    ! DELETE ENUMERATED LIST OF TRACER META INFORMATION
    DO jt=1, trset(zid)%ntrac
       NULLIFY(trset(zid)%ti(jt)%tp)
    END DO
    IF (ASSOCIATED(trset(zid)%ti)) DEALLOCATE(trset(zid)%ti)

    ! DELETE TRACER INFO STRUCT
    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       te => ti%next
       DEALLOCATE(ti, STAT=zstat)
       IF (zstat /= 0) THEN
          status = 1001 ! ERRRO: memory deallocation failed
          RETURN
       END IF
       ti => te
    END DO
    
    status = 0

  END SUBROUTINE clean_tracer_set
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE print_tracer_set

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, TRIM

    ! LOCAL
    INTEGER                        :: zid, ntrac, i, l
    TYPE(t_trinfo_list), POINTER   :: ti => NULL()
    TYPE(t_trinfo_list), POINTER   :: te => NULL()
    CHARACTER(LEN=200)                 :: form = ''
    CHARACTER(LEN=NAMES_CASK_I_STRLEN) :: tmpstr = ''
    CHARACTER(LEN=2*MAX_CASK_I+1)      :: oline = ''

    WRITE(form,*) '(1x,i3,1x,a15,1x,a8,1x,3(i1,1x),',MAX_CASK_I,'(1x,i1))'

    set_loop: DO zid = 1, NSETID

    ! NUMBER OF TRACERS
    ntrac = trset(zid)%ntrac

    WRITE(*,*) '============================================================='
    WRITE(*,*) 'TRACER-SET        : ', TRIM(trset(zid)%name)
    IF (trset(zid)%l_enable) THEN
       WRITE(*,*) 'STATUS            : enabled'
    ELSE
       WRITE(*,*) 'STATUS            : disabled'
    END IF
    WRITE(*,*) 'NUMBER OF TRACERS : ', ntrac
    WRITE(*,*)

    IF (ntrac > 0) THEN
    WRITE(*,*) &
         '----IDENT-------------------------  ----META--------------------- '

    DO l=1, NAMES_CASK_I_STRLEN
       DO i=1, MAX_CASK_I
          tmpstr = NAMES_CASK_I(i)
          WRITE(oline(2*i-1:2*i),'(a1,a1)') ' ',tmpstr(l:l)
       END DO
       WRITE(*,*) '                                   ', oline
    END DO
    WRITE(*,*) '                             M Q T                           '
    WRITE(*,*) 'idt fullname(1:15)  submodel                                 '
    WRITE(*,*) 
    END IF

    ! LOOP OVER TRACERS IN SET
    ti => trset(zid)%tilist
    tracer_loop: DO 
       IF (.NOT. ASSOCIATED(ti)) EXIT

       WRITE(*,form) &
          ti%info%ident%idx,             &
          ti%info%ident%fullname(1:15),  &
          ti%info%ident%submodel(1:8),   &
          ti%info%ident%medium,          &
          ti%info%ident%quantity,        &
          ti%info%ident%type,            &
          ti%info%meta%cask_i(:)

        ti => ti%next
    END DO tracer_loop

    WRITE(*,*) &
         '=================================================================='

    END DO set_loop

  END SUBROUTINE print_tracer_set
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE print_tracer_set_val

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, MINVAL, MAXVAL, TRIM, SIZE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: fstr1 = '(3x,a4,2x,e12.4,a5,e12.4,2x,a)'
    CHARACTER(LEN=*), PARAMETER  :: fstr2 = '(3x,i4.4,2x,e12.4,a5,e12.4,2x,a)'
    INTEGER                      :: zid, ntrac, jt
    TYPE(t_trinfo_list), POINTER :: ti => NULL()
    TYPE(t_trinfo_list), POINTER :: te => NULL()
    INTEGER                      :: n

    set_loop: DO zid=1, NSETID

    ! NUMBER OF TRACERS
    ntrac = trset(zid)%ntrac

    WRITE(*,*) '============================================================='
    WRITE(*,*) 'TRACER-SET        : ', TRIM(trset(zid)%name)
    IF (trset(zid)%l_enable) THEN
       WRITE(*,*) 'STATUS            : enabled'
    ELSE
       WRITE(*,*) 'STATUS            : disabled'
    END IF
    WRITE(*,*) 'NUMBER OF TRACERS : ', ntrac
    WRITE(*,*)

    ! LOOP OVER TRACERS IN SET
    ti => trset(zid)%tilist
    tracer_loop: DO 
       IF (.NOT. ASSOCIATED(ti)) EXIT

       jt = ti%info%ident%idx
       WRITE(*,*) jt &
            , ti%info%ident%basename,' - ', ti%info%ident%subname, ' - '  &
            , ti%info%ident%submodel

       WRITE(*,fstr1) 'XT  ', MINVAL(trset(zid)%xt(:,:,jt,:,:)), ' ... '  &
                            , MAXVAL(trset(zid)%xt(:,:,jt,:,:))           &
                            , ti%info%ident%unit

       IF (ASSOCIATED(trset(zid)%xtm1)) THEN
          WRITE(*,fstr1) 'XTM1', MINVAL(trset(zid)%xtm1(:,:,jt,:,:)), ' ... ' &
                               , MAXVAL(trset(zid)%xtm1(:,:,jt,:,:))          &
                               , ti%info%ident%unit
       END IF

       IF (ASSOCIATED(trset(zid)%xtte)) THEN
          WRITE(*,fstr1) 'XTTE', MINVAL(trset(zid)%xtte(:,:,jt,:,:)), ' ... ' &
                               , MAXVAL(trset(zid)%xtte(:,:,jt,:,:))          &
                               , ti%info%ident%unit

       END IF

       IF (ASSOCIATED(trset(zid)%xmem)) THEN
          DO n=1, SIZE(trset(zid)%xmem, 5)
             WRITE(*,fstr2) n, MINVAL(trset(zid)%xmem(:,:,jt,:,n)), ' ... ' &
                             , MAXVAL(trset(zid)%xmem(:,:,jt,:,n))          &
                             , ti%info%ident%unit

          END DO
       END IF

       ti => ti%next
    END DO tracer_loop

    END DO set_loop

    WRITE(*,*) '============================================================='

  END SUBROUTINE print_tracer_set_val
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_by_name(status, setname, basename &
       , subname, idx, fullname, longname, unit, submodel &
       , medium, quantity, type                           &
       , trinfo, pxt, pxtm1, pxtte, pxmem)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, INDEX, LEN, PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    !
    CHARACTER(LEN=*), INTENT(IN)            :: basename
    !
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: subname
    INTEGER,          INTENT(OUT), OPTIONAL :: idx
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: fullname
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: longname
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: unit
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: submodel
    INTEGER,          INTENT(OUT), OPTIONAL :: medium
    INTEGER,          INTENT(OUT), OPTIONAL :: quantity
    INTEGER,          INTENT(OUT), OPTIONAL :: type
    !
    TYPE(t_trinfo),   INTENT(OUT), OPTIONAL :: trinfo
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxt
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxtm1
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxtte
    REAL(DP), DIMENSION(:,:,:,:), POINTER, OPTIONAL :: pxmem

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    CHARACTER(LEN=STRLEN_FNAME)   :: zfullname
    INTEGER                       :: zid
    INTEGER                       :: zidx

    ! INIT
    IF (PRESENT(idx))        idx       = 0
    IF (PRESENT(fullname))   fullname  = ''
    IF (PRESENT(longname))   longname  = ''
    IF (PRESENT(unit))       unit      = ''
    IF (PRESENT(submodel))   submodel  = ''
    IF (PRESENT(medium))     medium    = 0
    IF (PRESENT(quantity))   quantity  = 0
    IF (PRESENT(type))       type      = 0
    !
    ! trinfo
    !
    IF (PRESENT(pxt))   NULLIFY(pxt)
    IF (PRESENT(pxtm1)) NULLIFY(pxtm1)
    IF (PRESENT(pxtte)) NULLIFY(pxtte)
    IF (PRESENT(pxmem)) NULLIFY(pxmem)

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    status = 208 ! ERROR: basename MUST NOT contain '_'
    IF (INDEX(basename, '_') /= 0) RETURN

    status = 200 ! ERROR: tracer basename too long
    IF (LEN(basename) > STRLEN_MEDIUM) RETURN
    
    IF (PRESENT(subname)) THEN
       status = 201 ! ERROR: tracer subname too long
       IF (LEN(subname) > STRLEN_MEDIUM) RETURN
       IF (TRIM(subname) /= '') THEN
          zfullname = TRIM(basename)//'_'//TRIM(subname)
       ELSE
          zfullname = TRIM(basename)
       END IF
    ELSE
       zfullname = TRIM(basename)
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (TRIM(ti%info%ident%fullname) == TRIM(zfullname)) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    ! RETURN TRACER INFORMATION
    ! - IDENTIFICATION
    zidx = ti%info%ident%idx
    IF (PRESENT(idx)) idx = zidx
    !
    IF (PRESENT(pxt)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xt)) THEN
          status = 600  ! ERROR: xt pointer not associated
          RETURN
       ELSE
          pxt => trset(zid)%xt(:,:,zidx,:,1)
       END IF
    END IF
    !
    IF (PRESENT(pxtm1)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xtm1)) THEN
          status = 601  ! ERROR: xtm1 pointer not associated
          RETURN
       ELSE
          pxtm1 => trset(zid)%xtm1(:,:,zidx,:,1)
       END IF
    END IF
    !
    IF (PRESENT(pxtte)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xtte)) THEN
          status = 602  ! ERROR: xtte pointer not associated
          RETURN
       ELSE
          pxtte => trset(zid)%xtte(:,:,zidx,:,1)
       END IF
    END IF
    !
    IF (PRESENT(pxmem)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xmem)) THEN
          status = 603  ! ERROR: xmem pointer not associated
          RETURN
       ELSE
          pxmem => trset(zid)%xmem(:,:,zidx,:,:)
       END IF
    END IF

    ! - OPTIONAL INFORMATION
    IF (PRESENT(fullname)) THEN
       IF (LEN(fullname) > STRLEN_FNAME) THEN
          status = 206 ! ERROR: fullname too long
          RETURN
       END IF
       fullname = TRIM(ti%info%ident%fullname)
    END IF
    !
    !                                    ! LONGNAME
    IF (PRESENT(longname)) THEN
       IF (LEN(longname) > STRLEN_LONG) THEN
          status = 207 ! ERROR: longname too long
          RETURN
       END IF
       longname = TRIM(ti%info%ident%longname)
    END IF
    !
    IF (PRESENT(unit)) THEN              ! UNIT
       IF (LEN(unit) > STRLEN_MEDIUM) THEN
          status = 203 ! ERROR: unit too long
          RETURN
       END IF
       unit = TRIM(ti%info%ident%unit)
    END IF
    !
    IF (PRESENT(submodel)) THEN          ! SUBMODEL
       IF (LEN(submodel) > STRLEN_MEDIUM) THEN
          status = 204 ! ERROR: submodel too long
          RETURN
       END IF
       submodel = TRIM(ti%info%ident%submodel)
    END IF
    !
    IF (PRESENT(medium)) THEN            ! MEDIUM
       IF ((ti%info%ident%medium < 0) .OR. &
            (ti%info%ident%medium > MAX_MEDIUM) ) THEN
          status = 300  ! ERROR: unknown medium
          RETURN
       END IF
       medium = ti%info%ident%medium
    END IF
    !
    IF (PRESENT(quantity)) THEN          ! QUANTITY
       IF ((ti%info%ident%quantity < 0) .OR. &
            (ti%info%ident%quantity > MAX_QUANTITY) ) &
            THEN
          status = 301  ! ERROR: unknown quantity
          RETURN
       END IF
       quantity = ti%info%ident%quantity
    END IF
    ! 
    IF (PRESENT(type)) THEN              ! TYPE
       IF ((ti%info%ident%type < 0) .OR. &
            (ti%info%ident%type > MAX_MEDIUM) ) THEN
          status = 302  ! ERROR: unknown type
          RETURN
       END IF
       type = ti%info%ident%type
    END IF

    IF (PRESENT(trinfo)) trinfo = ti%info

    status = 0

  END SUBROUTINE get_tracer_by_name
  ! -------------------------------------------------------------------

! mz_pj_20090507+
  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_by_id(status, setname, id, basename  &
       , subname, fullname, longname, unit, submodel         &
       , medium, quantity, type                              &
       , trinfo, pxt, pxtm1, pxtte, pxmem)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, LEN, PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: id
    !
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: basename
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: subname
    !
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: fullname
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: longname
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: unit
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: submodel
    INTEGER,          INTENT(OUT), OPTIONAL :: medium
    INTEGER,          INTENT(OUT), OPTIONAL :: quantity
    INTEGER,          INTENT(OUT), OPTIONAL :: type
    !
    TYPE(t_trinfo),   INTENT(OUT), OPTIONAL :: trinfo
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxt
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxtm1
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxtte
    REAL(DP), DIMENSION(:,:,:,:), POINTER, OPTIONAL :: pxmem

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! INIT
    IF (PRESENT(basename))   basename  = ''
    IF (PRESENT(subname))    subname   = ''

    IF (PRESENT(fullname))   fullname  = ''
    IF (PRESENT(longname))   longname  = ''
    IF (PRESENT(unit))       unit      = ''
    IF (PRESENT(submodel))   submodel  = ''
    IF (PRESENT(medium))     medium    = 0
    IF (PRESENT(quantity))   quantity  = 0
    IF (PRESENT(type))       type      = 0
    !
    ! trinfo
    !
    IF (PRESENT(pxt))   NULLIFY(pxt)
    IF (PRESENT(pxtm1)) NULLIFY(pxtm1)
    IF (PRESENT(pxtte)) NULLIFY(pxtte)
    IF (PRESENT(pxmem)) NULLIFY(pxmem)

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == id) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    ! RETURN TRACER INFORMATION
    ! - IDENTIFICATION
    IF (PRESENT(basename)) THEN
       IF (LEN(basename) > STRLEN_MEDIUM) THEN
          status = 200 ! ERROR: basename too long
          RETURN
       END IF
       basename = TRIM(ti%info%ident%basename)
    END IF

    IF (PRESENT(subname)) THEN
       IF (LEN(subname) > STRLEN_MEDIUM) THEN
          status = 201 ! ERROR: subname too long
          RETURN
       END IF
       subname = TRIM(ti%info%ident%subname)
    END IF

    IF (PRESENT(pxt)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xt)) THEN
          status = 600  ! ERROR: xt pointer not associated
          RETURN
       ELSE
          pxt => trset(zid)%xt(:,:,id,:,1)
       END IF
    END IF
    !
    IF (PRESENT(pxtm1)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xtm1)) THEN
          status = 601  ! ERROR: xtm1 pointer not associated
          RETURN
       ELSE
          pxtm1 => trset(zid)%xtm1(:,:,id,:,1)
       END IF
    END IF
    !
    IF (PRESENT(pxtte)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xtte)) THEN
          status = 602  ! ERROR: xtte pointer not associated
          RETURN
       ELSE
          pxtte => trset(zid)%xtte(:,:,id,:,1)
       END IF
    END IF
    !
    IF (PRESENT(pxmem)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xmem)) THEN
          status = 603  ! ERROR: xmem pointer not associated
          RETURN
       ELSE
          pxmem => trset(zid)%xmem(:,:,id,:,:)
       END IF
    END IF

    ! - OPTIONAL INFORMATION
    IF (PRESENT(fullname)) THEN
       IF (LEN(fullname) > STRLEN_FNAME) THEN
          status = 206 ! ERROR: fullname too long
          RETURN
       END IF
       fullname = TRIM(ti%info%ident%fullname)
    END IF
    !
    !                                    ! LONGNAME
    IF (PRESENT(longname)) THEN
       IF (LEN(longname) > STRLEN_LONG) THEN
          status = 207 ! ERROR: longname too long
          RETURN
       END IF
       longname = TRIM(ti%info%ident%longname)
    END IF
    !
    IF (PRESENT(unit)) THEN              ! UNIT
       IF (LEN(unit) > STRLEN_MEDIUM) THEN
          status = 203 ! ERROR: unit too long
          RETURN
       END IF
       unit = TRIM(ti%info%ident%unit)
    END IF
    !
    IF (PRESENT(submodel)) THEN          ! SUBMODEL
       IF (LEN(submodel) > STRLEN_MEDIUM) THEN
          status = 204 ! ERROR: submodel too long
          RETURN
       END IF
       submodel = TRIM(ti%info%ident%submodel)
    END IF
    !
    IF (PRESENT(medium)) THEN            ! MEDIUM
       IF ((ti%info%ident%medium < 0) .OR. &
            (ti%info%ident%medium > MAX_MEDIUM) ) THEN
          status = 300  ! ERROR: unknown medium
          RETURN
       END IF
       medium = ti%info%ident%medium
    END IF
    !
    IF (PRESENT(quantity)) THEN          ! QUANTITY
       IF ((ti%info%ident%quantity < 0) .OR. &
            (ti%info%ident%quantity > MAX_QUANTITY) ) &
            THEN
          status = 301  ! ERROR: unknown quantity
          RETURN
       END IF
       quantity = ti%info%ident%quantity
    END IF
    ! 
    IF (PRESENT(type)) THEN              ! TYPE
       IF ((ti%info%ident%type < 0) .OR. &
            (ti%info%ident%type > MAX_MEDIUM) ) THEN
          status = 302  ! ERROR: unknown type
          RETURN
       END IF
       type = ti%info%ident%type
    END IF

    IF (PRESENT(trinfo)) trinfo = ti%info

    status = 0

  END SUBROUTINE get_tracer_by_id
  ! -------------------------------------------------------------------
! mz_pj_20090507-

! op_pj_20100325+
  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_s(status, setname, idx, flag, s)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    INTEGER,          INTENT(IN)            :: flag
    CHARACTER(LEN=*), INTENT(OUT)           :: s

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF
    
    IF ( (flag < 1) .OR. (flag > MAX_CASK_S) ) THEN
       status = 702 ! STRING FLAG INDEX OUT OF RANGE
       RETURN
    END IF
    s = ti%info%meta%cask_s(flag)

    status = 0
    
  END SUBROUTINE get_tracer_s
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_i(status, setname, idx, flag, i)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    INTEGER,          INTENT(IN)            :: flag
    INTEGER,          INTENT(OUT)           :: i

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    IF ( (flag < 1) .OR. (flag > MAX_CASK_I) ) THEN
       status = 700 ! INTEGER FLAG INDEX OUT OF RANGE
       RETURN
    END IF
    i = ti%info%meta%cask_i(flag)

    status = 0

  END SUBROUTINE get_tracer_i
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_r(status, setname, idx, flag, r)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    INTEGER,          INTENT(IN)            :: flag
    REAL(DP),         INTENT(OUT)           :: r

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN
    
    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    IF ( (flag < 1) .OR. (flag > MAX_CASK_R) ) THEN
       status = 701 ! REAL FLAG INDEX OUT OF RANGE
       RETURN
    END IF
    r = ti%info%meta%cask_r(flag)

    status = 0

  END SUBROUTINE get_tracer_r
  ! -------------------------------------------------------------------
! op_pj_20100325-

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_list(status, setname, basename, idxs, subnames)


    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, INDEX, LEN, PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    CHARACTER(LEN=*), INTENT(IN)            :: basename
    INTEGER, DIMENSION(:), POINTER          :: idxs       ! INTENT(OUT)
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER, OPTIONAL :: subnames

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid
    INTEGER                       :: n

    ! INIT
    IF (ASSOCIATED(idxs)) THEN
       DEALLOCATE(idxs)
       NULLIFY(idxs)
    END IF

    IF (PRESENT(subnames)) THEN
       IF (ASSOCIATED(subnames)) THEN
          DEALLOCATE(subnames)
          NULLIFY(subnames)
       END IF
    END IF

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    status = 208 ! ERROR: basename MUST NOT contain '_'
    IF (INDEX(basename, '_') /= 0) RETURN

    status = 200 ! ERROR: tracer basename too long
    IF (LEN(basename) > STRLEN_MEDIUM) RETURN
    
    ! SEARCH FOR EXISITNG TRACERS IN LIST
    n = 0
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (TRIM(ti%info%ident%basename) == TRIM(basename)) n=n+1 ! FOUND
       ti => ti%next
    END DO

    ! ALLOCATE OUTPUT
    ALLOCATE(idxs(n))
    idxs(:) = 0
    IF (PRESENT(subnames)) THEN
       ALLOCATE(subnames(n))
       subnames(:) = ''
    END IF

    ! SEARCH FOR EXISITNG TRACERS IN LIST
    n = 0
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (TRIM(ti%info%ident%basename) == TRIM(basename)) THEN
          n=n+1
          idxs(n) = ti%info%ident%idx
          IF (PRESENT(subnames)) subnames(n) = TRIM(ti%info%ident%subname)
       END IF
       ti => ti%next
    END DO

    status = 0

  END SUBROUTINE get_tracer_list
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
!!$  SUBROUTINE tracer_iniflag(status, setname, id, linit)
  SUBROUTINE tracer_iniflag(status, setname, id, lset, lget)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: setname
    INTEGER,          INTENT(IN)  :: id
!!$    LOGICAL,          INTENT(IN)  :: linit
    LOGICAL,          INTENT(IN),  OPTIONAL  :: lset
    LOGICAL,          INTENT(OUT), OPTIONAL  :: lget

    ! LOCAL
    INTEGER                      :: zid
    TYPE(t_trinfo_list), POINTER :: ti => NULL()
    LOGICAL                      :: lfound

    IF (PRESENT(lget)) lget = .FALSE.

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK, IF SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    lfound = .FALSE.

    ! SEARCH FOR EXISITNG TRACER IN LIST WITH THIS ID
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == id) THEN
          lfound = .TRUE.
          EXIT   ! FOUND
       END IF
       ti => ti%next
    END DO

    IF (.NOT. lfound) THEN
       status = TR_NEXIST ! TRACER DOES NOT EXIST
       RETURN
    END IF

    ! SET / GET FLAG
!!$    ti%info%meta%linit = linit
    IF (PRESENT(lset)) ti%info%meta%linit = lset
    IF (PRESENT(lget)) lget = ti%info%meta%linit

  END SUBROUTINE tracer_iniflag
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  FUNCTION tracer_error_str(status)

    USE messy_main_tools, ONLY: int2str

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=STRLEN_VLONG)           :: tracer_error_str
    INTEGER,                   INTENT(IN) :: status

    ! LOCAL
    CHARACTER(LEN=4) :: echar 

    CALL int2str(echar, status, cpad='0', cerr='*')

    SELECT CASE(echar)
       ! NO ERROR
    CASE('0000')
       tracer_error_str = 'E'//echar//': NO ERROR'

       ! TRACER-SET ERRORS
    CASE('0100')
       tracer_error_str = 'E'//echar//': TRACER-SET NAME TOO LONG'
    CASE('0101')
       tracer_error_str = 'E'//echar//': TRACER-SET NAME EXISTS ALREADY'
    CASE('0102')
       tracer_error_str = 'E'//echar//': TRACER-SET LIST FULL'
    CASE('0103')
       tracer_error_str = 'E'//echar//': TRACER-SET NAME NOT EXISTENT'
    CASE('0104')
       tracer_error_str = 'E'//echar//': TRACER-SET ID NOT EXISTENT'
    CASE('0105')
       tracer_error_str = 'E'//echar//': WRONG NUMBER OF (TIME-) LEVELS'

       ! TRACER ERRORS
    CASE('0200')
       tracer_error_str = 'E'//echar//': TRACER BASENAME TOO LONG'
    CASE('0201')
       tracer_error_str = 'E'//echar//': TRACER SUBNAME TOO LONG'
    CASE('0202')
       tracer_error_str = 'E'//echar//': TRACER EXISTS ALREADY'
    CASE('0203')
       tracer_error_str = 'E'//echar//': TRACER UNIT TOO LONG'
    CASE('0204')
       tracer_error_str = 'E'//echar//': TRACER SUBMODEL NAME TOO LONG'
    CASE('0205')
       tracer_error_str = 'E'//echar//': TRACER DOES NOT EXIST'
    CASE('0206')
       tracer_error_str = 'E'//echar//': TRACER FULLNAME TOO LONG'
    CASE('0207')
       tracer_error_str = 'E'//echar//': TRACER LONGNAME TOO LONG'
    CASE('0208')
     tracer_error_str = 'E'//echar//': TRACER BASENAME MUST NOT CONTAIN ''_'''

       ! TYPE ERRORS
    CASE('0300')
       tracer_error_str = 'E'//echar//': UNKNOWN TRACER MEDIUM'
    CASE('0301')
       tracer_error_str = 'E'//echar//': UNKNOWN TRACER QUANTITY'
    CASE('0302')
       tracer_error_str = 'E'//echar//': UNKNOWN TRACER TYPE'

       ! POINTER ERRORS
    CASE('0600')
       tracer_error_str = 'E'//echar//': POINTER XT NOT ASSOCIATED'
    CASE('0601')
       tracer_error_str = 'E'//echar//': POINTER XTM1 NOT ASSOCIATED'
    CASE('0602')
       tracer_error_str = 'E'//echar//': POINTER XTTE NOT ASSOCIATED'
    CASE('0603')
       tracer_error_str = 'E'//echar//': POINTER XMEM NOT ASSOCIATED'

       ! META INFORMATION ERRORS
    CASE('0700')
       tracer_error_str = 'E'//echar//': INTEGER FLAG INDEX OUT OF RANGE'
    CASE('0701')
       tracer_error_str = 'E'//echar//': REAL FLAG INDEX OUT OF RANGE'
    CASE('0702')
       tracer_error_str = 'E'//echar//': STRING FLAG INDEX OUT OF RANGE'
       
       ! MEMORY ERRORS
    CASE('1000')
       tracer_error_str = 'E'//echar//': MEMORY ALLOCATION FAILED'
    CASE('1001')
       tracer_error_str = 'E'//echar//': MEMORY DEALLOCATION FAILED'

       ! TRACER FAMILY ERROR
    CASE('2000')
       tracer_error_str = 'E'//echar//': TRACER MODE IS ALREADY ACTIVE'
    CASE('2001')
       tracer_error_str = 'E'//echar//': FAMILY MODE IS ALREADY ACTIVE'
    CASE('2010')
       tracer_error_str = 'E'//echar//': UNKOWN CONVERSION FLAG'

    CASE DEFAULT
       tracer_error_str = 'E'//echar//': UNKONW ERROR STATUS'
    END SELECT

  END FUNCTION tracer_error_str
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  FUNCTION PARAM2STRING(i, mode)

    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, TRIM

    ! I/O
    CHARACTER(LEN=STRLEN_MEDIUM)               :: PARAM2STRING
    INTEGER,                     INTENT(IN)    :: i
    CHARACTER(LEN=*),            INTENT(IN)    :: mode

    PARAM2STRING = ''

    SELECT CASE(TRIM(ADJUSTL(mode)))
       !
       CASE('switch')
          SELECT CASE(i)
          CASE(ON)
             PARAM2STRING = 'ON' 
          CASE(OFF)
             PARAM2STRING = 'OFF'
          CASE(MODAL)
             PARAM2STRING = 'MODAL'
          CASE(BIN)
             PARAM2STRING = 'BIN'
          CASE DEFAULT
             PARAM2STRING = 'UNDEFINED'
          END SELECT

       CASE('type')
          SELECT CASE(i)
          CASE(SINGLE)
             PARAM2STRING = 'SINGLE'
          CASE(FAMILY)
             PARAM2STRING = 'FAMILY'
          CASE(ISOTOPE)
             PARAM2STRING = 'ISOTOPE'
          CASE DEFAULT
             PARAM2STRING = 'UNDEFINED'
          END SELECT

       CASE('medium')
          SELECT CASE(i)
          CASE(AIR)
             PARAM2STRING = 'AIR'
          CASE(AEROSOL)
             PARAM2STRING = 'AEROSOL'
          CASE(CLOUD)
             PARAM2STRING = 'CLOUD'
          CASE(OCEAN)
             PARAM2STRING = 'OCEAN'
          CASE(LAKE)
             PARAM2STRING = 'LAKE'
          CASE(RIVER)
             PARAM2STRING = 'RIVER'
          CASE(LANDICE)
             PARAM2STRING = 'LANDICE'
          CASE(SEAICE)
             PARAM2STRING = 'SEAICE'
          CASE(VEGETATION)
             PARAM2STRING = 'VEGETATION'
          CASE DEFAULT
             PARAM2STRING = 'UNDEFINED'
          END SELECT

       CASE('quantity')
          SELECT CASE(i)
          CASE(AMOUNTFRACTION)
             PARAM2STRING = 'AMOUNTFRACTION'
          CASE(NUMBERDENSITY)
             PARAM2STRING = 'NUMBERDENSITY'
          CASE(CONCENTRATION)
             PARAM2STRING = 'CONCENTRATION'
          CASE DEFAULT
             PARAM2STRING = 'UNDEFINED'
          END SELECT

       CASE DEFAULT
          PARAM2STRING = 'UNDEFINED'
    END SELECT

  END FUNCTION PARAM2STRING
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE full2base_sub(status, fullname, basename, subname)

    IMPLICIT NONE

    INTRINSIC :: LEN, INDEX, TRIM

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: fullname
    CHARACTER(LEN=*), INTENT(OUT) :: basename
    CHARACTER(LEN=*), INTENT(OUT) :: subname

    ! LOCAL
    INTEGER :: si

    status = 206 ! ERROR: fullname too long
    IF (LEN(fullname) > STRLEN_FNAME) RETURN

    status = 200 ! ERROR: tracer basename too long
    IF (LEN(basename) > STRLEN_MEDIUM) RETURN

    status = 201 ! ERROR: tracer subname too long
    IF (LEN(basename) > STRLEN_MEDIUM) RETURN

    ! INIT
    basename = ''
    subname = ''

    si = INDEX(fullname, '_')
    IF (si == 0) THEN
       basename = TRIM(fullname)
       subname  = ''
    ELSE
       basename = TRIM(fullname(:si-1))
       subname  = TRIM(fullname(si+1:))
    END IF

    status = 0

  END SUBROUTINE full2base_sub
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_read_nml_ctrl(status, iou)

    ! TRACER MODULE ROUTINE (CORE)
    !
    ! READ TRACER NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Jul 2003

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ l_family, l_pdef

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    IF (l_family) THEN
       WRITE(*,*) '  TRACER FAMILIES              : ON'
    ELSE
       WRITE(*,*) '  TRACER FAMILIES              : OFF'
    END IF

    IF (l_pdef) THEN
       WRITE(*,*) '  TRACER PDEF                  : ON'
    ELSE
       WRITE(*,*) '  TRACER PDEF                  : OFF'
    END IF

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_tracer_read_nml_ctrl
  ! -------------------------------------------------------------------
  
! **********************************************************************
END MODULE messy_main_tracer
! **********************************************************************
