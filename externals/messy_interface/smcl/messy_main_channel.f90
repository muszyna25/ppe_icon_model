! **********************************************************************
MODULE messy_main_channel
! **********************************************************************

  ! MESSY DATA TRANSFER AND EXPORT INTERFACE (MEMORY MANAGEMENT)
  !
  ! Author: Patrick Joeckel, MPICH, May 2005

!!$#if (defined MPI) && (defined PNETCDF)
!!$#define ZPNETCDF
!!$#endif

  USE messy_main_constants_mem,      ONLY: DP, SP, I8       &
       , STRLEN_MEDIUM, STRLEN_SHORT
  USE messy_main_channel_attributes, ONLY: t_attribute_list, add_attribute &
                                         , return_attribute, AF_NONE
  USE messy_main_channel_dimensions, ONLY: t_dimension_ptr, DIMID_UNDEF
  USE messy_main_channel_repr,       ONLY: IRANK, REPR_UNDEF &
                                         , t_representation &
                                         , get_representation
  USE messy_main_tools,              ONLY: PTR_4D_ARRAY

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: HUGE, NULL, RESHAPE

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'channel'
  CHARACTER(len=*), PARAMETER, PUBLIC :: modver = '2.1'

  PUBLIC :: DP, SP, I8, REPR_UNDEF, DIMID_UNDEF

  ! FILE-TYPES
  INTEGER, PARAMETER, PUBLIC :: FTYPE_UNDEFINED = -1
  INTEGER, PARAMETER, PUBLIC :: FTYPE_ASCII     = 1
!!$#ifdef ZPNETCDF
  INTEGER, PARAMETER, PUBLIC :: FTYPE_NETCDF    = 2
  INTEGER, PARAMETER, PUBLIC :: FTYPE_PNETCDF   = 3
!!$#else
!!$  INTEGER, PARAMETER, PUBLIC :: FTYPE_NETCDF    = 2
!!$!!$  INTEGER, PARAMETER, PUBLIC :: FTYPE_PNETCDF   = 3
!!$  INTEGER, PARAMETER, PUBLIC :: FTYPE_PNETCDF   = 2
!!$#endif
  INTEGER, PARAMETER, PUBLIC :: FTYPE_GRIB      = 4
  INTEGER, PARAMETER, PUBLIC :: FTYPE_HDF4      = 5
  INTEGER, PARAMETER, PUBLIC :: FTYPE_HDF5      = 6
  INTEGER, PARAMETER, PUBLIC :: FTYPE_MAXIMUM   = 6
  ! NOTE: DO NOT SET FTYPE_UNDEFINED AS DEFAULT
!!$#ifdef ZPNETCDF
!!$  INTEGER, PARAMETER, PUBLIC :: FTYPE_DEFAULT   = FTYPE_PNETCDF
!!$#else
  INTEGER, PARAMETER, PUBLIC :: FTYPE_DEFAULT   = FTYPE_NETCDF
!!$#endif
  CHARACTER(LEN=4), DIMENSION(FTYPE_MAXIMUM), PARAMETER, PUBLIC :: &
       FTYPE_EXT_TEXT = (/'.asc','.nc ','.nc ','.grb','.hdf','.hdf'/)

  ! INPUT/OUTPUT MODES
  INTEGER, PARAMETER, PUBLIC :: IOMODE_OUT = 1  ! OUTPUT
  INTEGER, PARAMETER, PUBLIC :: IOMODE_RST = 2  ! RESTART
  INTEGER, PARAMETER, PUBLIC :: IOMODE_MAX = 2  ! NUMBER OF O-MODES
  CHARACTER(LEN=*), DIMENSION(IOMODE_MAX), PARAMETER, PUBLIC :: &
       IOMODE_TEXT = (/'output ', 'restart'/)
  ! ACCESS MODES
  INTEGER, PARAMETER, PUBLIC :: AMODE_WRITE = 1
  INTEGER, PARAMETER, PUBLIC :: AMODE_READ  = 2
  INTEGER, PARAMETER, PUBLIC :: AMODE_INIT  = 3

  ! PRIMARY AND SECONDARY DATA TYPE
  INTEGER, PARAMETER, PUBLIC :: SND_UNDEF = 0
  INTEGER, PARAMETER, PUBLIC :: SND_INS = 1  !!! DO NOT CHANGE POSITION OF _INS
  INTEGER, PARAMETER, PUBLIC :: SND_AVE = 2
  INTEGER, PARAMETER, PUBLIC :: SND_STD = 3
  INTEGER, PARAMETER, PUBLIC :: SND_MIN = 4
  INTEGER, PARAMETER, PUBLIC :: SND_MAX = 5
  INTEGER, PARAMETER, PUBLIC :: SND_CNT = 6
  INTEGER, PARAMETER, PUBLIC :: SND_CAV = 7
  INTEGER, PARAMETER, PUBLIC :: SND_MAXLEN = 7
  CHARACTER(LEN=4), DIMENSION(SND_MAXLEN, IOMODE_MAX), PARAMETER, PUBLIC :: &
       SND_TEXT  = RESHAPE ( &
       (/ '    ', '_ave', '_std', '_min', '_max', '_cnt', '_cav',     &
          '    ', '_x1 ', '_x2 ', '_min', '_max', '_cnt', '_csm' /),  &
       (/ SND_MAXLEN, IOMODE_MAX /) )

  ! ATTRIBUTE FLAGS FOR RESTART HANDLING
  INTEGER, PARAMETER, PUBLIC :: AF_RST_NONE = AF_NONE ! NOT REQUIRED
  INTEGER, PARAMETER, PUBLIC :: AF_RST_CMP  = 1       ! COMPARE
  INTEGER, PARAMETER, PUBLIC :: AF_RST_INP  = 2       ! INPUT

  ! netCDF IDs
  INTEGER, PARAMETER, PUBLIC :: NC_ID_UNDEF = -1

  ! STRING LENGTHS
  ! op_bk_20130903+
#ifndef __ICON__
  INTEGER, PARAMETER, PUBLIC :: STRLEN_CHANNEL = 2*STRLEN_SHORT + 1
#else
  INTEGER, PARAMETER, PUBLIC :: STRLEN_CHANNEL = 128  
#endif
  ! op_bk_20130903-
  INTEGER, PARAMETER, PUBLIC :: STRLEN_OBJECT  = 2*STRLEN_MEDIUM + 1 + 4

  ! ====================================================================
  ! CHANNEL OBJECTS
  ! ====================================================================
  ! op_pj_20100827+
  TYPE t_chaobj_cpl
     CHARACTER(LEN=STRLEN_CHANNEL) :: cha  = ''
     CHARACTER(LEN=STRLEN_OBJECT ) :: obj  = ''
  END TYPE t_chaobj_cpl
  ! op_pj_20100827-

  ! MEMORY MANAGEMENT
  TYPE t_channel_object_mem
     !
     ! MEMORY USAGE
     INTEGER(I8) :: usage       = 0_I8  ! primary data section
     INTEGER(I8) :: usage_2nd   = 0_I8  ! secondary data section
     !
     ! FLAGS FOR INTERNAL USE
     LOGICAL     :: lalloc  = .FALSE.   ! AUTOMATIC MEMORY ALLOCATION
  END TYPE t_channel_object_mem

  ! I/O MANAGEMENT
  TYPE t_channel_object_io
     !
     ! RESTART HANDLING
     LOGICAL :: lrestart      = .FALSE. ! OUTPUT TO RESTART FILE ?
     LOGICAL :: lignore       = .FALSE. ! IGNORE lrestreq ?
     !
     ! OUTPUT FLAGS
     LOGICAL, DIMENSION(SND_MAXLEN) :: lout = .FALSE.
     ! SPECIAL FOR CONDITIONAL COUNTER (CNT) / AVAERAGE (CAV)
     REAL(DP), DIMENSION(2)         :: range = &
          (/ -HUGE(0.0_DP), HUGE(0.0_DP) /)
     !
  END TYPE t_channel_object_io

  ! netCDF I/O (internal use only !)
  TYPE t_channel_object_netcdf
     ! variable ID
     INTEGER                        :: varid        = NC_ID_UNDEF
     ! dimension IDs
     INTEGER                        :: dimid(IRANK) = NC_ID_UNDEF
     ! IDs OF SECONDARY VARIABLES
     INTEGER, DIMENSION(:), POINTER :: svarid => NULL()
  END TYPE t_channel_object_netcdf

  ! +++ ADD OTHER OUTPUT FORMATS HERE

  TYPE t_channel_object_int
     !
     ! OUTPUT AND RESTART
     LOGICAL :: lout = .FALSE.     ! ANY OUTPUT ?
     LOGICAL :: lrst = .FALSE.     ! ANY RESTART ?
     LOGICAL :: lign = .FALSE.     ! IGNORE lrestreq ?
     LOGICAL :: lref = .FALSE.     ! IS reference ?
     ! EXPORT DATA ?
     LOGICAL, DIMENSION(SND_MAXLEN, IOMODE_MAX) :: lexp = .FALSE.
     ! ... MEMORY MANAGEMENT FOR PRIMARY AND SECONDARY DATA
     INTEGER, DIMENSION(SND_MAXLEN) :: i2nd = 0  ! INDEX IN 2ndary DATA
     INTEGER :: n2nd = 0  ! SECONDARY DATA DIMENSION
     !
     ! MISC
     ! - FIELD HAS BEEN SET FROM RESTART FILE
     LOGICAL, DIMENSION(SND_MAXLEN) :: lrestart_read = .FALSE.
     !
     ! netCDF
     TYPE(t_channel_object_netcdf), DIMENSION(IOMODE_MAX) :: netcdf
     ! 
     ! +++ ADD OTHER OUTPUT FORMATS HERE
     !
  END TYPE t_channel_object_int

  TYPE t_channel_object
     CHARACTER(LEN=STRLEN_OBJECT)     :: name  = ''     ! NAME
     TYPE(t_attribute_list), POINTER  :: att => NULL()  ! OBJECT ATTRIBUTES
     TYPE(t_representation), POINTER  :: repr => NULL() ! REPRESENTATION
     TYPE(t_channel_object_mem)       :: memory         ! MEMORY MANAGEMENT
     TYPE(t_channel_object_io)        :: io             ! I/O MANAGEMENT
     ! ABSOLUTELY REQUIRED IN RESTART FILE ?
     LOGICAL                          :: lrestreq = .FALSE.
     TYPE(t_channel_object_int)       :: int            ! FOR INTERNAL USE
     REAL(DP), DIMENSION(:,:,:,:),     POINTER :: DATA => NULL() ! DATA
     TYPE(PTR_4D_ARRAY), DIMENSION(:), POINTER :: sdat => NULL() ! 2ndary DATA
     ! mz_ab_20090921+
     ! POINTER TO REGION WITHOUT BOUNDARIES FOR I/O
     REAL(DP), DIMENSION(:,:,:,:),     POINTER :: ioptr => NULL() 
     ! mz_ab_20090921-
  END TYPE t_channel_object

  TYPE t_channel_object_list
!     PRIVATE
     TYPE(t_channel_object)               :: this
     TYPE(t_channel_object_list), POINTER :: next => NULL()
  END TYPE t_channel_object_list
  ! ====================================================================

  ! ====================================================================
  ! CHANNELS
  ! ====================================================================
  TYPE t_channel_def
     !
     ! DEFAULT REPRESENTATION
     INTEGER                    :: reprid   = REPR_UNDEF
     ! DEFAULT: NOT REQUIRED IN RESTART
     LOGICAL                    :: lrestreq = .FALSE.
     ! DEFAULT OUTPUT FLAGS
     TYPE (t_channel_object_io) :: io
     !
  END TYPE t_channel_def

  ! netCDF I/O (internal use only !)
  TYPE t_channel_netcdf
     INTEGER  :: fileID     = NC_ID_UNDEF
     INTEGER  :: dimid_time = NC_ID_UNDEF
  END TYPE t_channel_netcdf

  ! +++ ADD OTHER OUTPUT FORMATS HERE

  TYPE t_channel_io
     ! OUTPUT FILE TYPE
     INTEGER, DIMENSION(IOMODE_MAX)  :: ftype  = FTYPE_DEFAULT
     ! NO. OF TIME STEPS PER FILE (OUT)
     INTEGER                         :: ntpf   = 0
     !
  END TYPE t_channel_io

  TYPE t_channel_int
     !
     ! OUTPUT AND RESTART
     LOGICAL  :: lout      = .FALSE.         ! OUTPUT ? ANY OBJECT ?
     LOGICAL  :: lrst      = .FALSE.         ! RESTART ? ANY OBJECT ?
     LOGICAL  :: lrestreq  = .FALSE.         ! ANY OBJECT REQUIRED IN RESTART
     LOGICAL  :: lign      = .FALSE.         ! IGNORE ALL (!) lrestreq ?
     ! - OUTPUT FILENAME : <EXPERIMENT (15)>_YYYYMMDD_HHMM_<CHANNEL>.<ext>
     ! - RESTART FILENAME: restart_<CHANNEL>.<ext>
     CHARACTER(LEN=30+STRLEN_CHANNEL+4), DIMENSION(IOMODE_MAX) :: fname = ''
     !
     ! TIMER
     LOGICAL  :: lout_now   = .FALSE.        ! TIME MANAGER
     INTEGER  :: ntpfcnt    = 0              ! COUNTER OF TIME STEPS PER FILE
     LOGICAL  :: lnew_file  = .FALSE.        ! OPEN NEW FILE ?
     REAL(DP) :: tslo       = 0.0_DP         ! time [s] since last output
     ! FORCE OUTPUT (to be set by set_channel_output)
     LOGICAL  :: lforce_out  = .FALSE.
     ! SUPPRESS OUTPUT (to be set by set_channel_output)
     LOGICAL  :: lsuppr_out  = .FALSE.
     ! FORCE NEW FILE
     LOGICAL  :: lforce_newfile = .FALSE.
     !
     ! op_pj_20100616+
     ! RE-INITIALIZE SECONDARY DATA
     LOGICAL :: l2ndreinit1 = .TRUE.
     LOGICAL :: l2ndreinit2 = .FALSE.
     ! op_pj_20100616-
     !
     ! netCDF
     TYPE(t_channel_netcdf), DIMENSION(IOMODE_MAX) :: netcdf
     !
     ! +++ ADD OTHER OUTPUT FORMATS HERE
     !
  END TYPE t_channel_int

  TYPE t_channel
     ! IDENTIFICATION
     CHARACTER(LEN=STRLEN_CHANNEL)    :: name  = ''       ! NAME
     INTEGER                          :: id    = 0        ! ID
     TYPE(t_attribute_list), POINTER  :: att   => NULL()  ! CHANNEL ATTRIBUTES
     TYPE(t_channel_object_mem)       :: memory           ! MEMORY MANAGEMENT
     TYPE(t_channel_io)               :: io               ! I/O
     TYPE(t_channel_def)              :: default          ! OBJECT DEFAULTS
     ! INTERNAL
     TYPE(t_channel_int)              :: int              ! FOR INTERNAL USE
     !
     ! CHANNEL OBJECTS
     TYPE(t_channel_object_list), POINTER :: list => NULL()
  END TYPE t_channel
  
  TYPE t_channel_list
!     PRIVATE
     TYPE(t_channel)               :: this
     TYPE(t_channel_list), POINTER :: next => NULL()
  END TYPE t_channel_list
  ! ====================================================================

  ! ====================================================================
  ! GLOBAL ATTRIBUTES COMMON TO ALL CHANNELS
  TYPE(t_attribute_list), POINTER, SAVE, PUBLIC :: GATT
  ! CONCT. LIST OF CHANNELS / OBJECTS
  TYPE(t_channel_list),   POINTER, SAVE, PUBLIC :: GCHANNELLIST  => NULL()
  INTEGER,                         SAVE, PUBLIC :: NCHANNEL = 0
  ! CHANNELS ALRADY FIXATED
  LOGICAL,                         SAVE         :: LFIXATE = .FALSE.
  ! ====================================================================

  ! ====================================================================
  ! FOR CTRL-NAMELIST
  TYPE t_channel_object_ctrl
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname = ''  ! CHANNEL NAME
     CHARACTER(LEN=STRLEN_OBJECT)  :: oname = ''  ! OBJECT NAME (TRACER !)
     TYPE(t_channel_object_io)     :: io
  END TYPE t_channel_object_ctrl
  !
  ! op_bk_20130905+
  ! Workaround for gfortran < 4.8.1, wrong handling of namelist input
  ! with more than one "sub-struct"
#ifdef __GFORTRAN__
  TYPE t_channel_ctrl_gf
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname = ''   ! CHANNEL NAME
     !
     ! (t_channel_io) cio
     ! OUTPUT FILE TYPE
     INTEGER, DIMENSION(IOMODE_MAX)  :: ftype  = FTYPE_DEFAULT
     ! NO. OF TIME STEPS PER FILE (OUT)
     INTEGER                         :: ntpf   = 0
     !
     ! (t_channel_object_io)
     !
     ! RESTART HANDLING
     LOGICAL :: lrestart      = .FALSE. ! OUTPUT TO RESTART FILE ?
     LOGICAL :: lignore       = .FALSE. ! IGNORE lrestreq ?
     !
     ! OUTPUT FLAGS
     LOGICAL, DIMENSION(SND_MAXLEN) :: lout = .FALSE.
     ! SPECIAL FOR CONDITIONAL COUNTER (CNT) / AVAERAGE (CAV)
     REAL(DP), DIMENSION(2)         :: range = &
          (/ -HUGE(0.0_DP), HUGE(0.0_DP) /)
     !
  END TYPE t_channel_ctrl_gf
#endif
  ! op_bk_20130905-
  TYPE t_channel_ctrl
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname = ''   ! CHANNEL NAME
     TYPE(t_channel_io)            :: cio
     TYPE(t_channel_object_io)     :: oio
  END TYPE t_channel_ctrl
  !
  TYPE t_channel_new
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname = ''   ! CHANNEL NAME
  END TYPE t_channel_new
  !
  TYPE t_channel_ref_new
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname1 = ''  ! CHANNEL NAME  SOURCE
     CHARACTER(LEN=STRLEN_OBJECT)  :: oname1 = ''  ! OBJECT NAME SOURCE
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname2 = ''  ! CHANNEL NAME  NEW
     CHARACTER(LEN=STRLEN_OBJECT)  :: oname2 = ''  ! OBJECT NAME NEW
  END TYPE t_channel_ref_new
  !
  CHARACTER(LEN=15), PUBLIC, SAVE :: EXP_NAME = ''
  !
  LOGICAL, PUBLIC, SAVE :: L_FLUSH_IOBUFFER = .TRUE.
  ! op_pj_20110803+
  ! 0: only error messages, 1: initialize and finalize, 2: 1 and time loop
  INTEGER, PUBLIC, SAVE :: I_VERBOSE_LEVEL  = 2
  ! op_pj_20110803-
  !
  ! op_bk_20130905+
  ! Workaround for gfortran < 4.8.1, wrong handling of namelist input
  ! with more than one "sub-struct"
#ifdef __GFORTRAN__
  TYPE(t_channel_ctrl_gf),                             SAVE  :: OUT_DEFAULT
#else
  TYPE(t_channel_ctrl),                                SAVE  :: OUT_DEFAULT
#endif
  ! op_bk_20130905-
  ! mz_pj_20080118+
  INTEGER,                    DIMENSION(FTYPE_MAXIMUM),SAVE  :: OUT_PREC = &
       (/1, 1, 1, 1, 1, 1/)
  ! mz_pj_20080118-
  INTEGER, PARAMETER  :: NMAXCHANNELS =  500
  ! op_bk_20130905+
  ! Workaround for gfortran < 4.8.1, wrong handling of namelist input
  ! with more than one "sub-struct"
#ifdef __GFORTRAN__
  TYPE(t_channel_ctrl_gf),    DIMENSION(NMAXCHANNELS), SAVE  :: OUT_CHANNEL
#else
  TYPE(t_channel_ctrl),       DIMENSION(NMAXCHANNELS), SAVE  :: OUT_CHANNEL
#endif
  ! op_bk_20130905-
  INTEGER, PARAMETER  :: NMAXOBJECTS  = 1000
  TYPE(t_channel_object_ctrl), DIMENSION(NMAXOBJECTS), SAVE  :: OUT_OBJECT
  INTEGER, PARAMETER  :: NMAXADDCHANNEL = 50
  TYPE(t_channel_new),       DIMENSION(NMAXADDCHANNEL),SAVE  :: ADD_CHANNEL
  INTEGER, PARAMETER  :: NMAXADDREF    = 500
  TYPE(t_channel_ref_new),    DIMENSION(NMAXADDREF),   SAVE  :: ADD_REF
  !
  PUBLIC :: NMAXCHANNELS, NMAXOBJECTS, NMAXADDCHANNEL, NMAXADDREF 
  PUBLIC :: OUT_DEFAULT, OUT_CHANNEL, OUT_OBJECT, ADD_CHANNEL, ADD_REF
  PUBLIC :: OUT_PREC ! mz_pj_20080118
  PUBLIC :: t_channel_list, t_channel                                 &
       , t_channel_int, t_channel_def, t_channel_io, t_channel_netcdf &
       , t_channel_new, t_channel_ref_new, t_channel_ctrl
  PUBLIC :: t_channel_object_list, t_channel_object                   &
       , t_channel_object_int, t_channel_object_io, t_channel_object_netcdf &
       , t_channel_object_mem, t_channel_object_ctrl
  PUBLIC :: t_chaobj_cpl ! op_pj_20100827
  ! ====================================================================

  ! PUBLIC INTERFACES ==================================================
  !
  ! - ATTRIBUTES
  INTERFACE new_attribute
     MODULE PROCEDURE add_attribute
     MODULE PROCEDURE new_global_attribute
     MODULE PROCEDURE new_channel_attribute
     MODULE PROCEDURE new_channel_object_attribute
  END INTERFACE
  PUBLIC :: new_attribute    ! add   GLOBAL, CHANNEL, or OBJECT ATTRIBUTES
  !
  INTERFACE get_attribute
     MODULE PROCEDURE return_attribute
     MODULE PROCEDURE get_global_attribute
     MODULE PROCEDURE get_channel_attribute
     MODULE PROCEDURE get_channel_object_attribute
  END INTERFACE
  PUBLIC :: get_attribute    ! add   GLOBAL, CHANNEL, or OBJECT ATTRIBUTES
  !
  INTERFACE write_attribute
     MODULE PROCEDURE write_global_attributes
     MODULE PROCEDURE write_channel_attributes
     MODULE PROCEDURE write_channel_object_attributes
  END INTERFACE
  PUBLIC :: write_attribute ! write GLOBAL, CHANNEL, or OBJECT ATTRIBUTES
  !
  ! - CHANNELS
  PUBLIC :: new_channel           ! define new CHANNEL
  INTERFACE write_channel
     MODULE PROCEDURE write_channel_by_name
     MODULE PROCEDURE write_channel_all
  END INTERFACE
  PUBLIC :: write_channel
  PUBLIC :: get_channel_info
  PUBLIC :: get_channel_name
  PUBLIC :: set_channel_output           ! force/suppress output
  PUBLIC :: set_channel_newfile          ! force new output file
  !
  ! - CHANNEL OBJECTS
  PUBLIC :: new_channel_object           ! define new channel OBJECT
  PUBLIC :: get_channel_object           ! get POINTER to OBJECT-memory
  PUBLIC :: get_channel_object_info      ! selected information
  PUBLIC :: new_channel_object_reference ! channel memory reference
  PUBLIC :: set_channel_object_restreq   ! set lrestreq = .TRUE.
  ! op_pj_20100827+
  PUBLIC :: get_channel_object_dimvar    ! get dimension variables and units
  ! op_pj_20100827-
  !
  ! - OVERALL SETUP
  PUBLIC :: main_channel_read_nml_ctrl
  !
  PUBLIC :: fixate_channels              ! consistency check and 2ndary memory
  !                                      ! (once after initialization)
  PUBLIC :: trigger_channel_output       ! trigger next output
  !                                      ! (once every time step)
  PUBLIC :: update_channels              ! update 2ndary variables
  !                                      ! (-> see flag)
  PUBLIC :: clean_channels               ! cleanup memory
  !
  ! ====================================================================

  ! PRIVATE INTERFACES =================================================
  !
  ! - ATTRIBUTES
  ! NONE
  !
  ! - CHANNELS
  INTERFACE loc_channel                    ! locate pointer to channel
     MODULE PROCEDURE loc_channel_by_name  ! locate pointer to channel by NAME
     MODULE PROCEDURE loc_channel_by_id    ! locate pointer to channel by ID
  END INTERFACE
  !PRIVATE :: write_channel_by_ptr
  !
  INTERFACE write_channel_object
     MODULE PROCEDURE write_channel_object_by_name
     MODULE PROCEDURE write_channel_object_by_ptr
  END INTERFACE
  !
  ! - CHANNEL OBJECTS
  !PRIVATE :: loc_channel_object ! locate pinter to channel-object (by NAME)
  !PRIVATE :: write_channel_object
  ! ====================================================================

CONTAINS

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES
  ! -------------------------------------------------------------------

  ! **********************************************************************
  ! ATTRIBUTES
  ! **********************************************************************

  ! -------------------------------------------------------------------
  SUBROUTINE new_global_attribute(status &
       , ganame, i, r, c, loverwrite, iflag)

    USE messy_main_channel_attributes, ONLY: add_attribute

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                   INTENT(OUT)          :: status
    CHARACTER(LEN=*),          INTENT(IN)           :: ganame
    INTEGER,                   INTENT(IN), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(IN), OPTIONAL :: c
    REAL(DP),                  INTENT(IN), OPTIONAL :: r
    LOGICAL,                   INTENT(IN), OPTIONAL :: loverwrite
    INTEGER,                   INTENT(IN), OPTIONAL :: iflag

    CALL add_attribute(status, GATT, TRIM(ganame), i, c, r, loverwrite, iflag)

  END SUBROUTINE new_global_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE new_channel_attribute(status, cname &
       , caname, i, r, c, loverwrite, iflag)

    USE messy_main_channel_attributes, ONLY: add_attribute

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                   INTENT(OUT)          :: status
    CHARACTER(LEN=*),          INTENT(IN)           :: cname
    CHARACTER(LEN=*),          INTENT(IN)           :: caname
    INTEGER,                   INTENT(IN), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(IN), OPTIONAL :: c
    REAL(DP),                  INTENT(IN), OPTIONAL :: r
    LOGICAL,                   INTENT(IN), OPTIONAL :: loverwrite
    INTEGER,                   INTENT(IN), OPTIONAL :: iflag 

    ! LOCAL
    TYPE(t_channel), POINTER :: channel => NULL()

    CALL loc_channel(status, GCHANNELLIST, cname, channel)
    IF (status /= 0) RETURN

    CALL add_attribute(status, channel%att, TRIM(caname), i, c, r &
         , loverwrite, iflag)

  END SUBROUTINE new_channel_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE new_channel_object_attribute(status, cname, oname &
       , oaname, i, r, c, loverwrite, iflag)

    USE messy_main_channel_attributes, ONLY: add_attribute

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                   INTENT(OUT)          :: status
    CHARACTER(LEN=*),          INTENT(IN)           :: cname
    CHARACTER(LEN=*),          INTENT(IN)           :: oname
    CHARACTER(LEN=*),          INTENT(IN)           :: oaname
    INTEGER,                   INTENT(IN), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(IN), OPTIONAL :: c
    REAL(DP),                  INTENT(IN), OPTIONAL :: r
    LOGICAL,                   INTENT(IN), OPTIONAL :: loverwrite
    INTEGER,                   INTENT(IN), OPTIONAL :: iflag 

    ! LOCAL
    TYPE(t_channel_object), POINTER :: object => NULL()

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object)
    IF (status /= 0) RETURN

    CALL add_attribute(status, object%att, TRIM(oaname), i, c, r &
         , loverwrite, iflag)

  END SUBROUTINE new_channel_object_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_global_attributes(status)

    USE messy_main_channel_attributes, ONLY: write_attribute

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status

    WRITE(*,*) '=== GLOBAL ATTRIBUTES: ============================='
    CALL write_attribute(status, GATT)
    IF (status /= 0) RETURN
    WRITE(*,*) '===================================================='

  END SUBROUTINE write_global_attributes
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_attributes(status, cname)

    USE messy_main_channel_attributes, ONLY: write_attribute

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: cname

    ! LOCAL
    TYPE(t_channel), POINTER :: channel => NULL()

    CALL loc_channel(status, GCHANNELLIST, cname, channel)
    IF (status /= 0) RETURN

    CALL write_attribute(status, channel%att)    

  END SUBROUTINE write_channel_attributes
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_object_attributes(status, cname, oname)

    USE messy_main_channel_attributes, ONLY: write_attribute

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: cname
    CHARACTER(LEN=*), INTENT(IN)  :: oname

    ! LOCAL
    TYPE(t_channel_object), POINTER :: object => NULL()

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object)
    IF (status /= 0) RETURN

    CALL write_attribute(status, object%att)    

  END SUBROUTINE write_channel_object_attributes
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_global_attribute(status &
       , ganame, i, r, c, iflag)

    USE messy_main_channel_attributes, ONLY: return_attribute

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                   INTENT(OUT)           :: status
    CHARACTER(LEN=*),          INTENT(IN)            :: ganame
    INTEGER,                   INTENT(OUT), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(OUT), OPTIONAL :: c
    REAL(DP),                  INTENT(OUT), OPTIONAL :: r
    INTEGER,                   INTENT(OUT), OPTIONAL :: iflag

    CALL return_attribute(status, GATT, TRIM(ganame), i, c, r, iflag)

  END SUBROUTINE get_global_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_attribute(status, cname &
       , caname, i, c, r, iflag)

    USE messy_main_channel_attributes, ONLY: return_attribute

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                   INTENT(OUT)           :: status
    CHARACTER(LEN=*),          INTENT(IN)            :: cname
    CHARACTER(LEN=*),          INTENT(IN)            :: caname
    INTEGER,                   INTENT(OUT), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(OUT), OPTIONAL :: c
    REAL(DP),                  INTENT(OUT), OPTIONAL :: r
    INTEGER,                   INTENT(OUT), OPTIONAL :: iflag 

    ! LOCAL
    TYPE(t_channel), POINTER :: channel => NULL()

    CALL loc_channel(status, GCHANNELLIST, cname, channel)
    IF (status /= 0) RETURN

    CALL return_attribute(status, channel%att, TRIM(caname), i, c, r, iflag)

  END SUBROUTINE get_channel_attribute
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_object_attribute(status, cname, oname &
       , oaname, i, r, c, iflag)

    USE messy_main_channel_attributes, ONLY: return_attribute

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                   INTENT(OUT)          :: status
    CHARACTER(LEN=*),          INTENT(IN)           :: cname
    CHARACTER(LEN=*),          INTENT(IN)           :: oname
    CHARACTER(LEN=*),          INTENT(IN)           :: oaname
    INTEGER,                   INTENT(OUT), OPTIONAL :: i
    CHARACTER(LEN=*),          INTENT(OUT), OPTIONAL :: c
    REAL(DP),                  INTENT(OUT), OPTIONAL :: r
    INTEGER,                   INTENT(OUT), OPTIONAL :: iflag 

    ! LOCAL
    TYPE(t_channel_object), POINTER :: object => NULL()

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object)
    IF (status /= 0) RETURN

    CALL return_attribute(status, object%att, TRIM(oaname), i, c, r, iflag)

  END SUBROUTINE get_channel_object_attribute
  ! -------------------------------------------------------------------

  ! **********************************************************************
  ! CHANNELS
  ! **********************************************************************

  ! -------------------------------------------------------------------
  SUBROUTINE new_channel(status, cname, reprid, lrestreq)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)          :: status
    CHARACTER(LEN=*), INTENT(IN)           :: cname    ! CHANNEL NAME
    INTEGER,          INTENT(IN), OPTIONAL :: reprid   ! REPRESENTATION ID
    LOGICAL,          INTENT(IN), OPTIONAL :: lrestreq ! REQUIRED IN RESTART

    ! LOCAL
    TYPE(t_channel_list), POINTER :: ai     => NULL()
    TYPE(t_channel_list), POINTER :: ae     => NULL()
    TYPE(t_channel),      POINTER :: channel => NULL()
    INTEGER                      :: zstat
    
    ! CHECKS
    IF (LFIXATE) THEN
       status = 3000
       RETURN
    END IF

    CALL loc_channel(zstat, GCHANNELLIST, cname, channel)
    IF (zstat /= 3003) THEN  ! CHANNEL (NAME) DOES NOT EXIST (IS OK HERE !)
       IF (zstat == 0) THEN  ! CHANNEL EXISTS ALREADY
          status = 3002      ! CHANNEL EXISTS ALREADY
       ELSE
          status = zstat    ! ERROR
       END IF
       RETURN
    END IF

    ! GOTO END OF LIST
    ai => GCHANNELLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       ae => ai
       ai => ai%next
    END DO

    ! ADD NEW
    ALLOCATE(ai)
    NULLIFY(ai%next)
    IF (.NOT. ASSOCIATED(GCHANNELLIST)) THEN
       GCHANNELLIST => ai              ! SET POINTER TO FIRST OBJECT
    ELSE
       ae%next => ai                   ! SET NEXT POINTER OF LAST OBJECT
       !                               ! TO NEW OBJECT
    END IF

    ! SET VALUES
    ai%this%name = TRIM(ADJUSTL(cname))

    ! SET DEFAULT REPRESENTATION
    IF (PRESENT(reprid))   ai%this%default%reprid   = reprid
    ! SET DEFAULT REQUIREMENT FOR RESTART
    IF (PRESENT(lrestreq)) ai%this%default%lrestreq = lrestreq

    ! COUNT AND SET ID
    NCHANNEL   = NCHANNEL + 1
    ai%this%id = NCHANNEL

    status = 0

  END SUBROUTINE new_channel
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_by_ptr(status, channel)

    USE messy_main_channel_attributes, ONLY: write_attribute

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,         INTENT(OUT) :: status
    TYPE(t_channel), POINTER     :: channel

    ! LOCAL
    TYPE(t_channel_object_list), POINTER :: le
    TYPE(t_channel_object),      POINTER :: object

    IF (.NOT. ASSOCIATED(channel)) THEN
       status = 3006 ! CHANNEL POINTER NOT ASSOCIATED
       RETURN
    END IF

    WRITE(*,*) '========================================================================'
    WRITE(*,*) ' NAME          : ',channel%name
    WRITE(*,*) ' ID            : ',channel%id
    WRITE(*,*) ' OUT-FILE-TYPE : ',channel%io%ftype(IOMODE_OUT)
    IF (channel%io%ntpf > 0) THEN                          ! mz_pj_20080905
       WRITE(*,*) ' STEPS PER FILE: ',channel%io%ntpf
    ELSE                                                   ! mz_pj_20080905
       WRITE(*,*) ' STEPS PER FILE: ','(event triggered)'  ! mz_pj_20080905
    ENDIF                                                  ! mz_pj_20080905
    WRITE(*,*) ' ANY OUTPUT  ? : ',channel%int%lout
    WRITE(*,*) ' RST-FILE-TYPE : ',channel%io%ftype(IOMODE_RST)
    WRITE(*,*) ' ANY RESTART ? : ',channel%int%lrst
    WRITE(*,*) ' MEMORY        : ',channel%memory%usage,' + ' &
         , channel%memory%usage_2nd
    WRITE(*,*) ' ATTRIBUTES    : '
    CALL write_attribute(status, channel%att)
    IF (status /= 0) RETURN

    WRITE(*,*) '------------------------------------------------------------------------'
    WRITE(*,'(1x,a24,1x,a7,1x,a13,1x,a4,1x,a20)') &
         '                        ','--RST--','---OUTPUT----','REPR','-------MEMORY-------'
    WRITE(*,'(1x,a24,1x,11(a1,1x))') &
         '                        ','R','R','I','U','I','A','S','M','M','C','C'
    WRITE(*,'(1x,a24,1x,11(a1,1x))') &
         '                        ','E','E','G','S','N','V','T','I','A','N','A'
    WRITE(*,'(1x,a24,1x,11(a1,1x))') &
         '                        ','S','S','N','E','S','E','D','N','X','T','V'
    WRITE(*,'(1x,a24,1x,11(a1,1x))') &
         '                        ','T','T','O','R',' ',' ',' ',' ',' ',' ',' '
    WRITE(*,'(1x,a24,1x,11(a1,1x))') &
         '                        ','A','R','R',' ',' ',' ',' ',' ',' ',' ',' '
    WRITE(*,'(1x,a24,1x,11(a1,1x))') &
         '                        ','R','E','E',' ',' ',' ',' ',' ',' ',' ',' '
    WRITE(*,'(1x,a24,1x,11(a1,1x),a4,1x,a1,1x,2(a8,1x))') &
         'NAME                    ','T','Q',' ',' ',' ',' ',' ',' ',' ',' ',' ','REPR','M','   MEM_1','  MEM_02'

    WRITE(*,*) '------------------------------------------------------------------------'

    le => channel%list
    object_loop: DO
       IF (.NOT. ASSOCIATED(le)) EXIT

       object => le%this

       CALL write_channel_object_by_ptr(status, object)
       IF (status /= 0) RETURN
       
       le => le%next
    END DO object_loop

    WRITE(*,*) '========================================================================'

    status = 0

  END SUBROUTINE write_channel_by_ptr
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_by_name(status, cname)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: cname

    ! LOCAL
    TYPE(t_channel), POINTER :: channel

    CALL loc_channel(status, GCHANNELLIST, cname, channel)
    IF (status /= 0) RETURN

    CALL write_channel_by_ptr(status, channel)

  END SUBROUTINE write_channel_by_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_all(status)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER, INTENT(OUT)  :: status

    ! LOCAL
    TYPE(t_channel_list),         POINTER :: ls
    TYPE(t_channel),              POINTER :: channel

    WRITE(*,*) '+++ CHANNELS: ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    ! EMPTY LIST ?
    IF (.NOT. ASSOCIATED(GCHANNELLIST)) THEN

       WRITE(*,*) '*** CHANNEL LIST IS EMPTY ***'

    ELSE

       ls => GCHANNELLIST
       channel_loop: DO
          IF (.NOT. ASSOCIATED(ls)) EXIT
          
          channel => ls%this
          
          CALL write_channel_by_ptr(status, channel)
          IF (status /= 0) RETURN
          
          ls => ls%next
       END DO channel_loop
    
    END IF

    WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    status = 0
    
  END SUBROUTINE write_channel_all
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_info(status, cname, LDIMS, LREPRS, ONAMES &
       , pick)

    USE messy_main_channel_dimensions,      ONLY: NDIM
    USE messy_main_channel_repr,            ONLY: NREP, IRANK

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, PRESENT, TRIM

    ! I/O
    INTEGER,               INTENT(OUT)          :: status
    CHARACTER(LEN=*),      INTENT(IN)           :: cname   ! CHANNEL NAME
    LOGICAL, DIMENSION(:), POINTER,    OPTIONAL :: LDIMS   ! DIMENSION FLAGS
    LOGICAL, DIMENSION(:), POINTER,    OPTIONAL :: LREPRS  ! REPR. FLAGS
    CHARACTER(LEN=STRLEN_OBJECT), DIMENSION(:), POINTER, OPTIONAL :: ONAMES
    CHARACTER(LEN=*),      INTENT(IN), OPTIONAL :: pick

    ! LOCAL
    TYPE(t_channel),             POINTER :: channel
    TYPE(t_channel_object_list), POINTER :: le
    TYPE(t_channel_object),      POINTER :: object
    INTEGER                              :: i
    INTEGER                              :: nvarcnt
    INTEGER                              :: ipick
    LOGICAL                              :: zl

    ! INIT
    IF (PRESENT(LDIMS)) THEN
       IF (ASSOCIATED(LDIMS)) DEALLOCATE(LDIMS)
       ALLOCATE(LDIMS(NDIM))
       LDIMS(:) = .FALSE.
    END IF
    !
    IF (PRESENT(LREPRS)) THEN
       IF (ASSOCIATED(LREPRS)) DEALLOCATE(LREPRS)
       ALLOCATE(LREPRS(NREP))
       LREPRS(:) = .FALSE.
    END IF
    !
    IF (PRESENT(ONAMES)) THEN
       IF (ASSOCIATED(ONAMES)) DEALLOCATE(ONAMES)
       NULLIFY(ONAMES)
    END IF
    !
    IF (PRESENT(pick)) THEN
       SELECT CASE(TRIM(ADJUSTL(pick)))
          CASE('restart')
             ipick = 0
          CASE('output')
             ipick = 1
          CASE('all')
             ipick = 2
          CASE DEFAULT
       END SELECT
    ELSE
       ipick = 2    ! DEFAULT: 'all'
    END IF

    CALL loc_channel(status, GCHANNELLIST, cname, channel)
    IF (status /= 0) RETURN

    nvarcnt = 0
    le => channel%list
    object_loop: DO
       IF (.NOT. ASSOCIATED(le)) EXIT
       
       object => le%this

       SELECT CASE(ipick)
          CASE(0) ! RESTART
             zl = object%int%lrst
          CASE(1) ! OUTPUT
             zl = object%int%lout
          CASE(2) ! ALL
             zl = .TRUE.
       END SELECT

       IF (PRESENT(LREPRS)) LREPRS(object%repr%id) = zl

       IF (PRESENT(LDIMS)) THEN
          DO i=1, IRANK
             IF (ASSOCIATED(object%repr%dim(i)%ptr)) &
                  LDIMS(object%repr%dim(i)%ptr%id) = zl
          END DO
       END IF

       IF (zl) nvarcnt = nvarcnt + 1

       le => le%next
    END DO object_loop

    IF (PRESENT(ONAMES)) THEN
       ALLOCATE(ONAMES(nvarcnt))
       ONAMES(:) = ''
       nvarcnt = 0
       le => channel%list
       object_loop2: DO
          IF (.NOT. ASSOCIATED(le)) EXIT

          object => le%this

          SELECT CASE(ipick)
          CASE(0) ! RESTART
             zl = object%int%lrst
          CASE(1) ! OUTPUT
             zl = object%int%lout
          CASE(2) ! ALL
             zl = .TRUE.
          END SELECT

          IF (zl) THEN
             nvarcnt = nvarcnt + 1
             ONAMES(nvarcnt) = TRIM(object%name)
          END IF

          le => le%next
       END DO object_loop2
    END IF

    status = 0

  END SUBROUTINE get_channel_info
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_name(status, id, cname)

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                       INTENT(OUT) :: status
    INTEGER,                       INTENT(IN)  :: id
    CHARACTER(LEN=STRLEN_CHANNEL), INTENT(OUT) :: cname

    ! LOCAL   
    TYPE(t_channel),                   POINTER :: channel  => NULL()

    ! INIT
    cname = ''

    CALL loc_channel_by_id(status, GCHANNELLIST, id, channel)
    IF (status /= 0) RETURN

    cname = TRIM(channel%name)

  END SUBROUTINE get_channel_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_channel_output(status, cname, lout)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: cname
    LOGICAL,          INTENT(IN)  :: lout

    ! LOCAL   
    TYPE(t_channel),                   POINTER :: channel  => NULL()

    CALL loc_channel(status, GCHANNELLIST, cname, channel)
    IF (status /= 0) RETURN

    IF (lout) THEN
       channel%int%lforce_out = .TRUE.
    ELSE
       channel%int%lsuppr_out = .TRUE.
    END IF

  END SUBROUTINE set_channel_output
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_channel_newfile(status, cname, lnew)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: cname
    LOGICAL,          INTENT(IN)  :: lnew

    ! LOCAL   
    TYPE(t_channel),                   POINTER :: channel  => NULL()

    CALL loc_channel(status, GCHANNELLIST, cname, channel)
    IF (status /= 0) RETURN

    channel%int%lforce_newfile = lnew

  END SUBROUTINE set_channel_newfile
  ! -------------------------------------------------------------------

  ! **********************************************************************
  ! CHANNEL OBJECTS
  ! **********************************************************************

  ! -------------------------------------------------------------------
  SUBROUTINE new_channel_object(status, cname, oname &
       , p0, p1, p2, p3, p4, mem                     &
       , reprid                                      &
       , lrestreq)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE, PRESENT, TRIM

    ! I/O
    INTEGER,                      INTENT(OUT)          :: status
    CHARACTER(LEN=*),             INTENT(IN)           :: cname ! CHANNEL NAME
    CHARACTER(LEN=*),             INTENT(IN)           :: oname ! OBJECT NAME
    REAL(DP),                     POINTER,    OPTIONAL :: p0    ! POINTER ...
    REAL(DP), DIMENSION(:),       POINTER,    OPTIONAL :: p1    ! ...
    REAL(DP), DIMENSION(:,:),     POINTER,    OPTIONAL :: p2    ! ...
    REAL(DP), DIMENSION(:,:,:),   POINTER,    OPTIONAL :: p3    ! ...
    REAL(DP), DIMENSION(:,:,:,:), POINTER,    OPTIONAL :: p4    ! ... TO MEMORY
    REAL(DP), DIMENSION(:,:,:,:), POINTER,    OPTIONAL :: mem   ! EXT. MEMORY
    !
    ! REPRESENTATION ID
    INTEGER,                      INTENT(IN), OPTIONAL :: reprid
    ! REQ. IN RESTART FILE ?
    LOGICAL,                      INTENT(IN), OPTIONAL :: lrestreq

    ! LOCAL
    TYPE(t_channel),             POINTER :: channel  => NULL()
    TYPE(t_channel_object),      POINTER :: object => NULL()
    TYPE(t_channel_object_list), POINTER :: ai      => NULL()
    TYPE(t_channel_object_list), POINTER :: ae      => NULL()
    INTEGER                              :: zstat   ! local status
    INTEGER                              :: i1,i2,i3,i4
    LOGICAL                              :: lshape

    ! CHECKS
    IF (LFIXATE) THEN
       status = 3000
       RETURN
    END IF

    CALL loc_channel(status, GCHANNELLIST, TRIM(cname), channel)
    IF (status /= 0) RETURN

    CALL get_channel_object(zstat, TRIM(cname), TRIM(oname))
    IF (zstat == 0) THEN
       status = 3102    ! CHANNEL OBJECT EXISTS ALREADY
       RETURN
    END IF

    ! GOTO END OF LIST
    ai => channel%list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       ae => ai
       ai => ai%next
    END DO    

    ! ADD NEW OBJECT TO LIST
    ALLOCATE(ai)
    NULLIFY(ai%next)
    IF (.NOT. ASSOCIATED(channel%list)) THEN
       channel%list => ai               ! SET POINTER TO FIRST OBJECT
    ELSE
       ae%next     => ai                ! SET NEXT POINTER OF LAST OBJECT
       !                                ! TO NEW OBJECT
    END IF

    object => ai%this

    ! SET NAME
    object%name = TRIM(oname)
    ! SET REPRESENTATION
    IF (PRESENT(reprid)) THEN
       CALL get_representation(status, reprid, object%repr)
       IF (status /= 0) RETURN
    ELSE
       ! CHECK DEFAULT REPRESENTATION OF CHANNEL
       IF (channel%default%reprid == REPR_UNDEF) THEN
          status = 3104 ! CHANNEL DEFAULT REPRESENTATION IS UNDEFINED
          RETURN
       ELSE
          CALL get_representation(status, channel%default%reprid, object%repr)
          IF (status /= 0) RETURN
       END IF
    END IF

    ! MEMORY MANAGEMENT
    object%memory%lalloc = (.NOT. PRESENT(mem))
    IF (object%memory%lalloc) THEN
       i1 = object%repr%ldimlen(1)
       i2 = object%repr%ldimlen(2)
       i3 = object%repr%ldimlen(3)
       i4 = object%repr%ldimlen(4)
       ALLOCATE(object%data(i1,i2,i3,i4), STAT=zstat)
       IF (zstat /= 0) THEN
          status = 1000 ! MEMORY ALLOCATION FAILED
          RETURN
       END IF
       object%data(:,:,:,:) = 0.0_DP
       ! OBJECT
       object%memory%usage = INT(SIZE(object%data), I8)
       ! CHANNEL
       channel%memory%usage = channel%memory%usage + object%memory%usage
    ELSE
       ! CHECK SIZE
       i1 = SIZE(mem, 1)
       i2 = SIZE(mem, 2)
       i3 = SIZE(mem, 3)
       i4 = SIZE(mem, 4)
       lshape = (i1 == object%repr%ldimlen(1)) .AND. &
                (i2 == object%repr%ldimlen(2)) .AND. &
                (i3 == object%repr%ldimlen(3)) .AND. &
                (i4 == object%repr%ldimlen(4))
       IF (.NOT. lshape) THEN
          status = 3105  ! SHAPE OF MEMORY NOT CONFORM WITH REPRESENTATION
          RETURN
       ELSE
          object%data => mem(:,:,:,:)
       END IF
    END IF

    ! SET DEFAULTS
    ! AVAILABILITY IN RESTART FILE IS MANDATORY ?
    IF (PRESENT(lrestreq)) THEN
       ! USER DEFINED PER OBJECT
       object%lrestreq = lrestreq
    ELSE
       ! USER DEFINED DEFAULT PER CHANNEL (DEFAULT: .FALSE.)
       object%lrestreq = channel%default%lrestreq
    END IF
    !
    object%io = channel%default%io

    ! mz_ab_20090921+
    ! POINTER TO MEMORY FOR I/O (WITHIN POTENTIAL BOUNDARIES)
    IF (object%repr%bounds%lbounds) THEN
       ! IN CASE OF BOUNDARIES
       object%ioptr => object%data( &
            object%repr%bounds%nbounds(1)+1: &
            object%repr%ldimlen(1)-object%repr%bounds%nbounds(1), &
            object%repr%bounds%nbounds(2)+1: &
            object%repr%ldimlen(2)-object%repr%bounds%nbounds(2), &
            object%repr%bounds%nbounds(3)+1: &
            object%repr%ldimlen(3)-object%repr%bounds%nbounds(3), &
            object%repr%bounds%nbounds(4)+1: &
            object%repr%ldimlen(4)-object%repr%bounds%nbounds(4)  &
            )
    ELSE
       object%ioptr => object%data
    ENDIF
    ! mz_ab_20090921-

    ! SET POINTER
    CALL get_channel_object(status, TRIM(cname), TRIM(oname) &
         , p0, p1, p2, p3, p4)

  END SUBROUTINE new_channel_object
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_object(status, cname, oname, p0, p1, p2, p3, p4 &
       , linner) ! mz_ab_20100309

    USE messy_main_channel_repr,  ONLY: repr_getptr

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                      INTENT(OUT)          :: status
    ! CHANNEL NAME
    CHARACTER(LEN=*),             INTENT(IN)           :: cname
    ! OBJECT NAME
    CHARACTER(LEN=*),             INTENT(IN)           :: oname
    REAL(DP),                     POINTER,    OPTIONAL :: p0
    REAL(DP), DIMENSION(:),       POINTER,    OPTIONAL :: p1
    REAL(DP), DIMENSION(:,:),     POINTER,    OPTIONAL :: p2
    REAL(DP), DIMENSION(:,:,:),   POINTER,    OPTIONAL :: p3
    REAL(DP), DIMENSION(:,:,:,:), POINTER,    OPTIONAL :: p4
    ! mz_ab_20100211+
    LOGICAL,                      INTENT(IN), OPTIONAL :: linner
    ! mz_ab_20100211-

    ! LOCAL
    TYPE(t_channel),              POINTER :: channel  => NULL()
    TYPE(t_channel_object),       POINTER :: object => NULL()

    ! mz_ab_20100211+
    LOGICAL :: zlinner 

    IF (PRESENT(linner)) THEN
       zlinner = linner
    ELSE
       zlinner = .FALSE.   ! default
    ENDIF   
    ! mz_ab_20100211-

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object)
    IF (status /= 0) RETURN

    ! mz_ab_20100211+ 
    ! return only inner part of object with bounds
    ! repr_getptr always returns the pointer p0...p4 to the entire input
    ! pointer (3rd argument). This means zlinner has no meaning if
    ! the reprsentation (2nd argument) has no boundaries. In case 
    ! boundaries are present
    ! - and linner is requested: We use the already present ioptr 
    !                            to the inner section.
    ! - and linner is not requested (i.e. we want the field incl bounds):
    !                            We use the full data pointer and force
    !                            repr_getptr to remap the indices to 
    !                            e.g. (-1,...)
    
    !!$CALL repr_getptr(status, object%repr, object%data, p0, p1, p2, p3, p4)

    IF (.NOT.zlinner) THEN
       CALL repr_getptr(status, object%repr, object%data,  p0, p1, p2, p3, p4 &
                      , .FALSE.)
    ELSE
       CALL repr_getptr(status, object%repr, object%ioptr, p0, p1, p2, p3, p4 &
                      , .TRUE.)
    END IF
    ! mz_ab_20100211-

    IF (status /= 0) RETURN

  END SUBROUTINE get_channel_object
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_object_info(status, cname, oname &
       , lrestart_read , reprid, axis, nbounds) ! um_ak_20091016 axis added
!                                               ! op_pj_20100121 nbounds added

    IMPLICIT NONE

    INTRINSIC :: PRESENT, TRIM

    ! I/O
    INTEGER,                      INTENT(OUT)           :: status
    ! CHANNEL NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: cname
    ! OBJECT NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: oname
    LOGICAL,                      INTENT(OUT), OPTIONAL :: lrestart_read
    INTEGER,                      INTENT(OUT), OPTIONAL :: reprid
    CHARACTER(LEN=IRANK),         INTENT(OUT), OPTIONAL :: axis ! um_ak_20091016
    ! op_pj_20100121+
    INTEGER, DIMENSION(IRANK),    INTENT(OUT), OPTIONAL :: nbounds
    ! op_pj_20100121-

    ! LOCAL
    TYPE(t_channel),              POINTER :: channel  => NULL()
    TYPE(t_channel_object),       POINTER :: object => NULL()

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object)
    IF (status /= 0) RETURN

    IF (PRESENT(lrestart_read)) &
         lrestart_read = object%int%lrestart_read(SND_INS)

    IF (PRESENT(reprid)) &
         reprid = object%repr%id

    IF (PRESENT(axis)) &             ! um_ak_20091016
         axis = object%repr%axis     ! um_ak_20091016

    IF (PRESENT(nbounds)) &                         ! op_pj_20100121
         nbounds(:) = object%repr%bounds%nbounds(:) ! op_pj_20100121

  END SUBROUTINE get_channel_object_info
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE new_channel_object_reference(status, &
       cname1, oname1, cname2, oname2, lcopyatt)

    USE messy_main_channel_attributes, ONLY: copy_attribute_list

    IMPLICIT NONE

    INTRINSIC :: PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)          :: status
    CHARACTER(LEN=*), INTENT(IN)           :: cname1   ! ORIGINAL CHANNEL NAME
    CHARACTER(LEN=*), INTENT(IN)           :: oname1   ! ORIGINAL OBJECT NAME
    CHARACTER(LEN=*), INTENT(IN)           :: cname2   ! NEW CHANNEL NAME
    CHARACTER(LEN=*), INTENT(IN)           :: oname2   ! NEW OBJECT NAME
    LOGICAL,          INTENT(IN), OPTIONAL :: lcopyatt ! COPY ALL ATTRIBUTES?

    ! LOCAL
    TYPE(t_channel_object),       POINTER  :: object1 => NULL()
    TYPE(t_channel_object),       POINTER  :: object2 => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER  :: mem
    LOGICAL                                :: zlcopyatt

    ! CHECKS
    IF (LFIXATE) THEN
       status = 3000
       RETURN
    END IF

    ! INIT
    IF (PRESENT(lcopyatt)) THEN
       zlcopyatt = lcopyatt
    ELSE
       zlcopyatt = .FALSE.   ! DEFAULT
    END IF

    CALL loc_channel_object(status, TRIM(cname1), TRIM(oname1), object1)
    IF (status /= 0) RETURN

    ! SHARE PRIMARY MEMORY
    ! SAME REPRESENTATAION
    mem => object1%data(:,:,:,:)
    CALL new_channel_object(status, TRIM(cname2), TRIM(oname2), mem=mem &
         ,reprid = object1%repr%id)
    IF (status /= 0) RETURN

    ! PRE-SET I/O SETTINGS
    CALL loc_channel_object(status, TRIM(cname2), TRIM(oname2), object2)
    IF (status /= 0) RETURN

    object2%lrestreq    = .FALSE.       ! SHARED MEMORY !!!
    object2%io          = object1%io
    object2%io%lrestart = .FALSE.       ! SHARED MEMORY !!!

    object2%int%lref    = .TRUE.        ! REFERENCE

    ! COPY ALL ATTRIBUTES
    IF (zlcopyatt) THEN
       CALL copy_attribute_list(status, object1%att, object2%att)
       IF (status /= 0) RETURN
    END IF

    ! ADD SPECIAL ATTRIBUTE
    CALL new_channel_object_attribute(status, TRIM(cname2), TRIM(oname2) &
         , 'REFERENCE_TO', c=TRIM(cname1)//': '//TRIM(oname1) &
         , loverwrite=.TRUE.  )

  END SUBROUTINE new_channel_object_reference
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_channel_object_restreq(status, cname, oname)

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                      INTENT(OUT)           :: status
    ! CHANNEL NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: cname
    ! OBJECT NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: oname

    ! LOCAL
    TYPE(t_channel_object),       POINTER :: object => NULL()

    ! CHECKS
    IF (LFIXATE) THEN
       status = 3000
       RETURN
    END IF

    CALL loc_channel_object(status, TRIM(cname), TRIM(oname), object)
    IF (status /= 0) RETURN
    
    object%lrestreq = .TRUE.

  END SUBROUTINE set_channel_object_restreq
  ! -------------------------------------------------------------------

  ! op_pj_20100827+
  ! -------------------------------------------------------------------
  SUBROUTINE get_channel_object_dimvar(status, cname, oname, dva, units) 

    USE messy_main_constants_mem,      ONLY: STRLEN_ULONG
    USE messy_main_tools,              ONLY: PTR_1D_ARRAY
    USE messy_main_channel_dimvar,     ONLY: t_dimvar, get_dimvar
    USE messy_main_channel_attributes, ONLY: return_attribute

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM, ADJUSTL

    ! I/O
    INTEGER,                      INTENT(OUT)           :: status
    ! CHANNEL NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: cname
    ! OBJECT NAME
    CHARACTER(LEN=*),             INTENT(IN)            :: oname
    ! POINTER TO DIMENSION VARIABLE DATA ARRAYS ! INTENT(OUT)
    TYPE (PTR_1D_ARRAY), DIMENSION(:), POINTER          :: dva
    ! UNITS  ! INTENT(OUT)
    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER, OPTIONAL :: units

    ! LOCAL
    INTEGER :: reprid
    TYPE(t_representation), POINTER  :: repr
    INTEGER :: i
    TYPE(t_dimvar),      POINTER :: dimvar => NULL()

    ! CLEAN FIRST
    IF (ASSOCIATED(dva)) DEALLOCATE(dva)
    NULLIFY(dva)

    IF (PRESENT(units)) THEN
       IF (ASSOCIATED(units)) DEALLOCATE(units)
       NULLIFY(units)
    END IF

    CALL get_channel_object_info(status, cname, oname, reprid=reprid)
    IF (status /= 0) RETURN

    CALL get_representation(status, reprid, repr)
    IF (status /= 0) RETURN

    ! ALLOCATE
    ALLOCATE(dva(repr%rank))
    DO i=1, repr%rank 
       NULLIFY(dva(i)%ptr)
    END DO

    IF (PRESENT(units)) THEN
       ALLOCATE(units(repr%rank))
       DO i=1, repr%rank 
          units(i) = ''
       END DO
    END IF

    dimension_loop: DO i=1, repr%rank

       CALL get_dimvar(status, repr%dim(i)%ptr%var &
            , TRIM(ADJUSTL(repr%dim(i)%ptr%name)), dimvar)
       ! 953: DIMENSION VARIABLE DOES NOT EXIST (OK HERE!)
       IF ((status /= 953) .AND. (status /= 0)) RETURN

       dva(i)%ptr => dimvar%val(:)

       IF (PRESENT(units)) THEN
          CALL return_attribute(status, dimvar%att, 'units', c=units(i))
          ! 805: ATTRIBUTE DOES NOT EXIST (OK HERE!)
          IF ((status /= 805) .AND. (status /= 0)) RETURN
       ENDIF

    END DO dimension_loop

    status = 0

  END SUBROUTINE get_channel_object_dimvar
  ! -------------------------------------------------------------------
  ! op_pj_20100827-

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_object_by_ptr(status, object)

    USE messy_main_channel_attributes, ONLY: write_attribute

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_channel_object), POINTER     :: object

    IF (.NOT. ASSOCIATED(object)) THEN
       status = 3107 ! CHANNEL OBJECT POINTER NOT ASSOCIATED
       RETURN
    END IF

    WRITE(*,'(1x,a24,1x,11(L1,1x),i4,1x,L1,1x,2(i8,1x))') &
         object%name    &
         , object%int%lrst                        & ! ANY RESTART ?
         , object%lrestreq                        & ! REQUIRED IN RESTART
         , object%int%lign                        & ! IGNORE lrestreq ?
         , object%io%lrestart                     & ! USER DEFINED
         , object%io%lout(:)                                 & ! WHICH OUTPUT
         , object%repr%id                                    &
         , object%memory%lalloc                              &
         , object%memory%usage, object%memory%usage_2nd

    CALL write_attribute(status, object%att)

    status = 0

  END SUBROUTINE write_channel_object_by_ptr
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_channel_object_by_name(status, cname, oname)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: cname
    CHARACTER(LEN=*), INTENT(IN)  :: oname

    ! LOCAL
    TYPE(t_channel_object), POINTER :: object

    CALL loc_channel_object(status, cname, oname, object)
    IF (status /= 0) RETURN

    CALL write_channel_object_by_ptr(status, object)

  END SUBROUTINE write_channel_object_by_name
  ! -------------------------------------------------------------------

  ! **********************************************************************
  ! ALL CHANNELS
  ! **********************************************************************

  ! -------------------------------------------------------------------
  SUBROUTINE fixate_channels(status)

    USE messy_main_tools,         ONLY: match_wild
    USE messy_main_channel_error, ONLY: channel_error_str
    USE messy_main_constants_mem, ONLY: STRLEN_VLONG

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, ANY, SIZE, TRIM, INDEX

    ! I/O
    INTEGER, INTENT(OUT) :: status

    ! LOCAL
    TYPE(t_channel_list),        POINTER :: ls
    TYPE(t_channel_object_list), POINTER :: le
    TYPE(t_channel),             POINTER :: channel
    TYPE(t_channel_object),      POINTER :: object
    INTEGER                              :: zstat
    INTEGER                              :: i,i1,i2,i3,i4
    INTEGER                              :: js, je
    INTEGER                              :: jsnd
    INTEGER                              :: m
    CHARACTER(LEN=STRLEN_CHANNEL)        :: s1 = ''
    CHARACTER(LEN=STRLEN_CHANNEL)        :: s2 = ''
    CHARACTER(LEN=STRLEN_OBJECT)         :: e1 = ''
    CHARACTER(LEN=STRLEN_OBJECT)         :: e2 = ''
    CHARACTER(LEN=STRLEN_VLONG)          :: errstr
    ! um_ak_20100305+
    ! search pointer for wildcards in source objects
    LOGICAL                              :: l_wild_exist
    TYPE(t_channel_list),        POINTER :: s_ls
    TYPE(t_channel_object_list), POINTER :: s_le
    TYPE(t_channel),             POINTER :: s_channel
    TYPE(t_channel_object),      POINTER :: s_object
    ! um_ak_20100305-

    ! CHECKS
    IF (LFIXATE) THEN
       status = 3000
       RETURN
    END IF

    ! CREATE NEW CHANNELS
    DO js=1, NMAXADDCHANNEL
       IF (TRIM(ADJUSTL(ADD_CHANNEL(js)%cname)) == '') CYCLE
       CALL new_channel(status, TRIM(ADJUSTL(ADD_CHANNEL(js)%cname)))
       IF (status /= 0) RETURN
    END DO

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT

       channel => ls%this

       ! -------------------------------------------
       ! CREATE NEW CHANNEL OBJECT REFERENCES
       ! -------------------------------------------
       DO je=1, NMAXADDREF
          s1 = TRIM(ADJUSTL(ADD_REF(je)%cname1))
          IF (TRIM(s1) == '') CYCLE
          e1 = TRIM(ADJUSTL(ADD_REF(je)%oname1))
          IF (TRIM(e1) == '') CYCLE

          s2 = TRIM(ADJUSTL(ADD_REF(je)%cname2))
          IF (TRIM(s2) == '') CYCLE
          IF ( TRIM(ADJUSTL(ADD_REF(je)%oname2)) == '' ) THEN
             e2 = TRIM(e1)                          ! USE SAME OBJECT NAME
          ELSE
             e2 = TRIM(ADJUSTL(ADD_REF(je)%oname2)) ! NEW NAME
          END IF          
          !
          ! WILDCARD-MATCH BUT NO SELF REFERENCE WITH IDENTICAL OBJECT NAME
          IF ( match_wild( TRIM(s2), TRIM(channel%name) )  .AND.    &
               ( .NOT. ( (TRIM(channel%name) == TRIM(s1)) .AND. &
               (TRIM(e2) == TRIM(e1)) ) ) ) THEN
             ! um_ak_20100305+
             ! enable wild matches in the source object name
             l_wild_exist = (INDEX(e1, '?') > 0) .OR. (INDEX(e1, '*') > 0)
             if_wc_ext: IF (l_wild_exist) THEN
                s_ls => GCHANNELLIST
                channel_loop2: DO
                   IF (.NOT. ASSOCIATED(s_ls)) EXIT

                   s_channel => s_ls%this
                   IF (TRIM(s_channel%name) == TRIM(s1)) THEN
                      s_le => s_channel%list
                      object_loop2: DO
                         IF (.NOT. ASSOCIATED(s_le)) EXIT

                         s_object => s_le%this
                         
                         IF (match_wild(TRIM(e1),TRIM(s_object%name))) THEN
                            !Note: when the source object name includes a
                            !      wildcard, the object names in the 
                            !      destination channel will be the same as no
                            !      specific name can be attributed for wildcard
                            !      strings. => e2 == s_object%name
                            CALL new_channel_object_reference(status       &
                                 , TRIM(s1), TRIM(s_object%name)           &
                                 , TRIM(channel%name), TRIM(s_object%name) &
                                 , lcopyatt = .TRUE.)
                            SELECT CASE(status)
                            CASE(0)
                               ! OK: CONTINUE
                            CASE(3003, 3103)
                               ! 3003: SOURCE CHANNEL DOES NOT EXIST
                               ! 3103: SOURCE CHANNEL OBJECT DOES NOT EXIST
                               WRITE(*,*) '**** WARNING CALLING', &
                                    'new_channel_object_reference: ' &
                                    , TRIM(s1),' / ', TRIM(s_object%name),' / '&
                                    , TRIM(channel%name), ' / '                &
                                    , TRIM(s_object%name)
                               errstr = channel_error_str(status)
                               WRITE(*,*) TRIM(errstr)
                            CASE DEFAULT
                               ! SEVERE ERROR
                               WRITE(*,*) '**** ERROR CALLING' &
                                    , 'new_channel_object_reference: ' &
                                    , TRIM(s1),' / ', TRIM(s_object%name),' / '&
                                    , TRIM(channel%name)                       &
                                    , ' / ', TRIM(s_object%name)
                               RETURN
                            END SELECT
                         END IF ! match wild source object
                         s_le => s_le%next
                     END DO object_loop2
                   END IF
                   s_ls => s_ls%next
                END DO channel_loop2
             ELSE
             ! um_ak_20100305-
                CALL new_channel_object_reference(status, &
                     TRIM(s1), TRIM(e1), TRIM(channel%name), TRIM(e2), &
                     lcopyatt = .TRUE.)
                SELECT CASE(status)
                CASE(0)
                   ! OK: CONTINUE
                CASE(3003, 3103)
                   ! 3003: SOURCE CHANNEL DOES NOT EXIST
                   ! 3103: SOURCE CHANNEL OBJECT DOES NOT EXIST
                   WRITE(*,*) &
                        '*** WARNING CALLING new_channel_object_reference: ' &
                        , TRIM(s1), ' / ', TRIM(e1), ' / '     &
                        , TRIM(channel%name), ' / ', TRIM(e2)
                   errstr = channel_error_str(status)
                   WRITE(*,*) TRIM(errstr)
                CASE DEFAULT
                   ! SEVERE ERROR
                   WRITE(*,*) &
                        '*** ERROR CALLING new_channel_object_reference: ' &
                        , TRIM(s1), ' / ', TRIM(e1), ' / '     &
                        , TRIM(channel%name), ' / ', TRIM(e2)
                   RETURN
                END SELECT
             ENDIF if_wc_ext ! wildcard exists um_ak_20100305
          ENDIF
       END DO

       ! -------------------------------------------
       ! (1) SET EVERYTHING TO DEFAULT VALUE
       ! -------------------------------------------
       ! op_bk_20130905+
       ! Workaround for gfortran < 4.8.1, wrong handling of namelist input
       ! with more than one "sub-struct"
#ifdef __GFORTRAN__
       channel%io%ftype    = OUT_DEFAULT%ftype
       channel%io%ntpf     = OUT_DEFAULT%ntpf
       channel%default%io%lrestart  = OUT_DEFAULT%lrestart
       channel%default%io%lignore   = OUT_DEFAULT%lignore
       channel%default%io%lout      = OUT_DEFAULT%lout
       channel%default%io%range     = OUT_DEFAULT%range
#else
       channel%io          = OUT_DEFAULT%cio
       channel%default%io  = OUT_DEFAULT%oio
#endif
       ! op_bk_20130905-
       ! -------------------------------------------
       ! (2) SET TO SPECIFIC, IF AVAIALABLE
       ! -------------------------------------------
       DO js=1, NMAXCHANNELS
          IF (TRIM(OUT_CHANNEL(js)%cname) == '') CYCLE
          IF (match_wild(TRIM(OUT_CHANNEL(js)%cname), TRIM(channel%name))) THEN
             ! op_bk_20130905+
             ! Workaround for gfortran < 4.8.1, wrong handling of namelist input
             ! with more than one "sub-struct"
#ifdef __GFORTRAN__
             channel%io%ftype    = OUT_CHANNEL(js)%ftype
             channel%io%ntpf     = OUT_CHANNEL(js)%ntpf
             channel%default%io%lrestart  = OUT_CHANNEL(js)%lrestart
             channel%default%io%lignore   = OUT_CHANNEL(js)%lignore
             channel%default%io%lout      = OUT_CHANNEL(js)%lout
             channel%default%io%range     = OUT_CHANNEL(js)%range
#else
             channel%io          = OUT_CHANNEL(js)%cio
             channel%default%io  = OUT_CHANNEL(js)%oio
#endif
       ! op_bk_20130905-
          END IF
       END DO

       ! -------------------------------------------
       ! (3) INITIALIZE INTERNAL SETUP
       ! -------------------------------------------
       channel%int%lign = .TRUE.  ! NEEDED FOR INITIALIZATION LOGICS BELOW

       le => channel%list
       object_loop: DO
          IF (.NOT. ASSOCIATED(le)) EXIT

          object => le%this

          ! -------------------------------------------
          ! (1) SET EVERYTHING TO DEFAULT VALUE
          ! -------------------------------------------
          object%io = channel%default%io

          ! -------------------------------------------
          ! (2) SET TO SPECIFIC, IF AVAIALABLE
          ! -------------------------------------------
          DO je=1, NMAXOBJECTS
             IF ( (TRIM(OUT_OBJECT(je)%cname) == '') .OR.      &
                  (TRIM(OUT_OBJECT(je)%oname) == '')) CYCLE
             IF ( match_wild(TRIM(OUT_OBJECT(je)%cname), TRIM(channel%name) ) &
                  .AND. &
                  match_wild(TRIM(OUT_OBJECT(je)%oname), TRIM(object%name)) ) &
                  THEN
                  object%io = OUT_OBJECT(je)%io
             END IF
          END DO

          ! -------------------------------------------
          ! (3) INTERNAL SETUP
          ! -------------------------------------------
          ! CHANNEL OBJECT EXPORT
          !
          ! - OUTPUT
          !   (USER DEFINED)
          object%int%lexp(:,IOMODE_OUT) = object%io%lout(:)
          ! - ANY OUTPUT ?
          object%int%lout = ANY(object%int%lexp(:, IOMODE_OUT))
          !
          ! - RESTART 
          !   (REQUESTED STATISTICS ALWAYS REQUIRES RESTART OUTPUT ...
          object%int%lexp(:,IOMODE_RST) = object%io%lout(:)
          !   ... INSTANTANEOUS VALUE ONLY ON REQUEST
          object%int%lexp(SND_INS,IOMODE_RST) = &
               object%io%lrestart .OR. & ! FORCED BY USER OR ...
               object%lrestreq           ! ... REQUIRED (CODE !)
          !   ... x2 (STD) REQUIRES x1 (AVE)
          object%int%lexp(SND_AVE,IOMODE_RST) = &
               object%int%lexp(SND_AVE,IOMODE_RST) .OR. &
               object%int%lexp(SND_STD,IOMODE_RST)
          !   ... csm (CAV) REQUIRES cnt (CNT)
          object%int%lexp(SND_CNT,IOMODE_RST) = &
               object%int%lexp(SND_CAV,IOMODE_RST) .OR. &
               object%int%lexp(SND_CNT,IOMODE_RST)
          ! - ANY RESTART ?
          object%int%lrst = ANY(object%int%lexp(:, IOMODE_RST))
          ! - IGNORE lrestreq
          object%int%lign = object%io%lignore

          ! - 2ndary MEMORY MANAGEMENT (1)
          DO jsnd=2, SND_MAXLEN
             IF (object%int%lexp(jsnd, IOMODE_RST)) THEN
                object%int%n2nd = object%int%n2nd + 1
                object%int%i2nd(jsnd) = object%int%n2nd
             END IF
          END DO
          !
          ! - 2ndary MEMORY MANAGEMENT (2)
          i1 = SIZE(object%data, 1)
          i2 = SIZE(object%data, 2)
          i3 = SIZE(object%data, 3)
          i4 = SIZE(object%data, 4)
          ALLOCATE(object%sdat(object%int%n2nd), STAT=zstat)
          IF (zstat /= 0) THEN
             status = 1000
             RETURN
          END IF
          DO i=1, object%int%n2nd
             ALLOCATE(object%sdat(i)%ptr(i1,i2,i3,i4), STAT=zstat)
             IF (zstat /= 0) THEN
                status = 1000
                RETURN
             END IF
             object%sdat(i)%ptr(:,:,:,:) = 0.0_DP
             object%memory%usage_2nd = &
                  object%memory%usage_2nd + &
                  INT( SIZE(object%sdat(i)%ptr), I8 )
          END DO

          ! CHANNEL
          channel%memory%usage_2nd = &
               channel%memory%usage_2nd + object%memory%usage_2nd
          ! - ANY OUTPUT ?
          channel%int%lout = channel%int%lout .OR. object%int%lout
          ! - ANY OBJECT REQUIRED IN RESTART ?
          channel%int%lrestreq = channel%int%lrestreq .OR. object%lrestreq
          ! - ANY RESTART ?
          channel%int%lrst = channel%int%lrst .OR. object%int%lrst
          ! - IGNORE ALL lrestreq ? (ONLY IF ALL ARE TRUE !!!)
          channel%int%lign = channel%int%lign .AND. object%int%lign
          ! -------------------------------------------

          ! -------------------------------------------
          ! OUTPUT I/O
          ! -------------------------------------------
          DO m = 1, IOMODE_MAX
             SELECT CASE(channel%io%ftype(m))
             CASE (FTYPE_NETCDF, FTYPE_PNETCDF)
                !
                ! netCDF
                !
                ALLOCATE(object%int%netcdf(m)%svarid(object%int%n2nd))
                object%int%netcdf(m)%svarid(:) = NC_ID_UNDEF
                !
                !
                !CASE(...) 
                ! +++ ADD OTHER OUTPUT FORMATS HERE
             CASE DEFAULT
                !
             END SELECT
          END DO

          le => le%next
       END DO object_loop

       ! -------------------------------------------
       ! OUTPUT TIMER
       ! -------------------------------------------
!!$       ! INITIALIZE COUNTER TO TRIGGER NEW FILE AT BEGINNING
!!$       channel%int%ntpfcnt = channel%io%ntpf
       ! TRIGGER NEW FILE AT BEGINNING
       channel%int%lnew_file = .TRUE.

       ls => ls%next
    END DO channel_loop

    ! SET FIXATED FLAG
    LFIXATE = .TRUE.

    status = 0

  END SUBROUTINE fixate_channels
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
! mz_pj_20080905+
!!$  SUBROUTINE trigger_channel_output(status, lnow, lforce_new)
  SUBROUTINE trigger_channel_output(status, lnow, ltnf, lforce_new)
! mz_pj_20080905-

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE

    ! I/O
    INTEGER,               INTENT(OUT) :: status
    LOGICAL, DIMENSION(:), INTENT(IN)  :: lnow       ! output now
    ! mz_pj_20080905+
    LOGICAL, DIMENSION(:), INTENT(IN)  :: ltnf       ! new file now
    ! mz_pj_20080905-
    LOGICAL,               INTENT(IN)  :: lforce_new ! new files for all now

    ! LOCAL
    TYPE(t_channel_list),         POINTER :: ls
    TYPE(t_channel),              POINTER :: channel
    LOGICAL                               :: lfull

    ! CHECKS
    IF (SIZE(lnow) /= NCHANNEL) THEN
       status = 3005
       RETURN
    END IF

    ! mz_pj_20080905+
    IF (SIZE(ltnf) /= NCHANNEL) THEN
       status = 3005
       RETURN
    END IF
    ! mz_pj_20080905-

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT

       channel => ls%this

       ! -------------------------------------------
!       channel%int%lout_now = lnow(channel%id) .AND. channel%int%lout
       channel%int%lout_now = channel%int%lout .AND. &
            (lnow(channel%id) .OR. channel%int%lforce_out) .AND. &
            (.NOT. channel%int%lsuppr_out)

       ! CHECK IF NEW FILE IS REQUIRED BECAUSE OLD IS FULL
       lfull = .FALSE.
       ! INCREMENT COUNTER
       IF (channel%int%lout_now) THEN
          channel%int%ntpfcnt   = channel%int%ntpfcnt + 1
       END IF
       IF (channel%io%ntpf >0) THEN          ! mz_pj_20080905
          ! COUNTER                          ! mz_pj_20080905
          lfull                 = (channel%int%ntpfcnt > channel%io%ntpf)
       ELSE                                  ! mz_pj_20080905
          ! EVENT TRIGGERED                  ! mz_pj_20080905
          lfull = ltnf(channel%id)           ! mz_pj_20080905
       END IF                                ! mz_pj_20080905
          ! SET TRIGGER FOR NEW FILE
          !  a) KEEP STATUS (RESET IN channel_finish_io (ONLY AFTER OUTPUT))
          !  b) OLD IS FULL
          !  c) FORCED (e.g., after restart) AND ANY OUTPUT REQUESTED
          !  d) FORCED BY USER (e.g. EVENT) AND ANY OUTPUT REQUESTED
       channel%int%lnew_file = channel%int%lnew_file .OR. &  ! a)
            lfull .OR.                                    &  ! b)
            (lforce_new .AND. channel%int%lout) .OR.      &  ! c)
            (channel%int%lforce_newfile .AND. channel%int%lout) ! d)
       ! RESET COUNTER IF NEW FILE IS TRIGGERED
       IF (channel%int%lnew_file) channel%int%ntpfcnt = 1

       ! RESET (only for one time step)
       channel%int%lforce_out = .FALSE.
       channel%int%lsuppr_out = .FALSE.
       channel%int%lforce_newfile = .FALSE.
       ! -------------------------------------------

       ls => ls%next
    END DO channel_loop

    status = 0

  END SUBROUTINE trigger_channel_output
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE update_channels(status, flag, dtime)

    USE messy_main_constants_mem, ONLY: FLAGGED_BAD

    IMPLICIT NONE

    INTRINSIC :: ABS, ASSOCIATED, SQRT, MIN, MAX

    ! I/O
    INTEGER,  INTENT(OUT) :: status
    INTEGER,  INTENT(IN)  :: flag
    REAL(DP), INTENT(IN)  :: dtime  ! time step length

    ! LOCAL
    TYPE(t_channel_list),        POINTER :: ls
    TYPE(t_channel_object_list), POINTER :: le
    TYPE(t_channel),             POINTER :: channel
    TYPE(t_channel_object),      POINTER :: object
    INTEGER                              :: n, m
    REAL(DP)                             :: quot
    LOGICAL                              :: zlout_now

! op_pj_20100616+
!!$    ! FIRST ACCUMULATION STEP
!!$    LOGICAL, SAVE                        :: lfirst = .TRUE.
!!$    LOGICAL                              :: zlfirst

!!$    ! INIT
!!$    SELECT CASE(flag)
!!$    CASE(1)
!!$       ! ACCUMULATE FIELDS --------------------------------
!!$       ! ... INITIALIZE MINIMUM AND MAXIMUM (VERY FIRST STEP
!!$       ! ... AND FIRST TIME STEP AFTER OUTPUT)
!!$       !     OR ACCUMULATE (LATER)
!!$       zlfirst = lfirst   ! SAVE (TO BE APPLIED BELOW ...)
!!$       lfirst = .FALSE.   ! RESET IN VERY FIRST CALL
!!$       !
!!$    CASE(2)
!!$       ! PREPARE FOR OUTPUT --------------------------------
!!$    CASE(3)
!!$       ! RESET AFTER OUTPUT --------------------------------
!!$    END SELECT
! op_pj_20100616-

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT

       channel => ls%this

       ! -------------------------------------------
       SELECT CASE(flag)
       CASE(1)
          ! ACCUMULATE FIELDS --------------------------------
          !
          ! TIME INTERVAL
          channel%int%tslo = channel%int%tslo + dtime
          !
! op_pj_20100616+
          ! ACCUMULATE FIELDS --------------------------------
          ! ... INITIALIZE MINIMUM AND MAXIMUM (VERY FIRST STEP
          !     AND FIRST TIME STEP AFTER OUTPUT)
          ! ... OR ACCUMULATE (LATER)
          ! SAVE (TO BE APPLIED BELOW ...)
          channel%int%l2ndreinit2 = channel%int%l2ndreinit1
          ! RESET IN VERY FIRST CALL
          channel%int%l2ndreinit1 = .FALSE.
! op_pj_20100616-
          !
       CASE(2)
          ! PREPARE FOR OUTPUT --------------------------------
          !
          IF (channel%int%lout_now) THEN
             IF (channel%int%tslo > 0.0_DP) THEN
                quot = 1.0_DP / channel%int%tslo 
             ELSE
                quot = 1.0_DP
             END IF
          END IF
          !
       CASE(3)
          ! RESET AFTER OUTPUT --------------------------------
          !
          zlout_now = channel%int%lout_now  ! -> SAVE FOR OBJECTS BELOW
          IF (channel%int%lout_now) THEN
             ! TIME INTERVAL
             channel%int%tslo = 0.0_DP
             ! OUTPUT FLAG
             channel%int%lout_now = .FALSE.
          END IF
          !
       CASE DEFAULT
          ! ERROR
          !
          STATUS = 3106
          RETURN
          !
       END SELECT
       ! -------------------------------------------

       le => channel%list
       object_loop: DO
          IF (.NOT. ASSOCIATED(le)) EXIT

          object => le%this

          ! -------------------------------------------
          SELECT CASE(flag)
          CASE(1)
             ! ACCUMULATE FIELDS -------------------------------- 
             !
             ! AVERAGE: SUM X
             n = object%int%i2nd(SND_AVE)
             IF (n /= 0) THEN
                object%sdat(n)%ptr(:,:,:,:) = object%sdat(n)%ptr(:,:,:,:) &
                     + object%data(:,:,:,:) * dtime
             END IF
             !
             ! STANDARD DEVIATION: SUM X^2
             n = object%int%i2nd(SND_STD)
             IF (n /= 0) THEN
                object%sdat(n)%ptr(:,:,:,:) = object%sdat(n)%ptr(:,:,:,:) &
                  + object%data(:,:,:,:)**2 * dtime
             END IF
             ! MINIMUM
             n = object%int%i2nd(SND_MIN)
             IF (n /= 0) THEN
! op_pj_20100616+
!!$                IF (zlfirst .AND. &
                IF (channel%int%l2ndreinit2 .AND. &
! op_pj_20100616-
                     (.NOT. object%int%lrestart_read(SND_MIN))) THEN
                   object%sdat(n)%ptr(:,:,:,:) = object%data(:,:,:,:)
                ELSE
                   object%sdat(n)%ptr(:,:,:,:) = &
                        MIN(object%sdat(n)%ptr(:,:,:,:), object%data(:,:,:,:))
                END IF
             END IF
             ! MAXIMUM
             n = object%int%i2nd(SND_MAX)
             IF (n /= 0) THEN
! op_pj_20100616+
!!$                IF (zlfirst .AND. &
                IF (channel%int%l2ndreinit2 .AND. &
! op_pj_20100616-
                     (.NOT. object%int%lrestart_read(SND_MAX))) THEN
                   object%sdat(n)%ptr(:,:,:,:) = object%data(:,:,:,:)
                ELSE
                   object%sdat(n)%ptr(:,:,:,:) = &
                        MAX(object%sdat(n)%ptr(:,:,:,:), object%data(:,:,:,:))
                END IF
             END IF
             !
             ! CONDITION COUNTER
             n = object%int%i2nd(SND_CNT)
             IF (n /= 0) THEN
                WHERE( (object%data(:,:,:,:) >= object%io%range(1)) .AND. &
                     (object%data(:,:,:,:) <= object%io%range(2)) )
                       object%sdat(n)%ptr = object%sdat(n)%ptr + 1.0_DP
                END WHERE
             END IF
             !
             ! CONDITIONAL AVERAGE: SUM
             n = object%int%i2nd(SND_CAV)
             IF (n /= 0) THEN
                 WHERE( (object%data(:,:,:,:) >= object%io%range(1)) .AND. &
                        (object%data(:,:,:,:) <= object%io%range(2)) )
                        object%sdat(n)%ptr = object%sdat(n)%ptr + object%data
                 END WHERE
             END IF
             !
          CASE (2)
             ! PREPARE FOR OUTPUT --------------------------------
             IF (channel%int%lout_now) THEN
                ! AVERAGE
                n = object%int%i2nd(SND_AVE)
                m = n
                IF (n > 0) THEN
                   object%sdat(n)%ptr(:,:,:,:) = &
                        object%sdat(n)%ptr(:,:,:,:) * quot
                END IF
                ! STANDARD DEVIATION
                n = object%int%i2nd(SND_STD)
                IF (n > 0) THEN
                   object%sdat(n)%ptr(:,:,:,:) = &
                        SQRT( ABS (object%sdat(n)%ptr(:,:,:,:)*quot   &
                                 - object%sdat(m)%ptr(:,:,:,:)**2 ) )
                END IF
                ! CONDITIONAL AVERAGE
                n = object%int%i2nd(SND_CAV)
                m = object%int%i2nd(SND_CNT)
                IF (n > 0) THEN
                   WHERE (object%sdat(m)%ptr(:,:,:,:) > 0.0_DP)
                      object%sdat(n)%ptr = &
                           object%sdat(n)%ptr / &
                           object%sdat(m)%ptr
                   ELSEWHERE
                      !object%sdat(n)%ptr = REAL(-HUGE(0.0_SP),DP)
                      object%sdat(n)%ptr = FLAGGED_BAD
                   END WHERE
                END IF
             END IF
             !
          CASE (3)
             ! RESET AFTER OUTPUT --------------------------------
             !
             IF (zlout_now) THEN
                ! 2ndary FIELDS
                IF (object%int%n2nd > 0) THEN
                   DO n=1, object%int%n2nd
                      object%sdat(n)%ptr(:,:,:,:) = 0.0_DP
                   END DO
! op_pj_20100616+
!!$                   lfirst = .TRUE. ! um_ak_20100616
                   channel%int%l2ndreinit1 = .TRUE.
! op_pj_20100616-
                END IF
             END IF
             !
          CASE DEFAULT
             ! ERROR
             !
             STATUS = 3106
             RETURN
             !
          END SELECT
          ! -------------------------------------------

          le => le%next
       END DO object_loop

       ls => ls%next
    END DO channel_loop

    status = 0

  END SUBROUTINE update_channels
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE clean_channels(status)

    USE messy_main_channel_attributes, ONLY: clean_attribute_list

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE

    ! I/O
    INTEGER,  INTENT(OUT) :: status

    ! LOCAL
    TYPE(t_channel_list),         POINTER :: lsai, lsae
    TYPE(t_channel_object_list),  POINTER :: leai, leae
    TYPE(t_channel),              POINTER :: channel
    TYPE(t_channel_object),       POINTER :: object
    INTEGER                               :: zstat
    INTEGER                               :: n, m

    status = 0
    zstat = 0

    lsai => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(lsai)) EXIT

       channel => lsai%this

       lsae => lsai
       lsai => lsai%next

       leai => channel%list

       object_loop: DO
          IF (.NOT. ASSOCIATED(leai)) EXIT

          object => leai%this

          leae => leai
          leai => leai%next

          ! -------------------------------------------
          ! ATTRIBUTES
          CALL clean_attribute_list(status, object%att)
          IF (status /= 0) RETURN
          !
          ! MEMORY: DATA
          !
          IF (object%memory%lalloc) THEN
             object%memory%usage = object%memory%usage &
                  - INT( SIZE(object%data), I8 )
             channel%memory%usage = channel%memory%usage &
                  - INT( SIZE(object%data), I8 )
             DEALLOCATE(object%data, STAT=zstat)
             IF (zstat /= 0) THEN
                status = 1001 ! MEMORY DEALLOCATION FAILED
                RETURN
             END IF
             IF (object%memory%usage /= 0_I8) THEN
                status = 3108 ! CHANNEL OBJECT PRIMARY MEMORY ERROR
                RETURN
             END IF
          END IF
          !
          ! MEMORY: SDAT
          DO n=1, object%int%n2nd
             object%memory%usage_2nd = object%memory%usage_2nd &
                  - INT(SIZE(object%sdat(n)%ptr), I8)
             channel%memory%usage_2nd = channel%memory%usage_2nd &
                  - INT(SIZE(object%sdat(n)%ptr), I8)
             DEALLOCATE(object%sdat(n)%ptr, STAT=zstat)
             IF (zstat /= 0) THEN
                status = 1001 ! MEMORY DEALLOCATION FAILED
                RETURN
             END IF
          END DO

          IF (object%int%n2nd > 0) THEN
             DEALLOCATE(object%sdat, STAT=zstat)
             NULLIFY(object%sdat)
          END IF
          IF (zstat /= 0) THEN
             status = 1001 ! MEMORY DEALLOCATION FAILED
             RETURN
          END IF

          IF (object%memory%usage_2nd /= 0_I8) THEN
             status = 3109 ! CHANNEL OBJECT SECONDARY MEMORY ERROR
             RETURN
          END IF
          !
          ! LIST OBJECT
          DEALLOCATE(leae)
          NULLIFY(leae)
          ! -------------------------------------------

          ! -------------------------------------------
          ! IO
          ! -------------------------------------------
          DO m=1, IOMODE_MAX
             ! netCDF
             IF (ASSOCIATED(object%int%netcdf(m)%svarid)) &
                  DEALLOCATE(object%int%netcdf(m)%svarid)
             !
             ! +++ ADD OTHER OUTPUT FORMATS HERE
          END DO

       END DO object_loop

       ! -------------------------------------------
       ! MEMORY: DATA
       IF (channel%memory%usage /= 0_I8) THEN
          status = 3108 ! CHANNEL SECONDARY MEMORY ERROR
          RETURN
       END IF
       ! MEMORY: SDAT
       IF (channel%memory%usage_2nd /= 0_I8) THEN
          status = 3109 ! CHANNEL SECONDARY MEMORY ERROR
          RETURN
       END IF
       !
       ! ATTRIBUTES
       CALL clean_attribute_list(status, channel%att)
       IF (status /= 0) RETURN
       !
       ! LIST OBJECT
       DEALLOCATE(lsae)
       NULLIFY(lsae)
       !
       ! COUNT
       NCHANNEL = NCHANNEL - 1
       ! -------------------------------------------

    END DO channel_loop

    NULLIFY(GCHANNELLIST)

    IF (NCHANNEL /= 0) THEN
       status = 1004
       RETURN
    END IF

    status = 0

  END SUBROUTINE clean_channels
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_read_nml_ctrl(status, iou)

    ! MODULE ROUTINE (CORE)
    !
    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2004

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! mz_pj_20080118: OUT_PREC added
    NAMELIST /CTRL/ EXP_NAME, L_FLUSH_IOBUFFER, I_VERBOSE_LEVEL &
         , ADD_CHANNEL, ADD_REF, OUT_DEFAULT, OUT_PREC, OUT_CHANNEL, OUT_OBJECT

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status
    INTEGER                     :: js           ! mz_pj_20081029

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE
! mz_pj_20081029+
! FOR SOME STRANGE REASONS THIS IS REQUIRED FOR LF8.1
    DO js=1, NMAXADDCHANNEL
       ADD_CHANNEL(js)%cname = ''
    END DO
    L_FLUSH_IOBUFFER = .TRUE. ! op_pj_20100819 this also
    I_VERBOSE_LEVEL = 1 ! op_pj_20110803
! mz_pj_20081029-

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_channel_read_nml_ctrl
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE loc_channel_by_name(status, list, name, channel)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, LEN_TRIM, TRIM

    ! I/O
    INTEGER,              INTENT(OUT) :: status
    TYPE(t_channel_list), POINTER     :: list
    CHARACTER(LEN=*),     INTENT(IN)  :: name
    TYPE(t_channel),      POINTER     :: channel

    ! LOCAL
    TYPE(t_channel_list), POINTER :: ai  => NULL()
    TYPE(t_channel_list), POINTER :: ae  => NULL()
    LOGICAL                       :: lex

    ! INIT
    lex = .FALSE.
    NULLIFY(channel)

    ! CHECKS
    IF (TRIM(name) == '') THEN
       status = 3010 ! CHANNEL NAME IS EMPTY
       RETURN
    END IF
    IF (LEN_TRIM(ADJUSTL(name)) > STRLEN_CHANNEL) THEN
       status = 3001  ! CHANNEL NAME TOO LONG
       RETURN
    END IF
    !
    IF (.NOT. ASSOCIATED(list)) THEN
       status = 3003  ! CHANNEL (NAME) DOES NOT EXIST
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
       channel => ai%this
    ELSE
       status = 3003  ! CHANNEL DOES NOT EXIST
       RETURN
    END IF

    status = 0


  END SUBROUTINE loc_channel_by_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE loc_channel_by_id(status, list, id, channel)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,              INTENT(OUT)          :: status
    TYPE(t_channel_list), POINTER              :: list
    INTEGER,              INTENT(IN)           :: id
    TYPE(t_channel),      POINTER              :: channel

    ! LOCAL
    TYPE(t_channel_list), POINTER :: ai  => NULL()
    TYPE(t_channel_list), POINTER :: ae  => NULL()
    LOGICAL                       :: lex

    ! INIT
    lex = .FALSE.
    NULLIFY(channel)

    ! CHECKS
    IF (id <= 0) THEN
       status = 3007  ! INVALID CHANNEL ID
       RETURN
    END IF
    !
    IF (.NOT. ASSOCIATED(list)) THEN
       status = 3004  ! CHANNEL (ID) DOES NOT EXIST
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
       channel => ai%this
    ELSE
       status = 3003  ! CHANNEL DOES NOT EXIST
       RETURN
    END IF

    status = 0

  END SUBROUTINE loc_channel_by_id
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE loc_channel_object(status, cname, oname, object)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, LEN_TRIM, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)   :: status
    CHARACTER(LEN=*), INTENT(IN)    :: cname
    CHARACTER(LEN=*), INTENT(IN)    :: oname
    TYPE(t_channel_object), POINTER :: object

    ! LOCAL
    TYPE(t_channel),             POINTER :: channel => NULL()
    TYPE(t_channel_object_list), POINTER :: ai     => NULL()
    TYPE(t_channel_object_list), POINTER :: ae     => NULL()
    LOGICAL                              :: lex

    ! INIT
    lex = .FALSE.
    NULLIFY(object)

    CALL loc_channel(status, GCHANNELLIST, TRIM(cname), channel)
    IF (status /= 0) RETURN

    ! CHECKS
    IF (TRIM(oname) == '') THEN
       status = 3110 ! CHANNEL OBJECT NAME IS EMPTY
    END IF
    IF (LEN_TRIM(ADJUSTL(oname)) > STRLEN_OBJECT) THEN
       status = 3101   ! CHANNEL OBJECT NAME TOO LONG
       RETURN
    END IF

    ! CHECK, IF NAME EXISTS
    ai => channel%list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       IF (TRIM(ADJUSTL(oname)) == TRIM(ai%this%name)) THEN
          lex = .TRUE.
          EXIT
       END IF
       ae => ai
       ai => ai%next
    END DO

    IF (lex) THEN
       object => ai%this
    ELSE
       status = 3103     ! CHANNEL OBJECT DOES NOT EXIST
       RETURN
    END IF

    status = 0

  END SUBROUTINE loc_channel_object
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_main_channel
! **********************************************************************
