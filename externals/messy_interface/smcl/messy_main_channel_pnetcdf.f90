! **********************************************************************
MODULE messy_main_channel_pnetcdf
! **********************************************************************

  USE messy_main_channel
#ifdef PNETCDF
  USE messy_main_constants_mem,      ONLY: STRLEN_MEDIUM ! op_pj_20091124
#endif

  IMPLICIT NONE
  PRIVATE
  SAVE

#ifdef PNETCDF
!  INCLUDE 'pnetcdf.inc'
#include <pnetcdf.inc>
 INCLUDE 'mpif.h'
#endif

  INTEGER :: p_pe
  INTEGER :: p_io
  INTEGER :: p_all_comm
  LOGICAL :: p_parallel_io

  PUBLIC :: ch_pnetcdf_init_pio

#ifdef PNETCDF
  INTEGER :: PNCOUT_PREC = NF_FLOAT ! mz_pj_20080118 default
  INTEGER :: MY_MPI_INFO            ! mz_pj_20081124

  PUBLIC :: ch_pnetcdf_init_rst
  PUBLIC :: ch_pnetcdf_init_io
  PUBLIC :: ch_pnetcdf_write_header
  PUBLIC :: ch_pnetcdf_write_time
  PUBLIC :: ch_pnetcdf_write_data
  PUBLIC :: ch_pnetcdf_finish_io
  !
  PUBLIC :: ch_pnetcdf_read_data
  !
! op_pj_20091124+
  PUBLIC :: ch_pnetcdf_read_nml_ctrl
! op_pj_20091124-
  !
  !PRIVATE :: pnetcdf_write_attribute_list
  !PRIVATE :: pnetcdf_define_dimvar_list
  !PRIVATE :: pnetcdf_write_dimvar_list
  !PRIVATE :: pnetcdf_check_attributes
  !PRIVATE :: pnf

! op_pj_20091124+
  ! FOR PERFORMANCE TUNING
  TYPE T_MPI_IO_HINTS
     CHARACTER(LEN=STRLEN_MEDIUM) :: hint  = ''
     CHARACTER(LEN=STRLEN_MEDIUM) :: value = ''
  END TYPE T_MPI_IO_HINTS
  PUBLIC :: T_MPI_IO_HINTS
  !
  INTEGER, PARAMETER, PUBLIC :: NMAX_MPI_IO_HINTS = 10
  !
  TYPE(T_MPI_IO_HINTS), DIMENSION(NMAX_MPI_IO_HINTS), PUBLIC :: MPI_IO_HINT
! op_pj_20091124-
#endif

CONTAINS

  ! ----------------------------------------------------------------

  SUBROUTINE ch_pnetcdf_init_pio(status, ex_p_pe, ex_p_io, ex_p_all_comm)

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: ex_p_pe
    INTEGER, INTENT(IN)  :: ex_p_io
    INTEGER, INTENT(IN)  :: ex_p_all_comm

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch_pnetcdf_init_pio'

    ! mz_pj_20081124+
#ifdef PNETCDF
    INTEGER :: ierror, ierror2
    CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: mpierrstr
    INTEGER :: errlen
    INTEGER :: i
#endif
    ! mz_pj_20081124-

    p_pe = ex_p_pe
    p_io = ex_p_io
    p_parallel_io = (p_pe == p_io)
    p_all_comm = ex_p_all_comm

    ! mz_pj_20081124+
#ifdef PNETCDF

    MY_MPI_INFO = MPI_INFO_NULL

    DO i=1, NMAX_MPI_IO_HINTS
       IF (TRIM(MPI_IO_HINT(i)%hint) /= "") THEN
          IF (p_parallel_io) &
               WRITE(*,*) substr,': setting MPI_IO_HINTS:'
          CALL MPI_INFO_CREATE(MY_MPI_INFO, ierror)
          IF (ierror /= MPI_SUCCESS) THEN
             CALL MPI_ERROR_STRING(ierror, mpierrstr, errlen, ierror2)
             WRITE(*,*) '*** MPI ERROR: ',mpierrstr(1:errlen)
             status = -1
             RETURN
          END IF
          EXIT
       END IF
    END DO
    
    DO i=1, NMAX_MPI_IO_HINTS
       IF (TRIM(MPI_IO_HINT(i)%hint) /= "") THEN
          IF (p_parallel_io) &
               WRITE(*,*) '  ',TRIM(MPI_IO_HINT(i)%hint),' = ' &
               ,TRIM(MPI_IO_HINT(i)%value)
          CALL MPI_INFO_SET(MY_MPI_INFO &
               , TRIM(MPI_IO_HINT(i)%hint), TRIM(MPI_IO_HINT(i)%value), ierror)
          IF (ierror /= MPI_SUCCESS) THEN
             CALL MPI_ERROR_STRING(ierror, mpierrstr, errlen, ierror2)
             WRITE(*,*) '*** MPI ERROR: ',mpierrstr(1:errlen)
             status = -1
             RETURN
          END IF
       END IF
    END DO

#endif

    status = 0

  END SUBROUTINE ch_pnetcdf_init_pio
  ! ----------------------------------------------------------------

#ifdef PNETCDF

  ! -------------------------------------------------------------------
  SUBROUTINE ch_pnetcdf_init_rst(status, fname, att)

    USE messy_main_channel_attributes, ONLY: t_attribute_list

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    CHARACTER(LEN=*),       INTENT(IN)  :: fname
    TYPE(t_attribute_list), POINTER     :: att 

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch_pnetcdf_init_rst'
    INTEGER                     :: ncid

    CALL pnf(NFMPI_OPEN( p_all_comm, TRIM(fname) &
         , NF_NOWRITE, MY_MPI_INFO, ncid), status, substr)
    IF (status /= 0) RETURN

    CALL pnetcdf_check_attributes(status, ncid, NF_GLOBAL, att) 
    IF (status /= 0) RETURN

    CALL pnf(NFMPI_CLOSE(ncid), status, substr)
    IF (status /= 0) RETURN

  END SUBROUTINE ch_pnetcdf_init_rst
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_pnetcdf_init_io(status, IOMODE, channel, AMODE, att)

    USE messy_main_channel_attributes, ONLY: t_attribute_list

    IMPLICIT NONE

    INTRINSIC :: TRIM, PRESENT, IOR

    ! I/O
    INTEGER,           INTENT(OUT) :: status
    INTEGER,           INTENT(IN)  :: IOMODE
    TYPE(t_channel),   POINTER     :: channel    ! INTENT(INOUT)    
    INTEGER,           INTENT(IN)  :: AMODE
    TYPE(t_attribute_list), POINTER, OPTIONAL :: att ! INTENT(IN)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'ch_pnetcdf_init_io'
    INTEGER                        :: ncid
    INTEGER                        :: ifmode
    LOGICAL                        :: lexist
    INTEGER                        :: nc_cmode
    LOGICAL, SAVE                  :: lfirst = .TRUE. ! mz_pj_20080118

    ! mz_pj_20080118+
    IF (lfirst) THEN
       SELECT CASE(OUT_PREC(FTYPE_PNETCDF))
       CASE(1)
          PNCOUT_PREC = NF_FLOAT
       CASE(2)
          PNCOUT_PREC = NF_DOUBLE
       CASE DEFAULT
          status = 4003 ! UNKNOWN PRECISION FLAG FOR (P)NETCDF
          RETURN
       END SELECT
       lfirst = .FALSE.
    END IF
    ! mz_pj_20080118-

    SELECT CASE(AMODE)
    CASE(AMODE_READ)
       !
       SELECT CASE(IOMODE)
       CASE(IOMODE_OUT)
          !
          status = 3205 ! NO INPUT OF OUTPUT FILES
          RETURN
          !
       CASE(IOMODE_RST)
          !
          IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
               WRITE(*,*) substr,': OPENING  FILE: ' &
               , TRIM(channel%int%fname(IOMODE))
          !
          INQUIRE(file = TRIM(channel%int%fname(IOMODE)), exist = lexist)
          IF (.NOT.lexist) THEN
             IF (channel%int%lrestreq) THEN
                IF (channel%int%lign) THEN
                   ! IGNORE ALL MISSING RESTART FIELDS
                   IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
                        WRITE(*,*) substr, &
                        ':       WARNING: REQUIRED RESTART ', &
                        'FILE NOT PRESENT ... ALL OBJECTS TO BE IGNORED'
                   ! RESET: PREVENT FROM CLOSING AND READING
                   channel%int%fname(IOMODE) = ''
                   status = 0
                ELSE
                   status = 3206 ! RESTART FILE REQUIRED BUT NOT PRESENT
                END IF
             ELSE
                IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
                     WRITE(*,*) substr,':       WARNING: RESTART ', &
                     'FILE NOT PRESENT ... NOT REQUIRED'
                ! RESET: PREVENT FROM CLOSING AND READING
                channel%int%fname(IOMODE) = ''
                status = 0
             END IF
             RETURN
          END IF
          !
          CALL pnf(NFMPI_OPEN( p_all_comm &
               , TRIM(channel%int%fname(IOMODE)) &
               , NF_NOWRITE, MY_MPI_INFO, ncid) &
               , status, substr)
          IF (status /= 0) RETURN
          channel%int%netcdf(IOMODE)%fileID = ncid
          !
          ! OPTIONAL SPECIAL ATTRIBUTES
          IF (PRESENT(att)) THEN
             CALL pnetcdf_check_attributes(status, ncid, NF_GLOBAL, att)
             IF (status /= 0) RETURN
          END IF
          CALL pnetcdf_check_attributes(status, ncid, NF_GLOBAL, GATT)
          IF (status /= 0) RETURN
          CALL pnetcdf_check_attributes(status, ncid, NF_GLOBAL, channel%att)
          IF (status /= 0) RETURN
          !
          IF (I_VERBOSE_LEVEL >= 1) WRITE(*,*)
          !
       END SELECT
       !
    CASE(AMODE_WRITE)
       !
       ! ---------------------------------------------------------
       ! CLOSE 'OLD' FILE FIRST, IF STILL OPEN
       ! ---------------------------------------------------------
       IF (channel%int%netcdf(IOMODE)%fileID /= NC_ID_UNDEF) THEN
          CALL pnf(NFMPI_CLOSE(channel%int%netcdf(IOMODE)%fileID) &
               , status, substr)
          IF (status /= 0) RETURN
          channel%int%netcdf(IOMODE)%fileID = NC_ID_UNDEF
       END IF
       !
       ! ---------------------------------------------------------
       ! OPEN 'NEW' FILE
       ! ---------------------------------------------------------
       IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 2)) &
            WRITE(*,*) substr,': CREATING FILE: ' &
            ,TRIM(channel%int%fname(IOMODE))
!qqq+
!       IF (IOMODE == IOMODE_RST) THEN
          nc_cmode = IOR(NF_CLOBBER,NF_64BIT_OFFSET)
!       ELSE
!          nc_cmode = NF_CLOBBER
!       END IF
!qqq-
       CALL pnf( &
            NFMPI_CREATE( p_all_comm &
            , TRIM(channel%int%fname(IOMODE)) &
            , nc_cmode &
            , MY_MPI_INFO &
            , ncid)    &
            , status, substr)
       IF (status /= 0) RETURN
       channel%int%netcdf(IOMODE)%fileID = ncid
       ! SWITCH OFF 'FILL MODE'
       CALL pnf(NFMPI_SET_FILL(ncid, NF_NOFILL, ifmode), status, substr)
       IF (status /= 0) RETURN
       !
    END SELECT

    status = 0

  END SUBROUTINE ch_pnetcdf_init_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_pnetcdf_write_header(status, IOMODE, channel, dim_time, att)

    USE messy_main_channel_dimensions, ONLY: NDIM, t_dimension
    USE messy_main_channel_attributes, ONLY: t_attribute_list
    USE messy_main_channel_repr,       ONLY: IRANK

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM, PRESENT

    ! I/O
    INTEGER,                INTENT(OUT)       :: status
    INTEGER,                INTENT(IN)        :: IOMODE
    TYPE(t_channel),        POINTER           :: channel    ! INTENT(INOUT)
    TYPE(t_dimension),      POINTER           :: dim_time   ! INTENT(IN)
    TYPE(t_attribute_list), POINTER, OPTIONAL :: att        ! INTENT(IN)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'ch_pnetcdf_write_header'
    INTEGER                              :: ncid
    INTEGER                              :: ifmode
    TYPE(t_channel_object_list), POINTER :: le
    TYPE(t_channel_object),      POINTER :: object
    INTEGER                              :: i
    INTEGER                              :: i2nd  ! INDEX IN SECONDARY DATA
    INTEGER                              :: jsnd  ! OUTPUT DATA TYPE
    INTEGER, DIMENSION(:), ALLOCATABLE   :: dimid
    LOGICAL, DIMENSION(:), ALLOCATABLE   :: ldim_out
    INTEGER                              :: varid
    INTEGER                              :: nrank
    INTEGER, DIMENSION(IRANK+1)          :: dimvec
    LOGICAL                              :: lskip
    INTEGER                              :: it
    INTEGER                              :: XNF90_PREC
    INTEGER                              :: PREC

    ! INIT
    ncid = channel%int%netcdf(IOMODE)%fileID

    ! mz_pj_20090115+
    CALL new_attribute(status, channel%att, 'channel_netcdf_lib' &
         , c = NFMPI_inq_libvers(), loverwrite=.TRUE.)
    IF (status /= 0) RETURN
    ! mz_pj_20090115-

    ! ---------------------------------------------------------
    ! WRITE GLOBAL ATTRIBUTES
    ! ---------------------------------------------------------
    ! - COMMON TO ALL CHANNELS
    CALL pnetcdf_write_attribute_list(status, ncid, NF_GLOBAL, GATT)
    IF (status /= 0) RETURN
    ! - CHANNEL SPECIFIC
    CALL pnetcdf_write_attribute_list(status, ncid, NF_GLOBAL, &
         channel%att)
    IF (status /= 0) RETURN
    ! - OPTIONAL SPECIAL ATRIBUTES
    IF (PRESENT(att)) THEN
       CALL pnetcdf_write_attribute_list(status, ncid, NF_GLOBAL, att)
       IF (status /= 0) RETURN
    END IF

    ! ---------------------------------------------------------
    ! DEFINE UNLIMITED DIMENSION WITH ATTRIBUTES
    ! ---------------------------------------------------------
    IF (IOMODE == IOMODE_OUT) THEN
       CALL pnf( &
            NFMPI_DEF_DIM( ncid, TRIM(dim_time%name) &
            , NFMPI_UNLIMITED                           &
            , channel%int%netcdf(IOMODE)%dimid_time) &
            , status, substr)
       IF (status /= 0) RETURN      
       CALL pnetcdf_define_dimvar_list(status       &
            , ncid                                  & ! netCDF file ID
            , channel%int%netcdf(IOMODE)%dimid_time & ! dimension ID
            , dim_time%var                          & ! list of variables
            , NF_DOUBLE                             & ! precision
            )
       IF (status /= 0) RETURN
    END IF
    
    ! ---------------------------------------------------------
    ! DEFINE DIMENSIONS AND DIMENSION VARIABLES WITH ATTRIBUTES
    ! ---------------------------------------------------------
    ! - INIT
    ALLOCATE(dimid(NDIM))
    dimid(:) = NC_ID_UNDEF
    ALLOCATE(ldim_out(NDIM))
    ldim_out(:) = .TRUE.

    !mz_pj_20060628+
    ! SET PRECISION
    IF (IOMODE == IOMODE_OUT) THEN
       ! mz_pj_20080118+
!!$       XNF90_PREC = NF_FLOAT
       XNF90_PREC = PNCOUT_PREC
       ! mz_pj_20080118-
    ELSE
       XNF90_PREC = NF_DOUBLE
    END IF
    !mz_pj_20060628-

    ! - LOOP OVER ALL DIMENSIONS IN ALL REPRESENTATIONS OF ALL OBJECTS
    le => channel%list
    object_loop1: DO
       IF (.NOT. ASSOCIATED(le)) EXIT
       object => le%this

       ! ANY OUTPUT FOR THIS OBJECT ???
       SELECT CASE(IOMODE)
       CASE(IOMODE_OUT)
          lskip = .NOT. object%int%lout
       CASE(IOMODE_RST)
          lskip = .NOT. object%int%lrst
       END SELECT
       !
       IF (lskip) THEN
          le => le%next          
          CYCLE
       END IF

       ! -------------------------------------------------------
       DO i=1, IRANK
          IF (ASSOCIATED(object%repr%dim(i)%ptr)) THEN
             IF (dimid(object%repr%dim(i)%ptr%id) == NC_ID_UNDEF) THEN
                !
                CALL pnf( &
                     NFMPI_DEF_DIM( ncid, TRIM(object%repr%dim(i)%ptr%name) &
                     , INT(object%repr%dim(i)%ptr%len, MPI_OFFSET_KIND)     &
                     , dimid(object%repr%dim(i)%ptr%id) )                  &
                     , status, substr)
                IF (status /= 0) RETURN
                !
                CALL pnetcdf_define_dimvar_list(status    &
                     , ncid                              & ! netCDF file ID
                     , dimid(object%repr%dim(i)%ptr%id)  & ! dimension ID
                     , object%repr%dim(i)%ptr%var        & ! list of variables
                     , XNF90_PREC                        & ! precision
                     )
                IF (status /= 0) RETURN
             END IF
             ! SAVE netCDF dimension IDs
             object%int%netcdf(IOMODE)%dimid(i) = &
                  dimid(object%repr%dim(i)%ptr%id)
          END IF
       END DO
       ! -------------------------------------------------------

       le => le%next
    END DO object_loop1

    ! ---------------------------------------------------------
    ! DEFINE VARIABLES (INCLUDING SECONDARY DATA)
    ! ---------------------------------------------------------
    ! - LOOP OVER ALL DIMENSIONS IN ALL REPRESENTATIONS OF ALL OBJECTS
    le => channel%list
    object_loop2: DO
       IF (.NOT. ASSOCIATED(le)) EXIT
       object => le%this

       ! ANY OUTPUT FOR THIS OBJECT ???
       SELECT CASE(IOMODE)
       CASE(IOMODE_OUT)
          lskip = .NOT. object%int%lout
       CASE(IOMODE_RST)
          lskip = .NOT. object%int%lrst
       END SELECT
       !
       IF (lskip) THEN
          le => le%next          
          CYCLE
       END IF

       ! -------------------------------------------------------
       ! SET DIMENSION VECTOR
       ! PERMUTATION ACCORDING TO REPRESENTATION
       nrank           = object%repr%rank
       DO i=1, nrank
          dimvec(i) = object%int%netcdf(IOMODE)%dimid( &
               object%repr%order_mem2out(i) )
       END DO
       IF (IOMODE == IOMODE_OUT) THEN
          it = 1
          dimvec(nrank+1) = channel%int%netcdf(IOMODE)%dimid_time
          ! mz_pj_20080118+
!!$          XNF90_PREC = NF_FLOAT
          XNF90_PREC = PNCOUT_PREC
          ! mz_pj_20080118-
       ELSE
          it = 0
          XNF90_PREC = NF_DOUBLE
       END IF

       ! DEFINE PRIMARY AND SECONDARY DATA
       DO jsnd=1, SND_MAXLEN

          IF (.NOT. object%int%lexp(jsnd, IOMODE)) CYCLE ! NO OUTPUT

          CALL pnf( &
               NFMPI_DEF_VAR(ncid           & ! netCDF file ID
               , TRIM(object%name)//&
               &TRIM(SND_TEXT(jsnd,IOMODE)) & ! variable name
               , XNF90_PREC                 & ! precision
               , nrank+it                   &
               , dimvec(1:nrank+it)         & ! netCDF dim-IDs
               , varid)                     & ! netCDF var-ID
               , status, substr )
          IF (status /= 0) RETURN
          IF (jsnd == SND_INS) then
             object%int%netcdf(IOMODE)%varid = varid        ! primary data
          ELSE
             i2nd = object%int%i2nd(jsnd)
             object%int%netcdf(IOMODE)%svarid(i2nd) = varid ! 2ndary data
          END IF
          ! WRITE VARIABLE ATTRIBUTES
          CALL pnetcdf_write_attribute_list(status, ncid, varid, &
               object%att)
          ! WRITE SPECIAL ATTRIBUTES FOR SECONDARY DATA
          IF ((jsnd == SND_CNT) .OR. (jsnd == SND_CAV)) THEN
             CALL pnf(NFMPI_PUT_ATT_INT(ncid, varid, 'range' &
                  , NF_INT, INT(2, MPI_OFFSET_KIND) &
                  , object%io%range) &
                  , status, substr)
             IF (status /= 0) RETURN
          END IF
       END DO
       ! -------------------------------------------------------

       le => le%next
    END DO object_loop2

    ! ---------------------------------------------------------
    ! SWITCH FROM DEFINE- TO DATA-OUTPUT-MODE
    ! ---------------------------------------------------------
    CALL pnf(NFMPI_ENDDEF(ncid), status, substr)
    IF (status /= 0) RETURN

    CALL pnf(NFMPI_BEGIN_INDEP_DATA(ncid), status, substr)
    IF (status /= 0) RETURN

    ! ---------------------------------------------------------
    ! WRITE DIMENSION-VARIABLES
    ! ---------------------------------------------------------
    ! - LOOP OVER ALL DIMENSIONS IN ALL REPRESENTATIONS OF ALL OBJECTS
    le => channel%list
    object_loop3: DO
       IF (.NOT. ASSOCIATED(le)) EXIT
       object => le%this

       ! ANY OUTPUT FOR THIS OBJECT ???
       SELECT CASE(IOMODE)
       CASE(IOMODE_OUT)
          lskip = .NOT. object%int%lout
          PREC  = SP
       CASE(IOMODE_RST)
          lskip = .NOT. object%int%lrst
          PREC  = DP
       END SELECT
       !
       IF (lskip) THEN
          le => le%next          
          CYCLE
       END IF

       ! -------------------------------------------------------
       DO i=1, IRANK
          IF (ASSOCIATED(object%repr%dim(i)%ptr)) THEN
             IF (ldim_out(object%repr%dim(i)%ptr%id)) THEN

                IF (p_parallel_io) THEN

                CALL pnetcdf_write_dimvar_list(status    &
                     , ncid                              & ! netCDF file ID
                     , object%repr%dim(i)%ptr%var        & ! list of variables
                     , PREC                              & ! precision
                     , 1                                 & ! start
                     )
                IF (status /= 0) RETURN

                END IF

                ldim_out(object%repr%dim(i)%ptr%id) = .FALSE. ! not twice
             END IF
          END IF
       END DO
       ! -------------------------------------------------------

       le => le%next
    END DO object_loop3

    CALL pnf(NFMPI_END_INDEP_DATA(ncid), status, substr)

    ! ---------------------------------------------------------
    ! FLUSH netCDF BUFFER
    ! ---------------------------------------------------------
    IF (L_FLUSH_IOBUFFER) THEN
       CALL pnf(NFMPI_SYNC(ncid), status, substr)
       IF (status /= 0) RETURN
    END IF

    ! ---------------------------------------------------------
    ! CLEAN UP
    ! ---------------------------------------------------------
    DEALLOCATE(dimid)
    DEALLOCATE(ldim_out)

  END SUBROUTINE ch_pnetcdf_write_header
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_pnetcdf_write_time(status, channel, dim_time)

    USE messy_main_channel_dimensions, ONLY: t_dimension

    IMPLICIT NONE
    
    ! I/O
    INTEGER,           INTENT(OUT) :: status
    TYPE(t_channel),   POINTER     :: channel ! INTENT(INOUT)
    TYPE(t_dimension), POINTER     :: dim_time

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'ch_pnetcdf_write_time'
    INTEGER :: ncid
    INTEGER :: outstep

    ncid    = channel%int%netcdf(IOMODE_OUT)%fileID
    IF (ncid == NC_ID_UNDEF) THEN
       status = 4001 ! UNDEFINED FILE-ID
       RETURN
    END IF

    outstep = channel%int%ntpfcnt

    CALL pnf(NFMPI_BEGIN_INDEP_DATA(ncid), status, substr)
    IF (status /= 0) RETURN

    IF (p_parallel_io) THEN
       CALL pnetcdf_write_dimvar_list(status  &
            , ncid                           & ! netCDF file ID
            , dim_time%var                   & ! list of variables
            , DP                             & ! precision
            , outstep                        & ! start
            )
    END IF

    CALL pnf(NFMPI_END_INDEP_DATA(ncid), status, substr)
    
  END SUBROUTINE ch_pnetcdf_write_time
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_pnetcdf_write_data(status, IOMODE, channel, object &
       , ptr, jsnd, i2nd)

    USE messy_main_channel_repr,  ONLY: IRANK, PIOTYPE_SGL, PIOTYPE_IND &
                                      , PIOTYPE_COL

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE, TRIM, REAL, NULL

    ! I/O
    INTEGER,                INTENT(OUT)   :: status
    INTEGER,                INTENT(IN)    :: IOMODE
    TYPE(t_channel),        POINTER       :: channel ! INTENT(INOUT)
    TYPE(t_channel_object), POINTER       :: object  ! INTENT(INOUT)
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: ptr
    INTEGER,                INTENT(IN)    :: jsnd    ! OUTPUT DATA TYPE
    INTEGER,                INTENT(IN)    :: i2nd    ! index of 2ndary DATA

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch_pnetcdf_write_data'
    INTEGER                               :: ncid
    INTEGER                               :: varid
    INTEGER                               :: nrank
    INTEGER(MPI_OFFSET_KIND), DIMENSION(IRANK+1) :: start, cnt
    REAL(DP),                     POINTER :: p0 => NULL()
    REAL(DP), DIMENSION(:),       POINTER :: p1 => NULL()
    REAL(DP), DIMENSION(:,:),     POINTER :: p2 => NULL()
    REAL(DP), DIMENSION(:,:,:),   POINTER :: p3 => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: p4 => NULL()
    INTEGER                               :: it
    INTEGER                               :: iseg, nseg
    INTEGER                               :: l1,l2,l3,l4
    INTEGER                               :: r1,r2,r3,r4

    IF (I_VERBOSE_LEVEL >= 2) THEN ! op_pj_20110803
       SELECT CASE(IOMODE)
       CASE(IOMODE_OUT)
          IF (p_parallel_io) &
               WRITE(*,*) '... parallel netCDF-output: '         &
               , TRIM(channel%int%fname(IOMODE))   & 
               , ' (',channel%int%ntpfcnt,') <- '  &
               , TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE))
       CASE(IOMODE_RST)
          IF (p_parallel_io) &
               WRITE(*,*) '... parallel netCDF-output: '         &
               , TRIM(channel%int%fname(IOMODE))   &
               , ' <- '                            &
               , TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE))
       END SELECT
    END IF ! op_pj_20110803

    ncid  = channel%int%netcdf(IOMODE)%fileID
    IF (ncid == NC_ID_UNDEF) THEN
       status = 4001 ! UNDEFINED FILE-ID
       RETURN
    END IF

    nrank = object%repr%rank

    SELECT CASE(jsnd)
    CASE(SND_UNDEF)
       status = 3203 ! SECONDARY DATA TYPE UNKNOWN
       RETURN
    CASE(SND_INS)
       varid = object%int%netcdf(IOMODE)%varid
    CASE DEFAULT ! SND_AVE, SND_STP, SND_MIN, SND_MAX, ...
       varid = object%int%netcdf(IOMODE)%svarid(i2nd)
    END SELECT

    start(:)  = 1
    cnt(:)    = 1

    SELECT CASE(IOMODE)
    CASE(IOMODE_OUT)
       it = 1
       start(nrank + 1) = channel%int%ntpfcnt   ! OUTPUT STEP IN netCDF FILE
       cnt(nrank + 1)   = 1                     ! 1 TIME STEP
    CASE(IOMODE_RST)
       it = 0                                   ! NO TIME AXIS IN RESTART FILE
    END SELECT

    ! LOOP OVER SEGMENTS
    nseg  = object%repr%pdecomp%nseg

    ! -----------------------------------------------------
    ! DATA MODE (independent or collective)
    SELECT CASE(object%repr%pdecomp%piotype)
    ! -----------------------------------------------------

    ! -----------------------------------------------------
    CASE(PIOTYPE_IND, PIOTYPE_SGL)
    ! -----------------------------------------------------

    CALL pnf(NFMPI_BEGIN_INDEP_DATA(ncid), status, substr)
    IF (status /= 0) RETURN

    SELECT CASE(IOMODE)
    CASE(IOMODE_OUT)

       DO iseg=1, nseg

          start(1:nrank) = object%repr%pdecomp%start(iseg,1:nrank)
          cnt(1:nrank)   = object%repr%pdecomp%cnt(iseg,1:nrank)
          
          l1 = object%repr%pdecomp%ml(iseg,1)
          l2 = object%repr%pdecomp%ml(iseg,2)
          l3 = object%repr%pdecomp%ml(iseg,3)
          l4 = object%repr%pdecomp%ml(iseg,4)
          
          r1 = object%repr%pdecomp%mu(iseg,1)
          r2 = object%repr%pdecomp%mu(iseg,2)
          r3 = object%repr%pdecomp%mu(iseg,3)
          r4 = object%repr%pdecomp%mu(iseg,4)

          SELECT CASE(nrank)
          CASE(0)
             p0 => ptr(1,1,1,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_REAL(ncid     & 
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , (/ REAL(p0, SP) /) ) &
                  , status, substr )
          CASE(1)
             p1 => ptr(l1:r1,1,1,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_REAL(ncid     & 
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , REAL(p1, SP) )      &
                  , status, substr )
             
          CASE(2)
             p2 => ptr(l1:r1,l2:r2,1,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_REAL(ncid     &
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , REAL(p2, SP) )      &
                  , status, substr )
             
          CASE(3)
             p3 => ptr(l1:r1,l2:r2,l3:r3,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_REAL(ncid     &
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , REAL(p3, SP) )      &
                  , status, substr )
             
          CASE(4)
             p4 => ptr(l1:r1,l2:r2,l3:r3,l4:r4)
             CALL pnf( &
                  NFMPI_PUT_VARA_REAL(ncid     &
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , REAL(p4, SP) )      &
                  , status, substr )
             
          END SELECT
          
       END DO
       
    CASE(IOMODE_RST)

       DO iseg=1, nseg

          start(1:nrank) = object%repr%pdecomp%start(iseg,1:nrank)
          cnt(1:nrank)   = object%repr%pdecomp%cnt(iseg,1:nrank)
          
          l1 = object%repr%pdecomp%ml(iseg,1)
          l2 = object%repr%pdecomp%ml(iseg,2)
          l3 = object%repr%pdecomp%ml(iseg,3)
          l4 = object%repr%pdecomp%ml(iseg,4)
          
          r1 = object%repr%pdecomp%mu(iseg,1)
          r2 = object%repr%pdecomp%mu(iseg,2)
          r3 = object%repr%pdecomp%mu(iseg,3)
          r4 = object%repr%pdecomp%mu(iseg,4)
          
          SELECT CASE(nrank)
          CASE(0)
             p0 => ptr(1,1,1,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_DOUBLE(ncid     & 
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , (/ p0 /) )          &
                  , status, substr )
          CASE(1)
             p1 => ptr(l1:r1,1,1,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_DOUBLE(ncid     & 
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , p1 )                &
                  , status, substr )
             
          CASE(2)
             p2 => ptr(l1:r1,l2:r2,1,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_DOUBLE(ncid     &
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , p2 )                &
                  , status, substr )
             
          CASE(3)
             p3 => ptr(l1:r1,l2:r2,l3:r3,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_DOUBLE(ncid     &
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , p3 )                &
                  , status, substr )
             
          CASE(4)
             p4 => ptr(l1:r1,l2:r2,l3:r3,l4:r4)
             CALL pnf( &
                  NFMPI_PUT_VARA_DOUBLE(ncid     &
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , p4 )                &
                  , status, substr )
             
          END SELECT
          
       END DO

    END SELECT

    CALL pnf(NFMPI_END_INDEP_DATA(ncid), status, substr)
    IF (status /= 0) RETURN

    ! -----------------------------------------------------
    CASE(PIOTYPE_COL)
    ! -----------------------------------------------------

    SELECT CASE(IOMODE)
    CASE(IOMODE_OUT)

       DO iseg=1, nseg

          start(1:nrank) = object%repr%pdecomp%start(iseg,1:nrank)
          cnt(1:nrank)   = object%repr%pdecomp%cnt(iseg,1:nrank)
          
          l1 = object%repr%pdecomp%ml(iseg,1)
          l2 = object%repr%pdecomp%ml(iseg,2)
          l3 = object%repr%pdecomp%ml(iseg,3)
          l4 = object%repr%pdecomp%ml(iseg,4)
          
          r1 = object%repr%pdecomp%mu(iseg,1)
          r2 = object%repr%pdecomp%mu(iseg,2)
          r3 = object%repr%pdecomp%mu(iseg,3)
          r4 = object%repr%pdecomp%mu(iseg,4)

          SELECT CASE(nrank)
          CASE(0)
             p0 => ptr(1,1,1,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_REAL_ALL(ncid     & 
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , (/ REAL(p0, SP) /) ) &
                  , status, substr )
          CASE(1)
             p1 => ptr(l1:r1,1,1,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_REAL_ALL(ncid     & 
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , REAL(p1, SP) )      &
                  , status, substr )
             
          CASE(2)
             p2 => ptr(l1:r1,l2:r2,1,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_REAL_ALL(ncid     &
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , REAL(p2, SP) )      &
                  , status, substr )
             
          CASE(3)
             p3 => ptr(l1:r1,l2:r2,l3:r3,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_REAL_ALL(ncid     &
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , REAL(p3, SP) )      &
                  , status, substr )
             
          CASE(4)
             p4 => ptr(l1:r1,l2:r2,l3:r3,l4:r4)
             CALL pnf( &
                  NFMPI_PUT_VARA_REAL_ALL(ncid     &
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , REAL(p4, SP) )      &
                  , status, substr )
             
          END SELECT
          
       END DO
       
    CASE(IOMODE_RST)

       DO iseg=1, nseg

          start(1:nrank) = object%repr%pdecomp%start(iseg,1:nrank)
          cnt(1:nrank)   = object%repr%pdecomp%cnt(iseg,1:nrank)
          
          l1 = object%repr%pdecomp%ml(iseg,1)
          l2 = object%repr%pdecomp%ml(iseg,2)
          l3 = object%repr%pdecomp%ml(iseg,3)
          l4 = object%repr%pdecomp%ml(iseg,4)
          
          r1 = object%repr%pdecomp%mu(iseg,1)
          r2 = object%repr%pdecomp%mu(iseg,2)
          r3 = object%repr%pdecomp%mu(iseg,3)
          r4 = object%repr%pdecomp%mu(iseg,4)
          
          SELECT CASE(nrank)
          CASE(0)
             p0 => ptr(1,1,1,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_DOUBLE_ALL(ncid     & 
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , (/ p0 /) )          &
                  , status, substr )
          CASE(1)
             p1 => ptr(l1:r1,1,1,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_DOUBLE_ALL(ncid     & 
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , p1 )                &
                  , status, substr )
             
          CASE(2)
             p2 => ptr(l1:r1,l2:r2,1,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_DOUBLE_ALL(ncid     &
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , p2 )                &
                  , status, substr )
             
          CASE(3)
             p3 => ptr(l1:r1,l2:r2,l3:r3,1)
             CALL pnf( &
                  NFMPI_PUT_VARA_DOUBLE_ALL(ncid     &
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , p3 )                &
                  , status, substr )
             
          CASE(4)
             p4 => ptr(l1:r1,l2:r2,l3:r3,l4:r4)
             CALL pnf( &
                  NFMPI_PUT_VARA_DOUBLE_ALL(ncid     &
                  , varid               &
                  , start(1:nrank+it)   &
                  , cnt(1:nrank+it)     &
                  , p4 )                &
                  , status, substr )
             
          END SELECT
          
       END DO

    END SELECT

    ! -----------------------------------------------------
    END SELECT
    ! -----------------------------------------------------

    status = 0

  END SUBROUTINE ch_pnetcdf_write_data
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_pnetcdf_finish_io(status, IOMODE, channel, lclose)

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    INTEGER,          INTENT(IN)  :: IOMODE
    TYPE(t_channel),  POINTER     :: channel ! INTENT(INOUT)
    LOGICAL,          INTENT(IN)  :: lclose

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'ch_pnetcdf_finish_io'
    INTEGER :: ncid

    IF (channel%int%fname(IOMODE) == '') THEN
       status = 0 ! O.K.: FILE NOT OPENED
       RETURN
    END IF

    ncid = channel%int%netcdf(IOMODE)%fileID
    IF (ncid == NC_ID_UNDEF) THEN
       status = 4001 ! UNDEFINED FILE-ID
       RETURN
    END IF

    IF (lclose) THEN
       IF (I_VERBOSE_LEVEL >= 2) THEN ! op_pj_20110803
          WRITE(*,*) substr,': CLOSING  FILE: ',TRIM(channel%int%fname(IOMODE))
       END IF ! op_pj_20110803
       CALL pnf(NFMPI_CLOSE(ncid), status, substr)
       channel%int%netcdf(IOMODE)%fileID = NC_ID_UNDEF
       channel%int%fname(IOMODE) = ''   ! mz_bk_20101119
    ELSE
       IF (L_FLUSH_IOBUFFER) THEN
          IF (I_VERBOSE_LEVEL >= 2) THEN ! op_pj_20110803
             WRITE(*,*) substr,': FLUSHING FILE: ' &
                  ,TRIM(channel%int%fname(IOMODE))
          END IF ! op_pj_20110803
          CALL pnf(NFMPI_SYNC(ncid), status, substr)
       ELSE           ! op_pj_20110803
          status = 0  ! op_pj_20110803
       END IF
    END IF

  END SUBROUTINE ch_pnetcdf_finish_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_pnetcdf_read_data(status, IOMODE, channel, object &
       , ptr, jsnd, i2nd)

    USE messy_main_channel_repr,  ONLY: IRANK, PIOTYPE_SGL, PIOTYPE_IND &
                                      , PIOTYPE_COL

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM

    ! I/O
    INTEGER,                INTENT(OUT)   :: status
    INTEGER,                INTENT(IN)    :: IOMODE
    TYPE(t_channel),        POINTER       :: channel ! INTENT(INOUT)
    TYPE(t_channel_object), POINTER       :: object  ! INTENT(INOUT)
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: ptr
    INTEGER,                INTENT(IN)    :: jsnd    ! OUTPUT DATA TYPE
    INTEGER,                INTENT(IN)    :: i2nd    ! index of 2ndary DATA

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch_pnetcdf_read_data'
    INTEGER                               :: ncid
    INTEGER                               :: varid
    INTEGER                               :: zstat
    INTEGER, DIMENSION(IRANK)             :: nsz
    INTEGER(MPI_OFFSET_KIND), DIMENSION(IRANK) :: start, cnt
    REAL(DP),                     POINTER :: p0 => NULL()
    REAL(DP), DIMENSION(:),       POINTER :: p1 => NULL()
    REAL(DP), DIMENSION(:,:),     POINTER :: p2 => NULL()
    REAL(DP), DIMENSION(:,:,:),   POINTER :: p3 => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: p4 => NULL()
    INTEGER                               :: nrank
    INTEGER                               :: iseg, nseg
    INTEGER                               :: l1,l2,l3,l4
    INTEGER                               :: r1,r2,r3,r4
    REAL(SP), DIMENSION(:,:,:,:), POINTER :: ptrsp

    IF (channel%int%fname(IOMODE) == '') THEN
       status = 0 ! O.K.: FILE NOT OPENED
       RETURN
    END IF

    ncid  = channel%int%netcdf(IOMODE)%fileID
    !
    IF (ncid == NC_ID_UNDEF) THEN
       status = 4001 ! UNDEFINED FILE-ID
       RETURN
    END IF

    IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
         WRITE(*,*) '... parallel netCDF-input: ' &
         , TRIM(channel%int%fname(IOMODE)), ' <- ' &
         , TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE))

    ! GET VARAIABLE ID
    zstat = NFMPI_INQ_VARID( ncid, &
         TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE)), varid)
    IF (zstat /= NF_NOERR) THEN
       ! netCDF VARIABLE DOES NOT EXIST -> CHECK
       SELECT CASE(jsnd)
       CASE(SND_INS)
          IF (object%lrestreq) THEN
             IF (object%int%lign) THEN
                IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
                     WRITE(*,*) substr &
                     ,'    WARNING: REQUIRED RESTART VARIABLE ', &
                     'NOT PRESENT! ' &
                     ,'HOWEVER: IGNORE = T FOR THIS OBJECT'
                status = 0
             ELSE
                status = 3211 ! RESTART VARIABLE REQUIRED BUT NOT PRESENT
             END IF
          ELSE
             IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
                  WRITE(*,*) '    WARNING: '//&
                  &'VARIABLE NOT PRESENT IN RESTART FILE'
             status = 0
          END IF
       CASE DEFAULT ! SND_AVE, SND_STP, SND_MIN, SND_MAX
          IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
               WRITE(*,*) '    WARNING: VARIABLE NOT PRESENT IN RESTART FILE'
          status = 0
       END SELECT
       RETURN ! RETURN IF VARIABLE NOT PRESENT
    END IF
    !
    SELECT CASE(jsnd)
    CASE(SND_UNDEF)
       status = 3203 ! SECONDARY DATA TYPE UNKNOWN
       RETURN
    CASE(SND_INS)
       object%int%netcdf(IOMODE)%varid = varid
    CASE DEFAULT ! SND_AVE, SND_STP, SND_MIN, SND_MAX, ...
       object%int%netcdf(IOMODE)%svarid(i2nd) = varid
    END SELECT

    ! CHECK VARIABLE ATTRIBUTES
    CALL pnetcdf_check_attributes(status, ncid, varid, object%att)
    IF (status /= 0) RETURN

    ! SETUP TEMPORARY MEMORY FOR IMPORT
    nsz(:) = object%repr%pdecomp%shape_out(:)
    ALLOCATE(ptr(nsz(1),nsz(2),nsz(3),nsz(4)), STAT=zstat)
    IF (zstat /= 0) THEN
       status = 1000 ! MEMORY ALLOCATION FAILED
       RETURN
    END IF

    ! LOOP OVER SEGMENTS
    nseg  = object%repr%pdecomp%nseg

    ! INIT
    start(:)  = 1
    cnt(:)    = 1

    nrank = object%repr%rank

    ! -----------------------------------------------------
    ! DATA MODE (independent or collective)
    SELECT CASE(object%repr%pdecomp%piotype)
    ! -----------------------------------------------------

    ! -----------------------------------------------------
    CASE(PIOTYPE_IND, PIOTYPE_SGL)
    ! -----------------------------------------------------

    CALL pnf(NFMPI_BEGIN_INDEP_DATA(ncid), status, substr)
    IF (status /= 0) RETURN

    SELECT CASE(IOMODE)
    CASE(IOMODE_OUT)

       DO iseg=1, nseg

          start(1:nrank) = object%repr%pdecomp%start(iseg,1:nrank)
          cnt(1:nrank)   = object%repr%pdecomp%cnt(iseg,1:nrank)
          
          l1 = object%repr%pdecomp%ml(iseg,1)
          l2 = object%repr%pdecomp%ml(iseg,2)
          l3 = object%repr%pdecomp%ml(iseg,3)
          l4 = object%repr%pdecomp%ml(iseg,4)
          
          r1 = object%repr%pdecomp%mu(iseg,1)
          r2 = object%repr%pdecomp%mu(iseg,2)
          r3 = object%repr%pdecomp%mu(iseg,3)
          r4 = object%repr%pdecomp%mu(iseg,4)

          ALLOCATE(ptrsp(1:(r1-l1+1),1:(r2-l2+1),1:(r3-l3+1),1:(r4-l4+1)))
          
          SELECT CASE(nrank)
          CASE(0)
             CALL pnf( &
                  NFMPI_GET_VARA_REAL(ncid     & 
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , ptrsp(1,1,1,1) ) &
                  , status, substr )
             ptr(1,1,1,1) = REAL(ptrsp(1,1,1,1),DP)

          CASE(1)
             CALL pnf( &
                  NFMPI_GET_VARA_REAL(ncid     & 
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , ptrsp(:,1,1,1) ) &
                  , status, substr )
             ptr(l1:r1,1,1,1) = REAL(ptrsp(:,1,1,1),DP)

          CASE(2)
             p2 => ptr(l1:r1,l2:r2,1,1)
             CALL pnf( &
                  NFMPI_GET_VARA_REAL(ncid     &
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , ptrsp(:,:,1,1) ) &
                  , status, substr )
             ptr(l1:r1,l2:r2,1,1) = REAL(ptrsp(:,:,1,1),DP)

          CASE(3)
             CALL pnf( &
                  NFMPI_GET_VARA_REAL(ncid     &
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , ptrsp(:,:,:,1) ) &
                  , status, substr )
             ptr(l1:r1,l2:r2,l3:r3,1) = REAL(ptrsp(:,:,:,1),DP)
             
          CASE(4)
             p4 => ptr(l1:r1,l2:r2,l3:r3,l4:r4)
             CALL pnf( &
                  NFMPI_GET_VARA_REAL(ncid     &
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , ptrsp(:,:,:,1) ) &
                  , status, substr )
             ptr(l1:r1,l2:r2,l3:r3,l4:r4) = REAL(ptrsp(:,:,:,:),DP)

          END SELECT
          
       END DO
       
    CASE(IOMODE_RST)

       DO iseg=1, nseg

          start(1:nrank) = object%repr%pdecomp%start(iseg,1:nrank)
          cnt(1:nrank)   = object%repr%pdecomp%cnt(iseg,1:nrank)
          
          l1 = object%repr%pdecomp%ml(iseg,1)
          l2 = object%repr%pdecomp%ml(iseg,2)
          l3 = object%repr%pdecomp%ml(iseg,3)
          l4 = object%repr%pdecomp%ml(iseg,4)
          
          r1 = object%repr%pdecomp%mu(iseg,1)
          r2 = object%repr%pdecomp%mu(iseg,2)
          r3 = object%repr%pdecomp%mu(iseg,3)
          r4 = object%repr%pdecomp%mu(iseg,4)
          
          SELECT CASE(nrank)
          CASE(0)
             p0 => ptr(1,1,1,1)
             CALL pnf( &
                  NFMPI_GET_VARA_DOUBLE(ncid     & 
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , (/ p0 /) )          &
                  , status, substr )
          CASE(1)
             p1 => ptr(l1:r1,1,1,1)
             CALL pnf( &
                  NFMPI_GET_VARA_DOUBLE(ncid     & 
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , p1 )                &
                  , status, substr )
             
          CASE(2)
             p2 => ptr(l1:r1,l2:r2,1,1)
             CALL pnf( &
                  NFMPI_GET_VARA_DOUBLE(ncid     &
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , p2 )                &
                  , status, substr )
             
          CASE(3)
             p3 => ptr(l1:r1,l2:r2,l3:r3,1)
             CALL pnf( &
                  NFMPI_GET_VARA_DOUBLE(ncid     &
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , p3 )                &
                  , status, substr )
             
          CASE(4)
             p4 => ptr(l1:r1,l2:r2,l3:r3,l4:r4)
             CALL pnf( &
                  NFMPI_GET_VARA_DOUBLE(ncid     &
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , p4 )                &
                  , status, substr )
             
          END SELECT
          
       END DO

    END SELECT
    
    CALL pnf(NFMPI_END_INDEP_DATA(ncid), status, substr)
    IF (status /= 0) RETURN

    ! -----------------------------------------------------
    CASE(PIOTYPE_COL)
    ! -----------------------------------------------------

    SELECT CASE(IOMODE)
    CASE(IOMODE_OUT)

       DO iseg=1, nseg

          start(1:nrank) = object%repr%pdecomp%start(iseg,1:nrank)
          cnt(1:nrank)   = object%repr%pdecomp%cnt(iseg,1:nrank)
          
          l1 = object%repr%pdecomp%ml(iseg,1)
          l2 = object%repr%pdecomp%ml(iseg,2)
          l3 = object%repr%pdecomp%ml(iseg,3)
          l4 = object%repr%pdecomp%ml(iseg,4)
          
          r1 = object%repr%pdecomp%mu(iseg,1)
          r2 = object%repr%pdecomp%mu(iseg,2)
          r3 = object%repr%pdecomp%mu(iseg,3)
          r4 = object%repr%pdecomp%mu(iseg,4)

          ALLOCATE(ptrsp(1:(r1-l1+1),1:(r2-l2+1),1:(r3-l3+1),1:(r4-l4+1)))
          
          SELECT CASE(nrank)
          CASE(0)
             CALL pnf( &
                  NFMPI_GET_VARA_REAL_ALL(ncid     & 
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , ptrsp(1,1,1,1) ) &
                  , status, substr )
             ptr(1,1,1,1) = REAL(ptrsp(1,1,1,1),DP)

          CASE(1)
             CALL pnf( &
                  NFMPI_GET_VARA_REAL_ALL(ncid     & 
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , ptrsp(:,1,1,1) ) &
                  , status, substr )
             ptr(l1:r1,1,1,1) = REAL(ptrsp(:,1,1,1),DP)

          CASE(2)
             p2 => ptr(l1:r1,l2:r2,1,1)
             CALL pnf( &
                  NFMPI_GET_VARA_REAL_ALL(ncid     &
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , ptrsp(:,:,1,1) ) &
                  , status, substr )
             ptr(l1:r1,l2:r2,1,1) = REAL(ptrsp(:,:,1,1),DP)

          CASE(3)
             CALL pnf( &
                  NFMPI_GET_VARA_REAL_ALL(ncid     &
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , ptrsp(:,:,:,1) ) &
                  , status, substr )
             ptr(l1:r1,l2:r2,l3:r3,1) = REAL(ptrsp(:,:,:,1),DP)
             
          CASE(4)
             p4 => ptr(l1:r1,l2:r2,l3:r3,l4:r4)
             CALL pnf( &
                  NFMPI_GET_VARA_REAL_ALL(ncid     &
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , ptrsp(:,:,:,1) ) &
                  , status, substr )
             ptr(l1:r1,l2:r2,l3:r3,l4:r4) = REAL(ptrsp(:,:,:,:),DP)

          END SELECT
          
       END DO
       
    CASE(IOMODE_RST)

       DO iseg=1, nseg

          start(1:nrank) = object%repr%pdecomp%start(iseg,1:nrank)
          cnt(1:nrank)   = object%repr%pdecomp%cnt(iseg,1:nrank)
          
          l1 = object%repr%pdecomp%ml(iseg,1)
          l2 = object%repr%pdecomp%ml(iseg,2)
          l3 = object%repr%pdecomp%ml(iseg,3)
          l4 = object%repr%pdecomp%ml(iseg,4)
          
          r1 = object%repr%pdecomp%mu(iseg,1)
          r2 = object%repr%pdecomp%mu(iseg,2)
          r3 = object%repr%pdecomp%mu(iseg,3)
          r4 = object%repr%pdecomp%mu(iseg,4)
          
          SELECT CASE(nrank)
          CASE(0)
             p0 => ptr(1,1,1,1)
             CALL pnf( &
                  NFMPI_GET_VARA_DOUBLE_ALL(ncid     & 
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , (/ p0 /) )          &
                  , status, substr )
          CASE(1)
             p1 => ptr(l1:r1,1,1,1)
             CALL pnf( &
                  NFMPI_GET_VARA_DOUBLE_ALL(ncid     & 
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , p1 )                &
                  , status, substr )
             
          CASE(2)
             p2 => ptr(l1:r1,l2:r2,1,1)
             CALL pnf( &
                  NFMPI_GET_VARA_DOUBLE_ALL(ncid     &
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , p2 )                &
                  , status, substr )
             
          CASE(3)
             p3 => ptr(l1:r1,l2:r2,l3:r3,1)
             CALL pnf( &
                  NFMPI_GET_VARA_DOUBLE_ALL(ncid     &
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , p3 )                &
                  , status, substr )
             
          CASE(4)
             p4 => ptr(l1:r1,l2:r2,l3:r3,l4:r4)
             CALL pnf( &
                  NFMPI_GET_VARA_DOUBLE_ALL(ncid     &
                  , varid               &
                  , start(1:nrank)   &
                  , cnt(1:nrank)     &
                  , p4 )                &
                  , status, substr )
             
          END SELECT
          
       END DO

    END SELECT
    
    ! -----------------------------------------------------
    END SELECT
    ! -----------------------------------------------------

    status = 0

  END SUBROUTINE ch_pnetcdf_read_data
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE pnetcdf_write_attribute_list(status, ncid, varid, list)

    USE messy_main_channel_attributes, ONLY: t_attribute_list, t_attribute &
         , TYPE_STRING, TYPE_INTEGER, TYPE_REAL_DP

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, NULL, TRIM, LEN_TRIM, INT

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    INTEGER,                INTENT(IN)  :: ncid
    INTEGER,                INTENT(IN)  :: varid
    TYPE(t_attribute_list), POINTER     :: list

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER     :: substr = 'pnetcdf_write_attribute_list'
    TYPE(t_attribute_list), POINTER :: ai  => NULL()
    TYPE(t_attribute),      POINTER :: att => NULL()

    ai => list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT

       att => ai%this

       SELECT CASE(ai%this%type)
       CASE(TYPE_STRING)
          CALL pnf(NFMPI_PUT_ATT_TEXT(ncid, varid &
               , TRIM(att%name), INT(LEN_TRIM(att%c),MPI_OFFSET_KIND) &
               , TRIM(att%c) ) &
               , status, substr)
       CASE(TYPE_INTEGER)
          CALL pnf(NFMPI_PUT_ATT_INT(ncid, varid &
               , TRIM(att%name), NF_INT, INT(1,MPI_OFFSET_KIND) &
               , att%i) &
               , status, substr)
       CASE(TYPE_REAL_DP)
          CALL pnf(NFMPI_PUT_ATT_DOUBLE(ncid, varid &
               , TRIM(att%name), NF_DOUBLE, INT(1,MPI_OFFSET_KIND) &
               , att%r) &
               , status, substr)
       CASE DEFAULT
          status = 806  ! UNKNOWN ATTRIBUTE TYPE
       END SELECT

       IF (status /= 0) RETURN

       ai => ai%next
    END DO    

    status = 0

  END SUBROUTINE pnetcdf_write_attribute_list
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE pnetcdf_define_dimvar_list(status, ncid, dimid, var, PREC)

    USE messy_main_channel_dimvar, ONLY: t_dimvar_list, t_dimvar

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM

    ! I/O
    INTEGER,             INTENT(OUT) :: status
    INTEGER,             INTENT(IN)  :: ncid
    INTEGER,             INTENT(IN)  :: dimid
    TYPE(t_dimvar_list), POINTER     :: var
    INTEGER,             INTENT(IN)  :: PREC

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'pnetcdf_define_dimvar_list'
    TYPE(t_dimvar_list), POINTER :: ai
    TYPE(t_dimvar),      POINTER :: dv
    INTEGER                      :: varid

    ai => var
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       dv => ai%this

       ! ---------------------------------------------
       CALL pnf( &
            NFMPI_DEF_VAR(ncid, TRIM(dv%name), PREC, 1, dimid, varid) &
            , status, substr)
       IF (status /= 0) RETURN

       ! DIMENSION VARIABLE ATTRIBUTES
       CALL pnetcdf_write_attribute_list(status, ncid, varid, dv%att)
       IF (status /= 0) RETURN
       ! ---------------------------------------------

       ai => ai%next
    END DO

    status = 0

  END SUBROUTINE pnetcdf_define_dimvar_list
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE pnetcdf_write_dimvar_list(status, ncid, var, PREC, start)

    USE messy_main_channel_dimvar, ONLY: t_dimvar_list, t_dimvar

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, REAL, TRIM, INT

    ! I/O
    INTEGER,             INTENT(OUT) :: status
    INTEGER,             INTENT(IN)  :: ncid
    TYPE(t_dimvar_list), POINTER     :: var
    INTEGER,             INTENT(IN)  :: PREC
    INTEGER,             INTENT(IN)  :: start

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'pnetcdf_write_dimvar_list'
    TYPE(t_dimvar_list), POINTER :: ai
    TYPE(t_dimvar),      POINTER :: dv
    INTEGER                      :: varid
    INTEGER(MPI_OFFSET_KIND)     :: zcount, zstart

    ai => var
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       dv => ai%this

       ! ---------------------------------------------
       CALL pnf( &
            NFMPI_INQ_VARID(ncid, TRIM(dv%name), varid) &
            , status, substr )
       IF (status /= 0) RETURN

       zstart = INT(start, MPI_OFFSET_KIND)
       zcount = SIZE(dv%val)

       SELECT CASE(PREC)
       CASE(SP)
          CALL pnf ( &
               NFMPI_PUT_VARA_REAL(ncid, varid &
               , (/zstart/), (/zcount/) &
               , REAL(dv%val, SP)) &
               , status, substr )
       CASE(DP)
          CALL pnf ( &
               NFMPI_PUT_VARA_DOUBLE(ncid, varid &
               , (/zstart/), (/zcount/) &
               , REAL(dv%val, DP)) &
               , status, substr )
       CASE DEFAULT
          status = 1 ! ERROR
       END SELECT

       IF (status /= 0) RETURN
       ! ---------------------------------------------

       ai => ai%next
    END DO

    status = 0

  END SUBROUTINE pnetcdf_write_dimvar_list
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE pnetcdf_check_attributes(status, ncid, varid, attlist)

    USE messy_main_channel_attributes, ONLY: t_attribute_list, t_attribute &
                                           , TYPE_INTEGER, TYPE_REAL_DP    &
                                           , TYPE_STRING
    USE messy_main_constants_mem,      ONLY: STRLEN_ULONG

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM, ABS, TINY

    ! I/O
    INTEGER,          INTENT(OUT)   :: status
    INTEGER,          INTENT(IN)    :: ncid
    INTEGER,          INTENT(IN)    :: varid
    TYPE(t_attribute_list), POINTER :: attlist  ! INTENT(IN)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER     :: substr = 'pnetcdf_check_attributes'
    TYPE(t_attribute_list), POINTER :: list
    TYPE(t_attribute),      POINTER :: att
    INTEGER                         :: zstat
    INTEGER                         :: xtype
    INTEGER(MPI_OFFSET_KIND)        :: len
    INTEGER                         :: ai
    CHARACTER(LEN=STRLEN_ULONG)     :: ac = ''
    REAL(DP)                        :: ar
    LOGICAL                         :: lok

    list => attlist
    DO
       IF (.NOT. ASSOCIATED(list)) EXIT
       att => list%this
       ! -------------------------
       SELECT CASE(att%iflag)
       CASE(AF_RST_NONE)
          list => list%next
          CYCLE
       CASE(AF_RST_CMP)
          IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
               WRITE(*,*) '   CHECKING ATTRIBUTE ''',TRIM(att%name),''' ...'
       CASE(AF_RST_INP)
          IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
               WRITE(*,*) '   READING  ATTRIBUTE ''',TRIM(att%name),''' ...'
       CASE DEFAULT
          list => list%next
          CYCLE
       END SELECT

       ! INQUIRE ATTRIBUTE
       zstat = NFMPI_INQ_ATT(ncid, varid &
            , TRIM(att%name), xtype, len)

       IF (zstat /= NF_NOERR) THEN
          IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
               WRITE(*,*) '   ... *** ERROR *** ATTRIBUTE NOT PRESENT'
          status = 3207  ! MISSING ATTRIBUTE IN RESTART FILE             
          RETURN
       END IF

       ! CHECK LENGTH
       IF (xtype /= NF_CHAR) THEN
          IF (len /= INT(1,MPI_OFFSET_KIND)) THEN
             status = 3208 ! RESTART ATTRIBUTE HAS NON-SCALAR RANK
             RETURN
          END IF
       ELSE
          IF (len > STRLEN_ULONG) THEN
             status = 3209 ! RESTART CHARACTER ATTRIBUTE TOO LONG
             RETURN
          END IF
       END IF
       
       ! GET ATTRIBUTE AND SET TYPE
       SELECT CASE(xtype)
       CASE(NF_BYTE)
          CALL pnf(  NFMPI_GET_ATT_INT1(ncid, varid, &
               TRIM(att%name), ai), status, substr)
       CASE(NF_CHAR)
          CALL pnf(  NFMPI_GET_ATT_TEXT(ncid, varid, &
               TRIM(att%name), ac(1:len+1)), status, substr)
       CASE(NF_SHORT)
          CALL pnf(  NFMPI_GET_ATT_INT2(ncid, varid, &
               TRIM(att%name), ai), status, substr)
       CASE(NF_INT)
          CALL pnf(  NFMPI_GET_ATT_INT(ncid, varid, &
               TRIM(att%name), ai), status, substr)
       CASE(NF_FLOAT)
          CALL pnf(  NFMPI_GET_ATT_REAL(ncid, varid, &
               TRIM(att%name), ar), status, substr)
       CASE(NF_DOUBLE)
          CALL pnf(  NFMPI_GET_ATT_DOUBLE(ncid, varid, &
               TRIM(att%name), ar), status, substr)
       END SELECT
       IF (status /= 0) RETURN
          
       SELECT CASE(att%iflag)
       CASE(AF_RST_CMP) ! COMPARE
          ! CHECK VALUE
          SELECT CASE(att%type)
          CASE(TYPE_INTEGER)
             lok = (att%i == ai)
             IF (.NOT. lok) THEN
                IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
                     WRITE(*,*) '      ... ATTRIBUTE MISMATCH: ',ai &
                     ,' /= ', att%i
             ELSE
                IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
                     WRITE(*,*) '      ... OK: ',ai
             END IF
          CASE(TYPE_REAL_DP)
             lok = (ABS(att%r - ar) < TINY(ar))
             IF (.NOT. lok) THEN
                IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
                     WRITE(*,*) '      ... ATTRIBUTE MISMATCH: ',ar &
                     ,' /= ', att%r
             ELSE
                IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
                     WRITE(*,*) '      ... OK: ',ar
             END IF
          CASE(TYPE_STRING)
             lok = (att%c(1:len) == ac(1:len))
             IF (.NOT. lok) THEN
                IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
                     WRITE(*,*) '      ... ATTRIBUTE MISMATCH: ',ac(1:len) &
                     , ' /= ',att%c(1:len)
             ELSE
                IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
                     WRITE(*,*) '      ... OK: ',TRIM(ac(1:len))
             END IF
          END SELECT
          !
          IF (.NOT. lok) THEN
             status = 3210  ! RESTART ATTRIBUTE MISMATCH
             RETURN
          END IF
          !
       CASE(AF_RST_INP) ! INPUT
          ! INPUT VALUE
          SELECT CASE(att%type)
          CASE(TYPE_INTEGER)
             att%i = ai
             IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
                  WRITE(*,*) '      ... ', att%i
          CASE(TYPE_REAL_DP)
             att%r = ar
             IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
                  WRITE(*,*) '      ... ', att%r
          CASE(TYPE_STRING)
             att%c = ac
             IF (p_parallel_io .AND. (I_VERBOSE_LEVEL >= 1)) &
                  WRITE(*,*) '      ... ', TRIM(att%c)
          END SELECT
          !
       CASE DEFAULT
          !
          !
       END SELECT
       
       ! -------------------------
       list => list%next
    END DO

    ! RETURN
    status = 0

  END SUBROUTINE pnetcdf_check_attributes
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE pnf(sin, sout, substr)

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,          INTENT(IN)  :: sin
    INTEGER,          INTENT(OUT) :: sout
    CHARACTER(LEN=*), INTENT(IN)  :: substr

    ! mz_pj_20080807+
    ! LOCAL
    INTEGER :: iou
    LOGICAL :: opened
    CHARACTER(LEN=22) :: errfname = ''
    ! mz_pj_20080807-

    IF (sin /= NF_NOERR) THEN
       WRITE(*,*) TRIM(substr),': *** parallel-netCDF ERROR on PE ',p_pe &
            ,' : ',NFMPI_STRERROR(sin)
       sout = 4002 ! PARALLEL NETCDF ERROR

       ! mz_pj_20080807+
       DO iou=100,300
          INQUIRE(unit=iou,opened=opened)
          IF (.NOT.opened) EXIT
       END DO
       WRITE(errfname,'(a13,a1,i4.4)') 'ERROR.pnetcdf','.',p_pe
       OPEN(iou, FILE=TRIM(errfname), STATUS='UNKNOWN')
       WRITE(iou,*) TRIM(substr),': *** parallel-netCDF ERROR on PE',p_pe &
            ,' : ',NFMPI_STRERROR(sin)
       CLOSE(iou)
       ! mz_pj_20080807-
    ELSE
       sout = 0
    ENDIF
    
  END SUBROUTINE pnf
  ! -------------------------------------------------------------------

! op_pj_20091124+
  ! -------------------------------------------------------------------
  SUBROUTINE ch_pnetcdf_read_nml_ctrl(status, iou)

    ! MODULE ROUTINE (CORE)
    !
    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, DLR-IPA, Nov 2009

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! mz_pj_20080118: OUT_PREC added
    NAMELIST /CTRL_PNETCDF/ MPI_IO_HINT

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch_pnetcdf_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CTRL_PNETCDF', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_PNETCDF, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_PNETCDF', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE ch_pnetcdf_read_nml_ctrl
  ! -------------------------------------------------------------------
! op_pj_20091124-

#endif
! **********************************************************************
END MODULE messy_main_channel_pnetcdf
! **********************************************************************
