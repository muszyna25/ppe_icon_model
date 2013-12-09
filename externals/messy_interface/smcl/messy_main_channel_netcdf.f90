! **********************************************************************
MODULE messy_main_channel_netcdf
! **********************************************************************

  ! MESSY DATA TRANSFER AND EXPORT INTERFACE (MEMORY MANAGEMENT)
  !
  ! Author: Patrick Joeckel, MPICH, May 2005

  USE messy_main_channel

  USE netcdf

  IMPLICIT NONE
  PRIVATE

#if defined (__SX__) || defined (__uxp__)
  INTEGER :: initialsize = 33554432      ! that's 32 MByte   
  INTEGER :: chunksize   = 33554432      ! too
#else
  ! mz_kk_20070305+
  !INTEGER :: initialsize =    32768      ! that's 32 kByte   
  !INTEGER :: chunksize   =    32768      ! too
  INTEGER :: initialsize =    1024*1024      ! that's 1 MByte   
  INTEGER :: chunksize   =    1024*1024      ! too
  ! mz_kk_20070305-
#endif

  INTEGER, SAVE :: NCOUT_PREC = NF90_FLOAT ! mz_pj_20080118 default

  PUBLIC :: ch_netcdf_init_rst
  PUBLIC :: ch_netcdf_init_io
  PUBLIC :: ch_netcdf_write_header
  PUBLIC :: ch_netcdf_write_time
  PUBLIC :: ch_netcdf_write_data
  PUBLIC :: ch_netcdf_finish_io
  !
  PUBLIC :: ch_netcdf_read_data
  !
  !PRIVATE :: netcdf_write_attribute_list
  !PRIVATE :: netcdf_define_dimvar_list
  !PRIVATE :: netcdf_write_dimvar_list
  !PRIVATE :: netcdf_check_attributes
  !PRIVATE :: nf

CONTAINS

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_netcdf_init_rst(status, fname, att)

    USE messy_main_channel_attributes, ONLY: t_attribute_list

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    CHARACTER(LEN=*),       INTENT(IN)  :: fname
    TYPE(t_attribute_list), POINTER     :: att 

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch_netcdf_init_rst'
    INTEGER                     :: ncid

    CALL nf(NF90_OPEN(TRIM(fname), NF90_NOWRITE, ncid), status, substr)
    IF (status /= 0) RETURN

    CALL netcdf_check_attributes(status, ncid, NF90_GLOBAL, att) 
    IF (status /= 0) RETURN

    CALL nf(NF90_CLOSE(ncid), status, substr)
    IF (status /= 0) RETURN

  END SUBROUTINE ch_netcdf_init_rst
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_netcdf_init_io(status, IOMODE, channel, AMODE, att)

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
    CHARACTER(LEN=*), PARAMETER    :: substr = 'ch_netcdf_init_io'
    INTEGER                        :: ncid
    INTEGER                        :: ifmode
    LOGICAL                        :: lexist
    INTEGER                        :: nc_cmode
    LOGICAL, SAVE                  :: lfirst = .TRUE. ! mz_pj_20080118

    ! mz_pj_20080118+
    IF (lfirst) THEN
       SELECT CASE(OUT_PREC(FTYPE_NETCDF))
       CASE(1)
          NCOUT_PREC = NF90_FLOAT
       CASE(2)
          NCOUT_PREC = NF90_DOUBLE
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
          IF (I_VERBOSE_LEVEL >= 1) THEN  ! op_pj_20110803
             WRITE(*,*) substr,': OPENING  FILE: ' &
                  , TRIM(channel%int%fname(IOMODE))
          ENDIF                           ! op_pj_20110803
          !
          INQUIRE(file = TRIM(channel%int%fname(IOMODE)), exist = lexist)
          IF (.NOT.lexist) THEN
             IF (channel%int%lrestreq) THEN
                IF (channel%int%lign) THEN
                   ! IGNORE ALL MISSING RESTART FIELDS
                   IF (I_VERBOSE_LEVEL >= 1) THEN  ! op_pj_20110803
                      WRITE(*,*) substr,':       WARNING: REQUIRED RESTART ', &
                           'FILE NOT PRESENT ... ALL OBJECTS TO BE IGNORED'
                   END IF                          ! op_pj_20110803
                   ! RESET: PREVENT FROM CLOSING AND READING
                   channel%int%fname(IOMODE) = ''
                   status = 0
                ELSE
                   status = 3206 ! RESTART FILE REQUIRED BUT NOT PRESENT
                END IF
             ELSE
                IF (I_VERBOSE_LEVEL >= 1) THEN  ! op_pj_20110803
                   WRITE(*,*) substr,':       WARNING: RESTART ', &
                        'FILE NOT PRESENT ... NOT REQUIRED'
                END IF                          ! op_pj_20110803
                ! RESET: PREVENT FROM CLOSING AND READING
                channel%int%fname(IOMODE) = ''
                status = 0
             END IF
             RETURN
          END IF
          !
          CALL nf(NF90_OPEN(TRIM(channel%int%fname(IOMODE)) &
               ,NF90_NOWRITE, ncid) &
               , status, substr)
          IF (status /= 0) RETURN
          channel%int%netcdf(IOMODE)%fileID = ncid
          !
          ! OPTIONAL SPECIAL ATTRIBUTES
          IF (PRESENT(att)) THEN
             CALL netcdf_check_attributes(status, ncid, NF90_GLOBAL, att)
             IF (status /= 0) RETURN
          END IF
          CALL netcdf_check_attributes(status, ncid, NF90_GLOBAL, GATT)
          IF (status /= 0) RETURN
          CALL netcdf_check_attributes(status, ncid, NF90_GLOBAL, channel%att)
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
          CALL nf(NF90_CLOSE(channel%int%netcdf(IOMODE)%fileID) &
               , status, substr)
          IF (status /= 0) RETURN
          channel%int%netcdf(IOMODE)%fileID = NC_ID_UNDEF
       END IF
       !
       ! ---------------------------------------------------------
       ! OPEN 'NEW' FILE
       ! ---------------------------------------------------------
       IF (I_VERBOSE_LEVEL >= 2) THEN ! op_pj_20110803
          WRITE(*,*) substr,': CREATING FILE: ',TRIM(channel%int%fname(IOMODE))
       END IF ! op_pj_20110803
!qqq+
!       IF (IOMODE == IOMODE_RST) THEN
          nc_cmode = IOR(NF90_CLOBBER,NF90_64BIT_OFFSET)
!       ELSE
!          nc_cmode = NF90_CLOBBER
!       END IF
!qqq-
       CALL nf( &
            NF90_CREATE( TRIM(channel%int%fname(IOMODE)) &
!            , NF90_CLOBBER                           &  ! overwrite
            , nc_cmode                               &
            , ncid                                   &
            , initialsize, chunksize )               &
            , status, substr)
       IF (status /= 0) RETURN
       channel%int%netcdf(IOMODE)%fileID = ncid
       ! SWITCH OFF 'FILL MODE'
       CALL nf(NF90_SET_FILL(ncid, NF90_NOFILL, ifmode), status, substr)
       IF (status /= 0) RETURN
       !
    END SELECT

    status = 0

  END SUBROUTINE ch_netcdf_init_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_netcdf_write_header(status, IOMODE, channel, dim_time, att)

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
    CHARACTER(LEN=*), PARAMETER    :: substr = 'ch_netcdf_write_header'
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
         , c = nf90_inq_libvers(), loverwrite=.TRUE.)
    IF (status /= 0) RETURN
    ! mz_pj_20090115-

    ! ---------------------------------------------------------
    ! WRITE GLOBAL ATTRIBUTES
    ! ---------------------------------------------------------
    ! - COMMON TO ALL CHANNELS
    CALL netcdf_write_attribute_list(status, ncid, NF90_GLOBAL, GATT)
    IF (status /= 0) RETURN
    ! - CHANNEL SPECIFIC
    CALL netcdf_write_attribute_list(status, ncid, NF90_GLOBAL, &
         channel%att)
    IF (status /= 0) RETURN
    ! - OPTIONAL SPECIAL ATRIBUTES
    IF (PRESENT(att)) THEN
       CALL netcdf_write_attribute_list(status, ncid, NF90_GLOBAL, att)
       IF (status /= 0) RETURN
    END IF

    ! ---------------------------------------------------------
    ! DEFINE UNLIMITED DIMENSION WITH ATTRIBUTES
    ! ---------------------------------------------------------
    IF (IOMODE == IOMODE_OUT) THEN
       CALL nf( &
            NF90_DEF_DIM( ncid, TRIM(dim_time%name) &
            , NF90_UNLIMITED                        &
            , channel%int%netcdf(IOMODE)%dimid_time)    &
            , status, substr)
       IF (status /= 0) RETURN      
       CALL netcdf_define_dimvar_list(status       &
            , ncid                                 & ! netCDF file ID
            , channel%int%netcdf(IOMODE)%dimid_time & ! dimension ID
            , dim_time%var                         & ! list of variables
            , NF90_DOUBLE                          & ! precision
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
!!$       XNF90_PREC = NF90_FLOAT
       XNF90_PREC = NCOUT_PREC
       ! mz_pj_20080118-
    ELSE
       XNF90_PREC = NF90_DOUBLE
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
                CALL nf( &
                     NF90_DEF_DIM( ncid, TRIM(object%repr%dim(i)%ptr%name) &
                     , object%repr%dim(i)%ptr%len                          &
                     , dimid(object%repr%dim(i)%ptr%id) )                  &
                     , status, substr)
                IF (status /= 0) RETURN
                !
                CALL netcdf_define_dimvar_list(status    &
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
!!$          XNF90_PREC = NF90_FLOAT
          XNF90_PREC = NCOUT_PREC
          ! mz_pj_20080118-
       ELSE
          it = 0
          XNF90_PREC = NF90_DOUBLE
       END IF

       ! DEFINE PRIMARY AND SECONDARY DATA
       DO jsnd=1, SND_MAXLEN

          IF (.NOT. object%int%lexp(jsnd, IOMODE)) CYCLE ! NO OUTPUT

          CALL nf( &
               NF90_DEF_VAR(ncid           & ! netCDF file ID
               , TRIM(object%name)//&
               &TRIM(SND_TEXT(jsnd,IOMODE)) & ! variable name
               , XNF90_PREC                & ! precision
               , dimvec(1:nrank+it)        & ! netCDF dim-IDs
               , varid)                    & ! netCDF var-ID
               , status, substr )
          IF (status /= 0) RETURN
          IF (jsnd == SND_INS) then
             object%int%netcdf(IOMODE)%varid = varid        ! primary data
          ELSE
             i2nd = object%int%i2nd(jsnd)
             object%int%netcdf(IOMODE)%svarid(i2nd) = varid ! 2ndary data
          END IF
          ! WRITE VARIABLE ATTRIBUTES
          CALL netcdf_write_attribute_list(status, ncid, varid, &
               object%att)
          ! WRITE SPECIAL ATTRIBUTES FOR SECONDARY DATA
          IF ((jsnd == SND_CNT) .OR. (jsnd == SND_CAV)) THEN
             CALL nf(NF90_PUT_ATT(ncid, varid, 'range', object%io%range) &
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
    CALL nf(NF90_ENDDEF(ncid), status, substr)
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

                CALL netcdf_write_dimvar_list(status     &
                     , ncid                              & ! netCDF file ID
                     , object%repr%dim(i)%ptr%var        & ! list of variables
                     , PREC                              & ! precision
                     , 1                                 & ! start
                     )
                IF (status /= 0) RETURN

                ldim_out(object%repr%dim(i)%ptr%id) = .FALSE. ! not twice
             END IF
          END IF
       END DO
       ! -------------------------------------------------------

       le => le%next
    END DO object_loop3

    ! ---------------------------------------------------------
    ! FLUSH netCDF BUFFER
    ! ---------------------------------------------------------
    IF (L_FLUSH_IOBUFFER) THEN
       CALL nf(NF90_SYNC(ncid), status, substr)
       IF (status /= 0) RETURN
    END IF

    ! ---------------------------------------------------------
    ! CLEAN UP
    ! ---------------------------------------------------------
    DEALLOCATE(dimid)
    DEALLOCATE(ldim_out)

  END SUBROUTINE ch_netcdf_write_header
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_netcdf_write_time(status, channel, dim_time)

    USE messy_main_channel_dimensions, ONLY: t_dimension

    IMPLICIT NONE
    
    ! I/O
    INTEGER,           INTENT(OUT) :: status
    TYPE(t_channel),   POINTER     :: channel ! INTENT(INOUT)
    TYPE(t_dimension), POINTER     :: dim_time

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER    :: substr = 'ch_netcdf_write_time'
    INTEGER :: ncid
    INTEGER :: outstep

    ncid    = channel%int%netcdf(IOMODE_OUT)%fileID
    IF (ncid == NC_ID_UNDEF) THEN
       status = 4001 ! UNDEFINED FILE-ID
       RETURN
    END IF

    outstep = channel%int%ntpfcnt

    CALL netcdf_write_dimvar_list(status  &
         , ncid                           & ! netCDF file ID
         , dim_time%var                   & ! list of variables
         , DP                             & ! precision
         , outstep                        & ! start
         )
    
  END SUBROUTINE ch_netcdf_write_time
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_netcdf_write_data(status, IOMODE, channel, object &
       , ptr, jsnd, i2nd)

    USE messy_main_channel_repr,  ONLY: IRANK

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE, TRIM

    ! I/O
    INTEGER,                INTENT(OUT)   :: status
    INTEGER,                INTENT(IN)    :: IOMODE
    TYPE(t_channel),        POINTER       :: channel ! INTENT(INOUT)
    TYPE(t_channel_object), POINTER       :: object  ! INTENT(INOUT)
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: ptr
    INTEGER,                INTENT(IN)    :: jsnd    ! OUTPUT DATA TYPE
    INTEGER,                INTENT(IN)    :: i2nd    ! index of 2ndary DATA

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch_netcdf_write_data'
    INTEGER                               :: ncid
    INTEGER                               :: varid
    INTEGER                               :: nrank
    INTEGER, DIMENSION(IRANK+1)           :: start, cnt
    INTEGER                               :: i
    REAL(DP),                     POINTER :: p0 => NULL()
    REAL(DP), DIMENSION(:),       POINTER :: p1 => NULL()
    REAL(DP), DIMENSION(:,:),     POINTER :: p2 => NULL()
    REAL(DP), DIMENSION(:,:,:),   POINTER :: p3 => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: p4 => NULL()
    INTEGER                               :: it

    IF (I_VERBOSE_LEVEL >= 2) THEN ! op_pj_20110803
       SELECT CASE(IOMODE)
       CASE(IOMODE_OUT)
          WRITE(*,*) '... netCDF-output: '         &
               , TRIM(channel%int%fname(IOMODE))   & 
               , ' (',channel%int%ntpfcnt,') <- '  &
               , TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE))

       CASE(IOMODE_RST)
          WRITE(*,*) '... netCDF-output: '         &
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

    start(:) = 1
    cnt(:)   = 1
    DO i=1, nrank
       cnt(i) = SIZE(ptr, i)
    END DO

    SELECT CASE(IOMODE)
    CASE(IOMODE_OUT)
       it = 1
       start(nrank + 1) = channel%int%ntpfcnt   ! OUTPUT STEP IN netCDF FILE
       cnt(nrank + 1)   = 1                     ! 1 TIME STEP
    CASE(IOMODE_RST)
       it = 0                                   ! NO TIME AXIS IN RESTART FILE
    END SELECT

    SELECT CASE(nrank)
    CASE(0)
       p0 => ptr(1,1,1,1)
       CALL nf( &
            NF90_PUT_VAR(ncid     & 
            , varid               &
            , (/ p0 /)            &
            , start(1:nrank+it)   &
            , cnt(1:nrank+it)   ) &
            , status, substr )
    CASE(1)
       p1 => ptr(:,1,1,1)
       CALL nf( &
            NF90_PUT_VAR(ncid     & 
            , varid               &
            , p1                  &
            , start(1:nrank+it)   &
            , cnt(1:nrank+it)   ) &
            , status, substr )

    CASE(2)
       p2 => ptr(:,:,1,1)
       CALL nf( &
            NF90_PUT_VAR(ncid     &
            , varid               &
            , p2                  &
            , start(1:nrank+it)   &
            , cnt(1:nrank+it)   ) &
            , status, substr )
       
    CASE(3)
       p3 => ptr(:,:,:,1)
       CALL nf( &
            NF90_PUT_VAR(ncid     &
            , varid               &
            , p3                  &
            , start(1:nrank+it)   &
            , cnt(1:nrank+it)   ) &
            , status, substr )

    CASE(4)
       p4 => ptr(:,:,:,:)
       CALL nf( &
            NF90_PUT_VAR(ncid     &
            , varid               &
            , p4                  &
            , start(1:nrank+it)   &
            , cnt(1:nrank+it)   ) &
            , status, substr )

    END SELECT

  END SUBROUTINE ch_netcdf_write_data
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_netcdf_finish_io(status, IOMODE, channel, lclose)

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    INTEGER,          INTENT(IN)  :: IOMODE
    TYPE(t_channel),  POINTER     :: channel ! INTENT(INOUT)
    LOGICAL,          INTENT(IN)  :: lclose

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'ch_netcdf_finish_io'
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
       CALL nf(NF90_CLOSE(ncid), status, substr)
       channel%int%netcdf(IOMODE)%fileID = NC_ID_UNDEF
       channel%int%fname(IOMODE) = '' ! mz_bk_20101119
    ELSE
       IF (L_FLUSH_IOBUFFER) THEN
          IF (I_VERBOSE_LEVEL >= 2) THEN ! op_pj_20110803
             WRITE(*,*) substr,': FLUSHING FILE: ' &
                  ,TRIM(channel%int%fname(IOMODE))
          END IF ! op_pj_20110803
          CALL nf(NF90_SYNC(ncid), status, substr)
       ELSE
          status = 0 ! mz_pj_20081030
       END IF
    END IF

  END SUBROUTINE ch_netcdf_finish_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_netcdf_read_data(status, IOMODE, channel, object &
       , ptr, jsnd, i2nd)

    USE messy_main_channel_repr,  ONLY: IRANK

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
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch_netcdf_read_data'
    INTEGER                               :: ncid
    INTEGER                               :: varid
    INTEGER                               :: zstat
    INTEGER, DIMENSION(IRANK)             :: nsz

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

    IF (I_VERBOSE_LEVEL >= 1) & ! op_pj_20110803
         WRITE(*,*) '... netCDF-input: ' &
         , TRIM(channel%int%fname(IOMODE)), ' <- ' &
         , TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE))

    ! GET VARAIABLE ID
    zstat = NF90_INQ_VARID( ncid, &
         TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE)), varid)
    IF (zstat /= NF90_NOERR) THEN
       ! netCDF VARIABLE DOES NOT EXIST -> CHECK
       SELECT CASE(jsnd)
       CASE(SND_INS)
          IF (object%lrestreq) THEN
             IF (object%int%lign) THEN
                IF (I_VERBOSE_LEVEL >= 1) & ! op_pj_20110803
                     WRITE(*,*) substr &
                     ,'    WARNING: REQUIRED RESTART VARIABLE ', &
                     'NOT PRESENT! ' &
                     ,'HOWEVER: IGNORE = T FOR THIS OBJECT'
                status = 0
             ELSE
                status = 3211 ! RESTART VARIABLE REQUIRED BUT NOT PRESENT
             END IF
          ELSE
             IF (I_VERBOSE_LEVEL >= 1) & ! op_pj_20110803
                  WRITE(*,*) '    WARNING: VARIABLE NOT PRESENT IN RESTART FILE'
             status = 0
          END IF
       CASE DEFAULT ! SND_AVE, SND_STP, SND_MIN, SND_MAX
          IF (I_VERBOSE_LEVEL >= 1) & ! op_pj_20110803
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
    CALL netcdf_check_attributes(status, ncid, varid, object%att)
    IF (status /= 0) RETURN

    ! SETUP TEMPORARY MEMORY FOR IMPORT
    nsz(:) = object%repr%shape_out(:)
    ALLOCATE(ptr(nsz(1),nsz(2),nsz(3),nsz(4)), STAT=zstat)
    IF (zstat /= 0) THEN
       status = 1000 ! MEMORY ALLOCATION FAILED
       RETURN
    END IF

    ! IMPORT DATA
    CALL nf( &
         NF90_GET_VAR(ncid, varid, ptr) &
         , status, substr)
    IF (status /= 0) RETURN
    
    status = 0

  END SUBROUTINE ch_netcdf_read_data
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE netcdf_write_attribute_list(status, ncid, varid, list)

    USE messy_main_channel_attributes, ONLY: t_attribute_list, t_attribute &
         , TYPE_STRING, TYPE_INTEGER, TYPE_REAL_DP

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, NULL, TRIM

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    INTEGER,                INTENT(IN)  :: ncid
    INTEGER,                INTENT(IN)  :: varid
    TYPE(t_attribute_list), POINTER     :: list

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER     :: substr = 'netcdf_write_attribute_list'
    TYPE(t_attribute_list), POINTER :: ai  => NULL()
    TYPE(t_attribute),      POINTER :: att => NULL()

    ai => list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT

       att => ai%this

       SELECT CASE(ai%this%type)
       CASE(TYPE_STRING)
          CALL nf(NF90_PUT_ATT(ncid, varid, TRIM(att%name), TRIM(att%c)) &
               , status, substr)
       CASE(TYPE_INTEGER)
          CALL nf(NF90_PUT_ATT(ncid, varid, TRIM(att%name), att%i) &
               , status, substr)
       CASE(TYPE_REAL_DP)
          CALL nf(NF90_PUT_ATT(ncid, varid, TRIM(att%name), att%r) &
               , status, substr)
       CASE DEFAULT
          status = 806  ! UNKNOWN ATTRIBUTE TYPE
       END SELECT

       IF (status /= 0) RETURN

       ai => ai%next
    END DO    

    status = 0

  END SUBROUTINE netcdf_write_attribute_list
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE netcdf_define_dimvar_list(status, ncid, dimid, var, PREC)

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
    CHARACTER(LEN=*), PARAMETER  :: substr = 'netcdf_define_dimvar_list'
    TYPE(t_dimvar_list), POINTER :: ai
    TYPE(t_dimvar),      POINTER :: dv
    INTEGER                      :: varid

    ai => var
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       dv => ai%this

       ! ---------------------------------------------
       CALL nf( &
            NF90_DEF_VAR(ncid, TRIM(dv%name), PREC, dimid, varid) &
            , status, substr)
       IF (status /= 0) RETURN

       ! DIMENSION VARIABLE ATTRIBUTES
       CALL netcdf_write_attribute_list(status, ncid, varid, dv%att)
       IF (status /= 0) RETURN
       ! ---------------------------------------------

       ai => ai%next
    END DO

    status = 0

  END SUBROUTINE netcdf_define_dimvar_list
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE netcdf_write_dimvar_list(status, ncid, var, PREC, start)

    USE messy_main_channel_dimvar, ONLY: t_dimvar_list, t_dimvar

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, REAL, TRIM

    ! I/O
    INTEGER,             INTENT(OUT) :: status
    INTEGER,             INTENT(IN)  :: ncid
    TYPE(t_dimvar_list), POINTER     :: var
    INTEGER,             INTENT(IN)  :: PREC
    INTEGER,             INTENT(IN)  :: start

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'netcdf_write_dimvar_list'
    TYPE(t_dimvar_list), POINTER :: ai
    TYPE(t_dimvar),      POINTER :: dv
    INTEGER                      :: varid

    ai => var
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       dv => ai%this

       ! ---------------------------------------------
       CALL nf( &
            NF90_INQ_VARID(ncid, TRIM(dv%name), varid) &
            , status, substr )
       IF (status /= 0) RETURN

       SELECT CASE(PREC)
       CASE(SP)
       CALL nf ( &
            NF90_PUT_VAR(ncid, varid, REAL(dv%val, SP), (/start/)) &
            , status, substr )
       CASE(DP)
       CALL nf ( &
            NF90_PUT_VAR(ncid, varid, REAL(dv%val, DP), (/start/)) &
            , status, substr )
       CASE DEFAULT
          status = 1 ! ERROR
       END SELECT

       IF (status /= 0) RETURN
       ! ---------------------------------------------

       ai => ai%next
    END DO

    status = 0

  END SUBROUTINE netcdf_write_dimvar_list
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE netcdf_check_attributes(status, ncid, varid, attlist)

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
    CHARACTER(LEN=*), PARAMETER     :: substr = 'netcdf_check_attributes'
    TYPE(t_attribute_list), POINTER :: list
    TYPE(t_attribute),      POINTER :: att
    INTEGER                         :: zstat
    INTEGER                         :: xtype, len, attnum
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
          IF (I_VERBOSE_LEVEL >= 1) & ! op_pj_20110803
               WRITE(*,*) '   CHECKING ATTRIBUTE ''',TRIM(att%name),''' ...'
       CASE(AF_RST_INP)
          IF (I_VERBOSE_LEVEL >= 1) & ! op_pj_20110803
               WRITE(*,*) '   READING  ATTRIBUTE ''',TRIM(att%name),''' ...'
       CASE DEFAULT
          list => list%next
          CYCLE
       END SELECT

       ! INQUIRE ATTRIBUTE
       zstat = NF90_INQUIRE_ATTRIBUTE(ncid, varid &
            , TRIM(att%name), xtype, len, attnum)

       IF (zstat /= NF90_NOERR) THEN
          IF (I_VERBOSE_LEVEL >= 0) & ! op_pj_20110803
               WRITE(*,*) '   ... *** ERROR *** ATTRIBUTE NOT PRESENT'
          status = 3207  ! MISSING ATTRIBUTE IN RESTART FILE             
          RETURN
       END IF

       ! CHECK LENGTH
       IF (xtype /= NF90_CHAR) THEN
          IF (len /= 1) THEN
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
       CASE(NF90_BYTE)
          CALL nf(  NF90_GET_ATT(ncid, varid, &
               TRIM(att%name), ai), status, substr)
       CASE(NF90_CHAR)
          CALL nf(  NF90_GET_ATT(ncid, varid, &
               TRIM(att%name), ac(1:len+1)), status, substr)
       CASE(NF90_SHORT)
          CALL nf(  NF90_GET_ATT(ncid, varid, &
               TRIM(att%name), ai), status, substr)
       CASE(NF90_INT)
          CALL nf(  NF90_GET_ATT(ncid, varid, &
               TRIM(att%name), ai), status, substr)
       CASE(NF90_FLOAT)
          CALL nf(  NF90_GET_ATT(ncid, varid, &
               TRIM(att%name), ar), status, substr)
       CASE(NF90_DOUBLE)
          CALL nf(  NF90_GET_ATT(ncid, varid, &
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
                IF (I_VERBOSE_LEVEL >= 1) &
                     WRITE(*,*) '      ... ATTRIBUTE MISMATCH: ',ai &
                     ,' /= ', att%i
             ELSE
                IF (I_VERBOSE_LEVEL >= 1) &
                     WRITE(*,*) '      ... OK: ',ai
             END IF
          CASE(TYPE_REAL_DP)
             lok = (ABS(att%r - ar) < TINY(ar))
             IF (.NOT. lok) THEN
                IF (I_VERBOSE_LEVEL >= 1) &
                     WRITE(*,*) '      ... ATTRIBUTE MISMATCH: ',ar &
                     ,' /= ', att%r
             ELSE
                IF (I_VERBOSE_LEVEL >= 1) &
                     WRITE(*,*) '      ... OK: ',ar
             END IF
          CASE(TYPE_STRING)
             lok = (att%c(1:len) == ac(1:len))
             IF (.NOT. lok) THEN
                IF (I_VERBOSE_LEVEL >= 1) &
                     WRITE(*,*) '      ... ATTRIBUTE MISMATCH: ',ac(1:len) &
                     , ' /= ',att%c(1:len)
             ELSE
                IF (I_VERBOSE_LEVEL >= 1) &
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
             IF (I_VERBOSE_LEVEL >= 1) WRITE(*,*) '      ... ', att%i
          CASE(TYPE_REAL_DP)
             att%r = ar
             IF (I_VERBOSE_LEVEL >= 1) WRITE(*,*) '      ... ', att%r
          CASE(TYPE_STRING)
             att%c = ac
             IF (I_VERBOSE_LEVEL >= 1) WRITE(*,*) '      ... ', TRIM(att%c)
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

  END SUBROUTINE netcdf_check_attributes
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE nf(sin, sout, substr)

!!$#ifdef NAG
!!$  USE f90_unix,         ONLY: flush
!!$#endif
!!$
!!$#ifdef __ibm__
!!$#define flush flush_
!!$  USE xlfutility,       ONLY: flush
!!$#endif

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
    ! mz_pj_20080807-

    IF (sin /= NF90_NOERR) THEN

       WRITE(*,*) TRIM(substr),': *** netCDF ERROR: ',NF90_STRERROR(sin)
       sout = 4000 ! NETCDF ERROR

       ! mz_pj_20080807+
       DO iou=100,300
          INQUIRE(unit=iou,opened=opened)
          IF (.NOT.opened) EXIT
       END DO
       OPEN(iou, FILE='ERROR.netcdf', STATUS='UNKNOWN')
       WRITE(iou,*) TRIM(substr),': *** netCDF ERROR: ',NF90_STRERROR(sin)
       CLOSE(iou)
       ! mz_pj_20080807-
    ELSE
       sout = 0
    ENDIF
    
  END SUBROUTINE nf
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_main_channel_netcdf
! **********************************************************************
 
