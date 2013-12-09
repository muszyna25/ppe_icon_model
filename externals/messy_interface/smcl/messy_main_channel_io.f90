! **********************************************************************
MODULE messy_main_channel_io
! **********************************************************************

  ! MESSY DATA TRANSFER AND EXPORT INTERFACE (MEMORY MANAGEMENT)
  !
  ! Author: Patrick Joeckel, MPICH, May 2005

  USE messy_main_channel
  
  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  PUBLIC :: initialize_parallel_io ! FOR PARALLEL I/O
  PUBLIC :: channel_init_restart ! INITIALIZE DATE FROM RESTART
  PUBLIC :: channel_init_io      ! OPEN FILE FOR READ/WRITE
  PUBLIC :: channel_write_header ! INITIALIZE OUTPUT FILE (WRITE HEADER)
  PUBLIC :: channel_write_time   ! WRITE TIME INFORMATION
  PUBLIC :: channel_write_data   ! WRITE DATA TO OUTPUT FILE
  PUBLIC :: channel_finish_io    ! FLUSH / CLOSE FILE
  !
  PUBLIC :: channel_read_data    ! READ DATA FROM RESTART FILE

CONTAINS

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE initialize_parallel_io(status, p_pe, p_io, p_all_comm)

    USE messy_main_channel_pnetcdf,     ONLY: ch_pnetcdf_init_pio

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: p_pe
    INTEGER, INTENT(IN)  :: p_io
    INTEGER, INTENT(IN)  :: p_all_comm

    CALL ch_pnetcdf_init_pio(status, p_pe, p_io, p_all_comm)

  END SUBROUTINE initialize_parallel_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_init_restart(status, lp, lp_io, fname_base, rstatt)

    USE messy_main_channel_attributes, ONLY: t_attribute_list
    USE messy_main_channel_netcdf,     ONLY: ch_netcdf_init_rst
#ifdef PNETCDF
    USE messy_main_channel_pnetcdf,    ONLY: ch_pnetcdf_init_rst
#endif
    USE messy_main_channel,            ONLY: new_attribute, get_attribute 

    IMPLICIT NONE

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    LOGICAL,                INTENT(OUT) :: lp
    LOGICAL,                INTENT(IN)  :: lp_io
    CHARACTER(LEN=*),       INTENT(IN)  :: fname_base
    TYPE(t_attribute_list), POINTER     :: rstatt     ! INTENT(IN)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'channel_init_restart'
    LOGICAL                     :: lex=.false.
    INTEGER                     :: i

    ! CHECK IF RESTART FILE IS PRESENT
    DO i=1,  FTYPE_MAXIMUM 
       IF (lp_io .AND. (I_VERBOSE_LEVEL >= 1)) &
            WRITE(*,*) substr,': checking for file '//&
            &TRIM(fname_base)//TRIM(FTYPE_EXT_TEXT(i))//' ...'
       INQUIRE(file=TRIM(fname_base)//TRIM(FTYPE_EXT_TEXT(i)), exist=lex)
       IF (lex) EXIT
    END DO

    ! READ ATTRIBUTES FROM FILE
    ! ------------------------------------------------
    SELECT CASE(i)
    CASE(FTYPE_UNDEFINED)
       status = 3200 ! OUTPUT FILE TYPE UNDEFINED
    CASE(FTYPE_ASCII)
       status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED 
#ifdef PNETCDF
    CASE(FTYPE_NETCDF)
#else
    CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
       !
       lp = .FALSE.
       !
       IF (lp_io) THEN
          CALL ch_netcdf_init_rst(status, &
               TRIM(fname_base)//TRIM(FTYPE_EXT_TEXT(i)), rstatt)
       ELSE
          status = 0
       END IF
       !
#ifdef PNETCDF
    CASE(FTYPE_PNETCDF)
       !
       lp = .TRUE.
       !
       CALL ch_pnetcdf_init_rst(status, &
            TRIM(fname_base)//TRIM(FTYPE_EXT_TEXT(i)), rstatt)
       !
#endif
    CASE(FTYPE_GRIB)
       status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED 
    CASE(FTYPE_HDF4)
       status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
    CASE(FTYPE_HDF5)
       status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
    CASE DEFAULT
       status = 3206 ! RESTART FILE REQUIRED BUT NOT PRESENT 
    END SELECT
    !
    IF (status /= 0) RETURN
    ! ------------------------------------------------

  END SUBROUTINE channel_init_restart
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
! um_ak_20110512+
!  SUBROUTINE channel_init_io(status, lp_io, IOMODE, fname, AMODE, att)
  SUBROUTINE channel_init_io(status, lp_io, IOMODE, fname, AMODE, att, chname)
! um_ak_20110512-

    USE messy_main_channel,         ONLY: new_attribute, AF_RST_CMP &
                                        , AF_RST_INP, get_attribute
    USE messy_main_channel_attributes, ONLY: t_attribute_list
    USE messy_main_channel_netcdf,  ONLY: ch_netcdf_init_io
#ifdef PNETCDF
    USE messy_main_channel_pnetcdf, ONLY: ch_pnetcdf_init_io
#endif

    IMPLICIT NONE
    
    INTRINSIC :: ASSOCIATED, TRIM

    ! I/O
    INTEGER,                INTENT(OUT)       :: status
    LOGICAL,                INTENT(IN)        :: lp_io
    INTEGER,                INTENT(IN)        :: IOMODE
    CHARACTER(LEN=*),       INTENT(IN)        :: fname
    INTEGER,                INTENT(IN)        :: AMODE
    TYPE(t_attribute_list), POINTER, OPTIONAL :: att ! INTENT(IN)
    CHARACTER(LEN=*),       INTENT(IN), OPTIONAL :: chname !um_ak_20110512

    ! LOCAL
    !CHARACTER(LEN=*),      PARAMETER :: substr = 'channel_init_io'
    TYPE(t_channel_list),  POINTER   :: ls
    TYPE(t_channel),       POINTER   :: channel
    LOGICAL                          :: lskip
    LOGICAL                          :: lp

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT
          
       channel => ls%this
          
       ! ------------------------------------------------
       SELECT CASE(AMODE)
       CASE(AMODE_READ)
          ! READ ALL AVAILABLE DATA
          lskip = .FALSE.
          ! um_ak_20110512+
          IF (PRESENT(chname)) THEN
             IF (TRIM(chname) /= TRIM(channel%name)) lskip = .TRUE.
          ENDIF
          ! um_ak_20110512-
       CASE(AMODE_WRITE)
          ! NEW OUTPUT FILE ? RESTART FILE ?
          SELECT CASE(IOMODE)
          CASE(IOMODE_OUT)
             lskip = .NOT. (channel%int%lnew_file .AND. channel%int%lout_now)
          CASE(IOMODE_RST)
             lskip = .NOT. channel%int%lrst
          END SELECT
       END SELECT
       IF (lskip) THEN
          ls => ls%next
          CYCLE
       END IF
       ! ------------------------------------------------

       ! ------------------------------------------------
       IF (channel%io%ftype(IOMODE) == FTYPE_UNDEFINED) THEN
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
          RETURN
       END IF

       IF (channel%io%ftype(IOMODE) > FTYPE_MAXIMUM) THEN
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
          RETURN
       END IF

       channel%int%fname(IOMODE)= &
            TRIM(fname)//TRIM(channel%name)//&
            &FTYPE_EXT_TEXT(channel%io%ftype(IOMODE))
       ! ------------------------------------------------

       ! ------------------------------------------------
       CALL new_attribute(status, channel%att &
            , 'channel_name', c=TRIM(channel%name)&
            , loverwrite=.TRUE., iflag = AF_RST_CMP)
       IF (status /= 0) RETURN

       CALL new_attribute(status, channel%att  &
            , 'channel_file_type', c=TRIM(IOMODE_TEXT(IOMODE)) &
            , loverwrite=.TRUE., iflag = AF_RST_CMP)
       IF (status /= 0) RETURN

       CALL new_attribute(status, channel%att &
            , 'channel_file_name', c=TRIM(channel%int%fname(IOMODE)) &
            , loverwrite=.TRUE.)
       IF (status /= 0) RETURN
  
       IF (IOMODE == IOMODE_RST) THEN
          !
          CALL new_attribute(status, channel%att &
               , 'channel_time_slo', r=channel%int%tslo &
               , loverwrite=.TRUE., iflag = AF_RST_INP)
          IF (status /= 0) RETURN
       END IF

       ! ------------------------------------------------

       ! ------------------------------------------------
       SELECT CASE(channel%io%ftype(IOMODE))
       CASE(FTYPE_UNDEFINED)
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
       CASE(FTYPE_ASCII)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED 
#ifdef PNETCDF
       CASE(FTYPE_NETCDF)
#else
       CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
          !
          lp = .FALSE.
          !
          IF (lp_io) THEN
             CALL ch_netcdf_init_io(status, IOMODE, channel, AMODE, att)
          ELSE 
             status = 0
          END IF
          !
#ifdef PNETCDF
       CASE(FTYPE_PNETCDF)
          !
          lp = .TRUE.
          !
          CALL ch_pnetcdf_init_io(status, IOMODE, channel, AMODE, att)
          !
#endif
       CASE(FTYPE_GRIB)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED 
       CASE(FTYPE_HDF4)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF5)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE DEFAULT
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
       END SELECT
       !
       IF (status /= 0) RETURN
       ! ------------------------------------------------

       ! ------------------------------------------------
       IF (lp_io .OR. lp) THEN
          IF (AMODE == AMODE_READ) THEN
             CALL get_attribute(status, channel%att, 'channel_time_slo' &
                  , r=channel%int%tslo)
             IF (status /= 0) RETURN
          END IF
       END IF
       ! ------------------------------------------------

       ls => ls%next
    END DO channel_loop
    
    status = 0

  END SUBROUTINE channel_init_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_write_header(status, lp_io, IOMODE, DIMID_TIME, att)

    USE messy_main_channel_dimensions, ONLY: get_dimension, t_dimension
    USE messy_main_channel_attributes, ONLY: t_attribute_list
    USE messy_main_channel_netcdf,     ONLY: ch_netcdf_write_header
#ifdef PNETCDF
    USE messy_main_channel_pnetcdf,    ONLY: ch_pnetcdf_write_header
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED
    
    ! I/O
    INTEGER,                INTENT(OUT)       :: status
    LOGICAL,                INTENT(IN)        :: lp_io
    INTEGER,                INTENT(IN)        :: IOMODE
    INTEGER,                INTENT(IN)        :: DIMID_TIME
    TYPE(t_attribute_list), POINTER, OPTIONAL :: att ! INTENT(IN)

    ! LOCAL
    !CHARACTER(LEN=*),    PARAMETER :: substr = 'channel_write_header'
    TYPE(t_dimension),     POINTER :: dim_time => NULL()
    TYPE(t_channel_list),  POINTER :: ls       => NULL()
    TYPE(t_channel),       POINTER :: channel  => NULL()
    LOGICAL                        :: lskip

    ! SET TIME DIMENSION
    CALL get_dimension(status, DIMID_TIME, dim_time)
    IF (status /= 0) RETURN

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT
          
       channel => ls%this
       ! ------------------------------------------------
       ! NEW OUTPUT FILE ? RESTART FILE ?
       SELECT CASE(IOMODE)
       CASE(IOMODE_OUT)
          lskip = .NOT. (channel%int%lnew_file .AND. channel%int%lout_now)
       CASE(IOMODE_RST)
          lskip = .NOT. channel%int%lrst
       END SELECT
       IF (lskip) THEN
          ls => ls%next
          CYCLE
       END IF
       ! ------------------------------------------------

       ! ------------------------------------------------
       SELECT CASE(channel%io%ftype(IOMODE))
       CASE(FTYPE_UNDEFINED)
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
       CASE(FTYPE_ASCII)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED 
#ifdef PNETCDF
       CASE(FTYPE_NETCDF)
#else
       CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
          !
          IF (lp_io) THEN
             ! op_bk_20130905+
             WRITE(*,*) "MESSY: ch_netcdf_write_header called"
             WRITE(*,'(A4, A128)') "CN:",TRIM(channel%name)
             ! op_bk_20130905-
             CALL ch_netcdf_write_header(status, IOMODE, channel &
                  , dim_time, att)
             WRITE(*,*) "MESSY: ch_netcdf_write_header DONE."
          ELSE
             status = 0
          END IF
          !
#ifdef PNETCDF
       CASE(FTYPE_PNETCDF)
          !
          CALL ch_pnetcdf_write_header(status, IOMODE, channel &
               , dim_time, att)
          !
#endif
       CASE(FTYPE_GRIB)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED 
       CASE(FTYPE_HDF4)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF5)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE DEFAULT
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
       END SELECT
       !
       IF (status /= 0) RETURN
       ! ------------------------------------------------

       ls => ls%next
    END DO channel_loop
    
    status = 0

  END SUBROUTINE channel_write_header
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_write_time(status, lp_io, IOMODE, DIMID_TIME)

    USE messy_main_channel_dimensions, ONLY: get_dimension, t_dimension
    USE messy_main_channel_netcdf,     ONLY: ch_netcdf_write_time
#ifdef PNETCDF
    USE messy_main_channel_pnetcdf,    ONLY: ch_pnetcdf_write_time
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED
    
    ! I/O
    INTEGER, INTENT(OUT) :: status
    LOGICAL, INTENT(IN)  :: lp_io
    INTEGER, INTENT(IN)  :: IOMODE
    INTEGER, INTENT(IN)  :: DIMID_TIME

    ! LOCAL
    !CHARACTER(LEN=*),   PARAMETER :: substr = 'channel_write_time'
    TYPE(t_dimension),    POINTER :: dim_time => NULL()
    TYPE(t_channel_list), POINTER :: ls       => NULL()
    TYPE(t_channel),      POINTER :: channel  => NULL()

    IF (IOMODE /= IOMODE_OUT) RETURN

    ! SET TIME DIMENSION
    CALL get_dimension(status, DIMID_TIME, dim_time)
    IF (status /= 0) RETURN

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT
          
       channel => ls%this

       ! ------------------------------------------------
       ! OUTPUT NOW ?
       IF (.NOT. channel%int%lout_now) THEN
          ls => ls%next
          CYCLE
       END IF
       ! ------------------------------------------------

       ! ------------------------------------------------
       SELECT CASE(channel%io%ftype(IOMODE))
       CASE(FTYPE_UNDEFINED)
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
       CASE(FTYPE_ASCII)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef PNETCDF
       CASE(FTYPE_NETCDF)
#else
       CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
          !
          IF (lp_io) THEN
             CALL ch_netcdf_write_time(status, channel, dim_time)
          ELSE
             status = 0
          END IF
          !
#ifdef PNETCDF
       CASE(FTYPE_PNETCDF)
          !
          CALL ch_pnetcdf_write_time(status, channel, dim_time)
          !
#endif
       CASE(FTYPE_GRIB)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF4)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF5)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE DEFAULT
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
       END SELECT
       !
       IF (status /= 0) RETURN
       ! ------------------------------------------------

       ls => ls%next
    END DO channel_loop
    
    status = 0

  END SUBROUTINE channel_write_time
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE channel_write_data(status, lp, lp_io, IOMODE, lexit, ptr, reprid)

    USE messy_main_channel_repr,   ONLY: REPR_UNDEF, repr_reorder
    USE messy_main_channel_netcdf, ONLY: ch_netcdf_write_data
#ifdef PNETCDF
    USE messy_main_channel_pnetcdf, ONLY: ch_pnetcdf_write_data
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                      INTENT(OUT)   :: status
    LOGICAL,                      INTENT(OUT)   :: lp
    LOGICAL,                      INTENT(IN)    :: lp_io
    INTEGER,                      INTENT(IN)    :: IOMODE
    LOGICAL,                      INTENT(OUT)   :: lexit
    REAL(DP), DIMENSION(:,:,:,:), POINTER       :: ptr
    INTEGER,                      INTENT(OUT)   :: reprid

    ! LOCAL
    INTEGER, PARAMETER :: MODE_INITIALIZE   = 0
    INTEGER, PARAMETER :: MODE_NEXT_CHANNEL = 1
    INTEGER, PARAMETER :: MODE_NEXT_OBJECT  = 2
    INTEGER, PARAMETER :: MODE_NEXT_DATA    = 3
    INTEGER, PARAMETER :: MODE_OUTPUT       = 4
    
    !CHARACTER(LEN=*),     PARAMETER     :: substr = 'channel_write_data'
    INTEGER,                       SAVE :: MODE = MODE_INITIALIZE  
    TYPE(t_channel_list), POINTER, SAVE :: ls
    TYPE(t_channel),      POINTER, SAVE :: channel
    TYPE(t_channel_object_list), POINTER, SAVE :: le
    TYPE(t_channel_object),      POINTER, SAVE :: object
    ! OUTPUT DATA TYPE 
    INTEGER,                              SAVE :: jsnd = 0
    ! INDEX IN SECONDARY DATA POINTER
    INTEGER,                              SAVE :: i2nd = 0
    REAL(DP), DIMENSION(:,:,:,:), POINTER      :: zptr => NULL()

    ! INIT
    lexit = .FALSE.

    DO

    SELECT CASE(mode)

    CASE(MODE_INITIALIZE)
       !
       ! INIT
       NULLIFY(ptr)
       !
       ls => GCHANNELLIST
       !
       reprid = REPR_UNDEF
       NULLIFY(le)
       NULLIFY(channel)
       NULLIFY(object)
       jsnd = 0
       i2nd = 0
       !
       MODE = MODE_NEXT_CHANNEL
       !
    CASE(MODE_NEXT_CHANNEL)
       !
       ! LOOK FOR NEXT CHANNEL WITH OUTPUT
       DO
          IF (.NOT. ASSOCIATED(ls)) THEN
             lexit = .TRUE.                  ! NO MORE CHANNEL
             lp    = .FALSE.                 ! ... always broadcast lexit ...
             MODE = MODE_INITIALIZE          ! NEXT OUTPUT
             status = 0
             RETURN
          ELSE
             channel => ls%this
             ! OUTPUT OR RESTART
             IF ( (IOMODE == IOMODE_OUT .AND. channel%int%lout_now) .OR. &
                  (IOMODE == IOMODE_RST .AND. channel%int%lrst) ) THEN
                ! O.K. CHANNEL FOR OUTPUT
                le => channel%list
                MODE = MODE_NEXT_OBJECT
                EXIT
             END IF
          END IF
          ls => ls%next
       END DO
       !
    CASE(MODE_NEXT_OBJECT)
       !
       ! LOOK FOR NEXT OBJECT WITH OUTPUT
       DO
          IF (.NOT. ASSOCIATED(le)) THEN  ! NO MORE OBJECT
             ls => ls%next
             MODE = MODE_NEXT_CHANNEL
             EXIT
          ELSE
             object => le%this
             IF ( ((IOMODE == IOMODE_OUT) .AND. object%int%lout) .OR. &
                  ((IOMODE == IOMODE_RST) .AND. object%int%lrst) ) THEN
                ! O.K. OBJECT FOR OUTPUT
                jsnd = 1
                MODE = MODE_NEXT_DATA
                EXIT
             END IF
          END IF
          le => le%next
       END DO
       !
    CASE(MODE_NEXT_DATA)
       !
       ! LOOK FOR NEXT DATA WITH OUTPUT
       DO
          IF (jsnd > SND_MAXLEN) THEN ! NO MORE DATA AVAILABLE
             le => le%next
             MODE = MODE_NEXT_OBJECT
             EXIT
          ELSE
             !
             SELECT CASE(jsnd)
             CASE(SND_INS)
                !!$ ptr => object%data(:,:,:,:) ! mz_ab_20090921
                ptr => object%ioptr(:,:,:,:) ! mz_ab_20090921
             CASE DEFAULT ! SND_AVE, SND_STP, SND_MIN, SND_MAX, ...
                i2nd = object%int%i2nd(jsnd)
                IF (i2nd > 0) THEN
                   ptr => object%sdat(i2nd)%ptr(:,:,:,:)
                ELSE
                   NULLIFY(ptr)
                END IF
             END SELECT
             !
             IF (object%int%lexp(jsnd, IOMODE)) THEN
                ! OUTPUT OF PRIMARY/SECONDARY DATA REQUESTED
                reprid = object%repr%id
                MODE = MODE_OUTPUT
                !
                SELECT CASE(channel%io%ftype(IOMODE))
                CASE(FTYPE_UNDEFINED)
                   status = 3200 ! OUTPUT FILE TYPE UNDEFINED
                CASE(FTYPE_ASCII)
                   status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef PNETCDF
                CASE(FTYPE_NETCDF)
#else
                CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
                   !
                   lp = .FALSE.
                   status = 0
                   !
#ifdef PNETCDF
                CASE(FTYPE_PNETCDF)
                   !
                   lp = .TRUE.
                   status = 0
                   !
#endif
                CASE(FTYPE_GRIB)
                   status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
                CASE(FTYPE_HDF4)
                   status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
                CASE(FTYPE_HDF5)
                   status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
                CASE DEFAULT
                   status = 3202 ! OUTPUT FILE TYPE UNKNOWN
                END SELECT
                !
                RETURN ! JUMP BACK TO GATHER ON ONE PE
             END IF
          END IF
          jsnd = jsnd + 1
       END DO
       !
    CASE(MODE_OUTPUT)
       !
!!$       IF (ASSOCIATED(ptr)) THEN  ! I/O - PE
!!$       END IF
       !
       SELECT CASE(channel%io%ftype(IOMODE))
       CASE(FTYPE_UNDEFINED)
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
       CASE(FTYPE_ASCII)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef PNETCDF
       CASE(FTYPE_NETCDF)
#else
       CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
          !
          lp = .FALSE.
          !
          IF (lp_io) THEN
! um_ak_20090120+
             ! CHECK, IF OBJECT HAS A "gather-able" SHAPE
             IF (ASSOCIATED(ptr)) THEN
! um_ak_20090120-
                CALL repr_reorder(status, 1, lp, object%repr, ptr, zptr)
                IF (status /= 0) RETURN
                !
                CALL ch_netcdf_write_data(status, IOMODE, channel &
                     , object, zptr, jsnd, i2nd)
! um_ak_20090120+
             ELSE
                status=0
             ENDIF
! um_ak_20090120-
          ELSE
             status = 0
          END IF
          !
#ifdef PNETCDF
       CASE(FTYPE_PNETCDF)
          !
          lp = .TRUE.
          !
! um_ak_20090120+
          ! OBJECT OUTPUT IS NOT IMPLEMENTED YET
          IF (ASSOCIATED(ptr)) THEN
! um_ak_20090120-
             CALL repr_reorder(status, 1, lp, object%repr, ptr, zptr)
             IF (status /= 0) RETURN
             !
             CALL ch_pnetcdf_write_data(status, IOMODE, channel &
                  , object, zptr, jsnd, i2nd)
! um_ak_20090120+
          ELSE
             status=0
          ENDIF
! um_ak_20090120-
             !
#endif
       CASE(FTYPE_GRIB)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF4)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF5)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE DEFAULT
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
       END SELECT
       !
       ! CLEAN UP
       IF (ASSOCIATED(zptr)) THEN
          DEALLOCATE(zptr)
          NULLIFY(zptr)
       END IF
       !
       jsnd = jsnd + 1
       MODE = MODE_NEXT_DATA
       RETURN
       !
    END SELECT

    END DO

  END SUBROUTINE channel_write_data
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
! um_ak_20110512+
!  SUBROUTINE channel_finish_io(status, lp_io, IOMODE, lclose)
  SUBROUTINE channel_finish_io(status, lp_io, IOMODE, lclose, chname)
! um_ak_20110512+

    USE messy_main_channel_netcdf, ONLY: ch_netcdf_finish_io
#ifdef PNETCDF
    USE messy_main_channel_pnetcdf, ONLY: ch_pnetcdf_finish_io
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED
    
    ! I/O
    INTEGER,  INTENT(OUT) :: status
    LOGICAL,  INTENT(IN)  :: lp_io
    INTEGER,  INTENT(IN)  :: IOMODE
    LOGICAL,  INTENT(IN)  :: lclose
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: chname !um_ak_20110512

    ! LOCAL
    !CHARACTER(LEN=*),     PARAMETER :: substr = 'channel_finish_io'
    TYPE(t_channel_list), POINTER   :: ls
    TYPE(t_channel),      POINTER   :: channel
    LOGICAL                         :: lskip
    LOGICAL                         :: lcycle ! um_ak_20110512

    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT
          
       channel => ls%this

       lcycle = .FALSE. !um_ak_20110512
       ! ------------------------------------------------
       SELECT CASE(IOMODE)
       CASE(IOMODE_OUT)
          lskip = .NOT. channel%int%lout_now
       CASE(IOMODE_RST)
          lskip = .NOT. channel%int%lrst   
          !um_ak_20110512+
          IF (PRESENT(chname)) THEN
             IF (TRIM(chname) /= TRIM(channel%name)) lcycle = .TRUE.
          ENDIF
          !um_ak_20110512-
       END SELECT
!      IF (lskip .AND. (.NOT.lclose)) THEN !um_ak_20110512
       IF ((lskip .AND. (.NOT.lclose)) .OR. lcycle) THEN !um_ak_20110512
          ls => ls%next
          CYCLE
       END IF
       ! ------------------------------------------------

       ! RESET TRIGGER FOR NEW FILE
       IF (IOMODE == IOMODE_OUT) channel%int%lnew_file = .FALSE.

       ! ------------------------------------------------
       ! OPEN FILE ? (= NAME ASSIGNED ?)
       !
       SELECT CASE(channel%io%ftype(IOMODE))
       CASE(FTYPE_UNDEFINED)
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
       CASE(FTYPE_ASCII)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef PNETCDF
       CASE(FTYPE_NETCDF)
#else
       CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
          !
          IF (lp_io) THEN
             CALL ch_netcdf_finish_io(status, IOMODE, channel, lclose)
          ELSE
             status = 0
          END IF
          !
#ifdef PNETCDF
       CASE(FTYPE_PNETCDF)
          !
          CALL ch_pnetcdf_finish_io(status, IOMODE, channel, lclose)
          !
#endif
       CASE(FTYPE_GRIB)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF4)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF5)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE DEFAULT
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
       END SELECT
       !
       IF (status /= 0) RETURN
       !
       ! ------------------------------------------------
       
       ls => ls%next
    END DO channel_loop
    
    status = 0

  END SUBROUTINE channel_finish_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
!um_ak_20110512+
!!$  SUBROUTINE channel_read_data(status, lp_io, IOMODE, lexit, ptr, reprid, lp)
  SUBROUTINE channel_read_data(status, lp_io, IOMODE, lexit, ptr, reprid, lp &
       , chname)
!um_ak_20110512-

    USE messy_main_channel_repr,   ONLY: REPR_UNDEF, repr_reorder
    USE messy_main_channel_netcdf, ONLY: ch_netcdf_read_data
#ifdef PNETCDF
    USE messy_main_channel_pnetcdf, ONLY: ch_pnetcdf_read_data
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                      INTENT(OUT)   :: status 
    LOGICAL,                      INTENT(IN)    :: lp_io
    INTEGER,                      INTENT(IN)    :: IOMODE
    LOGICAL,                      INTENT(OUT)   :: lexit
    REAL(DP), DIMENSION(:,:,:,:), POINTER       :: ptr
    INTEGER,                      INTENT(OUT)   :: reprid
    LOGICAL,                      INTENT(OUT)   :: lp
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL      :: chname !um_ak_20110512

    ! LOCAL
    INTEGER, PARAMETER :: MODE_INITIALIZE   = 0
    INTEGER, PARAMETER :: MODE_NEXT_CHANNEL = 1
    INTEGER, PARAMETER :: MODE_NEXT_OBJECT  = 2
    INTEGER, PARAMETER :: MODE_NEXT_DATA    = 3
    INTEGER, PARAMETER :: MODE_INPUT        = 4
    INTEGER, PARAMETER :: MODE_DISTRIBUTE   = 5
    
    !CHARACTER(LEN=*),     PARAMETER     :: substr = 'channel_read_data'
    INTEGER,                       SAVE :: MODE = MODE_INITIALIZE  
    TYPE(t_channel_list), POINTER, SAVE :: ls
    TYPE(t_channel),      POINTER, SAVE :: channel
    TYPE(t_channel_object_list), POINTER, SAVE :: le
    TYPE(t_channel_object),      POINTER, SAVE :: object
    ! OUTPT DATA TYPE 
    INTEGER,                              SAVE :: jsnd = 0
    ! INDEX IN SECONDARY DATA POINTER
    INTEGER,                              SAVE :: i2nd = 0
    REAL(DP), DIMENSION(:,:,:,:), POINTER      :: zptr  

    ! ONLY FOR RESTART FILES
    IF (IOMODE /= IOMODE_RST) THEN
       status = 3205 ! NO INPUT OF OUTPUT FILES
       RETURN
    END IF

    ! INIT
    lexit = .FALSE.

    DO

    SELECT CASE(mode)

    CASE(MODE_INITIALIZE)
       !
       ! INIT
       IF (ASSOCIATED(ptr)) DEALLOCATE(ptr)
       NULLIFY(ptr)
       !
       ls => GCHANNELLIST
       !
       reprid = REPR_UNDEF
       NULLIFY(le)
       NULLIFY(channel)
       NULLIFY(object)
       jsnd = 0
       i2nd = 0
       !
       MODE = MODE_NEXT_CHANNEL
       !
    CASE(MODE_NEXT_CHANNEL)
       !
       ! LOOK FOR NEXT CHANNEL
       IF (.NOT. ASSOCIATED(ls)) THEN
          lexit = .TRUE.                  ! NO MORE CHANNEL
          lp    = .FALSE.                 ! ... always broadcast lexit ...
          MODE = MODE_INITIALIZE          ! NEXT INPUT
          status = 0
          RETURN
       ELSE
          channel => ls%this
          le => channel%list
          MODE = MODE_NEXT_OBJECT             
! um_ak_20110512+
          IF (PRESENT(chname)) THEN
             IF (TRIM(chname) /= TRIM(channel%name) ) THEN
                ls => ls%next
                MODE = MODE_NEXT_CHANNEL
             ENDIF
          ENDIF
! um_ak_20110512-
       END IF
       !
    CASE(MODE_NEXT_OBJECT)
       !
       ! LOOK FOR NEXT OBJECT
       IF (.NOT. ASSOCIATED(le)) THEN  ! NO MORE OBJECT
          ls => ls%next
          MODE = MODE_NEXT_CHANNEL
       ELSE
          object => le%this
          jsnd = 1
          MODE = MODE_NEXT_DATA
       END IF
       !
    CASE(MODE_NEXT_DATA)
       !
       ! LOOK FOR NEXT DATA
       DO
          IF (jsnd > SND_MAXLEN) THEN ! NO MORE DATA AVAILABLE
             le => le%next
             MODE = MODE_NEXT_OBJECT
             EXIT
          ELSE
             reprid = object%repr%id
             MODE = MODE_INPUT
             !
             SELECT CASE(jsnd)
             CASE(SND_INS)
                ! PRIMARY DATA
                EXIT
             CASE DEFAULT ! SND_AVE, SND_STP, SND_MIN, SND_MAX
                ! SEONDARY DATA, IF PRESENT
                i2nd = object%int%i2nd(jsnd)
                IF (i2nd > 0) THEN
                   EXIT
                END IF
             END SELECT
             !
          END IF
          jsnd = jsnd + 1
       END DO
       !
    CASE(MODE_INPUT)
       !
       SELECT CASE(channel%io%ftype(IOMODE))
       CASE(FTYPE_UNDEFINED)
          status = 3200 ! OUTPUT FILE TYPE UNDEFINED
       CASE(FTYPE_ASCII)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
#ifdef PNETCDF
       CASE(FTYPE_NETCDF)
#else
       CASE(FTYPE_NETCDF, FTYPE_PNETCDF)
#endif
          !
          lp = .FALSE.
          !
          IF (lp_io) THEN  ! I/O - PE
             CALL ch_netcdf_read_data(status, IOMODE, channel &
                  , object, zptr, jsnd, i2nd)
          ELSE
             ! NOTHING TO DO FOR NON-I/O PE
             IF (ASSOCIATED(ptr)) DEALLOCATE(ptr)
             NULLIFY(ptr)
             status = 0
          END IF
          !
#ifdef PNETCDF
       CASE(FTYPE_PNETCDF)
          !
          lp = .TRUE.
          !
          CALL ch_pnetcdf_read_data(status, IOMODE, channel &
               , object, zptr, jsnd, i2nd)
          !
#endif
       CASE(FTYPE_GRIB)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF4)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE(FTYPE_HDF5)
          status = 3201 ! OUTPUT FILE TYPE NOT YET IMPLEMENTED
       CASE DEFAULT
          status = 3202 ! OUTPUT FILE TYPE UNKNOWN
       END SELECT
       !
       IF (status /= 0) RETURN
       !
       ! REORDER DATA
       IF (ASSOCIATED(zptr)) THEN
          CALL repr_reorder(status, -1, lp, object%repr, ptr, zptr)
          IF (status /= 0) RETURN
          ! CLEAN UP
          DEALLOCATE(zptr)
          NULLIFY(zptr)
       END IF
       !
       ! RETURN TO SCATTER ON ALL PEs AND DISTRIBUTE AFTERWARDS
       MODE = MODE_DISTRIBUTE
       RETURN
       !
    CASE(MODE_DISTRIBUTE)
       !
       IF (ASSOCIATED(ptr)) THEN
          ! DISTRIBUTE ON ALL PEs AFTER SCATTER
          SELECT CASE(jsnd)
          CASE(SND_INS)
          !!$ object%data(:,:,:,:) = ptr(:,:,:,:) ! mz_ab_20090921
             object%ioptr(:,:,:,:) = ptr(:,:,:,:) ! mz_ab_20090921
          CASE DEFAULT ! SND_AVE, SND_STP, SND_MIN, SND_MAX
             object%sdat(i2nd)%ptr(:,:,:,:) = ptr(:,:,:,:)
          END SELECT
          !
          ! FLAG DATA AS READ
          object%int%lrestart_read = .TRUE.
          !
       END IF
       !
       ! NEXT DATA
       jsnd = jsnd + 1
       MODE = MODE_NEXT_DATA
       RETURN
       !
    END SELECT

    END DO

  END SUBROUTINE channel_read_data
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_main_channel_io
! **********************************************************************
