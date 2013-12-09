! **********************************************************************
MODULE messy_main_channel_tracer
! **********************************************************************

  ! MESSY DATA TRANSFER AND EXPORT INTERFACE (MEMORY MANAGEMENT)
  !
  ! Author: Patrick Joeckel, MPICH, May 2005

  USE messy_main_constants_mem,     ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: create_tracer_channels
  PUBLIC :: set_channel_or_tracer

CONTAINS

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE create_tracer_channels(status, trsetname, channelname, reprid &
                                   ,nnexrst) ! um_ak_20101021

    USE messy_main_tracer,        ONLY: get_tracer_set, t_trinfo_list &
                                      , param2string &
                                      , MAX_CASK_I, NAMES_CASK_I &
                                      , MAX_CASK_R, NAMES_CASK_R &
                                      , MAX_CASK_S, NAMES_CASK_S
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
    USE messy_main_channel_repr,  ONLY: t_representation, get_representation
    USE messy_main_channel,       ONLY: new_channel, new_attribute &
                                      , STRLEN_OBJECT, new_channel_object

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, NULL, TRIM, SIZE
    
    ! I/O
    INTEGER,           INTENT(OUT) :: status
    CHARACTER(LEN=*),  INTENT(IN)  :: trsetname    ! TRACER SET NAME
    CHARACTER(LEN=*),  INTENT(IN)  :: channelname  ! CHANNEL NAME
    INTEGER,           INTENT(IN)  :: reprid       ! REPRESENTATION ID
    ! um_ak_20101021+
    INTEGER, OPTIONAL, INTENT(IN)  :: nnexrst      ! NUMBER OF EXTENDED MEMORY
                                                   ! IN RESTART
    ! um_ak_20101021-

    ! LOCAL
    TYPE(t_trinfo_list),             POINTER :: trlist
    INTEGER                                  :: ntrac
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER :: pxt
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER :: pxtte
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER :: pxtm1
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER :: pxmem
    !
    TYPE(t_representation),          POINTER :: repr
    !
    TYPE(t_trinfo_list), POINTER             :: til => NULL()
    INTEGER                                  :: i, j
    CHARACTER(LEN=STRLEN_OBJECT)             :: oname
    REAL(DP), DIMENSION(:,:,:,:),    POINTER :: mem
    CHARACTER(LEN=STRLEN_MEDIUM)             :: attstr
    INTEGER                                  :: n, nex
    CHARACTER(LEN=4)                         :: estr = ''
    ! um_ak_20101021+
    INTEGER                                  :: nexrst

    IF (PRESENT(nnexrst)) THEN
       nexrst = nnexrst
    ELSE
       nexrst = 0
    ENDIF
    ! um_ak_20101021-

    ! TRACER INFORMATION
    CALL get_tracer_set(status, trsetname, trlist=trlist, ntrac=ntrac &
         , xt=pxt, xtte=pxtte, xtm1=pxtm1, xmem=pxmem)
    IF (status /= 0) RETURN

    ! DO NOT CREATE AN EMPTY CHANNEL
    IF (ntrac == 0) THEN
       status = 0
       RETURN
    END IF

    ! REPRESENTATION
    CALL get_representation(status, reprid, repr)
    IF (status /= 0) RETURN

    ! CREATE NEW CHANNEL
    CALL new_channel(status, channelname, reprid)
    IF (status /= 0) RETURN
    ! ADD GLOBAL ATTRIBUTES
    CALL new_attribute(status, channelname &
         , 'TRACER_SET', c=TRIM(trsetname))
    IF (status /= 0) RETURN
    CALL new_attribute(status, channelname &
         , 'REPRESENTATION', c=TRIM(repr%name))
    IF (status /= 0) RETURN

    ! CREATE OPTIONAL CHANNELS
    IF (ASSOCIATED(pxtte)) THEN
       CALL new_channel(status, TRIM(channelname)//'_te', reprid)
       IF (status /= 0) RETURN
       ! ADD GLOBAL ATTRIBUTES
       CALL new_attribute(status, TRIM(channelname)//'_te' &
            , 'TRACER_SET', c=TRIM(trsetname))
       IF (status /= 0) RETURN
       CALL new_attribute(status, TRIM(channelname)//'_te' &
            , 'REPRESENTATION', c=TRIM(repr%name))
       IF (status /= 0) RETURN
    END IF
    !
    IF (ASSOCIATED(pxtm1)) THEN
       CALL new_channel(status, TRIM(channelname)//'_m1', reprid)
       IF (status /= 0) RETURN
       ! ADD GLOBAL ATTRIBUTES
       CALL new_attribute(status, TRIM(channelname)//'_m1' &
            , 'TRACER_SET', c=TRIM(trsetname))
       IF (status /= 0) RETURN
       CALL new_attribute(status, TRIM(channelname)//'_m1' &
            , 'REPRESENTATION', c=TRIM(repr%name))
       IF (status /= 0) RETURN
    END IF
    !
    IF (ASSOCIATED(pxmem)) THEN
       nex = SIZE(pxmem,5)
       DO n=1, nex
          WRITE(estr,'(a1,(i3.3))') 'x',n 
          CALL new_channel(status, TRIM(channelname)//'_'//estr, reprid)
          IF (status /= 0) RETURN
          ! ADD GLOBAL ATTRIBUTES
          CALL new_attribute(status, TRIM(channelname)//'_'//estr &
               , 'TRACER_SET', c=TRIM(trsetname))
          IF (status /= 0) RETURN
          CALL new_attribute(status, TRIM(channelname)//'_'//estr &
               , 'REPRESENTATION', c=TRIM(repr%name))
          IF (status /= 0) RETURN
       END DO
    END IF
       
    ! LOOP OVER ALL TRACERS AND ADD OBJECTS
    IF (ntrac > 0) THEN

       til => trlist
       tracer_loop: DO
          IF (.NOT. ASSOCIATED(til)) EXIT

          i = til%info%ident%idx
          oname = TRIM(til%info%ident%fullname)

          ! NOTE: TRACER MEMORY HAS BEEN SETUP, SUCH THAT NUMBER OF TRACERS
          !       IS ALWAYS AT POSITION 3, i.e.
          !       (:,:,ntrac,:,1)
          mem => pxt(:,:,i,:,:)
          CALL new_channel_object(status &
               , channelname, TRIM(oname), mem=mem &
               , lrestreq=.TRUE.)
          IF (status /= 0) RETURN

          IF (ASSOCIATED(pxtte)) THEN
             mem => pxtte(:,:,i,:,:)
             CALL new_channel_object(status &
                  , TRIM(channelname)//'_te', TRIM(oname), mem=mem)
             IF (status /= 0) RETURN             
          END IF

          IF (ASSOCIATED(pxtm1)) THEN
             mem => pxtm1(:,:,i,:,:)
             CALL new_channel_object(status &
                  , TRIM(channelname)//'_m1', TRIM(oname), mem=mem &
                  , lrestreq=.TRUE.)
             IF (status /= 0) RETURN             
          END IF
             
          IF (ASSOCIATED(pxmem)) THEN
             nex = SIZE(pxmem,5)
             DO n=1, nex
                WRITE(estr,'(a1,(i3.3))') 'x',n 
                mem => pxmem(:,:,i,:,n:n)
                ! um_ak_20101021+
!!$                CALL new_channel_object(status &
!!$                     , TRIM(channelname)//'_'//estr, TRIM(oname), mem=mem)
                IF (n<=nexrst) THEN
                   CALL new_channel_object(status &
                        , TRIM(channelname)//'_'//estr, TRIM(oname), mem=mem&
                        , lrestreq=.TRUE.)
                ELSE
                   CALL new_channel_object(status &
                        , TRIM(channelname)//'_'//estr, TRIM(oname), mem=mem)
                ENDIF
                ! um_ak_20101021-
                IF (status /= 0) RETURN             
             END DO
          END IF

          ! ADD ALL ATTRIBUTES
          CALL new_attribute(status, channelname, oname, &
               'long_name', c = TRIM(til%info%ident%longname) )

          CALL new_attribute(status, channelname, oname, &
               'units',     c = TRIM(til%info%ident%unit) )

          CALL new_attribute(status, channelname, oname, &
               'submodel',  c = TRIM(til%info%ident%submodel) )

          CALL new_attribute(status, channelname, oname, &
               'index',     i = til%info%ident%idx)

          attstr = param2string(til%info%ident%medium, 'medium')
          CALL new_attribute(status, channelname, oname, &
               'medium',    c = TRIM(attstr) )

          attstr = param2string(til%info%ident%quantity, 'quantity')
          CALL new_attribute(status, channelname, oname, &
               'quantity',  c = TRIM(attstr) )

          attstr = param2string(til%info%ident%type, 'type')
          CALL new_attribute(status, channelname, oname, &
               'type',      c = TRIM(attstr) )
          ! um_ak_20100604+
          ! This attribute normally indicates that a channel object
          ! contains the time-levels for a prognostic variable.
          ! However, this cannot be achieved for TRACERs, as the tendency
          ! is interrupting the continuity of the different time levels.
          ! Hae ? What?
          ! However, MMD coupling and indicator require that 
          ! tracers are "in general" time dependent. To indicate that this
          ! channel object itself does not contain the time levels, the integer
          ! is set to a negative value and not to the real number of 
          ! time levels the indicator needed for channel ele
           CALL new_attribute(status, channelname, oname, &
               'number_of_timelevels', i = -1 )
          ! um_ak_20100604-

          DO j=1, MAX_CASK_I
             attstr = param2string(til%info%meta%cask_i(j), 'switch')
             CALL new_attribute(status, channelname, oname, &
                  TRIM(NAMES_CASK_I(j)), c = TRIM(attstr) )            
          END DO

          DO j=1, MAX_CASK_R
             CALL new_attribute(status, channelname, oname, &
                  TRIM(NAMES_CASK_R(j)), r = til%info%meta%cask_r(j) )
          END DO

          DO j=1, MAX_CASK_S
             CALL new_attribute(status, channelname, oname, &
                  TRIM(NAMES_CASK_S(j)), c = TRIM(til%info%meta%cask_s(j)) )
          END DO

          til => til%next
       END DO tracer_loop
    END IF
   
    status = 0

  END SUBROUTINE create_tracer_channels
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_channel_or_tracer(status, TRSTR, chstr  &
       , cname, oname, pxt, pxtte)

    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
    USE messy_main_tracer,        ONLY: get_tracer
    USE messy_main_channel,       ONLY: get_channel_object

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, TRIM, INDEX, ASSOCIATED

    ! I/O
    INTEGER,                    INTENT(OUT)  :: status
    CHARACTER(LEN=*),           INTENT(IN)   :: TRSTR   ! TRACER SET NAME
    CHARACTER(LEN=*),           INTENT(IN)   :: chstr   ! TRACER CHANNEL NAME
    CHARACTER(LEN=*),           INTENT(IN)   :: cname, oname
    REAL(DP), DIMENSION(:,:,:), POINTER      :: pxt
    REAL(DP), DIMENSION(:,:,:), POINTER      :: pxtte

    ! LOCAL
    CHARACTER(LEN=2*STRLEN_MEDIUM+1) :: tracname = ''
    CHARACTER(LEN=STRLEN_MEDIUM)     :: basename = ''
    CHARACTER(LEN=STRLEN_MEDIUM)     :: subname  = ''
    INTEGER                          :: si

    ! INIT
    IF (ASSOCIATED(pxt)) THEN
       status = 1005
       RETURN
    END IF

    IF (ASSOCIATED(pxtte)) THEN
       status = 1006
       RETURN
    END IF

    NULLIFY(pxt)
    NULLIFY(pxtte)

    IF (TRIM(cname) == chstr) THEN
       ! TRACER
       tracname=ADJUSTL(TRIM(oname))
       si = INDEX(tracname, '_')
       IF (si == 0) THEN
          basename = TRIM(tracname)
          subname  = ''
       ELSE
          basename = TRIM(tracname(:si-1))
          subname  = TRIM(tracname(si+1:))
       END IF
       CALL get_tracer(status, TRSTR, basename, subname=subname &
            , pxt=pxt, pxtte=pxtte)
       IF (status /= 0) RETURN
    ELSE
       ! CHANNEL OBJECT
       CALL get_channel_object(status &
            , TRIM(ADJUSTL(cname)), TRIM(ADJUSTL(oname)) &
            , p3 = pxt )
       IF (status /= 0) RETURN
    END IF
    
    status = 0

  END SUBROUTINE set_channel_or_tracer
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_main_channel_tracer
! **********************************************************************
