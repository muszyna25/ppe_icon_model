! **********************************************************************+
MODULE messy_main_tracer_bi
! **********************************************************************+

  ! BM/MESSy
  USE messy_main_tracer_mem_bi
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  ! MESSy
  USE messy_main_tracer

  IMPLICIT NONE
  PRIVATE

  ! GLOBAL NAMELIST SWITCHES (CPL)
  LOGICAL, SAVE :: l_conv_lg2gp = .TRUE.
  LOGICAL, SAVE :: l_tracer_init = .TRUE.

! op_pj_20120120+ moved to messy_main_tracer_mem_bi.f90
!!$  ! GLOBAL SWITCHES FOR TRACER SETS
!!$  LOGICAL, SAVE :: L_GP   = .FALSE. ! GLOBAL SWITCH   ! mz_ab_20100226
!!$  LOGICAL, SAVE :: L_LG   = .FALSE. ! GLOBAL SWITCH
!!$  LOGICAL, SAVE :: L_OM   = .FALSE. ! GLOBAL SWITCH   ! mz_ap_20071023
!!$  LOGICAL, SAVE :: L_CMAT = .FALSE. ! GLOBAL SWITCH   ! um_ak_20100128
!!$  LOGICAL, SAVE :: L_EC   = .FALSE. ! GLOBAL SWITCH   ! um_ak_20100128
! op_pj_20120120-

  ! PUBLIC SUBROUTINES CALLED FROM BM(I)L
  ! ECHAM5/MESSy SPECIFIC
  PUBLIC :: main_tracer_initialize    ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_new_tracer    ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_init_memory   ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_init_coupling ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_init_tracer   ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_global_start  ! CALLED FROM messy_main_control_e5/c4
  ! SPECIAL: CALLED BEFORE OUTPUT
  PUBLIC :: main_tracer_write_output  ! CALLED FROM messy_main_control_e5/c4
  !
  ! SPECIAL: CALLED FROM BML, AFTER ADVECTION
  ! um_ak_20081030+
  PUBLIC :: main_tracer_beforeadv     ! CALLED DIRECTLY FROM BML
  ! um_ak_20081030-
  PUBLIC :: main_tracer_afteradv      ! CALLED DIRECTLY FROM BML
  !
  PUBLIC :: main_tracer_local_start   ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_global_end    ! CALLED FROM messy_main_control_e5/c4
  PUBLIC :: main_tracer_free_memory   ! CALLED FROM messy_main_control_e5/c4
  !
  !PRIVATE :: setup_tracer_set_gp
  !PRIVATE :: setup_tracer_set_lg
  !PRIVATE :: setup_tracer_set_om ! mz_ap_20071023
  !PRIVATE :: setup_tracer_set_cmat ! mz_ab_20090624
  !PRIVATE :: setup_tracer_set_ec   ! mz_ab_20090624
  !
  !PRIVATE :: main_tracer_read_nml_cpl

  ! SUBMODEL SMIL
  !
  PUBLIC :: tracer_init               ! CALLED FROM XXX_init_tracer
  PUBLIC :: tracer_halt
  !
  ! CONVERSION ROUTINES (FAMILIES) TO BE USED IN SUBMODEL SMIL
  PUBLIC :: main_tracer_fconv_loc     ! CONVERT FAMILIES <-> TRACERS
  PUBLIC :: main_tracer_fconv_glb     ! CONVERT FAMILIES <-> TRACERS
  !


CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_initialize

    ! TRACER MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, July 2002

    ! BM/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_initialize
    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_initialize
    USE messy_main_tools,            ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_tracer_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi( &
            'main_tracer_read_nml_ctrl reported an error', substr)
    END IF
    CALL p_bcast(l_family, p_io)
    CALL p_bcast(l_pdef, p_io)

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_tracer_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi( &
            'main_tracer_read_nml_cpl reported an error', substr)
    END IF
    CALL p_bcast(l_conv_lg2gp, p_io)
    CALL p_bcast(l_tracer_init, p_io)

    IF (l_family) CALL main_tracer_family_initialize

    IF (l_pdef) CALL main_tracer_pdef_initialize

  END SUBROUTINE main_tracer_initialize
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_new_tracer(flag)

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_new_tracer
    USE messy_main_channel_bi,       ONLY: REPR_UNDEF, GP_3D_MID &
                                         , GP_3D_MPIOM           &
                                         , REPRID_cmat3D, REPRID_ec3D

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_new_tracer'
    INTEGER                     :: status

    L_GP = (GP_3D_MID  /= REPR_UNDEF) ! mz_ab_20100226

!!$#if defined(ECHAM5) || defined(MBM_CMAT)
    ! INIT SWITCH
    L_LG = (NGCELL > 0)

    L_OM = (GP_3D_MPIOM /= REPR_UNDEF) ! mz_ap_20071023

    L_CMAT = (REPRID_cmat3D  /= REPR_UNDEF) ! um_ak_20100128
    L_EC   = (REPRID_ec3D    /= REPR_UNDEF) ! um_ak_20100128
!!$#endif

    SELECT CASE(flag)
    CASE(1)
       !
       CALL start_message_bi(modstr,'SETUP TRACER SETS',substr)
       !
       CALL new_tracer_set(status, GPTRSTR, L_GP)
       CALL tracer_halt(substr, status)
       !
       CALL new_tracer_set(status, LGTRSTR, L_LG)
       CALL tracer_halt(substr, status)
       !
       ! mz_ap_20071023+
       CALL new_tracer_set(status, OMTRSTR, L_OM)
       CALL tracer_halt(substr, status)
       ! mz_ap_20071023-
       !
       ! mz_ab_20090624+
       CALL new_tracer_set(status, CMATTRSTR, L_CMAT) ! um_ak_20100129
       CALL tracer_halt(substr, status)  
       !
       CALL new_tracer_set(status, ECTRSTR, L_EC) ! um_ak_20100129
       CALL tracer_halt(substr, status)       
       ! mz_ab_20090624-
       !
       CALL end_message_bi(modstr,'SETUP TRACER SETS',substr)
       !
    CASE (2)
       !
       IF (l_family) CALL main_tracer_family_new_tracer
       !
    CASE(3)
       !
       CALL start_message_bi(modstr,'SHOW TRACER SETS',substr)
       !
       IF (p_parallel_io) CALL print_tracer_set
       !
       CALL end_message_bi(modstr,'SHOW TRACER SETS',substr)
       !
    CASE DEFAULT
       !
       CALL error_bi( 'UNKNOWN FLAG !', substr)
       !
    END SELECT

  END SUBROUTINE main_tracer_new_tracer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_init_memory(flag)

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: error_bi
    ! MESSy
    USE messy_main_tools,            ONLY: remap_bounds ! mz_ab_20100610
    ! mz_ap_20071023: GP_3D_MPIOM added
    USE messy_main_channel_bi,       ONLY: channel_halt, GP_3D_MID &
                                         , LG_ATTILA, GP_3D_MPIOM
    ! mz_ab_20100118+ CMAT2 specific
    USE messy_main_channel_bi,       ONLY: REPRID_cmat3D     &
                                         , REPRID_cmat3D_BND & ! mz_ab_20100509
                                         , REPRID_ec3D
    ! mz_ab_20100118-
    USE messy_main_channel_tracer,   ONLY: create_tracer_channels
    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_init_mem
    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_init_mem
#ifdef COSMO
    USE messy_main_data_bi,          ONLY: l2tls
#endif


    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_init_memory'
    INTEGER :: status
    INTEGER :: iz     ! um_ak_20100729

    SELECT CASE(flag)
    CASE (1)
       !
       CALL start_message_bi(modstr,'SETUP TRACER MEMORY',substr)

#ifndef MBM_CMAT
       CALL setup_tracer_set_gp

       ! INITALIZE POINTER (messy_main_tracer_mem_bi.90)
       IF (L_GP) THEN
          xt    => pxt(:,:,:,:,1)
          xtm1  => pxtm1(:,:,:,:,1)
          xtte  => pxtte(:,:,:,:,1)
          xtf   => pxtf(:,:,:,:,1)
       ENDIF
#endif
#ifdef COSMO
       ! correct pointer from above and allocate pointer array xt_array
       IF (l2tls) THEN 
          NULLIFY(xtf)
          ALLOCATE(xt_array(2))
       ELSE
          xtf   => pxtf(:,:,:,:,3)
          ALLOCATE(xt_array(3))
       ENDIF
       DO iz = 1, SIZE(xt_array)
          NULLIFY(xt_array(iz)%PTR)
       END DO

       ! the boundary data is defined as place 1 and 2 of 
       ! extended tracer memory
       xt_bd => pxtf(:,:,:,:,1:2)
       !xt_bd(:,:,:,:,:) = 0._dp 
       xt_bd(:,:,:,:,:) = -999._dp 
#endif       
       ! ------------

#if defined(ECHAM5) || defined(MBM_CMAT)
! mz_ab_20090624+
       ! CMAT-DOMAIN-ONLY TRACER SET
       CALL setup_tracer_set_cmat

       ! INITALIZE POINTER (messy_main_tracer_mem_bi.90)

       ! ------------
#ifndef MBM_CMAT
       ! ECHAM+CMAT-DOMAIN TRACER SET
       CALL setup_tracer_set_ec

       ! INITALIZE POINTER (messy_main_tracer_mem_bi.90)
       IF (L_EC) THEN
       xt_ec    => pxt_ec(:,:,:,:,1)
!!$       xtm1_ec  => pxtm1_ec(:,:,:,:,1)
!!$       xtte_ec  => pxtte_ec(:,:,:,:,1)
!!$       xtf_ec   => pxtf_ec(:,:,:,:,1)    
!!$       ! INITALIZE POINTER TO SUBDOMAIN
!!$       xt_ecs    => pxt_ec(:,nlev-overlap:,:,:,1)
!!$       xtm1_ecs  => pxtm1_ec(:,:,:,:,1)
!!$       xtte_ecs  => pxtte_ec(:,:,:,:,1)
!!$       xtf_ecs   => pxtf_ec(:,:,:,:,1)
       ENDIF

       ! ------------
! mz_ab_20090624-

       IF ( l_conv_lg2gp .AND. (.NOT. L_LG ) ) THEN
          IF (p_parallel_io) THEN
             WRITE(*,*) ' L_CONV_LG2GP = T in namelist will be ignored,'
             WRITE(*,*) ' since no LAGRANGIAN scheme is running!'
          END IF
          l_conv_lg2gp = .FALSE.
       END IF
       
       CALL setup_tracer_set_lg
       
       ! INITALIZE POINTER (messy_main_tracer_mem_bi.90)
       IF (L_LG) THEN
          xt_a => pxt_a(:,:,:,:,1)
          xtm1_a => pxtm1_a(:,:,:,:,1)
          xtte_a => pxtte_a(:,:,:,:,1)
          xtf_a  => pxtf_a(:,:,:,:,1)
          !
          IF (l_conv_lg2gp) xt_lggp => pxt_lggp(:,:,:,:,1)
       END IF
       !
       ! mz_ap_20071023+
       CALL setup_tracer_set_om
       ! INITALIZE POINTER (messy_main_tracer_mem_bi.90)
       IF (L_OM) THEN
          xt_om => pxt_om(:,:,:,:,1)
       END IF
       ! mz_ap_20071023-
#endif
#endif
       !
       CALL end_message_bi(modstr,'SETUP TRACER MEMORY',substr)
       !
    CASE (2)
       !
       CALL start_message_bi(modstr,'TRACER MEMORY CHANNEL COUPLING',substr)
       !
       ! COUPLE TRACER MEMORY TO CHANNELS
#ifndef COSMO
       CALL create_tracer_channels(status, GPTRSTR, gp_channel, GP_3D_MID)
       CALL channel_halt(substr, status)
#else
       CALL create_tracer_channels(status, GPTRSTR, gp_channel, GP_3D_MID, 2)
       CALL channel_halt(substr, status)
#endif
       !
#if defined(ECHAM5)
       CALL create_tracer_channels(status, LGTRSTR, lg_channel, LG_ATTILA)
       CALL channel_halt(substr, status)
       !
       IF (l_conv_lg2gp) THEN
          CALL create_tracer_channels(status, LGGPTRSTR, &
               lggp_channel, GP_3D_MID)
          CALL channel_halt(substr, status)
       END IF
       !
       ! mz_ap_20071023+
       CALL create_tracer_channels(status, OMTRSTR, om_channel, GP_3D_MPIOM)
       CALL channel_halt(substr, status)
       ! mz_ap_20071023-

       ! mz_ab_20090625+ 
       CALL create_tracer_channels(status, ECTRSTR,   ec_channel,   &
            REPRID_ec3D)
       CALL channel_halt(substr, status)
       ! mz_ab_20090625-
#endif
       !
#if defined(ECHAM5) || defined(MBM_CMAT)
       ! mz_ab_20090625+ 
       CALL create_tracer_channels(status, CMATTRSTR, cmat_channel, &
!            REPRID_cmat3D) 
            REPRID_cmat3D_BND) ! mz_ab_20100509)
       CALL channel_halt(substr, status)       
       ! mz_ab_20090625-
#endif
       !
       !
       IF (l_pdef) CALL main_tracer_pdef_init_mem
       !
       ! setting meta information of family-members to fraction
       ! (for advection initialization)
       IF (l_family) CALL main_tracer_family_init_mem
       !
       CALL end_message_bi(modstr,'TRACER MEMORY CHANNEL COUPLING',substr)
       !
    CASE DEFAULT
       !
       CALL error_bi('UNKNOWN FLAG !', substr)
       !
    END SELECT

  END SUBROUTINE main_tracer_init_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_init_coupling

    ! MESSy
    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_init_cpl
    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_init_cpl

    IMPLICIT NONE

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_init_coupling'

    IF (l_family) CALL main_tracer_family_init_cpl

    IF (l_pdef) CALL main_tracer_pdef_init_cpl

  END SUBROUTINE main_tracer_init_coupling
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_init_tracer(flag)

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: error_bi, info_bi
    USE messy_main_data_bi,          ONLY: eps
#ifndef MESSYTIMER
    USE messy_main_data_bi,          ONLY: lresume
#else
    USE messy_main_timer,            ONLY: lresume
#endif
    USE messy_main_channel_bi,       ONLY: channel_halt
    ! MESSy
    USE messy_main_tracer,           ONLY: NSETID, STRLEN_TRSET &
                                         , get_tracer_set       &
                                         , I_MMD_INIT, ON
    USE messy_main_channel,          ONLY: get_channel_object_info

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, NULL, TRIM

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'main_tracer_init_tracer'
    INTEGER                      :: status
    TYPE(t_trinfo_list), POINTER :: til => NULL()
    CHARACTER(LEN=STRLEN_TRSET)  :: setname = ''
    REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: zpxt   => NULL()
    REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: zpxtm1 => NULL()
    ! NOTE: The order of names in chanle_str MUST be the same as
    !       in the call sequence of new_tracer_set (see CASE(1) in
    !       main_tracer_new_tracer !!!
    ! mz_ap_20071023: om_channel added
    ! mz_ab_20090625: cmat_channel and ec_channel added
    CHARACTER(LEN=*), PARAMETER, DIMENSION(5) :: channel_str = &
         (/ gp_channel, lg_channel, om_channel, cmat_channel, ec_channel/)
    LOGICAL :: linit, linit_m1
    LOGICAL :: l_init
    LOGICAL :: linit_mmd  ! um_ak_20090615
    INTEGER :: i, ntrac

    SELECT CASE(flag)
    CASE(1)
       !
       CALL start_message_bi(modstr,'CHECK TRACER INIT FROM RESTART',substr)
       !
! op_pj_20091129+
    CASE(4)
       !
       CALL start_message_bi(modstr,'INITIALISE TRACER VIA TRACER_INIT',substr)
       IF (l_tracer_init) THEN
          CALL tracer_init
       ELSE
          CALL info_bi(' ... skipped (l_tracer_init=F in CPL)', substr)
       END IF
       CALL end_message_bi(modstr,'INITIALISE TRACER VIA TRACER_INIT',substr)
       RETURN
       !
! op_pj_20091129-
    CASE(2)
       CALL start_message_bi(modstr,'CHECK TRACER INIT BY TRACER_INIT',substr)
       !
    CASE(3)
       !
       CALL start_message_bi(modstr,'DIAGNOSE TRACER INITIALIZATION',substr)
       !
    CASE DEFAULT
       !
       CALL error_bi('UNKNOWN FLAG !', substr)
       !
    END SELECT

    set_loop: DO i=1, NSETID
       !
       CALL get_tracer_set(status, i, setname=setname &
            , trlist=til, xt=zpxt, xtm1=zpxtm1 &
            , ntrac=ntrac, l_init=l_init)
       CALL tracer_halt(substr, status)
       ! mz_ab_20100124+
       IF (p_parallel_io) THEN
          WRITE(*,*) '*** TRACER SET '//TRIM(setname)//' ***'
       END IF       
       ! mz_ab_20100124-

       IF (.NOT. l_init) CYCLE
       IF (ntrac <= 0) CYCLE

! mz_ab_20100124+
! moved above
!!$       IF (flag < 3) THEN
!!$          IF (p_parallel_io) THEN
!!$             WRITE(*,*) '*** TRACER SET '//TRIM(setname)//' ***'
!!$          END IF
!!$       END IF
! mz_ab_20100124-
       !
       !
       tracer_loop: DO
          IF (.NOT.ASSOCIATED(til)) EXIT

          SELECT CASE(flag)
          CASE(1)
             !
             ! CHECK IF TRACER HAS BEEN INITIALIZED VIA RESTART
             ! -> CHANNEL OBJECTS FOR X and X_m1 have restart_read = .true.

             CALL get_channel_object_info(status, channel_str(i)        &
                  , TRIM(til%info%ident%fullname), lrestart_read=linit)
             CALL channel_halt(substr, status)
             ! mz_ap_20071023+
             ! bug fix for tracer sets with only 1 time level
             IF (ASSOCIATED(zpxtm1)) THEN
                CALL get_channel_object_info(status, channel_str(i)//'_m1' &
                     , TRIM(til%info%ident%fullname), lrestart_read=linit_m1)
                CALL channel_halt(substr, status)
             ELSE
                linit_m1 = .TRUE.
             END IF
             ! mz_ap_20071023-
             ! um_ak_20090615+
             ! CHECK IF gridpoint tracer was initialised by MMD
             linit_mmd = til%info%meta%cask_i(I_MMD_INIT) == ON
             til%info%meta%linit = (linit .AND. linit_m1) .OR. linit_mmd 
             ! um_ak_20090615-
             !
             IF (til%info%meta%linit) THEN
                IF (p_parallel_io) THEN
#ifndef I2CINC
                   WRITE(*,*) 'TRACER '//TRIM(til%info%ident%fullname)//' (', &
                        til%info%ident%idx,                  &
                        ') WAS INITIALIZED FROM RESTART-FILE!'
#else
                   IF (linit_mmd .AND. (.NOT. linit) ) THEN
                   WRITE(*,*) 'TRACER '//TRIM(til%info%ident%fullname)//' (', &
                        til%info%ident%idx,                  &
                        ') WAS INITIALIZED BY MMD!'
!                        ') WAS INITIALIZED FROM RESTART-FILE!'
                   ELSE
                   WRITE(*,*) 'TRACER '//TRIM(til%info%ident%fullname)//' (', &
                        til%info%ident%idx,                  &
                        ') WAS INITIALIZED FROM RESTART-FILE!'
!                        ') WAS INITIALIZED BY MMD!'
                   ENDIF
#endif
                END IF
                ! CHECK FOR lforce_init (FORCE TRACER_INIT OR vini)
                IF (til%info%meta%cask_i(I_force_init)==ON) THEN
                   IF (p_parallel_io) THEN
                      WRITE(*,*) ' force_init -> ', &
                           'RESTART INITIALIZATION WILL BE IGNORED'
                   END IF
                   til%info%meta%linit = .false.
                END IF
             ELSE
                IF (p_parallel_io) THEN
                   WRITE(*,*) 'TRACER '//TRIM(til%info%ident%fullname)//' (', &
                        til%info%ident%idx,                  &
                        ') WAS NOT INITIALIZED FROM RESTART-FILE ...'
                END IF
             END IF
             !
          CASE(2)
             !
             ! CHECK IF TRACER HAS BEEN INITIALIZED IN SOME WAY
             IF (til%info%meta%linit) THEN
                IF (p_parallel_io) THEN
                   WRITE(*,*) 'TRACER '//TRIM(til%info%ident%fullname)//' (', &
                        til%info%ident%idx,                  &
                        ') ALREADY INITIALIZED ... !'
                END IF
             ELSE
                IF (p_parallel_io) THEN
                   WRITE(*,*) 'TRACER '//TRIM(til%info%ident%fullname)//' (', &
                        til%info%ident%idx,                  &
                        ') WILL BE SET TO ', til%info%meta%cask_R(R_vini)
                END IF
                zpxt(:,:,til%info%ident%idx,:,:) = til%info%meta%cask_R(R_vini)
                ! mz_ap_20071023+
                ! bug fix for tracer sets with only 1 time level
                IF (ASSOCIATED(zpxtm1)) THEN
                   IF (lresume) &
                        zpxtm1(:,:,til%info%ident%idx,:,:) = (1._DP -eps) &
                        * zpxt(:,:,til%info%ident%idx,:,:)
                END IF
                ! mz_ap_20071023-
                til%info%meta%linit = .TRUE.
             END IF
             !
          CASE (3)
             !
             ! NOTHING TO DO
          END SELECT

          til => til%next
       END DO tracer_loop
       !
    END DO set_loop


    SELECT CASE(flag)
    CASE(1)
       !
       CALL end_message_bi(modstr,'CHECK TRACER INIT FROM RESTART',substr)
       !
    CASE(2)
       !
       CALL end_message_bi(modstr,'CHECK TRACER INIT BY TRACER_INIT',substr)
       !
    CASE(3)
       !
       ! --- DIAGNOSTIC OUTPUT
       !
       IF (p_parallel_io) CALL print_tracer_set_val
       !
       CALL end_message_bi(modstr,'DIAGNOSE TRACER INITIALIZATION',substr)
       !
    CASE DEFAULT
       !
       CALL error_bi('UNKNOWN FLAG !', substr)
       !
    END SELECT

  END SUBROUTINE main_tracer_init_tracer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_global_start

    USE messy_main_tracer_pdef_bi,     ONLY: main_tracer_pdef_global_start
    ! um_ak_20080709+
!!$     USE messy_main_tracer_family_bi,   ONLY: main_tracer_family_global_start
    ! um_ak_20080709-

    IMPLICIT NONE

    !  I/O
!!$    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_global_start'

!!$    SELECT CASE(flag)
!!$       !
!!$    CASE(1)
       !

#if defined(ECHAM5)
    ! LAGRANGIAN TRACERS
    IF (L_LG) THEN
       qxt_a    => xt_a(:,1,:,1)
       qxtte_a  => xtte_a(:,1,:,1)
       qxtm1_a  => xtm1_a(:,1,:,1)
       qxtf_a   => xtf_a(:,1,:,1)
    ENDIF
       !
#endif
       ! 
       IF (l_pdef) CALL main_tracer_pdef_global_start
       !
!!$       ! um_ak_20080709+
!!$       ! moved to main_tracer_beforeadv
!!$    CASE(2)
!!$       !
!!$       ! TYPE-1: t2f
!!$       ! TYPE-2: summation (GPTRSTR)
!!$       IF (l_family) CALL main_tracer_family_global_start
!!$       !
!!$       ! um_ak_20080709-
!!$    END SELECT

  END SUBROUTINE main_tracer_global_start
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! um_ak_20080709+
  SUBROUTINE main_tracer_beforeadv

    USE messy_main_tracer_family_bi,   ONLY: main_tracer_family_beforeadv
    
    IMPLICIT NONE

       ! TYPE-1: t2f
       ! TYPE-2: summation (GPTRSTR)
       IF (l_family) CALL main_tracer_family_beforeadv

  END SUBROUTINE main_tracer_beforeadv
  ! um_ak_20080709-
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_afteradv

    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_afteradv

    IMPLICIT NONE

    IF (l_family) CALL main_tracer_family_afteradv

  END SUBROUTINE main_tracer_afteradv
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_local_start

    ! ECHAM5/MESSy
    USE messy_main_data_bi,          ONLY: jrow

    IMPLICIT NONE

    ! SET POINTERS FOR WITHIN LOCAL LOOP

#if defined(ECHAM5)
    ! GRIDPOINT TRACERS
    IF (L_GP) THEN
       qxt    => xt(:,:,:,jrow)
       qxtte  => xtte(:,:,:,jrow)
       qxtm1  => xtm1(:,:,:,jrow)
       qxtf   => xtf(:,:,:,jrow)
    END IF

    IF (L_LG) THEN
       IF (l_conv_lg2gp) qxt_lggp => xt_lggp(:,:,:,jrow)
    END IF

!!$    ! mz_ab_20090624+
!!$    ! Probably not needed because CMAT runs only outside local loop
!!$    ! CMAT-DOMAIN-ONLY TRACERS
!!$    IF (L_CMAT) THEN
!!$    qxt_cmat    => xt_cmat(:,:,:,jrow)
!!$    qxtte_cmat  => xtte_cmat(:,:,:,jrow)
!!$    qxtm1_cmat  => xtm1_cmat(:,:,:,jrow)
!!$    qxtf_cmat   => xtf_cmat(:,:,:,jrow)   
!!$    END IF
!!$    IF (L_EC) THEN
!!$    ! ECHAM+CMAT-DOMAIN TRACERS
!!$    qxt_ec    => xt_ec(:,:,:,jrow)
!!$    qxtte_ec  => xtte_ec(:,:,:,jrow)
!!$    qxtm1_ec  => xtm1_ec(:,:,:,jrow)
!!$    qxtf_ec   => xtf_ec(:,:,:,jrow)   
!!$    ! INITALIZE POINTER TO CMAT SUBDOMAIN
!!$    qxt_ecs   => xt_ec(:,nlev-overlap:ht_dim,:,jrow)
!!$    qxtte_ecs => xtte_ec(:,:,:,jrow)
!!$    qxtm1_ecs => xtm1_ec(:,:,:,jrow)
!!$    qxtf_ecs  => xtf_ec(:,:,:,jrow)   
!!$    END IF
!!$    ! mz_ab_20090624-
#endif

#if defined(COSMO) || defined(BLANK)
    ! GRIDPOINT TRACERS
    IF (L_GP) THEN
       qxt    => xt(:,jrow,:,:)
       qxtte  => xtte(:,jrow,:,:)
       qxtm1  => xtm1(:,jrow,:,:)
       qxtf   => xtf(:,jrow,:,:)
    ENDIF
#endif

  END SUBROUTINE main_tracer_local_start
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_global_end
    
    ! BMIL
    ! op_bk_20130820+
#ifndef __ICON__
    ! op_bk_20130820-
#if ! (defined(BLANK) || defined(MBM_CMAT))
    USE messy_main_data_bi,          ONLY: grmass=>grmassdry
#endif
    ! op_bk_20130820+
#endif
    ! op_bk_20130820-
    ! BMIL
    USE messy_main_blather_bi,       ONLY: warning_bi
    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_global_end
    ! SMCL
    USE messy_main_tracer_pdef,      ONLY: tracpdef_airmass
#if defined(ECHAM5)
#endif

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_global_end'
    IF (l_pdef) THEN
    ! op_bk_20130820+
#ifndef __ICON__
    ! op_bk_20130820-
#if ! ( defined(BLANK) || defined(MBM_CMAT) )
       !
       CALL tracpdef_airmass(GPTRSTR, grmass)
       !
#if defined(ECHAM5)
       !
#endif
       !
       CALL main_tracer_pdef_global_end
#else
#ifdef BLANK
       CALL warning_bi(substr, 'pdef for BLANK not possible')
#endif
#ifdef MBM_CMAT
       CALL warning_bi(substr, 'pdef for CMAT not possible')
#endif
       CALL warning_bi(substr, 'define grmass first')
#endif
    ! op_bk_20130820+
#else
       CALL warning_bi(substr, 'pdef for ICON not possible')
#endif
    ! op_bk_20130820-
    END IF

  END SUBROUTINE main_tracer_global_end
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_free_memory

    USE messy_main_tracer_pdef_bi,   ONLY: main_tracer_pdef_free_mem
    USE messy_main_tracer_family_bi, ONLY: main_tracer_family_free_mem

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_free_memory'
    INTEGER :: status

    CALL start_message_bi(modstr,'FREE TRACER MEMORY',substr)

    CALL clean_tracer_set(status, GPTRSTR)
    CALL tracer_halt(substr, status)

#if defined(ECHAM5)
    CALL clean_tracer_set(status, LGTRSTR)
    CALL tracer_halt(substr, status)
    !
    IF (l_conv_lg2gp) THEN
       CALL clean_tracer_set(status, LGGPTRSTR)
       CALL tracer_halt(substr, status)
    END IF
#endif

    IF (l_pdef) CALL main_tracer_pdef_free_mem

    IF (l_family) CALL main_tracer_family_free_mem

    CALL end_message_bi(modstr,'FREE TRACER MEMORY',substr)

  END SUBROUTINE main_tracer_free_memory
  ! -------------------------------------------------------------------

  ! ------------------------------------------------------------------
#if defined(ECHAM5)
  SUBROUTINE TRACER_INIT !!$(modstr)

    ! INITIALIZES TRACER FIELDS FROM netCDF FILES
!!$    ! (USING NCREGRID) AS DEFINED IN NAMELIST modstr'_t.nml'
    ! (USING NCREGRID) AS DEFINED IN NAMELIST FILE tracer.nml'
    !
    ! Author: Patrick Joeckel, MPICH, Mainz, November 2002
    !         Patrick Joeckel, MPICH, Mainz, March    2004
    !         Patrick Joeckel, DLR-IPA, Oberpfaffenhofen, December 2009

    ! ECHAM5/MESSy
    ! ... MISC
    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast  &
                                      , dcg, dcl, scatter_gp
    ! mz_kk_20080218+
#ifdef PNCREGRID
    USE messy_main_mpi_bi,        ONLY: p_pe, p_nprocs, p_all_comm,   &
                                        p_send, p_recv
#endif
    ! mz_kk_20080218-
    USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
    USE messy_main_tools,         ONLY: find_next_free_unit
    ! ... GPTRSTR
    USE messy_main_data_bi,       ONLY: eps, nlev                  &
                                      , nlon, ngl, nproma, ngpblks
#ifndef MESSYTIMER
    USE messy_main_data_bi,       ONLY: lresume
#else
    USE messy_main_timer,         ONLY: lresume
#endif
    ! ... LGTRSTR
    ! ... TRACER SETS
    USE messy_main_tracer,        ONLY: NSETID, STRLEN_TRSET, get_tracer_set &
                                      , t_trinfo_list
    ! MESSy
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
    ! ... TRACER
    USE messy_main_tracer,        ONLY: TR_NEXIST
    ! ... NCREGRID
#ifndef PNCREGRID
    USE messy_ncregrid_control,   ONLY: RG_CTRL, RG_NML, RG_STATUS &
                                      , RG_PROC, RG_STOP, NML_NEXT &
                                      , RGSTAT_STOP, REGRID_CONTROL
#else
    USE messy_ncregrid_control,   ONLY: RG_CTRL, RG_NML             &
                                      , RG_PROC, RG_STOP, NML_NEXT  &
                                      , RGSTAT_STOP, REGRID_CONTROL &
                                      ! mz_kk_20080218+
                                      , RGSTAT_NULL, RGSTAT_CYCLE   &
                                      , t_mpi_def, INIT_PNCREGRID
                                      ! mz_kk_20080218-
#endif
    USE messy_ncregrid_base,      ONLY: RGMSG, ERRMSG, RGMLVL      &
                                      , RGMLW, RGMLWC
    USE messy_ncregrid_netcdf,    ONLY: GRD_MAXSTRLEN, ncvar, init_ncvar
    USE messy_ncregrid_geohyb,    ONLY: geohybgrid, init_geohybgrid
    USE messy_ncregrid_tools,     ONLY: RGTOOL_CONVERT
    !USE messy_ncregrid_diag,      ONLY: write_geohybgrid, write_ncvar

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, NULL, SIZE, TRIM

!!$    ! I/O
!!$    CHARACTER(LEN=*), INTENT(IN) :: modstr

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'tracer_init'
    INTEGER                      :: n       ! set counter
    INTEGER                      :: ntrac   ! number of tracers
    CHARACTER(LEN=STRLEN_TRSET)  :: setname ! tracer set name
    LOGICAL                      :: l_init  ! initialize this set ?
    TYPE(t_trinfo_list), POINTER :: til => NULL()
    TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti  => NULL()
    !
    INTEGER                      :: status  ! error status
    LOGICAL                      :: lskip   ! skip, if no uninitialized tracers
    INTEGER                      :: iunit   ! fortran unit for input file
    LOGICAL                      :: lex     ! file exists ?
    !
    ! FOR REGRIDDING TO GLOBAL FIELD (ON I/O PE)
    REAL(DP), DIMENSION(:,:,:), POINTER :: zin => NULL()
    !
    ! GPTRSTR (DECOMPOSITION)
    !
    ! LGTRSTR
    ! ... for GP -> LG TRAFO
    REAL(DP), DIMENSION(:,:,:), POINTER :: zinl => NULL()
    REAL(DP), DIMENSION(:),     POINTER :: zxt_lg => NULL()    
    !
    ! NCREGRID
    TYPE (ncvar), DIMENSION(:), POINTER   :: var => NULL() ! list of variables
    TYPE (geohybgrid)                     :: grid  ! struct with grid info
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: dat => NULL()
    CHARACTER(len=GRD_MAXSTRLEN)          :: cevar ! object of evar
    CHARACTER(len=STRLEN_MEDIUM)          :: basename
    CHARACTER(len=STRLEN_MEDIUM)          :: subname
    INTEGER                               :: isizev ! SIZE(evar)
    INTEGER                               :: i      ! species counter
    ! 
    INTEGER, DIMENSION(:), ALLOCATABLE    :: nu, nr, ni ! counter
    INTEGER                               :: jt         ! tracer index in set

#ifdef PNCREGRID
    ! mz_kk_20080218+
    TYPE(t_mpi_def)                       :: my_mpi
    INTEGER                               :: root_pe
    INTEGER                               :: my_RG_STATUS
    ! mz_kk_20080218-
#endif

!!$    CALL start_message_bi(modstr,'TRACER INITIALISATION',substr)
    CALL start_message_bi(substr,'TRACER INITIALISATION',substr)
    
    ! CHECK IF PROCEDURE CAN BE SKIPPED, E.G. AFTER RESTART WITHOUT
    ! ADDING NEW TRACERS ...
    ! NOTE: lforce_init HAS BEEN TESTED ALREADY !!!
    !       HERE IS TEST OF linit SUFFICIENT !!!
    lskip = .true.
    set_loop1: DO n=1, NSETID

       CALL get_tracer_set(status, n, setname, trlist=til, ntrac=ntrac &
            , l_init=l_init)
       CALL tracer_halt(substr, status)

       ! NO EMPTY SETS
       IF (ntrac == 0) CYCLE

       ! INITIALISATION MUST BE ALLOWED
       IF (.NOT. l_init) CYCLE

       ! CHECK FOR UNINITIALIZED TRACERS
       DO
          IF (.NOT. ASSOCIATED(til)) EXIT
          ! i = til%info%ident%idx
          IF (.NOT. til%info%meta%linit) lskip = .false.
          til => til%next
       END DO

    END DO set_loop1

    IF (lskip) THEN
!!$       CALL end_message_bi(modstr,'TRACER INITIALISATION (SKIPPED)',substr)
       CALL end_message_bi(substr,'TRACER INITIALISATION (SKIPPED)',substr)
       RETURN
    END IF

#ifndef PNCREGRID
    IF (p_parallel_io) THEN
#endif
!!$       INQUIRE(file=TRIM(modstr)//'_t.nml', exist=lex)
       INQUIRE(file=TRIM(modstr)//'.nml', exist=lex)  ! now tracer.nml
       IF (lex) THEN
          iunit = find_next_free_unit(100,200)
       ELSE
#ifdef PNCREGRID
          IF (p_parallel_io) THEN       ! mz_kk_20080218
#endif
          CALL RGMSG(substr, RGMLW, &
!!$               'NAMELIST FILE '''//TRIM(modstr)//'_t.nml'' NOT FOUND !')
               'NAMELIST FILE '''//TRIM(modstr)//'.nml'' NOT FOUND !')
          CALL RGMSG(substr, RGMLWC, &
               'NO TRACER INITIALIZATION POSSIBLE !')
       END IF
    END IF
#ifndef PNCREGRID
    CALL p_bcast(lex, p_io)
#endif

    IF (.NOT.lex) RETURN

    ! INIT
    ALLOCATE(ni(NSETID))
    ALLOCATE(nu(NSETID))
    ALLOCATE(nr(NSETID))
    nu(:) = 0
    nr(:) = 0
    ni(:) = 0

    ! ALLOCATE SPACE FOR I/O (global field)
    IF (p_parallel_io) THEN
       ALLOCATE(zin(nlon, nlev, ngl), STAT=status)
       CALL ERRMSG(substr,status,1)
    ELSE
       NULLIFY(zin)
    ENDIF

    ! ALLOCATE SPACE FOR GP -> LG TRAFO
    ALLOCATE(zinl(nproma, nlev, ngpblks), STAT=status)
    CALL ERRMSG(substr,status,2)

    ! START REGRIDDING
    ! EXAMPLE: REGRIDDING ALL VARIABLES IN ALL NAMELISTS
    RG_CTRL = RG_PROC   ! IMMEDIATE REGRIDDING
    RG_NML  = NML_NEXT  ! READ NEXT NAMELIST FROM FILE
    !
#ifdef PNCREGRID
    My_RG_STATUS = RGSTAT_NULL ! mz_kk_20080218
    CALL INIT_PNCREGRID
#endif

    regrid_loop: DO ! ENDLESS DO LOOP (MUST BE TERMINATED WITH EXIT)

#ifdef PNCREGRID
       ! mz_kk_20080218+
       ! Deallocate only after full CYCLE. 
       ! var is calculated in parallel at first iteration/cyle
       ! at passed to the BASEMODEL if root_pe == p_pe !
       IF (My_RG_STATUS /= RGSTAT_CYCLE) THEN
       ! mz_kk_20080218-
#endif
          ! INIT
          IF (ASSOCIATED(var)) THEN
             DO i=1, SIZE(var)
                CALL INIT_NCVAR(var(i))
             END DO
             DEALLOCATE(var, STAT=status)
             CALL ERRMSG(substr,status,3)
          END IF
          NULLIFY(var)
          !
          CALL INIT_GEOHYBGRID(grid)
          !
          IF (ASSOCIATED(dat)) THEN
             DEALLOCATE(dat, STAT=status)
             CALL ERRMSG(substr,status,4)
          END IF
          NULLIFY(dat)
#ifdef PNCREGRID
       ! mz_kk_20080218+
       END IF
       
       ! pass MPI Information to NCREGRID
       my_mpi%rank  = p_pe
       my_mpi%nproc = p_nprocs
       my_mpi%comm  = p_all_comm

       CALL REGRID_CONTROL(RG_CTRL, RG_NML, My_RG_STATUS &
       ! mz_kk_20080218-
#else
       IF (p_parallel_io) THEN
       CALL REGRID_CONTROL(RG_CTRL, RG_NML, RG_STATUS &
#endif

!!$               , TRIM(modstr)//'_t.nml' &
               , TRIM(modstr)//'.nml'   &
               , iounit = iunit         &
               , var = var              &
               , grid = grid            &
#ifdef PNCREGRID
               , my_mpi  = my_mpi       & ! mz_kk_20080218
               , root_pe = root_pe      & ! mz_kk_20080218
#endif
               )
#ifndef PNCREGRID
       ENDIF
       CALL p_bcast (RG_STATUS, p_io)

       ! LEAVE ENDLESS DO LOOP AFTER LAST NAMELIST IN FILE
       IF (RG_STATUS == RGSTAT_STOP) EXIT
#else
       ! mz_kk_20080218+
       CALL p_bcast (My_RG_STATUS, root_pe)

       ! LEAVE ENDLESS DO LOOP AFTER LAST NAMELIST IN FILE
       IF (My_RG_STATUS == RGSTAT_STOP) EXIT
       ! mz_kk_20080218-
#endif

#ifndef PNCREGRID
       IF (p_parallel_io) THEN
#else
       IF (p_pe == root_pe) THEN ! mz_kk_20080218       
#endif
          isizev = SIZE(var)
          CALL RGMSG(substr, RGMLVL, '')
          CALL RGMSG(substr, RGMLVL, 'SUBROUTINE '//TRIM(substr)//' REPORT:')
          CALL RGMSG(substr, RGMLVL, &
               'FOUND ',isizev,' FIELD(S); CHECKING FOR TRACER ...')
       END IF
#ifndef PNCREGRID
       CALL p_bcast (isizev, p_io)
#else
       CALL p_bcast (isizev, root_pe)  ! mz_kk_20080218
#endif

       set_loop2: DO n=1, NSETID

          CALL get_tracer_set(status, n, setname, ti=ti, ntrac=ntrac &
               , l_init=l_init)
          CALL tracer_halt(substr, status)

          ! NO EMPTY SETS
          IF (ntrac == 0) CYCLE
          
          ! INITIALISATION MUST BE ALLOWED
          IF (.NOT. l_init) CYCLE
         
          species_loop: DO i = 1, isizev ! LOOP OVER SPECIES

#ifndef PNCREGRID
             IF (p_parallel_io) THEN
#else
             IF (p_pe == root_pe) THEN   ! mz_kk_20080218
#endif
                CALL RGTOOL_CONVERT(var(i), dat, grid, order='xzyn')
#ifndef PNCREGRID
             END IF

             IF (p_parallel_io) THEN
#endif
                cevar = TRIM(var(i)%name)
             END IF
#ifndef PNCREGRID
             CALL p_bcast(cevar, p_io)
#else
             CALL p_bcast(cevar, root_pe)  ! mz_kk_20080218
#endif

             CALL full2base_sub(status, TRIM(cevar), basename, subname)
             CALL tracer_halt(substr, status)
             !
             CALL get_tracer(status, setname               &
                  , TRIM(basename), TRIM(subname), idx=jt)
             !
             tracer_exists: IF (status == TR_NEXIST) THEN
#ifndef PNCREGRID
                IF (p_parallel_io) THEN
#else
                IF (p_pe == root_pe) THEN   ! mz_kk_20080218
#endif
                   CALL RGMSG(substr, RGMLVL, &
                        '  TRACER '''//TRIM(cevar)//&
                        &''' (SET '//TRIM(setname)//') NOT DEFINED')
                END IF
                nu(n) = nu(n) + 1
             ELSE IF (status /= 0) THEN                 ! mz_pj_20060428
                CALL tracer_halt(substr, status)        ! mz_pj_20060428
             ELSE
                ! skip if already initialised (from restart-file)
                IF (ti(jt)%tp%meta%linit) THEN
#ifndef PNCREGRID
                   IF (p_parallel_io) THEN
#else
                   IF (p_pe == root_pe) THEN   ! mz_kk_20080218
#endif
                      CALL RGMSG(substr, RGMLVL, &
                           '  TRACER '''//TRIM(var(i)%name)&
                           &//''' (SET '//TRIM(setname)//&
                           &') ALREADY INITIALIZED (e.g. FROM RESTART)')
                   END IF
                   nr(n) = nr(n) + 1
                   CYCLE ! next species
                END IF
                !
#ifndef PNCREGRID
                IF (p_parallel_io) THEN
#else
                IF (p_pe == root_pe) THEN    ! mz_kk_20080218
#endif
                   CALL RGMSG(substr, RGMLVL, &
                        '  INITIALIZING TRACER '''//TRIM(var(i)%name)//&
                        &''' (SET '//TRIM(setname)//') ')
#ifndef PNCREGRID
                   zin(:,:,:) = dat(:,:,:,1)
#else
                   ! mz_kk_20080218+
                   IF(p_nprocs == 1 .OR. p_parallel_io)   THEN
                      zin(:,:,:) = dat(:,:,:,1)
                   ELSE
                      CALL p_send(dat(:,:,:,1), p_io, jt)
                   END IF
                ELSE IF (p_parallel_io .AND. p_pe /= root_pe)   THEN 
                   IF(p_nprocs > 1)   THEN
                      CALL p_recv(zin, root_pe, jt)
                   END IF
                   ! mz_kk_20080218-
#endif
                END IF

                SELECT CASE (TRIM(setname))
                   !
                CASE(GPTRSTR)
                   ! scatter the tracer over the processors
                   CALL scatter_gp (zin, xt(:,:,jt,:), dcg)
                   IF (lresume) xtm1(:,:,jt,:) = &
                        (1._DP - eps) * xt(:,:,jt,:)
                   ! set flag that shows that tracer is already initialized
                   ti(jt)%tp%meta%linit = .TRUE.
                   !
                CASE(LGTRSTR)
                   !
                   !
                END SELECT

                ni(n) = ni(n) + 1

             END IF tracer_exists

          END DO species_loop

       END DO set_loop2

       set_loop3: DO n=1, NSETID
       
          CALL get_tracer_set(status, n, setname &
               , ntrac=ntrac, l_init=l_init)
          CALL tracer_halt(substr, status)
          IF (ntrac == 0) CYCLE
          IF (.NOT. l_init) CYCLE
          
          IF (p_parallel_io) THEN
             CALL RGMSG(substr, RGMLVL, &
                  'TRACER SET : '//TRIM(setname))
             CALL RGMSG(substr, RGMLVL, &
                  '... ',nu(n), &
                  ' UNRECOGNIZED NAMES')
             CALL RGMSG(substr, RGMLVL, &
                  '... ',nr(n), &
                  ' TRACER(S) ALREADY INITIALIZED FROM RESTART')
             CALL RGMSG(substr, RGMLVL, &
                  '... ',ni(n),' TRACER(S) INITITALIZED')
          END IF

       END DO set_loop3

       IF (p_parallel_io) THEN
          CALL RGMSG(substr, RGMLVL, 'END OF SUBROUTINE '&
               &//TRIM(substr)//' REPORT!')
          CALL RGMSG(substr, RGMLVL, '')
       END IF

    END DO regrid_loop ! ENDLESS DO LOOP
    ! END REGRIDDING

    ! CLEAN
    IF (ASSOCIATED(zin)) THEN
       DEALLOCATE(zin, STAT=status)
       CALL ERRMSG(substr,status,5)
    END IF
    NULLIFY(zin)
 
    IF (ASSOCIATED(zinl)) THEN
       DEALLOCATE(zinl, STAT=status)
       CALL ERRMSG(substr,status,6)
    END IF
    NULLIFY(zinl)
 
    DEALLOCATE(nu)
    DEALLOCATE(nr)
    DEALLOCATE(ni)
 
    CALL end_message_bi(modstr,'TRACER INITIALISATION',substr)
  END SUBROUTINE TRACER_INIT
#else
  SUBROUTINE TRACER_INIT !!$(modstr)
    IMPLICIT NONE
!!$    CHARACTER(LEN=*) :: modstr
  END SUBROUTINE TRACER_INIT
#endif
! ------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE setup_tracer_set_gp

    ! SETUP MEMORY FOR TRACERS IN GRIDPOINT REPRESENTATION

    ! ECHAM5/MESSy
    USE messy_main_data_bi, ONLY: nproma, nlev, ngpblks, l2tls

    IMPLICIT NONE

    INTRINSIC :: NULL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'setup_tracer_set_gp'
    INTEGER                        :: status
    INTEGER, DIMENSION(3)          :: dims

#if defined(ECHAM5)
    dims(1) = nproma
    dims(2) = nlev
    dims(3) = ngpblks

    ! 4 time levels: xt, xtte, xtm1, xtf (the latter is 'extended')
    ! first 3 time levels are 'standard' (t, tendency, t-1) -> .TRUE.
    ! tracers can be initialiszed -> .TRUE.
    CALL setup_tracer_set(status, GPTRSTR, dims, 4, .TRUE., .TRUE.)
    CALL tracer_halt(substr, status)
#endif
#if defined(COSMO) || defined(BLANK)
    dims(1) = nproma
    dims(2) = ngpblks
    dims(3) = nlev
    ! 5/6 levels: xt, xtte, xtm1, xtf (the latter is 'extended')
    ! first 3 time levels are 'standard' (t, tendency, t-1) -> .TRUE.

    ! A) for leapfrog (.not. l2tls)
    ! pxtf contain the memory for xtf and the 2 timelevels of the 
    ! boundary data (xt_bd), with xt_bd level 1,2 and pxtf level 3
    ! B) RUNGE_KUTTA
    ! pxtf contain the memory only for 2 timelevels of the boundary data (xt_bd)

    IF (.NOT. l2tls) THEN
       ! tracers can be initialised -> .TRUE.
       CALL setup_tracer_set(status, GPTRSTR, dims, 6, .TRUE., .TRUE.)
       CALL tracer_halt(substr, status)
    ELSE
       ! tracers can be initialised -> .TRUE.
       CALL setup_tracer_set(status, GPTRSTR, dims, 5, .TRUE., .TRUE.)
       CALL tracer_halt(substr, status)
    ENDIF
#endif
    ! SETUP POINTERS TO GRIDPOINT TRACER MEMORY AND TRACER INFORMATION
    CALL get_tracer_set(status, GPTRSTR, trlist_gp, ti_gp, ntrac_gp &
         , xt=pxt, xtte=pxtte, xtm1=pxtm1, xmem=pxtf)
    CALL tracer_halt(substr, status)

  END SUBROUTINE setup_tracer_set_gp
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE setup_tracer_set_lg

    ! SETUP MEMORY FOR TRACERS IN LAGRANGIAN REPRESENTATION

    ! ECHAM5/MESSy
    USE messy_main_data_bi, ONLY: nproma, nlev, ngpblks

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, NULL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'setup_tracer_set_lg'
    INTEGER                        :: status
    INTEGER, DIMENSION(3)          :: dims
    TYPE(t_trinfo_list), POINTER   :: til => NULL()

    dims(1) = NCELL
    dims(2) = 1
    dims(3) = 1

    ! 4 time levels: xt_a, xtte_a, xtm1_a, xtf_a (the latter is 'extended')
    ! first 3 time levels are 'standard' (t, tendency, t-1) -> .TRUE.
    ! tracers can be initialiszed -> .TRUE.
    CALL setup_tracer_set(status, LGTRSTR, dims, 4, .TRUE., .TRUE.)
    CALL tracer_halt(substr, status)

    ! SETUP POINTERS TO LAGRANGIAN TRACER MEMORY AND TRACER INFORMATION
    CALL get_tracer_set(status, LGTRSTR, trlist_lg, ti_lg, ntrac_lg &
         , xt=pxt_a, xtte=pxtte_a, xtm1=pxtm1_a, xmem=pxtf_a)
    CALL tracer_halt(substr, status)

    til => trlist_lg
    DO
      IF (.NOT. ASSOCIATED(til)) EXIT
      IF (til%info%meta%cask_i(I_mix) == ON) THEN
        number_mix = number_mix + 1
      END IF
      til => til%next
    END DO

    conversion_lg2gp: IF (l_conv_lg2gp) THEN

       ! COPY META-INFORMATION
       CALL copy_tracer_set(status, LGTRSTR, LGGPTRSTR)
       CALL tracer_halt(substr, status)

       dims(1) = nproma
       dims(2) = nlev
       dims(3) = ngpblks
       
       ! only one time level, no standard, tracers shall not be externally
       ! initialized
       CALL setup_tracer_set(status, LGGPTRSTR, dims, 1, .FALSE., .FALSE.)
       CALL tracer_halt(substr, status)

       ! SETUP POINTERS TO MEMORY AND TRACER INFORMATION
       CALL get_tracer_set(status, LGGPTRSTR, trlist_lggp, ti_lggp &
            , ntrac_lggp, xt=pxt_lggp)
       CALL tracer_halt(substr, status)

    END IF conversion_lg2gp

  END SUBROUTINE setup_tracer_set_lg
  ! -------------------------------------------------------------------

! mz_ap_20071023+
  ! -------------------------------------------------------------------
  SUBROUTINE setup_tracer_set_om
#if defined(ECHAM5)
#endif
  END SUBROUTINE setup_tracer_set_om
  ! -------------------------------------------------------------------
! mz_ap_20071023-

! mz_ab_20090624+  
  ! -------------------------------------------------------------------
  SUBROUTINE setup_tracer_set_cmat
#if defined(ECHAM5) || defined(MBM_CMAT)
#endif
  END SUBROUTINE setup_tracer_set_cmat
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE setup_tracer_set_ec
#if defined(ECHAM5)
#endif
  END SUBROUTINE setup_tracer_set_ec
  ! -------------------------------------------------------------------
! mz_ab_20090624-

  ! -------------------------------------------------------------------
  SUBROUTINE tracer_halt(substr, status)

    ! ECHAM5/MESSy
    USE messy_main_blather_bi,    ONLY: error_bi
    ! MESSy
    USE messy_main_constants_mem, ONLY: STRLEN_VLONG

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: substr
    INTEGER,          INTENT(IN)  :: status
    ! LOCAL
    CHARACTER(LEN=STRLEN_VLONG)   :: errstr

    IF (status /= 0) THEN
       errstr = tracer_error_str(status)
       CALL error_bi(errstr, substr)
    END IF

  END SUBROUTINE tracer_halt
  ! -------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_fconv_loc(direction, callstr, TRSETSTR)

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_data_bi,          ONLY: jrow, kproma
#ifndef MESSYTIMER
    USE messy_main_data_bi,          ONLY: time_step_len
#else
    USE messy_main_timer,            ONLY: time_step_len
#endif
    !
    ! MESSy
    USE messy_main_tracer_family,    ONLY: tracfamily_1_f2t, tracfamily_1_t2f &
                                         , tracfamily_2_sum, tracfamily_2_rsc 

    IMPLICIT NONE

    ! I/O
    CHARACTER(len=3), INTENT(in) :: direction
    CHARACTER(len=*), INTENT(in) :: callstr
    CHARACTER(len=*), INTENT(in) :: TRSETSTR

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_fconv_loc'
    INTEGER                     :: status

    IF (.NOT. l_family) RETURN

    SELECT CASE (direction)
       !
    CASE ('f2t','F2t','f2T','F2T')
       !
       CALL tracfamily_1_f2t(status, callstr, p_pe, TRSETSTR, &
            time_step_len, jrow, kproma)
       CALL tracer_halt(substr, status)
       !
    CASE ('t2f','T2f','t2F','T2F')
       !
       CALL tracfamily_1_t2f(status, callstr, p_pe, TRSETSTR, &
            time_step_len, jrow, kproma)
       CALL tracer_halt(substr, status)
       !
    CASE ('sum','SUM')
       !
       CALL tracfamily_2_sum(TRSETSTR, jrow)
       !
    CASE ('rsc','RSC')
       !
       CALL tracfamily_2_rsc(TRSETSTR, time_step_len, jrow)
       !
    CASE default
       !
       status = 2010 ! UNKNOWN CONVERSION FLAG
       CALL tracer_halt(substr, status)
       !
    END SELECT

  END SUBROUTINE main_tracer_fconv_loc
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_fconv_glb(direction, callstr, TRSETSTR)

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_data_bi,          ONLY: jrow, nproma, npromz &
                                         , ngpblks
#ifndef MESSYTIMER
    USE messy_main_data_bi,          ONLY: time_step_len
#else
    USE messy_main_timer,            ONLY: time_step_len
#endif
    !
    ! MESSy
    USE messy_main_tracer_family,    ONLY: tracfamily_1_f2t, tracfamily_1_t2f &
                                         , tracfamily_2_sum, tracfamily_2_rsc 


    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(len=3), INTENT(in) :: direction
    CHARACTER(len=*), INTENT(in) :: callstr
    CHARACTER(len=*), INTENT(in) :: TRSETSTR

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_fconv_glb'
    INTEGER :: status
    INTEGER :: jjrow, jp

    IF (.NOT. l_family) RETURN

    SELECT CASE (direction)
       !
    CASE ('f2t','F2t','f2T','F2T')
       !
       DO jjrow = 1, ngpblks
          IF (jjrow == ngpblks) THEN
             jp = npromz
          ELSE
             jp = nproma
          END IF
          CALL tracfamily_1_f2t(status, callstr, p_pe, TRSETSTR, &
               time_step_len, jjrow, jp)
          CALL tracer_halt(substr, status)
       END DO
       !
    CASE ('t2f','T2f','t2F','T2F')
       !
       DO jjrow = 1, ngpblks
          IF (jjrow == ngpblks) THEN
             jp = npromz
          ELSE
             jp = nproma
          END IF
          CALL tracfamily_1_t2f(status, callstr, p_pe, TRSETSTR, &
               time_step_len, jjrow, jp)
          CALL tracer_halt(substr, status)
       END DO
       !
    CASE ('sum','SUM')
       !
       DO jjrow = 1, ngpblks
          CALL tracfamily_2_sum(TRSETSTR, jjrow)
       END DO
       !
    CASE ('rsc','RSC')
       !
       DO jjrow = 1, ngpblks
          CALL tracfamily_2_rsc(TRSETSTR, time_step_len, jjrow)
       END DO
       !
    CASE default
       !
       status = 2010 ! UNKNOWN CONVERSION FLAG
       CALL tracer_halt(substr, status)
       !
    END SELECT

  END SUBROUTINE main_tracer_fconv_glb
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_write_output(flag)

   ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_data_bi,          ONLY: jrow, nproma, npromz &
                                         , ngpblks
#ifndef MESSYTIMER
    USE messy_main_data_bi,          ONLY: time_step_len
#else
    USE messy_main_timer,            ONLY: time_step_len
#endif
    !
    ! MESSy
    USE messy_main_tracer_family,    ONLY: tracfamily_1_t2f

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_write_output'
    INTEGER :: status
    INTEGER :: jjrow, jp 

    IF ( (l_family) .AND. (flag==1) ) THEN

       ! UPDATE FAMILIES TO BE CONSISTENT WITH TRACERS (FOR OUTPUT)
       DO jjrow = 1, ngpblks
          IF (jjrow == ngpblks) THEN
             jp = npromz
          ELSE
             jp = nproma
          END IF
          CALL tracfamily_1_t2f(status, substr, p_pe, GPTRSTR, &
               time_step_len, jjrow, jp, .FALSE.)
          CALL tracer_halt(substr, status)
       END DO

    END IF

  END SUBROUTINE main_tracer_write_output
  ! ----------------------------------------------------------------------


  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_read_nml_cpl(status, iou)

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

    NAMELIST /CPL/ l_conv_lg2gp, l_tracer_init

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_read_nml_cpl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    IF (l_conv_lg2gp) THEN
       WRITE(*,*) '  CONVERSION OF LG TO GP       : ON '
    ELSE
       WRITE(*,*) '  CONVERSION OF LG TO GP       : OFF'
    END IF

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_tracer_read_nml_cpl
  ! -------------------------------------------------------------------
  
! **********************************************************************+
END MODULE messy_main_tracer_bi
! **********************************************************************+
