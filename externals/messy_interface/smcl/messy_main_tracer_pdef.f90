! ************************************************************************
MODULE messy_main_tracer_pdef
! ************************************************************************

  ! MODULE FOR FORCING POSITIVE TRACERS (MESSy-SMCL)
  ! 
  ! Authors:
  !    Patrick Joeckel, MPICH, March 2007
  !

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP, STRLEN_MEDIUM
  USE messy_main_tracer,        ONLY: NMAXSETID, NSETID, TRSET, STRLEN_TRSET &
                                    , modstr

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  PUBLIC :: DP

  PUBLIC :: tracpdef_initmem
  PUBLIC :: tracpdef_settings
  PUBLIC :: tracpdef_airmass
  PUBLIC :: tracpdef_integrate
  PUBLIC :: tracpdef_freemem
  PUBLIC :: tracpdef_print
  !
  PUBLIC :: tracer_pdef_read_nml_ctrl

  INTERFACE tracpdef_airmass
     MODULE PROCEDURE tracpdef_airmass_const
     MODULE PROCEDURE tracpdef_airmass_3d
  END INTERFACE

  ! ----------- <

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodstr = 'tracer_pdef'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodver = '1.3'

  TYPE IO_TRACPDEF
     CHARACTER(len=STRLEN_TRSET)  :: set=''      ! name of tracer set
     CHARACTER(len=STRLEN_MEDIUM) :: name=''     ! name of tracer family
     CHARACTER(len=STRLEN_MEDIUM) :: subname=''  ! subname of tracer family
     LOGICAL, DIMENSION(2)        :: lswitch = (/ .FALSE., .FALSE. /)
     REAL(DP)                     :: rtol = 1.0_dp
  END TYPE IO_TRACPDEF
  PUBLIC :: IO_TRACPDEF

  ! MAX NUMBER OF NAMELIST ENTRIES
  INTEGER, PARAMETER, PUBLIC :: NMAXTPDEF = 500 * NMAXSETID  ! 500 per set

  TYPE T_SET_WORKSPACE
     LOGICAL,  DIMENSION(:,:),       POINTER :: xlswitch => NULL()
     REAL(DP), DIMENSION(:),         POINTER :: xrtol    => NULL()
     !
     LOGICAL                                 :: ldiagonly
     !
     REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: mem_mass_p => NULL()
     REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: mem_mass_n => NULL()
     REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: mem_mass_c => NULL()
     REAL(DP), DIMENSION(:),         POINTER :: mass_p => NULL()
     REAL(DP), DIMENSION(:),         POINTER :: mass_n => NULL()
     REAL(DP), DIMENSION(:),         POINTER :: mass_c => NULL()
     ! ...
     REAL(DP), DIMENSION(:,:,:),     POINTER :: airmass => NULL()
     ! ...
     REAL(DP), DIMENSION(:,:,:),     POINTER :: mass_pe => NULL()
     LOGICAL,  DIMENSION(:),         POINTER :: lok     => NULL()
     REAL(DP), DIMENSION(:),         POINTER :: scalf   => NULL()
     CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: unit => NULL()
  END TYPE T_SET_WORKSPACE
  PUBLIC :: T_SET_WORKSPACE

  ! CTRL-NAMELIST
  LOGICAL,                                 PUBLIC, SAVE :: L_DIAGOUT = .FALSE.
  TYPE(IO_TRACPDEF), DIMENSION(NMAXSETID), PUBLIC, SAVE :: TPD_DEFAULT
  TYPE(IO_TRACPDEF), DIMENSION(NMAXTPDEF), PUBLIC, SAVE :: TPD

  ! GLOBAL VARIABLES
  TYPE(T_SET_WORKSPACE), DIMENSION(:), POINTER, PUBLIC, SAVE :: XWRK => NULL()

CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_initmem(nprocs)
    
    IMPLICIT NONE
    INTRINSIC :: SIZE, ASSOCIATED

    ! I/O
    INTEGER, INTENT(IN) :: nprocs

    ! LOCAL
    INTEGER :: n, s1, s2, s4, ntrac

    ! INIT
    ALLOCATE(XWRK(NSETID))

    set_loop: DO n=1, NSETID

       ! NUMBER OF TRACERS IN THIS SET
       ntrac = TRSET(n)%ntrac

       ! IF ONE OF TENDENCY OR VALUE OF T-1 IS NOT PRESENT
       ! MASS CAN ONLY BE DIAGNOSED BUT NOT CORRECTED
       XWRK(n)%ldiagonly = (.NOT. ASSOCIATED(TRSET(n)%xtte)) .OR. &
            (.NOT. ASSOCIATED(TRSET(n)%xtm1))
       
       ALLOCATE(XWRK(n)%xlswitch(ntrac,2))
       XWRK(n)%xlswitch(:,:) = .FALSE.       ! default

       ALLOCATE(XWRK(n)%xrtol(ntrac))
       XWRK(n)%xrtol(:) = 1.0_dp             ! default
       !
       ALLOCATE(XWRK(n)%mem_mass_n(ntrac,1,1,1,1))
       ALLOCATE(XWRK(n)%mem_mass_p(ntrac,1,1,1,1))
       ALLOCATE(XWRK(n)%mem_mass_c(ntrac,1,1,1,1))
       !
       XWRK(n)%mass_n => XWRK(n)%mem_mass_n(:,1,1,1,1)
       XWRK(n)%mass_p => XWRK(n)%mem_mass_p(:,1,1,1,1)
       XWRK(n)%mass_c => XWRK(n)%mem_mass_c(:,1,1,1,1)

       s1 = SIZE(TRSET(n)%xt, 1)
       s2 = SIZE(TRSET(n)%xt, 2)
       !s3 = SIZE(TRSET(n)%xt, 3) ! = ntrac
       s4 = SIZE(TRSET(n)%xt, 4)
       !s5 = SIZE(TRSET(n)%xt, 5) ! = 1
       ALLOCATE(XWRK(n)%airmass(s1, s2, s4))
       XWRK(n)%airmass(:,:,:) = 0.0_dp

       ALLOCATE(XWRK(n)%mass_pe(ntrac, 3, 0:nprocs-1))
       XWRK(n)%mass_pe(:,:,:) = 0.0_DP
       
       ALLOCATE(XWRK(n)%lok(ntrac))
       XWRK(n)%lok(:) = .TRUE.
       
       ALLOCATE(XWRK(n)%scalf(ntrac))
       XWRK(n)%scalf(:) = 0.0_DP

       ALLOCATE(XWRK(n)%unit(ntrac))
       XWRK(n)%unit(:) = ''

    END DO set_loop
    
  END SUBROUTINE tracpdef_initmem
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_settings(ldiagout)

    USE messy_main_tracer, ONLY: get_tracer

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    LOGICAL, INTENT(IN) :: ldiagout

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tracpdef_settings'
    INTEGER :: i, ierr, jt
    INTEGER :: n, nn
    INTEGER :: ntrac

    set_loop: DO n=1, NSETID

       ntrac = TRSET(n)%ntrac

       IF (ntrac == 0) CYCLE

       ! INIT: SET DEFAULT
       default_loop: DO nn=1, NMAXSETID
          IF (TRIM(TPD_DEFAULT(nn)%set) == '') CYCLE
          IF (TRIM(TPD_DEFAULT(nn)%set) /= TRIM(TRSET(n)%name)) CYCLE
          XWRK(n)%xlswitch(:,1) = TPD_DEFAULT(nn)%lswitch(1)
          XWRK(n)%xlswitch(:,2) = TPD_DEFAULT(nn)%lswitch(2)
          XWRK(n)%xrtol(:)      = TPD_DEFAULT(nn)%rtol
          IF (ldiagout) THEN
             WRITE(*,*) 'DEFAULT FOR SET '//TRIM(TRSET(n)%name)//': ' &
                  ,TPD_DEFAULT(nn)%lswitch(:), TPD_DEFAULT(nn)%rtol
          END IF
       END DO default_loop
    
       ! LOOP AND SET SPECIAL
       DO i=1, NMAXTPDEF
          IF (TRIM(TPD(i)%set)  == '') CYCLE
          IF (TRIM(TPD(i)%name) == '') CYCLE
          IF (TRIM(TPD(i)%set) /= TRIM(TRSET(n)%name)) CYCLE

          CALL get_tracer(ierr, TRIM(TPD(i)%set), TRIM(TPD(i)%name) &
               , subname=TRIM(TPD(i)%subname), idx=jt)
          IF (ierr == 0) THEN
             XWRK(n)%xlswitch(jt,:) = TPD(i)%lswitch(:)
             XWRK(n)%xrtol(jt)      = TPD(i)%rtol
          ELSE
             IF (ldiagout) THEN
                IF (TRIM(TPD(i)%subname) == '') THEN
                   WRITE(*,*) 'TRACER '''//TRIM(TPD(i)%name)//&
                        &''' not found in set '//TRIM(TPD(i)%set)//&
                        &' ... skipping ...'
                ELSE
                   WRITE(*,*) 'TRACER '''//TRIM(TPD(i)%name)//&
                        &'_'//TRIM(TPD(i)%subname)//&
                        &''' not found in set '//TRIM(TPD(i)%set)//&
                        &' ... skipping ...'
                END IF
             END IF
          END IF
       END DO

       IF (ldiagout) THEN
          WRITE(*,*) '======================================================='
          WRITE(*,*) 'TRACER SET: ',TRIM(TRSET(n)%name)
          WRITE(*,*) '-------------------------------------------------------'
          DO jt=1, ntrac
             WRITE(*,*) jt, TRSET(n)%ti(jt)%tp%ident%fullname &
                  , XWRK(n)%xlswitch(jt,:), XWRK(n)%xrtol(jt)
          END DO
          WRITE(*,*) '======================================================='
       END IF

    END DO set_loop

  END SUBROUTINE tracpdef_settings
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_airmass_3d(setname, airmass)

    IMPLICIT NONE
    INTRINSIC :: TRIM
    
    ! I/O
    CHARACTER(LEN=*),           INTENT(IN)  :: setname
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: airmass

    ! LOCAL
    LOGICAL :: l_found
    INTEGER :: n

    l_found = .FALSE.
    set_loop: DO n=1, NSETID
       IF (TRIM(TRSET(n)%name) == TRIM(setname)) THEN
          l_found = .TRUE.
          EXIT
       END IF
    END DO set_loop

    IF (.NOT. l_found) RETURN

    IF (TRSET(n)%ntrac == 0) RETURN

    XWRK(n)%airmass(:,:,:) = airmass(:,:,:)

  END SUBROUTINE tracpdef_airmass_3d
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_airmass_const(setname, airmass)

    IMPLICIT NONE
    INTRINSIC :: TRIM
    
    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: setname
    REAL(DP),         INTENT(IN)  :: airmass

    ! LOCAL
    LOGICAL :: l_found
    INTEGER :: n

    l_found = .FALSE.
    set_loop: DO n=1, NSETID
       IF (TRIM(TRSET(n)%name) == TRIM(setname)) THEN
          l_found = .TRUE.
          EXIT
       END IF
    END DO set_loop

    IF (.NOT. l_found) RETURN

    IF (TRSET(n)%ntrac == 0) RETURN

    XWRK(n)%airmass(:,:,:) = airmass

  END SUBROUTINE tracpdef_airmass_const
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! um_ak_20090317+
! um_ak_20090910+
!  SUBROUTINE tracpdef_integrate(status, flag, time_step_len, p_pe)
! um_ak-20110330 replace l2tls by lnewts (different meaning) correct
!                always xt for COSMO
  SUBROUTINE tracpdef_integrate(status, flag, time_step_len, p_pe, lnewtl &
                               ,ia1,ie1,ia2,ie2,ia3,ie3)
! um_ak_20090910-
! um_ak_20090317-

    USE messy_main_constants_mem, ONLY: FLAGGED_BAD, TINY_DP

    IMPLICIT NONE

    ! I/O
    INTEGER,   INTENT(OUT)      :: status
    INTEGER,   INTENT(IN)       :: flag
    REAL(DP),  INTENT(IN)       :: time_step_len
    INTEGER,   INTENT(IN)       :: p_pe
    ! um_ak_20110330 LOGICAL,   INTENT(IN), OPTIONAL :: l2tls ! um_ak_20090317
    LOGICAL,   INTENT(IN), OPTIONAL :: lnewtl ! um_ak_20110330
    INTEGER,   INTENT(IN), OPTIONAL :: ia1,ie1,ia2,ie2,ia3,ie3 ! um_ak_20090910

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tracpdef_integerate'
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zxt
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zxtte
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zxtp1
    REAL(DP), DIMENSION(:),     ALLOCATABLE :: ratio
    REAL(DP)                                :: div
    INTEGER                                 :: jt
    INTEGER                                 :: n
    INTEGER                                 :: s1, s2, s4
    INTEGER                                 :: ntrac
    ! um_ak_20110330 LOGICAL                :: l_2tls = .FALSE. !um_ak_20090317
    LOGICAL                                 :: l_newtl = .FALSE.! um_ak_20110330
    INTEGER                                 :: a1,e1,a2,e2,a3,e3

    INTRINSIC :: MERGE, SUM, MAX, ABS, TINY, TRIM, SIZE

    status = 0

    ! um_ak_20110330 IF (PRESENT(l2tls)) l_2tls = l2tls ! um_ak_20090317
    IF (PRESENT(lnewtl)) l_newtl = lnewtl ! um_ak_20110330

    SELECT CASE(flag)

    CASE(1)

       set_loop1: DO n=1, NSETID
          
          s1 = SIZE(TRSET(n)%xt(:,:,:,:,:), 1)
          s2 = SIZE(TRSET(n)%xt(:,:,:,:,:), 2)
          !s3 = SIZE(TRSET(n)%xt(:,:,:,:,:), 3) ! = ntrac
          s4 = SIZE(TRSET(n)%xt(:,:,:,:,:), 4)
          !s5 = SIZE(TRSET(n)%xt(:,:,:,:,:), 5) ! = 1
          ntrac = TRSET(n)%ntrac
          
          IF (ntrac == 0) CYCLE

          ! um_ak_20090910+
          a1 = 1
          a2 = 1
          a3 = 1
          e1 = s1
          e2 = s2
          e3 = s4
          IF (TRSET(n)%name == 'gp') THEN  ! um_ak_20090910 sorry, unfortunately
             ! I see no other way than to do this hardcoded (otherwise we have
             ! to do the tracer set loop in bi)
             IF (PRESENT(ia1)) a1 =ia1
             IF (PRESENT(ia2)) a2 =ia2
             IF (PRESENT(ia3)) a3 =ia3
             IF (PRESENT(ie1)) e1 =ie1
             IF (PRESENT(ie2)) e2 =ie2
             IF (PRESENT(ie3)) e3 =ie3
          ENDIF
          ! um_ak_20090910+

          ! INIT
          ALLOCATE(zxt(s1,s2,s4))
          ALLOCATE(zxtte(s1,s2,s4))
          ALLOCATE(zxtp1(s1,s2,s4))

          diag_only: IF (XWRK(n)%ldiagonly) THEN

             tracer_loop1: DO jt=1, ntrac

                IF (.NOT.XWRK(n)%lok(jt)) CYCLE

                ! INTEGRATE MASS (NEGATIVE)
                zxt(a1:e1,a2:e2,a3:e3) = &
                     MERGE(TRSET(n)%xt(a1:e1,a2:e2,jt,a3:e3,1), &
                     0.0_DP, TRSET(n)%xt(a1:e1,a2:e2,jt,a3:e3,1)<0.0_DP)
                zxt(a1:e1,a2:e2,a3:e3) = MERGE(0.0_dp, zxt(a1:e1,a2:e2,a3:e3), &
                     ABS(zxt(a1:e1,a2:e2,a3:e3)-FLAGGED_BAD)<TINY_DP)
                XWRK(n)%mass_pe(jt, 1, p_pe) = &
                     SUM( zxt(a1:e1,a2:e2,a3:e3) *          &
                     XWRK(n)%airmass(a1:e1,a2:e2,a3:e3) ) * &
                     XWRK(n)%scalf(jt)
                
                ! INTEGRATE MASS (POSITIVE)
                zxt(a1:e1,a2:e2,a3:e3) = &
                     MERGE(TRSET(n)%xt(a1:e1,a2:e2,jt,a3:e3,1), &
                     0.0_DP, TRSET(n)%xt(a1:e1,a2:e2,jt,a3:e3,1)>0.0_DP)
                zxt(a1:e1,a2:e2,a3:e3) = MERGE(0.0_dp, zxt(a1:e1,a2:e2,a3:e3), &
                     ABS(zxt(a1:e1,a2:e2,a3:e3)-FLAGGED_BAD)<TINY_DP)
                XWRK(n)%mass_pe(jt, 2, p_pe) =              &
                     SUM( zxt(a1:e1,a2:e2,a3:e3) *          &
                     XWRK(n)%airmass(a1:e1,a2:e2,a3:e3) ) * &
                     XWRK(n)%scalf(jt)
                
                ! NO CORRECTION POSSIBLE
                XWRK(n)%mass_pe(jt, 3, p_pe) = 0.0_DP

             END DO tracer_loop1

          ELSE
          
             tracer_loop2: DO jt=1, ntrac
             
                IF (.NOT.XWRK(n)%lok(jt)) CYCLE

                ! um_ak_20110330 if2tls: IF (.NOT. l_2tls) THEN ! um_ak_20090317
                ifnewtl: IF (.NOT. l_newtl) THEN ! um_ak_20110330
                   ! ACTUAL VALUE AT t+1
                   zxtp1(a1:e1,a2:e2,a3:e3) = &
                        TRSET(n)%xtm1(a1:e1,a2:e2,jt,a3:e3,1) + &
                        TRSET(n)%xtte(a1:e1,a2:e2,jt,a3:e3,1)*time_step_len

                   ! PSEUDO-TENDENCY TO FORCE RESULTS AT t+1 >= 0.0
                   zxtte(a1:e1,a2:e2,a3:e3) = &
                        MAX(-zxtp1(a1:e1,a2:e2,a3:e3)/time_step_len, 0.0_DP)

                   ! INTEGRATE MASS (NEGATIVE)
                   zxt(a1:e1,a2:e2,a3:e3) = MERGE(zxtp1(a1:e1,a2:e2,a3:e3) &
                        , 0.0_DP, zxtp1(a1:e1,a2:e2,a3:e3)<0.0_DP)
                   XWRK(n)%mass_pe(jt, 1, p_pe) =              &
                        SUM( zxt(a1:e1,a2:e2,a3:e3) *          &
                        XWRK(n)%airmass(a1:e1,a2:e2,a3:e3) ) * &
                        XWRK(n)%scalf(jt)
                   ! INTEGRATE MASS (POSITIVE)
                   zxt(a1:e1,a2:e2,a3:e3) = MERGE(zxtp1(a1:e1,a2:e2,a3:e3) &
                        , 0.0_DP, zxtp1(a1:e1,a2:e2,a3:e3)>0.0_DP)
                   XWRK(n)%mass_pe(jt, 2, p_pe) =              &
                        SUM( zxt(a1:e1,a2:e2,a3:e3) *          &
                        XWRK(n)%airmass(a1:e1,a2:e2,a3:e3) ) * &
                        XWRK(n)%scalf(jt)

                   ! RESET NEGATIVES ON REQUEST
                   XWRK(n)%mass_pe(jt, 3, p_pe) = 0.0_DP
                   IF (XWRK(n)%xlswitch(jt,1)) THEN
                      TRSET(n)%xtte(a1:e1,a2:e2,jt,a3:e3,1) =      &
                           TRSET(n)%xtte(a1:e1,a2:e2,jt,a3:e3,1) + &
                           zxtte(a1:e1,a2:e2,a3:e3)
                      XWRK(n)%mass_pe(jt, 3, p_pe) =              &
                           SUM(zxtte(a1:e1,a2:e2,a3:e3) *         &
                           XWRK(n)%airmass(a1:e1,a2:e2,a3:e3) ) * &
                           XWRK(n)%scalf(jt) * time_step_len
                   END IF
                ELSE ! l_newtl ! um_ak_20090317+

                   ! ACTUAL VALUE AT t+1
                   zxtp1(a1:e1,a2:e2,a3:e3) = &
                        TRSET(n)%xt(a1:e1,a2:e2,jt,a3:e3,1)

! um_ak_20110331+
!!$                   ! PSEUDO-TENDENCY TO FORCE RESULTS AT t+1 >= 0.0
!!$                   zxtte(a1:e1,a2:e2,a3:e3) = &
!!$                        MAX(-zxtp1(a1:e1,a2:e2,a3:e3)/time_step_len, 0.0_DP)
! um_ak_20110331-
                   ! INTEGRATE MASS (NEGATIVE)
                   zxt(a1:e1,a2:e2,a3:e3) = MERGE(zxtp1(a1:e1,a2:e2,a3:e3) &
                        , 0.0_DP,zxtp1(a1:e1,a2:e2,a3:e3)<0.0_DP)
                   XWRK(n)%mass_pe(jt, 1, p_pe) =              &
                        SUM( zxt(a1:e1,a2:e2,a3:e3) *          &
                        XWRK(n)%airmass(a1:e1,a2:e2,a3:e3) ) * &
                        XWRK(n)%scalf(jt)

                   ! INTEGRATE MASS (POSITIVE)
                   zxt(a1:e1,a2:e2,a3:e3) = MERGE(zxtp1(a1:e1,a2:e2,a3:e3) &
                        , 0.0_DP,zxtp1(a1:e1,a2:e2,a3:e3)>0.0_DP)
                   XWRK(n)%mass_pe(jt, 2, p_pe) =              &
                        SUM( zxt(a1:e1,a2:e2,a3:e3) *          &
                        XWRK(n)%airmass(a1:e1,a2:e2,a3:e3) ) * &
                        XWRK(n)%scalf(jt)

                   ! RESET NEGATIVES ON REQUEST
                   XWRK(n)%mass_pe(jt, 3, p_pe) = 0.0_DP
                   IF (XWRK(n)%xlswitch(jt,1)) THEN
! um_ak_20111121+
!!$                      TRSET(n)%xt(a1:e1,a2:e2,jt,a3:e3,1) = &
!!$                           zxt(a1:e1,a2:e2,a3:e3)
                      TRSET(n)%xt(:,:,jt,:,1) = &
                           MERGE(TRSET(n)%xt(:,:,jt,:,1) &
                        , 0.0_DP,TRSET(n)%xt(:,:,jt,:,1)>0.0_DP)
! um_ak_20111121-
                      XWRK(n)%mass_pe(jt, 3, p_pe) =              &
! um_ak_20110331+
                           XWRK(n)%mass_pe(jt, 1, p_pe)
!!$                           SUM(zxtte(a1:e1,a2:e2,a3:e3) *         &
!!$                           XWRK(n)%airmass(a1:e1,a2:e2,a3:e3) ) * &
!!$                           XWRK(n)%scalf(jt) * time_step_len
! um_ak_20110331-
                   END IF

                ENDIF ifnewtl ! um_ak_20090317-
             END DO tracer_loop2

          END IF diag_only

          ! CLEAN UP
          DEALLOCATE(zxt)
          DEALLOCATE(zxtte)
          DEALLOCATE(zxtp1)

       END DO set_loop1

    CASE(2)

       set_loop2: DO n=1, NSETID

          ntrac = TRSET(n)%ntrac

          IF (ntrac == 0) CYCLE

          ! INIT
          ALLOCATE(ratio(ntrac))
          ratio(:) = 0.0_DP

          ! INTEGRATE MASSES
          XWRK(n)%mass_n(:)  = SUM(XWRK(n)%mass_pe(:,1,:),2) ! SUM OVER p_pe
          XWRK(n)%mass_p(:)  = SUM(XWRK(n)%mass_pe(:,2,:),2) ! SUM OVER p_pe
          XWRK(n)%mass_c(:)  = SUM(XWRK(n)%mass_pe(:,3,:),2) ! SUM OVER p_pe

          ! CALCULATE RATIO
          DO jt=1, ntrac
             IF (.NOT.XWRK(n)%lok(jt)) THEN
                ratio(jt) = -777.0_DP
                CYCLE
             END IF
             !
             ! CONSERVED MASS
             div = ABS(XWRK(n)%mass_n(jt) + XWRK(n)%mass_p(jt))
             IF (div > TINY(0._DP)) THEN
                ratio(jt) = ABS(XWRK(n)%mass_n(jt))/div
             ELSE
                ratio(jt) = -999.0_DP
             END IF
          END DO

          ! APPLY TOLERANCE CRITERIUM
          tracer_loop3: DO jt=1, ntrac
             IF (.NOT.XWRK(n)%xlswitch(jt,2)) CYCLE 
             IF (.NOT.XWRK(n)%lok(jt)) CYCLE

             IF (ratio(jt) > XWRK(n)%xrtol(jt)) THEN
                
                WRITE(*,*) substr,' :NEGATIVE MASS OF TRACER '''//&
                     &TRIM(TRSET(n)%ti(jt)%tp%ident%fullname)//&
                     &''' IN SET '''//TRIM(TRSET(n)%name)//&
                     &''' EXCEEDS TOLERANCE: '            &
                     , XWRK(n)%mass_n(jt),XWRK(n)%mass_p(jt),ratio(jt)
                status = 1 ! terminate simulation
                RETURN
             END IF
          END DO tracer_loop3

          DEALLOCATE(ratio)

       END DO set_loop2

    END SELECT

  END SUBROUTINE tracpdef_integrate
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_freemem

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: n

    set_loop: DO n=1, NSETID

       DEALLOCATE(XWRK(n)%xlswitch)
       DEALLOCATE(XWRK(n)%xrtol)

       DEALLOCATE(XWRK(n)%mem_mass_n)
       DEALLOCATE(XWRK(n)%mem_mass_p)
       DEALLOCATE(XWRK(n)%mem_mass_c)
       
       DEALLOCATE(XWRK(n)%airmass)

       DEALLOCATE(XWRK(n)%mass_pe)
       DEALLOCATE(XWRK(n)%lok)
       DEALLOCATE(XWRK(n)%scalf)
       DEALLOCATE(XWRK(n)%unit)

    END DO set_loop

    DEALLOCATE(XWRK)

  END SUBROUTINE tracpdef_freemem
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_print(ldiagout)

    IMPLICIT NONE
    INTRINSIC :: SUM, TRIM

    ! I/O
    LOGICAL, INTENT(IN) :: ldiagout

    ! LOCAL
    INTEGER :: n, ntrac, jt

    IF (.NOT. L_DIAGOUT) RETURN
    IF (.NOT. ldiagout)  RETURN
    
    set_loop: DO n=1, NSETID

       ntrac = TRSET(n)%ntrac
       IF (ntrac == 0) CYCLE

       WRITE(*,*) '========================================================'
       WRITE(*,*) 'TRACER SET: '//TRIM(TRSET(n)%name)
       WRITE(*,*) 'AIR MASS  : ',SUM(XWRK(n)%airmass)
       WRITE(*,*) '--------------------------------------------------------'
       WRITE(*,'(a15,1x,3(a12,1x))') 'TRACER', 'MP', 'MN', 'MC'
       WRITE(*,*) '--------------------------------------------------------'

       tracer_loop1: DO jt=1, ntrac

          IF (.NOT.XWRK(n)%lok(jt)) CYCLE

          WRITE(*,'(a15,1x,3(e12.4,1x),a)') &
               TRSET(n)%ti(jt)%tp%ident%fullname  &
               , XWRK(n)%mass_p(jt), XWRK(n)%mass_n(jt), XWRK(n)%mass_c(jt) &
               , TRIM(XWRK(n)%unit(jt))

       END DO tracer_loop1

       WRITE(*,*) '========================================================'

    END DO set_loop

  END SUBROUTINE tracpdef_print
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracer_pdef_read_nml_ctrl(status, iou)

    ! Author: Patrick Joeckel, MPICH, Aug 2005

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'tracer_pdef_read_nml_ctrl'

    NAMELIST /CTRL_PDEF/ L_DIAGOUT, TPD_DEFAULT, TPD

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CTRL_PDEF', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_PDEF, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_PDEF', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE tracer_pdef_read_nml_ctrl
! ----------------------------------------------------------------------

! ************************************************************************
END MODULE messy_main_tracer_pdef
! ************************************************************************
