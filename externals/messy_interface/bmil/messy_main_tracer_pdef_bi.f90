#if defined(ECHAM5)
#define MESSYCHANNELS
#endif
#ifdef COSMO
#define MESSYCHANNELS
#endif
#ifdef BLANK
#define MESSYCHANNELS
#endif

!***********************************************************************
MODULE messy_main_tracer_pdef_bi
!***********************************************************************

#if defined(ECHAM5)
  ! BML -> BMIL
#ifndef MESSYTIMER
  USE messy_main_bmluse_bi,     ONLY: time_event, io_time_event
#else
  USE messy_main_timer_event,   ONLY: time_event, io_time_event
#endif
#endif

#if defined(COSMO) || defined(BLANK)
  USE messy_main_timer_event,   ONLY: time_event, io_time_event
#endif

  ! BMIL
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  ! SMCL
  USE messy_main_tracer,        ONLY: modstr, NMAXSETID
  USE messy_main_tracer_pdef

  IMPLICIT NONE
  PRIVATE

#if defined(ECHAM5) || defined(COSMO) || defined (BLANK)
  ! TIMER
  ! ... CPL_PDEF NAMELIST VARIABLE
  TYPE(io_time_event),                     SAVE :: TPD_TIMER = &
       io_time_event(1, 'steps','first',0)
  ! ... WORKSPACE
  TYPE(time_event),                        SAVE :: XTPD_TIMER
#endif
  LOGICAL,                                 SAVE :: lnow

  ! ROUTINES
  PUBLIC :: main_tracer_pdef_initialize
  PUBLIC :: main_tracer_pdef_init_mem
  PUBLIC :: main_tracer_pdef_init_cpl
  PUBLIC :: main_tracer_pdef_global_start
  PUBLIC :: main_tracer_pdef_global_end
  PUBLIC :: main_tracer_pdef_free_mem
  !
  !PRIVATE :: main_tracer_pdef_read_nml_cpl

CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_pdef_initialize

#if defined(ECHAM5)
    ! BML -> BMIL
#ifndef MESSYTIMER
    USE messy_main_bmluse_bi, ONLY: p_bcast_event &
                                  , timer_event_init => echam_ev_init 
#else
    USE messy_main_timer_bi,  ONLY: p_bcast_event, timer_event_init 
#endif
#endif

#if defined(COSMO) || defined(BLANK)
    USE messy_main_timer_bi,  ONLY: p_bcast_event, timer_event_init 
#endif
    ! BMIL
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_pdef_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i

    ! INITIALIZE CTRL_PDEF
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL tracer_pdef_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    CALL p_bcast(L_DIAGOUT, p_io)
    DO i=1, NMAXSETID
       CALL p_bcast(TPD_DEFAULT(i)%set, p_io)
       CALL p_bcast(TPD_DEFAULT(i)%name, p_io)
       CALL p_bcast(TPD_DEFAULT(i)%subname, p_io)
       CALL p_bcast(TPD_DEFAULT(i)%lswitch, p_io)
       CALL p_bcast(TPD_DEFAULT(i)%rtol, p_io)
    END DO
    DO i=1, NMAXTPDEF 
       CALL p_bcast(TPD(i)%set, p_io)
       CALL p_bcast(TPD(i)%name, p_io)
       CALL p_bcast(TPD(i)%subname, p_io)
       CALL p_bcast(TPD(i)%lswitch, p_io)
       CALL p_bcast(TPD(i)%rtol, p_io)
    END DO

#if defined(ECHAM5) || defined(COSMO) || defined(BLANK)
    ! INITIALIZE CPL_PDEF
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_tracer_pdef_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    ! INITALIZE EVENT
    CALL p_bcast_event(TPD_TIMER, p_io)
    ! um_ak_20080708+
    !CALL echam_ev_init(XTPD_TIMER, TPD_TIMER, 'tracer_pdef', 'next')
    CALL timer_event_init(XTPD_TIMER, TPD_TIMER, 'tracer_pdef', 'next')
    ! um_ak_20080708-
#endif

  END SUBROUTINE main_tracer_pdef_initialize
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_pdef_init_mem

#ifdef MESSYCHANNELS
    ! BMIL
    USE messy_main_channel_bi,    ONLY: channel_halt, SCALAR
    ! SMCL
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute
#endif
    ! BMIL
    USE messy_main_mpi_bi,        ONLY: p_nprocs
    USE messy_main_blather_bi,    ONLY: error_bi
    ! SMCL
    USE messy_main_constants_mem, ONLY: M_air
    USE messy_main_tracer,        ONLY: NSETID, TRSET, R_molarmass &
         , param2string, AIR, AEROSOL, OCEAN, CLOUD

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_pdef_init_mem'
    INTEGER                               :: jt
    INTEGER                               :: n
    INTEGER                               :: ntrac
#ifdef MESSYCHANNELS
    INTEGER                               :: status
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: mem
#endif

    CALL start_message_bi(submodstr, 'INITIALIZE MEMORY', substr)

    CALL tracpdef_initmem(p_nprocs)

    set_loop: DO n=1, NSETID

       ntrac  = TRSET(n)%ntrac
       IF (ntrac == 0) CYCLE

#ifdef MESSYCHANNELS
       ! CHANNEL AND OBJECTS
       CALL new_channel(status, submodstr//'_'//TRIM(TRSET(n)%name) &
            , reprid=SCALAR)
       CALL channel_halt(substr, status)
#endif       

       tracer_loop: DO jt=1, ntrac

          SELECT CASE(TRSET(n)%ti(jt)%tp%ident%medium)
             !
          CASE(AIR, AEROSOL, CLOUD)
             !
             ! WHICH UNIT ?
             SELECT CASE(TRIM(TRSET(n)%ti(jt)%tp%ident%unit))
             CASE('kg/kg')
                XWRK(n)%unit(jt) = 'kg'
                XWRK(n)%scalf(jt) = 1.0_DP
             CASE('mol/mol')  ! -> kg/kg
                XWRK(n)%unit(jt) = 'kg'
                XWRK(n)%scalf(jt) = &
                     TRSET(n)%ti(jt)%tp%meta%cask_r(R_molarmass) / M_air
             CASE('1/mol')    ! -> 1/kg
                XWRK(n)%unit(jt) = '1'
                XWRK(n)%scalf(jt) = 1.0_DP / M_air
             CASE DEFAULT
                XWRK(n)%unit(jt) = &
                     'kg*('//TRIM(TRSET(n)%ti(jt)%tp%ident%unit)//')'
                XWRK(n)%scalf(jt) = 1.0_DP
                XWRK(n)%lok(jt) = .FALSE.
             END SELECT
             !
          CASE(OCEAN)
             !
             SELECT CASE(TRIM(TRSET(n)%ti(jt)%tp%ident%unit))
             CASE('kmol/m^3')  
                XWRK(n)%unit(jt) = 'kg'
                XWRK(n)%scalf(jt) = &
                     TRSET(n)%ti(jt)%tp%meta%cask_r(R_molarmass) 
             CASE DEFAULT
                XWRK(n)%unit(jt) = &
                     'kg*('//TRIM(TRSET(n)%ti(jt)%tp%ident%unit)//')'
                XWRK(n)%scalf(jt) = 1.0_DP
                XWRK(n)%lok(jt) = .FALSE.
             END SELECT
             !
          CASE DEFAULT
             !
             CALL error_bi('integration not implemented for medium '//&
                param2string(TRSET(n)%ti(jt)%tp%ident%medium, 'medium'),substr)
             !
          END SELECT

#ifdef MESSYCHANNELS
          mem => XWRK(n)%mem_mass_p(jt,:,:,:,:)
          CALL new_channel_object(status, submodstr//'_'//TRIM(TRSET(n)%name) &
               , 'MP_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname)  &
               , mem=mem)
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, submodstr//'_'//TRIM(TRSET(n)%name)     &
               , 'MP_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname) &
               , 'long_name', c='positive mass of '&
               &//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname))
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, submodstr//'_'//TRIM(TRSET(n)%name)     &
               , 'MP_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname) &
               , 'units', c=TRIM(XWRK(n)%unit(jt)))
          CALL channel_halt(substr, status)
          
          mem => XWRK(n)%mem_mass_n(jt,:,:,:,:)
          CALL new_channel_object(status, submodstr//'_'//TRIM(TRSET(n)%name) &
               , 'MN_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname)  &
               , mem=mem)
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, submodstr//'_'//TRIM(TRSET(n)%name)     &
               , 'MN_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname) &
               , 'long_name', c='negative mass of '&
               &//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname))
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, submodstr//'_'//TRIM(TRSET(n)%name)  &
               , 'MN_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname) &
               , 'units', c=TRIM(XWRK(n)%unit(jt)))
          CALL channel_halt(substr, status)

          mem => XWRK(n)%mem_mass_c(jt,:,:,:,:)
          CALL new_channel_object(status, submodstr//'_'//TRIM(TRSET(n)%name) &
               , 'MC_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname)  &
               , mem=mem)
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, submodstr//'_'//TRIM(TRSET(n)%name)  &
               , 'MC_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname) &
               , 'long_name', c='mass correction of '&
               &//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname))
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, submodstr//'_'//TRIM(TRSET(n)%name)  &
               , 'MC_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname) &
               , 'units', c=TRIM(XWRK(n)%unit(jt)))
          CALL channel_halt(substr, status)
#endif

       END DO tracer_loop

    END DO set_loop

    CALL end_message_bi(submodstr, 'INITIALIZE MEMORY',substr)

  END SUBROUTINE main_tracer_pdef_init_mem
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_pdef_init_cpl

    ! BMIL
    USE messy_main_mpi_bi,        ONLY: p_parallel_io

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_pdef_init_cpl'

    CALL start_message_bi(submodstr &
         , 'SET TRACER SWITCHES/TOLERANCES', substr)

    CALL tracpdef_settings(p_parallel_io)

    CALL end_message_bi(submodstr &
         , 'SET TRACER SWITCHES/TOLERANCES', substr)

  END SUBROUTINE main_tracer_pdef_init_cpl
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_pdef_global_start

#if defined(ECHAM5)
    ! BML -> BMIL
#ifndef MESSYTIMER
    USE messy_main_bmluse_bi,     ONLY: next_date, event_state
#else
    USE messy_main_timer_bi,      ONLY: event_state 
    USE messy_main_timer,         ONLY: next_date
#endif
#endif

#if defined(COSMO) || defined(BLANK)
    USE messy_main_timer_bi,      ONLY: event_state 
    USE messy_main_timer,         ONLY: next_date
#endif

    IMPLICIT NONE

#if defined(ECHAM5) || defined(COSMO) || defined(BLANK)
    lnow = event_state(XTPD_TIMER, next_date)
#endif
#ifdef MBM_TRACER
    lnow = .TRUE.
#endif

  END SUBROUTINE main_tracer_pdef_global_start
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_pdef_global_end

#ifdef COSMO
    USE messy_main_data_bi,          ONLY: l2tls, istartpar,jstartpar &
                                         , iendpar, jendpar
#endif    

    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_pe, p_bcast, p_nprocs &
                                         , p_parallel_io
    USE messy_main_blather_bi,       ONLY: error_bi
#ifndef MESSYTIMER
    USE messy_main_data_bi,          ONLY: time_step_len
#else
    USE messy_main_timer,            ONLY: time_step_len
#endif
    ! SMCL
    USE messy_main_constants_mem,    ONLY: g
    USE messy_main_tracer,           ONLY: NSETID, TRSET

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_pdef_global_end'
    INTEGER :: status
    INTEGER :: jc
    INTEGER :: n
    INTEGER :: jrow, jk
    INTEGER :: js

    ! INIT
    DO n=1, NSETID
       XWRK(n)%mass_n(:) = 0.0_DP
       XWRK(n)%mass_p(:) = 0.0_DP
       XWRK(n)%mass_c(:) = 0.0_DP
    END DO

    IF (.NOT. lnow) RETURN

    ! 1st step: sum on each PE
#ifdef COSMO
    ! um_ak_20110330 l2tls replaced by TRUE (correct always xt = nnew)
    CALL tracpdef_integrate(status, 1, time_step_len, p_pe, .TRUE. &
        , ia1=istartpar,ie1=iendpar, ia2=jstartpar,ie2=jendpar)
#else
    CALL tracpdef_integrate(status, 1, time_step_len, p_pe, .FALSE.)
#endif
    ! status always 0 ...

    ! BROADCAST RESULT TO ALL PEs
    DO n=1, NSETID
       IF (TRSET(n)%ntrac == 0) CYCLE
       DO jc = 0, p_nprocs-1
#if defined(COSMO) || defined(BLANK)
          DO js=1,SIZE(XWRK(n)%mass_pe(:,:,jc),2)
             CALL p_bcast(XWRK(n)%mass_pe(:,js,jc), jc)
          END DO
#else
          CALL p_bcast(XWRK(n)%mass_pe(:,:,jc), jc)
#endif
       END DO
    END DO

    ! 2nd step: sum over all PEs and check tolerance criterium
    CALL tracpdef_integrate(status, 2, time_step_len, p_pe)
    IF (status /= 0) THEN
       CALL error_bi('negative mass of tracer exceeds tolerance', substr)
    END IF

    ! diagnostic output
    CALL tracpdef_print(p_parallel_io)

  END SUBROUTINE main_tracer_pdef_global_end
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_pdef_free_mem

    IMPLICIT NONE

    CALL tracpdef_freemem

  END SUBROUTINE main_tracer_pdef_free_mem
! ----------------------------------------------------------------------

#if defined(ECHAM5) || defined(COSMO) || defined(BLANK)
! **********************************************************************
! PRIVATE INTERFACE ROUTINES
! **********************************************************************

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_pdef_read_nml_cpl(status, iou)

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_pdef_read_nml_cpl'

    NAMELIST /CPL_PDEF/ TPD_TIMER

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL_PDEF', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL_PDEF, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL_PDEF', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_tracer_pdef_read_nml_cpl
! ----------------------------------------------------------------------
#endif
!***********************************************************************
END MODULE messy_main_tracer_pdef_bi
!***********************************************************************
