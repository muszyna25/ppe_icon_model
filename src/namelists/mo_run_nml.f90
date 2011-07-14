!>
!!        
!! @par Revision History
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!
MODULE mo_run_nml

  USE mo_run_config
  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: finish
  USE mo_impl_constants, ONLY: max_dom, max_ntracer, inoforcing, IHELDSUAREZ, &
                               INWP,IECHAM,ILDF_ECHAM,IMPIOM,INOFORCING,ILDF_DRY
  USE mo_io_units,       ONLY: nnml, nnml_output
  USE mo_namelist,       ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,            ONLY: p_pe, p_io
  USE mo_master_nml,     ONLY: lrestart

  USE mo_io_restart_attributes, ONLY: get_restart_attribute
  USE mo_io_restart_namelist,   ONLY: open_tmpfile, store_and_close_namelist,   &
                                    & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_run_namelist
  PUBLIC :: run_nml_setup, locean    !!!to be removed

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !------------------------------------------------------------------------
  ! Namelist variables
  !------------------------------------------------------------------------

  LOGICAL :: nml_ldump_states    ! Dump patch/interpolation/grid refinement state of every
                                 ! patch (after subdivision in case of a parallel run)
                                 ! end exit program
  LOGICAL :: nml_lrestore_states ! Restore patch/interpolation/grid refinement states
                                 ! from dump files instead of calculating them

  LOGICAL :: nml_ltestcase       ! if .TRUE. then
                                 ! - compute analytical initial state,
                                 !   depending on the specified test case,
                                 ! - compute analytical boundary conditions,
                                 ! - if applicable, compute analytical forcing

  LOGICAL :: nml_ldynamics       ! if .TRUE., switch on adiabatic dynamics
  INTEGER :: nml_iforcing        ! adiabatic forcing

  LOGICAL :: nml_ltransport      ! if .TRUE., switch on large-scale tracer transport
  INTEGER :: nml_ntracer         ! number of advected tracers

  INTEGER :: nml_nlev               ! number of full levels = number of layers
  LOGICAL :: nml_lvert_nest         ! if .TRUE., switch on vertical nesting
  INTEGER :: nml_num_lev(max_dom)   ! number of full levels for each domain
  INTEGER :: nml_nshift (max_dom)   ! half level of parent domain which coincides 
                                    ! with the upper boundary of the current domain jg

  INTEGER  :: nml_nsteps            ! number of time steps
  REAL(wp) :: nml_dtime             ! [s] length of a time step

  LOGICAL :: nml_ltimer        ! if .TRUE., wallclock timers are switched on
  INTEGER :: nml_timers_level  ! what level of timers to run

  INTEGER :: nml_msg_level     ! how much printout is generated during runtime

  INTEGER :: nml_inextra_2d    !> number of extra output fields for debugging
  INTEGER :: nml_inextra_3d    !> number of extra output fields for debugging


  NAMELIST /run_nml/ nml_ldump_states, nml_lrestore_states, &
                     nml_ltestcase,    nml_ldynamics,       &
                     nml_iforcing,     nml_ltransport,      &
                     nml_ntracer,                           &
                     nml_lvert_nest,   nml_nlev,            &
                     nml_num_lev,      nml_nshift,          &
                     nml_nsteps,       nml_dtime,           &
                     nml_ltimer,       nml_timers_level,    &
                     nml_msg_level,                         &
                     nml_inextra_2d,   nml_inextra_3d

  !-----------------------------------------------------------------------
  ! 2.0 Declaration of dependent control variables 
  !-----------------------------------------------------------------------
  LOGICAL  :: locean        ! if .TRUE., the model runs in ocean mode

CONTAINS
  !>
  !!
  SUBROUTINE read_run_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    CHARACTER(len=*), PARAMETER :: routine = 'mo_run_nml:read_run_namelist'

    !------------------------------------------------------------
    ! Default settings
    !------------------------------------------------------------
    nml_ldump_states    = .FALSE.
    nml_lrestore_states = .FALSE.

    nml_ltestcase       = .TRUE.
    nml_ldynamics       = .TRUE.
    nml_iforcing        = inoforcing

    nml_ltransport      = .FALSE.
    nml_ntracer         = 0

    nml_lvert_nest = .FALSE. ! no vertical nesting
    nml_nlev       = 31
    nml_num_lev(:) = nml_nlev  ! number of full levels for each domain
    nml_nshift(:)  = 0         ! please do not change the default.
                               ! otherwise the initialization of 
                               ! p_patch(jg)%nshift in "import patches" 
                               ! will not work properly.

    nml_nsteps = 0
    nml_dtime  = 600._wp   ! [s] for R2B04 + semi-implicit time steppping

    nml_ltimer       = .TRUE.
    nml_timers_level = 1
    nml_msg_level    = 10

    nml_inextra_2d   = 0           !> no extra output 2D fields
    nml_inextra_3d   = 0           !> no extra output 3D fields

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('run_nml')
      READ(funit,NML=run_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processors)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml('run_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, run_nml)
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! Sanity check
    !----------------------------------------------------
    SELECT CASE (nml_iforcing)                                                     
    CASE(INOFORCING,IHELDSUAREZ,INWP,IECHAM,ILDF_DRY,ILDF_ECHAM,IMPIOM)
    CASE DEFAULT
      CALL finish( TRIM(routine),'wrong value for nml_iforcing')
    END SELECT

    IF (nml_ltransport.AND.nml_ntracer<1) CALL finish(TRIM(routine), &
    'Tracer transport is switched on, but number of advected tracers is smaller than 1')

    IF (ANY(nml_num_lev < 0)) CALL finish(TRIM(routine),'"nml_num_lev" must be positive')
    IF (ANY(nml_nshift  < 0)) CALL finish(TRIM(routine),'"nml_nshift" must be positive')

    IF (nml_nsteps < 0) CALL finish(TRIM(routine),'"nsteps" must not be negative')
    IF (nml_dtime <= 0._wp) CALL finish(TRIM(routine),'"dtime" must be positive')

    IF ((nml_ntracer<0).OR.(nml_ntracer>max_ntracer)) CALL finish( TRIM(routine), &
    'wrong number of tracers. Valid range: 0<= ntracer <=20')

    !----------------------------------------------------
    ! Fill part of the configuration state
    !----------------------------------------------------
    ldump_states    = nml_ldump_states
    lrestore_states = nml_lrestore_states

    ltestcase       = nml_ltestcase 
    ldynamics       = nml_ldynamics 
    iforcing        = nml_iforcing 

    ltransport      = nml_ltransport 
    ntracer         = nml_ntracer 

    lvert_nest      = nml_lvert_nest
    nlev            = nml_nlev
    num_lev(:)      = nml_num_lev(:)
    nshift(:)       = nml_nshift(:)

    nsteps          = nml_nsteps  
    dtime           = nml_dtime 

    ltimer          = nml_ltimer
    timers_level    = nml_timers_level
    msg_level       = nml_msg_level

    inextra_2d      = nml_inextra_2d
    inextra_3d      = nml_inextra_3d

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=run_nml)
    CALL store_and_close_namelist(funit, 'run_nml')

    ! write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=run_nml)

  END SUBROUTINE read_run_namelist
  !-------------
  !>
  !!
  SUBROUTINE run_nml_setup
    !---------------------------------------------------------------
    ! 4. Check whether the namelist varibles have reasonable values
    !---------------------------------------------------------------

!     CASE (ihs_ocean)
!       locean         = .TRUE.
!     CASE default
!       CALL finish( TRIM(routine),'wrong equation specifier iequations')
!     END SELECT
! 
! 
  
    ! 5.3 End date and time, and length of integration
    !     Here we define "nsteps" as the number of time steps THIS integration
    !     will last, regardless of the restart status.

!     IF (nsteps/=0) THEN
! 
!     
!       length_sec   = REAL(nsteps,wp)*dtime
!       end_datetime = current_datetime
!       CALL add_time(length_sec,0,0,0,end_datetime)
!       !
!     ELSE IF (run_day/=0 .OR. run_hour/=0 .OR. run_minute/=0 .OR. run_second/=0.0_wp) THEN
!       IF (run_day    < 0    ) CALL finish(routine,'"run_day" must not be negative')
!       IF (run_hour   < 0    ) CALL finish(routine,'"run_hour" must not be negative')
!       IF (run_minute < 0    ) CALL finish(routine,'"run_minute" must not be negative')
!       IF (run_second < 0._wp) CALL finish(routine,'"run_second" must not be negative')
!       !
!       end_datetime = current_datetime
!       CALL add_time(run_second,run_minute,run_hour,run_day,end_datetime)
!       !
!       cur_datetime_calsec = (REAL(current_datetime%calday,wp)+current_datetime%caltime) &
!         &                   *REAL(current_datetime%daylen,wp)
!       end_datetime_calsec = (REAL(end_datetime%calday,wp)+end_datetime%caltime) &
!         &                   *REAL(end_datetime%daylen,wp)
!       nsteps=INT((end_datetime_calsec-cur_datetime_calsec)/dtime)
!       !
!     ELSE
!       ! compute nsteps from current_datetime, end_datetime and dtime
!       end_datetime%calendar = calendar
!       end_datetime%year     = end_year
!       end_datetime%month    = end_month
!       end_datetime%day      = end_day
!       end_datetime%hour     = end_hour
!       end_datetime%minute   = end_minute
!       end_datetime%second   = end_second
!       CALL date_to_time      (end_datetime) ! fill date time structure
!       !
!       cur_datetime_calsec = (REAL(current_datetime%calday,wp)+current_datetime%caltime) &
!         &                   *REAL(current_datetime%daylen,wp)
!       end_datetime_calsec = (REAL(end_datetime%calday,wp)+end_datetime%caltime) &
!         &                   *REAL(end_datetime%daylen,wp)
!       IF (end_datetime_calsec < cur_datetime_calsec) &
!         & CALL finish(routine,'The end date and time must not be '// &
!         &            'before the current date and time')
!       !
!       nsteps=INT((end_datetime_calsec-cur_datetime_calsec)/dtime)
!       !
!     END IF
! 
!     nsteps = MIN(nsteps,INT(dt_restart/dtime))

!   IF ( (dynamics_config(1)%iequations == IHS_ATM_TEMP) .OR. &
!        (dynamics_config(1)%iequations == IHS_ATM_THETA)     ) THEN

!    ! If running the HYDROSTATIC version,
!    ! let the model integrate one more step after the desired end of
!    ! simulation in order to get the proper output. This additional step is
!    ! necessary because the HYDROSTATIC model writes out values of step N
!    ! after the integration from N to N+1 is finished. Also note that
!    ! this additional step is done only for the regular output, and is 
!    ! ignored for restart.

!    nsteps = nsteps + 1

!    ! The additional step is not needed in the NON-hydrostatic version because
!    ! in this case the model writes out values of step N
!    ! after the integration from N-1 to N is finished.
!   ENDIF
    !..................................................................
    !
!     CALL message(' ',' ')
!     CALL message(routine,'End date and time')
!     CALL message(routine,'-----------------')
! !     CALL print_datetime_all(end_datetime)  ! print all date and time components
! 
!     CALL message(' ',' ')
!     CALL message(routine,'Length of restart cycle')
!     CALL message(routine,'-----------------------')
! !     WRITE(message_text,'(a,f10.2,a,f16.10,a)') &
! !          &'dt_restart :',dt_restart,' seconds =', dt_restart/86400._wp, ' days'
! !     CALL message(routine,message_text)
! 
!     CALL message(' ',' ')
!     CALL message(routine,'Length of this run')
!     CALL message(routine,'------------------')
!     WRITE (message_text,'(a,f7.2)') 'dtime [s] :',dtime
!     CALL message(routine,message_text)
!     WRITE (message_text,'(a,i7)')   'nsteps    :',nsteps
!     CALL message(routine,message_text)
!     CALL message(' ',' ')
    
 END SUBROUTINE run_nml_setup

END MODULE mo_run_nml
