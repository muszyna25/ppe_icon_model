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

  USE mo_run_config, ONLY: config_ldump_states    => ldump_states,    &
                         & config_lrestore_states => lrestore_states, &
                         & config_ltestcase       => ltestcase,       &
                         & config_ldynamics       => ldynamics,       &
                         & config_iforcing        => iforcing,        &
                         & config_ltransport      => ltransport,      &
                         & config_ntracer         => ntracer,         &
                         & config_lvert_nest      => lvert_nest,      &
                         & config_nlev            => nlev,            &
                         & config_num_lev         => num_lev,         &
                         & config_nshift          => nshift,          &
                         & config_nsteps          => nsteps,          &
                         & config_dtime           => dtime,           &
                         & config_ltimer          => ltimer,          &
                         & config_timers_level    => timers_level,    &
                         & config_msg_level       => msg_level,       &
                         & config_inextra_2d      => inextra_2d,      &
                         & config_inextra_3d      => inextra_3d

  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: finish
  USE mo_impl_constants, ONLY: max_dom, max_ntracer, inoforcing, IHELDSUAREZ, &
                               INWP,IECHAM,ILDF_ECHAM,IMPIOM,INOFORCING,ILDF_DRY
  USE mo_io_units,       ONLY: nnml, nnml_output
  USE mo_namelist,       ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,            ONLY: my_process_is_stdio 
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

  LOGICAL :: ldump_states    ! Dump patch/interpolation/grid refinement state of every
                             ! patch (after subdivision in case of a parallel run)
                             ! end exit program
  LOGICAL :: lrestore_states ! Restore patch/interpolation/grid refinement states
                             ! from dump files instead of calculating them

  LOGICAL :: ltestcase       ! if .TRUE. then
                             ! - compute analytical initial state,
                             !   depending on the specified test case,
                             ! - compute analytical boundary conditions,
                             ! - if applicable, compute analytical forcing

  LOGICAL :: ldynamics       ! if .TRUE., switch on adiabatic dynamics
  INTEGER :: iforcing        ! adiabatic forcing

  LOGICAL :: ltransport      ! if .TRUE., switch on large-scale tracer transport
  INTEGER :: ntracer         ! number of advected tracers

  INTEGER :: nlev               ! number of full levels = number of layers
  LOGICAL :: lvert_nest         ! if .TRUE., switch on vertical nesting
  INTEGER :: num_lev(max_dom)   ! number of full levels for each domain
  INTEGER :: nshift (max_dom)   ! half level of parent domain which coincides 
                                    ! with the upper boundary of the current domain jg

  INTEGER  :: nsteps            ! number of time steps
  REAL(wp) :: dtime             ! [s] length of a time step

  LOGICAL :: ltimer        ! if .TRUE., wallclock timers are switched on
  INTEGER :: timers_level  ! what level of timers to run

  INTEGER :: msg_level     ! how much printout is generated during runtime

  INTEGER :: inextra_2d    ! number of extra output fields for debugging
  INTEGER :: inextra_3d    ! number of extra output fields for debugging


  NAMELIST /run_nml/ ldump_states, lrestore_states, &
                     ltestcase,    ldynamics,       &
                     iforcing,     ltransport,      &
                     ntracer,                       &
                     lvert_nest,   nlev,            &
                     num_lev,      nshift,          &
                     nsteps,       dtime,           &
                     ltimer,       timers_level,    &
                     msg_level,                     &
                     inextra_2d,   inextra_3d

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
    ldump_states    = .FALSE.
    lrestore_states = .FALSE.

    ltestcase       = .TRUE.
    ldynamics       = .TRUE.
    iforcing        = inoforcing

    ltransport      = .FALSE.
    ntracer         = 0

    lvert_nest = .FALSE. ! no vertical nesting
    nlev       = 31
    num_lev(:) = nlev  ! number of full levels for each domain
    nshift(:)  = 0         ! please do not change the default.
                               ! otherwise the initialization of 
                               ! p_patch(jg)%nshift in "import patches" 
                               ! will not work properly.

    nsteps = 0
    dtime  = 600._wp   ! [s] for R2B04 + semi-implicit time steppping

    ltimer       = .TRUE.
    timers_level = 1
    msg_level    = 10

    inextra_2d   = 0           !> no extra output 2D fields
    inextra_3d   = 0           !> no extra output 3D fields

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('run_nml')
      READ(funit,NML=run_nml)
      CALL close_tmpfile(funit)
    END IF

    !----------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !----------------------------------------------------------------------
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
    SELECT CASE (iforcing)                                                     
    CASE(INOFORCING,IHELDSUAREZ,INWP,IECHAM,ILDF_DRY,ILDF_ECHAM,IMPIOM)
    CASE DEFAULT
      CALL finish( TRIM(routine),'wrong value for iforcing')
    END SELECT

    IF (ltransport.AND.ntracer<1) CALL finish(TRIM(routine), &
    'Tracer transport is switched on, but number of advected tracers is smaller than 1')

    IF (ANY(num_lev < 0)) CALL finish(TRIM(routine),'"num_lev" must be positive')
    IF (ANY(nshift  < 0)) CALL finish(TRIM(routine),'"nshift" must be positive')

    IF (nsteps < 0) CALL finish(TRIM(routine),'"nsteps" must not be negative')
    IF (dtime <= 0._wp) CALL finish(TRIM(routine),'"dtime" must be positive')

    IF ((ntracer<0).OR.(ntracer>max_ntracer)) CALL finish( TRIM(routine), &
    'wrong number of tracers. Valid range: 0<= ntracer <=20')

    !----------------------------------------------------
    ! Fill part of the configuration state
    !----------------------------------------------------
    config_ldump_states    = ldump_states
    config_lrestore_states = lrestore_states

    config_ltestcase       = ltestcase 
    config_ldynamics       = ldynamics 
    config_iforcing        = iforcing 

    config_ltransport      = ltransport 
    config_ntracer         = ntracer 

    config_lvert_nest      = lvert_nest
    config_nlev            = nlev
    config_num_lev(:)      = num_lev(:)
    config_nshift(:)       = nshift(:)

    config_nsteps          = nsteps  
    config_dtime           = dtime 

    config_ltimer          = ltimer
    config_timers_level    = timers_level
    config_msg_level       = msg_level

    config_inextra_2d      = inextra_2d
    config_inextra_3d      = inextra_3d

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=run_nml)
    CALL store_and_close_namelist(funit, 'run_nml')

    ! write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=run_nml)

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
