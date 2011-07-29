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
                         & config_msg_level       => msg_level

  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: finish
  USE mo_impl_constants, ONLY: max_dom, max_ntracer, inoforcing, IHELDSUAREZ, &
                               INWP,IECHAM,ILDF_ECHAM,IMPIOM,INOFORCING,ILDF_DRY
  USE mo_io_units,       ONLY: nnml, nnml_output
  USE mo_namelist,       ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,            ONLY: my_process_is_stdio 
  USE mo_master_control, ONLY: is_restart_run

  USE mo_io_restart_namelist,   ONLY: open_tmpfile, store_and_close_namelist,   &
                                    & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_run_namelist

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

  LOGICAL :: lvert_nest         ! if .TRUE., switch on vertical nesting
  INTEGER :: num_lev(max_dom)   ! number of full levels for each domain
  INTEGER :: nshift (max_dom)   ! half level of parent domain which coincides 
                                    ! with the upper boundary of the current domain jg

  INTEGER  :: nsteps            ! number of time steps
  REAL(wp) :: dtime             ! [s] length of a time step

  LOGICAL :: ltimer        ! if .TRUE., wallclock timers are switched on
  INTEGER :: timers_level  ! what level of timers to run

  INTEGER :: msg_level     ! how much printout is generated during runtime


  NAMELIST /run_nml/ ldump_states, lrestore_states, &
                     ltestcase,    ldynamics,       &
                     iforcing,     ltransport,      &
                     ntracer,                       &
                     lvert_nest,                    &
                     num_lev,      nshift,          &
                     nsteps,       dtime,           &
                     ltimer,       timers_level,    &
                     msg_level

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
    num_lev(:) = 31    ! number of full levels for each domain
    nshift(:)  = 0       ! please do not change the default.
                         ! otherwise the initialization of 
                         ! p_patch(jg)%nshift in "import patches" 
                         ! will not work properly.

    nsteps = 0
    dtime  = 600._wp     ! [s] for R2B04 + semi-implicit time steppping

    ltimer       = .TRUE.
    timers_level = 1
    msg_level    = 10

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
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

    ! If we choose to have NWP-forcing for the nonhydrostatic model, we want 
    ! to avoid the necessity of setting ntracer explicitly. Thus, a sanity 
    ! check for ntracer is triggered only, if iforcing /= INWP.
    IF (iforcing /= INWP) THEN
      IF ((ntracer<0).OR.(ntracer>max_ntracer)) CALL finish( TRIM(routine), &
      'wrong number of tracers. Valid range: 0<= ntracer <=20')

      IF (ltransport .AND. ntracer<1) CALL finish(TRIM(routine), &
      'Tracer transport is switched on, but number of advected tracers is smaller than 1')
    ENDIF

    IF (ANY(num_lev < 0)) CALL finish(TRIM(routine),'"num_lev" must be positive')
    IF (ANY(nshift  < 0)) CALL finish(TRIM(routine),'"nshift" must be positive')

    IF (nsteps < 0) CALL finish(TRIM(routine),'"nsteps" must not be negative')
    IF (dtime <= 0._wp) CALL finish(TRIM(routine),'"dtime" must be positive')


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
    config_nlev            = num_lev(1)
    config_num_lev(:)      = num_lev(:)
    config_nshift(:)       = nshift(:)

    config_nsteps          = nsteps  
    config_dtime           = dtime 

    config_ltimer          = ltimer
    config_timers_level    = timers_level
    config_msg_level       = msg_level

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=run_nml)
      CALL store_and_close_namelist(funit, 'run_nml')
    ENDIF
    ! write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=run_nml)

  END SUBROUTINE read_run_namelist
  !-------------

END MODULE mo_run_nml
